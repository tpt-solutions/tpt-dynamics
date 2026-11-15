//! Padded 4-atom cluster layout and the cluster-pair kernel.
//!
//! Atoms are grouped into clusters of [`CLUSTER_SIZE`] consecutive entries (the
//! simulation may permute atoms first). The last cluster is padded with
//! [`PAD`] indices and zero parameters. Each cluster stores its coordinates and
//! parameters as `[f64; 4]` lanes, so one cluster pair is a 4x4 block of atom
//! interactions. The kernel vectorises over the four atoms of the partner cluster.
//!
//! The cluster-pair list is built by [`build_cluster_pairs`]. That builder is a
//! brute-force O(clusters^2) bounding-sphere test, intended for tests and small
//! systems; production neighbour lists belong to `tpt-dynamics-sim`.

use tpt_simd::{F64x4, SimdMask};

use crate::simd_kernel::{hsum, lane_interaction, min_image_v, LANES};
use crate::types::{AtomParams, Energies, ForceBuffers, NonbondedConfig, Positions};

/// Atoms per cluster.
pub const CLUSTER_SIZE: usize = 4;
/// Atom index marking a padding lane.
pub const PAD: u32 = u32::MAX;

/// Padded structure-of-arrays cluster storage.
#[derive(Clone, Debug, PartialEq)]
pub struct ClusterSet {
    /// Original atom index per lane, or [`PAD`].
    pub atoms: Vec<[u32; CLUSTER_SIZE]>,
    /// x coordinates per lane (zero for padding).
    pub x: Vec<[f64; CLUSTER_SIZE]>,
    /// y coordinates per lane.
    pub y: Vec<[f64; CLUSTER_SIZE]>,
    /// z coordinates per lane.
    pub z: Vec<[f64; CLUSTER_SIZE]>,
    /// Lennard-Jones sigma per lane (zero for padding).
    pub sigma: Vec<[f64; CLUSTER_SIZE]>,
    /// Lennard-Jones epsilon per lane (zero for padding).
    pub epsilon: Vec<[f64; CLUSTER_SIZE]>,
    /// Charge per lane (zero for padding).
    pub charge: Vec<[f64; CLUSTER_SIZE]>,
}

/// Force accumulators in the same cluster-lane layout as [`ClusterSet`].
#[derive(Clone, Debug, PartialEq)]
pub struct ClusterForces {
    /// x force per lane.
    pub fx: Vec<[f64; CLUSTER_SIZE]>,
    /// y force per lane.
    pub fy: Vec<[f64; CLUSTER_SIZE]>,
    /// z force per lane.
    pub fz: Vec<[f64; CLUSTER_SIZE]>,
}

impl ClusterSet {
    /// Groups `pos.len()` atoms into consecutive clusters of four, padding the last.
    pub fn from_atoms(pos: Positions<'_>, params: AtomParams<'_>) -> Self {
        let n = pos.len();
        let nc = n.div_ceil(CLUSTER_SIZE);
        let mut set = Self {
            atoms: vec![[PAD; CLUSTER_SIZE]; nc],
            x: vec![[0.0; CLUSTER_SIZE]; nc],
            y: vec![[0.0; CLUSTER_SIZE]; nc],
            z: vec![[0.0; CLUSTER_SIZE]; nc],
            sigma: vec![[0.0; CLUSTER_SIZE]; nc],
            epsilon: vec![[0.0; CLUSTER_SIZE]; nc],
            charge: vec![[0.0; CLUSTER_SIZE]; nc],
        };
        for a in 0..n {
            let (c, l) = (a / CLUSTER_SIZE, a % CLUSTER_SIZE);
            set.atoms[c][l] = a as u32;
            set.x[c][l] = pos.x[a];
            set.y[c][l] = pos.y[a];
            set.z[c][l] = pos.z[a];
            set.sigma[c][l] = params.sigma[a];
            set.epsilon[c][l] = params.epsilon[a];
            set.charge[c][l] = params.charge[a];
        }
        set
    }

    /// Number of clusters.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// True when the set has no clusters.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Zeroed force accumulators matching this layout.
    pub fn zero_forces(&self) -> ClusterForces {
        let nc = self.len();
        ClusterForces {
            fx: vec![[0.0; CLUSTER_SIZE]; nc],
            fy: vec![[0.0; CLUSTER_SIZE]; nc],
            fz: vec![[0.0; CLUSTER_SIZE]; nc],
        }
    }

    /// Adds cluster-layout forces into per-atom force buffers, skipping padding.
    pub fn scatter_forces(&self, cf: &ClusterForces, forces: &mut ForceBuffers<'_>) {
        for c in 0..self.len() {
            for l in 0..CLUSTER_SIZE {
                let a = self.atoms[c][l];
                if a != PAD {
                    let a = a as usize;
                    forces.fx[a] += cf.fx[c][l];
                    forces.fy[a] += cf.fy[c][l];
                    forces.fz[a] += cf.fz[c][l];
                }
            }
        }
    }

    /// Bounding sphere `(center, radius)` of the valid atoms in cluster `c`.
    fn bounding_sphere(&self, c: usize) -> ([f64; 3], f64) {
        let mut n = 0.0;
        let mut ctr = [0.0; 3];
        for l in 0..CLUSTER_SIZE {
            if self.atoms[c][l] != PAD {
                ctr[0] += self.x[c][l];
                ctr[1] += self.y[c][l];
                ctr[2] += self.z[c][l];
                n += 1.0;
            }
        }
        if n == 0.0 {
            return (ctr, 0.0);
        }
        for v in &mut ctr {
            *v /= n;
        }
        let mut radius: f64 = 0.0;
        for l in 0..CLUSTER_SIZE {
            if self.atoms[c][l] != PAD {
                let d = [self.x[c][l] - ctr[0], self.y[c][l] - ctr[1], self.z[c][l] - ctr[2]];
                radius = radius.max((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());
            }
        }
        (ctr, radius)
    }
}

/// Builds the cluster-pair list `(ci, cj)` with `ci <= cj` that may hold an atom pair
/// within `cutoff`, using bounding spheres (with minimum image when `box_len` is set).
///
/// Brute force O(clusters^2). Intended for tests and small systems.
pub fn build_cluster_pairs(
    set: &ClusterSet,
    cutoff: f64,
    box_len: Option<[f64; 3]>,
) -> Vec<(u32, u32)> {
    let spheres: Vec<_> = (0..set.len()).map(|c| set.bounding_sphere(c)).collect();
    let mut pairs = Vec::new();
    for ci in 0..set.len() {
        for cj in ci..set.len() {
            let (c1, r1) = spheres[ci];
            let (c2, r2) = spheres[cj];
            let mut d2 = 0.0;
            for k in 0..3 {
                let mut d = c1[k] - c2[k];
                if let Some(l) = box_len {
                    d -= l[k] * (d / l[k]).round();
                }
                d2 += d * d;
            }
            let reach = cutoff + r1 + r2;
            if d2 <= reach * reach {
                pairs.push((ci as u32, cj as u32));
            }
        }
    }
    pairs
}

/// SIMD cluster-pair kernel.
///
/// For each `(ci, cj)` with `ci <= cj`, every valid atom of `ci` interacts with the
/// four lanes of `cj` in one vector operation. When `ci == cj` only lanes `j > i`
/// are used, so each intra-cluster pair is counted once. Forces are *added* to
/// `out` in cluster layout; call [`ClusterSet::scatter_forces`] to merge them.
pub fn cluster_pair_kernel(
    set: &ClusterSet,
    pairs: &[(u32, u32)],
    cfg: &NonbondedConfig,
    out: &mut ClusterForces,
) -> Energies {
    let rc2v = F64x4::splat(cfg.cutoff * cfg.cutoff);
    let zero = F64x4::splat(0.0);
    let mut en = Energies::default();

    for &(ci, cj) in pairs {
        let (ci, cj) = (ci as usize, cj as usize);
        let xj = F64x4::from(set.x[cj]);
        let yj = F64x4::from(set.y[cj]);
        let zj = F64x4::from(set.z[cj]);
        let sj = F64x4::from(set.sigma[cj]);
        let ej = F64x4::from(set.epsilon[cj]);
        let qj = F64x4::from(set.charge[cj]);
        let valid_j: [bool; LANES] = std::array::from_fn(|l| set.atoms[cj][l] != PAD);

        for il in 0..CLUSTER_SIZE {
            if set.atoms[ci][il] == PAD {
                continue;
            }
            let mut dx = F64x4::splat(set.x[ci][il]) - xj;
            let mut dy = F64x4::splat(set.y[ci][il]) - yj;
            let mut dz = F64x4::splat(set.z[ci][il]) - zj;
            if let Some(l) = cfg.box_len {
                dx = min_image_v(dx, l[0]);
                dy = min_image_v(dy, l[1]);
                dz = min_image_v(dz, l[2]);
            }
            let r2 = dx * dx + dy * dy + dz * dz;

            let ok: [bool; LANES] =
                std::array::from_fn(|jl| valid_j[jl] && (ci < cj || (ci == cj && jl > il)));
            let mask = r2.simd_lt(rc2v) & r2.simd_gt(zero) & SimdMask::<f64, LANES>::from_array(ok);

            // Lorentz-Berthelot combination per lane (padding lanes have zero params).
            let sig = (F64x4::splat(set.sigma[ci][il]) + sj) * F64x4::splat(0.5);
            let eps = (F64x4::splat(set.epsilon[ci][il]) * ej).sqrt();
            let qq = F64x4::splat(set.charge[ci][il]) * qj;

            let (e_lj, e_c, f) = lane_interaction(r2, mask, sig, eps, qq, cfg.alpha);
            en.lj += hsum(e_lj);
            en.coulomb += hsum(e_c);

            let fxa: [f64; LANES] = (f * dx).into();
            let fya: [f64; LANES] = (f * dy).into();
            let fza: [f64; LANES] = (f * dz).into();
            let (mut sx, mut sy, mut sz) = (0.0, 0.0, 0.0);
            for jl in 0..LANES {
                sx += fxa[jl];
                sy += fya[jl];
                sz += fza[jl];
                out.fx[cj][jl] -= fxa[jl];
                out.fy[cj][jl] -= fya[jl];
                out.fz[cj][jl] -= fza[jl];
            }
            out.fx[ci][il] += sx;
            out.fy[ci][il] += sy;
            out.fz[ci][il] += sz;
        }
    }
    en
}
