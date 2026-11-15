//! Scalar reference kernel over a pair list.
//!
//! This is the correctness reference for the SIMD and cluster kernels. It uses
//! only the formulas in [`crate::pairwise`].

use crate::pairwise::{coulomb_ewald_pair, excluded_ewald_pair, lj_pair, lorentz_berthelot};
use crate::types::{min_image, AtomParams, Energies, ForceBuffers, NonbondedConfig, Positions};

/// Minimum-image separation `r_i - r_j` for one pair, using the config's box.
#[inline]
pub(crate) fn separation(
    pos: &Positions<'_>,
    i: usize,
    j: usize,
    box_len: Option<[f64; 3]>,
) -> [f64; 3] {
    let mut d = [pos.x[i] - pos.x[j], pos.y[i] - pos.y[j], pos.z[i] - pos.z[j]];
    if let Some(l) = box_len {
        for (c, edge) in d.iter_mut().zip(l) {
            *c = min_image(*c, edge);
        }
    }
    d
}

/// Evaluates LJ and real-space Ewald Coulomb for every pair in `pairs`.
///
/// `pairs` holds `(i, j)` atom indices with `i != j`. Each pair is taken once; the
/// caller decides whether it lists both orders. Forces are *added* to `forces`
/// (equal and opposite on `i` and `j`). Pairs at or beyond `cfg.cutoff` are skipped.
pub fn pair_list_scalar(
    pos: Positions<'_>,
    params: AtomParams<'_>,
    pairs: &[(u32, u32)],
    cfg: &NonbondedConfig,
    forces: &mut ForceBuffers<'_>,
) -> Energies {
    let rc2 = cfg.cutoff * cfg.cutoff;
    let mut en = Energies::default();
    for &(i, j) in pairs {
        let (i, j) = (i as usize, j as usize);
        let d = separation(&pos, i, j, cfg.box_len);
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if r2 >= rc2 || r2 == 0.0 {
            continue;
        }
        let (sigma, eps) = lorentz_berthelot(
            params.sigma[i],
            params.epsilon[i],
            params.sigma[j],
            params.epsilon[j],
        );
        let (e_lj, f_lj) = lj_pair(r2, sigma, eps);
        let qq = params.charge[i] * params.charge[j];
        let (e_c, f_c) = coulomb_ewald_pair(r2, qq, cfg.alpha);
        en.lj += e_lj;
        en.coulomb += e_c;
        let f = f_lj + f_c;
        apply_pair_force(forces, i, j, f, d);
    }
    en
}

/// Adds `f_over_r * d` to atom `i` and subtracts it from atom `j`.
#[inline]
pub(crate) fn apply_pair_force(
    forces: &mut ForceBuffers<'_>,
    i: usize,
    j: usize,
    f_over_r: f64,
    d: [f64; 3],
) {
    let fx = f_over_r * d[0];
    let fy = f_over_r * d[1];
    let fz = f_over_r * d[2];
    forces.fx[i] += fx;
    forces.fy[i] += fy;
    forces.fz[i] += fz;
    forces.fx[j] -= fx;
    forces.fy[j] -= fy;
    forces.fz[j] -= fz;
}

/// Excluded-pair PME correction over a list of excluded `(i, j)` pairs.
///
/// Adds `-qq erf(alpha r)/r` for each pair and returns that energy. Uses the same
/// minimum-image separation as the real-space kernel.
pub fn excluded_pairs_correction(
    pos: Positions<'_>,
    charge: &[f64],
    excluded: &[(u32, u32)],
    alpha: f64,
    box_len: Option<[f64; 3]>,
    forces: &mut ForceBuffers<'_>,
) -> f64 {
    let mut energy = 0.0;
    for &(i, j) in excluded {
        let (i, j) = (i as usize, j as usize);
        let d = separation(&pos, i, j, box_len);
        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if r2 == 0.0 {
            continue;
        }
        let (e, f_over_r) = excluded_ewald_pair(r2, charge[i] * charge[j], alpha);
        energy += e;
        apply_pair_force(forces, i, j, f_over_r, d);
    }
    energy
}
