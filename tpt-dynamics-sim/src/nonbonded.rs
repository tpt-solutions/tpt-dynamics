//! Non-bonded interactions with exclusions, scaled 1-4 pairs and optional
//! particle-mesh Ewald (PME) electrostatics.
//!
//! Terms, for a pair at distance `r`:
//!
//! * Lennard-Jones: the [`LennardJones`] potential (cutoff with optional shift
//!   or switch), applied to every neighbour-list pair that is not excluded.
//! * Coulomb, real space (PME): `f·q_i q_j erfc(α r)/r` inside the cutoff.
//! * Coulomb, reciprocal (PME): the smooth-PME sum from `tpt-dynamics-compute`
//!   on a mesh whose grid spacing is at most 0.12 nm, plus the self-energy
//!   `−f α/√π Σ q²` and the neutralising-background term `−f π Q²/(2 α² V)`
//!   (for net charge `Q`).
//! * Exclusions: pairs listed in the topology are removed from the list. In PME
//!   runs their reciprocal-space contribution is removed with the correction
//!   `−f q_i q_j erf(α r)/r`.
//! * 1-4 pairs: LJ scaled by `fudge_lj` and Coulomb by `fudge_qq`, both without
//!   a cutoff (the GROMACS convention).
//!
//! Here `f` is the Coulomb factor: `138.935458 kJ mol⁻¹ nm e⁻²` in MD units and
//! `1` in reduced units. The PME solver works in units where `f = 1`, so its
//! outputs are multiplied by `f`.
//!
//! PME requires an orthorhombic box; [`Nonbonded::validate_box`] rejects
//! triclinic boxes when PME is enabled.

use std::sync::{Mutex, PoisonError};

use tpt_dynamics_compute::pairwise::{
    coulomb_ewald_pair, excluded_ewald_pair, lj_pair, lorentz_berthelot,
};
use tpt_dynamics_compute::pme::{alpha_for_tolerance, Pme, PmeParams};
use tpt_dynamics_compute::types::{ForceBuffers, Positions};
use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::forcefield::{ForceField, ForceOutput, LennardJones, LjType, Modifier};
use crate::neighbour::VerletList;
use crate::pbc::SimBox;

/// Coulomb factor `e²/(4πε₀)` in kJ mol⁻¹ nm e⁻² (CODATA 2018 value used by GROMACS).
pub const COULOMB_FACTOR_MD: f64 = 138.935_458;

/// Electrostatic treatment.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Electrostatics {
    /// No Coulomb interaction.
    None,
    /// Ewald real space plus smooth PME reciprocal space.
    Pme {
        /// Target real-space accuracy `erfc(α r_c) = tolerance`.
        tolerance: f64,
        /// B-spline interpolation order (2..=12; 4 is typical).
        order: usize,
    },
}

/// Settings for a [`Nonbonded`] force field.
#[derive(Clone, Debug, PartialEq)]
pub struct NonbondedSettings {
    /// Lennard-Jones parameters per atom type.
    pub lj_types: Vec<LjType>,
    /// Cutoff in nm.
    pub cutoff: f64,
    /// Cutoff treatment for Lennard-Jones.
    pub modifier: Modifier,
    /// Electrostatic treatment.
    pub electrostatics: Electrostatics,
    /// Coulomb factor (see the module docs).
    pub coulomb_factor: f64,
    /// Scale for 1-4 Lennard-Jones interactions.
    pub fudge_lj: f64,
    /// Scale for 1-4 Coulomb interactions.
    pub fudge_qq: f64,
}

/// Non-bonded force field with exclusions and 1-4 pairs.
pub struct Nonbonded {
    settings: NonbondedSettings,
    lj: LennardJones,
    /// Sorted excluded partners of each atom (both directions).
    exclusions: Vec<Vec<usize>>,
    /// 1-4 pairs `(i, j)`, `i < j`.
    pairs14: Vec<(usize, usize)>,
    /// Ewald splitting parameter, fixed at construction.
    alpha: f64,
    /// PME solver for the current box, rebuilt when the box lengths change.
    pme_cache: Mutex<Option<([f64; 3], Pme)>>,
}

impl Nonbonded {
    /// Builds the force field.
    ///
    /// `exclusions` must hold sorted partner lists for every atom (see
    /// [`crate::topology::Topology::exclusion_partners`]); `pairs14` lists the
    /// 1-4 pairs.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for invalid Lennard-Jones or electrostatic
    /// parameters or invalid fudge factors.
    pub fn new(
        settings: NonbondedSettings,
        exclusions: Vec<Vec<usize>>,
        pairs14: Vec<(usize, usize)>,
    ) -> Result<Self, SimError> {
        let lj = LennardJones::new(&settings.lj_types, settings.cutoff, settings.modifier)?;
        if !(settings.fudge_lj.is_finite() && settings.fudge_qq.is_finite()) {
            return Err(SimError::Config("1-4 fudge factors must be finite".into()));
        }
        if !(settings.coulomb_factor.is_finite() && settings.coulomb_factor > 0.0) {
            return Err(SimError::Config("Coulomb factor must be > 0".into()));
        }
        let alpha = match settings.electrostatics {
            Electrostatics::None => 0.0,
            Electrostatics::Pme { tolerance, order } => {
                if !(tolerance > 0.0 && tolerance < 1.0) {
                    return Err(SimError::Config("PME tolerance must lie in (0, 1)".into()));
                }
                if !(2..=12).contains(&order) {
                    return Err(SimError::Config("PME order must lie in 2..=12".into()));
                }
                alpha_for_tolerance(settings.cutoff, tolerance)
            }
        };
        Ok(Self { settings, lj, exclusions, pairs14, alpha, pme_cache: Mutex::new(None) })
    }

    /// The Ewald splitting parameter (0 when electrostatics are disabled).
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// The settings this force field was built from.
    pub fn settings(&self) -> &NonbondedSettings {
        &self.settings
    }

    fn is_excluded(&self, i: usize, j: usize) -> bool {
        self.exclusions[i].binary_search(&j).is_ok()
    }

    /// Reciprocal-space PME energy and forces (in MD units, already multiplied
    /// by the Coulomb factor) for the given configuration.
    fn reciprocal(
        &self,
        positions: &[Vec3],
        charges: &[f64],
        sim_box: &SimBox,
        forces: &mut [Vec3],
    ) -> f64 {
        let lengths = [sim_box.lengths.x, sim_box.lengths.y, sim_box.lengths.z];
        let Electrostatics::Pme { order, .. } = self.settings.electrostatics else {
            return 0.0;
        };
        let mut cache = self.pme_cache.lock().unwrap_or_else(PoisonError::into_inner);
        let stale = cache.as_ref().is_none_or(|(l, _)| *l != lengths);
        if stale {
            let rc = self.settings.cutoff;
            let tol = match self.settings.electrostatics {
                Electrostatics::Pme { tolerance, .. } => tolerance,
                Electrostatics::None => 0.0,
            };
            let mut params = PmeParams::from_accuracy(lengths, rc, tol, order);
            params.alpha = self.alpha;
            *cache = Some((lengths, Pme::new(lengths, params)));
        }
        let Some((_, pme)) = cache.as_ref() else {
            return 0.0;
        };
        let x: Vec<f64> = positions.iter().map(|p| p.x).collect();
        let y: Vec<f64> = positions.iter().map(|p| p.y).collect();
        let z: Vec<f64> = positions.iter().map(|p| p.z).collect();
        let mut fx = vec![0.0; positions.len()];
        let mut fy = vec![0.0; positions.len()];
        let mut fz = vec![0.0; positions.len()];
        let energy = {
            let mut buf = ForceBuffers::new(&mut fx, &mut fy, &mut fz);
            pme.reciprocal(Positions::new(&x, &y, &z), charges, &mut buf)
        };
        let f = self.settings.coulomb_factor;
        for (i, out) in forces.iter_mut().enumerate() {
            *out += Vec3::new(fx[i], fy[i], fz[i]) * f;
        }
        let self_e = pme.self_energy(charges);
        let q_total: f64 = charges.iter().sum();
        let background = -std::f64::consts::PI * q_total * q_total
            / (2.0 * self.alpha * self.alpha * sim_box.volume());
        f * (energy + self_e + background)
    }
}

impl ForceField for Nonbonded {
    fn name(&self) -> &'static str {
        "nonbonded"
    }

    fn cutoff(&self) -> f64 {
        self.settings.cutoff
    }

    fn validate_types(&self, atom_types: &[u16]) -> Result<(), SimError> {
        self.lj.validate_types(atom_types)
    }

    fn validate_box(&self, sim_box: &SimBox) -> Result<(), SimError> {
        if matches!(self.settings.electrostatics, Electrostatics::Pme { .. })
            && !sim_box.is_orthorhombic()
        {
            return Err(SimError::Invalid("PME electrostatics require an orthorhombic box".into()));
        }
        Ok(())
    }

    fn compute(
        &self,
        positions: &[Vec3],
        atom_types: &[u16],
        charges: &[f64],
        sim_box: &SimBox,
        neighbours: &VerletList,
        forces: &mut [Vec3],
    ) -> ForceOutput {
        let rc2 = self.settings.cutoff * self.settings.cutoff;
        let f = self.settings.coulomb_factor;
        let pme = matches!(self.settings.electrostatics, Electrostatics::Pme { .. });
        let alpha = self.alpha;
        let mut out = ForceOutput::default();

        // Short-range pairs from the neighbour list.
        for (i, j) in neighbours.pairs() {
            if self.is_excluded(i, j) {
                continue;
            }
            let dr = sim_box.min_image(positions[i] - positions[j]);
            let r2 = dr.norm_sq();
            if r2 >= rc2 || r2 == 0.0 {
                continue;
            }
            let r = r2.sqrt();
            let (v, dv) = self.lj.pair(atom_types[i], atom_types[j], r);
            let mut e = v;
            let mut f_over_r = -dv / r;
            if pme {
                let qq = f * charges[i] * charges[j];
                let (ec, fc) = coulomb_ewald_pair(r2, qq, alpha);
                e += ec;
                f_over_r += fc;
            }
            let force = dr * f_over_r;
            forces[i] += force;
            forces[j] -= force;
            out.potential += e;
            out.virial += dr.dot(force);
        }

        // Scaled 1-4 interactions (no cutoff).
        for &(i, j) in &self.pairs14 {
            let dr = sim_box.min_image(positions[i] - positions[j]);
            let r2 = dr.norm_sq();
            if r2 == 0.0 {
                continue;
            }
            let r = r2.sqrt();
            let ti = atom_types[i] as usize;
            let tj = atom_types[j] as usize;
            let (sigma, eps) = lorentz_berthelot(
                self.settings.lj_types[ti].sigma,
                self.settings.lj_types[ti].epsilon,
                self.settings.lj_types[tj].sigma,
                self.settings.lj_types[tj].epsilon,
            );
            let (v_lj, flj) = lj_pair(r2, sigma, eps);
            let mut e = self.settings.fudge_lj * v_lj;
            let mut f_over_r = self.settings.fudge_lj * flj;
            let qq = self.settings.fudge_qq * f * charges[i] * charges[j];
            e += qq / r;
            f_over_r += qq / (r2 * r);
            let force = dr * f_over_r;
            forces[i] += force;
            forces[j] -= force;
            out.potential += e;
            out.virial += dr.dot(force);
        }

        // Reciprocal-space PME and the exclusion correction.
        if pme {
            for (i, partners) in self.exclusions.iter().enumerate() {
                for &j in partners.iter().filter(|&&j| j > i) {
                    let dr = sim_box.min_image(positions[i] - positions[j]);
                    let r2 = dr.norm_sq();
                    if r2 == 0.0 {
                        continue;
                    }
                    let qq = f * charges[i] * charges[j];
                    let (e, f_over_r) = excluded_ewald_pair(r2, qq, alpha);
                    let force = dr * f_over_r;
                    forces[i] += force;
                    forces[j] -= force;
                    out.potential += e;
                    out.virial += dr.dot(force);
                }
            }
            out.potential += self.reciprocal(positions, charges, sim_box, forces);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bonded::{BondedForce, BondedTerms, HarmonicBond};
    use crate::topology::Topology;

    fn settings(electrostatics: Electrostatics, cutoff: f64) -> NonbondedSettings {
        NonbondedSettings {
            lj_types: vec![LjType { epsilon: 1.0, sigma: 1.0 }],
            cutoff,
            modifier: Modifier::None,
            electrostatics,
            coulomb_factor: 1.0,
            fudge_lj: 0.5,
            fudge_qq: 0.8,
        }
    }

    fn no_excl(n: usize) -> Vec<Vec<usize>> {
        vec![Vec::new(); n]
    }

    fn eval(
        ff: &dyn ForceField,
        pos: &[Vec3],
        q: &[f64],
        b: &SimBox,
        skin: f64,
    ) -> (ForceOutput, Vec<Vec3>) {
        let mut list = VerletList::new(ff.cutoff().max(0.5), skin).unwrap();
        list.update(pos, b).unwrap();
        let mut f = vec![Vec3::ZERO; pos.len()];
        let types = vec![0u16; pos.len()];
        let o = ff.compute(pos, &types, q, b, &list, &mut f);
        (o, f)
    }

    #[test]
    fn excluded_pairs_do_not_contribute_lj() {
        let b = SimBox::new(6.0, 6.0, 6.0).unwrap();
        let pos = [Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.1, 1.0, 1.0)];
        let q = [0.0, 0.0];
        let plain =
            Nonbonded::new(settings(Electrostatics::None, 2.0), no_excl(2), vec![]).unwrap();
        let excl = vec![vec![1], vec![0]];
        let excluded = Nonbonded::new(settings(Electrostatics::None, 2.0), excl, vec![]).unwrap();
        assert!(eval(&plain, &pos, &q, &b, 0.1).0.potential.abs() > 1.0);
        assert_eq!(eval(&excluded, &pos, &q, &b, 0.1).0.potential, 0.0);
    }

    #[test]
    fn one_four_pair_uses_fudge_factors() {
        let b = SimBox::new(6.0, 6.0, 6.0).unwrap();
        let r = 2f64.powf(1.0 / 6.0);
        let pos = [Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.0 + r, 1.0, 1.0)];
        let q = [0.4, -0.4];
        let ff = Nonbonded::new(
            settings(Electrostatics::None, 2.0),
            vec![vec![1], vec![0]],
            vec![(0, 1)],
        )
        .unwrap();
        let (o, _) = eval(&ff, &pos, &q, &b, 0.1);
        // LJ minimum depth scaled by fudge_lj, plus scaled Coulomb −0.8·0.16/r.
        let expected = -0.5 + 0.8 * (-0.16) / r;
        assert!((o.potential - expected).abs() < 1e-12, "{} vs {expected}", o.potential);
    }

    #[test]
    fn pme_forces_match_finite_differences() {
        let b = SimBox::new(3.0, 3.0, 3.0).unwrap();
        let pos = vec![
            Vec3::new(0.3, 0.4, 0.5),
            Vec3::new(1.1, 0.9, 1.4),
            Vec3::new(1.9, 1.7, 0.6),
            Vec3::new(0.8, 1.9, 2.0),
        ];
        let q = [1.0, -1.0, 0.5, -0.5];
        let ff = Nonbonded::new(
            NonbondedSettings {
                lj_types: vec![LjType { epsilon: 0.0, sigma: 1.0 }],
                ..settings(Electrostatics::Pme { tolerance: 1e-5, order: 4 }, 1.1)
            },
            vec![vec![]; 4],
            vec![],
        )
        .unwrap();
        let (_, f) = eval(&ff, &pos, &q, &b, 0.2);
        let h = 1e-5;
        for a in 0..4 {
            let mut p = pos.clone();
            p[a].x += h;
            let ep = eval(&ff, &p, &q, &b, 0.2).0.potential;
            p[a].x -= 2.0 * h;
            let em = eval(&ff, &p, &q, &b, 0.2).0.potential;
            let numeric = -(ep - em) / (2.0 * h);
            assert!(
                (f[a].x - numeric).abs() < 1e-3 * (1.0 + numeric.abs()),
                "atom {a}: analytic {} numeric {numeric}",
                f[a].x
            );
        }
    }

    /// Total Coulomb energy (real + reciprocal + self) against a brute-force Ewald
    /// sum from `tpt-dynamics-compute::reference`, for a small neutral system.
    #[test]
    fn pme_energy_matches_direct_ewald_sum() {
        use tpt_dynamics_compute::reference::{
            direct_real_images, direct_reciprocal, direct_self_energy,
        };
        use tpt_dynamics_compute::types::{ForceBuffers, Positions};

        let l = 2.0;
        let b = SimBox::new(l, l, l).unwrap();
        let pos = vec![
            Vec3::new(0.2, 0.3, 0.4),
            Vec3::new(1.1, 0.7, 1.5),
            Vec3::new(1.6, 1.6, 0.3),
            Vec3::new(0.5, 1.8, 1.2),
        ];
        let q = [1.0, -1.0, 1.0, -1.0];
        let rc = 0.9;
        let ff = Nonbonded::new(
            NonbondedSettings {
                lj_types: vec![LjType { epsilon: 0.0, sigma: 1.0 }],
                cutoff: rc,
                modifier: Modifier::None,
                electrostatics: Electrostatics::Pme { tolerance: 1e-6, order: 6 },
                coulomb_factor: 1.0,
                fudge_lj: 0.5,
                fudge_qq: 0.8,
            },
            vec![vec![]; 4],
            vec![],
        )
        .unwrap();
        let mut list = VerletList::new(rc, 0.05).unwrap();
        list.update(&pos, &b).unwrap();
        let mut forces = vec![Vec3::ZERO; 4];
        let ours = ff.compute(&pos, &[0; 4], &q, &b, &list, &mut forces).potential;

        let alpha = ff.alpha();
        let x: Vec<f64> = pos.iter().map(|p| p.x).collect();
        let y: Vec<f64> = pos.iter().map(|p| p.y).collect();
        let z: Vec<f64> = pos.iter().map(|p| p.z).collect();
        let (mut fx, mut fy, mut fz) = (vec![0.0; 4], vec![0.0; 4], vec![0.0; 4]);
        let real = direct_real_images(
            Positions::new(&x, &y, &z),
            &q,
            [l, l, l],
            alpha,
            2,
            &mut ForceBuffers::new(&mut fx, &mut fy, &mut fz),
        );
        let recip = direct_reciprocal(
            Positions::new(&x, &y, &z),
            &q,
            [l, l, l],
            alpha,
            20,
            &mut ForceBuffers::new(&mut fx, &mut fy, &mut fz),
        );
        let reference = real + recip + direct_self_energy(&q, alpha);
        assert!(
            (ours - reference).abs() < 1e-4 * reference.abs().max(1.0),
            "PME {ours} vs direct Ewald {reference}"
        );
    }

    #[test]
    fn pme_rejects_triclinic_box() {
        let ff = Nonbonded::new(
            settings(Electrostatics::Pme { tolerance: 1e-4, order: 4 }, 1.0),
            vec![],
            vec![],
        )
        .unwrap();
        let tri = SimBox::triclinic(3.0, 3.0, 3.0, 0.5, 0.0, 0.0).unwrap();
        assert!(ff.validate_box(&tri).is_err());
        assert!(ff.validate_box(&SimBox::new(3.0, 3.0, 3.0).unwrap()).is_ok());
    }

    #[test]
    fn bonded_plus_nonbonded_topology_is_consistent() {
        // A 3-atom chain built from bonds: exclusions and 1-4 pairs come from the topology.
        let topo = Topology::from_bonds(4, &[(0, 1), (1, 2), (2, 3)]).unwrap();
        assert_eq!(topo.pairs14, vec![(0, 3)]);
        let partners = topo.exclusion_partners(4);
        let nb =
            Nonbonded::new(settings(Electrostatics::None, 1.5), partners, topo.pairs14.clone())
                .unwrap();
        let bonded = BondedForce::new(
            BondedTerms {
                bonds: topo
                    .bonds
                    .iter()
                    .map(|&(i, j)| HarmonicBond { i, j, k: 100.0, r0: 0.15 })
                    .collect(),
                ..BondedTerms::default()
            },
            4,
        )
        .unwrap();
        let b = SimBox::new(8.0, 8.0, 8.0).unwrap();
        let pos = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.15, 1.0, 1.0),
            Vec3::new(1.15, 1.15, 1.0),
            Vec3::new(1.3, 1.15, 1.0),
        ];
        let combined = crate::bonded::CombinedForceField::new(vec![Box::new(nb), Box::new(bonded)]);
        let (o, f) = eval(&combined, &pos, &[0.0; 4], &b, 0.2);
        assert!(o.potential.is_finite());
        let net: Vec3 = f.iter().copied().fold(Vec3::ZERO, |a, x| a + x);
        assert!(net.norm() < 1e-10, "net force {net:?}");
    }
}
