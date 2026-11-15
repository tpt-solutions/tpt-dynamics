//! Force-field trait and the Lennard-Jones pair potential.
//!
//! The [`ForceField`] trait is the plug-in point for potentials: anything that
//! can add forces for the pairs in a [`VerletList`] and report its energy and
//! virial can drive the integrators. [`LennardJones`] is the first implementation.
//!
//! Energies are in kJ/mol when lengths are in nm (or in reduced units when
//! σ = ε = 1). The pair potential is
//! `V(r) = 4ε[(σ/r)¹² − (σ/r)⁶]` for `r < r_c`, with per-pair constants
//! combined by Lorentz-Berthelot rules `σ_ij = (σ_i + σ_j)/2`,
//! `ε_ij = √(ε_i ε_j)`.

use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::neighbour::VerletList;
use crate::pbc::SimBox;

/// Scalar results of a force evaluation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ForceOutput {
    /// Total potential energy (kJ/mol or reduced).
    pub potential: f64,
    /// Pair virial `Σ_pairs r_ij · f_ij`, where `f_ij` is the force on `i` due to `j`.
    /// Used for the pressure: `P = (2 KE + W) / (3 V)`.
    pub virial: f64,
}

/// A pairwise (or more general) potential that can be evaluated for a system.
pub trait ForceField: Send + Sync {
    /// Short identifier used in logs.
    fn name(&self) -> &'static str;

    /// Interaction cutoff in nm. Used to size the neighbour list.
    fn cutoff(&self) -> f64;

    /// Checks that the force field can be applied to a system with the given
    /// atom types. The default accepts everything.
    fn validate_types(&self, _atom_types: &[u16]) -> Result<(), SimError> {
        Ok(())
    }

    /// Checks that the force field supports the periodic box (for example, PME
    /// requires an orthorhombic box). The default accepts every box.
    fn validate_box(&self, _sim_box: &SimBox) -> Result<(), SimError> {
        Ok(())
    }

    /// Adds forces into `forces` (which the caller has zeroed) and returns the
    /// potential energy and virial. `neighbours` must cover all pairs closer
    /// than [`ForceField::cutoff`].
    fn compute(
        &self,
        positions: &[Vec3],
        atom_types: &[u16],
        charges: &[f64],
        sim_box: &SimBox,
        neighbours: &VerletList,
        forces: &mut [Vec3],
    ) -> ForceOutput;
}

/// Lennard-Jones parameters for one atom type.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LjType {
    /// Well depth ε (kJ/mol or reduced).
    pub epsilon: f64,
    /// Zero-crossing distance σ (nm or reduced).
    pub sigma: f64,
}

/// How the Lennard-Jones potential is modified at the cutoff.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Modifier {
    /// Plain truncation: energy is discontinuous at the cutoff.
    None,
    /// Energy shifted by `V(r_c)`: energy is continuous, force is not.
    Shift,
    /// Smooth polynomial switch between `r_switch` and `r_c`: energy and force
    /// both go to zero continuously at `r_c`.
    Switch {
        /// Distance (nm) where the switching function starts.
        r_switch: f64,
    },
}

#[derive(Clone, Copy, Debug)]
struct PairParams {
    c6: f64,
    c12: f64,
    /// `V(r_c)` of the unmodified potential, subtracted for `Modifier::Shift`.
    shift: f64,
}

/// Lennard-Jones 12-6 potential with a cutoff and an optional modifier.
#[derive(Clone, Debug)]
pub struct LennardJones {
    n_types: usize,
    table: Vec<PairParams>,
    cutoff: f64,
    modifier: Modifier,
}

impl LennardJones {
    /// Builds the force field from per-type parameters.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for empty tables, non-positive σ, negative ε,
    /// a non-positive cutoff, or a switch start outside `(0, r_c)`.
    pub fn new(types: &[LjType], cutoff: f64, modifier: Modifier) -> Result<Self, SimError> {
        if types.is_empty() {
            return Err(SimError::Config(
                "at least one Lennard-Jones atom type is required".into(),
            ));
        }
        if !(cutoff.is_finite() && cutoff > 0.0) {
            return Err(SimError::Config(format!("LJ cutoff must be > 0, got {cutoff}")));
        }
        for (k, t) in types.iter().enumerate() {
            if !(t.sigma.is_finite() && t.sigma > 0.0) {
                return Err(SimError::Config(format!("LJ type {k}: sigma must be > 0")));
            }
            if !(t.epsilon.is_finite() && t.epsilon >= 0.0) {
                return Err(SimError::Config(format!("LJ type {k}: epsilon must be >= 0")));
            }
        }
        if let Modifier::Switch { r_switch } = modifier {
            if !(r_switch > 0.0 && r_switch < cutoff) {
                return Err(SimError::Config("switch start must lie in (0, cutoff)".into()));
            }
        }
        let n = types.len();
        let mut table = Vec::with_capacity(n * n);
        for a in types {
            for b in types {
                let sigma = f64::midpoint(a.sigma, b.sigma);
                let epsilon = (a.epsilon * b.epsilon).sqrt();
                let s6 = sigma.powi(6);
                let c6 = 4.0 * epsilon * s6;
                let c12 = 4.0 * epsilon * s6 * s6;
                let shift = unmodified(c6, c12, cutoff).0;
                table.push(PairParams { c6, c12, shift });
            }
        }
        Ok(Self { n_types: n, table, cutoff, modifier })
    }

    /// Convenience constructor for a single atom type.
    pub fn single_type(
        epsilon: f64,
        sigma: f64,
        cutoff: f64,
        modifier: Modifier,
    ) -> Result<Self, SimError> {
        Self::new(&[LjType { epsilon, sigma }], cutoff, modifier)
    }

    /// Number of atom types in the parameter table.
    pub fn n_types(&self) -> usize {
        self.n_types
    }

    /// Pair energy and `dV/dr` for atoms of types `ti` and `tj` at distance `r`
    /// (zero beyond the cutoff). Exposed for tests and analysis.
    pub fn pair(&self, ti: u16, tj: u16, r: f64) -> (f64, f64) {
        let p = self.table[ti as usize * self.n_types + tj as usize];
        if r >= self.cutoff {
            return (0.0, 0.0);
        }
        let (v, dv) = unmodified(p.c6, p.c12, r);
        match self.modifier {
            Modifier::None => (v, dv),
            Modifier::Shift => (v - p.shift, dv),
            Modifier::Switch { r_switch } => {
                if r <= r_switch {
                    (v, dv)
                } else {
                    let width = self.cutoff - r_switch;
                    let t = (r - r_switch) / width;
                    let t2 = t * t;
                    let t3 = t2 * t;
                    let s = 1.0 - 10.0 * t3 + 15.0 * t2 * t2 - 6.0 * t3 * t2;
                    let ds_dr = (-30.0 * t2 + 60.0 * t3 - 30.0 * t2 * t2) / width;
                    (v * s, dv * s + v * ds_dr)
                }
            }
        }
    }
}

/// Unmodified LJ energy and derivative for constants `c6`, `c12` at distance `r`.
fn unmodified(c6: f64, c12: f64, r: f64) -> (f64, f64) {
    let inv_r2 = 1.0 / (r * r);
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;
    let v = c12 * inv_r6 * inv_r6 - c6 * inv_r6;
    let dv = (-12.0 * c12 * inv_r6 * inv_r6 + 6.0 * c6 * inv_r6) / r;
    (v, dv)
}

impl ForceField for LennardJones {
    fn name(&self) -> &'static str {
        "lennard-jones"
    }

    fn cutoff(&self) -> f64 {
        self.cutoff
    }

    fn validate_types(&self, atom_types: &[u16]) -> Result<(), SimError> {
        match atom_types.iter().find(|&&t| t as usize >= self.n_types) {
            Some(t) => Err(SimError::Invalid(format!(
                "atom type {t} has no Lennard-Jones parameters ({} types defined)",
                self.n_types
            ))),
            None => Ok(()),
        }
    }

    fn compute(
        &self,
        positions: &[Vec3],
        atom_types: &[u16],
        _charges: &[f64],
        sim_box: &SimBox,
        neighbours: &VerletList,
        forces: &mut [Vec3],
    ) -> ForceOutput {
        let rc2 = self.cutoff * self.cutoff;
        let mut out = ForceOutput::default();
        for (i, j) in neighbours.pairs() {
            let dr = sim_box.min_image(positions[i] - positions[j]);
            let r2 = dr.norm_sq();
            if r2 >= rc2 || r2 == 0.0 {
                continue;
            }
            let r = r2.sqrt();
            let (v, dv) = self.pair(atom_types[i], atom_types[j], r);
            // Force on i from j is -dV/dr * r̂_ij.
            let f = dr * (-dv / r);
            forces[i] += f;
            forces[j] -= f;
            out.potential += v;
            out.virial += dr.dot(f);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lj(modifier: Modifier) -> LennardJones {
        LennardJones::single_type(1.0, 1.0, 2.5, modifier).unwrap()
    }

    #[test]
    fn minimum_is_at_two_to_one_sixth_with_depth_epsilon() {
        let p = lj(Modifier::None);
        let rmin = 2f64.powf(1.0 / 6.0);
        let (v, dv) = p.pair(0, 0, rmin);
        assert!((v + 1.0).abs() < 1e-12);
        assert!(dv.abs() < 1e-10);
    }

    #[test]
    fn force_equals_negative_numerical_derivative() {
        for m in [Modifier::None, Modifier::Shift, Modifier::Switch { r_switch: 2.0 }] {
            let p = lj(m);
            for k in 0..200_usize {
                let r = 0.8 + 0.0085 * k as f64;
                let h = 1e-6;
                let (vp, _) = p.pair(0, 0, r + h);
                let (vm, _) = p.pair(0, 0, r - h);
                let numeric = (vp - vm) / (2.0 * h);
                let (_, analytic) = p.pair(0, 0, r);
                assert!(
                    (analytic - numeric).abs() < 1e-5 * (1.0 + numeric.abs()),
                    "{m:?} r={r}: analytic {analytic} numeric {numeric}"
                );
            }
        }
    }

    #[test]
    fn shift_makes_energy_vanish_at_cutoff() {
        let p = lj(Modifier::Shift);
        let (v, _) = p.pair(0, 0, 2.5 - 1e-9);
        assert!(v.abs() < 1e-7);
        assert_eq!(p.pair(0, 0, 2.5), (0.0, 0.0));
    }

    #[test]
    fn switch_is_continuous_in_energy_and_force_at_cutoff() {
        let p = lj(Modifier::Switch { r_switch: 2.0 });
        let (v, dv) = p.pair(0, 0, 2.5 - 1e-9);
        assert!(v.abs() < 1e-7 && dv.abs() < 1e-5, "v {v} dv {dv}");
        let (vs, _) = p.pair(0, 0, 2.0);
        let (vu, _) = unmodified_pair(2.0);
        assert!((vs - vu).abs() < 1e-12, "switch must be identity below r_switch");
    }

    fn unmodified_pair(r: f64) -> (f64, f64) {
        lj(Modifier::None).pair(0, 0, r)
    }

    #[test]
    fn lorentz_berthelot_combination() {
        let p = LennardJones::new(
            &[LjType { epsilon: 1.0, sigma: 1.0 }, LjType { epsilon: 4.0, sigma: 3.0 }],
            5.0,
            Modifier::None,
        )
        .unwrap();
        // sigma_01 = 2, epsilon_01 = 2: minimum depth -2 at r = 2^(1/6) * 2.
        let (v, _) = p.pair(0, 1, 2f64.powf(1.0 / 6.0) * 2.0);
        assert!((v + 2.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(LennardJones::single_type(1.0, 0.0, 2.5, Modifier::None).is_err());
        assert!(LennardJones::single_type(-1.0, 1.0, 2.5, Modifier::None).is_err());
        assert!(LennardJones::single_type(1.0, 1.0, 0.0, Modifier::None).is_err());
        assert!(
            LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Switch { r_switch: 3.0 }).is_err()
        );
        assert!(LennardJones::new(&[], 2.5, Modifier::None).is_err());
    }
}
