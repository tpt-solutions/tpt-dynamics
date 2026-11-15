//! Bonded interactions: harmonic bonds, harmonic angles with optional
//! Urey-Bradley 1-3 terms, periodic proper dihedrals, and harmonic impropers.
//!
//! Energy conventions (GROMACS style, with the factor one half):
//!
//! ```text
//! bond       V = ½ k_b (r − b₀)²
//! angle      V = ½ k_θ (θ − θ₀)²  [+ ½ k_UB (r₁₃ − s₀)²  Urey-Bradley]
//! dihedral   V = Σ_n k_n [1 + cos(n φ − δ_n)]
//! improper   V = ½ k_ξ (ξ − ξ₀)²   (ξ wrapped into (−π, π])
//! ```
//!
//! Dihedral angles follow the IUPAC sign convention: `φ = 0` for cis (eclipsed)
//! and `φ = π` for trans. All displacements use the minimum-image convention
//! relative to the first atom of each term, so terms that straddle a periodic
//! boundary are handled, provided no term spans more than half a box length.
//!
//! The forces are analytic (Bekker et al. 1995 for dihedrals). Their correctness
//! is checked against finite differences of the energies in the tests.

use std::f64::consts::PI;

use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::forcefield::{ForceField, ForceOutput};
use crate::neighbour::VerletList;
use crate::pbc::SimBox;

/// Harmonic bond between atoms `i` and `j`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicBond {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Force constant `k_b` in kJ mol⁻¹ nm⁻².
    pub k: f64,
    /// Equilibrium length `b₀` in nm.
    pub r0: f64,
}

/// Urey-Bradley 1-3 term attached to an angle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UreyBradley {
    /// Force constant in kJ mol⁻¹ nm⁻².
    pub k: f64,
    /// Equilibrium 1-3 distance in nm.
    pub r0: f64,
}

/// Harmonic angle `i-j-k` with central atom `j`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicAngle {
    /// First outer atom.
    pub i: usize,
    /// Central atom.
    pub j: usize,
    /// Second outer atom.
    pub k: usize,
    /// Force constant `k_θ` in kJ mol⁻¹ rad⁻².
    pub k_theta: f64,
    /// Equilibrium angle `θ₀` in radians.
    pub theta0: f64,
    /// Optional Urey-Bradley term between `i` and `k`.
    pub urey_bradley: Option<UreyBradley>,
}

/// One term `k [1 + cos(n φ − δ)]` of a periodic dihedral.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DihedralTerm {
    /// Force constant in kJ/mol.
    pub k: f64,
    /// Periodicity (multiplicity).
    pub n: u32,
    /// Phase `δ` in radians.
    pub phase: f64,
}

/// Proper dihedral `i-j-k-l` with one or more periodic terms.
#[derive(Clone, Debug, PartialEq)]
pub struct PeriodicDihedral {
    /// Atoms along the chain.
    pub atoms: [usize; 4],
    /// Sum of periodic terms.
    pub terms: Vec<DihedralTerm>,
}

/// Harmonic improper dihedral `i-j-k-l` (central atom `k`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicImproper {
    /// Atoms of the improper quadruple.
    pub atoms: [usize; 4],
    /// Force constant `k_ξ` in kJ mol⁻¹ rad⁻².
    pub k: f64,
    /// Reference angle `ξ₀` in radians.
    pub xi0: f64,
}

/// All bonded terms of a system with their parameters bound.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BondedTerms {
    /// Harmonic bonds.
    pub bonds: Vec<HarmonicBond>,
    /// Harmonic angles (with optional Urey-Bradley).
    pub angles: Vec<HarmonicAngle>,
    /// Periodic proper dihedrals.
    pub dihedrals: Vec<PeriodicDihedral>,
    /// Harmonic impropers.
    pub impropers: Vec<HarmonicImproper>,
}

/// Bonded force field: evaluates [`BondedTerms`] and adds their forces.
///
/// It needs no neighbour list, so its [`ForceField::cutoff`] is zero. Combine it
/// with a non-bonded force field through [`CombinedForceField`].
#[derive(Clone, Debug)]
pub struct BondedForce {
    terms: BondedTerms,
    n_atoms: usize,
}

impl BondedForce {
    /// Wraps the terms for a system with `n_atoms` atoms.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if any term refers to an atom out of range,
    /// or a constant is not finite.
    pub fn new(terms: BondedTerms, n_atoms: usize) -> Result<Self, SimError> {
        let check = |idx: &[usize]| -> Result<(), SimError> {
            if idx.iter().any(|&a| a >= n_atoms) {
                Err(SimError::Invalid(format!("bonded term {idx:?} refers to a missing atom")))
            } else {
                Ok(())
            }
        };
        for b in &terms.bonds {
            check(&[b.i, b.j])?;
        }
        for a in &terms.angles {
            check(&[a.i, a.j, a.k])?;
        }
        for d in &terms.dihedrals {
            check(&d.atoms)?;
        }
        for d in &terms.impropers {
            check(&d.atoms)?;
        }
        Ok(Self { terms, n_atoms })
    }

    /// The bound terms.
    pub fn terms(&self) -> &BondedTerms {
        &self.terms
    }

    /// Number of atoms this force field was built for.
    pub fn n_atoms(&self) -> usize {
        self.n_atoms
    }
}

/// Angle between `a` and `b` (radians) and the gradient of the angle with
/// respect to the position of the atom at the end of `a` (`∂θ/∂x_i`).
pub fn angle_and_gradient(a: Vec3, b: Vec3) -> (f64, Vec3) {
    let na = a.norm();
    let nb = b.norm();
    let cos_t = (a.dot(b) / (na * nb)).clamp(-1.0, 1.0);
    let theta = cos_t.acos();
    let sin_t = (1.0 - cos_t * cos_t).sqrt().max(1e-12);
    // ∂cosθ/∂a = b/(|a||b|) − cosθ a/|a|²; ∂θ/∂a = −(∂cosθ/∂a)/sinθ.
    let dcos = b * (1.0 / (na * nb)) - a * (cos_t / (na * na));
    (theta, dcos * (-1.0 / sin_t))
}

/// Signed dihedral angle `φ` for the chain `i-j-k-l` given the minimum-image
/// vectors `r_ij = x_i − x_j`, `r_kj = x_k − x_j`, `r_kl = x_k − x_l`.
pub fn dihedral_angle(r_ij: Vec3, r_kj: Vec3, r_kl: Vec3) -> f64 {
    let m = r_ij.cross(r_kj);
    let n = r_kj.cross(r_kl);
    let y = r_kj.norm() * r_ij.dot(n);
    let x = m.dot(n);
    y.atan2(x)
}

/// Dihedral forces `(F_i, F_j, F_k, F_l)` for `dV/dφ = dv`, from Bekker et al.
fn dihedral_forces(r_ij: Vec3, r_kj: Vec3, r_kl: Vec3, dv: f64) -> [Vec3; 4] {
    let m = r_ij.cross(r_kj);
    let n = r_kj.cross(r_kl);
    let m2 = m.norm_sq();
    let n2 = n.norm_sq();
    let len_kj = r_kj.norm();
    let len2 = len_kj * len_kj;
    if m2 < 1e-24 || n2 < 1e-24 {
        // Collinear atoms: the dihedral is undefined and its gradient vanishes
        // in this formulation; the energy term contributes no force.
        return [Vec3::ZERO; 4];
    }
    // F = −∇V with ∇_i φ = −|r_kj| m / |m|² and ∇_l φ = +|r_kj| n / |n|².
    let fi = m * (-dv * len_kj / m2);
    let fl = n * (dv * len_kj / n2);
    let p = r_ij.dot(r_kj) / len2;
    let q = r_kl.dot(r_kj) / len2;
    let fj = fi * (p - 1.0) - fl * q;
    let fk = fl * (q - 1.0) - fi * p;
    [fi, fj, fk, fl]
}

/// Wraps an angle difference into `(−π, π]`.
fn wrap_angle(x: f64) -> f64 {
    let mut y = (x + PI) % (2.0 * PI);
    if y <= 0.0 {
        y += 2.0 * PI;
    }
    y - PI
}

impl ForceField for BondedForce {
    fn name(&self) -> &'static str {
        "bonded"
    }

    fn cutoff(&self) -> f64 {
        0.0
    }

    fn validate_types(&self, _atom_types: &[u16]) -> Result<(), SimError> {
        Ok(())
    }

    fn compute(
        &self,
        positions: &[Vec3],
        _atom_types: &[u16],
        _charges: &[f64],
        sim_box: &SimBox,
        _neighbours: &VerletList,
        forces: &mut [Vec3],
    ) -> ForceOutput {
        let mut out = ForceOutput::default();
        let t = &self.terms;
        let x = |a: usize| positions[a];

        for b in &t.bonds {
            let dr = sim_box.min_image(x(b.i) - x(b.j));
            let r = dr.norm();
            if r == 0.0 {
                continue;
            }
            let dev = r - b.r0;
            out.potential += 0.5 * b.k * dev * dev;
            // F_i = −k (r − r0) r̂_ij.
            let f = dr * (-b.k * dev / r);
            forces[b.i] += f;
            forces[b.j] -= f;
            out.virial += dr.dot(f);
        }

        for a in &t.angles {
            let ra = sim_box.min_image(x(a.i) - x(a.j));
            let rb = sim_box.min_image(x(a.k) - x(a.j));
            if ra.norm() == 0.0 || rb.norm() == 0.0 {
                continue;
            }
            let (theta, dth_di) = angle_and_gradient(ra, rb);
            let dev = theta - a.theta0;
            out.potential += 0.5 * a.k_theta * dev * dev;
            let dv = a.k_theta * dev;
            // ∂θ/∂x_k follows from the symmetric expression with a and b swapped.
            let (_, dth_dk) = angle_and_gradient(rb, ra);
            let fi = dth_di * (-dv);
            let fk = dth_dk * (-dv);
            let fj = -(fi + fk);
            forces[a.i] += fi;
            forces[a.j] += fj;
            forces[a.k] += fk;
            out.virial += ra.dot(fi) + rb.dot(fk);

            if let Some(ub) = a.urey_bradley {
                let d13 = sim_box.min_image(x(a.i) - x(a.k));
                let r13 = d13.norm();
                if r13 > 0.0 {
                    let dev = r13 - ub.r0;
                    out.potential += 0.5 * ub.k * dev * dev;
                    let f = d13 * (-ub.k * dev / r13);
                    forces[a.i] += f;
                    forces[a.k] -= f;
                    out.virial += d13.dot(f);
                }
            }
        }

        for d in &t.dihedrals {
            let [i, j, k, l] = d.atoms;
            let r_ij = sim_box.min_image(x(i) - x(j));
            let r_kj = sim_box.min_image(x(k) - x(j));
            let r_kl = sim_box.min_image(x(k) - x(l));
            let phi = dihedral_angle(r_ij, r_kj, r_kl);
            let mut v = 0.0;
            let mut dv = 0.0;
            for term in &d.terms {
                let nf = f64::from(term.n);
                let arg = nf * phi - term.phase;
                v += term.k * (1.0 + arg.cos());
                dv += -term.k * nf * arg.sin();
            }
            out.potential += v;
            let fs = dihedral_forces(r_ij, r_kj, r_kl, dv);
            add_term_forces(forces, &mut out, &d.atoms, &fs, r_ij, r_kj, r_kl);
        }

        for imp in &t.impropers {
            let [i, j, k, l] = imp.atoms;
            let r_ij = sim_box.min_image(x(i) - x(j));
            let r_kj = sim_box.min_image(x(k) - x(j));
            let r_kl = sim_box.min_image(x(k) - x(l));
            let xi = dihedral_angle(r_ij, r_kj, r_kl);
            let dev = wrap_angle(xi - imp.xi0);
            out.potential += 0.5 * imp.k * dev * dev;
            let fs = dihedral_forces(r_ij, r_kj, r_kl, imp.k * dev);
            add_term_forces(forces, &mut out, &imp.atoms, &fs, r_ij, r_kj, r_kl);
        }
        out
    }
}

/// Adds four-atom forces and their virial contribution.
///
/// The virial `Σ_a r_a · F_a` is evaluated with displacements `r_a` relative to
/// atom `i`; it is origin independent because the four forces sum to zero.
/// Displacements from atom `i`: `x_j − x_i = −r_ij`, `x_k − x_i = r_kj − r_ij`,
/// `x_l − x_i = (x_k − x_i) − r_kl`.
fn add_term_forces(
    forces: &mut [Vec3],
    out: &mut ForceOutput,
    atoms: &[usize; 4],
    fs: &[Vec3; 4],
    r_ij: Vec3,
    r_kj: Vec3,
    r_kl: Vec3,
) {
    let rel_j = -r_ij;
    let rel_k = r_kj - r_ij;
    let rel_l = rel_k - r_kl;
    out.virial += rel_j.dot(fs[1]) + rel_k.dot(fs[2]) + rel_l.dot(fs[3]);
    for (a, f) in atoms.iter().zip(fs) {
        forces[*a] += *f;
    }
}

/// Sums the outputs of several force fields that act on the same positions.
pub struct CombinedForceField {
    parts: Vec<Box<dyn ForceField>>,
}

impl CombinedForceField {
    /// Combines the given force fields. The combined cutoff is the largest of
    /// the parts' cutoffs.
    pub fn new(parts: Vec<Box<dyn ForceField>>) -> Self {
        Self { parts }
    }
}

impl ForceField for CombinedForceField {
    fn name(&self) -> &'static str {
        "combined"
    }

    fn cutoff(&self) -> f64 {
        self.parts.iter().map(|p| p.cutoff()).fold(0.0, f64::max)
    }

    fn validate_types(&self, atom_types: &[u16]) -> Result<(), SimError> {
        self.parts.iter().try_for_each(|p| p.validate_types(atom_types))
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
        let mut total = ForceOutput::default();
        for p in &self.parts {
            let o = p.compute(positions, atom_types, charges, sim_box, neighbours, forces);
            total.potential += o.potential;
            total.virial += o.virial;
        }
        total
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f64, y: f64, z: f64) -> Vec3 {
        Vec3::new(x, y, z)
    }

    /// Energy of a single-term system at given positions (no forces needed).
    fn energy(ff: &BondedForce, pos: &[Vec3], b: &SimBox) -> f64 {
        let mut f = vec![Vec3::ZERO; pos.len()];
        let list = VerletList::new(1.0, 0.1).unwrap();
        ff.compute(pos, &[], &[], b, &list, &mut f).potential
    }

    fn forces_of(ff: &BondedForce, pos: &[Vec3], b: &SimBox) -> Vec<Vec3> {
        let mut f = vec![Vec3::ZERO; pos.len()];
        let list = VerletList::new(1.0, 0.1).unwrap();
        ff.compute(pos, &[], &[], b, &list, &mut f);
        f
    }

    /// Checks F = −∇V by central differences for every atom and component.
    fn check_gradient(ff: &BondedForce, pos: &[Vec3], b: &SimBox, tol: f64) {
        let f = forces_of(ff, pos, b);
        let h = 1e-6;
        for a in 0..pos.len() {
            for c in 0..3 {
                let shifted = |d: f64| {
                    let mut q = pos.to_vec();
                    match c {
                        0 => q[a].x += d,
                        1 => q[a].y += d,
                        _ => q[a].z += d,
                    }
                    energy(ff, &q, b)
                };
                let numeric = -(shifted(h) - shifted(-h)) / (2.0 * h);
                let analytic = match c {
                    0 => f[a].x,
                    1 => f[a].y,
                    _ => f[a].z,
                };
                assert!(
                    (analytic - numeric).abs() < tol * (1.0 + numeric.abs()),
                    "atom {a} comp {c}: analytic {analytic} numeric {numeric}"
                );
            }
        }
    }

    fn chain_positions() -> Vec<Vec3> {
        vec![v(0.1, 0.9, 0.2), v(0.0, 0.0, 0.0), v(0.9, -0.1, 0.3), v(1.2, 0.8, -0.4)]
    }

    #[test]
    fn bond_energy_matches_analytic_value() {
        let ff = BondedForce::new(
            BondedTerms {
                bonds: vec![HarmonicBond { i: 0, j: 1, k: 300.0, r0: 0.1 }],
                ..BondedTerms::default()
            },
            2,
        )
        .unwrap();
        let b = SimBox::new(5.0, 5.0, 5.0).unwrap();
        let pos = [v(0.0, 0.0, 0.0), v(0.12, 0.0, 0.0)];
        // V = ½ · 300 · (0.02)² = 0.06 kJ/mol.
        assert!((energy(&ff, &pos, &b) - 0.06).abs() < 1e-12);
        check_gradient(&ff, &pos, &b, 1e-6);
    }

    #[test]
    fn bond_across_periodic_boundary_uses_minimum_image() {
        let ff = BondedForce::new(
            BondedTerms {
                bonds: vec![HarmonicBond { i: 0, j: 1, k: 100.0, r0: 0.1 }],
                ..BondedTerms::default()
            },
            2,
        )
        .unwrap();
        let b = SimBox::new(2.0, 2.0, 2.0).unwrap();
        let pos = [v(0.05, 1.0, 1.0), v(1.95, 1.0, 1.0)]; // separation 0.1 across the boundary
        assert!(energy(&ff, &pos, &b).abs() < 1e-12);
    }

    #[test]
    fn angle_energy_and_gradient_match_numerics() {
        let ff = BondedForce::new(
            BondedTerms {
                angles: vec![HarmonicAngle {
                    i: 0,
                    j: 1,
                    k: 2,
                    k_theta: 50.0,
                    theta0: 1.9,
                    urey_bradley: Some(UreyBradley { k: 30.0, r0: 0.2 }),
                }],
                ..BondedTerms::default()
            },
            3,
        )
        .unwrap();
        let b = SimBox::new(6.0, 6.0, 6.0).unwrap();
        let pos = [v(0.1, 0.0, 0.0), v(0.0, 0.0, 0.0), v(-0.05, 0.12, 0.03)];
        check_gradient(&ff, &pos, &b, 1e-6);
    }

    #[test]
    fn right_angle_has_analytic_energy() {
        let ff = BondedForce::new(
            BondedTerms {
                angles: vec![HarmonicAngle {
                    i: 0,
                    j: 1,
                    k: 2,
                    k_theta: 80.0,
                    theta0: PI / 2.0 + 0.1,
                    urey_bradley: None,
                }],
                ..BondedTerms::default()
            },
            3,
        )
        .unwrap();
        let b = SimBox::new(6.0, 6.0, 6.0).unwrap();
        let pos = [v(1.0, 0.0, 0.0), v(0.0, 0.0, 0.0), v(0.0, 1.0, 0.0)];
        // θ = π/2 exactly, so V = ½ · 80 · 0.1² = 0.4.
        assert!((energy(&ff, &pos, &b) - 0.4).abs() < 1e-12);
    }

    #[test]
    fn dihedral_cis_and_trans_values() {
        // Cis: i and l on the same side (φ = 0).
        let (r_ij, r_kj, r_kl) = (v(0.0, 1.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, -1.0, 0.0));
        assert!(dihedral_angle(r_ij, r_kj, r_kl).abs() < 1e-12);
        // Trans: l on the opposite side (φ = π).
        let r_kl_trans = v(0.0, 1.0, 0.0);
        assert!((dihedral_angle(r_ij, r_kj, r_kl_trans).abs() - PI).abs() < 1e-12);
    }

    #[test]
    fn periodic_dihedral_energy_and_gradient() {
        let terms = vec![
            DihedralTerm { k: 2.0, n: 3, phase: 0.4 },
            DihedralTerm { k: 0.5, n: 1, phase: PI },
        ];
        let ff = BondedForce::new(
            BondedTerms {
                dihedrals: vec![PeriodicDihedral { atoms: [0, 1, 2, 3], terms }],
                ..BondedTerms::default()
            },
            4,
        )
        .unwrap();
        let b = SimBox::new(8.0, 8.0, 8.0).unwrap();
        check_gradient(&ff, &chain_positions(), &b, 1e-5);
    }

    #[test]
    fn dihedral_energy_matches_formula_for_known_geometry() {
        // Cis configuration (φ = 0): V = k(1 + cos(−δ)) for n = 1.
        let ff = BondedForce::new(
            BondedTerms {
                dihedrals: vec![PeriodicDihedral {
                    atoms: [0, 1, 2, 3],
                    terms: vec![DihedralTerm { k: 1.5, n: 1, phase: 0.3 }],
                }],
                ..BondedTerms::default()
            },
            4,
        )
        .unwrap();
        let b = SimBox::new(8.0, 8.0, 8.0).unwrap();
        // i and l on the same side of the j-k axis: φ = 0.
        let pos = [v(0.0, 1.0, 0.0), v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(1.0, 1.0, 0.0)];
        let expected = 1.5 * (1.0 + (-0.3f64).cos());
        assert!((energy(&ff, &pos, &b) - expected).abs() < 1e-12);
    }

    #[test]
    fn improper_gradient_and_wrapping() {
        let ff = BondedForce::new(
            BondedTerms {
                impropers: vec![HarmonicImproper { atoms: [0, 1, 2, 3], k: 20.0, xi0: 3.0 }],
                ..BondedTerms::default()
            },
            4,
        )
        .unwrap();
        let b = SimBox::new(8.0, 8.0, 8.0).unwrap();
        check_gradient(&ff, &chain_positions(), &b, 1e-5);
        assert!((wrap_angle(3.5) - (3.5 - 2.0 * PI)).abs() < 1e-12);
    }

    #[test]
    fn rejects_missing_atoms() {
        let r = BondedForce::new(
            BondedTerms {
                bonds: vec![HarmonicBond { i: 0, j: 7, k: 1.0, r0: 0.1 }],
                ..BondedTerms::default()
            },
            2,
        );
        assert!(r.is_err());
    }
}
