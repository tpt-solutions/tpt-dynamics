//! Holonomic constraints applied after each position update.
//!
//! Three algorithms are provided behind the [`Constraints`] trait:
//!
//! * [`Shake`]: iterative bond constraints (Ryckaert, Ciccotti & Berendsen 1977).
//!   Each constraint is relaxed in turn until every length is within the
//!   requested relative tolerance.
//! * [`Lincs`]: a LINCS-style solver (Hess et al. 1997). The coupled linear
//!   system for the constraint multipliers is solved with a truncated matrix
//!   expansion (`order` terms), followed by `iterations` centripetal-correction
//!   passes. This is the sequential form: GROMACS's parallel matrix-free variant
//!   is not implemented. The method is accurate for weakly coupled constraints
//!   (for example, isolated bonds and bonds on a chain with small coupling); the
//!   achieved accuracy is reported by the tests and it degrades for tightly
//!   coupled rings.
//! * [`Settle`]: rigid three-site water (O, H, H). The constrained geometry is
//!   the exact solution of the three coupled distance constraints. It is found
//!   by Newton iteration on the three multipliers rather than the closed-form
//!   Miyamoto-Kollman expressions; the two give the same configuration, and the
//!   tests compare the result with an independently converged SHAKE solution.
//!
//! Position constraints are applied through [`apply_constraints`], which also
//! corrects velocities by `Δx/dt` (the usual leap-frog/velocity-Verlet
//! correction). [`project_velocities`] removes velocity components along the
//! constraints, which is the velocity stage of RATTLE and is used for the
//! initial velocities.
//!
//! Constraint virials are not included in the pressure. Pressures reported for
//! constrained systems therefore differ from a full GROMACS-style virial.

use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::pbc::SimBox;
use crate::system::System;

/// A set of holonomic distance constraints.
pub trait Constraints: Send {
    /// Short identifier used in logs.
    fn name(&self) -> &'static str;

    /// False only for the no-op implementation, which lets integrators skip the
    /// reference-position copy.
    fn is_active(&self) -> bool {
        true
    }

    /// The constrained atom pairs `(i, j, d)` with target length `d` (nm).
    fn pairs(&self) -> Vec<(usize, usize, f64)>;

    /// Moves `positions` onto the constraint surface. `reference` holds the
    /// positions before the unconstrained update; the directions of the
    /// constraints are taken from it.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if the iterative solver does not converge.
    fn constrain_positions(
        &mut self,
        positions: &mut [Vec3],
        reference: &[Vec3],
        masses: &[f64],
        sim_box: &SimBox,
    ) -> Result<(), SimError>;
}

/// The no-op constraint set used when a system has no constraints.
#[derive(Clone, Copy, Debug, Default)]
pub struct NoConstraints;

impl Constraints for NoConstraints {
    fn name(&self) -> &'static str {
        "none"
    }

    fn is_active(&self) -> bool {
        false
    }

    fn pairs(&self) -> Vec<(usize, usize, f64)> {
        Vec::new()
    }

    fn constrain_positions(
        &mut self,
        _positions: &mut [Vec3],
        _reference: &[Vec3],
        _masses: &[f64],
        _sim_box: &SimBox,
    ) -> Result<(), SimError> {
        Ok(())
    }
}

/// Several constraint algorithms applied in sequence (for example, SHAKE for
/// solute bonds and SETTLE for water). The atom sets must not overlap.
#[derive(Default)]
pub struct ConstraintStack {
    parts: Vec<Box<dyn Constraints>>,
}

impl ConstraintStack {
    /// Creates an empty stack.
    pub fn new() -> Self {
        Self::default()
    }

    /// Appends a constraint set.
    pub fn push(&mut self, c: Box<dyn Constraints>) {
        self.parts.push(c);
    }
}

impl Constraints for ConstraintStack {
    fn name(&self) -> &'static str {
        "stack"
    }

    fn is_active(&self) -> bool {
        self.parts.iter().any(|p| p.is_active())
    }

    fn pairs(&self) -> Vec<(usize, usize, f64)> {
        self.parts.iter().flat_map(|p| p.pairs()).collect()
    }

    fn constrain_positions(
        &mut self,
        positions: &mut [Vec3],
        reference: &[Vec3],
        masses: &[f64],
        sim_box: &SimBox,
    ) -> Result<(), SimError> {
        for p in &mut self.parts {
            p.constrain_positions(positions, reference, masses, sim_box)?;
        }
        Ok(())
    }
}

/// Applies `constraints` to the system after an unconstrained update. Positions
/// are constrained, velocities are corrected by the displacement divided by
/// `dt`, and positions are wrapped back into the box.
///
/// # Errors
/// Propagates solver errors.
pub fn apply_constraints(
    constraints: &mut dyn Constraints,
    system: &mut System,
    reference: &[Vec3],
    dt: f64,
) -> Result<(), SimError> {
    if !constraints.is_active() {
        return Ok(());
    }
    let unconstrained = system.positions.clone();
    let sim_box = system.sim_box;
    constraints.constrain_positions(&mut system.positions, reference, &system.masses, &sim_box)?;
    for i in 0..system.len() {
        let shift = sim_box.min_image(system.positions[i] - unconstrained[i]);
        system.velocities[i] += shift * (1.0 / dt);
        system.positions[i] = sim_box.wrap(system.positions[i]);
    }
    Ok(())
}

/// Removes velocity components along the given constrained pairs (the RATTLE
/// velocity stage), iterating until the relative radial velocities are below
/// `tol` or `max_iter` sweeps have run.
///
/// # Errors
/// Returns [`SimError::Invalid`] if the sweeps do not converge.
pub fn project_velocities(
    pairs: &[(usize, usize, f64)],
    system: &mut System,
    tol: f64,
    max_iter: usize,
) -> Result<(), SimError> {
    let sim_box = system.sim_box;
    for _ in 0..max_iter {
        let mut worst: f64 = 0.0;
        for &(i, j, d) in pairs {
            let r = sim_box.min_image(system.positions[i] - system.positions[j]);
            let dv = system.velocities[i] - system.velocities[j];
            let rdv = r.dot(dv);
            worst = worst.max(rdv.abs() / d);
            let mu = 1.0 / system.masses[i] + 1.0 / system.masses[j];
            let lambda = -rdv / (r.norm_sq() * mu);
            system.velocities[i] += r * (lambda / system.masses[i]);
            system.velocities[j] -= r * (lambda / system.masses[j]);
        }
        if worst <= tol {
            return Ok(());
        }
    }
    Err(SimError::Invalid("velocity constraint projection did not converge".into()))
}

fn validate_pairs(pairs: &[(usize, usize, f64)], n_atoms: usize) -> Result<(), SimError> {
    for &(i, j, d) in pairs {
        if i == j || i >= n_atoms || j >= n_atoms {
            return Err(SimError::Invalid(format!("constraint ({i}, {j}) is invalid")));
        }
        if !(d.is_finite() && d > 0.0) {
            return Err(SimError::Invalid(format!("constraint length {d} must be > 0")));
        }
    }
    Ok(())
}

/// Iterative SHAKE for bond-length constraints.
#[derive(Clone, Debug)]
pub struct Shake {
    pairs: Vec<(usize, usize, f64)>,
    tolerance: f64,
    max_iter: usize,
}

impl Shake {
    /// Creates a SHAKE solver. `tolerance` is the relative length tolerance
    /// (`| |r| − d | / d`), `max_iter` the number of sweeps allowed.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] for invalid pairs or parameters.
    pub fn new(
        pairs: Vec<(usize, usize, f64)>,
        n_atoms: usize,
        tolerance: f64,
        max_iter: usize,
    ) -> Result<Self, SimError> {
        validate_pairs(&pairs, n_atoms)?;
        if !(tolerance.is_finite() && tolerance > 0.0) || max_iter == 0 {
            return Err(SimError::Invalid("SHAKE tolerance must be > 0 and max_iter >= 1".into()));
        }
        Ok(Self { pairs, tolerance, max_iter })
    }
}

impl Constraints for Shake {
    fn name(&self) -> &'static str {
        "shake"
    }

    fn pairs(&self) -> Vec<(usize, usize, f64)> {
        self.pairs.clone()
    }

    fn constrain_positions(
        &mut self,
        positions: &mut [Vec3],
        reference: &[Vec3],
        masses: &[f64],
        sim_box: &SimBox,
    ) -> Result<(), SimError> {
        for _ in 0..self.max_iter {
            let mut converged = true;
            for &(i, j, d) in &self.pairs {
                let s = sim_box.min_image(positions[i] - positions[j]);
                let s0 = sim_box.min_image(reference[i] - reference[j]);
                let diff = s.norm_sq() - d * d;
                if diff.abs() > 2.0 * self.tolerance * d * d {
                    converged = false;
                    let mu = 1.0 / masses[i] + 1.0 / masses[j];
                    let g = diff / (2.0 * s0.dot(s) * mu);
                    positions[i] -= s0 * (g / masses[i]);
                    positions[j] += s0 * (g / masses[j]);
                }
            }
            if converged {
                return Ok(());
            }
        }
        Err(SimError::Invalid(format!("SHAKE did not converge in {} iterations", self.max_iter)))
    }
}

/// LINCS-style constraint solver (sequential matrix expansion).
#[derive(Clone, Debug)]
pub struct Lincs {
    pairs: Vec<(usize, usize, f64)>,
    order: usize,
    iterations: usize,
}

impl Lincs {
    /// Creates a LINCS solver with `order` expansion terms (4 is typical) and
    /// `iterations` centripetal-correction passes (1 is typical).
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] for invalid pairs or parameters.
    pub fn new(
        pairs: Vec<(usize, usize, f64)>,
        n_atoms: usize,
        order: usize,
        iterations: usize,
    ) -> Result<Self, SimError> {
        validate_pairs(&pairs, n_atoms)?;
        if order == 0 {
            return Err(SimError::Invalid("LINCS expansion order must be >= 1".into()));
        }
        Ok(Self { pairs, order, iterations })
    }

    /// Solves `M λ = b` with `M` built from the old directions `u`, by the
    /// truncated expansion `λ = Σ_n Ã^n D⁻¹ b`.
    fn solve(&self, u: &[Vec3], masses: &[f64], b: &[f64]) -> Vec<f64> {
        let m = self.pairs.len();
        let diag: Vec<f64> =
            self.pairs.iter().map(|&(i, j, _)| 1.0 / masses[i] + 1.0 / masses[j]).collect();
        // Off-diagonal couplings between constraints sharing an atom.
        let mut couplings: Vec<Vec<(usize, f64)>> = vec![Vec::new(); m];
        for k in 0..m {
            let (ik, jk, _) = self.pairs[k];
            for l in 0..m {
                if l == k {
                    continue;
                }
                let (il, jl, _) = self.pairs[l];
                let mut mkl = 0.0;
                for (a, sa) in [(ik, 1.0), (jk, -1.0)] {
                    for (c, sc) in [(il, 1.0), (jl, -1.0)] {
                        if a == c {
                            mkl += sa * sc / masses[a];
                        }
                    }
                }
                if mkl != 0.0 {
                    couplings[k].push((l, -mkl * u[k].dot(u[l]) / diag[k]));
                }
            }
        }
        let mut term: Vec<f64> = (0..m).map(|k| b[k] / diag[k]).collect();
        let mut sum = term.clone();
        for _ in 0..self.order {
            let next: Vec<f64> =
                (0..m).map(|k| couplings[k].iter().map(|&(l, c)| c * term[l]).sum()).collect();
            for (s, t) in sum.iter_mut().zip(&next) {
                *s += t;
            }
            term = next;
        }
        sum
    }

    /// Applies the displacement `Σ_k λ_k σ u_k / m` to the positions, with the
    /// sign convention of the coupling matrix (`+` for the first atom of a pair).
    fn displace(&self, positions: &mut [Vec3], u: &[Vec3], masses: &[f64], lambda: &[f64]) {
        for (k, &(i, j, _)) in self.pairs.iter().enumerate() {
            let f = u[k] * lambda[k];
            positions[i] += f * (1.0 / masses[i]);
            positions[j] -= f * (1.0 / masses[j]);
        }
    }
}

impl Constraints for Lincs {
    fn name(&self) -> &'static str {
        "lincs"
    }

    fn pairs(&self) -> Vec<(usize, usize, f64)> {
        self.pairs.clone()
    }

    fn constrain_positions(
        &mut self,
        positions: &mut [Vec3],
        reference: &[Vec3],
        masses: &[f64],
        sim_box: &SimBox,
    ) -> Result<(), SimError> {
        let u: Vec<Vec3> = self
            .pairs
            .iter()
            .map(|&(i, j, _)| {
                let s0 = sim_box.min_image(reference[i] - reference[j]);
                s0 * (1.0 / s0.norm())
            })
            .collect();
        // Each pass solves for displacements along the old directions `u` that give
        // every constraint its target projection. Displacements along `u` leave the
        // component of the bond vector perpendicular to `u` unchanged, so the target
        // projection is sqrt(d² − |s_⊥|²): this captures the centripetal term.
        // The first pass is the projection step; the others are the correction passes.
        for _ in 0..=self.iterations {
            let b: Vec<f64> = self
                .pairs
                .iter()
                .enumerate()
                .map(|(k, &(i, j, d))| {
                    let s = sim_box.min_image(positions[i] - positions[j]);
                    let proj = u[k].dot(s);
                    let perp2 = (s.norm_sq() - proj * proj).max(0.0);
                    (d * d - perp2).max(0.0).sqrt() - proj
                })
                .collect();
            let lambda = self.solve(&u, masses, &b);
            self.displace(positions, &u, masses, &lambda);
        }
        Ok(())
    }
}

/// Rigid three-site water: atoms `[O, H1, H2]` with fixed O-H and H-H distances.
#[derive(Clone, Debug)]
pub struct Settle {
    molecules: Vec<[usize; 3]>,
    d_oh: f64,
    d_hh: f64,
    tolerance: f64,
}

/// Maximum Newton iterations for one water molecule.
const SETTLE_MAX_ITER: usize = 50;

impl Settle {
    /// Creates a SETTLE solver for the given molecules, each listed as
    /// `[oxygen, hydrogen1, hydrogen2]`.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] for invalid indices or non-positive distances.
    pub fn new(
        molecules: Vec<[usize; 3]>,
        n_atoms: usize,
        d_oh: f64,
        d_hh: f64,
        tolerance: f64,
    ) -> Result<Self, SimError> {
        if molecules.iter().flatten().any(|&a| a >= n_atoms) {
            return Err(SimError::Invalid("SETTLE molecule refers to a missing atom".into()));
        }
        if !(d_oh > 0.0 && d_hh > 0.0 && tolerance > 0.0) {
            return Err(SimError::Invalid("SETTLE distances and tolerance must be > 0".into()));
        }
        Ok(Self { molecules, d_oh, d_hh, tolerance })
    }

    /// The three constrained pairs of one molecule.
    fn molecule_pairs(&self, m: [usize; 3]) -> [(usize, usize, f64); 3] {
        [(m[0], m[1], self.d_oh), (m[0], m[2], self.d_oh), (m[1], m[2], self.d_hh)]
    }
}

impl Constraints for Settle {
    fn name(&self) -> &'static str {
        "settle"
    }

    fn pairs(&self) -> Vec<(usize, usize, f64)> {
        self.molecules.iter().flat_map(|&m| self.molecule_pairs(m)).collect()
    }

    fn constrain_positions(
        &mut self,
        positions: &mut [Vec3],
        reference: &[Vec3],
        masses: &[f64],
        sim_box: &SimBox,
    ) -> Result<(), SimError> {
        for &mol in &self.molecules {
            let pairs = self.molecule_pairs(mol);
            // Directions from the reference geometry, as in SHAKE.
            let s0: [Vec3; 3] =
                pairs.map(|(a, b, _)| sim_box.min_image(reference[a] - reference[b]));
            let mut eps = [0.0_f64; 3];
            let mut converged = false;
            for _ in 0..SETTLE_MAX_ITER {
                // Current positions for this multiplier guess.
                let mut x = [positions[mol[0]], positions[mol[1]], positions[mol[2]]];
                for (c, &(a, b, _)) in pairs.iter().enumerate() {
                    let ia = mol.iter().position(|&m| m == a).unwrap_or(0);
                    let ib = mol.iter().position(|&m| m == b).unwrap_or(0);
                    let f = s0[c] * eps[c];
                    x[ia] += f * (1.0 / masses[a]);
                    x[ib] -= f * (1.0 / masses[b]);
                }
                let mut resid = [0.0_f64; 3];
                let mut jac = [[0.0_f64; 3]; 3];
                // Atom `atom` moves by ε_c s0_c τ(atom, c), where τ = +1/m for the first
                // atom of constraint c, −1/m for the second and 0 otherwise.
                let tau = |atom: usize, c: usize| -> f64 {
                    let (ca, cb, _) = pairs[c];
                    if atom == ca {
                        1.0 / masses[atom]
                    } else if atom == cb {
                        -1.0 / masses[atom]
                    } else {
                        0.0
                    }
                };
                for (r, &(a, b, d)) in pairs.iter().enumerate() {
                    let ia = mol.iter().position(|&m| m == a).unwrap_or(0);
                    let ib = mol.iter().position(|&m| m == b).unwrap_or(0);
                    let s = sim_box.min_image(x[ia] - x[ib]);
                    resid[r] = s.norm_sq() - d * d;
                    for c in 0..3 {
                        // ∂(x_a − x_b)/∂ε_c = s0_c (τ(a,c) − τ(b,c)).
                        let coef = s0[c] * (tau(a, c) - tau(b, c));
                        jac[r][c] = 2.0 * s.dot(coef);
                    }
                }
                let worst = resid.iter().fold(0.0_f64, |m, &x| m.max(x.abs()));
                let scale = pairs.iter().map(|p| p.2 * p.2).fold(0.0_f64, f64::max);
                if worst <= 2.0 * self.tolerance * scale {
                    converged = true;
                    break;
                }
                let delta = solve3(jac, resid)
                    .ok_or_else(|| SimError::Invalid("SETTLE Jacobian is singular".into()))?;
                for c in 0..3 {
                    eps[c] -= delta[c];
                }
            }
            if !converged {
                return Err(SimError::Invalid("SETTLE did not converge".into()));
            }
            // Commit the displacement for this molecule.
            for (c, &(a, b, _)) in pairs.iter().enumerate() {
                let f = s0[c] * eps[c];
                positions[a] += f * (1.0 / masses[a]);
                positions[b] -= f * (1.0 / masses[b]);
            }
        }
        Ok(())
    }
}

/// Solves the 3×3 system `a x = b` by Gaussian elimination with partial
/// pivoting. Returns `None` for a singular matrix.
fn solve3(mut a: [[f64; 3]; 3], mut b: [f64; 3]) -> Option<[f64; 3]> {
    for col in 0..3 {
        let pivot = (col..3).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..3 {
            let f = a[row][col] / a[col][col];
            for k in col..3 {
                a[row][k] -= f * a[col][k];
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = [0.0; 3];
    for row in (0..3).rev() {
        let mut s = b[row];
        for k in row + 1..3 {
            s -= a[row][k] * x[k];
        }
        x[row] = s / a[row][row];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    /// Positions, masses, pairs and box of the four-atom test chain.
    type Chain = (Vec<Vec3>, Vec<f64>, Vec<(usize, usize, f64)>, SimBox);

    use super::*;
    use crate::system::NewAtom;

    fn atom(p: Vec3, m: f64) -> NewAtom {
        NewAtom {
            name: "X".into(),
            residue_name: "R".into(),
            residue_id: 1,
            element: "C".into(),
            position: p,
            velocity: Vec3::ZERO,
            mass: m,
            charge: 0.0,
            atom_type: 0,
        }
    }

    fn max_length_error(pairs: &[(usize, usize, f64)], pos: &[Vec3], b: &SimBox) -> f64 {
        pairs
            .iter()
            .map(|&(i, j, d)| (b.min_image(pos[i] - pos[j]).norm() - d).abs() / d)
            .fold(0.0, f64::max)
    }

    /// A small bonded chain 0-1-2-3 with masses 1, 2, 3, 4 and a deterministic kick.
    fn chain() -> Chain {
        let b = SimBox::new(10.0, 10.0, 10.0).unwrap();
        let reference = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.15, 1.0, 1.0),
            Vec3::new(1.25, 1.14, 1.0),
            Vec3::new(1.2, 1.22, 1.12),
        ];
        let pairs = vec![(0, 1, 0.15), (1, 2, 0.15), (2, 3, 0.15)];
        let masses = vec![1.0, 2.0, 3.0, 4.0];
        // Perturb every atom by a fixed, non-symmetric amount.
        let kicks = [
            Vec3::new(0.02, -0.01, 0.015),
            Vec3::new(-0.015, 0.02, -0.01),
            Vec3::new(0.01, 0.012, -0.02),
            Vec3::new(-0.02, -0.015, 0.01),
        ];
        let pos: Vec<Vec3> = reference.iter().zip(kicks).map(|(&p, k)| p + k).collect();
        (pos, masses, pairs, b)
    }

    #[test]
    fn shake_satisfies_all_bond_lengths() {
        let (mut pos, masses, pairs, b) = chain();
        let reference = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.15, 1.0, 1.0),
            Vec3::new(1.25, 1.14, 1.0),
            Vec3::new(1.2, 1.22, 1.12),
        ];
        let mut shake = Shake::new(pairs.clone(), 4, 1e-10, 500).unwrap();
        shake.constrain_positions(&mut pos, &reference, &masses, &b).unwrap();
        assert!(max_length_error(&pairs, &pos, &b) < 1e-9);
    }

    #[test]
    fn lincs_meets_its_documented_accuracy_on_a_chain() {
        let (mut pos, masses, pairs, b) = chain();
        let reference = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.15, 1.0, 1.0),
            Vec3::new(1.25, 1.14, 1.0),
            Vec3::new(1.2, 1.22, 1.12),
        ];
        let before = max_length_error(&pairs, &pos, &b);
        let mut lincs = Lincs::new(pairs.clone(), 4, 4, 2).unwrap();
        lincs.constrain_positions(&mut pos, &reference, &masses, &b).unwrap();
        let after = max_length_error(&pairs, &pos, &b);
        // LINCS is not exact for coupled constraints; require a large reduction.
        assert!(after < 1e-3 * before.max(1e-3), "before {before}, after {after}");
    }

    #[test]
    fn lincs_is_exact_for_isolated_bonds() {
        let b = SimBox::new(10.0, 10.0, 10.0).unwrap();
        let reference = [Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.1, 1.0, 1.0)];
        let mut pos = vec![Vec3::new(1.02, 1.01, 1.0), Vec3::new(1.13, 0.99, 1.02)];
        let pairs = vec![(0, 1, 0.1)];
        let masses = vec![1.0, 1.0];
        // A single constraint has no coupling, so only the centripetal correction
        // matters. Each pass reduces the length error; with a 10 % initial error
        // and four passes the residual is far below 1e-6.
        let mut lincs = Lincs::new(pairs.clone(), 2, 4, 4).unwrap();
        lincs.constrain_positions(&mut pos, &reference, &masses, &b).unwrap();
        assert!(
            max_length_error(&pairs, &pos, &b) < 1e-6,
            "{}",
            max_length_error(&pairs, &pos, &b)
        );
    }

    #[test]
    fn settle_matches_converged_shake_on_a_water_molecule() {
        // Equilibrium-like water in the xy plane, with an off-plane kick.
        let b = SimBox::new(10.0, 10.0, 10.0).unwrap();
        let d_oh = 0.1;
        let d_hh = 0.16330;
        let reference = [
            Vec3::new(5.0, 5.0, 5.0),
            Vec3::new(5.0 + d_oh * 0.8165, 5.0 + d_oh * 0.5774, 5.0),
            Vec3::new(5.0 - d_oh * 0.8165, 5.0 + d_oh * 0.5774, 5.0),
        ];
        let masses = vec![15.9994, 1.008, 1.008];
        let kicked: Vec<Vec3> = reference
            .iter()
            .enumerate()
            .map(|(i, &p)| {
                p + Vec3::new(0.004 * (i as f64 - 1.0), -0.003, 0.005 * (1.5 - i as f64))
            })
            .collect();
        let pairs = vec![(0, 1, d_oh), (0, 2, d_oh), (1, 2, d_hh)];

        let mut settle_pos = kicked.clone();
        let mut settle = Settle::new(vec![[0, 1, 2]], 3, d_oh, d_hh, 1e-12).unwrap();
        settle.constrain_positions(&mut settle_pos, &reference, &masses, &b).unwrap();

        let mut shake_pos = kicked;
        let mut shake = Shake::new(pairs.clone(), 3, 1e-13, 10_000).unwrap();
        shake.constrain_positions(&mut shake_pos, &reference, &masses, &b).unwrap();

        assert!(max_length_error(&pairs, &settle_pos, &b) < 1e-9);
        for (s, h) in settle_pos.iter().zip(&shake_pos) {
            assert!((*s - *h).norm() < 1e-9, "SETTLE and SHAKE differ: {s:?} vs {h:?}");
        }
    }

    #[test]
    fn settle_keeps_water_geometry_over_a_dynamics_like_loop() {
        let b = SimBox::new(10.0, 10.0, 10.0).unwrap();
        let (d_oh, d_hh) = (0.1, 0.16330);
        let mut pos = vec![
            Vec3::new(5.0, 5.0, 5.0),
            Vec3::new(5.0 + d_oh * 0.8165, 5.0 + d_oh * 0.5774, 5.0),
            Vec3::new(5.0 - d_oh * 0.8165, 5.0 + d_oh * 0.5774, 5.0),
        ];
        let masses = vec![15.9994, 1.008, 1.008];
        let mut settle = Settle::new(vec![[0, 1, 2]], 3, d_oh, d_hh, 1e-10).unwrap();
        let pairs = settle.pairs();
        let mut seed = 12345u64;
        let mut rnd = || {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((seed >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        for _ in 0..200 {
            let reference = pos.clone();
            for p in &mut pos {
                *p += Vec3::new(rnd(), rnd(), rnd()) * 0.01;
            }
            settle.constrain_positions(&mut pos, &reference, &masses, &b).unwrap();
            assert!(max_length_error(&pairs, &pos, &b) < 1e-8);
        }
    }

    #[test]
    fn velocity_projection_removes_bond_components() {
        let mut sys = System::new(SimBox::new(10.0, 10.0, 10.0).unwrap());
        sys.push_atom(atom(Vec3::new(1.0, 1.0, 1.0), 1.0));
        sys.push_atom(atom(Vec3::new(1.1, 1.0, 1.0), 2.0));
        sys.velocities[0] = Vec3::new(1.0, 0.5, 0.0);
        sys.velocities[1] = Vec3::new(-0.5, 0.2, 0.3);
        project_velocities(&[(0, 1, 0.1)], &mut sys, 1e-12, 50).unwrap();
        let r = sys.positions[0] - sys.positions[1];
        let rdv = r.dot(sys.velocities[0] - sys.velocities[1]);
        assert!(rdv.abs() < 1e-10);
    }

    #[test]
    fn apply_constraints_corrects_velocities_by_displacement() {
        let mut sys = System::new(SimBox::new(10.0, 10.0, 10.0).unwrap());
        sys.push_atom(atom(Vec3::new(1.0, 1.0, 1.0), 1.0));
        sys.push_atom(atom(Vec3::new(1.1, 1.0, 1.0), 1.0));
        let reference = sys.positions.clone();
        sys.positions[1] = Vec3::new(1.2, 1.0, 1.0); // stretched by 0.1 nm
        let mut shake = Shake::new(vec![(0, 1, 0.1)], 2, 1e-12, 100).unwrap();
        apply_constraints(&mut shake, &mut sys, &reference, 0.5).unwrap();
        let len = (sys.positions[1] - sys.positions[0]).norm();
        assert!((len - 0.1).abs() < 1e-10);
        // Atom 1 moved back by 0.05 nm total (split equally); velocity change = shift/dt.
        assert!((sys.velocities[1].x + 0.05 / 0.5).abs() < 1e-9);
    }

    #[test]
    fn solve3_matches_known_solution() {
        let a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 4.0]];
        let x = [1.0, -2.0, 0.5];
        let b = [2.0 * x[0] + x[1], x[0] + 3.0 * x[1] + x[2], x[1] + 4.0 * x[2]];
        let sol = solve3(a, b).unwrap();
        for k in 0..3 {
            assert!((sol[k] - x[k]).abs() < 1e-12);
        }
    }
}
