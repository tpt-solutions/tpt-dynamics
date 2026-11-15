//! Time integrators behind a runtime-swappable trait, and the force evaluator
//! that couples a [`ForceField`] to a [`VerletList`].
//!
//! Provided integrators:
//! - [`VelocityVerlet`]: positions, velocities and forces are all at the same
//!   time; kinetic energy is exact at each step.
//! - [`LeapFrog`]: velocities live at half steps. Kinetic energy at integer
//!   steps is the average of the two neighbouring half-step kinetic energies.
//!   The first call to [`Integrator::prepare`] shifts the initial velocities
//!   back by half a step.

use tpt_dynamics_core::Vec3;

use crate::constraints::{apply_constraints, Constraints};
use crate::error::SimError;
use crate::forcefield::{ForceField, ForceOutput};
use crate::neighbour::VerletList;
use crate::system::System;

/// Couples a force field with a neighbour list and remembers the latest result.
pub struct ForceEvaluator {
    forcefield: Box<dyn ForceField>,
    neighbours: VerletList,
    last: ForceOutput,
    /// One-thread pool used in deterministic mode; `None` uses the global pool.
    pool: Option<rayon::ThreadPool>,
}

impl ForceEvaluator {
    /// Creates an evaluator. The neighbour-list cutoff is the force field's
    /// cutoff; `skin` is the Verlet buffer width in nm.
    ///
    /// # Errors
    /// Returns an error if the skin or cutoff are invalid.
    pub fn new(forcefield: Box<dyn ForceField>, skin: f64) -> Result<Self, SimError> {
        let neighbours = VerletList::new(forcefield.cutoff(), skin)?;
        Ok(Self { forcefield, neighbours, last: ForceOutput::default(), pool: None })
    }

    /// Enables or disables reproducibility mode. When enabled, the neighbour
    /// list is built serially and every force evaluation runs inside a
    /// one-thread rayon pool, so that rayon-based kernels in dependencies (such
    /// as the PME spreading, which chunks by thread count) take one fixed path.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if the thread pool cannot be created.
    pub fn set_deterministic(&mut self, deterministic: bool) -> Result<(), SimError> {
        if deterministic {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(1)
                .build()
                .map_err(|e| SimError::Invalid(format!("cannot build thread pool: {e}")))?;
            self.pool = Some(pool);
            self.neighbours.set_parallel(false);
        } else {
            self.pool = None;
            self.neighbours.set_parallel(true);
        }
        Ok(())
    }

    /// True if reproducibility mode is active.
    pub fn is_deterministic(&self) -> bool {
        self.pool.is_some()
    }

    /// The force field in use.
    pub fn forcefield(&self) -> &dyn ForceField {
        self.forcefield.as_ref()
    }

    /// The Verlet list (for diagnostics such as rebuild counts).
    pub fn neighbours(&self) -> &VerletList {
        &self.neighbours
    }

    /// Potential energy and virial from the most recent evaluation.
    pub fn last(&self) -> ForceOutput {
        self.last
    }

    /// Forces the neighbour list to be rebuilt on the next evaluation (used after
    /// the box changes size, and before writing a checkpoint).
    pub fn invalidate_neighbours(&mut self) {
        self.neighbours.invalidate();
    }

    /// Updates the neighbour list if needed, recomputes `system.forces` and
    /// returns the energies for the current positions.
    ///
    /// # Errors
    /// Returns an error if the neighbour list cannot be built for this box.
    pub fn evaluate(&mut self, system: &mut System) -> Result<ForceOutput, SimError> {
        let Self { forcefield, neighbours, last, pool } = self;
        let out = match pool.as_ref() {
            Some(p) => p.install(|| evaluate_inner(neighbours, forcefield.as_ref(), system))?,
            None => evaluate_inner(neighbours, forcefield.as_ref(), system)?,
        };
        *last = out;
        Ok(out)
    }
}

fn evaluate_inner(
    neighbours: &mut VerletList,
    forcefield: &dyn ForceField,
    system: &mut System,
) -> Result<ForceOutput, SimError> {
    neighbours.update(&system.positions, &system.sim_box)?;
    system.forces.fill(Vec3::ZERO);
    Ok(forcefield.compute(
        &system.positions,
        &system.atom_types,
        &system.charges,
        &system.sim_box,
        neighbours,
        &mut system.forces,
    ))
}

/// A time-stepping scheme. Implementations are chosen at runtime (see
/// [`crate::config::IntegratorKind`]).
pub trait Integrator: Send {
    /// Short identifier used in logs.
    fn name(&self) -> &'static str;

    /// One-time setup before the first step. On entry `system.forces` are
    /// current for `system.positions`.
    fn prepare(&mut self, _system: &mut System, _dt: f64) {}

    /// Advances positions and velocities by `dt`, re-evaluating forces through
    /// `eval` at the new positions. `constraints` is applied to the positions
    /// after the unconstrained drift and before the force evaluation.
    ///
    /// # Errors
    /// Propagates neighbour-list and constraint-solver errors.
    fn step(
        &mut self,
        system: &mut System,
        eval: &mut ForceEvaluator,
        dt: f64,
        constraints: &mut dyn Constraints,
    ) -> Result<(), SimError>;

    /// Kinetic energy consistent with the current state, for energy output.
    fn kinetic_energy(&self, system: &System, _dt: f64) -> f64 {
        system.kinetic_energy()
    }
}

fn kick(system: &mut System, factor: f64) {
    for ((v, f), &m) in system.velocities.iter_mut().zip(&system.forces).zip(&system.masses) {
        *v += *f * (factor / m);
    }
}

fn drift(system: &mut System, dt: f64) {
    for (p, v) in system.positions.iter_mut().zip(&system.velocities) {
        *p += *v * dt;
    }
    for p in &mut system.positions {
        *p = system.sim_box.wrap(*p);
    }
}

/// Velocity Verlet: `v += dt/2 F/m; x += dt v; (new F); v += dt/2 F/m`.
#[derive(Clone, Copy, Debug, Default)]
pub struct VelocityVerlet;

impl Integrator for VelocityVerlet {
    fn name(&self) -> &'static str {
        "velocity-verlet"
    }

    fn step(
        &mut self,
        system: &mut System,
        eval: &mut ForceEvaluator,
        dt: f64,
        constraints: &mut dyn Constraints,
    ) -> Result<(), SimError> {
        kick(system, 0.5 * dt);
        let reference = constraints.is_active().then(|| system.positions.clone());
        drift(system, dt);
        if let Some(reference) = reference {
            apply_constraints(constraints, system, &reference, dt)?;
        }
        eval.evaluate(system)?;
        kick(system, 0.5 * dt);
        Ok(())
    }
}

/// Leap-frog: `v(t+dt/2) = v(t-dt/2) + dt F(t)/m; x(t+dt) = x(t) + dt v(t+dt/2)`.
#[derive(Clone, Copy, Debug, Default)]
pub struct LeapFrog;

impl Integrator for LeapFrog {
    fn name(&self) -> &'static str {
        "leap-frog"
    }

    fn prepare(&mut self, system: &mut System, dt: f64) {
        // Interpret the input velocities as v(0) and move them to v(-dt/2).
        kick(system, -0.5 * dt);
    }

    fn step(
        &mut self,
        system: &mut System,
        eval: &mut ForceEvaluator,
        dt: f64,
        constraints: &mut dyn Constraints,
    ) -> Result<(), SimError> {
        kick(system, dt);
        let reference = constraints.is_active().then(|| system.positions.clone());
        drift(system, dt);
        if let Some(reference) = reference {
            apply_constraints(constraints, system, &reference, dt)?;
        }
        eval.evaluate(system)?;
        Ok(())
    }

    fn kinetic_energy(&self, system: &System, dt: f64) -> f64 {
        // Half-step KE at t-dt/2 is the current velocities; the half-step KE at
        // t+dt/2 is exactly what the next kick will produce.
        let ke_before = system.kinetic_energy();
        let ke_after: f64 = system
            .masses
            .iter()
            .zip(&system.velocities)
            .zip(&system.forces)
            .map(|((&m, &v), &f)| {
                let vn = v + f * (dt / m);
                0.5 * m * vn.norm_sq()
            })
            .sum();
        f64::midpoint(ke_before, ke_after)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::constraints::NoConstraints;
    use crate::forcefield::{LennardJones, Modifier};
    use crate::pbc::SimBox;
    use crate::system::NewAtom;

    /// Simple cubic lattice of `n_side³` argon-like atoms with small random
    /// displacements and random velocities (zero total momentum).
    pub(crate) fn lj_gas(n_side: usize, density: f64, seed: u64) -> System {
        let a = (1.0 / density).cbrt();
        let l = a * n_side as f64;
        let mut sys = System::new(SimBox::new(l, l, l).unwrap());
        let mut s = seed | 1;
        let mut rnd = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for ix in 0..n_side {
            for iy in 0..n_side {
                for iz in 0..n_side {
                    let p = Vec3::new(
                        (ix as f64 + 0.5) * a + 0.05 * a * rnd(),
                        (iy as f64 + 0.5) * a + 0.05 * a * rnd(),
                        (iz as f64 + 0.5) * a + 0.05 * a * rnd(),
                    );
                    sys.push_atom(NewAtom {
                        name: "AR".into(),
                        residue_name: "LJ".into(),
                        residue_id: 1,
                        element: "Ar".into(),
                        position: p,
                        velocity: Vec3::new(rnd(), rnd(), rnd()) * 0.8,
                        mass: 1.0,
                        charge: 0.0,
                        atom_type: 0,
                    });
                }
            }
        }
        sys.remove_com_velocity();
        sys
    }

    fn evaluator(cutoff: f64, skin: f64) -> ForceEvaluator {
        let ff = LennardJones::single_type(1.0, 1.0, cutoff, Modifier::Shift).unwrap();
        ForceEvaluator::new(Box::new(ff), skin).unwrap()
    }

    fn run_nve(integ: &mut dyn Integrator, dt: f64, steps: usize) -> (f64, f64, f64) {
        let mut sys = lj_gas(6, 0.6, 11);
        let mut eval = evaluator(2.5, 0.3);
        let out = eval.evaluate(&mut sys).unwrap();
        integ.prepare(&mut sys, dt);
        let e0 = integ.kinetic_energy(&sys, dt) + out.potential;
        let mut max_dev: f64 = 0.0;
        for _ in 0..steps {
            integ.step(&mut sys, &mut eval, dt, &mut NoConstraints).unwrap();
            let e = integ.kinetic_energy(&sys, dt) + eval.last().potential;
            max_dev = max_dev.max((e - e0).abs());
            assert!(sys.total_momentum().norm() < 1e-9, "momentum drifted");
        }
        (e0, max_dev, sys.total_momentum().norm())
    }

    #[test]
    fn velocity_verlet_conserves_energy_and_momentum() {
        let (e0, dev, _) = run_nve(&mut VelocityVerlet, 0.002, 500);
        assert!(dev / e0.abs().max(1.0) < 1e-3, "energy deviation {dev} (E0 {e0})");
    }

    #[test]
    fn leap_frog_conserves_energy_and_momentum() {
        let (e0, dev, _) = run_nve(&mut LeapFrog, 0.002, 500);
        assert!(dev / e0.abs().max(1.0) < 1e-3, "energy deviation {dev} (E0 {e0})");
    }

    #[test]
    fn integrators_are_runtime_swappable() {
        let mut boxed: Vec<Box<dyn Integrator>> =
            vec![Box::new(VelocityVerlet), Box::new(LeapFrog)];
        let names: Vec<&str> = boxed.iter().map(|b| b.name()).collect();
        assert_eq!(names, ["velocity-verlet", "leap-frog"]);
        for b in &mut boxed {
            let _ = run_nve(b.as_mut(), 0.001, 10);
        }
    }
}
