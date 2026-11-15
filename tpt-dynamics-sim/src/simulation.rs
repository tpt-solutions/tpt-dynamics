//! The run loop: owns the system, the force evaluator, the integrator and the
//! optional thermostat; computes thermodynamic observables and dispatches
//! output to [`Observer`]s.

use std::time::Instant;

use tpt_dynamics_core::Vec3;

use crate::barostat::Barostat;
use crate::constraints::{Constraints, NoConstraints};
use crate::error::SimError;
use crate::integrator::{ForceEvaluator, Integrator};
use crate::output::{Energies, Observation, Observer};
use crate::rng::Rng;
use crate::system::{System, UnitSystem};
use crate::thermostat::Thermostat;

/// A configured simulation ready to advance.
pub struct Simulation {
    /// Unit system (sets the Boltzmann constant).
    pub units: UnitSystem,
    /// The particle state.
    pub system: System,
    /// Force evaluation and neighbour list.
    pub eval: ForceEvaluator,
    /// Time integrator.
    pub integrator: Box<dyn Integrator>,
    /// Optional temperature control.
    pub thermostat: Option<Box<dyn Thermostat>>,
    /// Optional isotropic pressure coupling.
    pub barostat: Option<Box<dyn Barostat>>,
    /// Holonomic constraints (none by default).
    pub constraints: Box<dyn Constraints>,
    /// Time step in ps.
    pub dt: f64,
    /// Number of completed steps.
    pub step: u64,
}

/// Summary of a [`Simulation::run`] call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RunSummary {
    /// Steps performed by this call.
    pub steps: u64,
    /// Wall-clock time of the call in seconds.
    pub wall_seconds: f64,
    /// Neighbour-list rebuilds during this call.
    pub neighbour_rebuilds: u64,
    /// Observables at the start of the call.
    pub initial: Energies,
    /// Observables at the end of the call.
    pub last: Energies,
}

impl RunSummary {
    /// Relative change of the total energy over the run, `(E_end − E_0) / |E_0|`.
    pub fn relative_energy_drift(&self) -> f64 {
        let e0 = self.initial.total;
        (self.last.total - e0) / e0.abs().max(f64::MIN_POSITIVE)
    }
}

impl Simulation {
    /// Assembles a simulation, checks the inputs, and evaluates initial forces.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] for an empty system, inconsistent arrays,
    /// non-positive masses, atom types without force-field parameters, a
    /// non-positive time step, or a box too small for the list radius.
    pub fn new(
        mut system: System,
        mut eval: ForceEvaluator,
        mut integrator: Box<dyn Integrator>,
        thermostat: Option<Box<dyn Thermostat>>,
        dt: f64,
        units: UnitSystem,
    ) -> Result<Self, SimError> {
        if system.is_empty() {
            return Err(SimError::Invalid("system has no atoms".into()));
        }
        system.validate()?;
        eval.forcefield().validate_types(&system.atom_types)?;
        eval.forcefield().validate_box(&system.sim_box)?;
        if !(dt.is_finite() && dt > 0.0) {
            return Err(SimError::Invalid(format!("time step must be > 0, got {dt}")));
        }
        eval.evaluate(&mut system)?;
        integrator.prepare(&mut system, dt);
        Ok(Self {
            units,
            system,
            eval,
            integrator,
            thermostat,
            barostat: None,
            constraints: Box::new(NoConstraints),
            dt,
            step: 0,
        })
    }

    /// Adds holonomic constraints. The initial positions are projected onto the
    /// constraint surface, the initial velocities are made tangent to it, and
    /// the constraint count is removed from the degrees of freedom.
    ///
    /// # Errors
    /// Propagates solver errors from the initial projection.
    pub fn with_constraints(
        mut self,
        mut constraints: Box<dyn Constraints>,
    ) -> Result<Self, SimError> {
        let pairs = constraints.pairs();
        let before = self.system.positions.clone();
        let sim_box = self.system.sim_box;
        if constraints.is_active() {
            constraints.constrain_positions(
                &mut self.system.positions,
                &before,
                &self.system.masses,
                &sim_box,
            )?;
            crate::constraints::project_velocities(&pairs, &mut self.system, 1e-10, 500)?;
        }
        self.system.n_constraints = pairs.len();
        self.constraints = constraints;
        self.eval.evaluate(&mut self.system)?;
        Ok(self)
    }

    /// Adds isotropic pressure coupling. Pressure is in the engine's unit
    /// (see [`crate::barostat`]).
    #[must_use]
    pub fn with_barostat(mut self, barostat: Box<dyn Barostat>) -> Self {
        self.barostat = Some(barostat);
        self
    }

    /// Name of the active integrator.
    pub fn integrator_name(&self) -> &'static str {
        self.integrator.name()
    }

    /// Current simulation time in ps.
    pub fn time(&self) -> f64 {
        self.step as f64 * self.dt
    }

    /// Thermodynamic observables for the current state.
    pub fn energies(&self) -> Energies {
        let out = self.eval.last();
        let kinetic = self.integrator.kinetic_energy(&self.system, self.dt);
        let dof = self.system.degrees_of_freedom();
        let temperature =
            if dof > 0.0 { 2.0 * kinetic / (dof * self.units.boltzmann()) } else { 0.0 };
        let volume = self.system.sim_box.volume();
        Energies {
            step: self.step,
            time: self.time(),
            kinetic,
            potential: out.potential,
            total: kinetic + out.potential,
            temperature,
            pressure: (2.0 * kinetic + out.virial) / (3.0 * volume),
        }
    }

    /// Advances the system by one time step (integrator, then thermostat).
    ///
    /// # Errors
    /// Propagates neighbour-list errors.
    pub fn advance(&mut self) -> Result<(), SimError> {
        self.integrator.step(
            &mut self.system,
            &mut self.eval,
            self.dt,
            self.constraints.as_mut(),
        )?;
        if let Some(th) = self.thermostat.as_mut() {
            let ke = self.integrator.kinetic_energy(&self.system, self.dt);
            th.apply(&mut self.system, ke, self.dt);
        }
        if self.barostat.is_some() {
            let e = self.energies();
            let kt = self.units.boltzmann() * e.temperature;
            let volume = self.system.sim_box.volume();
            let dt = self.dt;
            if let Some(bar) = self.barostat.as_mut() {
                let mu = bar.scale_factor(e.pressure, volume, kt, dt);
                self.scale_cell(mu);
            }
            // Forces refer to the old coordinates; refresh them for the next step.
            self.eval.evaluate(&mut self.system)?;
        }
        self.step += 1;
        Ok(())
    }

    /// Scales positions, box lengths and tilts by `mu` and invalidates the
    /// neighbour list.
    fn scale_cell(&mut self, mu: f64) {
        // `mu` is clamped to [0.97, 1.03] by the barostats, so the box stays valid
        // and the reduced-tilt condition is preserved.
        let sb = &mut self.system.sim_box;
        sb.lengths *= mu;
        for t in &mut sb.tilt {
            *t *= mu;
        }
        let sb = self.system.sim_box;
        for p in &mut self.system.positions {
            *p = sb.wrap(*p * mu);
        }
        self.eval.invalidate_neighbours();
    }

    /// Runs `nsteps` steps, calling each observer at its interval (including
    /// step 0 and the final step if they fall on the interval), then flushes
    /// the observers.
    ///
    /// # Errors
    /// Propagates simulation and output errors.
    pub fn run(
        &mut self,
        nsteps: u64,
        observers: &mut [Box<dyn Observer>],
    ) -> Result<RunSummary, SimError> {
        let started = Instant::now();
        let builds_before = self.eval.neighbours().builds();
        let initial = self.energies();
        self.observe(observers)?;
        for _ in 0..nsteps {
            self.advance()?;
            self.observe(observers)?;
        }
        for o in observers.iter_mut() {
            o.finish()?;
        }
        Ok(RunSummary {
            steps: nsteps,
            wall_seconds: started.elapsed().as_secs_f64(),
            neighbour_rebuilds: self.eval.neighbours().builds() - builds_before,
            initial,
            last: self.energies(),
        })
    }

    /// Serialises the full simulation state (see [`crate::checkpoint`]).
    ///
    /// The neighbour list is rebuilt and forces re-evaluated first, so that a
    /// run restored from the checkpoint follows exactly the same floating-point
    /// path as the run that wrote it.
    ///
    /// # Errors
    /// Propagates force-evaluation errors.
    pub fn checkpoint_bytes(&mut self) -> Result<Vec<u8>, SimError> {
        self.eval.invalidate_neighbours();
        self.eval.evaluate(&mut self.system)?;
        Ok(crate::checkpoint::encode(self))
    }

    /// Restores the state written by [`Simulation::checkpoint_bytes`] into this
    /// simulation, which must have been built from the same system and settings.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if the bytes are corrupt, or if the atom
    /// count, time step or thermostat/barostat kinds do not match.
    pub fn restore_checkpoint(&mut self, bytes: &[u8]) -> Result<(), SimError> {
        let ck = crate::checkpoint::Checkpoint::decode(bytes)?;
        ck.apply(self)
    }

    fn observe(&self, observers: &mut [Box<dyn Observer>]) -> Result<(), SimError> {
        let due = |o: &dyn Observer| self.step.is_multiple_of(o.interval().max(1));
        if !observers.iter().any(|o| due(o.as_ref())) {
            return Ok(());
        }
        let energies = self.energies();
        let obs = Observation { system: &self.system, energies };
        for o in observers.iter_mut() {
            if self.step.is_multiple_of(o.interval().max(1)) {
                o.observe(&obs)?;
            }
        }
        Ok(())
    }
}

/// Draws Maxwell-Boltzmann velocities at temperature `t` (K), removes the
/// centre-of-mass motion, and rescales so the instantaneous temperature is
/// exactly `t`. Deterministic for a given `seed`.
pub fn maxwell_boltzmann_velocities(system: &mut System, t: f64, seed: u64, units: UnitSystem) {
    let mut rng = Rng::seed_from_u64(seed);
    if t <= 0.0 || system.is_empty() {
        system.velocities.fill(Vec3::ZERO);
        return;
    }
    for (v, &m) in system.velocities.iter_mut().zip(&system.masses) {
        let sigma = (units.boltzmann() * t / m).sqrt();
        *v = Vec3::new(rng.gaussian(), rng.gaussian(), rng.gaussian()) * sigma;
    }
    system.remove_com_velocity();
    let t_now = system.temperature_in(units);
    if t_now > 0.0 {
        let scale = (t / t_now).sqrt();
        for v in &mut system.velocities {
            *v *= scale;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::{LennardJones, Modifier};
    use crate::integrator::{LeapFrog, VelocityVerlet};
    use crate::output::Observer;
    use crate::thermostat::Berendsen;

    fn make(
        integrator: Box<dyn Integrator>,
        thermostat: Option<Box<dyn Thermostat>>,
        dt: f64,
    ) -> Simulation {
        let mut sys = crate::integrator::tests::lj_gas(6, 0.7, 21);
        maxwell_boltzmann_velocities(&mut sys, 0.7, 5, UnitSystem::Reduced);
        let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap();
        let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
        Simulation::new(sys, eval, integrator, thermostat, dt, UnitSystem::Reduced).unwrap()
    }

    #[test]
    fn velocities_have_target_temperature_and_zero_momentum() {
        let mut sys = crate::integrator::tests::lj_gas(5, 0.5, 2);
        maxwell_boltzmann_velocities(&mut sys, 1.3, 77, UnitSystem::Reduced);
        assert!((sys.temperature_in(UnitSystem::Reduced) - 1.3).abs() < 1e-9);
        assert!(sys.total_momentum().norm() < 1e-9);
        // Same seed gives the same velocities.
        let mut again = crate::integrator::tests::lj_gas(5, 0.5, 2);
        maxwell_boltzmann_velocities(&mut again, 1.3, 77, UnitSystem::Reduced);
        assert_eq!(sys.velocities, again.velocities);
    }

    #[test]
    fn nve_run_conserves_momentum_and_energy() {
        for integ in [Box::new(VelocityVerlet) as Box<dyn Integrator>, Box::new(LeapFrog)] {
            let mut sim = make(integ, None, 0.001);
            let summary = sim.run(1000, &mut []).unwrap();
            assert!(sim.system.total_momentum().norm() < 1e-9);
            assert!(
                summary.relative_energy_drift().abs() < 1e-3,
                "{}: drift {}",
                sim.integrator_name(),
                summary.relative_energy_drift()
            );
            assert!(summary.last.temperature.is_finite());
        }
    }

    #[test]
    fn berendsen_run_reaches_target_temperature() {
        let th = Berendsen::new(0.4, 0.1, UnitSystem::Reduced).unwrap();
        let mut sim = make(Box::new(VelocityVerlet), Some(Box::new(th)), 0.002);
        sim.run(1500, &mut []).unwrap();
        let t = sim.energies().temperature;
        assert!((t - 0.4).abs() < 0.08, "T = {t}");
    }

    #[test]
    fn observers_are_called_at_their_interval() {
        struct Count(u64, std::sync::Arc<std::sync::Mutex<Vec<u64>>>);
        impl Observer for Count {
            fn interval(&self) -> u64 {
                self.0
            }
            fn observe(&mut self, o: &Observation<'_>) -> Result<(), SimError> {
                self.1.lock().unwrap().push(o.energies.step);
                Ok(())
            }
        }
        let seen = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let mut obs: Vec<Box<dyn Observer>> = vec![Box::new(Count(4, seen.clone()))];
        let mut sim = make(Box::new(VelocityVerlet), None, 0.001);
        sim.run(10, &mut obs).unwrap();
        assert_eq!(*seen.lock().unwrap(), vec![0, 4, 8]);
    }

    #[test]
    fn rejects_bad_inputs() {
        let sys = System::new(crate::pbc::SimBox::new(5.0, 5.0, 5.0).unwrap());
        let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap();
        let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
        assert!(Simulation::new(sys, eval, Box::new(VelocityVerlet), None, 0.001, UnitSystem::Md)
            .is_err());
    }
}
