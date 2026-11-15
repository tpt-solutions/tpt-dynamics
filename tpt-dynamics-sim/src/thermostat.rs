//! Temperature control. Currently the Berendsen weak-coupling thermostat.
//!
//! Berendsen rescales velocities by `λ = sqrt(1 + (dt/τ)(T₀/T − 1))` each step.
//! It does not sample the canonical ensemble exactly, but it drives the
//! temperature to `T₀` smoothly and is useful for equilibration.

use crate::error::SimError;
use crate::rng::Rng;
use crate::system::{System, UnitSystem};

/// A velocity-rescaling temperature control applied after each integrator step.
pub trait Thermostat: Send {
    /// Short identifier used in logs.
    fn name(&self) -> &'static str;

    /// Target temperature in K.
    fn target(&self) -> f64;

    /// Rescales `system.velocities` given the current kinetic energy.
    fn apply(&mut self, system: &mut System, kinetic: f64, dt: f64);

    /// Internal state as raw 64-bit words, for checkpoints. Stateless
    /// thermostats return an empty vector.
    fn save_state(&self) -> Vec<u64> {
        Vec::new()
    }

    /// Restores state produced by [`Thermostat::save_state`].
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if the state has the wrong length.
    fn load_state(&mut self, state: &[u64]) -> Result<(), SimError> {
        if state.is_empty() {
            Ok(())
        } else {
            Err(SimError::Invalid("thermostat has no state to restore".into()))
        }
    }
}

fn validate_thermostat(t_ref: f64, tau: f64) -> Result<(), SimError> {
    if !(t_ref.is_finite() && t_ref >= 0.0) {
        return Err(SimError::Config(format!("thermostat temperature must be >= 0, got {t_ref}")));
    }
    if !(tau.is_finite() && tau > 0.0) {
        return Err(SimError::Config(format!("thermostat tau must be > 0, got {tau}")));
    }
    Ok(())
}

/// Current temperature from kinetic energy with the given `k_B`.
fn temperature_from(system: &System, kinetic: f64, boltzmann: f64) -> Option<f64> {
    let dof = system.degrees_of_freedom();
    if dof <= 0.0 || kinetic <= 0.0 {
        None
    } else {
        Some(2.0 * kinetic / (dof * boltzmann))
    }
}

/// Stochastic velocity rescaling (Bussi, Donadio & Parrinello, J. Chem. Phys.
/// 126, 014101, 2007).
///
/// Each step the kinetic energy `K` is moved towards `K₀ = ½ N_f k_B T₀` and a
/// Gaussian noise term is added, so that the canonical distribution is sampled:
///
/// ```text
/// K' = K + (K₀ − K) Δt/τ + 2 √(K K₀ / N_f) √(Δt/τ) R,   R ~ N(0, 1)
/// ```
///
/// Velocities are scaled by `√(K'/K)`. This is the first-order
/// (Euler–Maruyama) discretisation of the paper's Langevin-type equation. The
/// exact per-step update in the paper's appendix is not implemented. The noise
/// comes from the crate's deterministic [`Rng`], seeded by the caller, so runs
/// are reproducible; the generator state is part of the checkpoint.
#[derive(Clone, Debug)]
pub struct VRescale {
    t_ref: f64,
    tau: f64,
    boltzmann: f64,
    rng: Rng,
}

impl VRescale {
    /// Creates the thermostat for target `t_ref`, relaxation time `tau` and a
    /// noise seed.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for a negative temperature or non-positive `tau`.
    pub fn new(t_ref: f64, tau: f64, units: UnitSystem, seed: u64) -> Result<Self, SimError> {
        validate_thermostat(t_ref, tau)?;
        Ok(Self { t_ref, tau, boltzmann: units.boltzmann(), rng: Rng::seed_from_u64(seed) })
    }
}

impl Thermostat for VRescale {
    fn name(&self) -> &'static str {
        "v-rescale"
    }

    fn target(&self) -> f64 {
        self.t_ref
    }

    fn apply(&mut self, system: &mut System, kinetic: f64, dt: f64) {
        // The noise is drawn every step, even when the kinetic energy is zero, so
        // the random stream does not depend on the trajectory.
        let r = self.rng.gaussian();
        let dof = system.degrees_of_freedom();
        if dof <= 0.0 || kinetic <= 0.0 {
            return;
        }
        let k0 = 0.5 * dof * self.boltzmann * self.t_ref;
        let h = dt / self.tau;
        let k_new = kinetic + (k0 - kinetic) * h + 2.0 * (kinetic * k0 / dof).sqrt() * h.sqrt() * r;
        let scale = (k_new.max(0.0) / kinetic).sqrt();
        for v in &mut system.velocities {
            *v *= scale;
        }
    }

    fn save_state(&self) -> Vec<u64> {
        self.rng.state().to_vec()
    }

    fn load_state(&mut self, state: &[u64]) -> Result<(), SimError> {
        match state {
            [a, b, c, d] => {
                self.rng = Rng::from_state([*a, *b, *c, *d]);
                Ok(())
            }
            _ => Err(SimError::Invalid("v-rescale state must have 4 words".into())),
        }
    }
}

/// Nosé–Hoover thermostat with a chain length of one.
///
/// The friction variable `ξ` obeys `dξ/dt = (T/T₀ − 1)·4π²/τ²` and the velocities
/// obey `dv/dt = −ξ v`, so the Nosé–Hoover equations sample the canonical
/// ensemble with period `τ`. The implementation applies an explicit first-order
/// update after each integrator step: `ξ` is advanced by the temperature
/// deviation, then velocities are scaled by `exp(−ξ Δt)`. The conserved
/// energy is not tracked, and only chain length one is provided. `ξ` is part
/// of the checkpoint state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NoseHoover {
    t_ref: f64,
    tau: f64,
    boltzmann: f64,
    xi: f64,
}

impl NoseHoover {
    /// Creates the thermostat for target `t_ref` and period `tau`.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for a negative temperature or non-positive `tau`.
    pub fn new(t_ref: f64, tau: f64, units: UnitSystem) -> Result<Self, SimError> {
        validate_thermostat(t_ref, tau)?;
        Ok(Self { t_ref, tau, boltzmann: units.boltzmann(), xi: 0.0 })
    }

    /// Current friction variable `ξ` (1/ps).
    pub fn xi(&self) -> f64 {
        self.xi
    }
}

impl Thermostat for NoseHoover {
    fn name(&self) -> &'static str {
        "nose-hoover"
    }

    fn target(&self) -> f64 {
        self.t_ref
    }

    fn apply(&mut self, system: &mut System, kinetic: f64, dt: f64) {
        let Some(t_now) = temperature_from(system, kinetic, self.boltzmann) else {
            return;
        };
        let omega2 = (2.0 * std::f64::consts::PI / self.tau).powi(2);
        self.xi += dt * omega2 * (t_now / self.t_ref - 1.0);
        let scale = (-self.xi * dt).exp();
        for v in &mut system.velocities {
            *v *= scale;
        }
    }

    fn save_state(&self) -> Vec<u64> {
        vec![self.xi.to_bits()]
    }

    fn load_state(&mut self, state: &[u64]) -> Result<(), SimError> {
        match state {
            [x] => {
                self.xi = f64::from_bits(*x);
                Ok(())
            }
            _ => Err(SimError::Invalid("Nose-Hoover state must have 1 word".into())),
        }
    }
}

/// Berendsen weak-coupling thermostat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Berendsen {
    t_ref: f64,
    tau: f64,
    boltzmann: f64,
}

impl Berendsen {
    /// Creates the thermostat for reference temperature `t_ref` and coupling
    /// time `tau` (in the time unit of the run); `units` sets the Boltzmann constant.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] if `t_ref` is negative or `tau` is not positive.
    pub fn new(t_ref: f64, tau: f64, units: UnitSystem) -> Result<Self, SimError> {
        if !(t_ref.is_finite() && t_ref >= 0.0) {
            return Err(SimError::Config(format!(
                "thermostat temperature must be >= 0, got {t_ref}"
            )));
        }
        if !(tau.is_finite() && tau > 0.0) {
            return Err(SimError::Config(format!("thermostat tau must be > 0, got {tau}")));
        }
        Ok(Self { t_ref, tau, boltzmann: units.boltzmann() })
    }
}

impl Thermostat for Berendsen {
    fn name(&self) -> &'static str {
        "berendsen"
    }

    fn target(&self) -> f64 {
        self.t_ref
    }

    fn apply(&mut self, system: &mut System, kinetic: f64, dt: f64) {
        let dof = system.degrees_of_freedom();
        if dof <= 0.0 || kinetic <= 0.0 {
            return;
        }
        let t_now = 2.0 * kinetic / (dof * self.boltzmann);
        if t_now <= 0.0 {
            return;
        }
        let lambda2 = 1.0 + (dt / self.tau) * (self.t_ref / t_now - 1.0);
        let lambda = lambda2.max(0.0).sqrt();
        for v in &mut system.velocities {
            *v *= lambda;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrator::{ForceEvaluator, Integrator, LeapFrog, VelocityVerlet};

    #[test]
    fn rejects_bad_parameters() {
        assert!(Berendsen::new(-1.0, 0.5, UnitSystem::Md).is_err());
        assert!(Berendsen::new(1.0, 0.0, UnitSystem::Md).is_err());
        assert!(VRescale::new(1.0, 0.0, UnitSystem::Md, 1).is_err());
        assert!(NoseHoover::new(-1.0, 0.1, UnitSystem::Md).is_err());
    }

    /// Runs `steps` of a fixed-temperature simulation and returns the time-averaged
    /// temperature over the second half of the run.
    fn average_temperature(th: Box<dyn Thermostat>, steps: usize) -> f64 {
        use crate::forcefield::{LennardJones, Modifier};
        use crate::simulation::Simulation;
        let sys = crate::integrator::tests::lj_gas(6, 0.7, 5);
        let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap();
        let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
        let mut sys = sys;
        crate::simulation::maxwell_boltzmann_velocities(&mut sys, 0.2, 9, UnitSystem::Reduced);
        let mut sim = Simulation::new(
            sys,
            eval,
            Box::new(VelocityVerlet),
            Some(th),
            0.002,
            UnitSystem::Reduced,
        )
        .unwrap();
        let mut acc = 0.0;
        for k in 0..steps {
            sim.advance().unwrap();
            if k >= steps / 2 {
                acc += sim.energies().temperature;
            }
        }
        acc / (steps - steps / 2) as f64
    }

    #[test]
    fn v_rescale_samples_target_temperature() {
        let th = VRescale::new(0.6, 0.1, UnitSystem::Reduced, 17).unwrap();
        let t = average_temperature(Box::new(th), 3000);
        assert!((t - 0.6).abs() < 0.06, "v-rescale average T = {t}");
    }

    #[test]
    fn nose_hoover_drives_temperature_to_target() {
        let th = NoseHoover::new(0.6, 0.1, UnitSystem::Reduced).unwrap();
        let t = average_temperature(Box::new(th), 3000);
        assert!((t - 0.6).abs() < 0.08, "Nose-Hoover average T = {t}");
    }

    #[test]
    fn thermostat_state_round_trips() {
        let mut a = VRescale::new(1.0, 0.1, UnitSystem::Reduced, 3).unwrap();
        let mut sys = crate::integrator::tests::lj_gas(3, 0.5, 1);
        a.apply(&mut sys, 1.0, 0.001);
        let state = a.save_state();
        let mut b = VRescale::new(1.0, 0.1, UnitSystem::Reduced, 999).unwrap();
        b.load_state(&state).unwrap();
        assert_eq!(a.rng.gaussian().to_bits(), b.rng.gaussian().to_bits());

        let mut n = NoseHoover::new(1.0, 0.1, UnitSystem::Reduced).unwrap();
        n.xi = 0.37;
        let mut m = NoseHoover::new(1.0, 0.1, UnitSystem::Reduced).unwrap();
        m.load_state(&n.save_state()).unwrap();
        assert_eq!(m.xi(), 0.37);
    }

    #[test]
    fn drives_temperature_towards_target() {
        use crate::forcefield::{LennardJones, Modifier};
        for integ in [&mut VelocityVerlet as &mut dyn Integrator, &mut LeapFrog] {
            let mut sys = crate::integrator::tests::lj_gas(6, 0.6, 3);
            let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap();
            let mut eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
            eval.evaluate(&mut sys).unwrap();
            integ.prepare(&mut sys, 0.002);
            let mut th = Berendsen::new(0.5, 0.05, UnitSystem::Reduced).unwrap();
            for _ in 0..2000 {
                let ke = integ.kinetic_energy(&sys, 0.002);
                integ
                    .step(&mut sys, &mut eval, 0.002, &mut crate::constraints::NoConstraints)
                    .unwrap();
                th.apply(&mut sys, ke, 0.002);
            }
            let t = sys.temperature_in(UnitSystem::Reduced);
            assert!((t - 0.5).abs() < 0.1, "{}: T = {t}", integ.name());
        }
    }
}
