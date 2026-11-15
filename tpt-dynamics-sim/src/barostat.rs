//! Pressure coupling. Both barostats here are isotropic: they scale all
//! coordinates and box vectors by one factor `μ` per step.
//!
//! * [`Berendsen`] (Berendsen et al., J. Chem. Phys. 81, 3684, 1984):
//!   `μ = [1 − (Δt/τ_p) β (P₀ − P)]^{1/3}`. Weak coupling, no fluctuations.
//! * [`CRescale`] (Bernetti & Bussi, J. Chem. Phys. 153, 114107, 2020), isotropic
//!   form: the Berendsen factor plus a stochastic term with variance
//!   `2 k_B T β Δt / (V τ_p)`, so the volume fluctuates with the correct
//!   isothermal-isobaric distribution.
//! * [`ParrinelloRahman`] (Parrinello & Rahman 1981), scalar (isotropic) form:
//!   the box strain obeys a second-order equation of motion driven by
//!   `3V(P − P₀)`. The full tensor form is not implemented.
//!
//! Units: `pressure` and `compressibility` must be in the engine's pressure
//! unit (kJ mol⁻¹ nm⁻³ for MD units, reduced units otherwise). Configuration
//! files take bar and bar⁻¹ and convert with [`BAR_PER_MD_PRESSURE`].

use crate::error::SimError;
use crate::rng::Rng;

/// Value of one kJ mol⁻¹ nm⁻³ expressed in bar.
pub const BAR_PER_MD_PRESSURE: f64 = 16.6054;

/// Isotropic pressure coupling.
pub trait Barostat: Send {
    /// Short identifier used in logs.
    fn name(&self) -> &'static str;

    /// Target pressure in the engine's pressure unit.
    fn target(&self) -> f64;

    /// Returns the isotropic scale factor `μ` for this step, given the
    /// instantaneous pressure, volume, `k_B T` (engine energy units) and
    /// time step.
    fn scale_factor(&mut self, pressure: f64, volume: f64, kt: f64, dt: f64) -> f64;

    /// Internal state as raw 64-bit words, for checkpoints.
    fn save_state(&self) -> Vec<u64> {
        Vec::new()
    }

    /// Restores state produced by [`Barostat::save_state`].
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if the state does not match.
    fn load_state(&mut self, state: &[u64]) -> Result<(), SimError> {
        if state.is_empty() {
            Ok(())
        } else {
            Err(SimError::Invalid("barostat has no state to restore".into()))
        }
    }
}

fn validate(tau_p: f64, compressibility: f64) -> Result<(), SimError> {
    if !(tau_p.is_finite() && tau_p > 0.0) {
        return Err(SimError::Config(format!("barostat tau must be > 0, got {tau_p}")));
    }
    if !(compressibility.is_finite() && compressibility > 0.0) {
        return Err(SimError::Config(format!(
            "barostat compressibility must be > 0, got {compressibility}"
        )));
    }
    Ok(())
}

/// Turns the bracketed factor into a cube root, clamped so that a single step
/// cannot shrink or expand the box by more than 1% in length.
fn cube_root_clamped(x: f64) -> f64 {
    x.max(0.97_f64.powi(3)).min(1.03_f64.powi(3)).cbrt()
}

/// Berendsen isotropic pressure coupling.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Berendsen {
    p0: f64,
    tau_p: f64,
    compressibility: f64,
}

impl Berendsen {
    /// Creates the barostat for target `p0`, coupling time `tau_p` and
    /// isothermal compressibility `compressibility`.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for a non-positive `tau_p` or compressibility.
    pub fn new(p0: f64, tau_p: f64, compressibility: f64) -> Result<Self, SimError> {
        validate(tau_p, compressibility)?;
        Ok(Self { p0, tau_p, compressibility })
    }
}

impl Barostat for Berendsen {
    fn name(&self) -> &'static str {
        "berendsen"
    }

    fn target(&self) -> f64 {
        self.p0
    }

    fn scale_factor(&mut self, pressure: f64, _volume: f64, _temperature: f64, dt: f64) -> f64 {
        cube_root_clamped(1.0 - dt / self.tau_p * self.compressibility * (self.p0 - pressure))
    }
}

/// Isotropic C-rescale barostat (stochastic cell rescaling).
#[derive(Clone, Debug)]
pub struct CRescale {
    p0: f64,
    tau_p: f64,
    compressibility: f64,
    rng: Rng,
}

impl CRescale {
    /// Creates the barostat. `seed` drives the stochastic term; the generator
    /// state is part of the checkpoint.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for a non-positive `tau_p` or compressibility.
    pub fn new(p0: f64, tau_p: f64, compressibility: f64, seed: u64) -> Result<Self, SimError> {
        validate(tau_p, compressibility)?;
        Ok(Self { p0, tau_p, compressibility, rng: Rng::seed_from_u64(seed) })
    }
}

impl Barostat for CRescale {
    fn name(&self) -> &'static str {
        "c-rescale"
    }

    fn target(&self) -> f64 {
        self.p0
    }

    fn scale_factor(&mut self, pressure: f64, volume: f64, temperature: f64, dt: f64) -> f64 {
        // `temperature` is in the engine's units, so k_B T is `temperature` times
        // the engine's Boltzmann constant. Callers pass k_B T directly.
        let noise = self.rng.gaussian();
        let beta = self.compressibility;
        let drift = -dt / self.tau_p * beta * (self.p0 - pressure);
        let sigma2 = 2.0 * temperature * beta * dt / (volume * self.tau_p);
        cube_root_clamped(1.0 + drift + sigma2.max(0.0).sqrt() * noise)
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
            _ => Err(SimError::Invalid("c-rescale state must have 4 words".into())),
        }
    }
}

/// Isotropic Parrinello–Rahman barostat (Parrinello & Rahman, J. Appl. Phys.
/// 52, 7182, 1981), scalar form.
///
/// The box is described by the strain `ε = ln L` of its edge length, with a
/// fictitious mass `W = 9 V τ_p² / (4π² β)` so that small oscillations of the
/// box have period `τ_p` for compressibility `β`. The equations of motion are
///
/// ```text
/// ε̇ = p_ε / W,    ṗ_ε = 3 V (P − P₀)
/// ```
///
/// integrated with a semi-implicit Euler step per MD step. The scale factor
/// applied to coordinates and box vectors is `μ = exp(Δt p_ε / W)`.
///
/// This is the isotropic (volume-only) variant. The full tensor form, which
/// lets the three box vectors respond separately to the pressure tensor, is
/// not implemented. The pressure is the scalar pressure from the virial, and
/// constraint virials are not included (see [`crate::constraints`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ParrinelloRahman {
    p0: f64,
    tau_p: f64,
    compressibility: f64,
    /// Momentum conjugate to the strain, `p_ε = W ε̇`.
    p_eps: f64,
}

impl ParrinelloRahman {
    /// Creates the barostat for target `p0`, period `tau_p` and compressibility
    /// `compressibility` (engine pressure units).
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for a non-positive `tau_p` or compressibility.
    pub fn new(p0: f64, tau_p: f64, compressibility: f64) -> Result<Self, SimError> {
        validate(tau_p, compressibility)?;
        Ok(Self { p0, tau_p, compressibility, p_eps: 0.0 })
    }

    /// Current box-strain momentum `p_ε`, in engine units times ps.
    pub fn strain_momentum(&self) -> f64 {
        self.p_eps
    }
}

impl Barostat for ParrinelloRahman {
    fn name(&self) -> &'static str {
        "parrinello-rahman"
    }

    fn target(&self) -> f64 {
        self.p0
    }

    fn scale_factor(&mut self, pressure: f64, volume: f64, _kt: f64, dt: f64) -> f64 {
        let w = 9.0 * volume * self.tau_p * self.tau_p
            / (4.0 * std::f64::consts::PI * std::f64::consts::PI * self.compressibility);
        self.p_eps += dt * 3.0 * volume * (pressure - self.p0);
        let mu = (dt * self.p_eps / w).exp();
        mu.clamp(0.97, 1.03)
    }

    fn save_state(&self) -> Vec<u64> {
        vec![self.p_eps.to_bits()]
    }

    fn load_state(&mut self, state: &[u64]) -> Result<(), SimError> {
        match state {
            [p] => {
                self.p_eps = f64::from_bits(*p);
                Ok(())
            }
            _ => Err(SimError::Invalid("Parrinello-Rahman state must have 1 word".into())),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parrinello_rahman_box_follows_the_pressure_difference() {
        let mut high = ParrinelloRahman::new(1.0, 1.0, 0.5).unwrap();
        let mut low = ParrinelloRahman::new(1.0, 1.0, 0.5).unwrap();
        let (mut mu_high, mut mu_low) = (1.0, 1.0);
        for _ in 0..50 {
            mu_high = high.scale_factor(3.0, 1.0, 0.0, 0.01);
            mu_low = low.scale_factor(-1.0, 1.0, 0.0, 0.01);
        }
        assert!(mu_high > 1.0 && high.strain_momentum() > 0.0, "high pressure must expand");
        assert!(mu_low < 1.0 && low.strain_momentum() < 0.0, "low pressure must shrink");
        // Balanced pressure leaves the box alone.
        let mut bal = ParrinelloRahman::new(1.0, 1.0, 0.5).unwrap();
        assert_eq!(bal.scale_factor(1.0, 1.0, 0.0, 0.01), 1.0);
    }

    #[test]
    fn parrinello_rahman_state_round_trips() {
        let mut a = ParrinelloRahman::new(1.0, 1.0, 0.5).unwrap();
        let _ = a.scale_factor(2.0, 1.0, 0.0, 0.01);
        let mut b = ParrinelloRahman::new(1.0, 1.0, 0.5).unwrap();
        b.load_state(&a.save_state()).unwrap();
        assert_eq!(a.strain_momentum().to_bits(), b.strain_momentum().to_bits());
        assert!(b.load_state(&[1, 2]).is_err());
    }

    #[test]
    fn berendsen_expands_when_pressure_is_too_high_and_shrinks_when_low() {
        let mut b = Berendsen::new(1.0, 1.0, 0.5).unwrap();
        // P above the target: the box expands.
        assert!(b.scale_factor(3.0, 1.0, 1.0, 0.01) > 1.0);
        // P below the target: the box shrinks.
        assert!(b.scale_factor(-1.0, 1.0, 1.0, 0.01) < 1.0);
        assert!((b.scale_factor(1.0, 1.0, 1.0, 0.01) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn c_rescale_noise_is_zero_for_zero_temperature() {
        let mut c = CRescale::new(1.0, 1.0, 0.5, 4).unwrap();
        // With T = 0 the stochastic term vanishes and the step equals Berendsen's.
        let mu = c.scale_factor(3.0, 1.0, 0.0, 0.01);
        let mut b = Berendsen::new(1.0, 1.0, 0.5).unwrap();
        assert!((mu - b.scale_factor(3.0, 1.0, 0.0, 0.01)).abs() < 1e-15);
    }

    #[test]
    fn rejects_invalid_parameters() {
        assert!(Berendsen::new(1.0, 0.0, 1.0).is_err());
        assert!(CRescale::new(1.0, 1.0, -1.0, 0).is_err());
    }

    #[test]
    fn c_rescale_state_round_trips() {
        let mut a = CRescale::new(1.0, 1.0, 0.5, 8).unwrap();
        let _ = a.scale_factor(1.0, 1.0, 1.0, 0.01);
        let mut b = CRescale::new(1.0, 1.0, 0.5, 1).unwrap();
        b.load_state(&a.save_state()).unwrap();
        assert_eq!(
            a.scale_factor(2.0, 1.0, 1.0, 0.01).to_bits(),
            b.scale_factor(2.0, 1.0, 1.0, 0.01).to_bits()
        );
    }
}
