//! Kinematic helpers for collider analysis.
//!
//! Particles are described in detector coordinates `(pt, eta, phi, m)`:
//! transverse momentum, pseudorapidity, azimuth in radians, and mass, in
//! consistent units (for example `GeV`). Use [`FourVector::from_pt_eta_phi_m`] to
//! convert to Cartesian components, then combine them with [`invariant_mass`].

use core::f64::consts::PI;

/// Energy-momentum four-vector in Cartesian components.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FourVector {
    /// Momentum x component.
    pub px: f64,
    /// Momentum y component.
    pub py: f64,
    /// Momentum z component.
    pub pz: f64,
    /// Energy.
    pub e: f64,
}

impl FourVector {
    /// Build from detector coordinates `(pt, eta, phi, m)`.
    ///
    /// `px = pt cos φ`, `py = pt sin φ`, `pz = pt sinh η`, `E = √(p² + m²)`.
    /// A negative mass is treated as its absolute value.
    pub fn from_pt_eta_phi_m(pt: f64, eta: f64, phi: f64, m: f64) -> Self {
        let pz = pt * eta.sinh();
        let p2 = pt * pt + pz * pz;
        Self { px: pt * phi.cos(), py: pt * phi.sin(), pz, e: (p2 + m * m).sqrt() }
    }

    /// Sum of two four-vectors.
    #[inline]
    #[must_use]
    pub fn plus(self, other: Self) -> Self {
        Self {
            px: self.px + other.px,
            py: self.py + other.py,
            pz: self.pz + other.pz,
            e: self.e + other.e,
        }
    }

    /// Transverse momentum `√(px² + py²)`.
    #[inline]
    pub fn pt(&self) -> f64 {
        self.px.hypot(self.py)
    }

    /// Invariant mass `√(E² − p²)`, clamped to zero where rounding makes the
    /// argument slightly negative.
    #[inline]
    pub fn mass(&self) -> f64 {
        let m2 = self.e * self.e - self.px * self.px - self.py * self.py - self.pz * self.pz;
        m2.max(0.0).sqrt()
    }

    /// Pseudorapidity `η = asinh(pz / pt)`. Returns `±∞` for momentum along the
    /// beam axis, with the sign of `pz`, and `0` for zero momentum.
    pub fn eta(&self) -> f64 {
        let pt = self.pt();
        if pt == 0.0 {
            if self.pz == 0.0 {
                0.0
            } else {
                f64::INFINITY.copysign(self.pz)
            }
        } else {
            (self.pz / pt).asinh()
        }
    }

    /// Azimuthal angle in `(-π, π]`.
    #[inline]
    pub fn phi(&self) -> f64 {
        self.py.atan2(self.px)
    }
}

/// Invariant mass of a set of particles given as four-vectors.
///
/// An empty slice gives `0.0`.
pub fn invariant_mass(particles: &[FourVector]) -> f64 {
    particles.iter().copied().reduce(FourVector::plus).map_or(0.0, |sum| sum.mass())
}

/// Invariant mass of two particles given in detector coordinates
/// `(pt, eta, phi, m)`.
pub fn pair_mass(a: (f64, f64, f64, f64), b: (f64, f64, f64, f64)) -> f64 {
    let va = FourVector::from_pt_eta_phi_m(a.0, a.1, a.2, a.3);
    let vb = FourVector::from_pt_eta_phi_m(b.0, b.1, b.2, b.3);
    va.plus(vb).mass()
}

/// Transverse momentum of a particle with components `(px, py)`.
#[inline]
pub fn transverse_momentum(px: f64, py: f64) -> f64 {
    px.hypot(py)
}

/// Azimuthal separation `φ₁ − φ₂` wrapped into `(-π, π]`.
#[inline]
pub fn delta_phi(phi1: f64, phi2: f64) -> f64 {
    wrap_pi(phi1 - phi2)
}

/// Angular separation `ΔR = √(Δη² + Δφ²)` in the (η, φ) plane, with `Δφ` wrapped
/// into `(-π, π]`.
#[inline]
pub fn delta_r(eta1: f64, phi1: f64, eta2: f64, phi2: f64) -> f64 {
    let de = eta1 - eta2;
    let dp = delta_phi(phi1, phi2);
    de.hypot(dp)
}

/// Wrap an angle into `(-π, π]`.
#[inline]
fn wrap_pi(x: f64) -> f64 {
    if !x.is_finite() {
        return f64::NAN;
    }
    let two_pi = 2.0 * PI;
    let mut y = x % two_pi;
    if y <= -PI {
        y += two_pi;
    } else if y > PI {
        y -= two_pi;
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn four_vector_round_trip_of_coordinates() {
        let v = FourVector::from_pt_eta_phi_m(30.0, 1.2, -0.7, 0.105);
        assert!(close(v.pt(), 30.0, 1e-12));
        assert!(close(v.eta(), 1.2, 1e-12));
        assert!(close(v.phi(), -0.7, 1e-12));
        assert!(close(v.mass(), 0.105, 1e-9));
    }

    #[test]
    fn invariant_mass_of_back_to_back_pair() {
        // Two massless particles at η = 0, Δφ = π, pt = 45 each: m = 2 · 45 = 90.
        let a = FourVector::from_pt_eta_phi_m(45.0, 0.0, 0.0, 0.0);
        let b = FourVector::from_pt_eta_phi_m(45.0, 0.0, PI, 0.0);
        assert!(close(invariant_mass(&[a, b]), 90.0, 1e-9));
        assert!(close(pair_mass((45.0, 0.0, 0.0, 0.0), (45.0, 0.0, PI, 0.0)), 90.0, 1e-9));
    }

    #[test]
    fn invariant_mass_of_collinear_massless_pair_is_zero() {
        let a = FourVector::from_pt_eta_phi_m(10.0, 0.5, 1.0, 0.0);
        let b = FourVector::from_pt_eta_phi_m(20.0, 0.5, 1.0, 0.0);
        // Rounding can make m² slightly negative; the result must still be ~0, not NaN.
        let m = invariant_mass(&[a, b]);
        assert!(m.is_finite() && m < 1e-6, "{m}");
    }

    #[test]
    fn empty_and_single_particle() {
        assert_eq!(invariant_mass(&[]), 0.0);
        let v = FourVector::from_pt_eta_phi_m(5.0, -2.0, 0.3, 4.0);
        assert!(close(invariant_mass(&[v]), 4.0, 1e-9));
    }

    #[test]
    fn mass_is_invariant_under_boost_along_z() {
        let a = FourVector::from_pt_eta_phi_m(12.0, 0.3, 0.2, 1.0);
        let b = FourVector::from_pt_eta_phi_m(7.0, -1.1, 2.5, 0.5);
        let m0 = invariant_mass(&[a, b]);
        // Boost by rapidity y: pz' = γ(pz + βE), E' = γ(E + β pz).
        let beta: f64 = 0.6;
        let gamma = 1.0 / (1.0 - beta * beta).sqrt();
        let boost = |v: FourVector| FourVector {
            pz: gamma * (v.pz + beta * v.e),
            e: gamma * (v.e + beta * v.pz),
            ..v
        };
        let m1 = invariant_mass(&[boost(a), boost(b)]);
        assert!(close(m0, m1, 1e-9));
    }

    #[test]
    fn transverse_momentum_helper() {
        assert!(close(transverse_momentum(3.0, 4.0), 5.0, 1e-12));
    }

    #[test]
    fn delta_phi_wraps() {
        assert!(close(delta_phi(3.0, -3.0), 6.0 - 2.0 * PI, 1e-12));
        // PI and -PI are the same angle; the result is 0 after wrapping.
        assert!(close(delta_phi(PI, -PI), 0.0, 1e-12));
        for k in -20..=20 {
            let d = delta_phi(f64::from(k) * 0.37, 0.0);
            assert!(d > -PI - 1e-12 && d <= PI + 1e-12, "{d}");
        }
        assert!(delta_phi(f64::NAN, 0.0).is_nan());
    }

    #[test]
    fn delta_r_basic() {
        // Δφ = -4 wraps to 2π - 4, so ΔR = sqrt(3² + (2π - 4)²).
        let expected = (9.0 + (2.0 * PI - 4.0).powi(2)).sqrt();
        assert!(close(delta_r(0.0, 0.0, 3.0, 4.0), expected, 1e-12));
        // Across the φ seam: Δφ = 0.2, not 2π − 0.2.
        assert!(close(delta_r(0.0, PI - 0.1, 0.0, -PI + 0.1), 0.2, 1e-12));
    }

    #[test]
    fn eta_on_beam_axis() {
        let fwd = FourVector { px: 0.0, py: 0.0, pz: 5.0, e: 5.0 };
        assert_eq!(fwd.eta(), f64::INFINITY);
        let bwd = FourVector { px: 0.0, py: 0.0, pz: -5.0, e: 5.0 };
        assert_eq!(bwd.eta(), f64::NEG_INFINITY);
        assert_eq!(FourVector::default().eta(), 0.0);
    }
}
