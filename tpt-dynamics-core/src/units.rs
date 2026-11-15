//! Scalar unit conversions.
//!
//! Every function takes and returns plain `f64`. For type-checked conversions
//! use the newtypes in [`crate::quantity`].
//!
//! Conventions: energies per mole use the *thermochemical* kcal (4.184 kJ), as
//! used by GROMACS and most MD force-field literature.

use crate::constants::{
    ATOMIC_MASS_UNIT, AVOGADRO, BOLTZMANN, ELECTRON_VOLT, THERMOCHEMICAL_CALORIE,
};

/// Joules per electron volt.
#[inline]
pub fn ev_to_joule(ev: f64) -> f64 {
    ev * ELECTRON_VOLT
}

/// Electron volts per joule.
#[inline]
pub fn joule_to_ev(joule: f64) -> f64 {
    joule / ELECTRON_VOLT
}

/// Ångström to nanometre (1 nm = 10 Å, exact).
#[inline]
pub fn angstrom_to_nm(angstrom: f64) -> f64 {
    angstrom * 0.1
}

/// Nanometre to Ångström.
#[inline]
pub fn nm_to_angstrom(nm: f64) -> f64 {
    nm * 10.0
}

/// Ångström to metre.
#[inline]
pub fn angstrom_to_metre(angstrom: f64) -> f64 {
    angstrom * 1e-10
}

/// Metre to Ångström.
#[inline]
pub fn metre_to_angstrom(metre: f64) -> f64 {
    metre * 1e10
}

/// kJ/mol to kcal/mol.
#[inline]
pub fn kj_per_mol_to_kcal_per_mol(kj: f64) -> f64 {
    kj / THERMOCHEMICAL_CALORIE
}

/// kcal/mol to kJ/mol.
#[inline]
pub fn kcal_per_mol_to_kj_per_mol(kcal: f64) -> f64 {
    kcal * THERMOCHEMICAL_CALORIE
}

/// Energy per particle in eV to molar energy in kJ/mol.
#[inline]
pub fn ev_to_kj_per_mol(ev: f64) -> f64 {
    ev * ELECTRON_VOLT * AVOGADRO * 1e-3
}

/// Molar energy in kJ/mol to energy per particle in eV.
#[inline]
pub fn kj_per_mol_to_ev(kj_per_mol: f64) -> f64 {
    kj_per_mol * 1e3 / (AVOGADRO * ELECTRON_VOLT)
}

/// Unified atomic mass units to kilograms.
#[inline]
pub fn amu_to_kg(amu: f64) -> f64 {
    amu * ATOMIC_MASS_UNIT
}

/// Kilograms to unified atomic mass units.
#[inline]
pub fn kg_to_amu(kg: f64) -> f64 {
    kg / ATOMIC_MASS_UNIT
}

/// Temperature (K) to thermal energy `k_B·T` (J).
#[inline]
pub fn kelvin_to_joule(kelvin: f64) -> f64 {
    kelvin * BOLTZMANN
}

/// Thermal energy (J) to temperature (K).
#[inline]
pub fn joule_to_kelvin(joule: f64) -> f64 {
    joule / BOLTZMANN
}

/// Temperature (K) to molar thermal energy R·T (kJ/mol).
#[inline]
pub fn kelvin_to_kj_per_mol(kelvin: f64) -> f64 {
    kelvin * BOLTZMANN * AVOGADRO * 1e-3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{ELEMENTARY_CHARGE, GAS_CONSTANT};

    fn close(a: f64, b: f64, rel: f64) -> bool {
        let scale = a.abs().max(b.abs()).max(f64::MIN_POSITIVE);
        (a - b).abs() <= rel * scale
    }

    #[test]
    fn ev_joule_known_value() {
        assert!(close(ev_to_joule(1.0), 1.602_176_634e-19, 1e-15));
        assert!(close(joule_to_ev(ev_to_joule(3.7)), 3.7, 1e-14));
    }

    #[test]
    fn length_conversions_are_exact() {
        assert_eq!(angstrom_to_nm(10.0), 1.0);
        assert_eq!(nm_to_angstrom(1.0), 10.0);
        assert!(close(metre_to_angstrom(angstrom_to_metre(2.5)), 2.5, 1e-15));
    }

    #[test]
    fn kcal_kj_ratio() {
        assert!(close(kcal_per_mol_to_kj_per_mol(1.0), 4.184, 1e-15));
        assert!(close(kj_per_mol_to_kcal_per_mol(4.184), 1.0, 1e-15));
    }

    #[test]
    fn ev_per_particle_to_kj_per_mol() {
        // 1 eV × N_A ≈ 96.485 332 12 kJ/mol
        assert!(close(ev_to_kj_per_mol(1.0), 96.485_332_12, 1e-8));
        assert!(close(kj_per_mol_to_ev(ev_to_kj_per_mol(0.123)), 0.123, 1e-14));
    }

    #[test]
    fn mass_round_trip() {
        assert!(close(kg_to_amu(amu_to_kg(12.011)), 12.011, 1e-15));
    }

    #[test]
    fn thermal_energy() {
        // k_B · 300 K ≈ 4.14 × 10⁻²¹ J
        assert!(close(kelvin_to_joule(300.0), 4.141_947e-21, 1e-6));
        // R · T in kJ/mol agrees with R (J/mol/K) × T / 1000
        assert!(close(kelvin_to_kj_per_mol(300.0), GAS_CONSTANT * 300.0 / 1000.0, 1e-12));
        assert!(close(joule_to_kelvin(kelvin_to_joule(77.0)), 77.0, 1e-14));
    }

    #[test]
    fn electron_volt_matches_charge() {
        assert_eq!(ev_to_joule(1.0), ELEMENTARY_CHARGE);
    }
}

#[cfg(test)]
mod property_tests {
    //! Property-style checks: round-trips over a deterministic sweep of inputs
    //! (a SplitMix64 stream), so failures are reproducible without a PRNG crate.
    use super::*;

    struct SplitMix64(u64);

    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        /// Uniform value in [lo, hi).
        fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
            let u = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
            lo + (hi - lo) * u
        }
    }

    fn rel_close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(f64::MIN_POSITIVE)
    }

    #[test]
    fn round_trips_over_sweep() {
        let mut rng = SplitMix64(0x5EED);
        for _ in 0..10_000 {
            let x = rng.uniform(1e-6, 1e6);
            assert!(rel_close(joule_to_ev(ev_to_joule(x)), x));
            assert!(rel_close(kg_to_amu(amu_to_kg(x)), x));
            assert!(rel_close(joule_to_kelvin(kelvin_to_joule(x)), x));
            assert!(rel_close(nm_to_angstrom(angstrom_to_nm(x)), x));
            assert!(rel_close(metre_to_angstrom(angstrom_to_metre(x)), x));
            assert!(rel_close(kj_per_mol_to_kcal_per_mol(kcal_per_mol_to_kj_per_mol(x)), x));
            assert!(rel_close(kj_per_mol_to_ev(ev_to_kj_per_mol(x)), x));
        }
    }

    #[test]
    fn conversions_are_monotonic() {
        let mut rng = SplitMix64(42);
        for _ in 0..10_000 {
            let a = rng.uniform(0.0, 100.0);
            let b = a + rng.uniform(1e-9, 10.0);
            assert!(ev_to_joule(a) < ev_to_joule(b));
            assert!(kelvin_to_joule(a) < kelvin_to_joule(b));
            assert!(ev_to_kj_per_mol(a) < ev_to_kj_per_mol(b));
        }
    }
}
