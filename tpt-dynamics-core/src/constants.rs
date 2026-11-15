//! Physical constants in SI units.
//!
//! Source: CODATA 2018 recommended values (the 2019 SI redefinition made several
//! of these exact). Constants marked *exact* are defined by the SI and have no
//! measurement uncertainty. Each constant is annotated with its origin.
//!
//! Values are stored at full published precision. Derived constants are computed
//! at compile time from the primary values where that is exact or standard.

/// Speed of light in vacuum, m s⁻¹. *Exact* (SI definition).
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;

/// Planck constant, J s. *Exact* (SI definition).
pub const PLANCK: f64 = 6.626_070_15e-34;

/// Reduced Planck constant ħ = h / 2π, J s.
pub const HBAR: f64 = PLANCK / (2.0 * core::f64::consts::PI);

/// Boltzmann constant, J K⁻¹. *Exact* (SI definition).
pub const BOLTZMANN: f64 = 1.380_649e-23;

/// Avogadro constant, mol⁻¹. *Exact* (SI definition).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Molar gas constant R = `N_A` · `k_B`, J mol⁻¹ K⁻¹ (derived).
pub const GAS_CONSTANT: f64 = AVOGADRO * BOLTZMANN;

/// Elementary charge, C. *Exact* (SI definition).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// Electron volt, J. *Exact*: 1 eV = e × 1 V.
pub const ELECTRON_VOLT: f64 = ELEMENTARY_CHARGE;

/// Vacuum permittivity ε₀, F m⁻¹. CODATA 2018 (measured).
pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_812_8e-12;

/// Vacuum permeability μ₀, N A⁻². CODATA 2018 (measured).
pub const VACUUM_PERMEABILITY: f64 = 1.256_637_062_12e-6;

/// Coulomb constant `k_e` = 1 / (4π `ε₀`), N m² C⁻² (derived).
pub const COULOMB_CONSTANT: f64 = 1.0 / (4.0 * core::f64::consts::PI * VACUUM_PERMITTIVITY);

/// Atomic mass unit (unified), kg. CODATA 2018 (measured).
pub const ATOMIC_MASS_UNIT: f64 = 1.660_539_066_60e-27;

/// Electron rest mass, kg. CODATA 2018 (measured).
pub const ELECTRON_MASS: f64 = 9.109_383_701_5e-31;

/// Proton rest mass, kg. CODATA 2018 (measured).
pub const PROTON_MASS: f64 = 1.672_621_923_69e-27;

/// Hartree energy `E_h`, J. CODATA 2018 (measured).
pub const HARTREE_ENERGY: f64 = 4.359_744_722_207_1e-18;

/// Bohr radius a₀, m. CODATA 2018 (measured).
pub const BOHR_RADIUS: f64 = 5.291_772_109_03e-11;

/// Thermochemical calorie, J. *Exact* by definition (International Table
/// calorie is not used; this is the thermochemical convention used by GROMACS).
pub const THERMOCHEMICAL_CALORIE: f64 = 4.184;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_si_values() {
        assert_eq!(SPEED_OF_LIGHT, 299_792_458.0);
        assert_eq!(PLANCK, 6.626_070_15e-34);
        assert_eq!(BOLTZMANN, 1.380_649e-23);
        assert_eq!(AVOGADRO, 6.022_140_76e23);
        assert_eq!(ELEMENTARY_CHARGE, 1.602_176_634e-19);
    }

    #[test]
    fn derived_constants_consistent() {
        // ħ = h / 2π
        assert!((HBAR * 2.0 * core::f64::consts::PI - PLANCK).abs() / PLANCK < 1e-15);
        // R = N_A k_B ≈ 8.314 462 618 J mol⁻¹ K⁻¹ (CODATA 2018)
        assert!((GAS_CONSTANT - 8.314_462_618).abs() < 1e-8);
        // k_e ≈ 8.987 551 792 e9 N m² C⁻²
        assert!((COULOMB_CONSTANT - 8.987_551_792e9).abs() / 8.987_551_792e9 < 1e-9);
    }

    #[test]
    fn electron_volt_is_elementary_charge() {
        assert_eq!(ELECTRON_VOLT, ELEMENTARY_CHARGE);
    }
}
