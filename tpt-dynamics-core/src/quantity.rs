//! Typed physical quantities.
//!
//! Each newtype stores its value in SI base units internally. Constructors and
//! accessors exist for the common units in each dimension, so that mixing e.g.
//! eV and kJ/mol is a compile-time-visible choice rather than a silent bug.
//!
//! The types deliberately implement only dimensionally safe arithmetic
//! (same-type addition and subtraction, scaling by a dimensionless `f64`).

use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::units;

macro_rules! quantity {
    ($(#[$meta:meta])* $name:ident, $si:literal) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug, Default, PartialEq, PartialOrd)]
        pub struct $name(f64);

        impl $name {
            /// Value in SI base units.
            #[doc = concat!("Unit: ", $si, ".")]
            #[inline]
            pub const fn si(self) -> f64 {
                self.0
            }

            /// Construct from a value already in SI base units.
            #[inline]
            pub const fn from_si(value: f64) -> Self {
                Self(value)
            }
        }

        impl Add for $name {
            type Output = Self;
            #[inline]
            fn add(self, rhs: Self) -> Self {
                Self(self.0 + rhs.0)
            }
        }

        impl Sub for $name {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: Self) -> Self {
                Self(self.0 - rhs.0)
            }
        }

        impl Neg for $name {
            type Output = Self;
            #[inline]
            fn neg(self) -> Self {
                Self(-self.0)
            }
        }

        impl Mul<f64> for $name {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: f64) -> Self {
                Self(self.0 * rhs)
            }
        }

        impl Div<f64> for $name {
            type Output = Self;
            #[inline]
            fn div(self, rhs: f64) -> Self {
                Self(self.0 / rhs)
            }
        }
    };
}

quantity!(
    /// An energy. Stored in joules.
    Energy,
    "J"
);
quantity!(
    /// A length. Stored in metres.
    Length,
    "m"
);
quantity!(
    /// A mass. Stored in kilograms.
    Mass,
    "kg"
);
quantity!(
    /// An absolute temperature. Stored in kelvin.
    Temperature,
    "K"
);

impl Energy {
    /// Construct from electron volts.
    #[inline]
    pub fn from_ev(ev: f64) -> Self {
        Self(units::ev_to_joule(ev))
    }
    /// Value in electron volts.
    #[inline]
    pub fn ev(self) -> f64 {
        units::joule_to_ev(self.0)
    }
    /// Construct from kJ/mol (molar energy).
    #[inline]
    pub fn from_kj_per_mol(kj: f64) -> Self {
        Self(kj * 1e3 / crate::constants::AVOGADRO)
    }
    /// Value in kJ/mol (molar energy).
    #[inline]
    pub fn kj_per_mol(self) -> f64 {
        self.0 * crate::constants::AVOGADRO * 1e-3
    }
    /// Construct from kcal/mol (molar energy).
    #[inline]
    pub fn from_kcal_per_mol(kcal: f64) -> Self {
        Self::from_kj_per_mol(units::kcal_per_mol_to_kj_per_mol(kcal))
    }
    /// Value in kcal/mol (molar energy).
    #[inline]
    pub fn kcal_per_mol(self) -> f64 {
        units::kj_per_mol_to_kcal_per_mol(self.kj_per_mol())
    }
}

impl Length {
    /// Construct from ångström.
    #[inline]
    pub fn from_angstrom(a: f64) -> Self {
        Self(units::angstrom_to_metre(a))
    }
    /// Value in ångström.
    #[inline]
    pub fn angstrom(self) -> f64 {
        units::metre_to_angstrom(self.0)
    }
    /// Construct from nanometres.
    #[inline]
    pub fn from_nm(nm: f64) -> Self {
        Self(nm * 1e-9)
    }
    /// Value in nanometres.
    #[inline]
    pub fn nm(self) -> f64 {
        self.0 * 1e9
    }
}

impl Mass {
    /// Construct from unified atomic mass units (Da).
    #[inline]
    pub fn from_amu(amu: f64) -> Self {
        Self(units::amu_to_kg(amu))
    }
    /// Value in unified atomic mass units.
    #[inline]
    pub fn amu(self) -> f64 {
        units::kg_to_amu(self.0)
    }
}

impl Temperature {
    /// Construct from kelvin (identity; provided for symmetry).
    #[inline]
    pub const fn from_kelvin(k: f64) -> Self {
        Self(k)
    }
    /// Value in kelvin.
    #[inline]
    pub const fn kelvin(self) -> f64 {
        self.0
    }
    /// Thermal energy `k_B·T`.
    #[inline]
    pub fn thermal_energy(self) -> Energy {
        Energy(units::kelvin_to_joule(self.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn energy_units_round_trip() {
        let e = Energy::from_ev(2.5);
        assert!((e.ev() - 2.5).abs() < 1e-14);
        let m = Energy::from_kcal_per_mol(1.0);
        assert!((m.kcal_per_mol() - 1.0).abs() < 1e-14);
        assert!((m.kj_per_mol() - 4.184).abs() < 1e-12);
    }

    #[test]
    fn length_units() {
        assert!((Length::from_angstrom(1.0).nm() - 0.1).abs() < 1e-15);
        assert!((Length::from_nm(0.5).angstrom() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn same_type_arithmetic() {
        let a = Energy::from_ev(1.0);
        let b = Energy::from_ev(0.25);
        assert!(((a - b).ev() - 0.75).abs() < 1e-14);
        assert!((((a + b) * 2.0).ev() - 2.5).abs() < 1e-14);
        assert!((((-a) / 2.0).ev() + 0.5).abs() < 1e-14);
    }

    #[test]
    fn temperature_thermal_energy() {
        let t = Temperature::from_kelvin(300.0);
        assert!((t.thermal_energy().ev() - 0.025_852).abs() < 1e-5);
    }
}
