//! Property tests for unit conversions: every forward/inverse pair must round-trip
//! to within floating-point rounding, over a wide range of magnitudes and signs.
//!
//! Inputs are generated from a fixed-seed SplitMix64 so failures are reproducible.

use tpt_dynamics_core::quantity::{Energy, Length, Mass, Temperature};
use tpt_dynamics_core::units::*;

const CASES: usize = 2_000;
const REL_TOL: f64 = 1e-12;

/// Deterministic generator over log-uniform magnitudes in `[1e-6, 1e6]` with random sign.
struct Gen(u64);

impl Gen {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    fn value(&mut self) -> f64 {
        let exp = -6.0 + 12.0 * self.unit();
        let mag = 10f64.powf(exp);
        if self.next_u64() & 1 == 0 {
            mag
        } else {
            -mag
        }
    }
}

fn close(a: f64, b: f64) -> bool {
    let scale = a.abs().max(b.abs());
    if scale == 0.0 {
        return a == b;
    }
    (a - b).abs() <= REL_TOL * scale
}

fn check_roundtrip(name: &str, forward: fn(f64) -> f64, inverse: fn(f64) -> f64) {
    let mut g = Gen(0xC0FF_EE00 ^ name.len() as u64);
    for _ in 0..CASES {
        let x = g.value();
        let back = inverse(forward(x));
        assert!(close(back, x), "{name}: x={x:e} round-tripped to {back:e}");
    }
}

#[test]
fn ev_joule_roundtrip() {
    check_roundtrip("ev<->joule", ev_to_joule, joule_to_ev);
}

#[test]
fn angstrom_nm_roundtrip() {
    check_roundtrip("angstrom<->nm", angstrom_to_nm, nm_to_angstrom);
}

#[test]
fn angstrom_metre_roundtrip() {
    check_roundtrip("angstrom<->metre", angstrom_to_metre, metre_to_angstrom);
}

#[test]
fn kj_kcal_roundtrip() {
    check_roundtrip("kj<->kcal", kj_per_mol_to_kcal_per_mol, kcal_per_mol_to_kj_per_mol);
}

#[test]
fn ev_kj_per_mol_roundtrip() {
    check_roundtrip("ev<->kj/mol", ev_to_kj_per_mol, kj_per_mol_to_ev);
}

#[test]
fn amu_kg_roundtrip() {
    check_roundtrip("amu<->kg", amu_to_kg, kg_to_amu);
}

#[test]
fn kelvin_joule_roundtrip() {
    check_roundtrip("kelvin<->joule", kelvin_to_joule, joule_to_kelvin);
}

#[test]
fn quantity_newtypes_roundtrip() {
    let mut g = Gen(0x1234_5678);
    for _ in 0..CASES {
        let x = g.value();
        assert!(close(Energy::from_ev(x).ev(), x), "ev {x:e}");
        assert!(close(Energy::from_kj_per_mol(x).kj_per_mol(), x), "kJ/mol {x:e}");
        assert!(close(Energy::from_kcal_per_mol(x).kcal_per_mol(), x), "kcal/mol {x:e}");
        assert!(close(Length::from_angstrom(x).angstrom(), x), "angstrom {x:e}");
        assert!(close(Length::from_nm(x).nm(), x), "nm {x:e}");
        assert!(close(Mass::from_amu(x).amu(), x), "amu {x:e}");
        assert!(close(Temperature::from_kelvin(x).kelvin(), x), "kelvin {x:e}");
    }
}

#[test]
fn angstrom_and_nm_quantities_agree_with_unit_functions() {
    let mut g = Gen(0xABCD);
    for _ in 0..CASES {
        let a = g.value();
        assert!(close(Length::from_angstrom(a).nm(), angstrom_to_nm(a)), "a={a:e}");
    }
}
