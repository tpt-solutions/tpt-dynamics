# tpt-dynamics-core

Foundation crate for [tpt-dynamics](../README.md): CODATA 2018 physical constants,
scalar unit conversions, typed physical quantities and a 3-vector type.

Everything here is shared by the rest of the workspace (`tpt-dynamics-data`,
`tpt-dynamics-sim`, `tpt-dynamics-compute`, the CLI and the Python bindings), so the
numbers are defined once.

- **`no_std`**: the crate does not link the standard library unless the `std` feature
  is enabled. It does not need `alloc`.
- **No `unsafe` code**, enforced with `#![forbid(unsafe_code)]`.
- **No FFI**. The only dependency is [`libm`](https://crates.io/crates/libm), which
  supplies the transcendental functions for targets without a math library.

## Contents

| Module | What it provides |
| --- | --- |
| `constants` | CODATA 2018 values: `SPEED_OF_LIGHT`, `PLANCK`, `HBAR`, `BOLTZMANN`, `AVOGADRO`, `GAS_CONSTANT`, `ELEMENTARY_CHARGE`, `ELECTRON_MASS`, `PROTON_MASS`, `ATOMIC_MASS_UNIT`, `ELECTRON_VOLT`, `HARTREE_ENERGY`, `BOHR_RADIUS`, `VACUUM_PERMITTIVITY`, `VACUUM_PERMEABILITY`, `COULOMB_CONSTANT`, `THERMOCHEMICAL_CALORIE`. |
| `units` | Scalar conversions between molecular-simulation and SI-style units: eV, J, Å, nm, m, amu, kg, K, kJ/mol and kcal/mol. |
| `quantity` | Newtype wrappers `Energy`, `Length`, `Mass` and `Temperature` around an `f64` in SI units, with named constructors and accessors (for example `Energy::from_ev` and `.kj_per_mol()`). |
| `vec3` | `Vec3`, a `Copy` 3-vector with arithmetic operators, `dot`, `cross`, `norm`, `norm_sq`, `normalized` (returns `Option`, so a zero vector is not silently normalised) and `hadamard`. |

The crate re-exports `Vec3` at the root, so `use tpt_dynamics_core::Vec3;` works.

## Quick start

```rust
use tpt_dynamics_core::units::{ev_to_kj_per_mol, kcal_per_mol_to_kj_per_mol};
use tpt_dynamics_core::Vec3;

// 1 eV per particle is about 96.485 kJ/mol.
let kj = ev_to_kj_per_mol(1.0);
assert!((kj - 96.485_332_12).abs() < 1e-6);
assert!((kcal_per_mol_to_kj_per_mol(1.0) - 4.184).abs() < 1e-12);

// Vector arithmetic.
let a = Vec3::new(1.0, 0.0, 0.0);
let b = Vec3::new(0.0, 1.0, 0.0);
assert_eq!(a.cross(b), Vec3::new(0.0, 0.0, 1.0));
assert!(Vec3::ZERO.normalized().is_none());
```

Typed quantities make the unit of a value explicit at the call site:

```rust
use tpt_dynamics_core::quantity::{Energy, Length};

let bond = Length::from_angstrom(1.09);
assert!((bond.nm() - 0.109).abs() < 1e-12);

let well_depth = Energy::from_kcal_per_mol(0.65);
assert!((well_depth.kj_per_mol() - 2.7196).abs() < 1e-3);
```

## Features

| Feature | Default | Effect |
| --- | --- | --- |
| `std` | no | Reserved for std-only conveniences. It is currently empty, so enabling it has no effect. |

## Design notes

- Units are not tracked by the type system except in `quantity`. Plain `f64`
  functions in `units` and `constants` are used by the numerical crates for speed
  and clarity. The conversion functions are the single place where factors live.
- `Vec3` is a small value type intended for particle-level code in
  `tpt-dynamics-sim`. The hot force kernels in `tpt-dynamics-compute` work on
  structure-of-arrays slices instead.
- The minimum supported Rust version is the workspace `rust-version` (see the root
  `Cargo.toml`).

## Licence

Licensed under either of [MIT](../LICENSE-MIT) or [Apache-2.0](../LICENSE-APACHE), at
your option. Part of [tpt-dynamics](../README.md).
