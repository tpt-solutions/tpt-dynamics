# Changelog

All notable changes to `tpt-dynamics-core` are documented in this file. The workspace
changelog is in the [root `CHANGELOG.md`](../CHANGELOG.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
crate uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0.

## [Unreleased]

### Added

- `constants`: CODATA 2018 values for the speed of light, Planck and reduced Planck
  constants, Boltzmann, Avogadro and gas constants, elementary charge, electron and
  proton masses, the atomic mass unit, the electronvolt, the Hartree energy, the Bohr
  radius, vacuum permittivity and permeability, the Coulomb constant and the
  thermochemical calorie.
- `units`: scalar conversions between eV, joules, Å, nm, metres, amu, kg, kelvin,
  kJ/mol and kcal/mol.
- `quantity`: `Energy`, `Length`, `Mass` and `Temperature` newtypes over SI values,
  with named constructors and accessors.
- `vec3`: `Vec3`, a `Copy` 3-vector with arithmetic operators, `dot`, `cross`,
  `norm`, `norm_sq`, `normalized` (returns `None` for a zero vector) and `hadamard`.
- `std` feature, off by default. It is reserved for std-only additions and is
  currently empty.
- `#![no_std]` and `#![forbid(unsafe_code)]`. The only dependency is `libm`.

### Known issues

- None.
