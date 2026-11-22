# Changelog

All notable changes to `tpt-dynamics-compute` are documented in this file. The
workspace changelog is in the [root `CHANGELOG.md`](../CHANGELOG.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
crate uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0.

## [Unreleased]

### Added

- `pairwise`: scalar Lennard-Jones, real-space Ewald Coulomb and excluded-pair
  formulas, Lorentz–Berthelot combining rules, and `alpha_for_tolerance`. These are
  the reference behaviour.
- `scalar`: `pair_list_scalar`, the scalar kernel over an `(i, j)` pair list.
- `simd_kernel`: `pair_list_simd`, the same kernel vectorised 4-wide with `tpt-simd`
  (ADR-003). It agrees with the scalar kernel to round-off.
- `listbuild`: SIMD cell-list pair search, with `build_pair_list` and
  `build_verlet_list`.
- `cluster`: padded 4-atom cluster layout (`ClusterSet`), a bounding-sphere
  cluster-pair builder and a SIMD cluster-pair kernel.
- `pme`: smooth particle-mesh Ewald. `Pme::from_accuracy` chooses the grid from an
  accuracy target, and `Pme::reciprocal` returns the reciprocal-space energy and
  adds forces.
- `fft`: radix-2 complex FFT in one and three dimensions, rayon-parallel.
- `reference`: brute-force Ewald sums for validating PME.
- `types`: `Positions`, `ForceBuffers`, `Energies`, `AtomParams` and
  `NonbondedConfig`.
- `examples/kernel_bench`: a timing harness for the kernels. Recorded results are in
  `docs/benchmarks/compute-kernels.md`.

### Changed

- The minimum supported Rust version is 1.95, because `tpt-simd` 0.1.1 requires it.

### Known issues

- Benchmark numbers depend strongly on the CPU and build flags, and are only
  comparable on the same machine.
