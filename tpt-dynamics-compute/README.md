# tpt-dynamics-compute

SIMD and mesh-based force kernels for [tpt-dynamics](../README.md). The kernels work
on plain structure-of-arrays slices, so `tpt-dynamics-sim` calls them directly and
this crate does not depend on any simulation types.

The crate has four jobs:

1. **Pair kernels** for Lennard-Jones and real-space Ewald Coulomb interactions, each
   as a scalar reference and as a 4-lane SIMD version.
2. **Neighbour search** with SIMD cell lists, Verlet lists and a padded 4-atom
   cluster layout with a cluster-pair kernel.
3. **Long-range electrostatics** with smooth particle-mesh Ewald (PME) on a pure-Rust
   radix-2 FFT, parallelised with rayon.
4. **Validation** against a brute-force Ewald sum, so that PME can be checked
   independently.

No `unsafe` code is used (`#![forbid(unsafe_code)]`). SIMD comes from
[`tpt-simd`](https://crates.io/crates/tpt-simd) (ADR-003), which is why the crate
requires Rust 1.95 or later.

## Modules

| Module | Contents |
| --- | --- |
| `pairwise` | Scalar formulas: `lj_pair`, `coulomb_ewald_pair`, `excluded_ewald_pair` and `lorentz_berthelot`. These define the reference behaviour. |
| `scalar` | `pair_list_scalar`: the scalar kernel over an `(i, j)` pair list. |
| `simd_kernel` | `pair_list_simd`: the same kernel, vectorised 4-wide. It agrees with the scalar kernel to round-off. |
| `listbuild` | SIMD cell-list pair search: `build_pair_list` and `build_verlet_list`. |
| `cluster` | Padded 4-atom clusters (`ClusterSet`), a bounding-sphere cluster-pair builder (`build_cluster_pairs`) and `cluster_pair_kernel`. |
| `pme` | Smooth PME: `Pme::from_accuracy` chooses the grid from a tolerance, `alpha_for_tolerance` chooses the Ewald splitting parameter, and `Pme::reciprocal` computes the reciprocal-space energy and forces. |
| `fft` | Radix-2 complex FFT in 1D (`Fft1d`) and 3D (`Fft3d`), rayon-parallel. |
| `reference` | Brute-force Ewald sums (`direct_real_images`, `direct_reciprocal`, `direct_self_energy`) used by the tests and for validation. |
| `types` | `Positions`, `ForceBuffers`, `Energies`, `AtomParams` and `NonbondedConfig`, the shared input and output types. |

## Conventions

- **Units.** The Coulomb constant is 1. The caller chooses the length and energy units
  and is responsible for any prefactor (for example the 138.935 of GROMACS-style
  units in kJ/mol).
- **Return values.** The pair formulas return `(energy, f_over_r)`: the energy of the
  pair and the force divided by the distance. The caller multiplies the separation
  vector by `f_over_r` to get the force on each atom.
- **Forces are added.** Kernels add into `ForceBuffers`. Zero the buffers before a
  step, and sum contributions from several kernels into the same buffers.
- **Positive `r2`.** Pair functions require `r2 > 0`. Coincident atoms are the
  caller's responsibility.

## Quick start

A Lennard-Jones pair at the potential minimum, `r = 2^(1/6) sigma`, has energy
`-epsilon`:

```rust
use tpt_dynamics_compute::pairwise::{coulomb_ewald_pair, lj_pair};
use tpt_dynamics_compute::pme::alpha_for_tolerance;

let r2 = 2.0_f64.powf(1.0 / 3.0); // (2^(1/6))^2
let (energy, f_over_r) = lj_pair(r2, 1.0, 1.0);
assert!((energy + 1.0).abs() < 1e-12);
assert!(f_over_r.abs() < 1e-12); // zero force at the minimum

// Real-space Ewald Coulomb for two opposite unit charges at r = 1 nm.
let alpha = alpha_for_tolerance(1.0, 1e-5); // erfc(alpha * rc) = tolerance at rc = 1
let (e_real, _) = coulomb_ewald_pair(1.0, -1.0, alpha);
assert!(e_real < 0.0);
```

Benchmarks:

```text
cargo run --release -p tpt-dynamics-compute --example kernel_bench
```

The harness reports the minimum and median wall-clock time per kernel. Numbers depend
strongly on the CPU and build flags; the recorded results and the method are in
[`docs/benchmarks/compute-kernels.md`](../docs/benchmarks/compute-kernels.md).

## Accuracy

`alpha_for_tolerance(rc, tol)` solves `erfc(alpha rc) = tol` for the Ewald splitting
parameter. `Pme::from_accuracy` uses the same target to choose the mesh, with spacing
at most `DEFAULT_MAX_SPACING` and dimensions rounded up to a power of two. The
reference Ewald sums in `reference` are the yardstick for both the real-space and
reciprocal-space paths.

## Licence

Licensed under either of [MIT](../LICENSE-MIT) or [Apache-2.0](../LICENSE-APACHE), at
your option. Part of [tpt-dynamics](../README.md).
