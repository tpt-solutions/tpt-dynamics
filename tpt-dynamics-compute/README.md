# tpt-dynamics-compute

SIMD and mesh-based force kernels for tpt-dynamics, operating on structure-of-arrays
slices so that `tpt-dynamics-sim` can call them directly.

* Lennard-Jones and real-space Ewald Coulomb pair kernels: scalar reference and
  `tpt-simd` 4-lane version.
* Padded 4-atom cluster layout with a SIMD cluster-pair kernel.
* Particle-mesh Ewald (smooth PME) with a pure-Rust radix-2 FFT, rayon-parallel.
* Brute-force Ewald reference for validation.

Benchmarks: `cargo run --release -p tpt-dynamics-compute --example kernel_bench`;
recorded numbers in [`docs/benchmarks/compute-kernels.md`](../docs/benchmarks/compute-kernels.md).

Part of [tpt-dynamics](../README.md). Licensed under MIT OR Apache-2.0.
