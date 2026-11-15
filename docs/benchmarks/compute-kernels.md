# tpt-dynamics-compute kernel benchmarks

Recorded 2026-11-16 with `cargo run --release -p tpt-dynamics-compute --example kernel_bench`.

**Machine:** Intel Core i7-12700 (12 cores, 20 logical), Windows 11, rayon pool of 20 threads.
Rust 1.100.0 stable, release profile (`lto = "thin"`, `codegen-units = 1`).
Target features are the default x86_64 baseline (no `-C target-cpu=native`), so the
tpt-simd lane code uses 128-bit SSE2 registers.

Timings are the minimum and median of 9 runs after one warm-up. Numbers are from a
single machine and are indicative only.

| Kernel | System | Min | Median | Notes |
|--------|--------|-----|--------|-------|
| Pair list, scalar reference | N=4096, rho=0.8, rc=2.5, 106,734 pairs | 4.60 ms | 4.62 ms | LJ + Ewald real-space Coulomb |
| Pair list, tpt-simd (4 lanes) | same | 4.15 ms | 4.28 ms | about 8% faster than scalar |
| Cluster pairs, atoms in input order | same, 524,798 cluster pairs | 142.7 ms | 144.8 ms | bounding spheres do not filter: clusters are spatially random |
| Cluster pairs, spatially sorted | same, 111,097 cluster pairs | 23.8 ms | 24.3 ms | atoms sorted by rc-sized cells before clustering |
| PME reciprocal | N=512, 8^3 box, grid 64^3, order 4 | 22.1 ms | 23.8 ms | spread, FFT, solve and gather, 20 threads |

## Interpretation

* The 4-lane kernel gains little because `erfc` and `exp` are evaluated per lane with
  `libm` (no vector `erfc` in tpt-simd). The arithmetic is vectorised; the transcendental
  calls dominate. A vector erfc would be the next step.
* The cluster kernel is slower than the pair-list kernel. Each cluster pair evaluates a
  4x4 block, so most lanes are outside the cutoff, and the cluster list in this benchmark
  holds about 10x more atom-pair work than the pair list. Atom ordering is decisive:
  sorting cuts the cluster list from 524,798 to 111,097 pairs and the time by about 6x.
  The list build is brute force (O(clusters^2)) and is not timed here; the simulation
  crate should provide the spatial sort and cell-based list build.
* PME spreading uses one private grid per thread and then sums them. An earlier version
  allocated a grid per rayon split and took about 640 ms for the same case.
