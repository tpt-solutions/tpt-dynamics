# tpt-dynamics — Project Checklist

Pure-Rust, zero-FFI replacement for ROOT (HEP analysis) and GROMACS (molecular dynamics).
**Owner:** TPT Solutions · **Licence:** MIT OR Apache-2.0 · © 2026 TPT Solutions
**Source of truth:** `spec.txt` (v0.1). Phases 0–4 are the path to v0.1.0; Phases 5–8 are post-0.1.

Rules for every dependency: MIT, BSD-2/3, Zlib, ISC, Unicode, or dual `MIT OR Apache-2.0`. **No Apache-2.0-only crates.** No GPL/AGPL/LGPL/MPL. No FFI/bindgen/`*-sys`.

Legend: `- [ ]` todo · `- [x]` done · **bold** = milestone

---

## Phase 0 — Repo, Licensing & Governance

### Repository setup
- [x] `git init`, `.gitignore`, default branch, remote at github.com/tpt-solutions
- [x] Root `Cargo.toml` workspace with members: `tpt-dynamics-core`, `-data`, `-sim`, `-compute`, `-cli`, `bindings/tpt-dynamics-py`
- [x] `[workspace.package]`: version, edition, `rust-version` (MSRV), `license = "MIT OR Apache-2.0"`, authors/repository
- [x] `rust-toolchain.toml`, `rustfmt.toml`, clippy config
- [x] Workspace lints (`unsafe_code` policy, `missing_docs`, clippy pedantic subset)
- [x] README (overview, status, licence section, build instructions)

### Licence & governance files
- [x] `LICENSE-MIT` (© 2026 TPT Solutions)
- [ ] `LICENSE-APACHE` (Apache-2.0 text, TPT Solutions notice)
- [x] `license = "MIT OR Apache-2.0"` and `license-file` correct in every crate manifest (incl. Python bindings / `pyproject.toml`)
- [x] README licence section with standard "contributions dual-licensed" clause
- [x] `CONTRIBUTING.md` (inbound = outbound licensing, dependency rules)
- [ ] `SECURITY.md`, `CODE_OF_CONDUCT.md`, issue/PR templates
- [x] `CHANGELOG.md` (Keep a Changelog)

### Dependency compliance (cargo-deny)
- [x] Write `deny.toml` using the **current** cargo-deny schema (spec's `unlicensed` / `copyleft` keys are deprecated)
- [x] Allow list: MIT, BSD-2-Clause, BSD-3-Clause, Zlib, ISC, Unicode-3.0, and Apache-2.0 **only** via `MIT OR Apache-2.0` expressions
- [x] Mechanism to reject Apache-only crates (e.g. `exceptions` / `licenses.clarify` / expression check script) — **deviation from spec's plain `Apache-2.0` allow entry**
- [x] Deny list: GPL-2.0, GPL-3.0, AGPL-3.0, LGPL-2.1, LGPL-3.0, MPL-2.0
- [x] `[bans]`: `root`, `gromacs`, `geant4`, `arrow`, `fftw-sys`, `mpi`/`rsmpi` (fix spec's invalid multi-`name` TOML entry — one `[[bans.deny]]` per crate)
- [ ] Ban/flag `*-sys`, `bindgen`, `cc`/`cmake` build deps in the dependency tree (zero-FFI policy)
- [x] `[sources]`: allow only crates.io + `github.com/tpt-solutions`
- [ ] Add `cargo deny check` as a CI gate and pre-release gate
- [ ] Test the gate: add a known GPL and a known Apache-only crate on a scratch branch and confirm the build fails

### Dependency audit (licence + FFI)
- [x] `tpt-math` (crates.io): licence, `no_std` support, FFT API, pin version
- [x] `tpt-simd` (crates.io): licence, supported targets, stable-Rust compatibility, pin version
- [x] `tpt-strata` (+ `tpt-strata-parquet`): licence of full tree via `cargo tree` / `cargo deny list`
- [x] Decide on `tpt-strata-parquet`: accept only if every Arrow/Parquet transitive dep is MIT/BSD/dual — otherwise avoid or isolate behind an optional feature
- [x] `glam` (no_std/libm feature), `clap`, `pyo3`, `rayon`, `maturin` licences
- [x] Record results in `docs/dependencies.md`

### CI
- [x] GitHub Actions: `fmt`, `clippy -D warnings`, `test`, `doc`
- [x] `cargo deny check` job
- [x] `no_std` build job for `tpt-dynamics-core` (bare-metal target, e.g. `thumbv7em-none-eabihf`)
- [x] MSRV job
- [x] OS matrix: Linux, Windows, macOS; arch: x86_64 + aarch64
- [x] Dependabot/Renovate with licence-aware review

### Docs & decisions
- [x] `docs/` skeleton (mdBook or rustdoc-only)
- [x] `docs/adr/` folder: ADR-001 (Pure-Rust MD), ADR-002 (Lazy data frames), ADR-003 (SIMD kernels)
- [x] New ADR-004: build data engine on tpt-strata (no polars/arrow)
- [x] New ADR-005: licence policy (no Apache-only deps, deviation from spec's deny.toml)
- [x] Unsafe-code policy: `#![forbid(unsafe_code)]` by default; any exception documented and reviewed

---

## Phase 1 — Core & Data Engine (spec months 1–3)

### tpt-dynamics-core (`no_std + alloc`, leaf crate)
- [x] Crate skeleton with `#![no_std]`, `extern crate alloc`, optional `std` feature
- [x] Physical constants (Planck, c, Boltzmann, Avogadro, elementary charge, etc.) with CODATA source/year noted
- [x] Unit conversions: eV↔J, Å↔nm, kJ/mol↔kcal/mol, amu↔kg, K↔energy
- [x] Optional typed units (newtypes) to prevent unit mix-ups
- [x] Vector math layer over `tpt-math` / `glam` (confirm `no_std`/`libm`)
- [ ] Unit tests + doc tests; property tests for conversion round-trips
- [x] Verify `no_std` build on embedded target in CI
- [x] Rustdoc with examples

### tpt-dynamics-data — foundation on tpt-strata
- [x] Add `tpt-strata` dependency; pin version
- [ ] API gap analysis: strata features vs spec needs (lazy frames, streaming, UINT64, binary columns, n-tuples)
- [ ] Open upstream issues/PRs on tpt-strata for gaps (lazy execution, chunked streaming, UINT64, binary/list columns)
- [x] Define physics-layer column model (jagged/variable-length arrays for per-event particles)

### Lazy data frame
- [x] `LazyFrame` API: scan → filter → select/with_column → group_by → aggregate → collect
- [x] Logical plan representation
- [ ] Optimiser: predicate pushdown, projection pushdown, expression fusion
- [x] Chunked/streaming executor (bounded memory)
- [x] Parallel execution across chunks (rayon)
- [x] Plan `explain()` output
- [x] Expression engine for physics quantities (invariant mass, pT, ΔR, etc.)
- [x] Tests: results equal eager evaluation; memory-bound tests on large inputs

### Physics histograms
- [x] 1D histogram: fixed, variable, and dynamic binning; under/overflow; weights; sum of squared weights
- [x] 2D histogram
- [x] Profile histograms (1D, 2D)
- [x] Efficient `merge` (associative, commutative) for distributed use
- [x] Serialisation (binary + JSON)
- [ ] Fill from lazy frame aggregation
- [x] Statistics: mean, std dev, integral, quantiles
- [ ] Property tests for merge correctness

### I/O
- [ ] Streaming CSV reader for very large files
- [ ] Parquet read/write via tpt-strata bridge (subject to Phase 0 licence decision)
- [ ] Custom binary columnar format for detector data (spec, versioning, checksums, compression)
- [ ] Zero-copy / memory-mapped readers where possible
- [ ] Fuzz tests for all parsers

### Benchmarks
- [ ] Generate synthetic collision-event dataset (10 GB CSV)
- [ ] Benchmark harness (criterion) + baseline PyROOT/RDataFrame comparison script
- [ ] Record results in `docs/benchmarks/`

- [ ] **Milestone: `tptdyn analyze` processes a 10 GB CSV of collision events and produces histograms faster than PyROOT**

---

## Phase 2 — Simulation Engine (spec months 4–6)

### System model & I/O
- [x] Core MD data structures: atoms, topology, box, velocities (SoA layout)
- [x] PDB reader/writer
- [x] GRO reader/writer
- [x] Simulation config format (TOML)
- [x] Mass/charge/type assignment

### Boundary conditions & neighbour search
- [x] Periodic boundary conditions (orthorhombic; triclinic planned)
- [x] Minimum-image convention
- [x] Cell-list algorithm
- [x] Verlet list with skin and rebuild criterion
- [x] Parallel neighbour-list build (rayon)
- [x] Tests vs brute-force O(N²)

### Forces & integrators
- [x] Lennard-Jones potential (cutoff, shift/switch)
- [x] Velocity Verlet integrator
- [x] Leap-frog integrator
- [x] Pluggable integrator / force-field traits (runtime swappable, per ADR-001)
- [x] Basic thermostat for stable runs (Berendsen or v-rescale)
- [x] Energy / momentum conservation tests (NVE)
- [ ] Reference validation: LJ fluid vs published thermodynamic data

### Output
- [x] Trajectory writer hook (plain binary initially)
- [x] Energy/log output

- [ ] **Milestone: `tptdyn simulate` runs a stable 10,000-atom Lennard-Jones fluid**

---

## Phase 3 — Scale & Acceleration (spec months 7–9)

### tpt-dynamics-compute: SIMD kernels
- [x] Integrate `tpt-simd`; scalar fallback path
- [x] SIMD Lennard-Jones pairwise kernel
- [x] SIMD Coulomb pairwise kernel (with Ewald real-space term)
- [x] Data layout for SIMD (SoA, padded cluster pair lists)
- [ ] SIMD cell-list / Verlet-list build
- [x] Correctness tests vs scalar reference
- [ ] CI coverage on x86_64 and aarch64
- [x] Kernel micro-benchmarks

### PME (Particle Mesh Ewald)
- [x] Charge spreading (B-splines)
- [ ] 3D FFT via `tpt-math`
- [x] Reciprocal-space solve, force interpolation
- [x] Self-energy and exclusion corrections
- [x] Parameter selection (grid size, order, Ewald coefficient) from accuracy target
- [x] Accuracy tests vs direct Ewald summation
- [x] Parallelise spread/gather/FFT

### Domain decomposition
- [ ] Spatial decomposition of the simulation box into sub-domains
- [ ] Halo / ghost-atom exchange
- [ ] Atom migration between domains
- [ ] Pure-Rust message passing layer (no `mpi` crate)
- [ ] rayon-based shared-memory execution
- [ ] Dynamic load balancing
- [ ] Reproducibility tests: decomposed vs single-domain results

### Scaling
- [ ] 1M-atom benchmark system
- [ ] Strong/weak scaling benchmarks up to 32+ cores
- [ ] Profile and remove bottlenecks (memory bandwidth, false sharing)

- [ ] **Milestone: `tptdyn simulate` scales linearly to 32+ CPU cores on a 1-million-atom system**

---

## Phase 4 — CLI, Python Bindings & v0.1.0 Release (spec months 10–12)

### tpt-dynamics-cli (`tptdyn`)
- [x] `clap` command skeleton and global options
- [x] `tptdyn analyze`
- [x] `tptdyn simulate`
- [x] `tptdyn histogram`
- [x] Config-file support, progress reporting, logging
- [x] Clear error messages and exit codes
- [x] Shell completions and man pages
- [x] CLI integration tests

### tpt-dynamics-py (PyO3)
- [x] PyO3 + maturin project setup (MIT OR Apache-2.0 metadata)
- [x] Bindings: lazy frames, histograms, simulation setup/run
- [x] NumPy interop (check licence of `numpy`/`rust-numpy` crate)
- [x] Type stubs (`.pyi`) and Python docs
- [x] PyROOT-style migration examples
- [ ] Wheels: Linux, Windows, macOS (x86_64 + arm64)
- [ ] Python test suite in CI

### Release
- [ ] Documentation site and tutorials
- [ ] Examples directory (HEP analysis, LJ fluid, water box)
- [ ] Benchmarks vs GROMACS and ROOT published
- [ ] Final `cargo deny check` and full-tree licence review
- [ ] Final `cargo publish --dry-run` for all crates
- [ ] Publish order: core → compute/data → sim → cli
- [ ] Publish to PyPI
- [ ] Tag v0.1.0, release notes, changelog

- [ ] **Milestone: v0.1.0 released on crates.io and PyPI with benchmarks vs GROMACS and ROOT**

---

## Phase 5 — Force Fields & MD Completeness (post-0.1)

- [x] Force-field parameter file formats and topology builder
- [x] Bonded terms: bonds, angles, proper/improper dihedrals, Urey-Bradley
- [x] 1-4 interactions and exclusions
- [ ] AMBER parameter sets
- [ ] CHARMM parameter sets
- [ ] GROMOS parameter sets
- [ ] Check redistribution licence of each parameter set before bundling
- [ ] Support for proteins, lipids, nucleic acids
- [x] SHAKE constraints
- [ ] LINCS constraints (spec writes "LINCOLN"; assumed to mean LINCS)
- [ ] SETTLE for rigid water models (TIP3P, SPC/E, TIP4P)
- [x] Thermostats: Nosé-Hoover, v-rescale
- [x] Barostats: Berendsen, Parrinello-Rahman / C-rescale
- [x] Triclinic boxes
- [x] Solvation and ion placement tools
- [x] Energy minimisation (steepest descent, L-BFGS)
- [x] Checkpoint / restart
- [ ] Validation vs reference energies/forces from established engines
- [ ] Demo: multi-million-atom viral capsid / ribosome system

---

## Phase 6 — Trajectory Output & Analysis

- [x] Compressed columnar trajectory format (Zarr-compatible, or custom binary fallback)
- [ ] Streaming trajectory writer during simulation
- [x] Trajectory reader via lazy engine (no full load into RAM)
- [x] RMSD
- [x] Radius of gyration
- [x] Hydrogen-bond counts
- [x] Radial distribution function
- [ ] Memory-bound tests on multi-GB trajectories
- [ ] `tptdyn analyze` support for trajectories

---

## Phase 7 — HEP Formats & Distributed Analysis

- [ ] RNTuple reader (clean-room from public spec, no FFI)
- [ ] ROOT TTree/file reader (clean-room)
- [x] Arrow-IPC-style columnar reader (own implementation, no `arrow` crate)
- [x] Jagged-array operations for per-event particle collections
- [ ] Distributed execution of lazy plans across nodes
- [x] Histogram merging across nodes
- [ ] Python histogram interop (boost-histogram/hist-style API)
- [ ] Network/object-store data sources

---

## Phase 8 — GPU, Embedded & Hardening

- [ ] GPU backend for `tpt-dynamics-compute` (evaluate wgpu etc.; licence-check first)
- [ ] GPU force kernels and PME
- [x] Embedded / FPGA DAQ demo using `no_std` core
- [ ] Fuzzing campaign for all parsers
- [x] Unsafe-code audit
- [x] Security review
- [x] Performance regression tracking in CI
- [ ] Reproducibility / determinism mode

---

## Ongoing / Cross-cutting

- [ ] Licence review on every new or updated dependency
- [ ] Keep `deny.toml` and ADRs current
- [ ] Semver and MSRV policy
- [ ] Issue triage and roadmap review each phase
- [ ] Keep docs and benchmarks current
