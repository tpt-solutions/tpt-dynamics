# tpt-dynamics — Project Checklist

Pure-Rust, zero-FFI replacement for ROOT (HEP analysis) and GROMACS (molecular dynamics).
**Owner:** TPT Solutions · **Licence:** MIT OR Apache-2.0 · © 2026 TPT Solutions
**Source of truth:** `spec.txt` (v0.1). Phases 0–4 are the path to v0.1.0; Phases 5–8 are post-0.1.

Rules for every dependency: MIT, BSD-2/3, Zlib, ISC, Unicode, or dual `MIT OR Apache-2.0`. **No Apache-2.0-only crates.** No GPL/AGPL/LGPL/MPL. No FFI/bindgen/`*-sys`.

Legend: `- [ ]` todo · `- [x]` done · **bold** = milestone

---

## Phase 0 — Repo, Licensing & Governance

### Repository setup
- [ ] `git init`, `.gitignore`, default branch, remote at github.com/tpt-solutions
- [ ] Root `Cargo.toml` workspace with members: `tpt-dynamics-core`, `-data`, `-sim`, `-compute`, `-cli`, `bindings/tpt-dynamics-py`
- [ ] `[workspace.package]`: version, edition, `rust-version` (MSRV), `license = "MIT OR Apache-2.0"`, authors/repository
- [ ] `rust-toolchain.toml`, `rustfmt.toml`, clippy config
- [ ] Workspace lints (`unsafe_code` policy, `missing_docs`, clippy pedantic subset)
- [ ] README (overview, status, licence section, build instructions)

### Licence & governance files
- [ ] `LICENSE-MIT` (© 2026 TPT Solutions)
- [ ] `LICENSE-APACHE` (Apache-2.0 text, TPT Solutions notice)
- [ ] `license = "MIT OR Apache-2.0"` and `license-file` correct in every crate manifest (incl. Python bindings / `pyproject.toml`)
- [ ] README licence section with standard "contributions dual-licensed" clause
- [ ] `CONTRIBUTING.md` (inbound = outbound licensing, dependency rules)
- [ ] `SECURITY.md`, `CODE_OF_CONDUCT.md`, issue/PR templates
- [ ] `CHANGELOG.md` (Keep a Changelog)

### Dependency compliance (cargo-deny)
- [ ] Write `deny.toml` using the **current** cargo-deny schema (spec's `unlicensed` / `copyleft` keys are deprecated)
- [ ] Allow list: MIT, BSD-2-Clause, BSD-3-Clause, Zlib, ISC, Unicode-3.0, and Apache-2.0 **only** via `MIT OR Apache-2.0` expressions
- [ ] Mechanism to reject Apache-only crates (e.g. `exceptions` / `licenses.clarify` / expression check script) — **deviation from spec's plain `Apache-2.0` allow entry**
- [ ] Deny list: GPL-2.0, GPL-3.0, AGPL-3.0, LGPL-2.1, LGPL-3.0, MPL-2.0
- [ ] `[bans]`: `root`, `gromacs`, `geant4`, `arrow`, `fftw-sys`, `mpi`/`rsmpi` (fix spec's invalid multi-`name` TOML entry — one `[[bans.deny]]` per crate)
- [ ] Ban/flag `*-sys`, `bindgen`, `cc`/`cmake` build deps in the dependency tree (zero-FFI policy)
- [ ] `[sources]`: allow only crates.io + `github.com/tpt-solutions`
- [ ] Add `cargo deny check` as a CI gate and pre-release gate
- [ ] Test the gate: add a known GPL and a known Apache-only crate on a scratch branch and confirm the build fails

### Dependency audit (licence + FFI)
- [ ] `tpt-math` (crates.io): licence, `no_std` support, FFT API, pin version
- [ ] `tpt-simd` (crates.io): licence, supported targets, stable-Rust compatibility, pin version
- [ ] `tpt-strata` (+ `tpt-strata-parquet`): licence of full tree via `cargo tree` / `cargo deny list`
- [ ] Decide on `tpt-strata-parquet`: accept only if every Arrow/Parquet transitive dep is MIT/BSD/dual — otherwise avoid or isolate behind an optional feature
- [ ] `glam` (no_std/libm feature), `clap`, `pyo3`, `rayon`, `maturin` licences
- [ ] Record results in `docs/dependencies.md`

### CI
- [ ] GitHub Actions: `fmt`, `clippy -D warnings`, `test`, `doc`
- [ ] `cargo deny check` job
- [ ] `no_std` build job for `tpt-dynamics-core` (bare-metal target, e.g. `thumbv7em-none-eabihf`)
- [ ] MSRV job
- [ ] OS matrix: Linux, Windows, macOS; arch: x86_64 + aarch64
- [ ] Dependabot/Renovate with licence-aware review

### Docs & decisions
- [ ] `docs/` skeleton (mdBook or rustdoc-only)
- [ ] `docs/adr/` folder: ADR-001 (Pure-Rust MD), ADR-002 (Lazy data frames), ADR-003 (SIMD kernels)
- [ ] New ADR-004: build data engine on tpt-strata (no polars/arrow)
- [ ] New ADR-005: licence policy (no Apache-only deps, deviation from spec's deny.toml)
- [ ] Unsafe-code policy: `#![forbid(unsafe_code)]` by default; any exception documented and reviewed

---

## Phase 1 — Core & Data Engine (spec months 1–3)

### tpt-dynamics-core (`no_std + alloc`, leaf crate)
- [ ] Crate skeleton with `#![no_std]`, `extern crate alloc`, optional `std` feature
- [ ] Physical constants (Planck, c, Boltzmann, Avogadro, elementary charge, etc.) with CODATA source/year noted
- [ ] Unit conversions: eV↔J, Å↔nm, kJ/mol↔kcal/mol, amu↔kg, K↔energy
- [ ] Optional typed units (newtypes) to prevent unit mix-ups
- [ ] Vector math layer over `tpt-math` / `glam` (confirm `no_std`/`libm`)
- [ ] Unit tests + doc tests; property tests for conversion round-trips
- [ ] Verify `no_std` build on embedded target in CI
- [ ] Rustdoc with examples

### tpt-dynamics-data — foundation on tpt-strata
- [ ] Add `tpt-strata` dependency; pin version
- [ ] API gap analysis: strata features vs spec needs (lazy frames, streaming, UINT64, binary columns, n-tuples)
- [ ] Open upstream issues/PRs on tpt-strata for gaps (lazy execution, chunked streaming, UINT64, binary/list columns)
- [ ] Define physics-layer column model (jagged/variable-length arrays for per-event particles)

### Lazy data frame
- [ ] `LazyFrame` API: scan → filter → select/with_column → group_by → aggregate → collect
- [ ] Logical plan representation
- [ ] Optimiser: predicate pushdown, projection pushdown, expression fusion
- [ ] Chunked/streaming executor (bounded memory)
- [ ] Parallel execution across chunks (rayon)
- [ ] Plan `explain()` output
- [ ] Expression engine for physics quantities (invariant mass, pT, ΔR, etc.)
- [ ] Tests: results equal eager evaluation; memory-bound tests on large inputs

### Physics histograms
- [ ] 1D histogram: fixed, variable, and dynamic binning; under/overflow; weights; sum of squared weights
- [ ] 2D histogram
- [ ] Profile histograms (1D, 2D)
- [ ] Efficient `merge` (associative, commutative) for distributed use
- [ ] Serialisation (binary + JSON)
- [ ] Fill from lazy frame aggregation
- [ ] Statistics: mean, std dev, integral, quantiles
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
- [ ] Core MD data structures: atoms, topology, box, velocities (SoA layout)
- [ ] PDB reader/writer
- [ ] GRO reader/writer
- [ ] Simulation config format (TOML)
- [ ] Mass/charge/type assignment

### Boundary conditions & neighbour search
- [ ] Periodic boundary conditions (orthorhombic; triclinic planned)
- [ ] Minimum-image convention
- [ ] Cell-list algorithm
- [ ] Verlet list with skin and rebuild criterion
- [ ] Parallel neighbour-list build (rayon)
- [ ] Tests vs brute-force O(N²)

### Forces & integrators
- [ ] Lennard-Jones potential (cutoff, shift/switch)
- [ ] Velocity Verlet integrator
- [ ] Leap-frog integrator
- [ ] Pluggable integrator / force-field traits (runtime swappable, per ADR-001)
- [ ] Basic thermostat for stable runs (Berendsen or v-rescale)
- [ ] Energy / momentum conservation tests (NVE)
- [ ] Reference validation: LJ fluid vs published thermodynamic data

### Output
- [ ] Trajectory writer hook (plain binary initially)
- [ ] Energy/log output

- [ ] **Milestone: `tptdyn simulate` runs a stable 10,000-atom Lennard-Jones fluid**

---

## Phase 3 — Scale & Acceleration (spec months 7–9)

### tpt-dynamics-compute: SIMD kernels
- [ ] Integrate `tpt-simd`; scalar fallback path
- [ ] SIMD Lennard-Jones pairwise kernel
- [ ] SIMD Coulomb pairwise kernel (with Ewald real-space term)
- [ ] Data layout for SIMD (SoA, padded cluster pair lists)
- [ ] SIMD cell-list / Verlet-list build
- [ ] Correctness tests vs scalar reference
- [ ] CI coverage on x86_64 and aarch64
- [ ] Kernel micro-benchmarks

### PME (Particle Mesh Ewald)
- [ ] Charge spreading (B-splines)
- [ ] 3D FFT via `tpt-math`
- [ ] Reciprocal-space solve, force interpolation
- [ ] Self-energy and exclusion corrections
- [ ] Parameter selection (grid size, order, Ewald coefficient) from accuracy target
- [ ] Accuracy tests vs direct Ewald summation
- [ ] Parallelise spread/gather/FFT

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
- [ ] `clap` command skeleton and global options
- [ ] `tptdyn analyze`
- [ ] `tptdyn simulate`
- [ ] `tptdyn histogram`
- [ ] Config-file support, progress reporting, logging
- [ ] Clear error messages and exit codes
- [ ] Shell completions and man pages
- [ ] CLI integration tests

### tpt-dynamics-py (PyO3)
- [ ] PyO3 + maturin project setup (MIT OR Apache-2.0 metadata)
- [ ] Bindings: lazy frames, histograms, simulation setup/run
- [ ] NumPy interop (check licence of `numpy`/`rust-numpy` crate)
- [ ] Type stubs (`.pyi`) and Python docs
- [ ] PyROOT-style migration examples
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

- [ ] Force-field parameter file formats and topology builder
- [ ] Bonded terms: bonds, angles, proper/improper dihedrals, Urey-Bradley
- [ ] 1-4 interactions and exclusions
- [ ] AMBER parameter sets
- [ ] CHARMM parameter sets
- [ ] GROMOS parameter sets
- [ ] Check redistribution licence of each parameter set before bundling
- [ ] Support for proteins, lipids, nucleic acids
- [ ] SHAKE constraints
- [ ] LINCS constraints (spec writes "LINCOLN"; assumed to mean LINCS)
- [ ] SETTLE for rigid water models (TIP3P, SPC/E, TIP4P)
- [ ] Thermostats: Nosé-Hoover, v-rescale
- [ ] Barostats: Berendsen, Parrinello-Rahman / C-rescale
- [ ] Triclinic boxes
- [ ] Solvation and ion placement tools
- [ ] Energy minimisation (steepest descent, L-BFGS)
- [ ] Checkpoint / restart
- [ ] Validation vs reference energies/forces from established engines
- [ ] Demo: multi-million-atom viral capsid / ribosome system

---

## Phase 6 — Trajectory Output & Analysis

- [ ] Compressed columnar trajectory format (Zarr-compatible, or custom binary fallback)
- [ ] Streaming trajectory writer during simulation
- [ ] Trajectory reader via lazy engine (no full load into RAM)
- [ ] RMSD
- [ ] Radius of gyration
- [ ] Hydrogen-bond counts
- [ ] Radial distribution function
- [ ] Memory-bound tests on multi-GB trajectories
- [ ] `tptdyn analyze` support for trajectories

---

## Phase 7 — HEP Formats & Distributed Analysis

- [ ] RNTuple reader (clean-room from public spec, no FFI)
- [ ] ROOT TTree/file reader (clean-room)
- [ ] Arrow-IPC-style columnar reader (own implementation, no `arrow` crate)
- [ ] Jagged-array operations for per-event particle collections
- [ ] Distributed execution of lazy plans across nodes
- [ ] Histogram merging across nodes
- [ ] Python histogram interop (boost-histogram/hist-style API)
- [ ] Network/object-store data sources

---

## Phase 8 — GPU, Embedded & Hardening

- [ ] GPU backend for `tpt-dynamics-compute` (evaluate wgpu etc.; licence-check first)
- [ ] GPU force kernels and PME
- [ ] Embedded / FPGA DAQ demo using `no_std` core
- [ ] Fuzzing campaign for all parsers
- [ ] Unsafe-code audit
- [ ] Security review
- [ ] Performance regression tracking in CI
- [ ] Reproducibility / determinism mode

---

## Ongoing / Cross-cutting

- [ ] Licence review on every new or updated dependency
- [ ] Keep `deny.toml` and ADRs current
- [ ] Semver and MSRV policy
- [ ] Issue triage and roadmap review each phase
- [ ] Keep docs and benchmarks current
