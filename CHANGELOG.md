# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0.

## [Unreleased]

### Added

- Cargo workspace with six members: `tpt-dynamics-core`, `-data`, `-sim`,
  `-compute`, `-cli` and `bindings/tpt-dynamics-py`.
- `tpt-dynamics-core` (`no_std`): CODATA 2018 physical constants, scalar unit
  conversions (eV, Å, nm, kJ/mol, kcal/mol, amu, K), typed quantities
  (`Energy`, `Length`, `Mass`, `Temperature`) and a `Vec3` type.
- `tpt-dynamics-data`: `Histogram1D` with fixed and variable binning, under- and
  overflow, weights, sum of squared weights, mean, standard deviation, quantiles,
  and merging.
- `tptdyn` binary with `--help` and `--version` (other subcommands are not yet
  implemented).
- Licence files (MIT and Apache-2.0), `deny.toml` with the project licence policy,
  and the dependency audit in `docs/dependencies.md`.
- Architecture decision records ADR-001 to ADR-006 in `docs/adr/`.
- Continuous integration: formatting, clippy, tests on three operating systems,
  documentation build, `cargo deny`, `no_std` build, and an MSRV job.

### Added (Phases 1–8 implementation pass, unreleased)

- `tpt-dynamics-data`: 2D, profile histograms; JSON and checksummed binary serialisation; streaming CSV reader; checksummed columnar format; physics helpers; a lazy data frame with predicate and projection pushdown and chunked, parallel execution; jagged arrays; a read-only Arrow-IPC-style reader (not Arrow-compatible); compressed chunked trajectory format; distributed-merge and fuzz-style tests.
- `tpt-dynamics-sim`: PDB/GRO I/O, TOML config, cell and Verlet lists, Lennard-Jones, bonded terms, 1-4 and exclusions, SHAKE, SETTLE, Nosé–Hoover, v-rescale, Berendsen, C-rescale, triclinic boxes, solvation and ion placement, steepest descent and L-BFGS minimisation, checkpoint and restart, RMSD, radius of gyration, hydrogen bonds, RDF.
- `tpt-dynamics-compute`: scalar and tpt-simd pair kernels, Ewald real-space terms, PME with a pure-Rust FFT, kernel benchmarks.
- `tptdyn`: `simulate`, `analyze`, `histogram`, completions and man pages, exit-code table, integration tests.
- `bindings/tpt-dynamics-py`: PyO3 bindings (lazy frames, histograms, simulation), NumPy interop, type stubs, pytest suite, and a PyROOT migration guide. Requires ADR-007.
- Embedded demo (`embedded-demo/`, separate workspace) building for `thumbv7em-none-eabihf`.
- Audit and documentation: unsafe-code audit, security review, fuzzing guide, tutorials, benchmark regression workflow.

### Known issues

- None open on licences: `target-lexicon` is replaced by the MIT-licensed `target-lexicon-mit` via `[patch.crates-io]`, and the LLVM-exception allowance is removed (ADR-005 amendment, later 2026-11-16).
