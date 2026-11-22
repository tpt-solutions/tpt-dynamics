# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
the project uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0.

## [Unreleased]

### Added

- `tptdyn validate` prints an estimate of non-bonded pair evaluations per step
  and over the run, from the uniform-density cutoff-sphere count.
- Python `run_simulation` writes the `energy_file` and `trajectory_file` named in
  the TOML `[output]` table, at their configured intervals, as `tptdyn simulate` does.
- `[system]` molecule grid for `tptdyn simulate` and `validate`: copies of one
  ffparams molecule template on a cubic grid, so bonded force fields (with PME and
  SETTLE) run from a TOML file. Requires `[forcefield] kind = "bonded"`.
- `tptdyn init water-box`: a scaffold of 512 rigid water molecules. Its parameter file
  is supplied by the user, because no parameter set is bundled.
- Security policy (`SECURITY.md`) and Code of Conduct (`CODE_OF_CONDUCT.md`).
- Cargo workspace with six members: `tpt-dynamics-core`, `-data`, `-sim`,
  `-compute`, `-cli` and `bindings/tpt-dynamics-py`.
- `tpt-dynamics-core` (`no_std`): CODATA 2018 physical constants, scalar unit
  conversions (eV, Å, nm, kJ/mol, kcal/mol, amu, K), typed quantities
  (`Energy`, `Length`, `Mass`, `Temperature`) and a `Vec3` type.
- `tpt-dynamics-data`: `Histogram1D` with fixed and variable binning, under- and
  overflow, weights, sum of squared weights, mean, standard deviation, quantiles,
  and merging.
- `tptdyn` binary with `--help`, `--version`, `simulate`, `analyze`, `histogram`,
  shell completions and man pages.
- Licence files (MIT and Apache-2.0), `deny.toml` with the project licence policy,
  and the dependency audit in `docs/dependencies.md`.
- Architecture decision records ADR-001 to ADR-006 in `docs/adr/`.
- Continuous integration: formatting, clippy, tests on three operating systems,
  documentation build, `cargo deny`, `no_std` build, and an MSRV job.

### Added (Phase 9 adoption and fixes, unreleased)

- `release.yml`: tag-triggered builds of `tptdyn` for four targets, packaged with checksums into a draft GitHub Release. The crates.io and PyPI publish jobs are manual and approval-gated.
- `install.sh` and `install.ps1` (verify SHA-256 before installing), and draft Homebrew and winget manifests in `packaging/`. The installers work after the first release.
- Root `README.md` rewritten, with a quickstart and `examples/` (a small LJ run, a CSV generator and a Python histogram script).
- `.devcontainer/`, a multi-stage `Dockerfile` and `.dockerignore`. CI builds `embedded-demo` and the Docker image.
- `docs-site.yml`: checks that `docs/site` is up to date and deploys it to GitHub Pages from `master`.
- Criterion benchmark harness for `tpt-dynamics-data` (`cargo bench -p tpt-dynamics-data`), covering 1D histogram filling.
- `SECURITY.md` and `CODE_OF_CONDUCT.md` added. The debug `.pdb` file is no longer tracked.

### Added (Phases 1–8 implementation pass, unreleased)

- `tpt-dynamics-data`: 2D, profile histograms; JSON and checksummed binary serialisation; streaming CSV reader; checksummed columnar format; physics helpers; a lazy data frame with predicate and projection pushdown and chunked, parallel execution; jagged arrays; a read-only Arrow-IPC-style reader (not Arrow-compatible); compressed chunked trajectory format; distributed-merge and fuzz-style tests.
- `tpt-dynamics-sim`: PDB/GRO I/O, TOML config, cell and Verlet lists, Lennard-Jones, bonded terms, 1-4 and exclusions, SHAKE, SETTLE, Nosé–Hoover, v-rescale, Berendsen, C-rescale, triclinic boxes, solvation and ion placement, steepest descent and L-BFGS minimisation, checkpoint and restart, RMSD, radius of gyration, hydrogen bonds, RDF.
- `tpt-dynamics-compute`: scalar and tpt-simd pair kernels, Ewald real-space terms, PME with a pure-Rust FFT, kernel benchmarks.
- `tptdyn`: `simulate`, `analyze`, `histogram`, completions and man pages, exit-code table, integration tests.
- `bindings/tpt-dynamics-py`: PyO3 bindings (lazy frames, histograms, simulation), NumPy interop, type stubs, pytest suite, and a PyROOT migration guide. Requires ADR-007.
- Embedded demo (`embedded-demo/`, separate workspace) building for `thumbv7em-none-eabihf`.
- Audit and documentation: unsafe-code audit, security review, fuzzing guide, tutorials, benchmark regression workflow.

### Added (checklist reconciliation, unreleased)

- `SECURITY.md` (private reporting, scope), `CODE_OF_CONDUCT.md`, and a pull request template.
- `deny.toml` now denies `cc`; the `*-sys` rule is enforced by the CI `zero-ffi` job.
- `LazyFrame::fill_histogram`: fills a `Histogram1D` from a lazy plan, optionally weighted.
- `parquet` feature (off by default): read-only Parquet import with an own reader, for flat schemas, SNAPPY/GZIP/uncompressed codecs and dictionary encoding (ADR-008).

### Added (large data and distributed execution, unreleased)

- `tptdyn analyze` streams the CSV in chunks of 1,048,576 rows. Memory no longer grows with the file: a 10.49 GB file peaks at 554 MB. Column statistics are running (Welford) values; a column that is integer in some chunks and fractional in others is reported as `f64`.
- `lazy::distributed`: a coordinator splits a source's rows across TCP workers, which run the row-wise stages and return partial aggregates or rows. Results are merged in worker order. Plans are `DistributedPlan` (filter, select, with-column, group-by aggregate). The protocol has no authentication, so run workers only on a trusted network.
- `rowstore`: a fixed-width binary column file that serves any row range with one seek per column. It implements `lazy::Source` and can be read over plain HTTP `Range` requests (`RowStore::open_http`, no TLS).
- `examples/generate_trajectory.rs`: writes a synthetic trajectory one frame at a time.
- Measurements and reproduction commands are in `docs/benchmarks/large-data.md`.

### Known issues

- None open on licences: `target-lexicon` is replaced by the MIT-licensed `target-lexicon-mit` via `[patch.crates-io]`, and the LLVM-exception allowance is removed (ADR-005 amendment, later 2026-11-16).
