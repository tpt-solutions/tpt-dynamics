# Changelog

All notable changes to `tpt-dynamics-cli` (the `tptdyn` binary) are documented in this
file. The workspace changelog is in the [root `CHANGELOG.md`](../CHANGELOG.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
crate uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0.

## [Unreleased]

### Added

- `tptdyn` binary with `--help` and `--version`, global `--log-level`, `-q/--quiet`,
  `--threads` and a `--config` file for CLI defaults. Unknown keys are rejected.
- `tptdyn init TEMPLATE [DIR]`: writes a project from the `lj-fluid`, `hep-analysis` or
  `water-box` template. The `water-box` template is a scaffold without a force field.
  Nothing is overwritten without `--force`.
- `tptdyn validate CONFIG.toml`: parses a simulation configuration and prints the atom
  count, box, step count, expected energy-log rows, estimated trajectory size and an
  estimate of non-bonded pair evaluations per step and over the run. Unknown keys get
  a "did you mean" hint.
- `tptdyn simulate CONFIG.toml [--output-dir DIR] [--report]`: runs the MD engine from a
  TOML configuration and a `[system]` table. `--report` writes an HTML report of the run.
- `tptdyn report`: builds an HTML report from a finished run.
- `tptdyn analyze`: summaries, filters, expressions and histograms of CSV, Parquet and
  `tptdyn` binary tables. CSV is streamed in chunks of 1,048,576 rows, so memory does
  not grow with file size; a 10.49 GB file peaks at 554 MB. Column statistics are
  running (Welford) values.
- `tptdyn histogram INPUT.csv --column NAME|INDEX --bins N --range LO,HI --out FILE`:
  streams one column into a uniform histogram, written as JSON or binary.
- `tptdyn completions SHELL` (bash, zsh, fish, PowerShell, elvish) and `tptdyn man`.
- Exit codes: 0 success, 1 runtime error, 2 usage error, 3 input or configuration
  error. Errors go to stderr as `tptdyn: error: <message>`.
- A progress display for long runs, which `--quiet` suppresses.
- cargo-binstall metadata, so prebuilt archives from the GitHub release can be
  installed with `cargo binstall tptdyn`.
- Integration tests in `tests/cli.rs` for help, exit codes, histogram bin counts,
  analyze filters and expressions, binary inputs, init, validate, report, completions,
  man pages and short simulations.
