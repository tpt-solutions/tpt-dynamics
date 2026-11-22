# Changelog

All notable changes to the `tpt-dynamics-py` bindings (the `tpt_dynamics` Python
package) are documented in this file. The workspace changelog is in the
[root `CHANGELOG.md`](../../CHANGELOG.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0. The
crate is not published to crates.io (`publish = false`); the Python package ships to
PyPI.

## [Unreleased]

### Added

- PyO3 bindings, packaged as an abi3 wheel for Python 3.11 and later, with NumPy as
  the only runtime dependency.
- `Histogram1D` and `Histogram2D`: numpy fill with optional weights, `merge`, bin
  contents and errors, edges, under- and overflow, and dictionary export.
- Hist-style interop with the boost-histogram and `hist` naming conventions.
- Lazy data frames: scan from numpy arrays or CSV, `filter`, `select`,
  `with_column`, `group_by` with `agg`, and `collect` to numpy arrays.
- Four-vector helpers for invariant mass, transverse momentum and ΔR.
- `run_simulation`: molecular dynamics runs configured from a TOML string. It writes
  the `energy_file` and `trajectory_file` named in the `[output]` table, at their
  configured intervals, as `tptdyn simulate` does.
- Python type stubs (`_tpt_dynamics.pyi`) for the extension module.
- A pytest suite (`tests/test_basic.py`, `tests/test_histogram_interop.py`).
- A PyROOT migration guide in `docs/site/migration/pyroot.html`.
- Python exceptions mapped from the Rust error types. See the README for the mapping.

### Notes

- The crate contains no hand-written `unsafe` code (`#![deny(unsafe_code)]`). The
  generated PyO3 and rust-numpy glue is covered by ADR-007.
