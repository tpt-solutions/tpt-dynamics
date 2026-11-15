# ADR-007: Unsafe-code exception for the Python bindings (`tpt-dynamics-py`)

- **Status:** Accepted
- **Date:** 2026-11-16

## Context

ADR-006 forbids `unsafe` across the workspace (`unsafe_code = "forbid"`). It
anticipated that `pyo3` generates `unsafe` glue and that the Python bindings
would need a documented exception, either per item or by relaxing the lint for
that crate only, recorded in a new ADR. This is that ADR.

The bindings crate `bindings/tpt-dynamics-py` depends on:

- `pyo3` 0.30 (MIT OR Apache-2.0). Its `#[pyclass]`, `#[pymethods]`,
  `#[pyfunction]` and `#[pymodule]` macros expand to `unsafe` code that calls the
  CPython C API (reference counting, GIL token handling, type slots, module
  initialisation).
- `numpy` 0.30 (rust-numpy, BSD-2-Clause). Its array types wrap the NumPy C API
  and generate `unsafe` for buffer access and dtype checks.

Per-item `#[allow(unsafe_code)]` cannot be placed on macro output, since the
attributes would have to be emitted by the macro, so the exception is expressed
at crate level.

## Decision

1. In `bindings/tpt-dynamics-py/Cargo.toml`, `[lints.rust]` sets
   `unsafe_code = "allow"` instead of inheriting the workspace table. All other
   workspace lints (`missing_docs`, `unreachable_pub`, the clippy pedantic set
   and its allow-list) are restated in that crate's manifest unchanged.
2. The crate root sets `#![deny(unsafe_code)]`. Hand-written code in the crate
   therefore cannot contain `unsafe` (a source-level `deny` overrides the
   manifest `allow`). The only `unsafe` the compiler accepts is what the macros
   generate.
3. No other crate in the workspace is affected. `forbid` remains in force
   everywhere else.

## Scope and justification

- The exception covers generated glue only. Hand-written `unsafe` is not
  permitted in this crate; adding any would require a revision of this ADR
  naming the block, and a `// SAFETY:` comment per ADR-006.
- Safe code is insufficient: the Python extension module must export a C
  entry point (`PyInit_*`) and interact with the CPython object model, which
  requires FFI. The Rust-side code that does the physics and data work is
  entirely in the safe crates `tpt-dynamics-data` and `tpt-dynamics-sim`.
- The FFI boundary is the pinned `pyo3` and `numpy` versions. Upgrading them is
  an event that requires re-reading their changelogs for soundness fixes.

## Consequences

- Phase 8 (unsafe-code audit) must include this crate's generated glue by
  reviewing the pinned `pyo3` and `numpy` versions, not the expansion itself.
- The workspace `zero-ffi` CI check (no `*-sys` crate in the normal dependency
  tree) is not affected by this ADR. `pyo3` depends on `pyo3-ffi` (not a `-sys`
  crate name). The check must be re-run after the dependency is added.
- The Python wheel build and tests are not covered by `cargo deny`; they are
  verified by the wheel workflow (`.github/workflows/wheels.yml`).
