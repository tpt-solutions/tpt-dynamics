# Unsafe-code audit

- **Date:** 2026-11-16
- **Scope:** every `.rs` and `.toml` file in the workspace, excluding `target/`
  and `.git/`. Searched with `rg -w unsafe`.
- **Policy:** ADR-006. Every crate forbids `unsafe_code`. An exception needs an
  ADR, a local `#[allow(unsafe_code)]`, and a `// SAFETY:` comment.

## Result

No `unsafe` block, `unsafe fn`, `unsafe impl` or `unsafe trait` exists in
hand-written code. The only matches are documentation comments and lint
configuration.

| Location | Occurrence | Justification |
|---|---|---|
| `Cargo.toml` (root) | `unsafe_code = "forbid"` | Workspace default (ADR-006). |
| `bindings/tpt-dynamics-py/Cargo.toml` | `unsafe_code = "allow"` | Crate-level override for PyO3 and rust-numpy generated glue. Hand-written code stays under `#![deny(unsafe_code)]`. Covered by ADR-007. |
| `bindings/tpt-dynamics-py/src/lib.rs` (lines 12–13) | Doc comment and `#![deny(unsafe_code)]` | Documents the override. `deny` (not `forbid`) lets the generated code compile. |
| `tpt-dynamics-compute/src/lib.rs` (line 22) | Doc comment | States the crate has no `unsafe`. |
| `tpt-dynamics-core/src/lib.rs` (line 5) | Doc comment | States the crate has no `unsafe`. |

## Crate-level attributes

| Crate | Attribute | Status |
|---|---|---|
| tpt-dynamics-core | `#![forbid(unsafe_code)]` | present |
| tpt-dynamics-data | `#![forbid(unsafe_code)]` | present |
| tpt-dynamics-sim | `#![forbid(unsafe_code)]` | present |
| tpt-dynamics-compute | `#![forbid(unsafe_code)]` | present |
| tpt-dynamics-cli | `#![forbid(unsafe_code)]` | present |
| tpt-dynamics-py | `#![deny(unsafe_code)]` | deliberate deviation (see above) |

## Open items

1. **ADR-007 is missing.** `bindings/tpt-dynamics-py/Cargo.toml` and its
   `lib.rs` both cite ADR-007, but `docs/adr/` contains only ADR-001 to ADR-006.
   The exception is not yet reviewed. This must be written before release.
2. **Dependencies.** This audit covers workspace code only. Transitive
   crates such as `pyo3`, `numpy`, `rayon` and `libm` contain `unsafe`. They
   are out of scope here. Re-check with `cargo geiger` (MIT OR Apache-2.0) in a
   later pass.
3. **Generated code.** The PyO3 glue has not been inspected. ADR-007 should
   record how it was reviewed.
