# ADR-006: Unsafe-code policy

- **Status:** Accepted
- **Date:** 2026-11-16

## Decision

- Every crate has `#![forbid(unsafe_code)]`, enforced workspace-wide by the
  `unsafe_code = "forbid"` lint.
- An exception requires: (a) a new ADR that names the crate, the block, and why
  safe code is insufficient; (b) a local `#[allow(unsafe_code)]` on the item; and
  (c) a `// SAFETY:` comment giving the invariant that makes it sound.
- Exceptions are reviewed again in Phase 8 (unsafe-code audit).

## Consequences

- `pyo3` generates `unsafe` glue. When the Python bindings are added, the crate
  will need a documented exception (or the lint is relaxed for that crate only,
  via a new ADR).
