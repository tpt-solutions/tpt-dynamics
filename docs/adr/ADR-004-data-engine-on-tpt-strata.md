# ADR-004: Build the data engine on tpt-strata (no polars or arrow)

- **Status:** Accepted (adoption pending; see `todo.md` Phase 1)
- **Date:** 2026-11-16

## Context

`polars` and `arrow` are candidate foundations for the data engine. `arrow` is
Apache-2.0 only and is banned by ADR-005. `polars` is a large dependency tree.
`tpt-strata` is a dependency-free crate under the owner's control.

## Decision

Build the data engine on `tpt-strata` (`tpt-strata` 0.1.0 has no dependencies).
`tpt-strata-parquet` is **not** adopted: it depends on `arrow` and `parquet`, which
are Apache-2.0 only. Parquet support, if needed, must come from an own reader or
be placed behind a feature that is excluded from the default build and
documented in ADR-005.

## Consequences

- Full control of the columnar model, including jagged arrays for per-event
  particles.
- Gaps in `tpt-strata` (lazy execution, streaming, UINT64, binary and list columns)
  must be filled in the crate or upstream.
