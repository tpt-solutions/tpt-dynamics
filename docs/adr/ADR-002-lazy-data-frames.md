# ADR-002: Lazy data frames for HEP analysis

- **Status:** Accepted (implementation pending; see `todo.md` Phase 1)
- **Date:** 2026-11-16

## Context

LHC-scale datasets are too large to hold in RAM. Eager evaluation materialises
every intermediate result.

## Decision

Analysis is expressed as a logical plan (scan, filter, project, group, aggregate)
that is optimised and then executed in chunks. Results must equal eager
evaluation.

## Consequences

- Memory use is bounded by chunk size, not dataset size.
- Requires a plan optimiser and an executor; more work than an eager API but
  required by the 10 GB milestone.
