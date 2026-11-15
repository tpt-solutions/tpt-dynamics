# ADR-003: SIMD-accelerated N-body kernels

- **Status:** Accepted (implementation pending; see `todo.md` Phase 3)
- **Date:** 2026-11-16

## Context

Most MD run time is spent evaluating pairwise non-bonded forces. Scalar loops are
not enough to approach GROMACS performance.

## Decision

Write the non-bonded inner loops with `tpt-simd`, with a scalar reference path.
Correctness is tested against the scalar path. No architecture-specific assembly.

## Consequences

- Portable across x86_64 and aarch64 once `tpt-simd` is adopted.
- `tpt-simd` declares a high MSRV (1.95 for v0.1.1); this may force the workspace
  MSRV up. Decide when the dependency is added (see `docs/dependencies.md`).
