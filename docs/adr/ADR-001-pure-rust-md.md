# ADR-001: Pure-Rust molecular dynamics over FFI to GROMACS

- **Status:** Accepted
- **Date:** 2026-11-16

## Context

GROMACS is GPL/LGPL C++. Linking it would contaminate the workspace licence and
make the engine hard to embed or swap components in.

## Decision

Implement classical MD integrators and force fields in Rust in `tpt-dynamics-sim`.
Do not wrap GROMACS or any other C/C++ MD engine.

## Consequences

- Full control over the licence and the API; components can be swapped at runtime
  through traits.
- Performance must be earned with SIMD (ADR-003) and parallelism (rayon), and
  validated against published reference data.
