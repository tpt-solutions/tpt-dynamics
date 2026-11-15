# ADR-009: Keep the own data engine; supersede ADR-004

- **Status:** Accepted
- **Date:** 2026-11-16
- **Supersedes:** ADR-004 (build the data engine on tpt-strata)
- **Related:** ADR-004 addendum (gap analysis), ADR-002 (lazy data frames), ADR-008 (own Parquet reader)

## Context

ADR-004 chose tpt-strata as the data-engine foundation, on the premise that its
gaps were small. The gap analysis in the ADR-004 addendum found that premise does
not hold. The needs that define the engine are missing from tpt-strata 0.1.0:
a logical lazy plan with expressions, chunked file streaming, unsigned integers,
binary columns, and list (jagged) columns.

By the time of this decision, `tpt-dynamics-data` already implements those layers
itself, and none of them depends on tpt-strata:

- lazy plans with expressions, predicate and projection pushdown, constant
  folding, a chunked and parallel executor, and histogram fill (`src/lazy/`);
- chunked CSV reading (`csv.rs`), including parallel batch parsing;
- the columnar and row-store formats, and an own read-only Parquet reader (ADR-008);
- jagged arrays (`jagged.rs`) and a trajectory reader with bounded memory;
- distributed execution of lazy plans over TCP, and HTTP range reads.

## Decision

1. The data engine remains the in-house `tpt-dynamics-data` crate. It does not
   depend on tpt-strata, and tpt-strata is not added as a dependency.
2. ADR-004 is superseded by this record.
3. If ad-hoc SQL over tables is wanted later, tpt-strata may be added as an
   optional feature for that purpose only. That would need its own ADR.
4. The drafted upstream issues in `docs/upstream/tpt-strata-issues.md` are not
   filed. The owner decides separately whether to raise them.

## Consequences

- No new dependency and no licence review for tpt-strata.
- The engine's behaviour is covered by the project's own tests, including the
  fuzz and distributed tests.
- Any future gap in the engine is fixed here, not upstream.
- The cost is maintaining the data layers in this repo. That cost already exists,
  because the layers are implemented.
