# ADR-004 addendum: tpt-strata gap analysis

- **Status:** Analysis only. ADR-004 is unchanged. No dependency added.
- **Date:** 2026-11-16
- **Checked against:** `tpt-strata` 0.1.0 (MIT OR Apache-2.0, zero dependencies),
  read from the crates.io tarball. Nothing was built or executed from it.

## Question

ADR-004 accepts building the data engine on `tpt-strata` and says gaps are to be
filled "in the crate or upstream". This addendum checks which needs `tpt-strata`
covers today, and which it does not.

## What tpt-strata 0.1.0 provides

| Area | Surface |
|------|---------|
| Column types | `Boolean`, `Int32`, `Int64`, `Float64`, `Utf8`, `Date` (i32 days), `Timestamp` (i64 us). Nullable via validity. |
| In-memory tables | `Table`, `Column`, `Batch`, typed arrays (`F64Array`, `I64Array`, ...). |
| Physical operators | `FilterExec`, `ProjectExec`, `JoinExec`, `SortExec`, `LimitExec`, `CondFilterExec`, aggregation. |
| Streaming | `PhysicalPlan::open()` returns a pull-based `BatchStream`. Operators that cannot stream fall back to materialising through `execute()`. |
| Explain | Each operator implements `explain()` (indented tree, root first). |
| Parallelism | Optional `parallel` feature (internal pool). |
| Query surface | `QueryBuilder` (filter by column/op/value, project, group_by, aggregate, order_by, limit), a SQL parser, and CSV/JSON parsing from `&str`. |
| Aggregates | `Sum`, `Count`, `Min`, `Max`, `Avg`. |

## Gaps against the data-engine needs

| Need (source) | tpt-strata 0.1.0 | Gap |
|---------------|------------------|-----|
| Lazy frames: logical plan, optimiser, `collect` (todo Phase 1, ADR-002) | Physical operators only. No logical plan, no optimiser. | **Large.** Predicate and projection pushdown must be written on top. The crate's own pruning covers only some filter cases. |
| Expression engine for physics quantities (invariant mass, pT, dR) (todo Phase 1) | `ProjectExec` takes `Expr { column: usize }`: a column reference only. Filters compare a column to a constant. | **Large.** No computed expressions, no `with_column`. |
| Chunked, bounded-memory streaming from disk (todo Phase 1, 10 GB CSV milestone) | `open()` streams between operators, but every source is an in-memory `Table`. CSV and JSON parse from `&str`. No file scan, no chunked reader. | **Large.** A file source with chunked reads is needed. |
| UINT64 (ADR/spec) | No unsigned types. `DataType` has no `UInt64` or `UInt32`. | **Medium.** Adding a variant touches the type enum, arrays, comparisons, and key encoding. |
| Binary columns (spec, detector payloads) | No binary type. `Utf8` only. | **Medium.** |
| List / jagged columns, per-event particles (todo Phase 1 "jagged arrays") | No nested types. | **Large.** No offsets-plus-values representation. The data crate already has `jagged.rs`. |
| N-tuples, multi-column event records | Column-oriented tables cover the flat case. Nested records are not covered. | Follows from list support. |
| Histogram fill from lazy frames (todo Phase 1, unchecked) | No binning operator. | Medium. Can be a custom operator. |
| Quantiles, variance, sum of squared weights (histogram stats) | No such aggregates. | Small to medium. |
| Parquet (todo Phase 0 decision) | Not in this crate. `tpt-strata-parquet` is rejected (Arrow, Apache-only). | Out of scope. |

## Existing overlap in tpt-dynamics-data

`tpt-dynamics-data` already has its own lazy module (`src/lazy/`, with `Plan`,
`Expr`, an optimiser, and a chunked executor with a window and rayon). It also has
its own columnar format (`columnar.rs`), an Arrow-IPC-style reader
(`arrow_ipc.rs`), jagged arrays (`jagged.rs`), and a trajectory reader. None of
these depends on tpt-strata. The crate's `DType` has `F64`, `I64`, and `Bool`;
columnar adds `U32`.

## Conclusion

The gap analysis does **not** support adopting tpt-strata as the data engine. It
covers the in-memory relational core (filter, project, join, sort, limit, group-by
aggregate) and has a streaming operator trait. The needs that define the data
engine (logical lazy plan with expressions, chunked file streaming, UINT64, binary,
lists/jagged) are missing. Adopting it would mean building most of those on top of
it while the crate already owns the same layers. Because the data crate is
already self-contained, `tpt-strata` is not added as a dependency.

ADR-004's premise ("tpt-strata covers the core, gaps are small") does not hold for
the physics-layer needs. The owner should decide whether to:

1. Keep the own engine and mark ADR-004 as superseded, or
2. Keep ADR-004 and accept that the own `lazy/` module is the implementation, with
   tpt-strata used only for SQL-style ad hoc queries if that is wanted later.

Upstream issues and PRs on tpt-strata were not opened. The user did not authorise
that.
