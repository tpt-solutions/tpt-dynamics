# Draft upstream issues for tpt-strata

- **Status:** DRAFTS. Not filed. Filing a GitHub issue publishes text outside the
  project, so each draft needs owner approval first.
- **Source:** gap list in [ADR-004 addendum](../adr/ADR-004-addendum-strata-gap-analysis.md),
  checked against tpt-strata 0.1.0 (MIT OR Apache-2.0, zero dependencies).
- **Scope:** the addendum recommends not adopting tpt-strata as the data engine.
  These drafts are for the owner to decide whether to raise upstream at all. If the
  owner chooses option 1 (keep the own engine, supersede ADR-004), none of them need filing.

---

## 1. Logical plan and expression evaluation

**Title:** Logical plan with computed expressions (`with_column`, arithmetic, function calls)

**Body:** `ProjectExec` takes `Expr { column: usize }`, a column reference only, and
filters compare a column to a constant. Analysis pipelines need computed columns such
as `pt1 * 2`, `invariant_mass(...)`, `delta_r(...)`, plus a logical plan that a
pushdown optimiser can rewrite before execution. Could a logical layer sit above the
physical operators, with an expression type that supports arithmetic and user functions?

## 2. Chunked file scan source

**Title:** File-backed scan source with chunked, bounded-memory reads

**Body:** Every source in tpt-strata is an in-memory `Table`, and CSV/JSON parse from
`&str`. Streaming between operators works through `PhysicalPlan::open()`, but reading
a multi-GB CSV needs a scan that yields batches from disk without loading the file.
Is a file source with a configurable chunk size planned?

## 3. Unsigned integer types

**Title:** Add `UInt64` (and `UInt32`) column types

**Body:** `DataType` has `Boolean`, `Int32`, `Int64`, `Float64`, `Utf8`, `Date`,
`Timestamp`, but no unsigned types. Event IDs and run numbers are unsigned 64-bit.
Would comparisons, arrays and key encoding be extended to cover them?

## 4. Binary columns

**Title:** Add a binary column type

**Body:** Only `Utf8` exists for variable-length data. Detector payloads are raw bytes,
not valid UTF-8. A `Binary` type with offsets and values would cover this.

## 5. List (jagged) columns

**Title:** Nested list columns for per-event collections

**Body:** No nested types exist. Per-event particle lists (variable-length arrays of
momenta and IDs) need an offsets-plus-values representation with operators that work
per list (length, element access, filter inside lists). Is this on the roadmap?

---

Not included: histogram binning operator and quantile/variance aggregates. These are
more specific to the physics layer and would be better as local operators than as
upstream requests.
