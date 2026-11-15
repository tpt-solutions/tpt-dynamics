# ADR-008: Own read-only Parquet reader (import only)

- **Status:** Accepted
- **Date:** 2026-11-16

## Context

Parquet is needed only to import existing data. The Parquet crates
(`parquet`, `arrow`) are Apache-2.0 only and banned by ADR-005, and
`tpt-strata-parquet` depends on them (ADR-004). So the Phase 0 licence decision
closed the Parquet item.

Parquet is a public format with a public specification. Its codec crates are
permissive: `snap` (BSD-3-Clause) and `flate2` (MIT OR Apache-2.0).

## Decision

Implement a clean-room, read-only Parquet reader in `tpt-dynamics-data`, behind
the optional `parquet` feature (off by default):

- **Scope:** flat schemas; `BOOLEAN`, `INT32`, `INT64`, `FLOAT`, `DOUBLE`;
  `PLAIN`, `RLE_DICTIONARY` and `PLAIN_DICTIONARY` encodings; `UNCOMPRESSED`,
  `SNAPPY` and `GZIP` codecs; data page v1.
- **Not supported (explicit errors):** nested and repeated fields, byte-array,
  `INT96` and fixed-length binary columns, data page v2, `ZSTD`/`LZ4`/`BROTLI`,
  and any null value, because `Column` has no null representation.
- **Dependencies:** `snap` and `flate2` (default features off, `rust_backend`).
  No `*-sys` crate and no C code.
- **Writing Parquet** is out of scope.
- **Testing:** fixtures written by pyarrow are committed under
  `tests/fixtures/parquet/`. Tests cover every codec and encoding combination,
  projection, lazy-engine integration, error paths, and a byte-flip corruption
  sweep that must not panic.

## Consequences

- Corrupt input is bounded: chunk reads are checked against the file length,
  page sizes are capped at 1 GiB, and RLE allocations use fallible reserve.
- Each column chunk is decoded in memory, so peak memory is one chunk plus the
  output. Streaming decode is not implemented.
- Unsigned and timestamp logical types are read as raw signed integers.
- Extending the scope (nested types, ZSTD, data page v2) needs a new ADR entry
  and tests against pyarrow fixtures.
