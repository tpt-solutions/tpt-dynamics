# Changelog

All notable changes to `tpt-dynamics-data` are documented in this file. The workspace
changelog is in the [root `CHANGELOG.md`](../CHANGELOG.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
crate uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0.

## [Unreleased]

### Added

- `hist1d`: `Histogram1D` with uniform and variable binning, under- and overflow,
  weights, sum of squared weights, bin errors, mean, standard deviation, quantiles,
  integral, and a commutative, associative `merge`.
- `hist2d`: `Histogram2D` with the same weighting and merge semantics.
- `profile`: `Profile1D` and `Profile2D`.
- `serialise`: binary (CRC-32 checksummed) and JSON round-trips for every histogram
  type.
- `json`: dependency-free JSON writer and reader.
- `csv`: `CsvReader`, a streaming reader with bounded memory (`max_line_bytes`),
  batched output via `next_batch`, and line-numbered parse errors.
- `columnar`: a versioned, checksummed binary columnar table, with `write_table`,
  `read_table` and `ColumnarWriter`.
- `rowstore`: a fixed-width binary column file with one seek per column per range,
  implementing `lazy::Source`.
- `http`: `RowStore::open_http` reads row stores over plain HTTP `Range` requests.
  No TLS.
- `trajectory`: a compressed, chunked trajectory format used for simulation output.
- `jagged`: variable-length per-event collections.
- `lazy`: `LazyFrame` with predicate and projection pushdown, chunked and parallel
  execution, `group_by` and `agg`, `fill_histogram` for filling a `Histogram1D` from a
  plan, and `lazy::distributed`, which runs plans over loopback TCP workers.
- `arrow_ipc`: read-only reader for an Arrow-IPC-style format. It is not
  Arrow-compatible.
- `physics`: `invariant_mass`, `pair_mass`, `transverse_momentum`, `delta_phi` and
  `delta_r`.
- `parquet` feature (off by default): read-only Parquet import with an own
  clean-room reader for flat schemas, with SNAPPY, GZIP and uncompressed codecs and
  dictionary encoding (ADR-008).
- `mmap` feature (off by default): zero-copy memory-mapped reader for the columnar
  format (ADR-011).
- Criterion benchmark for 1D histogram filling (`cargo bench -p tpt-dynamics-data`).
- Examples: `generate_events`, `generate_trajectory`, `tutorial_hep` and `bench_json`.

### Security

- The distributed execution protocol has no authentication. Run workers only on a
  trusted network.

### Known issues

- ROOT and RNTuple files cannot be read.
- Writing Parquet is not planned.
