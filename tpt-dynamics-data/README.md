# tpt-dynamics-data

Data analysis for tpt-dynamics: the ROOT replacement. It provides physics histograms,
lazy data frames, file readers and writers, and a distributed-merge path. Everything
is pure Rust with no C++ or FFI dependencies.

Modules:

- `hist1d`, `hist2d`, `profile`: 1D and 2D histograms and profiles with fixed or
  variable binning, under- and overflow, weights, statistics, quantiles, and a
  commutative `merge` for distributed filling.
- `serialise`, `json`: binary (checksummed) and JSON round-trips for histograms.
- `csv`: streaming CSV reader with bounded memory and line-numbered errors.
- `columnar`: versioned, checksummed binary columnar format.
- `rowstore`, `trajectory`: row-oriented storage and a compressed chunked trajectory
  format.
- `jagged`: variable-length per-event collections.
- `lazy`: lazy data frames with predicate and projection pushdown, chunked parallel
  execution, histogram filling from a plan (`fill_histogram`), and distributed
  execution of plans over loopback TCP (no authentication; trusted networks only).
- `arrow_ipc`: read-only reader for an Arrow-IPC-style format. It is not
  Arrow-compatible.
- `physics`: invariant mass, transverse momentum and delta R helpers.
- `http`: reads row stores over plain HTTP `Range` requests (no TLS).
- `rng`, `crc32`: small building blocks used by the modules above.

Optional features:

- `parquet`: read-only Parquet import, a clean-room reader with no Arrow or Parquet
  crates (ADR-008). Writing Parquet is not planned.
- `mmap`: zero-copy memory-mapped reader for the columnar format (ADR-011).

ROOT and RNTuple files cannot be read.

Part of [tpt-dynamics](../README.md). Licensed under MIT OR Apache-2.0.
