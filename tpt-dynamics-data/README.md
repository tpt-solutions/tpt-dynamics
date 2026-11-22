# tpt-dynamics-data

Data analysis for [tpt-dynamics](../README.md): the ROOT replacement. It provides
physics histograms, lazy data frames, streaming file readers, a checksummed binary
columnar format and a distributed-merge path.

Everything is pure Rust. There are no C++ libraries, no ROOT and no FFI, and the
crate contains no `unsafe` code outside one documented memory-mapping call behind
the `mmap` feature (ADR-011).

## Contents

| Module | What it provides |
| --- | --- |
| `hist1d` | `Histogram1D`: fixed or variable binning, under- and overflow, weights, sum of squared weights, bin errors, mean, standard deviation, quantiles, integral and a commutative `merge` for distributed filling. |
| `hist2d` | `Histogram2D`: two-dimensional histograms with the same weighting and merge semantics. |
| `profile` | `Profile1D` and `Profile2D`: the mean of a value in bins of one or two variables, with per-bin spread. |
| `serialise` | Binary (checksummed) and JSON round-trips for every histogram type. |
| `json` | Dependency-free JSON writer and reader used by `serialise`. |
| `csv` | `CsvReader`: a streaming CSV reader with bounded memory (`max_line_bytes`) and line-numbered errors. Yields batches with `next_batch`. |
| `columnar` | A versioned, checksummed binary columnar table (`write_table`, `read_table`, `ColumnarWriter`). |
| `rowstore` | A fixed-width binary column file that serves any row range with one seek per column (`write_rowstore`, `RowStore::open`). |
| `http` | Reads row stores over plain HTTP `Range` requests (`RowStore::open_http`). No TLS. |
| `trajectory` | A compressed, chunked trajectory format used by `tpt-dynamics-sim` output. |
| `jagged` | Variable-length per-event collections (jagged arrays). |
| `lazy` | Lazy data frames: `LazyFrame` plans with predicate and projection pushdown, chunked parallel execution, `group_by` aggregation, histogram filling from a plan (`fill_histogram`), and distributed execution over loopback TCP (`lazy::distributed`). |
| `arrow_ipc` | Read-only reader for an Arrow-IPC-style format. It is **not** Arrow-compatible. |
| `physics` | Invariant mass, transverse momentum and ΔR helpers for four-vectors. |
| `mapped` | Zero-copy memory-mapped reader for the columnar format (feature `mmap`). |
| `parquet` | Read-only Parquet import (feature `parquet`). |
| `rng`, `crc32` | Small building blocks used by the modules above. |

## Quick start

Filling and summarising a histogram:

```rust
use tpt_dynamics_data::hist1d::Histogram1D;

// 50 uniform bins over [0, 200).
let mut h = Histogram1D::uniform(50, 0.0, 200.0).expect("valid binning");

for &mass in &[88.1, 91.2, 91.4, 125.0] {
    h.fill(mass);
}
h.fill_weighted(91.3, 0.5);

assert_eq!(h.nbins(), 50);
assert_eq!(h.entries(), 5);
// Bin 22 covers [88, 92): 88.1, 91.2, 91.4 and the weighted 91.3 (weight 0.5).
assert!((h.bin(22) - 3.5).abs() < 1e-12);
println!("mean = {:?}", h.mean());
println!("median = {:?}", h.quantile(0.5));
```

Invariant mass of two particles in detector coordinates `(pt, eta, phi, m)`, and
transverse momentum from its components:

```rust
use tpt_dynamics_data::physics::{pair_mass, transverse_momentum};

// Two massless 50 GeV particles, back to back in phi at eta = 0: m = 2 * pt.
let mass = pair_mass((50.0, 0.0, 0.0, 0.0), (50.0, 0.0, std::f64::consts::PI, 0.0));
assert!((mass - 100.0).abs() < 1e-9);
assert!((transverse_momentum(3.0, 4.0) - 5.0).abs() < 1e-12);
```

Lazy plans fill a histogram directly. Only the selected columns are read, and
filters are applied before the fill:

```rust
use tpt_dynamics_data::hist1d::Histogram1D;
use tpt_dynamics_data::lazy::{col, Column, DataFrame, LazyFrame};

let df = DataFrame::new(vec![
    ("pt", Column::F64(vec![5.0, 15.0, 25.0, 35.0])),
    ("keep", Column::Bool(vec![true, true, false, true])),
])
.unwrap();

let mut h = Histogram1D::uniform(4, 0.0, 40.0).unwrap();
LazyFrame::scan(df)
    .filter(col("keep"))
    .fill_histogram(&mut h, col("pt"), None)
    .unwrap();
assert_eq!(h.integral(), 3.0);
```

Any type implementing `lazy::Source` can be scanned. `DataFrame` and `RowStore`
(including row stores opened over HTTP) do. The CSV reader is used directly, as
`tptdyn analyze` does, and `examples/` shows complete runs.

Distributed merging needs no special code for histograms: each worker fills its own
`Histogram1D` and the coordinator calls `merge`. The operation is associative and
commutative, so the result does not depend on the order in which parts arrive.

## Features

| Feature | Default | Effect |
| --- | --- | --- |
| `parquet` | no | Read-only Parquet import. A clean-room reader for flat schemas with SNAPPY, GZIP or uncompressed pages and dictionary encoding. Adds the `snap` (BSD-3-Clause) and `flate2` (MIT OR Apache-2.0) dependencies. Writing Parquet is not planned (ADR-008). |
| `mmap` | no | Zero-copy memory-mapped reader for the columnar format. Adds `memmap2`. The mapping call is the only `unsafe` block in the crate (ADR-011). |

## Formats and limits

- **Binary histograms and columnar tables are checksummed** (CRC-32). A corrupted
  file is reported as an error rather than read silently.
- **JSON histograms** round-trip through `serialise`, and are the human-readable
  option used by the CLI.
- **CSV is read, not written by this crate.** The reader streams, so memory use does
  not grow with file size; `tptdyn analyze` processes multi-gigabyte files in
  1,048,576-row chunks (see `docs/benchmarks/large-data.md`).
- **Not supported**: ROOT and RNTuple files, Arrow-compatible IPC, Parquet writing,
  and TLS on the HTTP reader.
- **Networking is trusted-network only.** The distributed execution protocol has no
  authentication. Run workers only on a network you control.

## Benchmarks

The criterion benchmark for 1D histogram filling runs with:

```text
cargo bench -p tpt-dynamics-data
```

Measurements are recorded in `docs/benchmarks/`.

## Licence

Licensed under either of [MIT](../LICENSE-MIT) or [Apache-2.0](../LICENSE-APACHE), at
your option. Part of [tpt-dynamics](../README.md).
