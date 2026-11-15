# Tutorial: dilepton mass analysis with the data crate

This tutorial walks through a small HEP analysis: read per-particle rows, pick the
opposite-charge muon pair in each event, histogram the invariant mass, merge
partial histograms, and save the result. The complete program is
[`tpt-dynamics-data/examples/tutorial_hep.rs`](../../tpt-dynamics-data/examples/tutorial_hep.rs).
The snippets below are copied from that file, so they compile with the crate.

Run it:

```text
cargo run --release -p tpt-dynamics-data --example tutorial_hep
```

Expected output (the synthetic sample is seeded, so the numbers are reproducible):

```text
events with a dilepton pair: 10010
entries: 10010, peak bin centre: 91.2 GeV, mean: 91.02 GeV
serialised histogram: 2995 bytes
```

## 1. Input format

The reader expects one CSV row per particle, with a header:

```text
event,pt,eta,phi,mass,charge
```

Rows of one event must be contiguous. The example generator writes them that way.
`synthetic_csv` produces the text in memory. For a real file, pass a `BufReader`
over a `std::fs::File` instead of the `&str`. The reader then holds one line at a
time.

## 2. Streaming read

```rust
let reader = CsvReader::new(BufReader::new(csv.as_bytes()), CsvOptions::default())?;
for record in reader {
    let rec = record?;
    let event = rec.u64(0)?;
    let charge = rec.i64(5)?;
    // ...
}
```

`CsvReader` is an iterator of `Result<CsvRecord, CsvError>`. The typed accessors
(`f64`, `u64`, `i64`) report the line number on failure, so malformed input is easy to
locate. `CsvOptions::default()` uses a comma, a header row, and a 1 MiB line limit.

## 3. Event selection

Each event contributes at most one pair: the first positive and the first negative
lepton. `pair_mass` takes detector coordinates `(pt, eta, phi, m)` and returns the
invariant mass.

```rust
fn push_pair_mass(leptons: &[(i64, (f64, f64, f64, f64))], masses: &mut Vec<f64>) {
    let positive = leptons.iter().find(|(q, _)| *q > 0);
    let negative = leptons.iter().find(|(q, _)| *q < 0);
    if let (Some((_, a)), Some((_, b))) = (positive, negative) {
        masses.push(pair_mass(*a, *b));
    }
}
```

For more than two particles, `physics::invariant_mass(&[FourVector])` sums the
four-vectors. `FourVector::from_pt_eta_phi_m` builds one from detector coordinates.

## 4. Histogramming

```rust
let mut h = Histogram1D::uniform(120, 60.0, 120.0)?;
h.fill(m);
```

Bins are half-open, `[edge[i], edge[i+1])`. Values outside the range go to the
underflow and overflow counters, which are kept apart from the bins.
`fill_weighted(x, w)` takes a weight. `contents()`, `bin_error(i)`, `mean()`,
`std_dev()` and `quantile(p)` give the summaries.

## 5. Merging partial results

Workers fill their own histograms and the results are merged. `merge` requires
identical binning and returns `HistogramError::BinningMismatch` otherwise.

```rust
let mut total = Histogram1D::uniform(120, 60.0, 120.0)?;
total.merge(&second)?;
total.merge(&first)?;
```

Merging is associative and commutative, so the order does not matter. The
distributed-merge test (`tpt-dynamics-data/tests/distributed_merge.rs`) checks this
across threads and serialised partial results. It runs in-process only; no network
is involved.

## 6. Saving and loading

```rust
let bytes = total.to_bytes();
let restored = Histogram1D::from_bytes(&bytes)?;
assert_eq!(restored, total);
```

The binary format is versioned and carries a CRC-32. `from_bytes` rejects truncated,
corrupted, or wrong-type input with a `SerialError`. `to_json()` and `from_json()`
give a readable form with the same schema checks.

## Related modules and current limits

- `tpt_dynamics_data::jagged` holds variable-length per-event collections. It can
  replace the per-row loop above with per-event slices, reductions and filters.
- `tpt_dynamics_data::arrow_ipc` reads a documented Arrow-IPC-style format. It is
  **not** Arrow-compatible, so files written by Arrow libraries will not load.
- ROOT and RNTuple files cannot be read yet. Input is CSV, the project columnar format,
  or the Arrow-IPC-style format above.
- The lazy data frame (`tpt_dynamics_data::lazy`) exists. Filling a histogram directly
  from a lazy frame aggregation is not yet connected.
