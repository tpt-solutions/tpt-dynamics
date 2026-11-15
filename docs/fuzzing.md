# Fuzzing

Two layers exist. Both target the untrusted-input readers:

- `tpt-dynamics-data`: CSV, columnar, serialise/JSON, Arrow-IPC-style, Parquet
  (feature `parquet`), trajectory, row store, the distributed wire decoder,
  the HTTP response parser, and a live worker receiving garbage over TCP.
- `tpt-dynamics-sim`: PDB, GRO, the simulation config TOML, and the
  force-field parameter TOML (`ffparams`).

Section 4 records the long campaign runs.

## 1. Stable-toolchain randomised tests (run in CI)

`tpt-dynamics-data/tests/fuzz_parsers.rs` runs with `cargo test`. It uses the
deterministic `SplitMix64` generator, so every failure reproduces.

Each reader gets four kinds of input:

1. random bytes;
2. random text from a parser-relevant alphabet;
3. every truncation of a valid input;
4. random mutations of a valid input (bit flips, overwrites, insertions,
   deletions, and large values written over 4-byte windows to hit length fields).

The contract is that every input either parses or returns `Err`. No input may
panic, hang, or allocate in proportion to an unchecked length.

Iteration count defaults to 10,000 per generator. Override it with `FUZZ_ITERS`:

```text
cargo test -p tpt-dynamics-data --test fuzz_parsers
FUZZ_ITERS=200000 cargo test -p tpt-dynamics-data --release --test fuzz_parsers
```

## 2. cargo-fuzz (longer campaigns, nightly)

`cargo-fuzz` uses libFuzzer and needs a nightly toolchain. It is not part of the
workspace, and it is not run in CI. Use it for long campaigns on a dedicated machine.

Setup (one time, on a machine with nightly):

```text
rustup toolchain install nightly
cargo install cargo-fuzz
cd tpt-dynamics-data
cargo +nightly fuzz init
```

Then add one target per reader under `tpt-dynamics-data/fuzz/fuzz_targets/`. Each target
is a few lines that call the public API on the input bytes and ignore the result:

```rust
#![no_main]
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let _ = tpt_dynamics_data::serialise::Histogram1D::from_bytes(data);
});
```

Suggested targets: `csv` (drive `CsvReader` over `data`), `columnar`
(`columnar::read_table`), `serialise_hist1d` (`Histogram1D::from_bytes`),
`serialise_json` (`Histogram1D::from_json` on `str::from_utf8(data)`), and
`arrow_ipc` (`arrow_ipc::read_all`).

Run one target:

```text
cargo +nightly fuzz run arrow_ipc -- -max_total_time=3600
```

Crashes are saved under `fuzz/artifacts/<target>/`. Reproduce with
`cargo +nightly fuzz run <target> <artifact-file>`, fix the bug, and add the artifact
as a regression case in `fuzz_parsers.rs` or as a unit test.

`fuzz/` is a separate package with its own `Cargo.toml`. Its dependencies (`libfuzzer-sys`,
MIT OR Apache-2.0) must pass `cargo deny check` before the directory is committed.
Do not add it to the workspace `members` list.

## 3. Sanitizers (optional)

Nightly `-Zsanitizer=address` catches out-of-bounds reads in any code that uses
`unsafe`. The workspace currently has none in hand-written code (see
`docs/audit/unsafe-audit.md`), so sanitizers mainly add value for dependencies.

## 4. Campaign runs (stable toolchain, release build)

Each run uses the deterministic generators, so a failure reproduces from the same
seed. Results from 16 November 2026 on the development machine:

| Target | Command | Iterations | Result |
|---|---|---|---|
| Data-crate parsers | `FUZZ_ITERS=200000 cargo test -p tpt-dynamics-data --release --features parquet --test fuzz_parsers` | 200,000 per generator | 12 passed, no panics |
| Network readers and live worker | `FUZZ_ITERS=20000 cargo test -p tpt-dynamics-data --release --test fuzz_network` | 20,000 | 4 passed, worker kept serving valid jobs |
| Simulation PDB/GRO | `FUZZ_ITERS=200000 cargo test -p tpt-dynamics-sim --release --test fuzz_sim_parsers` | 200,000 | 5 passed |
| Config and force-field TOML | `FUZZ_ITERS=200000 cargo test -p tpt-dynamics-sim --release --test fuzz_config` | 200,000 | 5 passed |

No input caused a panic, hang, or unbounded allocation. The campaign did not
find a defect, so no regression artifacts were added. Coverage-guided runs
(section 2) are still the stronger check and need a nightly toolchain.
