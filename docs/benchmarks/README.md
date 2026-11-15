# Performance regression tracking

CI job: `.github/workflows/bench.yml`.

- The probe is `tpt-dynamics-data/examples/bench_json.rs`. It times four hot
  paths (histogram fill, histogram serialise round trip, CSV parse, jagged
  sums). Each is the median of five runs.
- `docs/benchmarks/baseline.json` holds the reference timings in nanoseconds.
- A metric fails when it is more than 25% slower than the baseline. Faster
  results always pass.

## Running locally

```text
cargo run --release -p tpt-dynamics-data --example bench_json -- --check docs/benchmarks/baseline.json
cargo run --release -p tpt-dynamics-data --example bench_json -- --out docs/benchmarks/baseline.json   # refresh
```

## Caveats

- The baseline was recorded on one Windows development machine. It is not
  comparable to GitHub-hosted runners. Before relying on the CI gate, record a
  baseline from the CI runner (download the `bench-current` artifact from a
  green run and commit it as `baseline.json`).
- The histogram-fill metric varied by about 35% between two local runs. The
  25% threshold can therefore flag noise. If this happens, re-run first.
- Only the data crate is covered. Kernel micro-benchmarks in
  `tpt-dynamics-compute` are tracked separately (`docs/benchmarks/compute-kernels.md`).
