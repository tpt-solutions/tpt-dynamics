//! Timing probe for performance-regression tracking in CI.
//!
//! Measures a few hot paths of `tpt-dynamics-data`, takes the median of several
//! runs for each, and writes the results as JSON. With `--check`, it compares
//! the results against a stored baseline and exits with status 1 if any metric
//! is more than [`MAX_REGRESSION`] slower.
//!
//! ```text
//! cargo run --release -p tpt-dynamics-data --example bench_json -- --out current.json
//! cargo run --release -p tpt-dynamics-data --example bench_json -- --check docs/benchmarks/baseline.json
//! ```
//!
//! Timings are wall-clock and machine-dependent. Compare only against a baseline
//! recorded on the same runner class. See `docs/benchmarks/README.md`.

use std::fmt::Write as _;
use std::hint::black_box;
use std::process::ExitCode;
use std::time::Instant;

use tpt_dynamics_data::csv::{CsvOptions, CsvReader};
use tpt_dynamics_data::hist1d::Histogram1D;
use tpt_dynamics_data::jagged::Jagged;
use tpt_dynamics_data::json::{self, Value};
use tpt_dynamics_data::rng::SplitMix64;

/// A metric fails the check when it is more than this fraction slower.
const MAX_REGRESSION: f64 = 0.25;
const RUNS: usize = 5;

fn median_ns(mut f: impl FnMut()) -> f64 {
    let mut samples: Vec<f64> = (0..RUNS)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_nanos() as f64
        })
        .collect();
    samples.sort_by(f64::total_cmp);
    samples[RUNS / 2]
}

fn bench_hist_fill() -> f64 {
    let mut rng = SplitMix64::new(1);
    let xs: Vec<f64> = (0..2_000_000).map(|_| rng.gaussian()).collect();
    let ns = median_ns(|| {
        let mut h = Histogram1D::uniform(200, -5.0, 5.0).unwrap();
        for &x in &xs {
            h.fill(x);
        }
        black_box(&h);
    });
    ns / xs.len() as f64
}

fn bench_serialise() -> f64 {
    let mut h = Histogram1D::uniform(1000, -5.0, 5.0).unwrap();
    for i in 0..10_000 {
        h.fill(f64::from(i) * 0.001 - 5.0);
    }
    let ns = median_ns(|| {
        for _ in 0..100 {
            let bytes = h.to_bytes();
            black_box(Histogram1D::from_bytes(&bytes).unwrap());
        }
    });
    ns / 100.0
}

fn bench_csv() -> f64 {
    let mut text = String::from("a,b,c\n");
    let mut rng = SplitMix64::new(2);
    let nrows: i32 = 200_000;
    for _ in 0..nrows {
        let _ = writeln!(text, "{:.4},{:.4},{:.4}", rng.gaussian(), rng.uniform(0.0, 1.0), 3.0);
    }
    let ns = median_ns(|| {
        let reader = CsvReader::new(text.as_bytes(), CsvOptions::default()).unwrap();
        let mut sum = 0.0;
        for rec in reader {
            sum += rec.unwrap().f64(0).unwrap();
        }
        black_box(sum);
    });
    ns / f64::from(nrows)
}

fn bench_jagged() -> f64 {
    let mut rng = SplitMix64::new(3);
    let mut j: Jagged<f64> = Jagged::new();
    let nevents: i32 = 1_000_000;
    for _ in 0..nevents {
        let n = rng.below(6);
        j.push_event((0..n).map(|_| rng.uniform(0.0, 100.0)));
    }
    let ns = median_ns(|| {
        black_box(j.sums());
    });
    ns / f64::from(nevents)
}

fn run_all() -> Vec<(&'static str, f64)> {
    vec![
        ("hist1d_fill_ns_per_entry", bench_hist_fill()),
        ("serialise_roundtrip_ns_per_hist", bench_serialise()),
        ("csv_parse_ns_per_record", bench_csv()),
        ("jagged_sum_ns_per_event", bench_jagged()),
    ]
}

fn to_json(results: &[(&str, f64)]) -> String {
    let mut out = String::from("{\n  \"metrics\": {\n");
    for (i, (name, value)) in results.iter().enumerate() {
        let comma = if i + 1 < results.len() { "," } else { "" };
        let _ = writeln!(out, "    \"{name}\": {value:.3}{comma}");
    }
    out.push_str("  }\n}\n");
    out
}

fn read_baseline(path: &str) -> Result<Vec<(String, f64)>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("reading {path}: {e}"))?;
    let root = json::parse(&text).map_err(|e| format!("parsing {path}: {e}"))?;
    let Value::Object(top) = root else {
        return Err(format!("{path}: top level must be an object"));
    };
    let metrics = top
        .into_iter()
        .find(|(k, _)| k == "metrics")
        .map(|(_, v)| v)
        .ok_or_else(|| format!("{path}: missing \"metrics\""))?;
    let Value::Object(fields) = metrics else {
        return Err(format!("{path}: \"metrics\" must be an object"));
    };
    fields
        .into_iter()
        .map(|(k, v)| match v {
            Value::Number(n) => Ok((k, n)),
            _ => Err(format!("{path}: metric {k} is not a number")),
        })
        .collect()
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut out_path = None;
    let mut check_path = None;
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--out" => out_path = it.next().cloned(),
            "--check" => check_path = it.next().cloned(),
            other => {
                eprintln!("unknown argument {other}");
                return ExitCode::from(2);
            }
        }
    }

    let results = run_all();
    let json_text = to_json(&results);
    print!("{json_text}");
    if let Some(path) = out_path {
        if let Err(e) = std::fs::write(&path, &json_text) {
            eprintln!("writing {path}: {e}");
            return ExitCode::from(2);
        }
    }

    let Some(path) = check_path else {
        return ExitCode::SUCCESS;
    };
    let baseline = match read_baseline(&path) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("{e}");
            return ExitCode::from(2);
        }
    };
    let mut failed = false;
    for (name, now) in &results {
        let Some((_, base)) = baseline.iter().find(|(k, _)| k == name) else {
            eprintln!("{name}: no baseline entry, skipped");
            continue;
        };
        let ratio = now / base;
        let status = if ratio > 1.0 + MAX_REGRESSION { "REGRESSION" } else { "ok" };
        println!("{name}: baseline {base:.3} now {now:.3} ratio {ratio:.3} {status}");
        if ratio > 1.0 + MAX_REGRESSION {
            failed = true;
        }
    }
    if failed {
        eprintln!("performance regression above {:.0}%", MAX_REGRESSION * 100.0);
        ExitCode::from(1)
    } else {
        ExitCode::SUCCESS
    }
}
