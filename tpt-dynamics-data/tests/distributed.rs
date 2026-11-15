//! Distributed execution over loopback TCP: workers run in threads, the
//! coordinator talks to them over real sockets, and every result is checked
//! against a local run of the same plan.

use std::net::{SocketAddr, TcpListener};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use tpt_dynamics_data::lazy::distributed::{collect, serve, DistributedPlan, Resolver};
use tpt_dynamics_data::lazy::{col, lit, AggExpr, Column, DataFrame, LazyError, Source};
use tpt_dynamics_data::rowstore::{write_rowstore, RowStore};

/// A deterministic event table: integer-valued `pt` (so float sums are exact
/// whatever the summation order), a fractional `eta`, an integer `charge` and a
/// boolean `isolated`.
fn events(n: usize) -> DataFrame {
    let mut state = 12345u64;
    let mut next = move || {
        state =
            state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        state >> 11
    };
    let mut pt = Vec::with_capacity(n);
    let mut eta = Vec::with_capacity(n);
    let mut charge = Vec::with_capacity(n);
    let mut isolated = Vec::with_capacity(n);
    for _ in 0..n {
        pt.push((next() % 100) as f64);
        eta.push((next() % 1000) as f64 / 500.0 - 1.0);
        charge.push(if next() % 2 == 0 { 1 } else { -1 });
        isolated.push(next() % 3 == 0);
    }
    DataFrame::new(vec![
        ("pt", Column::F64(pt)),
        ("eta", Column::F64(eta)),
        ("charge", Column::I64(charge)),
        ("isolated", Column::Bool(isolated)),
    ])
    .expect("valid frame")
}

fn temp_store(name: &str, df: &DataFrame) -> PathBuf {
    let path = std::env::temp_dir().join(format!("tpt-dist-{}-{name}.rows", std::process::id()));
    write_rowstore(&path, df).expect("write row store");
    path
}

/// Starts `n` worker threads, each serving `jobs` connections, and returns
/// their addresses. Each worker resolves `"events"` to the row store at `path`.
fn start_workers(n: usize, jobs: usize, path: &Path) -> Vec<SocketAddr> {
    (0..n)
        .map(|_| {
            let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
            let addr = listener.local_addr().expect("addr");
            let path = path.to_path_buf();
            let resolve: Resolver = Arc::new(move |name| {
                if name != "events" {
                    return Err(LazyError::Invalid(format!("unknown source '{name}'")));
                }
                Ok(Arc::new(RowStore::open(&path)?))
            });
            std::thread::spawn(move || serve(&listener, &resolve, Some(jobs)));
            addr
        })
        .collect()
}

fn assert_frames_equal(a: &DataFrame, b: &DataFrame) {
    assert_eq!(a.names(), b.names(), "column names differ");
    assert_eq!(a.len(), b.len(), "row counts differ");
    for name in a.names() {
        assert_eq!(a.column(name), b.column(name), "column '{name}' differs");
    }
}

#[test]
fn row_plan_matches_local_execution() {
    let df = events(10_000);
    let path = temp_store("rows", &df);
    let store = RowStore::open(&path).unwrap();
    let plan = DistributedPlan::new()
        .filter(col("pt").gt(30.0))
        .with_column((col("pt") * lit(2.0)).alias("pt2"))
        .select([col("pt2"), col("charge"), col("isolated")]);
    let local = plan.to_lazy(Arc::new(store.clone())).collect().unwrap();
    let workers = start_workers(4, 1, &path);
    let remote = collect(&plan, "events", &store, &workers).unwrap();
    assert_frames_equal(&remote, &local);
    assert!(!remote.is_empty() && remote.len() < 10_000);
    let _ = std::fs::remove_file(path);
}

#[test]
fn group_by_matches_local_execution() {
    let df = events(20_000);
    let path = temp_store("agg", &df);
    let store = RowStore::open(&path).unwrap();
    let plan = DistributedPlan::new().filter(col("isolated").eq_to(true)).group_by_agg(
        ["charge"],
        [
            AggExpr::count(),
            AggExpr::sum(col("pt")).alias("pt_sum"),
            AggExpr::min(col("pt")).alias("pt_min"),
            AggExpr::max(col("pt")).alias("pt_max"),
            AggExpr::mean(col("eta")).alias("eta_mean"),
        ],
    );
    let local = plan.to_lazy(Arc::new(store.clone())).collect().unwrap();
    let workers = start_workers(3, 1, &path);
    let remote = collect(&plan, "events", &store, &workers).unwrap();

    assert_eq!(remote.names(), local.names());
    assert_eq!(remote.len(), local.len());
    // Counts, sums of integers, min and max are exact.
    for name in ["count", "pt_sum", "pt_min", "pt_max"] {
        assert_eq!(remote.column(name), local.column(name), "column '{name}' differs");
    }
    // Mean of a fractional column depends on summation order only in the last bits.
    let (Some(Column::F64(remote_eta)), Some(Column::F64(local_eta))) =
        (remote.column("eta_mean"), local.column("eta_mean"))
    else {
        panic!("eta_mean must be f64");
    };
    for (x, y) in remote_eta.iter().zip(local_eta) {
        assert!((x - y).abs() <= 1e-12 * y.abs().max(1.0), "eta_mean {x} vs {y}");
    }
    let _ = std::fs::remove_file(path);
}

#[test]
fn more_workers_than_rows_still_matches() {
    let df = events(3);
    let path = temp_store("tiny", &df);
    let store = RowStore::open(&path).unwrap();
    let plan = DistributedPlan::new().group_by_agg(["charge"], [AggExpr::count()]);
    let local = plan.to_lazy(Arc::new(store.clone())).collect().unwrap();
    let workers = start_workers(5, 1, &path);
    let remote = collect(&plan, "events", &store, &workers).unwrap();
    assert_frames_equal(&remote, &local);
    let _ = std::fs::remove_file(path);
}

#[test]
fn empty_source_gives_empty_result() {
    let df = events(0);
    let path = temp_store("empty", &df);
    let store = RowStore::open(&path).unwrap();
    let plan = DistributedPlan::new().group_by_agg(["charge"], [AggExpr::count()]);
    let workers = start_workers(2, 1, &path);
    let remote = collect(&plan, "events", &store, &workers).unwrap();
    assert_eq!(remote.len(), 0);
    let _ = std::fs::remove_file(path);
}

#[test]
fn worker_error_is_reported_to_the_coordinator() {
    let df = events(100);
    let path = temp_store("err", &df);
    let store = RowStore::open(&path).unwrap();
    let plan = DistributedPlan::new().filter(col("pt").gt(1.0));
    let workers = start_workers(2, 1, &path);
    let err = collect(&plan, "not-events", &store, &workers).unwrap_err();
    assert!(err.to_string().contains("unknown source"), "got: {err}");
    let _ = std::fs::remove_file(path);
}

#[test]
fn invalid_plan_is_rejected_before_any_job_is_sent() {
    let df = events(10);
    let path = temp_store("badplan", &df);
    let store = RowStore::open(&path).unwrap();
    // No workers are listening here, so a connection attempt would fail with a
    // different error: validation must come first.
    let plan = DistributedPlan::new().filter(col("no_such_column").gt(1.0));
    let nowhere: SocketAddr = "127.0.0.1:9".parse().unwrap();
    let err = collect(&plan, "events", &store, &[nowhere]).unwrap_err();
    assert!(!err.to_string().contains("worker"), "validation should run first: {err}");
    let _ = std::fs::remove_file(path);
}

#[test]
fn no_workers_is_an_error() {
    let df = events(10);
    let path = temp_store("noworkers", &df);
    let store = RowStore::open(&path).unwrap();
    assert!(collect(&DistributedPlan::new(), "events", &store, &[]).is_err());
    let _ = std::fs::remove_file(path);
}

#[test]
fn row_store_serves_arbitrary_ranges() {
    let df = events(1000);
    let path = temp_store("ranges", &df);
    let store = RowStore::open(&path).unwrap();
    assert_eq!(store.num_rows(), 1000);
    let names: Vec<String> = df.names().to_vec();
    let slice = store.read(250..400, &names).unwrap();
    assert_eq!(slice.len(), 150);
    assert_eq!(slice.column("pt"), df.column("pt").map(|c| slice_of(c, 250, 400)).as_ref());
    // Ranges past the end are clamped.
    assert_eq!(store.read(990..2000, &names).unwrap().len(), 10);
    let _ = std::fs::remove_file(path);
}

fn slice_of(c: &Column, a: usize, b: usize) -> Column {
    match c {
        Column::F64(v) => Column::F64(v[a..b].to_vec()),
        Column::I64(v) => Column::I64(v[a..b].to_vec()),
        Column::Bool(v) => Column::Bool(v[a..b].to_vec()),
    }
}

#[test]
fn truncated_row_store_is_rejected() {
    let df = events(100);
    let path = temp_store("trunc", &df);
    let bytes = std::fs::read(&path).unwrap();
    std::fs::write(&path, &bytes[..bytes.len() - 1]).unwrap();
    assert!(RowStore::open(&path).is_err());
    std::fs::write(&path, b"not a row store at all, really").unwrap();
    assert!(RowStore::open(&path).is_err());
    let _ = std::fs::remove_file(path);
}
