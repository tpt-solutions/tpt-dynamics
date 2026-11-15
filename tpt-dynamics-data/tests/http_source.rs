//! Row stores served over plain HTTP `Range` requests, checked against a tiny
//! loopback server. Also runs a distributed plan whose workers read their
//! ranges over HTTP.

use std::io::{Read, Write};
use std::net::{SocketAddr, TcpListener, TcpStream};
use std::path::PathBuf;
use std::sync::Arc;

use tpt_dynamics_data::lazy::distributed::{collect, serve, DistributedPlan, Resolver};
use tpt_dynamics_data::lazy::{col, AggExpr, Column, DataFrame, LazyError, Source};
use tpt_dynamics_data::rowstore::{write_rowstore, RowStore};

fn events(n: usize) -> DataFrame {
    let pt: Vec<f64> = (0..n).map(|i| (i * 37 % 101) as f64).collect();
    let charge: Vec<i64> = (0..n).map(|i| if i % 3 == 0 { -1 } else { 1 }).collect();
    let flag: Vec<bool> = (0..n).map(|i| i % 5 == 0).collect();
    DataFrame::new(vec![
        ("pt", Column::F64(pt)),
        ("charge", Column::I64(charge)),
        ("flag", Column::Bool(flag)),
    ])
    .expect("valid frame")
}

/// Answers HTTP/1.1 `GET` requests with optional `Range` headers from `data`.
/// `ignore_ranges` makes it answer every request with `200` and the full body.
fn serve_bytes(data: Arc<Vec<u8>>, ignore_ranges: bool) -> SocketAddr {
    let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
    let addr = listener.local_addr().expect("addr");
    std::thread::spawn(move || {
        for stream in listener.incoming() {
            let Ok(mut stream) = stream else { continue };
            let data = Arc::clone(&data);
            std::thread::spawn(move || handle(&mut stream, &data, ignore_ranges));
        }
    });
    addr
}

fn handle(stream: &mut TcpStream, data: &[u8], ignore_ranges: bool) {
    let mut req = Vec::new();
    let mut byte = [0u8; 1];
    while !req.ends_with(b"\r\n\r\n") {
        match stream.read(&mut byte) {
            Ok(1) => req.push(byte[0]),
            _ => return,
        }
    }
    let text = String::from_utf8_lossy(&req).to_string();
    let range = text
        .lines()
        .find_map(|l| l.strip_prefix("Range: bytes="))
        .and_then(|r| r.trim().split_once('-'))
        .and_then(|(a, b)| Some((a.parse::<usize>().ok()?, b.parse::<usize>().ok()?)));
    let response = match range {
        Some((a, b)) if !ignore_ranges => {
            if a >= data.len() {
                b"HTTP/1.1 416 Range Not Satisfiable\r\nContent-Length: 0\r\nConnection: close\r\n\r\n".to_vec()
            } else {
                let b = b.min(data.len() - 1);
                let body = &data[a..=b];
                let mut r = format!(
                    "HTTP/1.1 206 Partial Content\r\nContent-Range: bytes {a}-{b}/{}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    data.len(),
                    body.len()
                )
                .into_bytes();
                r.extend_from_slice(body);
                r
            }
        }
        _ => {
            let mut r = format!(
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                data.len()
            )
            .into_bytes();
            r.extend_from_slice(data);
            r
        }
    };
    let _ = stream.write_all(&response);
}

fn store_bytes(df: &DataFrame, name: &str) -> (PathBuf, Arc<Vec<u8>>) {
    let path = std::env::temp_dir().join(format!("tpt-http-{}-{name}.rows", std::process::id()));
    write_rowstore(&path, df).expect("write");
    let bytes = std::fs::read(&path).expect("read back");
    (path, Arc::new(bytes))
}

#[test]
fn http_store_matches_file_store() {
    let df = events(5_000);
    let (path, bytes) = store_bytes(&df, "match");
    let addr = serve_bytes(bytes, false);
    let url = format!("http://{addr}/events.rows");

    let local = RowStore::open(&path).unwrap();
    let remote = RowStore::open_http(&url).unwrap();
    assert_eq!(remote.num_rows(), local.num_rows());
    assert_eq!(remote.schema(), local.schema());
    let names: Vec<String> = df.names().to_vec();
    for range in [0..5_000, 0..1, 4_999..5_000, 1_234..2_345] {
        let a = local.read(range.clone(), &names).unwrap();
        let b = remote.read(range, &names).unwrap();
        for n in &names {
            assert_eq!(a.column(n), b.column(n), "column {n}");
        }
    }
    let _ = std::fs::remove_file(path);
}

#[test]
fn distributed_plan_with_http_backed_workers() {
    let df = events(12_000);
    let (path, bytes) = store_bytes(&df, "dist");
    let addr = serve_bytes(bytes, false);
    let url = format!("http://{addr}/events.rows");

    let local_store = RowStore::open(&path).unwrap();
    let plan = DistributedPlan::new()
        .filter(col("pt").gt(10.0))
        .group_by_agg(["charge"], [AggExpr::count(), AggExpr::sum(col("pt")).alias("pt_sum")]);
    let expected = plan.to_lazy(Arc::new(local_store.clone())).collect().unwrap();

    // Each worker opens the store over HTTP for every job it serves.
    let mut workers = Vec::new();
    for _ in 0..3 {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        workers.push(listener.local_addr().unwrap());
        let url = url.clone();
        let resolve: Resolver = Arc::new(move |name| {
            if name != "events" {
                return Err(LazyError::Invalid(format!("unknown source '{name}'")));
            }
            Ok(Arc::new(RowStore::open_http(&url)?))
        });
        std::thread::spawn(move || serve(&listener, &resolve, Some(1)));
    }
    let got = collect(&plan, "events", &local_store, &workers).unwrap();
    assert_eq!(got.column("count"), expected.column("count"));
    assert_eq!(got.column("pt_sum"), expected.column("pt_sum"));
    let _ = std::fs::remove_file(path);
}

#[test]
fn server_without_range_support_is_rejected() {
    let (path, bytes) = store_bytes(&events(100), "norange");
    let addr = serve_bytes(bytes, true);
    let err = RowStore::open_http(&format!("http://{addr}/x.rows")).unwrap_err();
    assert!(err.to_string().contains("status 200"), "got: {err}");
    let _ = std::fs::remove_file(path);
}

#[test]
fn clamped_reads_over_http() {
    let df = events(50);
    let (path, bytes) = store_bytes(&df, "past");
    let addr = serve_bytes(bytes, false);
    let remote = RowStore::open_http(&format!("http://{addr}/x.rows")).unwrap();
    // Asking for rows beyond the end is clamped by Source::read, so the store itself
    // stays valid; a raw read past the end of the resource must fail.
    assert_eq!(remote.read(40..100, &["pt".to_string()]).unwrap().len(), 10);
    let _ = std::fs::remove_file(path);
}

#[test]
fn https_and_unreachable_hosts_are_errors() {
    assert!(RowStore::open_http("https://example.org/x.rows").is_err());
    // Port 9 (discard) on loopback is almost never open; the connect must fail cleanly.
    assert!(RowStore::open_http("http://127.0.0.1:9/x.rows").is_err());
}
