//! Robustness tests for the network-facing readers added with the distributed
//! executor: the row-store header, the HTTP response parser, and the worker's
//! wire decoder. Garbage must produce errors. A worker must keep serving valid
//! jobs after receiving it.
//!
//! `FUZZ_ITERS` sets the iteration count (default 300; TCP tests are slower).

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::Duration;

use tpt_dynamics_data::http::HttpRange;
use tpt_dynamics_data::lazy::distributed::{collect, serve, DistributedPlan, Resolver};
use tpt_dynamics_data::lazy::{col, AggExpr, Column, DataFrame, LazyError, Source};
use tpt_dynamics_data::rng::SplitMix64;
use tpt_dynamics_data::rowstore::{write_rowstore, RangeRead, RowStore};

fn iterations() -> usize {
    std::env::var("FUZZ_ITERS").ok().and_then(|v| v.parse().ok()).unwrap_or(300)
}

/// A row store held in memory, for header tests without a file.
struct Bytes(Vec<u8>);

impl RangeRead for Bytes {
    fn byte_len(&self) -> Result<u64, LazyError> {
        Ok(self.0.len() as u64)
    }

    fn read_at(&self, offset: u64, buf: &mut [u8]) -> Result<(), LazyError> {
        let start = usize::try_from(offset).map_err(|_| LazyError::Invalid("offset".into()))?;
        let end = start.checked_add(buf.len()).filter(|&e| e <= self.0.len());
        let end = end.ok_or_else(|| LazyError::Invalid("read past end".into()))?;
        buf.copy_from_slice(&self.0[start..end]);
        Ok(())
    }
}

fn sample_frame() -> DataFrame {
    DataFrame::new(vec![
        ("pt", Column::F64(vec![1.0, 2.0, 3.0, 4.0])),
        ("charge", Column::I64(vec![1, -1, 1, -1])),
        ("flag", Column::Bool(vec![true, false, true, true])),
    ])
    .expect("frame")
}

fn sample_store_bytes() -> Vec<u8> {
    // Tests run on parallel threads in one process, so the name must be unique
    // per call; otherwise one test's `remove_file` can delete another's file.
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let n = NEXT.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("tpt-fuzz-net-{}-{n}.rows", std::process::id()));
    write_rowstore(&path, &sample_frame()).expect("write");
    let bytes = std::fs::read(&path).expect("read");
    let _ = std::fs::remove_file(path);
    bytes
}

fn mutate(rng: &mut SplitMix64, base: &[u8]) -> Vec<u8> {
    let mut out = base.to_vec();
    for _ in 0..=(rng.next_u64() % 5) {
        if out.is_empty() {
            break;
        }
        let at = (rng.next_u64() % out.len() as u64) as usize;
        match rng.next_u64() % 3 {
            0 => out[at] = rng.next_u64() as u8,
            1 => out[at] ^= 1 << (rng.next_u64() % 8),
            // Large values over a length field, to hit allocation paths.
            _ => {
                let big = [0xFFu8; 8];
                let end = (at + big.len()).min(out.len());
                out.splice(at..end, big[..end - at].iter().copied());
            }
        }
    }
    out
}

#[test]
fn row_store_headers_never_panic() {
    let valid = sample_store_bytes();
    // The valid file must open and read back.
    let store = RowStore::open_range(Arc::new(Bytes(valid.clone())), "valid").expect("valid store");
    assert_eq!(store.num_rows(), 4);

    let mut rng = SplitMix64::new(0x5EED);
    for _ in 0..iterations() * 20 {
        let garbage: Vec<u8> =
            (0..(rng.next_u64() % 64) as usize).map(|_| rng.next_u64() as u8).collect();
        let _ = RowStore::open_range(Arc::new(Bytes(garbage)), "garbage");
        let mutated = mutate(&mut rng, &valid);
        if let Ok(store) = RowStore::open_range(Arc::new(Bytes(mutated)), "mutated") {
            // An accepted header must also read without panicking.
            let names: Vec<String> = store.schema().into_iter().map(|f| f.name).collect();
            let _ = store.read(0..store.num_rows(), &names);
        }
    }
}

#[test]
fn huge_row_count_in_header_is_rejected_quickly() {
    let mut bytes = sample_store_bytes();
    // Claim u64::MAX rows: the size check must fail before any allocation.
    bytes[8..16].copy_from_slice(&u64::MAX.to_le_bytes());
    assert!(RowStore::open_range(Arc::new(Bytes(bytes)), "huge").is_err());
}

/// Serves one response per connection, taken from `next` at the time of the request.
fn scripted_http_server(next: Arc<Mutex<Vec<u8>>>) -> std::net::SocketAddr {
    let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
    let addr = listener.local_addr().expect("addr");
    std::thread::spawn(move || {
        for stream in listener.incoming() {
            let Ok(mut stream) = stream else { continue };
            let mut req = [0u8; 2048];
            let _ = stream.read(&mut req);
            let response = next.lock().expect("lock").clone();
            let _ = stream.write_all(&response);
        }
    });
    addr
}

#[test]
fn http_response_parser_never_panics() {
    let script = Arc::new(Mutex::new(Vec::new()));
    let addr = scripted_http_server(Arc::clone(&script));
    let url = format!("http://{addr}/x.rows");
    let valid = b"HTTP/1.1 206 Partial Content\r\nContent-Range: bytes 0-0/10\r\nContent-Length: 1\r\n\r\nX".to_vec();
    *script.lock().unwrap() = valid.clone();
    assert!(HttpRange::new(&url).is_ok());

    let mut rng = SplitMix64::new(0xABCD);
    for _ in 0..iterations() {
        let response = if rng.next_u64().is_multiple_of(2) {
            (0..(rng.next_u64() % 200) as usize).map(|_| rng.next_u64() as u8).collect()
        } else {
            mutate(&mut rng, &valid)
        };
        *script.lock().unwrap() = response;
        // Any outcome is acceptable except a panic.
        let _ = HttpRange::new(&url);
    }
}

fn send_frame(stream: &mut TcpStream, payload: &[u8]) -> std::io::Result<()> {
    stream.write_all(&(payload.len() as u32).to_le_bytes())?;
    stream.write_all(payload)
}

#[test]
fn worker_survives_garbage_and_keeps_serving() {
    let df = DataFrame::new(vec![
        ("pt", Column::F64((0..500).map(|i| f64::from(i % 50)).collect())),
        ("charge", Column::I64((0..500).map(|i| if i % 2 == 0 { 1 } else { -1 }).collect())),
    ])
    .expect("frame");
    let df = Arc::new(df);
    let listener = TcpListener::bind("127.0.0.1:0").expect("bind");
    let addr = listener.local_addr().expect("addr");
    let resolve: Resolver = {
        let df = Arc::clone(&df);
        Arc::new(move |name| {
            if name == "events" {
                Ok(Arc::clone(&df) as Arc<dyn Source>)
            } else {
                Err(LazyError::Invalid("unknown source".into()))
            }
        })
    };
    std::thread::spawn(move || serve(&listener, &resolve, None));

    let mut rng = SplitMix64::new(0x7070);
    for _ in 0..iterations() {
        let Ok(mut stream) = TcpStream::connect(addr) else { panic!("worker stopped listening") };
        stream.set_read_timeout(Some(Duration::from_secs(5))).expect("timeout");
        let len = (rng.next_u64() % 120) as usize;
        let payload: Vec<u8> = (0..len).map(|_| rng.next_u64() as u8).collect();
        let _ = send_frame(&mut stream, &payload);
        // Read whatever comes back (an error reply or a closed connection).
        let mut sink = Vec::new();
        let _ = stream.read_to_end(&mut sink);
    }

    // After the garbage, a real plan must still run through this worker.
    let plan = DistributedPlan::new()
        .filter(col("pt").gt(10.0))
        .group_by_agg(["charge"], [AggExpr::count()]);
    let local = plan.to_lazy(Arc::clone(&df) as Arc<dyn Source>).collect().expect("local");
    let remote = collect(&plan, "events", df.as_ref(), &[addr]).expect("distributed after garbage");
    assert_eq!(remote.column("count"), local.column("count"));
}
