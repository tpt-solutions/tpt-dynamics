//! Randomised robustness tests for the data parsers (Phase 8 fuzzing campaign,
//! stable-toolchain half).
//!
//! Each target is fed (a) random bytes, (b) random text from a parser-relevant
//! alphabet, (c) every truncation of a valid input, and (d) random mutations of a
//! valid input (bit flips, byte overwrites, insertions, deletions). The contract
//! is that every input either parses or returns `Err`; none may panic, hang, or
//! allocate in proportion to a length field it has not checked.
//!
//! The iteration count is `FUZZ_ITERS` (default `10_000` per generator). Run with
//! `--release` for a larger campaign:
//! `FUZZ_ITERS=200000 cargo test -p tpt-dynamics-data --release --test fuzz_parsers`.
//! See `docs/fuzzing.md` for the cargo-fuzz setup used for longer runs.

use std::io::BufReader;

use tpt_dynamics_data::arrow_ipc::{self, DataType, Values};
use tpt_dynamics_data::columnar::{self, ColumnData};
use tpt_dynamics_data::csv::{CsvOptions, CsvReader};
use tpt_dynamics_data::rng::SplitMix64;
use tpt_dynamics_data::{Histogram1D, Histogram2D, Profile1D, Profile2D};

fn iters() -> usize {
    std::env::var("FUZZ_ITERS").ok().and_then(|v| v.parse().ok()).unwrap_or(10_000)
}

/// Random bytes, or random bytes drawn from a small alphabet of interesting
/// characters, with lengths up to `max_len`.
fn random_input(rng: &mut SplitMix64, alphabet: &[u8], max_len: u64) -> Vec<u8> {
    let len = rng.below(max_len + 1) as usize;
    let biased = rng.below(2) == 0;
    (0..len)
        .map(|_| {
            if biased {
                alphabet[rng.below(alphabet.len() as u64) as usize]
            } else {
                rng.next_u64() as u8
            }
        })
        .collect()
}

/// Apply 1..=4 random structural mutations to a valid input.
fn mutate(rng: &mut SplitMix64, valid: &[u8]) -> Vec<u8> {
    let mut out = valid.to_vec();
    let ops = 1 + rng.below(4);
    for _ in 0..ops {
        if out.is_empty() {
            out.push(rng.next_u64() as u8);
            continue;
        }
        let i = rng.below(out.len() as u64) as usize;
        match rng.below(5) {
            0 => out[i] ^= 1 << rng.below(8),
            1 => out[i] = rng.next_u64() as u8,
            2 => out.insert(i, rng.next_u64() as u8),
            3 => {
                out.remove(i);
            }
            _ => {
                // Overwrite a 4-byte window with a large little-endian value, which is
                // how length fields are most likely to be corrupted.
                let big = [0xFF, 0xFF, 0xFF, 0x7F];
                for (k, b) in big.iter().enumerate() {
                    if let Some(slot) = out.get_mut(i + k) {
                        *slot = *b;
                    }
                }
            }
        }
    }
    out
}

fn drain_csv(data: &[u8], opts: CsvOptions) {
    let Ok(reader) = CsvReader::new(BufReader::new(data), opts) else {
        return;
    };
    let mut records = 0usize;
    for rec in reader {
        records += 1;
        assert!(records <= data.len() + 1, "CSV reader did not terminate");
        if let Ok(rec) = rec {
            for i in 0..rec.len() {
                let _ = rec.f64(i);
                let _ = rec.i64(i);
                let _ = rec.u64(i);
            }
        }
    }
}

#[test]
fn csv_random_inputs_never_panic() {
    let alphabet = b"0123456789.,-+eE\"\n\r\t ab\xff\x00;#";
    let mut rng = SplitMix64::new(0xC5F0_0001);
    for i in 0..iters() {
        let data = random_input(&mut rng, alphabet, 400);
        let opts = CsvOptions {
            delimiter: if i % 3 == 0 { b';' } else { b',' },
            has_header: i % 4 != 0,
            max_line_bytes: 1 + (i % 96),
        };
        drain_csv(&data, opts);
    }
}

#[test]
fn csv_truncations_and_mutations_never_panic() {
    let valid = b"run,pt,eta\n1,45.5,0.25\n1,12.0,-1.5\n2,7.25,2.0\n".to_vec();
    for n in 0..=valid.len() {
        drain_csv(&valid[..n], CsvOptions::default());
    }
    let mut rng = SplitMix64::new(0xC5F0_0002);
    for _ in 0..iters() {
        let data = mutate(&mut rng, &valid);
        drain_csv(&data, CsvOptions::default());
    }
}

fn valid_columnar() -> Vec<u8> {
    let cols = vec![
        ("pt".to_owned(), ColumnData::F64(vec![1.0, 2.5, -3.0])),
        ("n".to_owned(), ColumnData::I64(vec![1, -2, 3])),
        ("id".to_owned(), ColumnData::U32(vec![7, 8, 9])),
    ];
    columnar::write_table(Vec::new(), &cols).unwrap()
}

#[test]
fn columnar_random_inputs_never_panic() {
    let alphabet = b"TPTC\x01\x00\x02\x00pt\xff\x00\x80";
    let mut rng = SplitMix64::new(0xC01A_0001);
    for _ in 0..iters() {
        let data = random_input(&mut rng, alphabet, 300);
        let _ = columnar::read_table(&data[..]);
    }
}

#[test]
fn columnar_truncations_and_mutations_never_panic() {
    let valid = valid_columnar();
    for n in 0..valid.len() {
        assert!(columnar::read_table(&valid[..n]).is_err(), "prefix {n} accepted");
    }
    let mut rng = SplitMix64::new(0xC01A_0002);
    for _ in 0..iters() {
        let data = mutate(&mut rng, &valid);
        let _ = columnar::read_table(&data[..]);
    }
}

fn valid_histogram_bytes() -> Vec<Vec<u8>> {
    let mut h1 = Histogram1D::uniform(8, -2.0, 2.0).unwrap();
    let mut h2 = Histogram2D::uniform(4, 0.0, 1.0, 3, 0.0, 1.0).unwrap();
    let mut p1 = Profile1D::uniform(5, 0.0, 1.0).unwrap();
    let mut p2 = Profile2D::uniform(2, 0.0, 1.0, 2, 0.0, 1.0).unwrap();
    for k in 0..40 {
        let x = f64::from(k) * 0.1 - 2.0;
        h1.fill_weighted(x, 0.5);
        h2.fill(x.abs() % 1.0, (x * 1.7).abs() % 1.0);
        p1.fill(x.abs() % 1.0, x);
        p2.fill(x.abs() % 1.0, (x * 2.0).abs() % 1.0, x);
    }
    vec![h1.to_bytes(), h2.to_bytes(), p1.to_bytes(), p2.to_bytes()]
}

/// Decode with every serialise reader; the result is ignored, only panics matter.
fn decode_all_binary(data: &[u8]) {
    let _ = Histogram1D::from_bytes(data);
    let _ = Histogram2D::from_bytes(data);
    let _ = Profile1D::from_bytes(data);
    let _ = Profile2D::from_bytes(data);
}

#[test]
fn serialise_random_binary_never_panics() {
    let alphabet = b"TPTH\x01\x00\x00\x01\xff\x7f";
    let mut rng = SplitMix64::new(0x5E1A_0001);
    for _ in 0..iters() {
        let data = random_input(&mut rng, alphabet, 256);
        decode_all_binary(&data);
    }
}

#[test]
fn serialise_truncations_and_mutations_never_panic() {
    let mut rng = SplitMix64::new(0x5E1A_0002);
    for valid in valid_histogram_bytes() {
        for n in 0..valid.len() {
            decode_all_binary(&valid[..n]);
        }
        for _ in 0..iters() / 4 {
            let data = mutate(&mut rng, &valid);
            decode_all_binary(&data);
        }
    }
}

#[test]
fn serialise_random_json_never_panics() {
    let alphabet = br#"{}[]:,"0123456789.-+eE truefalsnul\\u"#;
    let mut rng = SplitMix64::new(0x5E1A_0003);
    for _ in 0..iters() {
        let len = rng.below(200) as usize;
        let text: String = (0..len)
            .map(|_| char::from(alphabet[rng.below(alphabet.len() as u64) as usize]))
            .collect();
        let _ = Histogram1D::from_json(&text);
        let _ = Histogram2D::from_json(&text);
        let _ = Profile1D::from_json(&text);
        let _ = Profile2D::from_json(&text);
    }
}

#[test]
fn arrow_ipc_random_and_mutated_inputs_never_panic() {
    let valid = arrow_ipc::test_support::encode(
        &[("a", DataType::Int64, true), ("b", DataType::Bool, false)],
        &[(
            5,
            vec![
                (Values::Int64(vec![1, 2, 3, 4, 5]), Some(vec![true, false, true, true, true])),
                (Values::Bool(vec![true, false, true, false, true]), None),
            ],
        )],
    );
    let alphabet = b"TPTA\x01\x00\x00\x00\x01\x03\xff";
    let mut rng = SplitMix64::new(0xA770_0001);
    for _ in 0..iters() {
        let data = random_input(&mut rng, alphabet, 300);
        let _ = arrow_ipc::read_all(&data);
    }
    for n in 0..valid.len() {
        assert!(arrow_ipc::read_all(&valid[..n]).is_err());
    }
    for _ in 0..iters() {
        let data = mutate(&mut rng, &valid);
        let _ = arrow_ipc::read_all(&data);
    }
}

/// Parquet import (feature `parquet`): random, truncated and mutated inputs must
/// produce `Ok` or `Err`, never a panic. Seeds are the committed fixtures.
#[cfg(feature = "parquet")]
mod parquet_targets {
    use super::*;
    use tpt_dynamics_data::parquet::ParquetFile;

    fn seeds() -> Vec<Vec<u8>> {
        ["flat_snappy_plain.parquet", "flat_gzip_dict.parquet", "flat_none_dict.parquet"]
            .iter()
            .map(|name| {
                let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                    .join("tests/fixtures/parquet")
                    .join(name);
                std::fs::read(path).unwrap()
            })
            .collect()
    }

    fn decode(data: &[u8]) {
        if let Ok(mut f) = ParquetFile::open(std::io::Cursor::new(data.to_vec())) {
            let _ = f.read_frame(None);
        }
    }

    #[test]
    fn parquet_random_inputs_never_panic() {
        let mut rng = SplitMix64::new(0x9A4B_0001);
        let alphabet = b"PAR1\x00\x01\x02\xff";
        for _ in 0..iters() {
            decode(&random_input(&mut rng, alphabet, 512));
        }
    }

    #[test]
    fn parquet_truncations_and_mutations_never_panic() {
        let mut rng = SplitMix64::new(0x9A4B_0002);
        for valid in seeds() {
            for n in 0..valid.len() {
                decode(&valid[..n]);
            }
            for _ in 0..iters() / 4 {
                decode(&mutate(&mut rng, &valid));
            }
        }
    }
}

/// Chunked trajectory format: random, truncated and mutated inputs must be
/// rejected or decoded without a panic, including random frame access.
mod trajectory_targets {
    use super::*;
    use tpt_dynamics_data::trajectory::{Frame, FrameSink, Precision, TrajectoryReader, TrajectoryWriter};

    fn valid_trajectory() -> Vec<u8> {
        let mut w = TrajectoryWriter::new(Vec::new(), 3, Precision::F32, 4).unwrap();
        for step in 0..10u64 {
            let f = Frame {
                step,
                time: 0.002 * step as f64,
                box_lengths: [3.0, 3.0, 3.0],
                positions: vec![[0.1 * step as f64, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 0.5, 0.25]],
            };
            w.push_frame(&f).unwrap();
        }
        w.finish().unwrap()
    }

    fn decode(data: &[u8], rng: &mut SplitMix64) {
        let Ok(mut r) = TrajectoryReader::open(std::io::Cursor::new(data.to_vec())) else {
            return;
        };
        // Bounded: a corrupt header must not be able to stall the test.
        for item in r.frames().take(64) {
            let _ = item;
        }
        let n = r.n_frames();
        for _ in 0..4 {
            let i = rng.below(n + 2);
            let _ = r.frame(i);
        }
    }

    #[test]
    fn trajectory_random_inputs_never_panic() {
        let mut rng = SplitMix64::new(0x7A7_0001);
        let alphabet = b"TRJ\x01\x00\x00\x00\xff";
        for _ in 0..iters() {
            let data = random_input(&mut rng, alphabet, 512);
            decode(&data, &mut rng);
        }
    }

    #[test]
    fn trajectory_truncations_and_mutations_never_panic() {
        let mut rng = SplitMix64::new(0x7A7_0002);
        let valid = valid_trajectory();
        for n in 0..valid.len() {
            decode(&valid[..n], &mut rng);
        }
        for _ in 0..iters() {
            let data = mutate(&mut rng, &valid);
            decode(&data, &mut rng);
        }
    }
}
