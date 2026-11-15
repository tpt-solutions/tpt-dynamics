//! Parquet import against files written by pyarrow (see
//! `tests/fixtures/parquet/`; regenerate with pyarrow, `data_page_version="1.0"`).
#![cfg(feature = "parquet")]

use std::fs::File;
use std::path::PathBuf;

use tpt_dynamics_data::lazy::{col, Column, LazyFrame};
use tpt_dynamics_data::parquet::{read_frame, ParquetError, ParquetFile};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/parquet").join(name)
}

fn assert_flat(df: &tpt_dynamics_data::lazy::DataFrame) {
    assert_eq!(df.len(), 1000);
    let Some(Column::I64(i64s)) = df.column("i64") else {
        panic!("i64 column");
    };
    let Some(Column::I64(i32s)) = df.column("i32") else {
        panic!("i32 column");
    };
    let Some(Column::F64(f64s)) = df.column("f64") else {
        panic!("f64 column");
    };
    let Some(Column::F64(f32s)) = df.column("f32") else {
        panic!("f32 column");
    };
    let Some(Column::Bool(flags)) = df.column("flag") else {
        panic!("flag column");
    };
    let Some(Column::I64(cats)) = df.column("cat") else {
        panic!("cat column");
    };
    for i in 0..1000usize {
        let k = i64::try_from(i).unwrap();
        assert_eq!(i64s[i], k * 7 - 300, "i64 row {i}");
        assert_eq!(i32s[i], k - 500, "i32 row {i}");
        assert_eq!(f64s[i], i as f64 * 0.25 - 10.0, "f64 row {i}");
        assert_eq!(f32s[i], i as f64 * 0.5, "f32 row {i}");
        assert_eq!(flags[i], i % 3 == 0, "flag row {i}");
        assert_eq!(cats[i], k % 5, "cat row {i}");
    }
}

#[test]
fn reads_every_codec_and_encoding_combination() {
    for name in [
        "flat_none_plain.parquet",
        "flat_none_dict.parquet",
        "flat_snappy_plain.parquet",
        "flat_snappy_dict.parquet",
        "flat_gzip_plain.parquet",
        "flat_gzip_dict.parquet",
        "default_writer.parquet",
    ] {
        let df = read_frame(fixture(name), None).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_flat(&df);
    }
}

#[test]
fn metadata_reports_schema_and_row_groups() {
    let file = ParquetFile::open(File::open(fixture("flat_none_plain.parquet")).unwrap()).unwrap();
    assert_eq!(file.num_rows(), 1000);
    assert_eq!(file.num_row_groups(), 4); // row_group_size=300 -> 4 groups
    let names: Vec<&str> = file.columns().iter().map(|c| c.name.as_str()).collect();
    assert_eq!(names, ["i64", "i32", "f64", "f32", "flag", "cat"]);
}

#[test]
fn column_projection_reads_only_requested_columns() {
    let df = read_frame(fixture("flat_snappy_dict.parquet"), Some(&["f64", "cat"])).unwrap();
    assert_eq!(df.names(), ["f64", "cat"]);
    assert_eq!(df.len(), 1000);
}

#[test]
fn imported_frame_feeds_the_lazy_engine() {
    let df = read_frame(fixture("flat_gzip_dict.parquet"), Some(&["cat", "f64"])).unwrap();
    let out = LazyFrame::scan(df).filter(col("cat").eq_to(0.0)).collect().unwrap();
    assert_eq!(out.len(), 200);
}

#[test]
fn unknown_column_is_reported() {
    let mut f = ParquetFile::open(File::open(fixture("flat_none_plain.parquet")).unwrap()).unwrap();
    assert!(matches!(f.read_column("nope"), Err(ParquetError::UnknownColumn(_))));
}

#[test]
fn optional_column_without_nulls_is_read() {
    let df = read_frame(fixture("optional_no_nulls.parquet"), None).unwrap();
    let Some(Column::F64(x)) = df.column("x") else {
        panic!("x column");
    };
    assert_eq!(x.len(), 20);
    assert_eq!(x[19], 19.0);
}

#[test]
fn null_values_are_rejected_not_silently_dropped() {
    let err = read_frame(fixture("nulls_plain.parquet"), None).unwrap_err();
    assert!(matches!(err, ParquetError::Unsupported(_)), "{err}");
    assert!(err.to_string().contains("null"), "{err}");
}

#[test]
fn string_columns_are_rejected_with_a_clear_error() {
    let err = read_frame(fixture("unsupported_string.parquet"), None).unwrap_err();
    assert!(matches!(err, ParquetError::Unsupported(_)), "{err}");
}

#[test]
fn truncated_file_is_a_format_error() {
    let bytes = std::fs::read(fixture("flat_none_plain.parquet")).unwrap();
    let cut = &bytes[..bytes.len() - 100];
    let err = ParquetFile::open(std::io::Cursor::new(cut.to_vec())).unwrap_err();
    assert!(matches!(err, ParquetError::Format(_)), "{err}");
}

#[test]
fn corrupted_magic_is_a_format_error() {
    let mut bytes = std::fs::read(fixture("flat_none_plain.parquet")).unwrap();
    bytes[0] = b'X';
    let err = ParquetFile::open(std::io::Cursor::new(bytes)).unwrap_err();
    assert!(matches!(err, ParquetError::Format(_)), "{err}");
}

#[test]
fn corrupted_page_data_never_panics() {
    // Flip every byte of the data region in turn; each result must be Ok or Err.
    let original = std::fs::read(fixture("flat_snappy_plain.parquet")).unwrap();
    for i in 4..original.len() - 8 {
        let mut bytes = original.clone();
        bytes[i] ^= 0xA5;
        let _ = ParquetFile::open(std::io::Cursor::new(bytes)).and_then(|mut f| f.read_frame(None));
    }
}
