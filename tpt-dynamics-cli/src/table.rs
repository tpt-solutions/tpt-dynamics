//! Input formats for `tptdyn analyze`, chosen by the file's magic bytes.
//!
//! * `PAR1`: Parquet, read with `tpt_dynamics_data::parquet` (feature `parquet`).
//! * `TPTC`: the columnar table format, read whole with `tpt_dynamics_data::columnar`.
//! * `TPTROWS1`: the row-store format, read whole through `tpt_dynamics_data::rowstore`.
//! * anything else: CSV, streamed by the caller.
//!
//! The binary formats are loaded into memory as one frame. CSV is streamed in chunks.

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use tpt_dynamics_data::columnar::{self, ColumnData};
use tpt_dynamics_data::lazy::{Column, DataFrame, LazyFrame};
use tpt_dynamics_data::rowstore::RowStore;

use crate::error::{input, CliError};

/// Magic bytes of the Parquet format.
const PARQUET_MAGIC: &[u8] = b"PAR1";
/// Magic bytes of the columnar table format (`tpt_dynamics_data::columnar::MAGIC`).
const COLUMNAR_MAGIC: &[u8] = b"TPTC";
/// Magic bytes of the row-store format.
const ROWSTORE_MAGIC: &[u8] = b"TPTROWS1";
/// Magic bytes of the chunked trajectory format, which starts with `TPTC` too.
const CHUNKED_TRAJECTORY_MAGIC: &[u8] = b"TPTCTRJ1";
/// Magic bytes of the plain trajectory format.
const PLAIN_TRAJECTORY_MAGIC: &[u8] = b"TPTTRJ01";

/// The input format of an event table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Kind {
    /// Comma (or other delimiter) separated text, streamed.
    Csv,
    /// Apache Parquet.
    Parquet,
    /// The `tpt-dynamics-data` columnar table format.
    Columnar,
    /// The `tpt-dynamics-data` row-store format.
    RowStore,
}

/// Decides the format of `path` from its first bytes.
///
/// # Errors
/// Input error if the file cannot be read, or if it is a trajectory (which `analyze`
/// reads only with `--trajectory`).
pub(crate) fn detect(path: &Path) -> Result<Kind, CliError> {
    let mut file =
        File::open(path).map_err(|e| input(format!("cannot open '{}': {e}", path.display())))?;
    let mut head = [0u8; 8];
    let mut len = 0;
    while len < head.len() {
        match file.read(&mut head[len..]) {
            Ok(0) => break,
            Ok(n) => len += n,
            Err(e) => return Err(input(format!("cannot read '{}': {e}", path.display()))),
        }
    }
    let head = &head[..len];
    if head.starts_with(CHUNKED_TRAJECTORY_MAGIC) || head.starts_with(PLAIN_TRAJECTORY_MAGIC) {
        return Err(input(format!(
            "'{}' is a trajectory; use `analyze --trajectory` for it",
            path.display()
        )));
    }
    Ok(if head.starts_with(PARQUET_MAGIC) {
        Kind::Parquet
    } else if head.starts_with(ROWSTORE_MAGIC) {
        Kind::RowStore
    } else if head.starts_with(COLUMNAR_MAGIC) {
        Kind::Columnar
    } else {
        Kind::Csv
    })
}

/// Loads a whole binary table into a frame. `Kind::Csv` is not handled here.
///
/// # Errors
/// Input error if the file is missing, truncated, corrupt or uses an unsupported feature.
pub(crate) fn load(path: &Path, kind: Kind) -> Result<DataFrame, CliError> {
    let shown = path.display();
    match kind {
        Kind::Csv => Err(input(format!("'{shown}' is CSV; it is streamed, not loaded"))),
        Kind::Parquet => tpt_dynamics_data::parquet::read_frame(path, None)
            .map_err(|e| input(format!("'{shown}': {e}"))),
        Kind::Columnar => {
            let file =
                File::open(path).map_err(|e| input(format!("cannot open '{shown}': {e}")))?;
            let table = columnar::read_table(BufReader::new(file))
                .map_err(|e| input(format!("'{shown}': {e}")))?;
            let rows = usize::try_from(table.nrows())
                .map_err(|_| input(format!("'{shown}': too many rows for this platform")))?;
            let columns = table.columns().iter().map(|c| {
                let data = match &c.data {
                    ColumnData::F64(v) => Column::F64(v.clone()),
                    ColumnData::I64(v) => Column::I64(v.clone()),
                    ColumnData::U32(v) => Column::I64(v.iter().map(|&x| i64::from(x)).collect()),
                };
                (c.name.clone(), data)
            });
            DataFrame::with_rows(rows, columns).map_err(|e| input(format!("'{shown}': {e}")))
        }
        Kind::RowStore => {
            let store = RowStore::open(path).map_err(|e| input(format!("'{shown}': {e}")))?;
            LazyFrame::scan(store).collect().map_err(|e| input(format!("'{shown}': {e}")))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn scratch_file(name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("tptdyn-table-test-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        fs::write(&path, bytes).unwrap();
        path
    }

    #[test]
    fn detects_formats_by_magic_bytes() {
        assert_eq!(detect(&scratch_file("a.bin", b"PAR1....")).unwrap(), Kind::Parquet);
        assert_eq!(detect(&scratch_file("b.bin", b"TPTROWS1xxxx")).unwrap(), Kind::RowStore);
        assert_eq!(detect(&scratch_file("c.bin", b"TPTC\x01\x00")).unwrap(), Kind::Columnar);
        assert_eq!(detect(&scratch_file("d.csv", b"x,y\n1,2\n")).unwrap(), Kind::Csv);
        assert_eq!(detect(&scratch_file("e.csv", b"")).unwrap(), Kind::Csv);
    }

    #[test]
    fn trajectories_are_redirected() {
        let err = detect(&scratch_file("t1.bin", b"TPTCTRJ1rest")).unwrap_err();
        assert!(matches!(err, CliError::Input(ref m) if m.contains("--trajectory")), "{err}");
        assert!(detect(&scratch_file("t2.bin", b"TPTTRJ01rest")).is_err());
    }

    #[test]
    fn columnar_round_trip_loads_as_frame() {
        let path =
            std::env::temp_dir().join(format!("tptdyn-table-columnar-{}.tpc", std::process::id()));
        let file = File::create(&path).unwrap();
        let cols = vec![
            ("x".to_string(), ColumnData::F64(vec![0.5, 1.5, 2.5])),
            ("n".to_string(), ColumnData::U32(vec![1, 2, 3])),
        ];
        columnar::write_table(file, &cols).unwrap();
        assert_eq!(detect(&path).unwrap(), Kind::Columnar);
        let frame = load(&path, Kind::Columnar).unwrap();
        assert_eq!(frame.len(), 3);
        assert_eq!(frame.column("n"), Some(&Column::I64(vec![1, 2, 3])));
        assert_eq!(frame.column("x"), Some(&Column::F64(vec![0.5, 1.5, 2.5])));
    }
}
