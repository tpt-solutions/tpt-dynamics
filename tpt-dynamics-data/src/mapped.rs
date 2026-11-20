//! Zero-copy, memory-mapped reader for the version-1 columnar format.
//!
//! [`MappedTable::open`] maps the file read-only and verifies every header and
//! column checksum, without copying column data into the heap. Column values are
//! then decoded from the mapped bytes as they are iterated, so a scan touches the
//! pages it needs and nothing else.
//!
//! This module is the only place in the crate that uses `unsafe` code. The single
//! unsafe call is the memory mapping itself. See ADR-011 for the exception and its
//! invariant, which is that the file must not be modified while it is mapped.
//!
//! ```no_run
//! use tpt_dynamics_data::mapped::MappedTable;
//!
//! let table = MappedTable::open("events.tptc")?;
//! let total_pt: f64 = table.f64s("pt").expect("pt column").sum();
//! # Ok::<(), tpt_dynamics_data::columnar::ColumnarError>(())
//! ```

use std::fs::File;
use std::io;
use std::path::Path;

use memmap2::Mmap;

use crate::columnar::{ColumnarError, DType, MAGIC, VERSION};
use crate::crc32::checksum;

/// Where one column sits inside the mapped file.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Span {
    name: String,
    dtype: DType,
    start: usize,
    len: usize,
}

/// A columnar file opened as a read-only memory mapping.
#[derive(Debug)]
pub struct MappedTable {
    map: Mmap,
    nrows: u64,
    spans: Vec<Span>,
}

impl MappedTable {
    /// Map `path` and verify all header and column checksums.
    ///
    /// Verification reads every byte once but allocates nothing for column data.
    /// Any corruption or truncation is reported as a [`ColumnarError`].
    pub fn open(path: impl AsRef<Path>) -> Result<Self, ColumnarError> {
        let file = File::open(path)?;
        if file.metadata()?.len() == 0 {
            return Err(ColumnarError::Truncated);
        }
        let map = map_file(&file)?;
        let parsed = parse(&map)?;
        Ok(Self { map, nrows: parsed.nrows, spans: parsed.spans })
    }

    /// Number of rows in every column.
    pub fn nrows(&self) -> u64 {
        self.nrows
    }

    /// Column names and types, in file order.
    pub fn columns(&self) -> impl Iterator<Item = (&str, DType)> + '_ {
        self.spans.iter().map(|s| (s.name.as_str(), s.dtype))
    }

    /// Type of the named column, if present.
    pub fn dtype(&self, name: &str) -> Option<DType> {
        self.span(name).map(|s| s.dtype)
    }

    /// Raw little-endian bytes of the named column, borrowed from the mapping.
    pub fn bytes(&self, name: &str) -> Option<&[u8]> {
        self.span(name).map(|s| &self.map[s.start..s.start + s.len])
    }

    /// Values of an `f64` column, decoded lazily from the mapping.
    ///
    /// Returns `None` if the column is missing or has a different type.
    pub fn f64s(&self, name: &str) -> Option<impl Iterator<Item = f64> + '_> {
        self.typed(name, DType::F64)
            .map(|b| b.as_chunks::<8>().0.iter().map(|c| f64::from_le_bytes(*c)))
    }

    /// Values of an `i64` column, decoded lazily from the mapping.
    ///
    /// Returns `None` if the column is missing or has a different type.
    pub fn i64s(&self, name: &str) -> Option<impl Iterator<Item = i64> + '_> {
        self.typed(name, DType::I64)
            .map(|b| b.as_chunks::<8>().0.iter().map(|c| i64::from_le_bytes(*c)))
    }

    /// Values of a `u32` column, decoded lazily from the mapping.
    ///
    /// Returns `None` if the column is missing or has a different type.
    pub fn u32s(&self, name: &str) -> Option<impl Iterator<Item = u32> + '_> {
        self.typed(name, DType::U32)
            .map(|b| b.as_chunks::<4>().0.iter().map(|c| u32::from_le_bytes(*c)))
    }

    fn span(&self, name: &str) -> Option<&Span> {
        self.spans.iter().find(|s| s.name == name)
    }

    fn typed(&self, name: &str, want: DType) -> Option<&[u8]> {
        let s = self.span(name)?;
        (s.dtype == want).then(|| &self.map[s.start..s.start + s.len])
    }
}

/// Map a file read-only.
///
/// SAFETY: the mapping is never written through and is only read via shared slices.
/// Rust's soundness here depends on the file not being truncated or rewritten by any
/// process while it is mapped. Violating that is undefined behaviour, so this reader
/// must only be pointed at files that are complete and not being written. ADR-011 records
/// the exception and the policy for callers.
#[allow(unsafe_code)]
fn map_file(file: &File) -> io::Result<Mmap> {
    // SAFETY: see the function documentation above.
    unsafe { Mmap::map(file) }
}

struct Parsed {
    nrows: u64,
    spans: Vec<Span>,
}

/// A bounds-checked cursor over the mapped bytes.
struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], ColumnarError> {
        let end = self.pos.checked_add(n).ok_or(ColumnarError::Truncated)?;
        let out = self.buf.get(self.pos..end).ok_or(ColumnarError::Truncated)?;
        self.pos = end;
        Ok(out)
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N], ColumnarError> {
        self.take(N)?.try_into().map_err(|_| ColumnarError::Truncated)
    }
}

/// Parse and verify a whole file held in `buf`, recording where each column lies.
fn parse(buf: &[u8]) -> Result<Parsed, ColumnarError> {
    let mut c = Cursor { buf, pos: 0 };
    if c.take(4)? != MAGIC {
        return Err(ColumnarError::BadMagic);
    }
    let version = u16::from_le_bytes(c.array()?);
    if version != VERSION {
        return Err(ColumnarError::UnsupportedVersion(version));
    }
    let flags = u16::from_le_bytes(c.array()?);
    if flags != 0 {
        return Err(ColumnarError::ReservedFlags(flags));
    }
    let ncols = u32::from_le_bytes(c.array()?) as usize;
    let nrows = u64::from_le_bytes(c.array()?);

    let mut specs: Vec<(String, DType)> = Vec::with_capacity(ncols.min(1024));
    for _ in 0..ncols {
        let len = u16::from_le_bytes(c.array()?) as usize;
        let name_bytes = c.take(len)?;
        let name = std::str::from_utf8(name_bytes)
            .map_err(|_| ColumnarError::BadColumnName(String::from_utf8_lossy(name_bytes).into()))?
            .to_owned();
        let code = c.array::<1>()?[0];
        let dtype = DType::from_code(code).ok_or(ColumnarError::UnknownDtype(code))?;
        specs.push((name, dtype));
    }

    let header_end = c.pos;
    let stored_header_crc = u32::from_le_bytes(c.array()?);
    if checksum(&buf[..header_end]) != stored_header_crc {
        return Err(ColumnarError::ChecksumMismatch { what: "header".into() });
    }
    for (i, (name, _)) in specs.iter().enumerate() {
        if name.is_empty() || specs[..i].iter().any(|(other, _)| other == name) {
            return Err(ColumnarError::BadColumnName(name.clone()));
        }
    }

    let mut spans = Vec::with_capacity(specs.len().min(1024));
    for (name, dtype) in specs {
        let nbytes = nrows
            .checked_mul(dtype.width() as u64)
            .and_then(|n| usize::try_from(n).ok())
            .ok_or(ColumnarError::Truncated)?;
        let start = c.pos;
        let data = c.take(nbytes)?;
        let stored = u32::from_le_bytes(c.array()?);
        if checksum(data) != stored {
            return Err(ColumnarError::ChecksumMismatch { what: format!("column \"{name}\"") });
        }
        spans.push(Span { name, dtype, start, len: nbytes });
    }

    if c.pos != buf.len() {
        return Err(ColumnarError::TrailingData);
    }
    Ok(Parsed { nrows, spans })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::columnar::{write_table, ColumnData};
    use std::path::PathBuf;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("tptdyn-mapped-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch dir");
        dir.join(name)
    }

    fn sample_file(name: &str) -> (PathBuf, Vec<(String, ColumnData)>) {
        let cols = vec![
            (
                "pt".to_string(),
                ColumnData::F64((0..50).map(|i| f64::from(i) * 0.5 - 3.0).collect()),
            ),
            (
                "event".to_string(),
                ColumnData::I64((0..50).map(|i| i64::from(i) * 1000 - 7).collect()),
            ),
            ("ntrk".to_string(), ColumnData::U32((0..50).map(|i| i % 9).collect())),
        ];
        let path = scratch(name);
        let file = File::create(&path).expect("create");
        write_table(file, &cols).expect("write");
        (path, cols)
    }

    #[test]
    fn mapped_values_match_the_written_columns() {
        let (path, cols) = sample_file("match.tptc");
        let t = MappedTable::open(&path).expect("open");
        assert_eq!(t.nrows(), 50);
        let pt: Vec<f64> = t.f64s("pt").expect("pt").collect();
        assert_eq!(Some(&pt[..]), cols[0].1.as_f64());
        let ev: Vec<i64> = t.i64s("event").expect("event").collect();
        assert_eq!(Some(&ev[..]), cols[1].1.as_i64());
        let nt: Vec<u32> = t.u32s("ntrk").expect("ntrk").collect();
        assert_eq!(Some(&nt[..]), cols[2].1.as_u32());
        assert_eq!(t.bytes("pt").map(<[u8]>::len), Some(50 * 8));
    }

    #[test]
    fn wrong_type_or_missing_column_gives_none() {
        let (path, _) = sample_file("types.tptc");
        let t = MappedTable::open(&path).expect("open");
        assert!(t.f64s("event").is_none());
        assert!(t.i64s("pt").is_none());
        assert!(t.f64s("absent").is_none());
        assert_eq!(t.dtype("ntrk"), Some(DType::U32));
    }

    #[test]
    fn columns_are_listed_in_file_order() {
        let (path, _) = sample_file("order.tptc");
        let t = MappedTable::open(&path).expect("open");
        let names: Vec<&str> = t.columns().map(|(n, _)| n).collect();
        assert_eq!(names, ["pt", "event", "ntrk"]);
    }

    #[test]
    fn corrupted_column_is_rejected() {
        let (path, _) = sample_file("corrupt.tptc");
        let mut bytes = std::fs::read(&path).expect("read");
        let mid = bytes.len() / 2;
        bytes[mid] ^= 0x55;
        std::fs::write(&path, &bytes).expect("write back");
        assert!(matches!(MappedTable::open(&path), Err(ColumnarError::ChecksumMismatch { .. })));
    }

    #[test]
    fn truncated_file_is_rejected() {
        let (path, _) = sample_file("trunc.tptc");
        let bytes = std::fs::read(&path).expect("read");
        std::fs::write(&path, &bytes[..bytes.len() - 7]).expect("write");
        assert!(matches!(MappedTable::open(&path), Err(ColumnarError::Truncated)));
    }

    #[test]
    fn empty_and_foreign_files_are_rejected() {
        let empty = scratch("empty.tptc");
        std::fs::write(&empty, b"").expect("write");
        assert!(matches!(MappedTable::open(&empty), Err(ColumnarError::Truncated)));
        let foreign = scratch("foreign.tptc");
        std::fs::write(&foreign, b"NOT A COLUMNAR FILE AT ALL").expect("write");
        assert!(matches!(MappedTable::open(&foreign), Err(ColumnarError::BadMagic)));
    }

    #[test]
    fn trailing_bytes_are_rejected() {
        let (path, _) = sample_file("trail.tptc");
        let mut bytes = std::fs::read(&path).expect("read");
        bytes.extend_from_slice(b"junk");
        std::fs::write(&path, &bytes).expect("write");
        assert!(matches!(MappedTable::open(&path), Err(ColumnarError::TrailingData)));
    }
}
