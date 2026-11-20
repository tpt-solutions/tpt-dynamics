//! Custom binary columnar format for detector data (version 1).
//!
//! Layout, all integers little-endian:
//!
//! ```text
//! header:
//!   magic      4 bytes  b"TPTC"
//!   version    u16      1
//!   flags      u16      0 (reserved; non-zero is rejected)
//!   ncols      u32
//!   nrows      u64
//!   per column:
//!     name_len u16, name (UTF-8, unique within the file)
//!     dtype    u8      1 = f64, 2 = i64, 3 = u32
//!   header_crc u32      CRC-32 of every header byte before it
//! per column, in header order:
//!   data       nrows × element width, little-endian values
//!   crc        u32      CRC-32 of the data bytes
//! ```
//!
//! Each column is one contiguous block, so a reader that wants one column can
//! skip the others by their known sizes. Compression is not part of version 1.
//! Every read path bounds its allocations by the bytes actually present, so a
//! corrupted row count cannot trigger a huge allocation. Corruption and truncation
//! return a [`ColumnarError`]; nothing panics.
//!
//! ```
//! use tpt_dynamics_data::columnar::{write_table, read_table, ColumnData, DType, Table};
//!
//! let cols = vec![
//!     ("pt".to_string(), ColumnData::F64(vec![1.5, 2.5])),
//!     ("event".to_string(), ColumnData::I64(vec![7, 8])),
//! ];
//! let mut bytes = Vec::new();
//! write_table(&mut bytes, &cols).unwrap();
//! let table: Table = read_table(&bytes[..]).unwrap();
//! assert_eq!(table.nrows(), 2);
//! assert_eq!(table.column("event").unwrap().dtype(), DType::I64);
//! ```

use crate::crc32::checksum;
use core::fmt;
use std::io::{self, Read, Write};

/// File magic.
pub const MAGIC: &[u8; 4] = b"TPTC";
/// Format version written by this crate.
pub const VERSION: u16 = 1;

/// Element type of a column.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DType {
    /// 64-bit IEEE float.
    F64,
    /// 64-bit signed integer.
    I64,
    /// 32-bit unsigned integer.
    U32,
}

impl DType {
    /// Size of one element in bytes.
    pub fn width(self) -> usize {
        match self {
            Self::F64 | Self::I64 => 8,
            Self::U32 => 4,
        }
    }

    fn code(self) -> u8 {
        match self {
            Self::F64 => 1,
            Self::I64 => 2,
            Self::U32 => 3,
        }
    }

    pub(crate) fn from_code(c: u8) -> Option<Self> {
        match c {
            1 => Some(Self::F64),
            2 => Some(Self::I64),
            3 => Some(Self::U32),
            _ => None,
        }
    }
}

/// One column of values.
#[derive(Clone, Debug, PartialEq)]
pub enum ColumnData {
    /// `f64` values.
    F64(Vec<f64>),
    /// `i64` values.
    I64(Vec<i64>),
    /// `u32` values.
    U32(Vec<u32>),
}

impl ColumnData {
    /// Element type.
    pub fn dtype(&self) -> DType {
        match self {
            Self::F64(_) => DType::F64,
            Self::I64(_) => DType::I64,
            Self::U32(_) => DType::U32,
        }
    }

    /// Number of values.
    pub fn len(&self) -> usize {
        match self {
            Self::F64(v) => v.len(),
            Self::I64(v) => v.len(),
            Self::U32(v) => v.len(),
        }
    }

    /// Whether the column is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Values as `f64` if the column is `F64`.
    pub fn as_f64(&self) -> Option<&[f64]> {
        match self {
            Self::F64(v) => Some(v),
            _ => None,
        }
    }

    /// Values as `i64` if the column is `I64`.
    pub fn as_i64(&self) -> Option<&[i64]> {
        match self {
            Self::I64(v) => Some(v),
            _ => None,
        }
    }

    /// Values as `u32` if the column is `U32`.
    pub fn as_u32(&self) -> Option<&[u32]> {
        match self {
            Self::U32(v) => Some(v),
            _ => None,
        }
    }

    fn encode_into(&self, out: &mut Vec<u8>) {
        match self {
            Self::F64(v) => v.iter().for_each(|x| out.extend_from_slice(&x.to_le_bytes())),
            Self::I64(v) => v.iter().for_each(|x| out.extend_from_slice(&x.to_le_bytes())),
            Self::U32(v) => v.iter().for_each(|x| out.extend_from_slice(&x.to_le_bytes())),
        }
    }

    fn decode(dtype: DType, bytes: &[u8]) -> Self {
        // `as_chunks` yields exact-width arrays; a trailing partial chunk is ignored.
        match dtype {
            DType::F64 => {
                Self::F64(bytes.as_chunks::<8>().0.iter().map(|c| f64::from_le_bytes(*c)).collect())
            }
            DType::I64 => {
                Self::I64(bytes.as_chunks::<8>().0.iter().map(|c| i64::from_le_bytes(*c)).collect())
            }
            DType::U32 => {
                Self::U32(bytes.as_chunks::<4>().0.iter().map(|c| u32::from_le_bytes(*c)).collect())
            }
        }
    }
}

/// Errors from writing or reading the columnar format.
#[derive(Debug)]
pub enum ColumnarError {
    /// Underlying I/O failure.
    Io(io::Error),
    /// The input ended before the expected data.
    Truncated,
    /// The file does not start with [`MAGIC`].
    BadMagic,
    /// The version is not supported by this build.
    UnsupportedVersion(u16),
    /// Reserved header flags are set.
    ReservedFlags(u16),
    /// A column has an unknown type code.
    UnknownDtype(u8),
    /// A column name is not valid UTF-8, is empty, or is duplicated.
    BadColumnName(String),
    /// A checksum did not match. `what` names the header or column.
    ChecksumMismatch {
        /// Which part failed verification.
        what: String,
    },
    /// Bytes follow the last column block.
    TrailingData,
    /// The writer was misused (wrong column order, type, or length).
    InvalidInput(&'static str),
}

impl fmt::Display for ColumnarError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Truncated => write!(f, "columnar data is truncated"),
            Self::BadMagic => write!(f, "not a tpt-dynamics columnar file (bad magic)"),
            Self::UnsupportedVersion(v) => write!(f, "unsupported columnar format version {v}"),
            Self::ReservedFlags(x) => write!(f, "reserved header flags set: {x:#06x}"),
            Self::UnknownDtype(c) => write!(f, "unknown column type code {c}"),
            Self::BadColumnName(n) => write!(f, "invalid or duplicate column name {n:?}"),
            Self::ChecksumMismatch { what } => {
                write!(f, "checksum mismatch in {what} (data corrupted)")
            }
            Self::TrailingData => write!(f, "unexpected data after the last column"),
            Self::InvalidInput(m) => write!(f, "invalid columnar write: {m}"),
        }
    }
}

impl std::error::Error for ColumnarError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for ColumnarError {
    fn from(e: io::Error) -> Self {
        if e.kind() == io::ErrorKind::UnexpectedEof {
            Self::Truncated
        } else {
            Self::Io(e)
        }
    }
}

/// Schema entry for one column: name and type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ColumnSpec {
    /// Column name.
    pub name: String,
    /// Element type.
    pub dtype: DType,
}

/// Streaming writer. The header is written by [`ColumnarWriter::new`]; columns
/// must then be written in schema order with [`ColumnarWriter::write_column`].
#[derive(Debug)]
pub struct ColumnarWriter<W: Write> {
    inner: W,
    specs: Vec<ColumnSpec>,
    nrows: u64,
    next: usize,
}

impl<W: Write> ColumnarWriter<W> {
    /// Validate the schema and write the header.
    pub fn new(mut inner: W, specs: Vec<ColumnSpec>, nrows: u64) -> Result<Self, ColumnarError> {
        if specs.len() > u32::MAX as usize {
            return Err(ColumnarError::InvalidInput("too many columns"));
        }
        let mut seen = std::collections::HashSet::new();
        for s in &specs {
            if s.name.is_empty()
                || s.name.len() > u16::MAX as usize
                || !seen.insert(s.name.as_str())
            {
                return Err(ColumnarError::BadColumnName(s.name.clone()));
            }
        }
        let mut header = Vec::new();
        header.extend_from_slice(MAGIC);
        header.extend_from_slice(&VERSION.to_le_bytes());
        header.extend_from_slice(&0u16.to_le_bytes());
        header.extend_from_slice(&(specs.len() as u32).to_le_bytes());
        header.extend_from_slice(&nrows.to_le_bytes());
        for s in &specs {
            header.extend_from_slice(&(s.name.len() as u16).to_le_bytes());
            header.extend_from_slice(s.name.as_bytes());
            header.push(s.dtype.code());
        }
        let crc = checksum(&header);
        header.extend_from_slice(&crc.to_le_bytes());
        inner.write_all(&header)?;
        Ok(Self { inner, specs, nrows, next: 0 })
    }

    /// Write the next column block. Its type and length must match the schema.
    pub fn write_column(&mut self, data: &ColumnData) -> Result<(), ColumnarError> {
        let spec = self
            .specs
            .get(self.next)
            .ok_or(ColumnarError::InvalidInput("all columns already written"))?;
        if spec.dtype != data.dtype() {
            return Err(ColumnarError::InvalidInput("column type does not match schema"));
        }
        if data.len() as u64 != self.nrows {
            return Err(ColumnarError::InvalidInput("column length does not match row count"));
        }
        let mut bytes = Vec::with_capacity(data.len() * data.dtype().width());
        data.encode_into(&mut bytes);
        self.inner.write_all(&bytes)?;
        self.inner.write_all(&checksum(&bytes).to_le_bytes())?;
        self.next += 1;
        Ok(())
    }

    /// Check that every column was written, flush, and return the underlying writer.
    pub fn finish(mut self) -> Result<W, ColumnarError> {
        if self.next != self.specs.len() {
            return Err(ColumnarError::InvalidInput("not all columns were written"));
        }
        self.inner.flush()?;
        Ok(self.inner)
    }
}

/// Write a whole table in one call. Column names must be unique and non-empty;
/// all columns must have the same length.
pub fn write_table<W: Write>(w: W, columns: &[(String, ColumnData)]) -> Result<W, ColumnarError> {
    let nrows = columns.first().map_or(0, |(_, c)| c.len() as u64);
    let specs = columns
        .iter()
        .map(|(name, c)| ColumnSpec { name: name.clone(), dtype: c.dtype() })
        .collect();
    let mut writer = ColumnarWriter::new(w, specs, nrows)?;
    for (_, c) in columns {
        writer.write_column(c)?;
    }
    writer.finish()
}

/// A named column read from a file.
#[derive(Clone, Debug, PartialEq)]
pub struct Column {
    /// Column name.
    pub name: String,
    /// Values.
    pub data: ColumnData,
}

/// A fully read and verified table.
#[derive(Clone, Debug, PartialEq)]
pub struct Table {
    nrows: u64,
    columns: Vec<Column>,
}

impl Table {
    /// Number of rows.
    pub fn nrows(&self) -> u64 {
        self.nrows
    }

    /// Columns in file order.
    pub fn columns(&self) -> &[Column] {
        &self.columns
    }

    /// Look up a column by name.
    pub fn column(&self, name: &str) -> Option<&ColumnData> {
        self.columns.iter().find(|c| c.name == name).map(|c| &c.data)
    }
}

fn read_exact_arr<const N: usize, R: Read>(r: &mut R) -> Result<[u8; N], ColumnarError> {
    let mut buf = [0u8; N];
    r.read_exact(&mut buf)?;
    Ok(buf)
}

/// Read and verify a complete columnar file.
///
/// Every header and column checksum is verified before the data is returned.
/// Allocation for a column is proportional to the bytes actually present in the
/// input, so a corrupt row count cannot exhaust memory.
pub fn read_table<R: Read>(mut r: R) -> Result<Table, ColumnarError> {
    let mut header = Vec::new();
    let magic: [u8; 4] = read_exact_arr(&mut r)?;
    header.extend_from_slice(&magic);
    if &magic != MAGIC {
        return Err(ColumnarError::BadMagic);
    }
    let version_bytes: [u8; 2] = read_exact_arr(&mut r)?;
    header.extend_from_slice(&version_bytes);
    let version = u16::from_le_bytes(version_bytes);
    if version != VERSION {
        return Err(ColumnarError::UnsupportedVersion(version));
    }
    let flags_bytes: [u8; 2] = read_exact_arr(&mut r)?;
    header.extend_from_slice(&flags_bytes);
    let flags = u16::from_le_bytes(flags_bytes);
    if flags != 0 {
        return Err(ColumnarError::ReservedFlags(flags));
    }
    let ncols_bytes: [u8; 4] = read_exact_arr(&mut r)?;
    header.extend_from_slice(&ncols_bytes);
    let ncols = u32::from_le_bytes(ncols_bytes) as usize;
    let nrows_bytes: [u8; 8] = read_exact_arr(&mut r)?;
    header.extend_from_slice(&nrows_bytes);
    let nrows = u64::from_le_bytes(nrows_bytes);

    let mut specs = Vec::with_capacity(ncols.min(1024));
    for _ in 0..ncols {
        let len_bytes: [u8; 2] = read_exact_arr(&mut r)?;
        header.extend_from_slice(&len_bytes);
        let len = u16::from_le_bytes(len_bytes) as usize;
        let mut name_bytes = vec![0u8; len];
        r.read_exact(&mut name_bytes)?;
        header.extend_from_slice(&name_bytes);
        let name = String::from_utf8(name_bytes.clone()).map_err(|_| {
            ColumnarError::BadColumnName(String::from_utf8_lossy(&name_bytes).into_owned())
        })?;
        let code: [u8; 1] = read_exact_arr(&mut r)?;
        header.extend_from_slice(&code);
        let dtype = DType::from_code(code[0]).ok_or(ColumnarError::UnknownDtype(code[0]))?;
        specs.push(ColumnSpec { name, dtype });
    }
    let stored_header_crc = u32::from_le_bytes(read_exact_arr(&mut r)?);
    if checksum(&header) != stored_header_crc {
        return Err(ColumnarError::ChecksumMismatch { what: "header".into() });
    }
    for (i, a) in specs.iter().enumerate() {
        if a.name.is_empty() || specs[..i].iter().any(|b| b.name == a.name) {
            return Err(ColumnarError::BadColumnName(a.name.clone()));
        }
    }

    let mut columns = Vec::with_capacity(specs.len().min(1024));
    for spec in specs {
        let nbytes =
            nrows.checked_mul(spec.dtype.width() as u64).ok_or(ColumnarError::Truncated)?;
        let mut data = Vec::new();
        // `take` + `read_to_end` grows the buffer only as bytes arrive.
        (&mut r).take(nbytes).read_to_end(&mut data)?;
        if data.len() as u64 != nbytes {
            return Err(ColumnarError::Truncated);
        }
        let stored = u32::from_le_bytes(read_exact_arr(&mut r)?);
        if checksum(&data) != stored {
            return Err(ColumnarError::ChecksumMismatch {
                what: format!("column \"{}\"", spec.name),
            });
        }
        columns.push(Column { name: spec.name, data: ColumnData::decode(spec.dtype, &data) });
    }

    let mut probe = [0u8; 1];
    if r.read(&mut probe)? != 0 {
        return Err(ColumnarError::TrailingData);
    }
    Ok(Table { nrows, columns })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::SplitMix64;

    fn sample() -> Vec<(String, ColumnData)> {
        vec![
            ("pt".into(), ColumnData::F64((0..50).map(|i| f64::from(i) * 0.5 - 3.0).collect())),
            ("event".into(), ColumnData::I64((0..50).map(|i| i64::from(i) * 1000 - 7).collect())),
            ("ntrk".into(), ColumnData::U32((0..50).map(|i| i % 9).collect())),
        ]
    }

    fn bytes_of(cols: &[(String, ColumnData)]) -> Vec<u8> {
        let mut out = Vec::new();
        write_table(&mut out, cols).unwrap();
        out
    }

    #[test]
    fn round_trip() {
        let cols = sample();
        let table = read_table(&bytes_of(&cols)[..]).unwrap();
        assert_eq!(table.nrows(), 50);
        assert_eq!(table.columns().len(), 3);
        for (name, data) in &cols {
            assert_eq!(table.column(name).unwrap(), data);
        }
    }

    #[test]
    fn empty_table_round_trip() {
        let cols = vec![("x".to_string(), ColumnData::F64(vec![]))];
        let table = read_table(&bytes_of(&cols)[..]).unwrap();
        assert_eq!(table.nrows(), 0);
        assert_eq!(table.column("x"), Some(&ColumnData::F64(vec![])));
    }

    #[test]
    fn nan_and_extremes_survive() {
        let cols =
            vec![("x".to_string(), ColumnData::F64(vec![f64::NAN, f64::MIN_POSITIVE, -0.0]))];
        let table = read_table(&bytes_of(&cols)[..]).unwrap();
        let v = table.column("x").unwrap().as_f64().unwrap();
        assert!(v[0].is_nan());
        assert_eq!(v[1], f64::MIN_POSITIVE);
        assert!(v[2].is_sign_negative());
    }

    #[test]
    fn writer_enforces_schema() {
        let specs = vec![ColumnSpec { name: "a".into(), dtype: DType::F64 }];
        let mut w = ColumnarWriter::new(Vec::new(), specs.clone(), 2).unwrap();
        assert!(w.write_column(&ColumnData::I64(vec![1, 2])).is_err());
        assert!(w.write_column(&ColumnData::F64(vec![1.0])).is_err());
        let w2 = ColumnarWriter::new(Vec::new(), specs.clone(), 2).unwrap();
        assert!(w2.finish().is_err());
        assert!(
            ColumnarWriter::new(Vec::new(), vec![specs[0].clone(), specs[0].clone()], 0).is_err()
        );
    }

    #[test]
    fn every_single_byte_corruption_is_detected() {
        let bytes = bytes_of(&sample());
        for i in 0..bytes.len() {
            let mut bad = bytes.clone();
            bad[i] ^= 0x01;
            assert!(read_table(&bad[..]).is_err(), "corruption at byte {i} accepted");
        }
    }

    #[test]
    fn every_truncation_is_rejected() {
        let bytes = bytes_of(&sample());
        for n in 0..bytes.len() {
            assert!(read_table(&bytes[..n]).is_err(), "prefix of {n} bytes accepted");
        }
    }

    #[test]
    fn trailing_data_is_rejected() {
        let mut bytes = bytes_of(&sample());
        bytes.push(0);
        assert!(matches!(read_table(&bytes[..]), Err(ColumnarError::TrailingData)));
    }

    #[test]
    fn huge_row_count_does_not_allocate() {
        // A header claiming 2^60 rows with a valid header checksum must fail with
        // Truncated, not attempt a giant allocation.
        let mut header = Vec::new();
        header.extend_from_slice(MAGIC);
        header.extend_from_slice(&VERSION.to_le_bytes());
        header.extend_from_slice(&0u16.to_le_bytes());
        header.extend_from_slice(&1u32.to_le_bytes());
        header.extend_from_slice(&(1u64 << 60).to_le_bytes());
        header.extend_from_slice(&1u16.to_le_bytes());
        header.push(b'x');
        header.push(DType::F64.code());
        let crc = checksum(&header);
        header.extend_from_slice(&crc.to_le_bytes());
        header.extend_from_slice(&[1, 2, 3]);
        assert!(matches!(read_table(&header[..]), Err(ColumnarError::Truncated)));
    }

    #[test]
    fn random_garbage_never_panics() {
        let mut rng = SplitMix64::new(99);
        for _ in 0..2000 {
            let len = rng.below(200) as usize;
            let mut data: Vec<u8> = (0..len).map(|_| rng.next_u64() as u8).collect();
            // Sometimes give the garbage a valid magic and version to go deeper.
            if len >= 6 && rng.below(2) == 0 {
                data[..4].copy_from_slice(MAGIC);
                data[4..6].copy_from_slice(&VERSION.to_le_bytes());
            }
            let _ = read_table(&data[..]);
        }
    }
}
