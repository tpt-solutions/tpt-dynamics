//! Arrow-IPC-style columnar stream reader (own implementation, read-only).
//!
//! **This format is NOT Arrow-compatible.** It borrows the Arrow *shape*
//! (a schema followed by record batches, each column stored as a values buffer
//! plus an optional validity bitmap) but uses its own byte layout with no
//! `FlatBuffers` metadata. Files written by `pyarrow` or any other Arrow library
//! will be rejected with [`ArrowIpcError::BadMagic`]. The decision and its
//! reasons are recorded in `todo.md` (Phase 7) and `docs/audit/security-review.md`.
//!
//! The format supports three primitive column types: `int64`, `float64` and
//! `bool`. Each column may be nullable.
//!
//! # Layout (all integers little-endian)
//!
//! ```text
//! header:   b"TPTA" | u16 version (= 1) | u16 flags (must be 0)
//! schema:   u16 ncols | ncols * (u16 name_len | name (UTF-8) | u8 dtype | u8 nullable)
//!           | u32 crc32 of the schema bytes (from `ncols` through the last field)
//!           dtype: 1 = int64, 2 = float64, 3 = bool
//! messages: repeated, each starting with a u8 tag:
//!           1 = record batch, 0 = end of stream (must be last)
//! batch:    u32 nrows | per column:
//!               [validity bitmap, ceil(nrows/8) bytes, bit i = row i valid]  if nullable
//!               values: nrows * 8 bytes (int64/float64), or ceil(nrows/8) bytes bitmap (bool)
//!           | u32 crc32 of the batch bytes (after the tag, excluding this checksum)
//! ```
//!
//! Limits: at most [`MAX_COLUMNS`] columns, [`MAX_NAME_BYTES`] bytes per name, and
//! [`MAX_BATCH_ROWS`] rows per batch. Every length is checked against the bytes
//! actually remaining in the input before any allocation, so a corrupt length
//! field cannot trigger a large allocation.
//!
//! ```
//! # use tpt_dynamics_data::arrow_ipc::{read_all, ArrowIpcError};
//! // An empty stream (header + empty schema + end marker) has zero batches.
//! # let bytes = tpt_dynamics_data::arrow_ipc::test_support::empty_stream();
//! let batches = read_all(&bytes)?;
//! assert!(batches.is_empty());
//! # Ok::<(), ArrowIpcError>(())
//! ```

use core::fmt;

use crate::crc32;

/// Magic bytes at the start of every stream.
pub const MAGIC: &[u8; 4] = b"TPTA";
/// Format version understood by this reader.
pub const VERSION: u16 = 1;
/// Maximum number of columns in a schema.
pub const MAX_COLUMNS: usize = 1024;
/// Maximum UTF-8 byte length of a column name.
pub const MAX_NAME_BYTES: usize = 256;
/// Maximum number of rows in one record batch.
pub const MAX_BATCH_ROWS: usize = 1 << 24;

const TAG_END: u8 = 0;
const TAG_BATCH: u8 = 1;

/// Primitive column types supported by this format.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DataType {
    /// Signed 64-bit integer.
    Int64,
    /// IEEE 754 double.
    Float64,
    /// Boolean (bit-packed on disk).
    Bool,
}

impl DataType {
    fn from_code(code: u8) -> Option<Self> {
        match code {
            1 => Some(Self::Int64),
            2 => Some(Self::Float64),
            3 => Some(Self::Bool),
            _ => None,
        }
    }

    /// The on-disk type code.
    pub fn code(self) -> u8 {
        match self {
            Self::Int64 => 1,
            Self::Float64 => 2,
            Self::Bool => 3,
        }
    }
}

/// One field of the schema.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Field {
    /// Column name.
    pub name: String,
    /// Column type.
    pub dtype: DataType,
    /// Whether the column carries a validity bitmap.
    pub nullable: bool,
}

/// Values of one column in one batch.
#[derive(Clone, Debug, PartialEq)]
pub enum Values {
    /// `int64` values.
    Int64(Vec<i64>),
    /// `float64` values.
    Float64(Vec<f64>),
    /// `bool` values.
    Bool(Vec<bool>),
}

impl Values {
    /// Number of values.
    pub fn len(&self) -> usize {
        match self {
            Self::Int64(v) => v.len(),
            Self::Float64(v) => v.len(),
            Self::Bool(v) => v.len(),
        }
    }

    /// Whether there are no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// One column of a record batch.
#[derive(Clone, Debug, PartialEq)]
pub struct ColumnArray {
    /// Values; for null slots the value is unspecified (zero on write).
    pub values: Values,
    /// `Some(bits)` for nullable columns: `bits[i]` is true when row `i` is valid.
    pub validity: Option<Vec<bool>>,
}

impl ColumnArray {
    /// Number of null rows (zero for non-nullable columns).
    pub fn null_count(&self) -> usize {
        self.validity.as_ref().map_or(0, |v| v.iter().filter(|b| !**b).count())
    }
}

/// A record batch: one array per schema field, all of equal length.
#[derive(Clone, Debug, PartialEq)]
pub struct RecordBatch {
    nrows: usize,
    columns: Vec<ColumnArray>,
}

impl RecordBatch {
    /// Number of rows.
    pub fn nrows(&self) -> usize {
        self.nrows
    }

    /// Columns in schema order.
    pub fn columns(&self) -> &[ColumnArray] {
        &self.columns
    }
}

/// Errors from reading an Arrow-IPC-style stream.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ArrowIpcError {
    /// The input ended before the expected data, or the end marker is missing.
    Truncated,
    /// The input does not start with [`MAGIC`].
    BadMagic,
    /// The version is not supported.
    UnsupportedVersion(u16),
    /// Reserved header flags are set.
    ReservedFlags(u16),
    /// A column has an unknown type code.
    UnknownDtype(u8),
    /// A column name is empty, too long, not UTF-8, or duplicated.
    BadName(String),
    /// A field exceeds a documented limit.
    TooLarge {
        /// What is too large.
        what: &'static str,
        /// Value found.
        value: usize,
        /// Limit.
        limit: usize,
    },
    /// A checksum did not match.
    ChecksumMismatch(&'static str),
    /// An unknown message tag was found.
    BadTag(u8),
    /// Bytes follow the end-of-stream marker.
    TrailingData,
}

impl fmt::Display for ArrowIpcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Truncated => write!(f, "stream is truncated"),
            Self::BadMagic => write!(f, "not a TPTA stream (bad magic)"),
            Self::UnsupportedVersion(v) => write!(f, "unsupported stream version {v}"),
            Self::ReservedFlags(v) => write!(f, "reserved header flags set: {v:#06x}"),
            Self::UnknownDtype(c) => write!(f, "unknown column type code {c}"),
            Self::BadName(n) => write!(f, "invalid column name {n:?}"),
            Self::TooLarge { what, value, limit } => {
                write!(f, "{what} {value} exceeds limit {limit}")
            }
            Self::ChecksumMismatch(what) => write!(f, "checksum mismatch in {what}"),
            Self::BadTag(t) => write!(f, "unknown message tag {t}"),
            Self::TrailingData => write!(f, "data follows end-of-stream marker"),
        }
    }
}

impl std::error::Error for ArrowIpcError {}

/// Streaming reader over a byte slice. Yields one [`RecordBatch`] at a time.
#[derive(Debug)]
pub struct ArrowReader<'a> {
    data: &'a [u8],
    pos: usize,
    schema: Vec<Field>,
    done: bool,
}

impl<'a> ArrowReader<'a> {
    /// Parse the header and schema. Batches are read lazily by [`Iterator::next`].
    pub fn new(data: &'a [u8]) -> Result<Self, ArrowIpcError> {
        let mut cur = Cursor { data, pos: 0 };
        let magic = cur.take(4)?;
        if magic != MAGIC {
            return Err(ArrowIpcError::BadMagic);
        }
        let version = cur.u16()?;
        if version != VERSION {
            return Err(ArrowIpcError::UnsupportedVersion(version));
        }
        let flags = cur.u16()?;
        if flags != 0 {
            return Err(ArrowIpcError::ReservedFlags(flags));
        }
        let schema_start = cur.pos;
        let ncols = usize::from(cur.u16()?);
        if ncols > MAX_COLUMNS {
            return Err(ArrowIpcError::TooLarge {
                what: "column count",
                value: ncols,
                limit: MAX_COLUMNS,
            });
        }
        let mut schema: Vec<Field> = Vec::with_capacity(ncols);
        for _ in 0..ncols {
            let name_len = usize::from(cur.u16()?);
            if name_len == 0 || name_len > MAX_NAME_BYTES {
                return Err(ArrowIpcError::BadName(format!("length {name_len}")));
            }
            let raw = cur.take(name_len)?;
            let name = core::str::from_utf8(raw)
                .map_err(|_| ArrowIpcError::BadName("not UTF-8".to_owned()))?
                .to_owned();
            if schema.iter().any(|f| f.name == name) {
                return Err(ArrowIpcError::BadName(format!("duplicate {name}")));
            }
            let code = cur.u8()?;
            let dtype = DataType::from_code(code).ok_or(ArrowIpcError::UnknownDtype(code))?;
            let nullable = cur.u8()? != 0;
            schema.push(Field { name, dtype, nullable });
        }
        let schema_end = cur.pos;
        let stored = cur.u32()?;
        if crc32::checksum(&data[schema_start..schema_end]) != stored {
            return Err(ArrowIpcError::ChecksumMismatch("schema"));
        }
        Ok(Self { data, pos: cur.pos, schema, done: false })
    }

    /// The schema, in column order.
    pub fn schema(&self) -> &[Field] {
        &self.schema
    }

    fn read_batch(&mut self) -> Result<RecordBatch, ArrowIpcError> {
        let body_start = self.pos;
        let mut cur = Cursor { data: self.data, pos: self.pos };
        let nrows = cur.u32()? as usize;
        if nrows > MAX_BATCH_ROWS {
            return Err(ArrowIpcError::TooLarge {
                what: "batch rows",
                value: nrows,
                limit: MAX_BATCH_ROWS,
            });
        }
        let bitmap_len = nrows.div_ceil(8);
        let mut columns = Vec::with_capacity(self.schema.len());
        for field in &self.schema {
            let validity = if field.nullable {
                let bits = unpack_bits(cur.take(bitmap_len)?, nrows);
                Some(bits)
            } else {
                None
            };
            let values = match field.dtype {
                DataType::Int64 => {
                    let raw = cur.take(nrows * 8)?;
                    Values::Int64(
                        raw.as_chunks::<8>().0.iter().map(|c| i64::from_le_bytes(*c)).collect(),
                    )
                }
                DataType::Float64 => {
                    let raw = cur.take(nrows * 8)?;
                    Values::Float64(
                        raw.as_chunks::<8>().0.iter().map(|c| f64::from_le_bytes(*c)).collect(),
                    )
                }
                DataType::Bool => Values::Bool(unpack_bits(cur.take(bitmap_len)?, nrows)),
            };
            columns.push(ColumnArray { values, validity });
        }
        let body_end = cur.pos;
        let stored = cur.u32()?;
        if crc32::checksum(&self.data[body_start..body_end]) != stored {
            return Err(ArrowIpcError::ChecksumMismatch("record batch"));
        }
        self.pos = cur.pos;
        Ok(RecordBatch { nrows, columns })
    }
}

impl Iterator for ArrowReader<'_> {
    type Item = Result<RecordBatch, ArrowIpcError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let result = (|| {
            let mut cur = Cursor { data: self.data, pos: self.pos };
            let tag = cur.u8()?;
            match tag {
                TAG_END => {
                    if cur.pos != self.data.len() {
                        return Err(ArrowIpcError::TrailingData);
                    }
                    Ok(None)
                }
                TAG_BATCH => {
                    self.pos = cur.pos;
                    self.read_batch().map(Some)
                }
                other => Err(ArrowIpcError::BadTag(other)),
            }
        })();
        match result {
            Ok(Some(batch)) => Some(Ok(batch)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// Read every batch in a stream into memory.
pub fn read_all(data: &[u8]) -> Result<Vec<RecordBatch>, ArrowIpcError> {
    ArrowReader::new(data)?.collect()
}

fn unpack_bits(bytes: &[u8], nbits: usize) -> Vec<bool> {
    (0..nbits).map(|i| bytes[i / 8] >> (i % 8) & 1 == 1).collect()
}

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], ArrowIpcError> {
        let end = self.pos.checked_add(n).ok_or(ArrowIpcError::Truncated)?;
        if end > self.data.len() {
            return Err(ArrowIpcError::Truncated);
        }
        let out = &self.data[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn u8(&mut self) -> Result<u8, ArrowIpcError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, ArrowIpcError> {
        let b = self.take(2)?;
        Ok(u16::from_le_bytes([b[0], b[1]]))
    }

    fn u32(&mut self) -> Result<u32, ArrowIpcError> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }
}

/// Helpers for writing streams, used by tests and by callers that need to
/// produce this format. Not an Arrow writer.
pub mod test_support {
    use super::{DataType, Values, MAGIC, MAX_NAME_BYTES, TAG_BATCH, TAG_END, VERSION};
    use crate::crc32;

    /// One column of a batch: its values and optional validity bits.
    pub type EncodedColumn = (Values, Option<Vec<bool>>);

    /// Encode a schema followed by batches. Each batch is a list of
    /// `(values, validity)` pairs in schema order.
    pub fn encode(
        schema: &[(&str, DataType, bool)],
        batches: &[(usize, Vec<EncodedColumn>)],
    ) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&VERSION.to_le_bytes());
        out.extend_from_slice(&0u16.to_le_bytes());
        let mut schema_bytes = Vec::new();
        schema_bytes.extend_from_slice(&(schema.len() as u16).to_le_bytes());
        for (name, dtype, nullable) in schema {
            assert!(name.len() <= MAX_NAME_BYTES);
            schema_bytes.extend_from_slice(&(name.len() as u16).to_le_bytes());
            schema_bytes.extend_from_slice(name.as_bytes());
            schema_bytes.push(dtype.code());
            schema_bytes.push(u8::from(*nullable));
        }
        out.extend_from_slice(&schema_bytes);
        out.extend_from_slice(&crc32::checksum(&schema_bytes).to_le_bytes());
        for (nrows, cols) in batches {
            out.push(TAG_BATCH);
            let body_start = out.len();
            out.extend_from_slice(&(*nrows as u32).to_le_bytes());
            for ((values, validity), (_, dtype, nullable)) in cols.iter().zip(schema) {
                if *nullable {
                    let bits = validity.clone().unwrap_or_else(|| vec![true; *nrows]);
                    out.extend_from_slice(&pack_bits(&bits));
                }
                match (values, dtype) {
                    (Values::Int64(v), DataType::Int64) => {
                        for x in v {
                            out.extend_from_slice(&x.to_le_bytes());
                        }
                    }
                    (Values::Float64(v), DataType::Float64) => {
                        for x in v {
                            out.extend_from_slice(&x.to_le_bytes());
                        }
                    }
                    (Values::Bool(v), DataType::Bool) => out.extend_from_slice(&pack_bits(v)),
                    _ => panic!("value type does not match schema"),
                }
            }
            let crc = crc32::checksum(&out[body_start..]);
            out.extend_from_slice(&crc.to_le_bytes());
        }
        out.push(TAG_END);
        out
    }

    /// Bit-pack booleans, least significant bit first.
    pub fn pack_bits(bits: &[bool]) -> Vec<u8> {
        let mut out = vec![0u8; bits.len().div_ceil(8)];
        for (i, &b) in bits.iter().enumerate() {
            if b {
                out[i / 8] |= 1 << (i % 8);
            }
        }
        out
    }

    /// An empty stream with an empty schema.
    pub fn empty_stream() -> Vec<u8> {
        encode(&[], &[])
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::{empty_stream, encode};
    use super::*;

    fn sample_bytes() -> Vec<u8> {
        encode(
            &[
                ("run", DataType::Int64, false),
                ("pt", DataType::Float64, true),
                ("trig", DataType::Bool, false),
            ],
            &[
                (
                    3,
                    vec![
                        (Values::Int64(vec![1, 1, 2]), None),
                        (Values::Float64(vec![10.5, 0.0, 3.25]), Some(vec![true, false, true])),
                        (Values::Bool(vec![true, false, true]), None),
                    ],
                ),
                (
                    1,
                    vec![
                        (Values::Int64(vec![-7]), None),
                        (Values::Float64(vec![f64::MIN_POSITIVE]), None),
                        (Values::Bool(vec![false]), None),
                    ],
                ),
            ],
        )
    }

    #[test]
    fn round_trip_two_batches() {
        let bytes = sample_bytes();
        let batches = read_all(&bytes).unwrap();
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].nrows(), 3);
        assert_eq!(batches[0].columns()[0].values, Values::Int64(vec![1, 1, 2]));
        assert_eq!(batches[0].columns()[1].null_count(), 1);
        assert_eq!(batches[0].columns()[1].validity, Some(vec![true, false, true]));
        assert_eq!(batches[0].columns()[2].values, Values::Bool(vec![true, false, true]));
        assert_eq!(batches[1].columns()[0].values, Values::Int64(vec![-7]));
        assert_eq!(batches[1].columns()[1].values, Values::Float64(vec![f64::MIN_POSITIVE]));
    }

    #[test]
    fn schema_is_exposed() {
        let bytes = sample_bytes();
        let r = ArrowReader::new(&bytes).unwrap();
        let names: Vec<&str> = r.schema().iter().map(|f| f.name.as_str()).collect();
        assert_eq!(names, ["run", "pt", "trig"]);
        assert!(r.schema()[1].nullable);
    }

    #[test]
    fn empty_stream_has_no_batches() {
        assert_eq!(read_all(&empty_stream()).unwrap().len(), 0);
    }

    #[test]
    fn zero_row_batch() {
        let bytes = encode(
            &[("x", DataType::Float64, true)],
            &[(0, vec![(Values::Float64(vec![]), Some(vec![]))])],
        );
        let b = read_all(&bytes).unwrap();
        assert_eq!(b[0].nrows(), 0);
        assert_eq!(b[0].columns()[0].values.len(), 0);
    }

    #[test]
    fn bad_magic_and_version() {
        let mut bytes = sample_bytes();
        bytes[0] = b'X';
        assert_eq!(read_all(&bytes), Err(ArrowIpcError::BadMagic));
        let mut bytes = sample_bytes();
        bytes[4] = 9;
        assert_eq!(read_all(&bytes), Err(ArrowIpcError::UnsupportedVersion(9)));
    }

    #[test]
    fn every_truncation_is_an_error() {
        let bytes = sample_bytes();
        for n in 0..bytes.len() {
            assert!(read_all(&bytes[..n]).is_err(), "prefix of {n} bytes accepted");
        }
    }

    #[test]
    fn corrupted_payload_fails_checksum() {
        let mut bytes = sample_bytes();
        // Flip a byte inside the first batch's int64 data.
        let idx = bytes.len() / 2;
        bytes[idx] ^= 0x40;
        assert!(read_all(&bytes).is_err());
    }

    #[test]
    fn trailing_data_rejected() {
        let mut bytes = sample_bytes();
        bytes.push(0);
        assert_eq!(read_all(&bytes), Err(ArrowIpcError::TrailingData));
    }

    #[test]
    fn huge_row_count_rejected_without_allocating() {
        let mut bytes = encode(&[("x", DataType::Int64, false)], &[]);
        // Replace the end marker with a batch claiming u32::MAX rows.
        bytes.pop();
        bytes.push(TAG_BATCH);
        bytes.extend_from_slice(&u32::MAX.to_le_bytes());
        let err = read_all(&bytes).unwrap_err();
        assert!(matches!(err, ArrowIpcError::TooLarge { what: "batch rows", .. }));
    }

    #[test]
    fn oversized_name_and_duplicate_names_rejected() {
        // Duplicate: build a two-column schema with the same name manually.
        let dup = encode(&[("a", DataType::Int64, false), ("a", DataType::Int64, false)], &[]);
        assert!(matches!(read_all(&dup), Err(ArrowIpcError::BadName(_))));
        // A zero-length name is rejected before the bytes are read.
        let mut empty_name = encode(&[("a", DataType::Int64, false)], &[]);
        empty_name[10] = 0;
        assert!(matches!(read_all(&empty_name), Err(ArrowIpcError::BadName(_))));
    }

    #[test]
    fn unknown_dtype_rejected() {
        let mut bytes = encode(&[("a", DataType::Int64, false)], &[]);
        // Byte layout: magic(4) version(2) flags(2) ncols(2) namelen(2) name(1) dtype(1)
        bytes[4 + 2 + 2 + 2 + 2 + 1] = 42;
        assert_eq!(read_all(&bytes), Err(ArrowIpcError::UnknownDtype(42)));
    }

    #[test]
    fn bit_packing_round_trip() {
        let bits: Vec<bool> = (0..19).map(|i| i % 3 == 0).collect();
        let packed = test_support::pack_bits(&bits);
        assert_eq!(packed.len(), 3);
        assert_eq!(unpack_bits(&packed, 19), bits);
    }
}
