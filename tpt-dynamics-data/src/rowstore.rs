//! Fixed-width row store: a columnar file that can serve any row range with one
//! seek per column.
//!
//! The lazy engine reads sources by row range, so a distributed worker needs a
//! file where a range costs one seek and one read per column. The columnar
//! format in [`crate::columnar`] stores whole tables and cannot do that.
//!
//! Layout (little endian):
//!
//! ```text
//! magic      8 bytes  b"TPTROWS1"
//! nrows      u64
//! ncols      u32
//! per column: dtype u8 (0 = i64, 1 = f64, 2 = bool), name_len u32, name bytes
//! data       column 0 values, then column 1 values, ...
//!            (8 bytes per i64 or f64 value, 1 byte per bool value)
//! ```
//!
//! Each column's data is contiguous, so the value of column `c` at row `r` sits
//! at `offset(c) + r * width(c)`. Reading rows `a..b` of one column is a single
//! contiguous read.
//!
//! A row store can be read from a local file ([`RowStore::open`]), from any
//! [`RangeRead`] implementation, or over plain HTTP ([`RowStore::open_http`]).

use std::fmt;
use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::http::HttpRange;
use crate::lazy::{Column, DType, DataFrame, Field, LazyError, Schema, Source};

const MAGIC: &[u8; 8] = b"TPTROWS1";
/// Bytes in the fixed part of the header: magic, row count, column count.
const FIXED_HEADER: u64 = 20;

/// Random access to a byte sequence of known length.
///
/// Implementations are called from several threads at once, so reads must not
/// depend on a shared cursor.
pub trait RangeRead: Send + Sync {
    /// Total length in bytes.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Invalid`] if the length cannot be determined.
    fn byte_len(&self) -> Result<u64, LazyError>;

    /// Fills `buf` with the bytes starting at `offset`.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Invalid`] if the read fails or runs past the end.
    fn read_at(&self, offset: u64, buf: &mut [u8]) -> Result<(), LazyError>;
}

/// A local file. Each read opens the file, so the value is cheap to clone and
/// share.
#[derive(Clone, Debug)]
pub struct FileRange {
    path: PathBuf,
}

impl FileRange {
    /// A range over the file at `path`. The file is not opened until a read.
    pub fn new(path: impl Into<PathBuf>) -> Self {
        Self { path: path.into() }
    }
}

impl RangeRead for FileRange {
    fn byte_len(&self) -> Result<u64, LazyError> {
        fs::metadata(&self.path).map(|m| m.len()).map_err(|e| io_err(&self.path, "stat", &e))
    }

    fn read_at(&self, offset: u64, buf: &mut [u8]) -> Result<(), LazyError> {
        let mut file = File::open(&self.path).map_err(|e| io_err(&self.path, "open", &e))?;
        file.seek(SeekFrom::Start(offset)).map_err(|e| io_err(&self.path, "seek", &e))?;
        file.read_exact(buf).map_err(|e| io_err(&self.path, "read", &e))
    }
}

fn width(dtype: DType) -> u64 {
    match dtype {
        DType::I64 | DType::F64 => 8,
        DType::Bool => 1,
    }
}

fn dtype_code(dtype: DType) -> u8 {
    match dtype {
        DType::I64 => 0,
        DType::F64 => 1,
        DType::Bool => 2,
    }
}

fn dtype_from_code(code: u8) -> Result<DType, LazyError> {
    match code {
        0 => Ok(DType::I64),
        1 => Ok(DType::F64),
        2 => Ok(DType::Bool),
        other => Err(LazyError::Invalid(format!("row store: unknown dtype code {other}"))),
    }
}

fn io_err(path: &Path, what: &str, e: &std::io::Error) -> LazyError {
    LazyError::Invalid(format!("row store '{}': {what}: {e}", path.display()))
}

/// Reads `len` bytes at `offset`, refusing ranges that run past `total`.
fn read_vec(
    source: &dyn RangeRead,
    offset: u64,
    len: u64,
    total: u64,
) -> Result<Vec<u8>, LazyError> {
    let end = offset.checked_add(len).filter(|&e| e <= total);
    if end.is_none() {
        return Err(LazyError::Invalid("row store: header runs past the end of the data".into()));
    }
    let len = usize::try_from(len)
        .map_err(|_| LazyError::Invalid("row store: header too large".into()))?;
    let mut buf = vec![0u8; len];
    source.read_at(offset, &mut buf)?;
    Ok(buf)
}

/// Writes `df` to `path` in row-store format, replacing any existing file.
///
/// # Errors
///
/// Returns [`LazyError::Invalid`] if the file cannot be written.
pub fn write_rowstore(path: &Path, df: &DataFrame) -> Result<(), LazyError> {
    let mut out = Vec::new();
    out.extend_from_slice(MAGIC);
    out.extend_from_slice(&(df.len() as u64).to_le_bytes());
    let schema = df.schema();
    out.extend_from_slice(&(schema.len() as u32).to_le_bytes());
    for field in &schema {
        out.push(dtype_code(field.dtype));
        out.extend_from_slice(&(field.name.len() as u32).to_le_bytes());
        out.extend_from_slice(field.name.as_bytes());
    }
    for field in &schema {
        let column =
            df.column(&field.name).ok_or_else(|| LazyError::UnknownColumn(field.name.clone()))?;
        match column {
            Column::I64(v) => {
                for x in v {
                    out.extend_from_slice(&x.to_le_bytes());
                }
            }
            Column::F64(v) => {
                for x in v {
                    out.extend_from_slice(&x.to_bits().to_le_bytes());
                }
            }
            Column::Bool(v) => out.extend(v.iter().map(|&x| u8::from(x))),
        }
    }
    let mut file = File::create(path).map_err(|e| io_err(path, "create", &e))?;
    file.write_all(&out).map_err(|e| io_err(path, "write", &e))?;
    file.sync_all().map_err(|e| io_err(path, "sync", &e))
}

/// A row store opened for reading. Implements [`Source`], so it can back a
/// [`crate::lazy::LazyFrame`] or be served to distributed workers.
#[derive(Clone)]
pub struct RowStore {
    source: Arc<dyn RangeRead>,
    label: String,
    schema: Schema,
    nrows: usize,
    /// Absolute byte offset of each column's data.
    offsets: Vec<u64>,
}

impl fmt::Debug for RowStore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RowStore")
            .field("label", &self.label)
            .field("rows", &self.nrows)
            .field("columns", &self.schema.len())
            .finish_non_exhaustive()
    }
}

impl RowStore {
    /// Opens a row store file and checks its header against the file length.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Invalid`] if the file is missing, truncated, or not
    /// in row-store format.
    pub fn open(path: impl Into<PathBuf>) -> Result<Self, LazyError> {
        let range = FileRange::new(path);
        let label = range.path.display().to_string();
        Self::open_range(Arc::new(range), label)
    }

    /// Opens a row store reachable through `source`. `label` names it in error
    /// messages and in `Debug` output.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Invalid`] if the header cannot be read or does not
    /// match the length of `source`.
    pub fn open_range(
        source: Arc<dyn RangeRead>,
        label: impl Into<String>,
    ) -> Result<Self, LazyError> {
        let label = label.into();
        let total = source.byte_len()?;
        let head = read_vec(source.as_ref(), 0, FIXED_HEADER.min(total), total)?;
        if head.len() < FIXED_HEADER as usize || &head[..8] != MAGIC {
            return Err(LazyError::Invalid(format!("'{label}' is not a row store")));
        }
        let nrows = u64::from_le_bytes(head[8..16].try_into().expect("8 bytes"));
        let ncols = u32::from_le_bytes(head[16..20].try_into().expect("4 bytes")) as usize;
        let mut pos = FIXED_HEADER;
        let mut schema = Vec::new();
        let mut dtypes = Vec::new();
        for _ in 0..ncols {
            let fixed = read_vec(source.as_ref(), pos, 5, total)?;
            let dtype = dtype_from_code(fixed[0])?;
            let name_len = u64::from(u32::from_le_bytes(fixed[1..5].try_into().expect("4 bytes")));
            if name_len > total {
                return Err(LazyError::Invalid("row store: column name too long".into()));
            }
            let name = read_vec(source.as_ref(), pos + 5, name_len, total)?;
            let name = String::from_utf8(name)
                .map_err(|_| LazyError::Invalid("row store: column name is not UTF-8".into()))?;
            pos += 5 + name_len;
            schema.push(Field { name, dtype });
            dtypes.push(dtype);
        }
        let mut offsets = Vec::with_capacity(ncols);
        let mut at = pos;
        for dtype in &dtypes {
            offsets.push(at);
            let size = nrows
                .checked_mul(width(*dtype))
                .ok_or_else(|| LazyError::Invalid("row store: size overflow".into()))?;
            at = at
                .checked_add(size)
                .ok_or_else(|| LazyError::Invalid("row store: size overflow".into()))?;
        }
        if at != total {
            return Err(LazyError::Invalid(format!(
                "row store '{label}' has {total} bytes, header says {at}"
            )));
        }
        let nrows = usize::try_from(nrows)
            .map_err(|_| LazyError::Invalid("row store: row count too large".into()))?;
        Ok(Self { source, label, schema, nrows, offsets })
    }

    /// Opens a row store served over plain HTTP (`http://host[:port]/path`).
    /// The server must answer `Range` requests with `206 Partial Content`.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Invalid`] if the URL is not plain HTTP, the server
    /// does not support range requests, or the header is invalid.
    pub fn open_http(url: &str) -> Result<Self, LazyError> {
        let range = HttpRange::new(url)?;
        Self::open_range(Arc::new(range), url)
    }

    /// Name of the store: the file path or URL it was opened from.
    pub fn label(&self) -> &str {
        &self.label
    }

    fn decode(dtype: DType, bytes: &[u8]) -> Column {
        match dtype {
            DType::I64 => Column::I64(
                bytes.as_chunks::<8>().0.iter().map(|c| i64::from_le_bytes(*c)).collect(),
            ),
            DType::F64 => Column::F64(
                bytes
                    .as_chunks::<8>()
                    .0
                    .iter()
                    .map(|c| f64::from_bits(u64::from_le_bytes(*c)))
                    .collect(),
            ),
            DType::Bool => Column::Bool(bytes.iter().map(|&b| b != 0).collect()),
        }
    }
}

impl Source for RowStore {
    fn schema(&self) -> Schema {
        self.schema.clone()
    }

    fn num_rows(&self) -> usize {
        self.nrows
    }

    fn read(&self, rows: Range<usize>, projection: &[String]) -> Result<DataFrame, LazyError> {
        let end = rows.end.min(self.nrows);
        let start = rows.start.min(end);
        let n = end - start;
        let mut columns = Vec::with_capacity(projection.len());
        for name in projection {
            let idx = self
                .schema
                .iter()
                .position(|f| &f.name == name)
                .ok_or_else(|| LazyError::UnknownColumn(name.clone()))?;
            let dtype = self.schema[idx].dtype;
            let w = width(dtype);
            let mut buf = vec![0u8; n * w as usize];
            self.source.read_at(self.offsets[idx] + start as u64 * w, &mut buf)?;
            columns.push((name.clone(), Self::decode(dtype, &buf)));
        }
        DataFrame::with_rows(n, columns)
    }
}
