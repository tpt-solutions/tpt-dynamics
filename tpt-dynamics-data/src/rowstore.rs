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

use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use crate::lazy::{Column, DType, DataFrame, Field, LazyError, Schema, Source};

const MAGIC: &[u8; 8] = b"TPTROWS1";

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
#[derive(Clone, Debug)]
pub struct RowStore {
    path: PathBuf,
    schema: Schema,
    nrows: usize,
    /// Absolute byte offset of each column's data.
    offsets: Vec<u64>,
}

impl RowStore {
    /// Opens a row store and checks its header against the file length.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Invalid`] if the file is missing, truncated, or not
    /// in row-store format.
    pub fn open(path: impl Into<PathBuf>) -> Result<Self, LazyError> {
        let path = path.into();
        let meta = fs::metadata(&path).map_err(|e| io_err(&path, "stat", &e))?;
        let mut file = File::open(&path).map_err(|e| io_err(&path, "open", &e))?;
        let mut head = [0u8; 20];
        file.read_exact(&mut head).map_err(|e| io_err(&path, "read header", &e))?;
        if &head[..8] != MAGIC {
            return Err(LazyError::Invalid(format!("'{}' is not a row store", path.display())));
        }
        let nrows = u64::from_le_bytes(head[8..16].try_into().expect("8 bytes"));
        let ncols = u32::from_le_bytes(head[16..20].try_into().expect("4 bytes")) as usize;
        let mut pos = 20u64;
        let mut schema = Vec::new();
        let mut dtypes = Vec::new();
        for _ in 0..ncols {
            let mut fixed = [0u8; 5];
            file.read_exact(&mut fixed).map_err(|e| io_err(&path, "read column", &e))?;
            let dtype = dtype_from_code(fixed[0])?;
            let name_len = u32::from_le_bytes(fixed[1..5].try_into().expect("4 bytes")) as usize;
            if name_len > meta.len() as usize {
                return Err(LazyError::Invalid("row store: column name too long".into()));
            }
            let mut name = vec![0u8; name_len];
            file.read_exact(&mut name).map_err(|e| io_err(&path, "read name", &e))?;
            let name = String::from_utf8(name)
                .map_err(|_| LazyError::Invalid("row store: column name is not UTF-8".into()))?;
            pos += 5 + name_len as u64;
            schema.push(Field { name, dtype });
            dtypes.push(dtype);
        }
        let mut offsets = Vec::with_capacity(ncols);
        let mut at = pos;
        for dtype in &dtypes {
            offsets.push(at);
            at = at
                .checked_add(
                    nrows
                        .checked_mul(width(*dtype))
                        .ok_or_else(|| LazyError::Invalid("row store: size overflow".into()))?,
                )
                .ok_or_else(|| LazyError::Invalid("row store: size overflow".into()))?;
        }
        if at != meta.len() {
            return Err(LazyError::Invalid(format!(
                "row store '{}' has {} bytes, header says {at}",
                path.display(),
                meta.len()
            )));
        }
        Ok(Self { path, schema, nrows: nrows as usize, offsets })
    }

    /// Path of the underlying file.
    pub fn path(&self) -> &Path {
        &self.path
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

    fn read(
        &self,
        rows: std::ops::Range<usize>,
        projection: &[String],
    ) -> Result<DataFrame, LazyError> {
        let end = rows.end.min(self.nrows);
        let start = rows.start.min(end);
        let n = end - start;
        let mut file = File::open(&self.path).map_err(|e| io_err(&self.path, "open", &e))?;
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
            file.seek(SeekFrom::Start(self.offsets[idx] + start as u64 * w))
                .map_err(|e| io_err(&self.path, "seek", &e))?;
            file.read_exact(&mut buf).map_err(|e| io_err(&self.path, "read", &e))?;
            columns.push((name.clone(), Self::decode(dtype, &buf)));
        }
        DataFrame::with_rows(n, columns)
    }
}
