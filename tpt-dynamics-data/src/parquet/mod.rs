//! Read-only Parquet import (ADR-008).
//!
//! This is a clean-room reader written from the public Parquet format
//! specification. It uses no Arrow or Parquet crates, which are Apache-2.0 only
//! and rejected by ADR-005.
//!
//! Supported:
//!
//! * flat schemas (a root with leaf children only);
//! * physical types `BOOLEAN`, `INT32`, `INT64`, `FLOAT`, `DOUBLE`, read as
//!   [`Column::Bool`], [`Column::I64`] and [`Column::F64`]. Logical annotations
//!   are not interpreted, so unsigned and timestamp columns come back as their
//!   raw signed integers;
//! * encodings `PLAIN`, `RLE_DICTIONARY` and `PLAIN_DICTIONARY`;
//! * codecs `UNCOMPRESSED`, `SNAPPY` and `GZIP`;
//! * data page v1 and dictionary pages;
//! * `REQUIRED` and `OPTIONAL` columns, provided an optional column holds no nulls.
//!
//! Not supported (returns [`ParquetError::Unsupported`]): nested and repeated
//! fields, `BYTE_ARRAY`, `INT96` and fixed-length binary columns, data page v2,
//! the `ZSTD`/`LZ4`/`BROTLI` codecs, and any null value (the column types have no
//! null representation).
//!
//! Memory: the footer is read once, and each column chunk is read and decoded
//! on its own. Peak memory is one chunk plus the decoded output.
//!
//! # Examples
//!
//! ```no_run
//! use tpt_dynamics_data::lazy::{col, LazyFrame};
//! use tpt_dynamics_data::parquet::read_frame;
//!
//! let df = read_frame("events.parquet", None).unwrap();
//! let out = LazyFrame::scan(df).filter(col("pt").gt(20.0)).collect().unwrap();
//! println!("{} rows pass", out.len());
//! ```

mod rle;
mod thrift;

use std::fmt;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use crate::lazy::{Column, DataFrame, LazyError};
use rle::{decode_bools, decode_hybrid};
use thrift::Thrift;

/// Errors from reading a Parquet file.
#[derive(Debug)]
pub enum ParquetError {
    /// Reading the file failed.
    Io(std::io::Error),
    /// The file is not valid Parquet, or its metadata is inconsistent.
    Format(String),
    /// The file uses a feature this reader does not implement.
    Unsupported(String),
    /// No column with the requested name exists.
    UnknownColumn(String),
}

impl fmt::Display for ParquetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Format(m) => write!(f, "invalid Parquet file: {m}"),
            Self::Unsupported(m) => write!(f, "unsupported Parquet feature: {m}"),
            Self::UnknownColumn(c) => write!(f, "unknown column `{c}`"),
        }
    }
}

impl std::error::Error for ParquetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for ParquetError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<LazyError> for ParquetError {
    fn from(e: LazyError) -> Self {
        Self::Format(e.to_string())
    }
}

const MAGIC: &[u8; 4] = b"PAR1";

/// Upper bound on one page's uncompressed size (1 GiB), so corrupt headers
/// cannot request absurd allocations.
const MAX_PAGE_BYTES: usize = 1 << 30;

// Physical types.
const PT_BOOLEAN: i32 = 0;
const PT_INT32: i32 = 1;
const PT_INT64: i32 = 2;
const PT_FLOAT: i32 = 4;
const PT_DOUBLE: i32 = 5;

// Repetition.
const REP_REQUIRED: i32 = 0;
const REP_OPTIONAL: i32 = 1;

// Encodings.
const ENC_PLAIN: i32 = 0;
const ENC_PLAIN_DICTIONARY: i32 = 2;
const ENC_RLE_DICTIONARY: i32 = 8;

// Codecs.
const CODEC_UNCOMPRESSED: i32 = 0;
const CODEC_SNAPPY: i32 = 1;
const CODEC_GZIP: i32 = 2;

// Page types.
const PAGE_DATA: i32 = 0;
const PAGE_DICTIONARY: i32 = 2;
const PAGE_DATA_V2: i32 = 3;

#[derive(Debug)]
struct SchemaElem {
    name: String,
    ptype: Option<i32>,
    repetition: Option<i32>,
    num_children: Option<i32>,
}

#[derive(Debug)]
struct ColumnMeta {
    ptype: i32,
    codec: i32,
    num_values: i64,
    total_compressed_size: i64,
    data_page_offset: i64,
    dictionary_page_offset: Option<i64>,
}

#[derive(Debug)]
struct RowGroup {
    num_rows: i64,
    columns: Vec<ColumnMeta>,
}

#[derive(Debug)]
struct Footer {
    num_rows: i64,
    schema: Vec<SchemaElem>,
    row_groups: Vec<RowGroup>,
}

/// A leaf (top-level) column in the schema.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LeafInfo {
    /// Column name.
    pub name: String,
    /// Physical type code (`BOOLEAN` = 0, `INT32` = 1, `INT64` = 2, `FLOAT` = 4,
    /// `DOUBLE` = 5, `BYTE_ARRAY` = 6, and so on).
    pub physical_type: i32,
    /// Whether the column is `OPTIONAL`.
    pub optional: bool,
}

fn schema_leaves(schema: &[SchemaElem]) -> Result<Vec<LeafInfo>, ParquetError> {
    let root = schema.first().ok_or_else(|| ParquetError::Format("empty schema".into()))?;
    let n = root.num_children.unwrap_or(0).max(0) as usize;
    if schema.len() != n + 1 {
        return Err(ParquetError::Unsupported("nested schema (group fields)".into()));
    }
    let mut out = Vec::with_capacity(n);
    for el in &schema[1..] {
        let physical_type = el
            .ptype
            .ok_or_else(|| ParquetError::Unsupported("nested schema (group fields)".into()))?;
        let optional = match el.repetition {
            Some(REP_OPTIONAL) => true,
            Some(REP_REQUIRED) | None => false,
            Some(_) => {
                return Err(ParquetError::Unsupported(format!("repeated field `{}`", el.name)))
            }
        };
        out.push(LeafInfo { name: el.name.clone(), physical_type, optional });
    }
    Ok(out)
}

fn read_list<T>(
    th: &mut Thrift<'_>,
    mut f: impl FnMut(&mut Thrift<'_>) -> Result<T, ParquetError>,
) -> Result<Vec<T>, ParquetError> {
    let (n, _) = th.list()?;
    (0..n).map(|_| f(th)).collect()
}

fn parse_schema_elem(th: &mut Thrift<'_>) -> Result<SchemaElem, ParquetError> {
    let mut el =
        SchemaElem { name: String::new(), ptype: None, repetition: None, num_children: None };
    th.structure(|th, id, _| match id {
        1 => {
            el.ptype = Some(th.zigzag()? as i32);
            Ok(true)
        }
        3 => {
            el.repetition = Some(th.zigzag()? as i32);
            Ok(true)
        }
        4 => {
            el.name = th.string()?;
            Ok(true)
        }
        5 => {
            el.num_children = Some(th.zigzag()? as i32);
            Ok(true)
        }
        _ => Ok(false),
    })?;
    Ok(el)
}

fn parse_column_meta(th: &mut Thrift<'_>) -> Result<ColumnMeta, ParquetError> {
    let mut m = ColumnMeta {
        ptype: -1,
        codec: CODEC_UNCOMPRESSED,
        num_values: 0,
        total_compressed_size: 0,
        data_page_offset: 0,
        dictionary_page_offset: None,
    };
    th.structure(|th, id, _| match id {
        1 => {
            m.ptype = th.zigzag()? as i32;
            Ok(true)
        }
        4 => {
            m.codec = th.zigzag()? as i32;
            Ok(true)
        }
        5 => {
            m.num_values = th.zigzag()?;
            Ok(true)
        }
        7 => {
            m.total_compressed_size = th.zigzag()?;
            Ok(true)
        }
        9 => {
            m.data_page_offset = th.zigzag()?;
            Ok(true)
        }
        11 => {
            m.dictionary_page_offset = Some(th.zigzag()?);
            Ok(true)
        }
        _ => Ok(false),
    })?;
    if m.ptype < 0 {
        return Err(ParquetError::Format("column metadata has no type".into()));
    }
    Ok(m)
}

fn parse_footer(bytes: &[u8]) -> Result<Footer, ParquetError> {
    let mut f = Footer { num_rows: 0, schema: Vec::new(), row_groups: Vec::new() };
    let mut th = Thrift::new(bytes);
    th.structure(|th, id, _| match id {
        2 => {
            f.schema = read_list(th, parse_schema_elem)?;
            Ok(true)
        }
        3 => {
            f.num_rows = th.zigzag()?;
            Ok(true)
        }
        4 => {
            f.row_groups = read_list(th, |th| {
                let mut rg = RowGroup { num_rows: 0, columns: Vec::new() };
                th.structure(|th, id, _| match id {
                    1 => {
                        rg.columns = read_list(th, |th| {
                            let mut meta = None;
                            th.structure(|th, id, _| {
                                if id == 3 {
                                    meta = Some(parse_column_meta(th)?);
                                    Ok(true)
                                } else {
                                    Ok(false)
                                }
                            })?;
                            meta.ok_or_else(|| {
                                ParquetError::Format("column chunk has no metadata".into())
                            })
                        })?;
                        Ok(true)
                    }
                    3 => {
                        rg.num_rows = th.zigzag()?;
                        Ok(true)
                    }
                    _ => Ok(false),
                })?;
                Ok(rg)
            })?;
            Ok(true)
        }
        _ => Ok(false),
    })?;
    Ok(f)
}

/// Page header fields needed to decode a page.
#[derive(Debug)]
struct PageHeader {
    kind: i32,
    uncompressed_size: usize,
    compressed_size: usize,
    num_values: usize,
    encoding: i32,
}

fn parse_page_header(th: &mut Thrift<'_>) -> Result<PageHeader, ParquetError> {
    let mut h = PageHeader {
        kind: -1,
        uncompressed_size: 0,
        compressed_size: 0,
        num_values: 0,
        encoding: -1,
    };
    th.structure(|th, id, _| match id {
        1 => {
            h.kind = th.zigzag()? as i32;
            Ok(true)
        }
        2 => {
            h.uncompressed_size = th.zigzag()? as usize;
            Ok(true)
        }
        3 => {
            h.compressed_size = th.zigzag()? as usize;
            Ok(true)
        }
        // DataPageHeader (5), DictionaryPageHeader (7): both carry num_values
        // (field 1) and encoding (field 2).
        5 | 7 => {
            th.structure(|th, id, _| match id {
                1 => {
                    h.num_values = th.zigzag()? as usize;
                    Ok(true)
                }
                2 => {
                    h.encoding = th.zigzag()? as i32;
                    Ok(true)
                }
                _ => Ok(false),
            })?;
            Ok(true)
        }
        _ => Ok(false),
    })?;
    if h.kind < 0 {
        return Err(ParquetError::Format("page header has no type".into()));
    }
    Ok(h)
}

/// Decoded values of one physical type.
#[derive(Clone, Debug, PartialEq)]
enum Values {
    I64(Vec<i64>),
    F64(Vec<f64>),
    Bool(Vec<bool>),
}

impl Values {
    fn empty(physical_type: i32) -> Result<Self, ParquetError> {
        match physical_type {
            PT_BOOLEAN => Ok(Self::Bool(Vec::new())),
            PT_INT32 | PT_INT64 => Ok(Self::I64(Vec::new())),
            PT_FLOAT | PT_DOUBLE => Ok(Self::F64(Vec::new())),
            other => Err(ParquetError::Unsupported(physical_type_name(other))),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::I64(v) => v.len(),
            Self::F64(v) => v.len(),
            Self::Bool(v) => v.len(),
        }
    }

    fn append(&mut self, other: Self) -> Result<(), ParquetError> {
        match (self, other) {
            (Self::I64(a), Self::I64(b)) => a.extend(b),
            (Self::F64(a), Self::F64(b)) => a.extend(b),
            (Self::Bool(a), Self::Bool(b)) => a.extend(b),
            _ => return Err(ParquetError::Format("mixed value types in a column".into())),
        }
        Ok(())
    }

    /// Looks up dictionary entries by index.
    fn gather(&self, indices: &[u32]) -> Result<Self, ParquetError> {
        let oob = || ParquetError::Format("dictionary index out of range".into());
        Ok(match self {
            Self::I64(d) => Self::I64(
                indices
                    .iter()
                    .map(|&i| d.get(i as usize).copied().ok_or_else(oob))
                    .collect::<Result<_, _>>()?,
            ),
            Self::F64(d) => Self::F64(
                indices
                    .iter()
                    .map(|&i| d.get(i as usize).copied().ok_or_else(oob))
                    .collect::<Result<_, _>>()?,
            ),
            Self::Bool(d) => Self::Bool(
                indices
                    .iter()
                    .map(|&i| d.get(i as usize).copied().ok_or_else(oob))
                    .collect::<Result<_, _>>()?,
            ),
        })
    }

    fn into_column(self) -> Column {
        match self {
            Self::I64(v) => Column::I64(v),
            Self::F64(v) => Column::F64(v),
            Self::Bool(v) => Column::Bool(v),
        }
    }
}

fn physical_type_name(pt: i32) -> String {
    match pt {
        3 => "INT96 columns".into(),
        6 => "BYTE_ARRAY (string/binary) columns".into(),
        7 => "FIXED_LEN_BYTE_ARRAY columns".into(),
        other => format!("physical type {other}"),
    }
}

/// Decodes `count` PLAIN values of physical type `pt`.
fn decode_plain(data: &[u8], pt: i32, count: usize) -> Result<Values, ParquetError> {
    let need = |width: usize| -> Result<&[u8], ParquetError> {
        count
            .checked_mul(width)
            .and_then(|n| data.get(..n))
            .ok_or_else(|| ParquetError::Format("truncated data page".into()))
    };
    Ok(match pt {
        PT_BOOLEAN => Values::Bool(decode_bools(data, count)?),
        PT_INT32 => Values::I64(
            need(4)?.as_chunks::<4>().0.iter().map(|c| i64::from(i32::from_le_bytes(*c))).collect(),
        ),
        PT_INT64 => Values::I64(
            need(8)?.as_chunks::<8>().0.iter().map(|c| i64::from_le_bytes(*c)).collect(),
        ),
        PT_FLOAT => Values::F64(
            need(4)?.as_chunks::<4>().0.iter().map(|c| f64::from(f32::from_le_bytes(*c))).collect(),
        ),
        PT_DOUBLE => Values::F64(
            need(8)?.as_chunks::<8>().0.iter().map(|c| f64::from_le_bytes(*c)).collect(),
        ),
        other => return Err(ParquetError::Unsupported(physical_type_name(other))),
    })
}

fn decompress(codec: i32, body: &[u8], uncompressed: usize) -> Result<Vec<u8>, ParquetError> {
    match codec {
        CODEC_UNCOMPRESSED => Ok(body.to_vec()),
        CODEC_SNAPPY => snap::raw::Decoder::new()
            .decompress_vec(body)
            .map_err(|e| ParquetError::Format(format!("snappy page: {e}"))),
        CODEC_GZIP => {
            let mut out = Vec::with_capacity(uncompressed);
            flate2::read::GzDecoder::new(body).read_to_end(&mut out)?;
            Ok(out)
        }
        other => Err(ParquetError::Unsupported(format!("compression codec {other}"))),
    }
}

/// Reads and decodes one column chunk (all of its pages).
fn read_chunk<R: Read + Seek>(
    source: &mut R,
    file_len: u64,
    leaf: &LeafInfo,
    meta: &ColumnMeta,
) -> Result<Values, ParquetError> {
    let start =
        meta.dictionary_page_offset.map_or(meta.data_page_offset, |d| d.min(meta.data_page_offset));
    let start = u64::try_from(start).map_err(|_| ParquetError::Format("negative offset".into()))?;
    let len = usize::try_from(meta.total_compressed_size)
        .map_err(|_| ParquetError::Format("negative chunk size".into()))?;
    if start.checked_add(len as u64).is_none_or(|end| end > file_len) {
        return Err(ParquetError::Format(format!(
            "column `{}` chunk lies outside the file",
            leaf.name
        )));
    }
    source.seek(SeekFrom::Start(start))?;
    let mut buf = vec![0u8; len];
    source.read_exact(&mut buf)?;

    let mut out = Values::empty(meta.ptype)?;
    let mut dict: Option<Values> = None;
    let mut pos = 0usize;
    let mut seen = 0i64;

    while seen < meta.num_values {
        let rest = buf
            .get(pos..)
            .ok_or_else(|| ParquetError::Format("chunk ended before all values".into()))?;
        let mut th = Thrift::new(rest);
        let hdr = parse_page_header(&mut th)?;
        if i64::try_from(hdr.num_values).unwrap_or(i64::MAX) > meta.num_values - seen {
            return Err(ParquetError::Format(format!(
                "column `{}` page holds more values than the chunk declares",
                leaf.name
            )));
        }
        let body_start = pos + th.position();
        let body_end = body_start
            .checked_add(hdr.compressed_size)
            .filter(|&e| e <= buf.len())
            .ok_or_else(|| ParquetError::Format("page extends past its chunk".into()))?;
        if hdr.uncompressed_size > MAX_PAGE_BYTES {
            return Err(ParquetError::Format("page uncompressed size is implausibly large".into()));
        }
        let body = decompress(meta.codec, &buf[body_start..body_end], hdr.uncompressed_size)?;
        pos = body_end;

        match hdr.kind {
            PAGE_DICTIONARY => {
                if hdr.encoding != ENC_PLAIN && hdr.encoding != ENC_PLAIN_DICTIONARY {
                    return Err(ParquetError::Unsupported(format!(
                        "dictionary page encoding {}",
                        hdr.encoding
                    )));
                }
                dict = Some(decode_plain(&body, meta.ptype, hdr.num_values)?);
            }
            PAGE_DATA => {
                let page = decode_data_page(&body, leaf, &hdr, dict.as_ref(), meta.ptype)?;
                seen += i64::try_from(hdr.num_values).unwrap_or(i64::MAX);
                out.append(page)?;
            }
            PAGE_DATA_V2 => {
                return Err(ParquetError::Unsupported("data page v2".into()));
            }
            other => {
                return Err(ParquetError::Format(format!("unknown page type {other}")));
            }
        }
    }

    if i64::try_from(out.len()).unwrap_or(i64::MAX) != meta.num_values {
        return Err(ParquetError::Format(format!(
            "column `{}` decoded {} values, metadata says {}",
            leaf.name,
            out.len(),
            meta.num_values
        )));
    }
    Ok(out)
}

fn decode_data_page(
    body: &[u8],
    leaf: &LeafInfo,
    hdr: &PageHeader,
    dict: Option<&Values>,
    ptype: i32,
) -> Result<Values, ParquetError> {
    let mut p = 0usize;
    let mut present = hdr.num_values;
    if leaf.optional {
        let len_bytes = body
            .get(..4)
            .ok_or_else(|| ParquetError::Format("truncated definition levels".into()))?;
        let len =
            u32::from_le_bytes([len_bytes[0], len_bytes[1], len_bytes[2], len_bytes[3]]) as usize;
        let levels = body
            .get(4..4 + len)
            .ok_or_else(|| ParquetError::Format("truncated definition levels".into()))?;
        let defs = decode_hybrid(levels, 1, hdr.num_values)?;
        let nulls = defs.iter().filter(|&&d| d == 0).count();
        if nulls > 0 {
            return Err(ParquetError::Unsupported(format!(
                "column `{}` contains {nulls} null value(s); null values have no representation in Column",
                leaf.name
            )));
        }
        present = defs.len();
        p = 4 + len;
    }
    let values = &body[p..];
    match hdr.encoding {
        ENC_PLAIN => decode_plain(values, ptype, present),
        ENC_PLAIN_DICTIONARY | ENC_RLE_DICTIONARY => {
            let dict = dict.ok_or_else(|| {
                ParquetError::Format(format!("column `{}` has no dictionary page", leaf.name))
            })?;
            let (&bit_width, rest) = values
                .split_first()
                .ok_or_else(|| ParquetError::Format("truncated dictionary indices".into()))?;
            let idx = decode_hybrid(rest, u32::from(bit_width), present)?;
            dict.gather(&idx)
        }
        other => Err(ParquetError::Unsupported(format!("encoding {other}"))),
    }
}

/// An open Parquet file. Column data is read on demand.
pub struct ParquetFile<R> {
    source: R,
    file_len: u64,
    footer: Footer,
    leaves: Vec<LeafInfo>,
}

impl<R> fmt::Debug for ParquetFile<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ParquetFile")
            .field("rows", &self.footer.num_rows)
            .field("row_groups", &self.footer.row_groups.len())
            .field("columns", &self.leaves)
            .finish_non_exhaustive()
    }
}

impl<R: Read + Seek> ParquetFile<R> {
    /// Reads and validates the footer of `source`.
    ///
    /// # Errors
    ///
    /// Returns [`ParquetError::Format`] if the magic bytes or footer are invalid,
    /// and [`ParquetError::Unsupported`] for nested schemas.
    pub fn open(mut source: R) -> Result<Self, ParquetError> {
        let len = source.seek(SeekFrom::End(0))?;
        if len < 12 {
            return Err(ParquetError::Format("file too short".into()));
        }
        source.seek(SeekFrom::Start(0))?;
        let mut head = [0u8; 4];
        source.read_exact(&mut head)?;
        if &head != MAGIC {
            return Err(ParquetError::Format("missing PAR1 header".into()));
        }
        source.seek(SeekFrom::End(-8))?;
        let mut tail = [0u8; 8];
        source.read_exact(&mut tail)?;
        if &tail[4..] != MAGIC {
            return Err(ParquetError::Format("missing PAR1 footer marker".into()));
        }
        let footer_len = u64::from(u32::from_le_bytes([tail[0], tail[1], tail[2], tail[3]]));
        if footer_len + 12 > len {
            return Err(ParquetError::Format("footer length exceeds file".into()));
        }
        source.seek(SeekFrom::Start(len - 8 - footer_len))?;
        let mut bytes = vec![0u8; footer_len as usize];
        source.read_exact(&mut bytes)?;
        let footer = parse_footer(&bytes)?;
        let leaves = schema_leaves(&footer.schema)?;
        for rg in &footer.row_groups {
            if rg.num_rows < 0 || rg.num_rows > footer.num_rows.max(0) {
                return Err(ParquetError::Format("row group row count is inconsistent".into()));
            }
        }
        Ok(Self { source, file_len: len, footer, leaves })
    }

    /// Total number of rows.
    pub fn num_rows(&self) -> u64 {
        self.footer.num_rows.max(0) as u64
    }

    /// Number of row groups.
    pub fn num_row_groups(&self) -> usize {
        self.footer.row_groups.len()
    }

    /// The top-level columns, in file order.
    pub fn columns(&self) -> &[LeafInfo] {
        &self.leaves
    }

    /// Reads one column across all row groups.
    ///
    /// # Errors
    ///
    /// Returns [`ParquetError::UnknownColumn`] if the name is not in the schema,
    /// and [`ParquetError::Unsupported`] or [`ParquetError::Format`] for data
    /// this reader cannot decode.
    pub fn read_column(&mut self, name: &str) -> Result<Column, ParquetError> {
        let index = self
            .leaves
            .iter()
            .position(|l| l.name == name)
            .ok_or_else(|| ParquetError::UnknownColumn(name.to_owned()))?;
        let leaf = self.leaves[index].clone();
        let mut acc = Values::empty(leaf.physical_type)?;
        for rg in &self.footer.row_groups {
            let meta = rg.columns.get(index).ok_or_else(|| {
                ParquetError::Format("row group has fewer columns than the schema".into())
            })?;
            let chunk = read_chunk(&mut self.source, self.file_len, &leaf, meta)?;
            acc.append(chunk)?;
        }
        Ok(acc.into_column())
    }

    /// Reads several columns (all of them if `names` is `None`) into a frame.
    ///
    /// # Errors
    ///
    /// As [`Self::read_column`], and [`ParquetError::Format`] if the columns
    /// disagree in length.
    pub fn read_frame(&mut self, names: Option<&[&str]>) -> Result<DataFrame, ParquetError> {
        let wanted: Vec<String> = match names {
            Some(ns) => ns.iter().map(|s| (*s).to_owned()).collect(),
            None => self.leaves.iter().map(|l| l.name.clone()).collect(),
        };
        let mut cols = Vec::with_capacity(wanted.len());
        for name in &wanted {
            cols.push((name.clone(), self.read_column(name)?));
        }
        let rows = self.num_rows() as usize;
        Ok(DataFrame::with_rows(rows, cols)?)
    }
}

/// Opens a Parquet file and reads the requested columns (all if `columns` is
/// `None`) into a [`DataFrame`].
///
/// # Errors
///
/// As [`ParquetFile::open`] and [`ParquetFile::read_frame`].
pub fn read_frame(
    path: impl AsRef<Path>,
    columns: Option<&[&str]>,
) -> Result<DataFrame, ParquetError> {
    let file = BufReader::new(File::open(path)?);
    ParquetFile::open(file)?.read_frame(columns)
}
