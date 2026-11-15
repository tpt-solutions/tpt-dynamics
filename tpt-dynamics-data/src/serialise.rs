//! Serialisation of histograms and profiles.
//!
//! Every histogram type offers two forms, both exact round-trips:
//!
//! * **Binary** (`to_bytes` / `from_bytes`): little-endian, with a 7-byte header
//!   (`b"TPTH"`, format version `u16`, kind `u8`), the payload, and a trailing
//!   CRC-32 over everything before it. Corrupted, truncated or mis-typed input is
//!   rejected with a [`SerialError`]; the decoder never panics and never allocates
//!   more than the input length.
//! * **JSON** (`to_json` / `from_json`): a self-describing object with
//!   `"format"`, `"version"` and `"type"` keys. Integers such as entry counts must be
//!   below 2^53 to survive the round-trip through JSON numbers.
//!
//! Histograms with non-finite contents cannot be written as JSON; their numbers
//! become `null`, which the reader rejects.

use crate::crc32::checksum;
use crate::hist1d::{Histogram1D, HistogramError};
use crate::hist2d::Histogram2D;
use crate::json::{self, JsonError, Value};
use crate::profile::{Profile1D, Profile2D, ProfileCell};
use core::fmt;

/// Magic bytes at the start of every binary histogram.
pub const MAGIC: &[u8; 4] = b"TPTH";
/// Current binary format version.
pub const FORMAT_VERSION: u16 = 1;
const JSON_FORMAT: &str = "tpt-dynamics-histogram";
const HEADER_LEN: usize = 7;
const CRC_LEN: usize = 4;

/// Error from decoding a serialised histogram.
#[derive(Debug)]
pub enum SerialError {
    /// The input ended before the expected data.
    Truncated,
    /// The input does not start with [`MAGIC`].
    BadMagic,
    /// The format version is not supported by this build.
    UnsupportedVersion(u16),
    /// The payload is a different histogram type than requested.
    WrongKind {
        /// Kind the caller asked for.
        expected: &'static str,
        /// Kind found in the input.
        found: u8,
    },
    /// The stored CRC-32 does not match the data.
    ChecksumMismatch,
    /// Extra bytes follow the payload.
    TrailingData,
    /// The decoded values do not form a valid histogram.
    Invalid(HistogramError),
    /// The JSON text could not be parsed.
    Json(JsonError),
    /// The JSON parsed but does not match the histogram schema.
    Schema(String),
}

impl fmt::Display for SerialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Truncated => write!(f, "histogram data is truncated"),
            Self::BadMagic => write!(f, "not a tpt-dynamics histogram (bad magic)"),
            Self::UnsupportedVersion(v) => write!(f, "unsupported histogram format version {v}"),
            Self::WrongKind { expected, found } => {
                write!(f, "expected a {expected} but found histogram kind {found}")
            }
            Self::ChecksumMismatch => write!(f, "histogram checksum mismatch (data corrupted)"),
            Self::TrailingData => write!(f, "unexpected data after histogram payload"),
            Self::Invalid(e) => write!(f, "invalid histogram: {e}"),
            Self::Json(e) => write!(f, "{e}"),
            Self::Schema(msg) => write!(f, "invalid histogram JSON: {msg}"),
        }
    }
}

impl std::error::Error for SerialError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Invalid(e) => Some(e),
            Self::Json(e) => Some(e),
            _ => None,
        }
    }
}

impl From<HistogramError> for SerialError {
    fn from(e: HistogramError) -> Self {
        Self::Invalid(e)
    }
}

impl From<JsonError> for SerialError {
    fn from(e: JsonError) -> Self {
        Self::Json(e)
    }
}

// ---------------------------------------------------------------------------
// Shared codec plumbing
// ---------------------------------------------------------------------------

/// Implemented by each serialisable type. Private: the public surface is the
/// inherent `to_bytes` / `from_bytes` / `to_json` / `from_json` methods.
trait Codec: Sized {
    const KIND: u8;
    const NAME: &'static str;
    fn enc(&self, e: &mut Enc);
    fn dec(d: &mut Dec<'_>) -> Result<Self, SerialError>;
    fn to_json_fields(&self) -> Vec<(String, Value)>;
    fn from_json_obj(v: &Value) -> Result<Self, SerialError>;
}

struct Enc {
    buf: Vec<u8>,
}

impl Enc {
    fn new(kind: u8) -> Self {
        let mut buf = Vec::with_capacity(64);
        buf.extend_from_slice(MAGIC);
        buf.extend_from_slice(&FORMAT_VERSION.to_le_bytes());
        buf.push(kind);
        Self { buf }
    }

    fn u64(&mut self, v: u64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn f64(&mut self, v: f64) {
        self.buf.extend_from_slice(&v.to_le_bytes());
    }

    fn f64s(&mut self, xs: &[f64]) {
        self.u64(xs.len() as u64);
        for &x in xs {
            self.f64(x);
        }
    }

    fn cells(&mut self, cells: &[ProfileCell]) {
        self.u64(cells.len() as u64);
        for c in cells {
            self.f64(c.sumw);
            self.f64(c.sumwy);
            self.f64(c.sumwy2);
            self.f64(c.sumw2);
            self.u64(c.entries);
        }
    }

    fn finish(mut self) -> Vec<u8> {
        let crc = checksum(&self.buf);
        self.buf.extend_from_slice(&crc.to_le_bytes());
        self.buf
    }
}

struct Dec<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Dec<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], SerialError> {
        let end = self.pos.checked_add(n).ok_or(SerialError::Truncated)?;
        if end > self.data.len() {
            return Err(SerialError::Truncated);
        }
        let s = &self.data[self.pos..end];
        self.pos = end;
        Ok(s)
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N], SerialError> {
        let s = self.take(N)?;
        let mut out = [0u8; N];
        out.copy_from_slice(s);
        Ok(out)
    }

    fn u64(&mut self) -> Result<u64, SerialError> {
        Ok(u64::from_le_bytes(self.array()?))
    }

    fn f64(&mut self) -> Result<f64, SerialError> {
        Ok(f64::from_le_bytes(self.array()?))
    }

    /// Read a length prefix, rejecting lengths that cannot fit in the remaining input.
    fn len(&mut self, elem_bytes: usize) -> Result<usize, SerialError> {
        let n = self.u64()?;
        let remaining = (self.data.len() - self.pos) / elem_bytes;
        if n > remaining as u64 {
            return Err(SerialError::Truncated);
        }
        Ok(n as usize)
    }

    fn f64s(&mut self) -> Result<Vec<f64>, SerialError> {
        let n = self.len(8)?;
        (0..n).map(|_| self.f64()).collect()
    }

    fn cells(&mut self) -> Result<Vec<ProfileCell>, SerialError> {
        let n = self.len(40)?;
        (0..n)
            .map(|_| {
                Ok(ProfileCell {
                    sumw: self.f64()?,
                    sumwy: self.f64()?,
                    sumwy2: self.f64()?,
                    sumw2: self.f64()?,
                    entries: self.u64()?,
                })
            })
            .collect()
    }
}

fn encode<T: Codec>(t: &T) -> Vec<u8> {
    let mut e = Enc::new(T::KIND);
    t.enc(&mut e);
    e.finish()
}

fn decode<T: Codec>(data: &[u8]) -> Result<T, SerialError> {
    if data.len() < HEADER_LEN + CRC_LEN {
        return Err(SerialError::Truncated);
    }
    if &data[..4] != MAGIC {
        return Err(SerialError::BadMagic);
    }
    let version = u16::from_le_bytes([data[4], data[5]]);
    if version != FORMAT_VERSION {
        return Err(SerialError::UnsupportedVersion(version));
    }
    let body_end = data.len() - CRC_LEN;
    let stored = u32::from_le_bytes([
        data[body_end],
        data[body_end + 1],
        data[body_end + 2],
        data[body_end + 3],
    ]);
    if checksum(&data[..body_end]) != stored {
        return Err(SerialError::ChecksumMismatch);
    }
    let kind = data[6];
    if kind != T::KIND {
        return Err(SerialError::WrongKind { expected: T::NAME, found: kind });
    }
    let mut d = Dec { data: &data[..body_end], pos: HEADER_LEN };
    let value = T::dec(&mut d)?;
    if d.pos != body_end {
        return Err(SerialError::TrailingData);
    }
    Ok(value)
}

fn to_json_text<T: Codec>(t: &T) -> String {
    let mut fields = vec![
        ("format".to_owned(), Value::String(JSON_FORMAT.to_owned())),
        ("version".to_owned(), Value::Number(f64::from(FORMAT_VERSION))),
        ("type".to_owned(), Value::String(T::NAME.to_owned())),
    ];
    fields.extend(t.to_json_fields());
    Value::Object(fields).write()
}

fn from_json_text<T: Codec>(text: &str) -> Result<T, SerialError> {
    let v = json::parse(text)?;
    if str_field(&v, "format")? != JSON_FORMAT {
        return Err(SerialError::Schema("format is not tpt-dynamics-histogram".into()));
    }
    if u64_field(&v, "version")? != u64::from(FORMAT_VERSION) {
        return Err(SerialError::Schema("unsupported version".into()));
    }
    let found = str_field(&v, "type")?;
    if found != T::NAME {
        return Err(SerialError::Schema(format!("expected type {}, found {found}", T::NAME)));
    }
    T::from_json_obj(&v)
}

// JSON field helpers.

fn field<'a>(v: &'a Value, key: &str) -> Result<&'a Value, SerialError> {
    v.get(key).ok_or_else(|| SerialError::Schema(format!("missing field \"{key}\"")))
}

fn str_field<'a>(v: &'a Value, key: &str) -> Result<&'a str, SerialError> {
    match field(v, key)? {
        Value::String(s) => Ok(s),
        _ => Err(SerialError::Schema(format!("\"{key}\" must be a string"))),
    }
}

fn num_field(v: &Value, key: &str) -> Result<f64, SerialError> {
    match field(v, key)? {
        Value::Number(x) => Ok(*x),
        _ => Err(SerialError::Schema(format!("\"{key}\" must be a number"))),
    }
}

fn u64_field(v: &Value, key: &str) -> Result<u64, SerialError> {
    let x = num_field(v, key)?;
    if x.fract() != 0.0 || !(0.0..=9_007_199_254_740_992.0).contains(&x) {
        return Err(SerialError::Schema(format!("\"{key}\" must be a non-negative integer")));
    }
    Ok(x as u64)
}

fn f64s_field(v: &Value, key: &str) -> Result<Vec<f64>, SerialError> {
    match field(v, key)? {
        Value::Array(items) => items
            .iter()
            .map(|item| match item {
                Value::Number(x) => Ok(*x),
                _ => Err(SerialError::Schema(format!("\"{key}\" must contain only numbers"))),
            })
            .collect(),
        _ => Err(SerialError::Schema(format!("\"{key}\" must be an array"))),
    }
}

fn u64s_field(v: &Value, key: &str, n: usize) -> Result<Vec<u64>, SerialError> {
    let arr = f64s_field(v, key)?;
    if arr.len() != n {
        return Err(SerialError::Schema(format!("\"{key}\" has the wrong length")));
    }
    arr.into_iter()
        .map(|x| {
            if x.fract() == 0.0 && (0.0..=9_007_199_254_740_992.0).contains(&x) {
                Ok(x as u64)
            } else {
                Err(SerialError::Schema(format!("\"{key}\" must contain integers")))
            }
        })
        .collect()
}

fn nums(xs: &[f64]) -> Value {
    Value::Array(xs.iter().map(|&x| Value::Number(x)).collect())
}

fn ints(xs: impl Iterator<Item = u64>) -> Value {
    Value::Array(xs.map(|x| Value::Number(x as f64)).collect())
}

fn cells_from_json(v: &Value) -> Result<Vec<ProfileCell>, SerialError> {
    let n = match field(v, "sumw")? {
        Value::Array(a) => a.len(),
        _ => return Err(SerialError::Schema("\"weights\" must be an array".into())),
    };
    let weights = f64s_field(v, "sumw")?;
    let weighted_y = f64s_field(v, "sumwy")?;
    let weighted_y_sq = f64s_field(v, "sumwy2")?;
    let weights_sq = f64s_field(v, "sumw2")?;
    let entries = u64s_field(v, "entries", n)?;
    if weighted_y.len() != n || weighted_y_sq.len() != n || weights_sq.len() != n {
        return Err(SerialError::Schema("profile arrays have mismatched lengths".into()));
    }
    Ok((0..n)
        .map(|i| ProfileCell {
            sumw: weights[i],
            sumwy: weighted_y[i],
            sumwy2: weighted_y_sq[i],
            sumw2: weights_sq[i],
            entries: entries[i],
        })
        .collect())
}

fn cells_to_json(cells: &[ProfileCell]) -> Vec<(String, Value)> {
    vec![
        ("sumw".into(), nums(&cells.iter().map(|c| c.sumw).collect::<Vec<_>>())),
        ("sumwy".into(), nums(&cells.iter().map(|c| c.sumwy).collect::<Vec<_>>())),
        ("sumwy2".into(), nums(&cells.iter().map(|c| c.sumwy2).collect::<Vec<_>>())),
        ("sumw2".into(), nums(&cells.iter().map(|c| c.sumw2).collect::<Vec<_>>())),
        ("entries".into(), ints(cells.iter().map(|c| c.entries))),
    ]
}

// ---------------------------------------------------------------------------
// Histogram1D
// ---------------------------------------------------------------------------

impl Codec for Histogram1D {
    const KIND: u8 = 1;
    const NAME: &'static str = "histogram1d";

    fn enc(&self, e: &mut Enc) {
        e.f64s(&self.edges);
        e.f64s(&self.sumw);
        e.f64s(&self.sumw2);
        for x in self.flow() {
            e.f64(x);
        }
        e.u64(self.entries);
        e.u64(self.nan_fills);
        for x in self.moments() {
            e.f64(x);
        }
    }

    fn dec(d: &mut Dec<'_>) -> Result<Self, SerialError> {
        let edges = d.f64s()?;
        let sumw = d.f64s()?;
        let sumw2 = d.f64s()?;
        let flow = [d.f64()?, d.f64()?, d.f64()?, d.f64()?];
        let entries = d.u64()?;
        let nan_fills = d.u64()?;
        let moments = [d.f64()?, d.f64()?, d.f64()?];
        Ok(Self::from_parts(edges, sumw, sumw2, flow, entries, nan_fills, moments)?)
    }

    fn to_json_fields(&self) -> Vec<(String, Value)> {
        let flow = self.flow();
        let m = self.moments();
        vec![
            ("edges".into(), nums(&self.edges)),
            ("contents".into(), nums(&self.sumw)),
            ("sumw2".into(), nums(&self.sumw2)),
            ("underflow".into(), Value::Number(flow[0])),
            ("overflow".into(), Value::Number(flow[1])),
            ("underflow_sumw2".into(), Value::Number(flow[2])),
            ("overflow_sumw2".into(), Value::Number(flow[3])),
            ("entries".into(), Value::Number(self.entries as f64)),
            ("nan_fills".into(), Value::Number(self.nan_fills as f64)),
            ("moments".into(), nums(&m)),
        ]
    }

    fn from_json_obj(v: &Value) -> Result<Self, SerialError> {
        let edges = f64s_field(v, "edges")?;
        let sumw = f64s_field(v, "contents")?;
        let sumw2 = f64s_field(v, "sumw2")?;
        let flow = [
            num_field(v, "underflow")?,
            num_field(v, "overflow")?,
            num_field(v, "underflow_sumw2")?,
            num_field(v, "overflow_sumw2")?,
        ];
        let entries = u64_field(v, "entries")?;
        let nan_fills = u64_field(v, "nan_fills")?;
        let m = f64s_field(v, "moments")?;
        if m.len() != 3 {
            return Err(SerialError::Schema("\"moments\" must have 3 values".into()));
        }
        Ok(Self::from_parts(edges, sumw, sumw2, flow, entries, nan_fills, [m[0], m[1], m[2]])?)
    }
}

impl Histogram1D {
    /// Encode as the binary format (see the [module docs](crate::serialise)).
    pub fn to_bytes(&self) -> Vec<u8> {
        encode(self)
    }

    /// Decode from the binary format, verifying magic, version, kind and checksum.
    pub fn from_bytes(data: &[u8]) -> Result<Self, SerialError> {
        decode(data)
    }

    /// Encode as a JSON document.
    pub fn to_json(&self) -> String {
        to_json_text(self)
    }

    /// Decode from a JSON document produced by [`Self::to_json`].
    pub fn from_json(text: &str) -> Result<Self, SerialError> {
        from_json_text(text)
    }
}

// ---------------------------------------------------------------------------
// Histogram2D
// ---------------------------------------------------------------------------

impl Codec for Histogram2D {
    const KIND: u8 = 2;
    const NAME: &'static str = "histogram2d";

    fn enc(&self, e: &mut Enc) {
        e.f64s(&self.xedges);
        e.f64s(&self.yedges);
        e.f64s(&self.sumw);
        e.f64s(&self.sumw2);
        e.u64(self.entries);
        e.u64(self.nan_fills);
    }

    fn dec(d: &mut Dec<'_>) -> Result<Self, SerialError> {
        let xedges = d.f64s()?;
        let yedges = d.f64s()?;
        let sumw = d.f64s()?;
        let sumw2 = d.f64s()?;
        let entries = d.u64()?;
        let nan_fills = d.u64()?;
        Ok(Self::from_parts(xedges, yedges, sumw, sumw2, entries, nan_fills)?)
    }

    fn to_json_fields(&self) -> Vec<(String, Value)> {
        vec![
            ("x_edges".into(), nums(&self.xedges)),
            ("y_edges".into(), nums(&self.yedges)),
            ("contents".into(), nums(&self.sumw)),
            ("sumw2".into(), nums(&self.sumw2)),
            ("entries".into(), Value::Number(self.entries as f64)),
            ("nan_fills".into(), Value::Number(self.nan_fills as f64)),
        ]
    }

    fn from_json_obj(v: &Value) -> Result<Self, SerialError> {
        Ok(Self::from_parts(
            f64s_field(v, "x_edges")?,
            f64s_field(v, "y_edges")?,
            f64s_field(v, "contents")?,
            f64s_field(v, "sumw2")?,
            u64_field(v, "entries")?,
            u64_field(v, "nan_fills")?,
        )?)
    }
}

impl Histogram2D {
    /// Encode as the binary format (see the [module docs](crate::serialise)).
    pub fn to_bytes(&self) -> Vec<u8> {
        encode(self)
    }

    /// Decode from the binary format, verifying magic, version, kind and checksum.
    pub fn from_bytes(data: &[u8]) -> Result<Self, SerialError> {
        decode(data)
    }

    /// Encode as a JSON document.
    pub fn to_json(&self) -> String {
        to_json_text(self)
    }

    /// Decode from a JSON document produced by [`Self::to_json`].
    pub fn from_json(text: &str) -> Result<Self, SerialError> {
        from_json_text(text)
    }
}

// ---------------------------------------------------------------------------
// Profile1D / Profile2D
// ---------------------------------------------------------------------------

impl Codec for Profile1D {
    const KIND: u8 = 3;
    const NAME: &'static str = "profile1d";

    fn enc(&self, e: &mut Enc) {
        e.f64s(&self.edges);
        e.cells(&self.cells);
        e.f64(self.underflow);
        e.f64(self.overflow);
        e.u64(self.nan_fills);
    }

    fn dec(d: &mut Dec<'_>) -> Result<Self, SerialError> {
        let edges = d.f64s()?;
        let cells = d.cells()?;
        let underflow = d.f64()?;
        let overflow = d.f64()?;
        let nan_fills = d.u64()?;
        Ok(Self::from_parts(edges, cells, underflow, overflow, nan_fills)?)
    }

    fn to_json_fields(&self) -> Vec<(String, Value)> {
        let mut fields = vec![("edges".into(), nums(&self.edges))];
        fields.extend(cells_to_json(&self.cells));
        fields.push(("underflow".into(), Value::Number(self.underflow)));
        fields.push(("overflow".into(), Value::Number(self.overflow)));
        fields.push(("nan_fills".into(), Value::Number(self.nan_fills as f64)));
        fields
    }

    fn from_json_obj(v: &Value) -> Result<Self, SerialError> {
        Ok(Self::from_parts(
            f64s_field(v, "edges")?,
            cells_from_json(v)?,
            num_field(v, "underflow")?,
            num_field(v, "overflow")?,
            u64_field(v, "nan_fills")?,
        )?)
    }
}

impl Profile1D {
    /// Encode as the binary format (see the [module docs](crate::serialise)).
    pub fn to_bytes(&self) -> Vec<u8> {
        encode(self)
    }

    /// Decode from the binary format, verifying magic, version, kind and checksum.
    pub fn from_bytes(data: &[u8]) -> Result<Self, SerialError> {
        decode(data)
    }

    /// Encode as a JSON document.
    pub fn to_json(&self) -> String {
        to_json_text(self)
    }

    /// Decode from a JSON document produced by [`Self::to_json`].
    pub fn from_json(text: &str) -> Result<Self, SerialError> {
        from_json_text(text)
    }
}

impl Codec for Profile2D {
    const KIND: u8 = 4;
    const NAME: &'static str = "profile2d";

    fn enc(&self, e: &mut Enc) {
        e.f64s(&self.xedges);
        e.f64s(&self.yedges);
        e.cells(&self.cells);
        e.f64(self.outside);
        e.u64(self.nan_fills);
    }

    fn dec(d: &mut Dec<'_>) -> Result<Self, SerialError> {
        let xedges = d.f64s()?;
        let yedges = d.f64s()?;
        let cells = d.cells()?;
        let outside = d.f64()?;
        let nan_fills = d.u64()?;
        Ok(Self::from_parts(xedges, yedges, cells, outside, nan_fills)?)
    }

    fn to_json_fields(&self) -> Vec<(String, Value)> {
        let mut fields =
            vec![("x_edges".into(), nums(&self.xedges)), ("y_edges".into(), nums(&self.yedges))];
        fields.extend(cells_to_json(&self.cells));
        fields.push(("outside".into(), Value::Number(self.outside)));
        fields.push(("nan_fills".into(), Value::Number(self.nan_fills as f64)));
        fields
    }

    fn from_json_obj(v: &Value) -> Result<Self, SerialError> {
        Ok(Self::from_parts(
            f64s_field(v, "x_edges")?,
            f64s_field(v, "y_edges")?,
            cells_from_json(v)?,
            num_field(v, "outside")?,
            u64_field(v, "nan_fills")?,
        )?)
    }
}

impl Profile2D {
    /// Encode as the binary format (see the [module docs](crate::serialise)).
    pub fn to_bytes(&self) -> Vec<u8> {
        encode(self)
    }

    /// Decode from the binary format, verifying magic, version, kind and checksum.
    pub fn from_bytes(data: &[u8]) -> Result<Self, SerialError> {
        decode(data)
    }

    /// Encode as a JSON document.
    pub fn to_json(&self) -> String {
        to_json_text(self)
    }

    /// Decode from a JSON document produced by [`Self::to_json`].
    pub fn from_json(text: &str) -> Result<Self, SerialError> {
        from_json_text(text)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_h1() -> Histogram1D {
        let mut h = Histogram1D::variable(&[0.0, 0.1, 1.0, 10.0]).unwrap();
        for i in 0..200 {
            h.fill_weighted(f64::from(i) * 0.06 - 0.5, 1.0 + f64::from(i % 7) * 0.3);
        }
        h.fill(f64::NAN);
        h
    }

    fn sample_h2() -> Histogram2D {
        let mut h = Histogram2D::uniform(3, 0.0, 3.0, 2, -1.0, 1.0).unwrap();
        for i in 0..100 {
            h.fill_weighted(
                f64::from(i % 5) * 0.8 - 0.5,
                f64::from(i % 3) - 1.2,
                0.5 + f64::from(i),
            );
        }
        h
    }

    fn sample_p1() -> Profile1D {
        let mut p = Profile1D::uniform(4, 0.0, 4.0).unwrap();
        for i in 0..50 {
            p.fill_weighted(f64::from(i % 5), f64::from(i) * 0.3, 1.0 + f64::from(i % 2));
        }
        p
    }

    fn sample_p2() -> Profile2D {
        let mut p = Profile2D::uniform(2, 0.0, 2.0, 2, 0.0, 2.0).unwrap();
        for i in 0..30 {
            p.fill(f64::from(i % 3) * 0.9, f64::from(i % 2) + 0.5, f64::from(i) * 1.1);
        }
        p
    }

    #[test]
    fn binary_round_trip_all_types() {
        let h1 = sample_h1();
        assert_eq!(Histogram1D::from_bytes(&h1.to_bytes()).unwrap(), h1);
        let h2 = sample_h2();
        assert_eq!(Histogram2D::from_bytes(&h2.to_bytes()).unwrap(), h2);
        let p1 = sample_p1();
        assert_eq!(Profile1D::from_bytes(&p1.to_bytes()).unwrap(), p1);
        let p2 = sample_p2();
        assert_eq!(Profile2D::from_bytes(&p2.to_bytes()).unwrap(), p2);
    }

    #[test]
    fn json_round_trip_all_types() {
        let h1 = sample_h1();
        assert_eq!(Histogram1D::from_json(&h1.to_json()).unwrap(), h1);
        let h2 = sample_h2();
        assert_eq!(Histogram2D::from_json(&h2.to_json()).unwrap(), h2);
        let p1 = sample_p1();
        assert_eq!(Profile1D::from_json(&p1.to_json()).unwrap(), p1);
        let p2 = sample_p2();
        assert_eq!(Profile2D::from_json(&p2.to_json()).unwrap(), p2);
    }

    #[test]
    fn wrong_kind_is_rejected() {
        let bytes = sample_h1().to_bytes();
        assert!(matches!(
            Histogram2D::from_bytes(&bytes),
            Err(SerialError::WrongKind { found: 1, .. })
        ));
        let json = sample_h1().to_json();
        assert!(Profile1D::from_json(&json).is_err());
    }

    #[test]
    fn every_single_byte_corruption_is_detected() {
        let bytes = sample_h2().to_bytes();
        for i in 0..bytes.len() {
            let mut bad = bytes.clone();
            bad[i] ^= 0x40;
            assert!(Histogram2D::from_bytes(&bad).is_err(), "corruption at byte {i} accepted");
        }
    }

    #[test]
    fn every_truncation_is_rejected() {
        let bytes = sample_p2().to_bytes();
        for n in 0..bytes.len() {
            assert!(Profile2D::from_bytes(&bytes[..n]).is_err(), "prefix {n} accepted");
        }
    }

    #[test]
    fn bad_magic_and_version_are_rejected() {
        let mut bytes = sample_h1().to_bytes();
        bytes[0] = b'X';
        assert!(matches!(Histogram1D::from_bytes(&bytes), Err(SerialError::BadMagic)));
        let mut bytes = sample_h1().to_bytes();
        bytes[4] = 9;
        // Fix up the CRC so that the version check, not the checksum, is what fails.
        let body = bytes.len() - 4;
        let crc = checksum(&bytes[..body]).to_le_bytes();
        bytes[body..].copy_from_slice(&crc);
        assert!(matches!(Histogram1D::from_bytes(&bytes), Err(SerialError::UnsupportedVersion(9))));
    }

    #[test]
    fn json_schema_violations_are_rejected() {
        let bad = [
            "",
            "{}",
            r#"{"format":"other","version":1,"type":"histogram1d"}"#,
            r#"{"format":"tpt-dynamics-histogram","version":1,"type":"histogram1d"}"#,
            r#"{"format":"tpt-dynamics-histogram","version":1,"type":"histogram2d","x_edges":[0,1],"y_edges":[0,1],"contents":[1],"sumw2":[1],"entries":1.5,"nan_fills":0}"#,
            r#"{"format":"tpt-dynamics-histogram","version":1,"type":"histogram2d","x_edges":[1,0],"y_edges":[0,1],"contents":[1],"sumw2":[1],"entries":1,"nan_fills":0}"#,
        ];
        for text in bad {
            assert!(Histogram2D::from_json(text).is_err(), "accepted {text}");
        }
    }
}
