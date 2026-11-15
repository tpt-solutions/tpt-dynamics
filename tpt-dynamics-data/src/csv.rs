//! Streaming CSV reader for files too large to load into memory.
//!
//! [`CsvReader`] pulls one line at a time from any [`BufRead`]. The line buffer is
//! reused and capped at [`CsvOptions::max_line_bytes`], so memory use does not
//! grow with file size. Each record is returned as owned fields; the caller
//! chooses which columns to parse, by index or by header name.
//!
//! Supported syntax: a single-byte delimiter, RFC 4180 double-quoted fields with
//! `""` as an escaped quote, `\n` or `\r\n` line endings, and blank lines
//! (skipped). A quoted field may not span lines. Such a line is reported as an
//! error with its line number, which keeps the reader's buffer bounded.
//!
//! ```
//! use std::io::BufReader;
//! use tpt_dynamics_data::csv::{CsvOptions, CsvReader};
//!
//! let data = "pt,eta\n12.5,0.3\n7.0,-1.2\n";
//! let mut rdr = CsvReader::new(BufReader::new(data.as_bytes()), CsvOptions::default()).unwrap();
//! let cols = rdr.columns(&["pt", "eta"]).unwrap();
//! let mut total = 0.0;
//! for rec in rdr {
//!     total += rec.unwrap().f64(cols[0]).unwrap();
//! }
//! assert_eq!(total, 19.5);
//! ```

use core::fmt;
use std::io::{self, BufRead, Read};

/// Options controlling how a CSV stream is parsed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CsvOptions {
    /// Field separator byte. Must not be `"`, `\n` or `\r`.
    pub delimiter: u8,
    /// Whether the first non-blank line is a header.
    pub has_header: bool,
    /// Maximum length of one physical line, including its terminator.
    pub max_line_bytes: usize,
}

impl Default for CsvOptions {
    fn default() -> Self {
        Self { delimiter: b',', has_header: true, max_line_bytes: 1 << 20 }
    }
}

/// Errors from reading or parsing CSV data. Parse errors carry the 1-based
/// line number of the offending line.
#[derive(Debug)]
pub enum CsvError {
    /// The underlying reader failed.
    Io(io::Error),
    /// The options were unusable (for example, a quote used as delimiter).
    InvalidOptions(&'static str),
    /// A line exceeded `max_line_bytes`.
    LineTooLong {
        /// Line number.
        line: u64,
        /// Configured limit.
        limit: usize,
    },
    /// A line was not valid UTF-8.
    InvalidUtf8 {
        /// Line number.
        line: u64,
    },
    /// Quoting was malformed (unterminated quote, or text after a closing quote).
    Malformed {
        /// Line number.
        line: u64,
        /// What went wrong.
        message: &'static str,
    },
    /// A record had a different number of fields than the header or first record.
    FieldCount {
        /// Line number.
        line: u64,
        /// Expected number of fields.
        expected: usize,
        /// Fields found.
        found: usize,
    },
    /// A requested column name is not in the header.
    UnknownColumn(String),
    /// A requested field index does not exist in the record.
    NoSuchField {
        /// Line number.
        line: u64,
        /// Requested index.
        index: usize,
    },
    /// A field could not be parsed as the requested type.
    BadValue {
        /// Line number.
        line: u64,
        /// Field index.
        index: usize,
        /// Raw field text (truncated to 64 bytes).
        value: String,
        /// Type name that was expected.
        expected: &'static str,
    },
}

impl fmt::Display for CsvError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error reading CSV: {e}"),
            Self::InvalidOptions(m) => write!(f, "invalid CSV options: {m}"),
            Self::LineTooLong { line, limit } => {
                write!(f, "line {line}: longer than {limit} bytes")
            }
            Self::InvalidUtf8 { line } => write!(f, "line {line}: not valid UTF-8"),
            Self::Malformed { line, message } => write!(f, "line {line}: {message}"),
            Self::FieldCount { line, expected, found } => {
                write!(f, "line {line}: expected {expected} fields, found {found}")
            }
            Self::UnknownColumn(name) => write!(f, "no column named \"{name}\""),
            Self::NoSuchField { line, index } => {
                write!(f, "line {line}: no field at index {index}")
            }
            Self::BadValue { line, index, value, expected } => {
                write!(f, "line {line}, field {index}: cannot parse \"{value}\" as {expected}")
            }
        }
    }
}

impl std::error::Error for CsvError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for CsvError {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

/// One parsed record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CsvRecord {
    line: u64,
    fields: Vec<String>,
}

impl CsvRecord {
    /// 1-based line number this record came from.
    #[inline]
    pub fn line(&self) -> u64 {
        self.line
    }

    /// Number of fields.
    #[inline]
    pub fn len(&self) -> usize {
        self.fields.len()
    }

    /// Whether the record has no fields.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }

    /// Raw text of field `i`, if present.
    #[inline]
    pub fn get(&self, i: usize) -> Option<&str> {
        self.fields.get(i).map(String::as_str)
    }

    /// All fields as text.
    #[inline]
    pub fn fields(&self) -> &[String] {
        &self.fields
    }

    fn raw(&self, i: usize) -> Result<&str, CsvError> {
        self.get(i).ok_or(CsvError::NoSuchField { line: self.line, index: i })
    }

    fn bad(&self, i: usize, value: &str, expected: &'static str) -> CsvError {
        let mut shown = value.to_owned();
        if shown.len() > 64 {
            let mut cut = 64;
            while !shown.is_char_boundary(cut) {
                cut -= 1;
            }
            shown.truncate(cut);
        }
        CsvError::BadValue { line: self.line, index: i, value: shown, expected }
    }

    /// Parse field `i` as a finite `f64`. Surrounding whitespace is ignored.
    pub fn f64(&self, i: usize) -> Result<f64, CsvError> {
        let text = self.raw(i)?.trim();
        match text.parse::<f64>() {
            Ok(x) if x.is_finite() => Ok(x),
            _ => Err(self.bad(i, text, "finite f64")),
        }
    }

    /// Parse field `i` as a `u64`.
    pub fn u64(&self, i: usize) -> Result<u64, CsvError> {
        let text = self.raw(i)?.trim();
        text.parse::<u64>().map_err(|_| self.bad(i, text, "u64"))
    }

    /// Parse field `i` as an `i64`.
    pub fn i64(&self, i: usize) -> Result<i64, CsvError> {
        let text = self.raw(i)?.trim();
        text.parse::<i64>().map_err(|_| self.bad(i, text, "i64"))
    }
}

/// Iterator-based CSV reader over any buffered source.
#[derive(Debug)]
pub struct CsvReader<R: BufRead> {
    reader: R,
    opts: CsvOptions,
    header: Option<Vec<String>>,
    expected_fields: Option<usize>,
    line_no: u64,
    buf: Vec<u8>,
    finished: bool,
}

impl<R: BufRead> CsvReader<R> {
    /// Create a reader. If [`CsvOptions::has_header`] is set, the header line is
    /// read immediately so that [`Self::columns`] is available.
    pub fn new(reader: R, opts: CsvOptions) -> Result<Self, CsvError> {
        if matches!(opts.delimiter, b'"' | b'\n' | b'\r') {
            return Err(CsvError::InvalidOptions("delimiter may not be a quote or newline"));
        }
        if opts.max_line_bytes == 0 {
            return Err(CsvError::InvalidOptions("max_line_bytes must be positive"));
        }
        let mut rdr = Self {
            reader,
            opts,
            header: None,
            expected_fields: None,
            line_no: 0,
            buf: Vec::with_capacity(256),
            finished: false,
        };
        if rdr.opts.has_header {
            match rdr.next_fields()? {
                Some((_, fields)) => {
                    rdr.expected_fields = Some(fields.len());
                    rdr.header = Some(fields);
                }
                None => rdr.finished = true,
            }
        }
        Ok(rdr)
    }

    /// Header names, if the file has a header.
    pub fn header(&self) -> Option<&[String]> {
        self.header.as_deref()
    }

    /// Index of the column called `name`, if the header contains it.
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.header.as_ref()?.iter().position(|h| h == name)
    }

    /// Resolve several column names to indices in one call.
    pub fn columns(&self, names: &[&str]) -> Result<Vec<usize>, CsvError> {
        names
            .iter()
            .map(|n| self.column_index(n).ok_or_else(|| CsvError::UnknownColumn((*n).to_owned())))
            .collect()
    }

    /// Number of physical lines consumed so far (1-based line of the last record).
    pub fn lines_read(&self) -> u64 {
        self.line_no
    }

    /// Reads up to `max_rows` records, splitting their fields on the rayon pool.
    ///
    /// Lines are read in order on the calling thread, and the records come back
    /// in file order. The result is the same as collecting the iterator, with
    /// the same field-count check. Returns an empty vector at end of input. An
    /// error discards the records of this batch: the caller should stop.
    ///
    /// # Errors
    ///
    /// Returns the first read, UTF-8, quoting or field-count error in the batch.
    pub fn next_batch(&mut self, max_rows: usize) -> Result<Vec<CsvRecord>, CsvError> {
        use rayon::prelude::*;

        if self.finished || max_rows == 0 {
            return Ok(Vec::new());
        }
        let mut lines: Vec<(u64, String)> = Vec::with_capacity(max_rows.min(1 << 20));
        while lines.len() < max_rows {
            if !self.read_line()? {
                self.finished = true;
                break;
            }
            if self.buf.is_empty() {
                continue;
            }
            let line = self.line_no;
            // Move the bytes out rather than copying them; read_line clears the buffer.
            let bytes = std::mem::take(&mut self.buf);
            let text = String::from_utf8(bytes).map_err(|_| CsvError::InvalidUtf8 { line })?;
            lines.push((line, text));
        }
        let delim = self.opts.delimiter;
        let split: Vec<Result<CsvRecord, CsvError>> = lines
            .into_par_iter()
            .map(|(line, text)| Ok(CsvRecord { line, fields: split_fields(&text, delim, line)? }))
            .collect();
        let mut out = Vec::with_capacity(split.len());
        for record in split {
            let record = record?;
            match self.expected_fields {
                None => self.expected_fields = Some(record.fields.len()),
                Some(expected) if expected != record.fields.len() => {
                    self.finished = true;
                    return Err(CsvError::FieldCount {
                        line: record.line,
                        expected,
                        found: record.fields.len(),
                    });
                }
                Some(_) => {}
            }
            out.push(record);
        }
        Ok(out)
    }

    /// Read the next line into `self.buf`, without its terminator.
    /// Returns `false` at end of input.
    fn read_line(&mut self) -> Result<bool, CsvError> {
        self.buf.clear();
        let limit = self.opts.max_line_bytes;
        // Reading through `take` caps how much one line can pull into the buffer.
        let mut limited = Read::take(&mut self.reader, limit as u64 + 1);
        let n = limited.read_until(b'\n', &mut self.buf)?;
        if n == 0 {
            return Ok(false);
        }
        self.line_no += 1;
        if self.buf.len() > limit {
            return Err(CsvError::LineTooLong { line: self.line_no, limit });
        }
        if self.buf.last() == Some(&b'\n') {
            self.buf.pop();
            if self.buf.last() == Some(&b'\r') {
                self.buf.pop();
            }
        }
        Ok(true)
    }

    /// Read and split the next non-blank line.
    fn next_fields(&mut self) -> Result<Option<(u64, Vec<String>)>, CsvError> {
        loop {
            if !self.read_line()? {
                return Ok(None);
            }
            if self.buf.is_empty() {
                continue;
            }
            let line = self.line_no;
            let text =
                core::str::from_utf8(&self.buf).map_err(|_| CsvError::InvalidUtf8 { line })?;
            let fields = split_fields(text, self.opts.delimiter, line)?;
            return Ok(Some((line, fields)));
        }
    }
}

impl<R: BufRead> Iterator for CsvReader<R> {
    type Item = Result<CsvRecord, CsvError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        match self.next_fields() {
            Ok(Some((line, fields))) => {
                match self.expected_fields {
                    None => self.expected_fields = Some(fields.len()),
                    Some(expected) if expected != fields.len() => {
                        self.finished = true;
                        return Some(Err(CsvError::FieldCount {
                            line,
                            expected,
                            found: fields.len(),
                        }));
                    }
                    Some(_) => {}
                }
                Some(Ok(CsvRecord { line, fields }))
            }
            Ok(None) => {
                self.finished = true;
                None
            }
            Err(e) => {
                // A read error or malformed line ends the stream. Continuing could
                // yield misaligned data.
                self.finished = true;
                Some(Err(e))
            }
        }
    }
}

/// Split one line into fields, honouring RFC 4180 quoting.
fn split_fields(line: &str, delim: u8, line_no: u64) -> Result<Vec<String>, CsvError> {
    let bytes = line.as_bytes();
    let mut fields = Vec::new();
    let mut start = 0;
    loop {
        let (field, end) = parse_field(bytes, start, delim, line_no)?;
        fields.push(field);
        if end < bytes.len() {
            // `end` is the delimiter. Another field follows, possibly empty.
            start = end + 1;
        } else {
            return Ok(fields);
        }
    }
}

/// Parse the field starting at `start`. Returns its text and the index of the
/// delimiter (or the line length) that terminates it.
fn parse_field(
    bytes: &[u8],
    start: usize,
    delim: u8,
    line_no: u64,
) -> Result<(String, usize), CsvError> {
    if bytes.get(start) == Some(&b'"') {
        let mut field: Vec<u8> = Vec::new();
        let mut i = start + 1;
        loop {
            match bytes.get(i) {
                None => {
                    return Err(CsvError::Malformed {
                        line: line_no,
                        message: "unterminated quoted field",
                    })
                }
                Some(&b'"') if bytes.get(i + 1) == Some(&b'"') => {
                    field.push(b'"');
                    i += 2;
                }
                Some(&b'"') => {
                    i += 1;
                    break;
                }
                Some(&c) => {
                    field.push(c);
                    i += 1;
                }
            }
        }
        match bytes.get(i) {
            None => {}
            Some(&c) if c == delim => {}
            Some(_) => {
                return Err(CsvError::Malformed {
                    line: line_no,
                    message: "text after closing quote",
                })
            }
        }
        // Only ASCII bytes were removed from a valid UTF-8 slice, so this holds.
        let text = String::from_utf8(field).map_err(|_| CsvError::InvalidUtf8 { line: line_no })?;
        Ok((text, i))
    } else {
        let mut i = start;
        while i < bytes.len() && bytes[i] != delim {
            if bytes[i] == b'"' {
                return Err(CsvError::Malformed {
                    line: line_no,
                    message: "quote inside unquoted field",
                });
            }
            i += 1;
        }
        // Boundaries are at ASCII bytes, so they fall on char boundaries.
        let text = String::from_utf8_lossy(&bytes[start..i]).into_owned();
        Ok((text, i))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::SplitMix64;
    use std::io::BufReader;

    fn reader(text: &str) -> CsvReader<BufReader<&[u8]>> {
        CsvReader::new(BufReader::new(text.as_bytes()), CsvOptions::default()).unwrap()
    }

    #[test]
    fn header_and_typed_access() {
        let rdr = reader("event,pt,eta,label\n1,12.5,0.25,\"a,b\"\n2,7,-1e-1,plain\n");
        assert_eq!(rdr.header().unwrap(), ["event", "pt", "eta", "label"]);
        let cols = rdr.columns(&["pt", "eta"]).unwrap();
        assert_eq!(cols, vec![1, 2]);
        let recs: Vec<CsvRecord> = rdr.map(Result::unwrap).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].f64(1).unwrap(), 12.5);
        assert_eq!(recs[0].get(3), Some("a,b"));
        assert_eq!(recs[1].u64(0).unwrap(), 2);
        assert_eq!(recs[1].f64(2).unwrap(), -0.1);
        assert_eq!(recs[1].line(), 3);
    }

    #[test]
    fn crlf_blank_lines_and_escaped_quotes() {
        let rdr = CsvReader::new(
            BufReader::new("a,b\r\n\r\n\"x\"\"y\",2\r\n".as_bytes()),
            CsvOptions::default(),
        )
        .unwrap();
        let recs: Vec<CsvRecord> = rdr.map(Result::unwrap).collect();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].get(0), Some("x\"y"));
        assert_eq!(recs[0].line(), 3);
    }

    #[test]
    fn trailing_delimiter_gives_empty_field() {
        let mut rdr =
            CsvReader::new(BufReader::new("a,b\n1,\n".as_bytes()), CsvOptions::default()).unwrap();
        let rec = rdr.next().unwrap().unwrap();
        assert_eq!(rec.fields(), ["1", ""]);
        assert!(matches!(rec.f64(1), Err(CsvError::BadValue { .. })));
    }

    #[test]
    fn errors_carry_line_numbers() {
        let mut rdr = reader("a,b\n1,2\n3,oops\n");
        assert!(rdr.next().unwrap().is_ok());
        let rec = rdr.next().unwrap().unwrap();
        let err = rec.f64(1).unwrap_err();
        assert!(err.to_string().contains("line 3"), "{err}");
        assert!(rdr.next().is_none());
    }

    #[test]
    fn field_count_mismatch_is_reported() {
        let mut rdr = reader("a,b\n1,2\n1,2,3\n");
        assert!(rdr.next().unwrap().is_ok());
        match rdr.next().unwrap() {
            Err(CsvError::FieldCount { line: 3, expected: 2, found: 3 }) => {}
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn overlong_line_is_rejected() {
        let opts = CsvOptions { has_header: false, max_line_bytes: 16, ..CsvOptions::default() };
        let text = format!("{}\n", "9".repeat(40));
        let mut rdr = CsvReader::new(BufReader::new(text.as_bytes()), opts).unwrap();
        assert!(matches!(rdr.next(), Some(Err(CsvError::LineTooLong { line: 1, limit: 16 }))));
    }

    #[test]
    fn malformed_quotes_are_rejected() {
        for text in ["a\n\"open\n", "a\n\"x\"y\n", "a\nx\"y\n"] {
            let mut rdr = reader(text);
            assert!(matches!(rdr.next(), Some(Err(CsvError::Malformed { .. }))), "{text:?}");
        }
    }

    #[test]
    fn invalid_options_rejected() {
        let bad = CsvOptions { delimiter: b'"', ..CsvOptions::default() };
        assert!(CsvReader::new(BufReader::new(&b""[..]), bad).is_err());
        let zero = CsvOptions { max_line_bytes: 0, ..CsvOptions::default() };
        assert!(CsvReader::new(BufReader::new(&b""[..]), zero).is_err());
    }

    #[test]
    fn unknown_column_and_index() {
        let rdr = reader("a,b\n1,2\n");
        assert!(matches!(rdr.columns(&["zz"]), Err(CsvError::UnknownColumn(_))));
        let mut rdr = reader("a\n1\n");
        let rec = rdr.next().unwrap().unwrap();
        assert!(matches!(rec.f64(5), Err(CsvError::NoSuchField { index: 5, .. })));
    }

    /// Feed random and structured-random bytes through the reader. The only
    /// requirement is that nothing panics and iteration terminates.
    #[test]
    fn fuzz_random_bytes_never_panic() {
        let alphabet: &[u8] = b"0123456789.,-+eE\"\n\r\t ab\xff\x00;";
        let mut rng = SplitMix64::new(0x00C0_FFEE);
        for round in 0..2000 {
            let len = rng.below(300) as usize;
            let data: Vec<u8> = (0..len)
                .map(|_| {
                    if round % 2 == 0 {
                        rng.next_u64() as u8
                    } else {
                        alphabet[rng.below(alphabet.len() as u64) as usize]
                    }
                })
                .collect();
            let opts = CsvOptions {
                delimiter: if round % 3 == 0 { b';' } else { b',' },
                has_header: round % 4 != 0,
                max_line_bytes: 64,
            };
            let Ok(rdr) = CsvReader::new(BufReader::new(&data[..]), opts) else {
                continue;
            };
            let mut count = 0usize;
            for rec in rdr {
                count += 1;
                assert!(count <= data.len() + 1, "reader did not terminate");
                if let Ok(rec) = rec {
                    for i in 0..=rec.len() {
                        let _ = rec.f64(i);
                        let _ = rec.u64(i);
                        let _ = rec.i64(i);
                    }
                }
            }
        }
    }
}
