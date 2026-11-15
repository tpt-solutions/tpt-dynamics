//! Bounded line reading for the fixed-width structure readers.
//!
//! [`std::io::BufRead::lines`] allocates a line of any length, so a file with
//! no newline can exhaust memory before parsing starts (security finding S-1).
//! [`BoundedLines`] reads at most [`ReadLimits::max_line_bytes`] bytes per line
//! and [`ReadLimits::max_atoms`] atoms per file, and reports either limit as a
//! [`SimError::Parse`] that names the line.

use std::io::BufRead;

use crate::error::SimError;

/// Largest accepted line length in bytes (1 MiB). Real PDB and GRO lines are
/// at most a few hundred bytes.
pub const MAX_LINE_BYTES: usize = 1 << 20;

/// Largest accepted atom count per structure file (10 million).
pub const MAX_ATOMS: usize = 10_000_000;

/// Size limits applied while reading a structure file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Maximum bytes in one line, excluding the line terminator.
    pub max_line_bytes: usize,
    /// Maximum number of atoms in one file.
    pub max_atoms: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self { max_line_bytes: MAX_LINE_BYTES, max_atoms: MAX_ATOMS }
    }
}

/// Iterator over the lines of a reader, each at most `max_line_bytes` long.
///
/// Lines end at `\n`; a trailing `\r` is removed, as [`BufRead::lines`] does.
pub struct BoundedLines<R> {
    inner: R,
    format: &'static str,
    max_line_bytes: usize,
    line_no: usize,
}

impl<R: BufRead> BoundedLines<R> {
    /// Wraps `inner`. `format` names the file format in error messages.
    pub fn new(inner: R, format: &'static str, max_line_bytes: usize) -> Self {
        Self { inner, format, max_line_bytes, line_no: 0 }
    }

    fn read_raw(&mut self) -> Result<Option<Vec<u8>>, SimError> {
        let mut buf: Vec<u8> = Vec::new();
        let line_no = self.line_no + 1;
        let mut saw_any = false;
        loop {
            let available = self.inner.fill_buf()?;
            if available.is_empty() {
                return Ok(saw_any.then_some(buf));
            }
            saw_any = true;
            let (take, done) = match available.iter().position(|&b| b == b'\n') {
                Some(pos) => (pos + 1, true),
                None => (available.len(), false),
            };
            let content = if done { &available[..take - 1] } else { &available[..take] };
            buf.extend_from_slice(content);
            self.inner.consume(take);
            if buf.len() > self.max_line_bytes {
                return Err(SimError::parse(
                    self.format,
                    line_no,
                    format!("line longer than {} bytes", self.max_line_bytes),
                ));
            }
            if done {
                return Ok(Some(buf));
            }
        }
    }
}

impl<R: BufRead> Iterator for BoundedLines<R> {
    type Item = Result<String, SimError>;

    fn next(&mut self) -> Option<Self::Item> {
        let raw = match self.read_raw() {
            Ok(Some(raw)) => raw,
            Ok(None) => return None,
            Err(e) => return Some(Err(e)),
        };
        self.line_no += 1;
        let mut raw = raw;
        if raw.last() == Some(&b'\r') {
            raw.pop();
        }
        Some(
            String::from_utf8(raw)
                .map_err(|_| SimError::parse(self.format, self.line_no, "line is not valid UTF-8")),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(text: &[u8], max: usize) -> Vec<Result<String, SimError>> {
        BoundedLines::new(text, "TST", max).collect()
    }

    #[test]
    fn splits_lines_like_bufread_lines() {
        let got: Vec<String> =
            collect(b"a\r\nbb\n\nccc", 100).into_iter().map(Result::unwrap).collect();
        assert_eq!(got, ["a", "bb", "", "ccc"]);
        assert!(collect(b"", 10).is_empty());
        let trailing: Vec<String> = collect(b"x\n", 10).into_iter().map(Result::unwrap).collect();
        assert_eq!(trailing, ["x"]);
    }

    #[test]
    fn rejects_overlong_line_without_a_newline() {
        let big = vec![b'a'; 5000];
        let res = collect(&big, 1024);
        assert_eq!(res.len(), 1);
        let err = res[0].as_ref().unwrap_err().to_string();
        assert!(err.contains("line 1") && err.contains("1024"), "{err}");
    }

    #[test]
    fn rejects_invalid_utf8() {
        let res = collect(&[0xFF, 0xFE, b'\n'], 10);
        assert!(res[0].is_err());
    }
}
