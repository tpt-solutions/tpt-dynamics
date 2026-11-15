//! Shared CSV input helpers for `analyze` and `histogram`.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use tpt_dynamics_data::csv::{CsvError, CsvOptions, CsvReader};

use crate::error::{input, usage, CliError};

/// Converts a single-character delimiter option into a byte.
pub(crate) fn delimiter_byte(c: char) -> Result<u8, CliError> {
    match u8::try_from(c) {
        Ok(b) if b.is_ascii() && !matches!(b, b'"' | b'\n' | b'\r') => Ok(b),
        _ => Err(usage(format!(
            "invalid delimiter {c:?}: use one ASCII character other than quote or newline"
        ))),
    }
}

/// Opens a CSV file for streaming.
pub(crate) fn open_csv(
    path: &Path,
    delimiter: u8,
    has_header: bool,
) -> Result<CsvReader<BufReader<File>>, CliError> {
    let file =
        File::open(path).map_err(|e| input(format!("cannot open '{}': {e}", path.display())))?;
    let opts = CsvOptions { delimiter, has_header, ..CsvOptions::default() };
    CsvReader::new(BufReader::new(file), opts).map_err(|e| csv_err(&e))
}

/// Maps a CSV parse or read error to an input error.
pub(crate) fn csv_err(e: &CsvError) -> CliError {
    input(e.to_string())
}
