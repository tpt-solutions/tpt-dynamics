//! `tptdyn histogram`: streams one CSV column into a uniform 1D histogram.

use std::fs;
use std::path::Path;

use tpt_dynamics_data::hist1d::Histogram1D;

use crate::cli::{Format, HistogramArgs};
use crate::csvin::{csv_err, delimiter_byte, open_csv};
use crate::error::{input, runtime, usage, CliError};

/// Runs the `histogram` subcommand.
pub(crate) fn run(args: &HistogramArgs, quiet: bool) -> Result<(), CliError> {
    let (lo, hi) = parse_range(&args.range)?;
    let mut hist = Histogram1D::uniform(args.bins, lo, hi)
        .map_err(|e| usage(format!("invalid histogram: {e}")))?;
    let delimiter = delimiter_byte(args.delimiter)?;
    let mut reader = open_csv(&args.input, delimiter, !args.no_header)?;
    let idx = resolve_column(&reader, &args.column)?;

    let mut rows = 0_u64;
    for record in &mut reader {
        let record = record.map_err(|e| csv_err(&e))?;
        let x = record.f64(idx).map_err(|e| csv_err(&e))?;
        hist.fill(x);
        rows += 1;
    }
    let format = args.format.unwrap_or_else(|| infer_format(&args.out));
    write_histogram(&hist, &args.out, format)?;
    log::info!("histogram of '{}' from {rows} rows written to {}", args.column, args.out.display());
    if !quiet {
        eprintln!(
            "filled {} entries ({} underflow, {} overflow) -> {}",
            hist.entries(),
            hist.underflow(),
            hist.overflow(),
            args.out.display()
        );
    }
    Ok(())
}

/// Parses `LO,HI`.
pub(crate) fn parse_range(text: &str) -> Result<(f64, f64), CliError> {
    let parts: Vec<&str> = text.split(',').map(str::trim).collect();
    let [lo, hi] = parts.as_slice() else {
        return Err(usage(format!("--range expects LO,HI, got '{text}'")));
    };
    let lo: f64 = lo.parse().map_err(|_| usage(format!("invalid range bound '{lo}'")))?;
    let hi: f64 = hi.parse().map_err(|_| usage(format!("invalid range bound '{hi}'")))?;
    Ok((lo, hi))
}

/// Resolves a column given by name (with a header) or by zero-based index.
pub(crate) fn resolve_column<R: std::io::BufRead>(
    reader: &tpt_dynamics_data::csv::CsvReader<R>,
    spec: &str,
) -> Result<usize, CliError> {
    if let Some(i) = reader.column_index(spec) {
        return Ok(i);
    }
    if reader.header().is_none() {
        if let Ok(i) = spec.parse::<usize>() {
            return Ok(i);
        }
    }
    Err(input(format!("no column named '{spec}' in the input")))
}

/// Picks the format from the file extension when `--format` is not given.
pub(crate) fn infer_format(path: &Path) -> Format {
    match path.extension().and_then(|e| e.to_str()) {
        Some(e) if e.eq_ignore_ascii_case("bin") => Format::Bin,
        _ => Format::Json,
    }
}

/// Writes a histogram in the requested format, creating parent directories.
pub(crate) fn write_histogram(
    hist: &Histogram1D,
    path: &Path,
    format: Format,
) -> Result<(), CliError> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)
            .map_err(|e| runtime(format!("cannot create '{}': {e}", parent.display())))?;
    }
    let bytes = match format {
        Format::Json => hist.to_json().into_bytes(),
        Format::Bin => hist.to_bytes(),
    };
    fs::write(path, bytes).map_err(|e| runtime(format!("cannot write '{}': {e}", path.display())))
}
