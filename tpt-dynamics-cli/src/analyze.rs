//! `tptdyn analyze`: loads a CSV event file, applies lazy-frame filters, prints
//! per-column summaries and writes uniform histograms.

use std::fs;

use tpt_dynamics_data::csv::CsvRecord;
use tpt_dynamics_data::hist1d::Histogram1D;
use tpt_dynamics_data::lazy::{col, Column, DType, DataFrame, Expr, LazyFrame};

use crate::cli::{AnalyzeArgs, Format};
use crate::csvin::{csv_err, delimiter_byte, open_csv};
use crate::error::{input, runtime, usage, CliError};
use crate::histogram::{parse_range, write_histogram};

/// Comparison operators accepted by `--filter`.
#[derive(Clone, Copy, Debug)]
enum Op {
    Gt,
    Ge,
    Lt,
    Le,
    Eq,
    Ne,
}

impl Op {
    fn parse(text: &str) -> Result<Self, CliError> {
        Ok(match text {
            ">" => Self::Gt,
            ">=" => Self::Ge,
            "<" => Self::Lt,
            "<=" => Self::Le,
            "==" => Self::Eq,
            "!=" => Self::Ne,
            other => return Err(usage(format!("unknown operator '{other}' in --filter"))),
        })
    }
}

/// A parsed `--filter`: column, operator and the raw value text.
type FilterSpec = (String, Op, String);

/// A parsed `--hist`: column, bins, lo, hi.
type HistSpec = (String, usize, f64, f64);

/// Runs the `analyze` subcommand.
pub(crate) fn run(args: &AnalyzeArgs, quiet: bool) -> Result<(), CliError> {
    if args.trajectory.is_some() {
        return crate::traj_analysis::run(args, quiet);
    }
    // Parse every argument before touching the input, so usage errors come first.
    let delimiter = delimiter_byte(args.delimiter)?;
    let filters = args.filters.iter().map(|f| parse_filter(f)).collect::<Result<Vec<_>, _>>()?;
    let hists = args.hists.iter().map(|h| parse_hist_spec(h)).collect::<Result<Vec<_>, _>>()?;

    let path = args.input.as_deref().ok_or_else(|| usage("analyze needs INPUT.csv"))?;
    let mut reader = open_csv(path, delimiter, !args.no_header)?;
    // Without a header, column names are `c0, c1, ...`, taken from the first chunk.
    let mut names: Vec<String> = reader.header().map(<[String]>::to_vec).unwrap_or_default();
    let mut hist_objs = hists
        .iter()
        .map(|(column, bins, lo, hi)| {
            let hist = Histogram1D::uniform(*bins, *lo, *hi)
                .map_err(|e| usage(format!("invalid --hist for '{column}': {e}")))?;
            Ok((column.clone(), hist))
        })
        .collect::<Result<Vec<_>, CliError>>()?;

    // Stream the file in chunks: memory is bounded by CHUNK_ROWS, not by the file size.
    // Each chunk is filtered, then folded into the summary and the histograms.
    let mut summary = Summary::default();
    let mut total_rows = 0usize;
    loop {
        let records: Vec<CsvRecord> =
            reader.by_ref().take(CHUNK_ROWS).collect::<Result<_, _>>().map_err(|e| csv_err(&e))?;
        if records.is_empty() {
            break;
        }
        if names.is_empty() {
            names = (0..records[0].len()).map(|i| format!("c{i}")).collect();
        }
        total_rows += records.len();
        let frame = apply_filters(frame_from_records(&names, &records)?, &filters)?;
        summary.add(&frame);
        for (column, hist) in &mut hist_objs {
            for x in column_f64(&frame, column)? {
                hist.fill(x);
            }
        }
    }
    print_summary(&summary, total_rows, quiet);

    if let Some((column, _)) = hist_objs.iter().find(|(c, _)| !names.contains(c)) {
        return Err(input(format!("no column named '{column}' in the input")));
    }
    if !hist_objs.is_empty() {
        fs::create_dir_all(&args.out_dir)
            .map_err(|e| runtime(format!("cannot create '{}': {e}", args.out_dir.display())))?;
    }
    let ext = if args.format == Format::Bin { "bin" } else { "json" };
    for (column, hist) in &hist_objs {
        let path = args.out_dir.join(format!("{}.{ext}", sanitise(column)));
        write_histogram(hist, &path, args.format)?;
        log::info!("histogram of '{column}' written to {}", path.display());
    }
    Ok(())
}

/// Rows per chunk when streaming a CSV file.
const CHUNK_ROWS: usize = 1 << 20;

/// Per-column statistics over every filtered row seen so far.
#[derive(Default)]
struct Summary {
    rows: usize,
    columns: Vec<(String, DType, RunningStats)>,
}

impl Summary {
    /// Folds one filtered chunk into the summary. A column that is `i64` in some
    /// chunks and `f64` in others is reported as `f64`.
    fn add(&mut self, frame: &DataFrame) {
        self.rows += frame.len();
        for name in frame.names() {
            let Some(c) = frame.column(name) else { continue };
            let idx = if let Some(i) = self.columns.iter().position(|(n, _, _)| n == name) {
                i
            } else {
                self.columns.push((name.clone(), c.dtype(), RunningStats::new()));
                self.columns.len() - 1
            };
            let entry = &mut self.columns[idx];
            if entry.1 != c.dtype() {
                entry.1 = DType::F64;
            }
            let stats = &mut entry.2;
            match c {
                Column::F64(v) => stats.push_all(v.iter().copied()),
                Column::I64(v) => stats.push_all(v.iter().map(|&x| x as f64)),
                Column::Bool(_) => {}
            }
        }
    }
}

/// Welford running mean and variance, with min and max.
#[derive(Clone, Copy)]
struct RunningStats {
    n: u64,
    mean: f64,
    m2: f64,
    min: f64,
    max: f64,
}

impl RunningStats {
    fn new() -> Self {
        Self { n: 0, mean: 0.0, m2: 0.0, min: f64::INFINITY, max: f64::NEG_INFINITY }
    }

    fn push_all(&mut self, values: impl Iterator<Item = f64>) {
        for x in values {
            self.n += 1;
            let delta = x - self.mean;
            self.mean += delta / self.n as f64;
            self.m2 += delta * (x - self.mean);
            self.min = self.min.min(x);
            self.max = self.max.max(x);
        }
    }
}

/// Builds a typed frame from one chunk of records. Every column must be numeric:
/// a column whose fields all parse as integers in this chunk becomes `i64`,
/// otherwise `f64`.
fn frame_from_records(names: &[String], records: &[CsvRecord]) -> Result<DataFrame, CliError> {
    let mut columns = Vec::with_capacity(names.len());
    for (i, name) in names.iter().enumerate() {
        let ints: Result<Vec<i64>, _> = records.iter().map(|r| r.i64(i)).collect();
        let column = if let Ok(v) = ints {
            Column::I64(v)
        } else {
            let floats: Result<Vec<f64>, _> = records.iter().map(|r| r.f64(i)).collect();
            let v = floats.map_err(|e| input(format!("column '{name}' is not numeric: {e}")))?;
            Column::F64(v)
        };
        columns.push((name.clone(), column));
    }
    DataFrame::with_rows(records.len(), columns)
        .map_err(|e| input(format!("cannot build table: {e}")))
}

/// Keeps the rows that satisfy every filter (logical AND), using the lazy
/// engine. Returns the input unchanged when there are no filters.
fn apply_filters(frame: DataFrame, filters: &[FilterSpec]) -> Result<DataFrame, CliError> {
    if filters.is_empty() {
        return Ok(frame);
    }
    let exprs = filters
        .iter()
        .map(|(column, op, value)| filter_expr(&frame, column, *op, value))
        .collect::<Result<Vec<_>, _>>()?;
    let mut lf = LazyFrame::scan(frame);
    for e in exprs {
        lf = lf.filter(e);
    }
    lf.collect().map_err(|e| input(format!("filter failed: {e}")))
}

/// Splits `COLUMN OP VALUE` (whitespace separated) into its three parts.
fn parse_filter(text: &str) -> Result<FilterSpec, CliError> {
    let parts: Vec<&str> = text.split_whitespace().collect();
    match parts.as_slice() {
        [column, op, value] => Ok(((*column).to_string(), Op::parse(op)?, (*value).to_string())),
        _ => Err(usage(format!("--filter expects 'COLUMN OP VALUE', got '{text}'"))),
    }
}

/// Builds the lazy-frame predicate for one filter, checking that the column exists.
fn filter_expr(frame: &DataFrame, column: &str, op: Op, value: &str) -> Result<Expr, CliError> {
    let dtype = frame
        .column(column)
        .map(Column::dtype)
        .ok_or_else(|| input(format!("no column named '{column}' in the input")))?;
    let rhs: Expr = match dtype {
        DType::I64 => match value.parse::<i64>() {
            Ok(v) => Expr::from(v),
            Err(_) => Expr::from(parse_number(value)?),
        },
        DType::F64 => Expr::from(parse_number(value)?),
        DType::Bool => return Err(usage(format!("cannot filter boolean column '{column}'"))),
    };
    let c = col(column);
    Ok(match op {
        Op::Gt => c.gt(rhs),
        Op::Ge => c.ge(rhs),
        Op::Lt => c.lt(rhs),
        Op::Le => c.le(rhs),
        Op::Eq => c.eq_to(rhs),
        Op::Ne => c.ne_to(rhs),
    })
}

fn parse_number(text: &str) -> Result<f64, CliError> {
    text.parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| usage(format!("'{text}' is not a finite number")))
}

/// Parses `COLUMN,BINS,LO,HI`.
fn parse_hist_spec(text: &str) -> Result<HistSpec, CliError> {
    let parts: Vec<&str> = text.splitn(3, ',').collect();
    let [column, bins, range] = parts.as_slice() else {
        return Err(usage(format!("--hist expects COLUMN,BINS,LO,HI, got '{text}'")));
    };
    let bins: usize = bins
        .trim()
        .parse()
        .ok()
        .filter(|&n| n >= 1)
        .ok_or_else(|| usage(format!("invalid bin count in --hist '{text}'")))?;
    let (lo, hi) = parse_range(range)?;
    Ok(((*column).to_string(), bins, lo, hi))
}

/// Values of a numeric column as `f64`.
fn column_f64(frame: &DataFrame, name: &str) -> Result<Vec<f64>, CliError> {
    match frame.column(name) {
        Some(Column::F64(v)) => Ok(v.clone()),
        Some(Column::I64(v)) => Ok(v.iter().map(|&x| x as f64).collect()),
        Some(Column::Bool(_)) => Err(usage(format!("cannot histogram boolean column '{name}'"))),
        None => Err(input(format!("no column named '{name}' in the input"))),
    }
}

/// Replaces characters that are unsafe in file names.
fn sanitise(name: &str) -> String {
    name.chars()
        .map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
        .collect()
}

/// Prints the row counts and a per-column table to stdout.
fn print_summary(summary: &Summary, total: usize, quiet: bool) {
    if quiet {
        return;
    }
    println!("rows: {} (of {total} read)", summary.rows);
    println!(
        "{:<20} {:>5} {:>14} {:>14} {:>14} {:>14}",
        "column", "type", "min", "max", "mean", "stddev"
    );
    for (name, dtype, stats) in &summary.columns {
        let (min, max, mean, sd) = stats.formatted();
        println!("{name:<20} {dtype:>5} {min:>14} {max:>14} {mean:>14} {sd:>14}");
    }
}

impl RunningStats {
    /// Min, max, mean and population standard deviation, formatted; `-` when empty.
    fn formatted(&self) -> (String, String, String, String) {
        if self.n == 0 {
            return ("-".into(), "-".into(), "-".into(), "-".into());
        }
        (
            format!("{:.6e}", self.min),
            format!("{:.6e}", self.max),
            format!("{:.6e}", self.mean),
            format!("{:.6e}", (self.m2 / self.n as f64).sqrt()),
        )
    }
}
