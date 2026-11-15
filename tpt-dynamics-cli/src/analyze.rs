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

    let (loaded, total_rows) = load(args, delimiter)?;
    let frame = apply_filters(loaded, &filters)?;
    print_summary(&frame, total_rows, quiet);

    if !hists.is_empty() {
        fs::create_dir_all(&args.out_dir)
            .map_err(|e| runtime(format!("cannot create '{}': {e}", args.out_dir.display())))?;
    }
    let ext = if args.format == Format::Bin { "bin" } else { "json" };
    for (column, bins, lo, hi) in &hists {
        let mut hist = Histogram1D::uniform(*bins, *lo, *hi)
            .map_err(|e| usage(format!("invalid --hist for '{column}': {e}")))?;
        for x in column_f64(&frame, column)? {
            hist.fill(x);
        }
        let path = args.out_dir.join(format!("{}.{ext}", sanitise(column)));
        write_histogram(&hist, &path, args.format)?;
        log::info!("histogram of '{column}' written to {}", path.display());
    }
    Ok(())
}

/// Reads the whole CSV into a typed frame and returns it with the number of
/// rows read. Every column must be numeric: a column whose fields all parse as
/// integers becomes `i64`, otherwise `f64`.
fn load(args: &AnalyzeArgs, delimiter: u8) -> Result<(DataFrame, usize), CliError> {
    let path = args.input.as_deref().ok_or_else(|| usage("analyze needs INPUT.csv"))?;
    let reader = open_csv(path, delimiter, !args.no_header)?;
    let header: Option<Vec<String>> = reader.header().map(<[String]>::to_vec);
    let records: Vec<CsvRecord> = reader.collect::<Result<_, _>>().map_err(|e| csv_err(&e))?;
    let ncols = match (&header, records.first()) {
        (Some(h), _) => h.len(),
        (None, Some(r)) => r.len(),
        (None, None) => 0,
    };
    let names: Vec<String> =
        header.unwrap_or_else(|| (0..ncols).map(|i| format!("c{i}")).collect());
    let mut columns = Vec::with_capacity(ncols);
    for (i, name) in names.into_iter().enumerate() {
        let ints: Result<Vec<i64>, _> = records.iter().map(|r| r.i64(i)).collect();
        let column = if let Ok(v) = ints {
            Column::I64(v)
        } else {
            let floats: Result<Vec<f64>, _> = records.iter().map(|r| r.f64(i)).collect();
            let v = floats.map_err(|e| input(format!("column '{name}' is not numeric: {e}")))?;
            Column::F64(v)
        };
        columns.push((name, column));
    }
    let n = records.len();
    let frame =
        DataFrame::with_rows(n, columns).map_err(|e| input(format!("cannot build table: {e}")))?;
    Ok((frame, n))
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
fn print_summary(frame: &DataFrame, total: usize, quiet: bool) {
    if quiet {
        return;
    }
    println!("rows: {} (of {total} read)", frame.len());
    println!(
        "{:<20} {:>5} {:>14} {:>14} {:>14} {:>14}",
        "column", "type", "min", "max", "mean", "stddev"
    );
    for name in frame.names() {
        let Some(c) = frame.column(name) else { continue };
        let values = match c {
            Column::F64(v) => v.clone(),
            Column::I64(v) => v.iter().map(|&x| x as f64).collect(),
            Column::Bool(_) => continue,
        };
        let (min, max, mean, sd) = stats(&values);
        println!("{:<20} {:>5} {:>14} {:>14} {:>14} {:>14}", name, c.dtype(), min, max, mean, sd);
    }
}

/// Min, max, mean and population standard deviation, formatted; `-` when empty.
fn stats(values: &[f64]) -> (String, String, String, String) {
    if values.is_empty() {
        return ("-".into(), "-".into(), "-".into(), "-".into());
    }
    let n = values.len() as f64;
    let min = values.iter().copied().fold(f64::INFINITY, f64::min);
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    (
        format!("{min:.6e}"),
        format!("{max:.6e}"),
        format!("{mean:.6e}"),
        format!("{:.6e}", var.sqrt()),
    )
}
