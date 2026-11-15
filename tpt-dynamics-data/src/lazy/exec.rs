//! Chunked executor.
//!
//! A scan is read in ranges of `chunk_size` rows. Up to `window` chunks are held
//! in memory at once; within a window the row-wise operators run on each chunk
//! (in parallel with rayon when `parallel` is set). Results are merged in chunk
//! order, so the output does not depend on thread scheduling.
//!
//! Plans that are not a scan followed by row-wise operators (for example a
//! filter over an aggregate) materialise their input first and then apply the
//! operator eagerly.

use std::borrow::Cow;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::ops::Range;

use rayon::prelude::*;

use super::expr::eval;
use super::frame::{Column, DataFrame, LazyError, Schema, Source};
use super::plan::{validate, AggExpr, AggKind, Plan};

/// Controls how a plan is executed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CollectOptions {
    /// Rows per chunk. Peak memory is about `window * chunk_size` rows of the
    /// columns the plan reads.
    pub chunk_size: usize,
    /// Maximum number of chunks held in memory at once.
    pub window: usize,
    /// Process the chunks of each window on the rayon thread pool.
    pub parallel: bool,
}

impl Default for CollectOptions {
    fn default() -> Self {
        Self { chunk_size: 65_536, window: 16, parallel: true }
    }
}

pub(crate) fn execute(plan: &Plan, opts: &CollectOptions) -> Result<DataFrame, LazyError> {
    if opts.chunk_size == 0 || opts.window == 0 {
        return Err(LazyError::Invalid("chunk_size and window must both be at least 1".into()));
    }
    exec_node(plan, opts)
}

/// A scan plus the row-wise operators above it, bottom-up.
struct Stream<'a> {
    source: &'a dyn Source,
    projection: Vec<String>,
    ops: Vec<&'a Plan>,
}

fn streamable(plan: &Plan) -> Option<Stream<'_>> {
    let mut ops = Vec::new();
    let mut node = plan;
    loop {
        match node {
            Plan::Scan { source, projection } => {
                ops.reverse();
                let projection = match projection {
                    Some(p) => p.clone(),
                    None => source.schema().into_iter().map(|f| f.name).collect(),
                };
                return Some(Stream { source: source.as_ref(), projection, ops });
            }
            Plan::Filter { input, .. }
            | Plan::Project { input, .. }
            | Plan::WithColumn { input, .. } => {
                ops.push(node);
                node = input.as_ref();
            }
            Plan::Aggregate { .. } => return None,
        }
    }
}

impl Stream<'_> {
    /// Reads every chunk, applies the operators, runs `f` on each result and
    /// hands the outputs to `sink` in chunk order.
    fn for_each_chunk<T: Send>(
        &self,
        opts: &CollectOptions,
        f: impl Fn(DataFrame) -> Result<T, LazyError> + Sync,
        mut sink: impl FnMut(T) -> Result<(), LazyError>,
    ) -> Result<(), LazyError> {
        let rows = self.source.num_rows();
        let ranges: Vec<Range<usize>> = (0..rows)
            .step_by(opts.chunk_size)
            .map(|s| s..s.saturating_add(opts.chunk_size).min(rows))
            .collect();
        let process = |df: DataFrame| -> Result<T, LazyError> {
            let df = apply_ops(df, &self.ops)?;
            f(df)
        };
        for window in ranges.chunks(opts.window) {
            let batch = window
                .iter()
                .map(|r| self.source.read(r.clone(), &self.projection))
                .collect::<Result<Vec<_>, _>>()?;
            let results: Vec<Result<T, LazyError>> = if opts.parallel {
                batch.into_par_iter().map(process).collect()
            } else {
                batch.into_iter().map(process).collect()
            };
            for r in results {
                sink(r?)?;
            }
        }
        Ok(())
    }
}

fn exec_node(plan: &Plan, opts: &CollectOptions) -> Result<DataFrame, LazyError> {
    if let Plan::Aggregate { input, keys, aggs } = plan {
        let schema = validate(plan)?;
        let mut groups = Partial::new();
        if let Some(stream) = streamable(input) {
            stream.for_each_chunk(
                opts,
                |df| partial_agg(&df, keys, aggs),
                |p| {
                    merge(&mut groups, p);
                    Ok(())
                },
            )?;
        } else {
            let df = exec_node(input, opts)?;
            merge(&mut groups, partial_agg(&df, keys, aggs)?);
        }
        return finalize(groups, &schema);
    }

    if let Some(stream) = streamable(plan) {
        let schema = validate(plan)?;
        let mut parts = Vec::new();
        stream.for_each_chunk(opts, Ok, |df| {
            parts.push(df);
            Ok(())
        })?;
        return concat(&schema, &parts);
    }

    // A row-wise operator above a materialised input (for example a filter
    // over an aggregate): run the input, then apply the operator once.
    let input =
        plan.input().ok_or_else(|| LazyError::Invalid("internal: scan without stream".into()))?;
    let df = exec_node(input, opts)?;
    apply_op(plan, df)
}

fn apply_ops(df: DataFrame, ops: &[&Plan]) -> Result<DataFrame, LazyError> {
    let mut df = df;
    for op in ops {
        df = apply_op(op, df)?;
    }
    Ok(df)
}

fn apply_op(op: &Plan, df: DataFrame) -> Result<DataFrame, LazyError> {
    match op {
        Plan::Filter { predicate, .. } => {
            let mask_col = eval(predicate, &df)?;
            let Column::Bool(mask) = mask_col.as_ref() else {
                return Err(LazyError::Type("filter predicate must be boolean".into()));
            };
            df.filter_rows(mask)
        }
        Plan::Project { exprs, .. } => {
            let mut cols = Vec::with_capacity(exprs.len());
            for e in exprs {
                let name = e
                    .output_name()
                    .ok_or_else(|| LazyError::Invalid(format!("unnamed expression `{e}`")))?;
                let col = eval(e, &df)?.into_owned();
                cols.push((name.to_string(), col));
            }
            df.derive(df.len(), cols)
        }
        Plan::WithColumn { expr, .. } => {
            let name = expr
                .output_name()
                .ok_or_else(|| LazyError::Invalid(format!("unnamed expression `{expr}`")))?
                .to_string();
            let col = eval(expr, &df)?.into_owned();
            df.set_column(name, col)
        }
        _ => Err(LazyError::Invalid("internal: not a row-wise operator".into())),
    }
}

/// Group key element. Ordering gives a deterministic output order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum KeyVal {
    I64(i64),
    Bool(bool),
}

fn key_at(col: &Column, row: usize) -> KeyVal {
    match col {
        Column::I64(v) => KeyVal::I64(v[row]),
        Column::Bool(v) => KeyVal::Bool(v[row]),
        // Rejected by validation; kept total for safety.
        Column::F64(v) => KeyVal::I64(v[row] as i64),
    }
}

/// Running state of one aggregate within one group. Partial states merge
/// associatively.
#[derive(Clone, Copy, Debug)]
enum Acc {
    Count(i64),
    SumI(i64),
    SumF(f64),
    Mean { sum: f64, n: i64 },
    MinF(f64),
    MaxF(f64),
    MinI(i64),
    MaxI(i64),
}

fn f64_at(col: Option<&Column>, row: usize) -> f64 {
    match col {
        Some(Column::F64(v)) => v[row],
        Some(Column::I64(v)) => v[row] as f64,
        Some(Column::Bool(v)) => f64::from(u8::from(v[row])),
        None => 0.0,
    }
}

fn i64_at(col: Option<&Column>, row: usize) -> i64 {
    match col {
        Some(Column::I64(v)) => v[row],
        Some(Column::F64(v)) => v[row] as i64,
        Some(Column::Bool(v)) => i64::from(v[row]),
        None => 0,
    }
}

impl Acc {
    fn new(kind: AggKind, input: Option<&Column>) -> Self {
        let int = matches!(input, Some(Column::I64(_)));
        match kind {
            AggKind::Count => Self::Count(0),
            AggKind::Sum if int => Self::SumI(0),
            AggKind::Sum => Self::SumF(0.0),
            AggKind::Mean => Self::Mean { sum: 0.0, n: 0 },
            AggKind::Min if int => Self::MinI(i64::MAX),
            AggKind::Min => Self::MinF(f64::INFINITY),
            AggKind::Max if int => Self::MaxI(i64::MIN),
            AggKind::Max => Self::MaxF(f64::NEG_INFINITY),
        }
    }

    fn update(&mut self, col: Option<&Column>, row: usize) {
        match self {
            Self::Count(n) => *n += 1,
            Self::SumI(s) => *s = s.wrapping_add(i64_at(col, row)),
            Self::SumF(s) => *s += f64_at(col, row),
            Self::Mean { sum, n } => {
                *sum += f64_at(col, row);
                *n += 1;
            }
            Self::MinF(m) => *m = m.min(f64_at(col, row)),
            Self::MaxF(m) => *m = m.max(f64_at(col, row)),
            Self::MinI(m) => *m = (*m).min(i64_at(col, row)),
            Self::MaxI(m) => *m = (*m).max(i64_at(col, row)),
        }
    }

    fn merge(&mut self, other: Acc) {
        match (self, other) {
            (Self::Count(a), Self::Count(b)) => *a += b,
            (Self::SumI(a), Self::SumI(b)) => *a = a.wrapping_add(b),
            (Self::SumF(a), Self::SumF(b)) => *a += b,
            (Self::Mean { sum: a, n: an }, Self::Mean { sum: b, n: bn }) => {
                *a += b;
                *an += bn;
            }
            (Self::MinF(a), Self::MinF(b)) => *a = a.min(b),
            (Self::MaxF(a), Self::MaxF(b)) => *a = a.max(b),
            (Self::MinI(a), Self::MinI(b)) => *a = (*a).min(b),
            (Self::MaxI(a), Self::MaxI(b)) => *a = (*a).max(b),
            // Partials are built from the same aggregate list, so kinds match.
            _ => {}
        }
    }

    fn as_i64(&self) -> i64 {
        match *self {
            Self::Count(n) | Self::SumI(n) | Self::MinI(n) | Self::MaxI(n) => n,
            Self::SumF(x) | Self::MinF(x) | Self::MaxF(x) => x as i64,
            Self::Mean { sum, n } => (sum / n as f64) as i64,
        }
    }

    fn as_f64(&self) -> f64 {
        match *self {
            Self::Count(n) | Self::SumI(n) | Self::MinI(n) | Self::MaxI(n) => n as f64,
            Self::SumF(x) | Self::MinF(x) | Self::MaxF(x) => x,
            Self::Mean { sum, n } => sum / n as f64,
        }
    }
}

type Partial = HashMap<Vec<KeyVal>, Vec<Acc>>;

fn partial_agg(df: &DataFrame, keys: &[String], aggs: &[AggExpr]) -> Result<Partial, LazyError> {
    let key_cols = keys
        .iter()
        .map(|k| df.column(k).ok_or_else(|| LazyError::UnknownColumn(k.clone())))
        .collect::<Result<Vec<_>, _>>()?;
    let inputs: Vec<Option<Cow<'_, Column>>> = aggs
        .iter()
        .map(|a| a.input().map(|e| eval(e, df)).transpose())
        .collect::<Result<_, _>>()?;
    let mut map = Partial::new();
    for row in 0..df.len() {
        let key: Vec<KeyVal> = key_cols.iter().map(|c| key_at(c, row)).collect();
        let accs = map.entry(key).or_insert_with(|| {
            aggs.iter().zip(&inputs).map(|(a, inp)| Acc::new(a.kind(), inp.as_deref())).collect()
        });
        for (acc, inp) in accs.iter_mut().zip(&inputs) {
            acc.update(inp.as_deref(), row);
        }
    }
    Ok(map)
}

fn merge(into: &mut Partial, part: Partial) {
    for (key, accs) in part {
        match into.entry(key) {
            Entry::Occupied(mut e) => {
                for (a, b) in e.get_mut().iter_mut().zip(accs) {
                    a.merge(b);
                }
            }
            Entry::Vacant(v) => {
                v.insert(accs);
            }
        }
    }
}

fn push_key(col: &mut Column, kv: KeyVal) -> Result<(), LazyError> {
    match (col, kv) {
        (Column::I64(v), KeyVal::I64(x)) => v.push(x),
        (Column::Bool(v), KeyVal::Bool(x)) => v.push(x),
        _ => return Err(LazyError::Type("group key type mismatch".into())),
    }
    Ok(())
}

fn push_acc(col: &mut Column, acc: &Acc) -> Result<(), LazyError> {
    match col {
        Column::I64(v) => v.push(acc.as_i64()),
        Column::F64(v) => v.push(acc.as_f64()),
        Column::Bool(_) => return Err(LazyError::Type("aggregate cannot produce bool".into())),
    }
    Ok(())
}

fn finalize(groups: Partial, schema: &Schema) -> Result<DataFrame, LazyError> {
    let mut entries: Vec<(Vec<KeyVal>, Vec<Acc>)> = groups.into_iter().collect();
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    let mut cols: Vec<Column> = schema.iter().map(|f| Column::empty(f.dtype)).collect();
    for (key, accs) in &entries {
        for (i, kv) in key.iter().enumerate() {
            push_key(&mut cols[i], *kv)?;
        }
        for (j, acc) in accs.iter().enumerate() {
            push_acc(&mut cols[key.len() + j], acc)?;
        }
    }
    DataFrame::with_rows(
        entries.len(),
        schema.iter().map(|f| f.name.clone()).zip(cols).collect::<Vec<_>>(),
    )
}

fn concat(schema: &Schema, parts: &[DataFrame]) -> Result<DataFrame, LazyError> {
    let mut cols: Vec<Column> = schema.iter().map(|f| Column::empty(f.dtype)).collect();
    let mut rows = 0;
    for part in parts {
        rows += part.len();
        for (f, col) in schema.iter().zip(cols.iter_mut()) {
            let src =
                part.column(&f.name).ok_or_else(|| LazyError::UnknownColumn(f.name.clone()))?;
            col.append(src)?;
        }
    }
    DataFrame::with_rows(rows, schema.iter().map(|f| f.name.clone()).zip(cols).collect::<Vec<_>>())
}
