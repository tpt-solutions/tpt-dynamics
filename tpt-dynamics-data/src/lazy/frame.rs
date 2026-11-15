//! Columnar in-memory frames, the [`Source`] trait and the [`LazyFrame`] builder.

use std::any::Any;
use std::fmt;
use std::ops::Range;
use std::sync::Arc;

use super::exec::{execute, CollectOptions};
use super::expr::Expr;
use super::optimize::optimize;
use super::plan::{format_plan, validate, AggExpr, Plan};

/// Errors raised while building, validating or executing a lazy plan.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LazyError {
    /// A referenced column does not exist.
    UnknownColumn(String),
    /// Two columns or outputs share a name.
    DuplicateColumn(String),
    /// An operation was applied to a value of the wrong type.
    Type(String),
    /// The plan or an option is invalid, or columns have mismatched lengths.
    Invalid(String),
}

impl fmt::Display for LazyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownColumn(c) => write!(f, "unknown column `{c}`"),
            Self::DuplicateColumn(c) => write!(f, "duplicate column `{c}`"),
            Self::Type(m) => write!(f, "type error: {m}"),
            Self::Invalid(m) => write!(f, "invalid: {m}"),
        }
    }
}

impl std::error::Error for LazyError {}

/// Element type of a column.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DType {
    /// 64-bit float.
    F64,
    /// 64-bit signed integer.
    I64,
    /// Boolean.
    Bool,
}

impl fmt::Display for DType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::F64 => "f64",
            Self::I64 => "i64",
            Self::Bool => "bool",
        })
    }
}

/// A typed column of values. Columns have no nulls.
#[derive(Clone, Debug, PartialEq)]
pub enum Column {
    /// Floating-point values.
    F64(Vec<f64>),
    /// Integer values.
    I64(Vec<i64>),
    /// Boolean values.
    Bool(Vec<bool>),
}

impl Column {
    /// Element type of the column.
    pub fn dtype(&self) -> DType {
        match self {
            Self::F64(_) => DType::F64,
            Self::I64(_) => DType::I64,
            Self::Bool(_) => DType::Bool,
        }
    }

    /// Number of rows.
    pub fn len(&self) -> usize {
        match self {
            Self::F64(v) => v.len(),
            Self::I64(v) => v.len(),
            Self::Bool(v) => v.len(),
        }
    }

    /// Whether the column has no rows.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// An empty column of the given element type.
    pub fn empty(dtype: DType) -> Self {
        match dtype {
            DType::F64 => Self::F64(Vec::new()),
            DType::I64 => Self::I64(Vec::new()),
            DType::Bool => Self::Bool(Vec::new()),
        }
    }

    pub(crate) fn slice(&self, range: Range<usize>) -> Self {
        match self {
            Self::F64(v) => Self::F64(v[range].to_vec()),
            Self::I64(v) => Self::I64(v[range].to_vec()),
            Self::Bool(v) => Self::Bool(v[range].to_vec()),
        }
    }

    pub(crate) fn filter(&self, mask: &[bool]) -> Self {
        fn keep<T: Copy>(v: &[T], mask: &[bool]) -> Vec<T> {
            v.iter().zip(mask).filter(|(_, &m)| m).map(|(&x, _)| x).collect()
        }
        match self {
            Self::F64(v) => Self::F64(keep(v, mask)),
            Self::I64(v) => Self::I64(keep(v, mask)),
            Self::Bool(v) => Self::Bool(keep(v, mask)),
        }
    }

    pub(crate) fn append(&mut self, other: &Self) -> Result<(), LazyError> {
        match (self, other) {
            (Self::F64(a), Self::F64(b)) => a.extend_from_slice(b),
            (Self::I64(a), Self::I64(b)) => a.extend_from_slice(b),
            (Self::Bool(a), Self::Bool(b)) => a.extend_from_slice(b),
            _ => return Err(LazyError::Type("appending columns of different types".into())),
        }
        Ok(())
    }
}

/// A named, typed field of a schema.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Field {
    /// Column name.
    pub name: String,
    /// Element type.
    pub dtype: DType,
}

/// Ordered list of fields.
pub type Schema = Vec<Field>;

/// An in-memory columnar table.
///
/// A frame may carry a *lease*: an arbitrary value kept alive as long as the
/// frame (or any frame derived from it by the executor). Sources use this to
/// tie a chunk to a buffer or a resource-accounting token.
#[derive(Clone)]
pub struct DataFrame {
    names: Vec<String>,
    columns: Vec<Column>,
    rows: usize,
    lease: Option<Arc<dyn Any + Send + Sync>>,
}

impl fmt::Debug for DataFrame {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("DataFrame")
            .field("rows", &self.rows)
            .field("columns", &self.names)
            .finish_non_exhaustive()
    }
}

impl DataFrame {
    /// Builds a frame whose row count is the length of its first column
    /// (zero if there are no columns).
    pub fn new<N: Into<String>>(
        columns: impl IntoIterator<Item = (N, Column)>,
    ) -> Result<Self, LazyError> {
        let columns: Vec<(N, Column)> = columns.into_iter().collect();
        let rows = columns.first().map_or(0, |(_, c)| c.len());
        Self::with_rows(rows, columns)
    }

    /// Builds a frame with an explicit row count (needed when no columns are
    /// selected but the rows still matter).
    pub fn with_rows<N: Into<String>>(
        rows: usize,
        columns: impl IntoIterator<Item = (N, Column)>,
    ) -> Result<Self, LazyError> {
        let mut names = Vec::new();
        let mut cols = Vec::new();
        for (name, col) in columns {
            let name = name.into();
            if col.len() != rows {
                return Err(LazyError::Invalid(format!(
                    "column `{name}` has {} rows, expected {rows}",
                    col.len()
                )));
            }
            if names.contains(&name) {
                return Err(LazyError::DuplicateColumn(name));
            }
            names.push(name);
            cols.push(col);
        }
        Ok(Self { names, columns: cols, rows, lease: None })
    }

    /// Attaches a lease that lives as long as this frame and its derivatives.
    #[must_use]
    pub fn with_lease(mut self, lease: impl Any + Send + Sync) -> Self {
        self.lease = Some(Arc::new(lease));
        self
    }

    /// Number of rows.
    pub fn len(&self) -> usize {
        self.rows
    }

    /// Whether the frame has no rows.
    pub fn is_empty(&self) -> bool {
        self.rows == 0
    }

    /// Column names in order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Looks up a column by name.
    pub fn column(&self, name: &str) -> Option<&Column> {
        self.names.iter().position(|n| n == name).map(|i| &self.columns[i])
    }

    /// The schema of this frame.
    pub fn schema(&self) -> Schema {
        self.names
            .iter()
            .zip(&self.columns)
            .map(|(n, c)| Field { name: n.clone(), dtype: c.dtype() })
            .collect()
    }

    /// A new frame with the same lease, the given rows and columns.
    pub(crate) fn derive(
        &self,
        rows: usize,
        columns: Vec<(String, Column)>,
    ) -> Result<Self, LazyError> {
        let mut out = Self::with_rows(rows, columns)?;
        out.lease.clone_from(&self.lease);
        Ok(out)
    }

    /// Replaces or appends a column, keeping the lease.
    pub(crate) fn set_column(mut self, name: String, col: Column) -> Result<Self, LazyError> {
        if col.len() != self.rows {
            return Err(LazyError::Invalid(format!("column `{name}` has the wrong length")));
        }
        if let Some(i) = self.names.iter().position(|n| *n == name) {
            self.columns[i] = col;
        } else {
            self.names.push(name);
            self.columns.push(col);
        }
        Ok(self)
    }

    /// Keeps the rows where `mask` is true.
    pub(crate) fn filter_rows(&self, mask: &[bool]) -> Result<Self, LazyError> {
        let rows = mask.iter().filter(|&&m| m).count();
        let cols = self
            .names
            .iter()
            .cloned()
            .zip(self.columns.iter().map(|c| c.filter(mask)))
            .collect::<Vec<_>>();
        self.derive(rows, cols)
    }
}

/// A source of rows that can be read in ranges, so the executor never needs the
/// whole dataset in memory.
pub trait Source: Send + Sync {
    /// Schema of every column the source can provide.
    fn schema(&self) -> Schema;

    /// Total number of rows.
    fn num_rows(&self) -> usize;

    /// Reads the rows in `rows` (clamped to the source) restricted to the named
    /// `projection` columns. The result must have `rows.len()` rows, or fewer
    /// when the range runs past the end of the source.
    fn read(&self, rows: Range<usize>, projection: &[String]) -> Result<DataFrame, LazyError>;
}

impl Source for DataFrame {
    fn schema(&self) -> Schema {
        DataFrame::schema(self)
    }

    fn num_rows(&self) -> usize {
        self.rows
    }

    fn read(&self, rows: Range<usize>, projection: &[String]) -> Result<DataFrame, LazyError> {
        let end = rows.end.min(self.rows);
        let start = rows.start.min(end);
        let mut cols = Vec::with_capacity(projection.len());
        for name in projection {
            let col = self.column(name).ok_or_else(|| LazyError::UnknownColumn(name.clone()))?;
            cols.push((name.clone(), col.slice(start..end)));
        }
        DataFrame::with_rows(end - start, cols)
    }
}

/// A deferred computation over a source. Nothing runs until [`LazyFrame::collect`].
#[derive(Clone)]
pub struct LazyFrame {
    plan: Plan,
    options: CollectOptions,
}

impl LazyFrame {
    /// Starts a plan that scans `source`.
    pub fn scan<S: Source + 'static>(source: S) -> Self {
        Self::scan_arc(Arc::new(source))
    }

    /// Starts a plan that scans a shared source.
    pub fn scan_arc(source: Arc<dyn Source>) -> Self {
        Self { plan: Plan::Scan { source, projection: None }, options: CollectOptions::default() }
    }

    /// Keeps the rows for which `predicate` is true.
    #[must_use]
    pub fn filter(self, predicate: Expr) -> Self {
        self.map_plan(|input| Plan::Filter { input: Box::new(input), predicate })
    }

    /// Keeps only the given outputs. Each expression needs a name: a bare
    /// column or an expression created with [`Expr::alias`].
    #[must_use]
    pub fn select(self, exprs: impl IntoIterator<Item = Expr>) -> Self {
        let exprs: Vec<Expr> = exprs.into_iter().collect();
        self.map_plan(|input| Plan::Project { input: Box::new(input), exprs })
    }

    /// Adds a column, or replaces the column with the same name.
    #[must_use]
    pub fn with_column(self, expr: Expr) -> Self {
        self.map_plan(|input| Plan::WithColumn { input: Box::new(input), expr })
    }

    /// Starts a group-by on the given key columns (`i64` or `bool`).
    pub fn group_by<S: Into<String>>(self, keys: impl IntoIterator<Item = S>) -> GroupBy {
        GroupBy { frame: self, keys: keys.into_iter().map(Into::into).collect() }
    }

    /// Replaces the execution options.
    #[must_use]
    pub fn with_options(mut self, options: CollectOptions) -> Self {
        self.options = options;
        self
    }

    /// Sets the number of rows per chunk.
    #[must_use]
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.options.chunk_size = chunk_size;
        self
    }

    /// The logical plan as built, before optimisation.
    pub fn plan(&self) -> &Plan {
        &self.plan
    }

    /// Renders the optimised plan as an indented tree, one operator per line.
    pub fn explain(&self) -> String {
        format_plan(&optimize(self.plan.clone()))
    }

    /// Renders the plan as built, without optimisation.
    pub fn explain_unoptimized(&self) -> String {
        format_plan(&self.plan)
    }

    /// Validates, optimises and executes the plan.
    pub fn collect(self) -> Result<DataFrame, LazyError> {
        validate(&self.plan)?;
        let plan = optimize(self.plan);
        execute(&plan, &self.options)
    }

    /// Validates and executes the plan without optimisation. Results must equal
    /// [`LazyFrame::collect`]; this exists for testing and debugging.
    pub fn collect_unoptimized(self) -> Result<DataFrame, LazyError> {
        validate(&self.plan)?;
        execute(&self.plan, &self.options)
    }

    fn map_plan(self, f: impl FnOnce(Plan) -> Plan) -> Self {
        Self { plan: f(self.plan), options: self.options }
    }
}

/// Intermediate state of a group-by, waiting for [`GroupBy::agg`].
#[derive(Clone)]
pub struct GroupBy {
    frame: LazyFrame,
    keys: Vec<String>,
}

impl GroupBy {
    /// Aggregates each group. The result has the key columns followed by one
    /// column per aggregate.
    #[must_use]
    pub fn agg(self, aggs: impl IntoIterator<Item = AggExpr>) -> LazyFrame {
        let aggs: Vec<AggExpr> = aggs.into_iter().collect();
        let keys = self.keys;
        self.frame.map_plan(|input| Plan::Aggregate { input: Box::new(input), keys, aggs })
    }
}
