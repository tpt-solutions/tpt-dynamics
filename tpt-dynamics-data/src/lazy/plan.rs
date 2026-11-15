//! Logical plan, schema validation and the `explain` renderer.

use std::fmt;
use std::sync::Arc;

use super::expr::{eval, Expr};
use super::frame::{Column, DType, DataFrame, Field, LazyError, Schema, Source};

/// A logical operator. Plans are trees: every node has at most one input.
#[derive(Clone)]
pub enum Plan {
    /// Reads from a source. `projection: None` means every column.
    Scan {
        /// Data source.
        source: Arc<dyn Source>,
        /// Columns to read, or `None` for all of them.
        projection: Option<Vec<String>>,
    },
    /// Keeps the rows where `predicate` is true.
    Filter {
        /// Input plan.
        input: Box<Plan>,
        /// Boolean expression.
        predicate: Expr,
    },
    /// Keeps only the listed outputs (each named).
    Project {
        /// Input plan.
        input: Box<Plan>,
        /// Named output expressions.
        exprs: Vec<Expr>,
    },
    /// Adds or replaces one column.
    WithColumn {
        /// Input plan.
        input: Box<Plan>,
        /// Named expression giving the new column.
        expr: Expr,
    },
    /// Groups rows by key columns and aggregates each group.
    Aggregate {
        /// Input plan.
        input: Box<Plan>,
        /// Key column names.
        keys: Vec<String>,
        /// Aggregates to compute per group.
        aggs: Vec<AggExpr>,
    },
}

impl Plan {
    /// The single input of this node, if any.
    pub(crate) fn input(&self) -> Option<&Plan> {
        match self {
            Self::Scan { .. } => None,
            Self::Filter { input, .. }
            | Self::Project { input, .. }
            | Self::WithColumn { input, .. }
            | Self::Aggregate { input, .. } => Some(input),
        }
    }
}

/// Aggregate function.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AggKind {
    /// Number of rows in the group.
    Count,
    /// Sum of values.
    Sum,
    /// Arithmetic mean (always `f64`).
    Mean,
    /// Minimum value.
    Min,
    /// Maximum value.
    Max,
}

impl AggKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Count => "count",
            Self::Sum => "sum",
            Self::Mean => "mean",
            Self::Min => "min",
            Self::Max => "max",
        }
    }
}

/// One aggregate of a group-by, with its output name.
#[derive(Clone, Debug, PartialEq)]
pub struct AggExpr {
    kind: AggKind,
    input: Option<Expr>,
    name: String,
}

impl AggExpr {
    /// Builds an aggregate from its parts. Used when decoding a shipped plan.
    pub(crate) fn from_parts(kind: AggKind, input: Option<Expr>, name: String) -> Self {
        Self { kind, input, name }
    }

    /// Number of rows per group.
    pub fn count() -> Self {
        Self { kind: AggKind::Count, input: None, name: "count".into() }
    }

    /// Sum of `input` per group.
    pub fn sum(input: Expr) -> Self {
        Self::with_input(AggKind::Sum, input)
    }

    /// Mean of `input` per group.
    pub fn mean(input: Expr) -> Self {
        Self::with_input(AggKind::Mean, input)
    }

    /// Minimum of `input` per group.
    pub fn min(input: Expr) -> Self {
        Self::with_input(AggKind::Min, input)
    }

    /// Maximum of `input` per group.
    pub fn max(input: Expr) -> Self {
        Self::with_input(AggKind::Max, input)
    }

    fn with_input(kind: AggKind, input: Expr) -> Self {
        let name = match input.output_name() {
            Some(n) => format!("{n}_{}", kind.as_str()),
            None => format!("{}_{input}", kind.as_str()),
        };
        Self { kind, input: Some(input), name }
    }

    /// Sets the output column name.
    #[must_use]
    pub fn alias(mut self, name: impl Into<String>) -> Self {
        self.name = name.into();
        self
    }

    /// Output column name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Aggregate function.
    pub fn kind(&self) -> AggKind {
        self.kind
    }

    /// Input expression (`None` for `count`).
    pub fn input(&self) -> Option<&Expr> {
        self.input.as_ref()
    }

    pub(crate) fn result_dtype(&self, input: Option<DType>) -> Result<DType, LazyError> {
        let numeric = |d: Option<DType>| match d {
            Some(DType::F64) => Ok(DType::F64),
            Some(DType::I64) => Ok(DType::I64),
            _ => Err(LazyError::Type(format!("{} needs a numeric input", self.kind.as_str()))),
        };
        match self.kind {
            AggKind::Count => Ok(DType::I64),
            AggKind::Mean => numeric(input).map(|_| DType::F64),
            AggKind::Sum | AggKind::Min | AggKind::Max => numeric(input),
        }
    }
}

impl fmt::Display for AggExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.input {
            Some(e) => write!(f, "{}={}({e})", self.name, self.kind.as_str()),
            None => write!(f, "{}={}()", self.name, self.kind.as_str()),
        }
    }
}

fn push_field(out: &mut Schema, name: String, dtype: DType) -> Result<(), LazyError> {
    if out.iter().any(|f| f.name == name) {
        return Err(LazyError::DuplicateColumn(name));
    }
    out.push(Field { name, dtype });
    Ok(())
}

fn output_name(expr: &Expr) -> Result<String, LazyError> {
    expr.output_name().map(str::to_string).ok_or_else(|| {
        LazyError::Invalid(format!("expression `{expr}` needs a name; use .alias(...)"))
    })
}

fn empty_frame(schema: &Schema) -> Result<DataFrame, LazyError> {
    DataFrame::with_rows(
        0,
        schema.iter().map(|f| (f.name.clone(), Column::empty(f.dtype))).collect::<Vec<_>>(),
    )
}

/// Computes the output schema of `plan`, checking every column reference and
/// operand type. Works on zero-row frames, so it also validates plans over empty
/// sources.
pub(crate) fn validate(plan: &Plan) -> Result<Schema, LazyError> {
    match plan {
        Plan::Scan { source, projection } => {
            let schema = source.schema();
            match projection {
                None => Ok(schema),
                Some(names) => names
                    .iter()
                    .map(|n| {
                        schema
                            .iter()
                            .find(|f| &f.name == n)
                            .cloned()
                            .ok_or_else(|| LazyError::UnknownColumn(n.clone()))
                    })
                    .collect(),
            }
        }
        Plan::Filter { input, predicate } => {
            let schema = validate(input)?;
            let frame = empty_frame(&schema)?;
            let mask = eval(predicate, &frame)?;
            if mask.dtype() != DType::Bool {
                return Err(LazyError::Type(format!(
                    "filter predicate `{predicate}` must be boolean"
                )));
            }
            Ok(schema)
        }
        Plan::Project { input, exprs } => {
            let frame = empty_frame(&validate(input)?)?;
            let mut out = Vec::with_capacity(exprs.len());
            for e in exprs {
                let name = output_name(e)?;
                let dtype = eval(e, &frame)?.dtype();
                push_field(&mut out, name, dtype)?;
            }
            Ok(out)
        }
        Plan::WithColumn { input, expr } => {
            let mut schema = validate(input)?;
            let name = output_name(expr)?;
            let dtype = eval(expr, &empty_frame(&schema)?)?.dtype();
            match schema.iter_mut().find(|f| f.name == name) {
                Some(f) => f.dtype = dtype,
                None => schema.push(Field { name, dtype }),
            }
            Ok(schema)
        }
        Plan::Aggregate { input, keys, aggs } => {
            let frame = empty_frame(&validate(input)?)?;
            if keys.is_empty() {
                return Err(LazyError::Invalid("group_by needs at least one key column".into()));
            }
            let mut out = Vec::new();
            for k in keys {
                let col = frame.column(k).ok_or_else(|| LazyError::UnknownColumn(k.clone()))?;
                if col.dtype() == DType::F64 {
                    return Err(LazyError::Type(format!("cannot group by f64 column `{k}`")));
                }
                push_field(&mut out, k.clone(), col.dtype())?;
            }
            for a in aggs {
                let input_dtype = match a.input() {
                    Some(e) => Some(eval(e, &frame)?.dtype()),
                    None => None,
                };
                let dtype = a.result_dtype(input_dtype)?;
                push_field(&mut out, a.name().to_string(), dtype)?;
            }
            Ok(out)
        }
    }
}

/// Renders a plan as an indented tree, root first.
pub(crate) fn format_plan(plan: &Plan) -> String {
    let mut out = String::new();
    write_tree(plan, 0, &mut out);
    out
}

fn write_tree(plan: &Plan, depth: usize, out: &mut String) {
    let line = match plan {
        Plan::Scan { source, projection } => {
            let cols = match projection {
                Some(p) => format!("[{}]", p.join(", ")),
                None => "[*]".to_string(),
            };
            format!("Scan {cols} rows={}", source.num_rows())
        }
        Plan::Filter { predicate, .. } => format!("Filter {predicate}"),
        Plan::Project { exprs, .. } => {
            let items: Vec<String> = exprs.iter().map(ToString::to_string).collect();
            format!("Select [{}]", items.join(", "))
        }
        Plan::WithColumn { expr, .. } => format!("WithColumn {expr}"),
        Plan::Aggregate { keys, aggs, .. } => {
            let items: Vec<String> = aggs.iter().map(ToString::to_string).collect();
            format!("Aggregate by=[{}] aggs=[{}]", keys.join(", "), items.join(", "))
        }
    };
    out.push_str(&"  ".repeat(depth));
    out.push_str(&line);
    out.push('\n');
    if let Some(input) = plan.input() {
        write_tree(input, depth + 1, out);
    }
}
