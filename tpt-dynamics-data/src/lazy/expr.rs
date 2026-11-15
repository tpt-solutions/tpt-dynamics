//! Expression trees and their column-wise (vectorised) evaluator.
//!
//! Every operator consumes whole columns, so a chunk is processed with one
//! tight loop per operator rather than one interpreted step per value.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Not, Sub};

use super::frame::{Column, DataFrame, LazyError};

/// A constant value.
#[derive(Clone, Debug, PartialEq)]
pub enum Scalar {
    /// Floating-point constant.
    F64(f64),
    /// Integer constant.
    I64(i64),
    /// Boolean constant.
    Bool(bool),
}

impl From<f64> for Scalar {
    fn from(v: f64) -> Self {
        Self::F64(v)
    }
}

impl From<i64> for Scalar {
    fn from(v: i64) -> Self {
        Self::I64(v)
    }
}

impl From<i32> for Scalar {
    fn from(v: i32) -> Self {
        Self::I64(i64::from(v))
    }
}

impl From<bool> for Scalar {
    fn from(v: bool) -> Self {
        Self::Bool(v)
    }
}

impl fmt::Display for Scalar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::F64(v) => write!(f, "{v}"),
            Self::I64(v) => write!(f, "{v}"),
            Self::Bool(v) => write!(f, "{v}"),
        }
    }
}

/// Binary operators.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinOp {
    /// `+`
    Add,
    /// `-`
    Sub,
    /// `*`
    Mul,
    /// `/` (always produces `f64`)
    Div,
    /// `==`
    Eq,
    /// `!=`
    Ne,
    /// `<`
    Lt,
    /// `<=`
    Le,
    /// `>`
    Gt,
    /// `>=`
    Ge,
    /// Logical and.
    And,
    /// Logical or.
    Or,
}

impl BinOp {
    fn symbol(self) -> &'static str {
        match self {
            Self::Add => "+",
            Self::Sub => "-",
            Self::Mul => "*",
            Self::Div => "/",
            Self::Eq => "==",
            Self::Ne => "!=",
            Self::Lt => "<",
            Self::Le => "<=",
            Self::Gt => ">",
            Self::Ge => ">=",
            Self::And => "&",
            Self::Or => "|",
        }
    }
}

/// Unary operators.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnOp {
    /// Arithmetic negation.
    Neg,
    /// Logical not.
    Not,
}

impl UnOp {
    fn symbol(self) -> &'static str {
        match self {
            Self::Neg => "-",
            Self::Not => "!",
        }
    }
}

/// Built-in physics functions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Func {
    /// Invariant mass of two particles from `(pt, eta, phi, m)` each.
    InvariantMass,
    /// Angular distance `sqrt(d_eta^2 + d_phi^2)` with `d_phi` wrapped to `[-pi, pi)`.
    DeltaR,
    /// Transverse momentum `sqrt(px^2 + py^2)`.
    TransverseMomentum,
}

impl Func {
    fn name(self) -> &'static str {
        match self {
            Self::InvariantMass => "invariant_mass",
            Self::DeltaR => "delta_r",
            Self::TransverseMomentum => "transverse_momentum",
        }
    }

    fn arity(self) -> usize {
        match self {
            Self::InvariantMass => 8,
            Self::DeltaR => 4,
            Self::TransverseMomentum => 2,
        }
    }
}

/// An expression over the columns of a frame.
#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    /// Reference to a column by name.
    Col(String),
    /// Constant broadcast to every row.
    Lit(Scalar),
    /// Binary operation.
    Binary {
        /// Operator.
        op: BinOp,
        /// Left operand.
        lhs: Box<Expr>,
        /// Right operand.
        rhs: Box<Expr>,
    },
    /// Unary operation.
    Unary {
        /// Operator.
        op: UnOp,
        /// Operand.
        arg: Box<Expr>,
    },
    /// Physics function call.
    Call {
        /// Function.
        func: Func,
        /// Arguments, in the order of the constructor.
        args: Vec<Expr>,
    },
    /// Names the output of an expression.
    Alias {
        /// Aliased expression.
        expr: Box<Expr>,
        /// Output name.
        name: String,
    },
}

/// Column reference.
pub fn col(name: impl Into<String>) -> Expr {
    Expr::Col(name.into())
}

/// Constant.
pub fn lit(value: impl Into<Scalar>) -> Expr {
    Expr::Lit(value.into())
}

/// Invariant mass of two particles, each given as `(pt, eta, phi, m)`.
/// Returns `sqrt(max(0, (E1+E2)^2 - |p1+p2|^2))`.
pub fn invariant_mass(
    pt1: Expr,
    eta1: Expr,
    phi1: Expr,
    m1: Expr,
    pt2: Expr,
    eta2: Expr,
    phi2: Expr,
    m2: Expr,
) -> Expr {
    func_call(Func::InvariantMass, vec![pt1, eta1, phi1, m1, pt2, eta2, phi2, m2])
}

/// Angular distance between two objects given by `(eta, phi)`.
pub fn delta_r(eta1: Expr, phi1: Expr, eta2: Expr, phi2: Expr) -> Expr {
    func_call(Func::DeltaR, vec![eta1, phi1, eta2, phi2])
}

/// Transverse momentum `sqrt(px^2 + py^2)`.
pub fn transverse_momentum(px: Expr, py: Expr) -> Expr {
    func_call(Func::TransverseMomentum, vec![px, py])
}

fn func_call(func: Func, args: Vec<Expr>) -> Expr {
    Expr::Call { func, args }
}

impl From<f64> for Expr {
    fn from(v: f64) -> Self {
        Expr::Lit(Scalar::F64(v))
    }
}

impl From<i64> for Expr {
    fn from(v: i64) -> Self {
        Expr::Lit(Scalar::I64(v))
    }
}

impl From<bool> for Expr {
    fn from(v: bool) -> Self {
        Expr::Lit(Scalar::Bool(v))
    }
}

impl Expr {
    fn binary(op: BinOp, lhs: Expr, rhs: Expr) -> Expr {
        Expr::Binary { op, lhs: Box::new(lhs), rhs: Box::new(rhs) }
    }

    /// Names this expression's output.
    #[must_use]
    pub fn alias(self, name: impl Into<String>) -> Expr {
        Expr::Alias { expr: Box::new(self), name: name.into() }
    }

    /// `self > rhs`
    #[must_use]
    pub fn gt(self, rhs: impl Into<Expr>) -> Expr {
        Self::binary(BinOp::Gt, self, rhs.into())
    }

    /// `self >= rhs`
    #[must_use]
    pub fn ge(self, rhs: impl Into<Expr>) -> Expr {
        Self::binary(BinOp::Ge, self, rhs.into())
    }

    /// `self < rhs`
    #[must_use]
    pub fn lt(self, rhs: impl Into<Expr>) -> Expr {
        Self::binary(BinOp::Lt, self, rhs.into())
    }

    /// `self <= rhs`
    #[must_use]
    pub fn le(self, rhs: impl Into<Expr>) -> Expr {
        Self::binary(BinOp::Le, self, rhs.into())
    }

    /// `self == rhs`
    #[must_use]
    pub fn eq_to(self, rhs: impl Into<Expr>) -> Expr {
        Self::binary(BinOp::Eq, self, rhs.into())
    }

    /// `self != rhs`
    #[must_use]
    pub fn ne_to(self, rhs: impl Into<Expr>) -> Expr {
        Self::binary(BinOp::Ne, self, rhs.into())
    }

    /// Logical and.
    #[must_use]
    pub fn and(self, rhs: Expr) -> Expr {
        Self::binary(BinOp::And, self, rhs)
    }

    /// Logical or.
    #[must_use]
    pub fn or(self, rhs: Expr) -> Expr {
        Self::binary(BinOp::Or, self, rhs)
    }

    /// The output name, if the expression has one (a column or an alias).
    pub fn output_name(&self) -> Option<&str> {
        match self {
            Self::Col(n) => Some(n),
            Self::Alias { name, .. } => Some(name),
            _ => None,
        }
    }

    /// The expression without an outer alias.
    pub(crate) fn unaliased(&self) -> &Expr {
        match self {
            Self::Alias { expr, .. } => expr,
            other => other,
        }
    }

    /// All column names the expression reads.
    pub(crate) fn referenced_columns(&self) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        self.collect_refs(&mut out);
        out
    }

    fn collect_refs(&self, out: &mut BTreeSet<String>) {
        match self {
            Self::Col(n) => {
                out.insert(n.clone());
            }
            Self::Lit(_) => {}
            Self::Binary { lhs, rhs, .. } => {
                lhs.collect_refs(out);
                rhs.collect_refs(out);
            }
            Self::Unary { arg, .. } => arg.collect_refs(out),
            Self::Call { args, .. } => args.iter().for_each(|a| a.collect_refs(out)),
            Self::Alias { expr, .. } => expr.collect_refs(out),
        }
    }
}

macro_rules! bin_op {
    ($trait:ident, $method:ident, $op:expr) => {
        impl $trait<Expr> for Expr {
            type Output = Expr;
            fn $method(self, rhs: Expr) -> Expr {
                Expr::binary($op, self, rhs)
            }
        }
    };
}

bin_op!(Add, add, BinOp::Add);
bin_op!(Sub, sub, BinOp::Sub);
bin_op!(Mul, mul, BinOp::Mul);
bin_op!(Div, div, BinOp::Div);

impl Neg for Expr {
    type Output = Expr;
    fn neg(self) -> Expr {
        Expr::Unary { op: UnOp::Neg, arg: Box::new(self) }
    }
}

impl Not for Expr {
    type Output = Expr;
    fn not(self) -> Expr {
        Expr::Unary { op: UnOp::Not, arg: Box::new(self) }
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Col(n) => write!(f, "{n}"),
            Self::Lit(s) => write!(f, "{s}"),
            Self::Binary { op, lhs, rhs } => write!(f, "({lhs} {} {rhs})", op.symbol()),
            Self::Unary { op, arg } => write!(f, "{}{arg}", op.symbol()),
            Self::Call { func, args } => {
                write!(f, "{}(", func.name())?;
                for (i, a) in args.iter().enumerate() {
                    if i > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{a}")?;
                }
                write!(f, ")")
            }
            Self::Alias { expr, name } => write!(f, "{expr} AS {name}"),
        }
    }
}

/// Evaluates `expr` over every row of `df`. Column references are borrowed;
/// every other node produces a new column of `df.len()` values.
pub(crate) fn eval<'a>(expr: &Expr, df: &'a DataFrame) -> Result<Cow<'a, Column>, LazyError> {
    match expr {
        Expr::Col(name) => {
            df.column(name).map(Cow::Borrowed).ok_or_else(|| LazyError::UnknownColumn(name.clone()))
        }
        Expr::Lit(s) => Ok(Cow::Owned(broadcast(s, df.len()))),
        Expr::Alias { expr, .. } => eval(expr, df),
        Expr::Unary { op, arg } => {
            let a = eval(arg, df)?;
            unary(*op, &a).map(Cow::Owned)
        }
        Expr::Binary { op, lhs, rhs } => {
            let l = eval(lhs, df)?;
            let r = eval(rhs, df)?;
            binary(*op, &l, &r).map(Cow::Owned)
        }
        Expr::Call { func, args } => {
            if args.len() != func.arity() {
                return Err(LazyError::Invalid(format!(
                    "{} takes {} arguments, got {}",
                    func.name(),
                    func.arity(),
                    args.len()
                )));
            }
            let vals = args.iter().map(|a| eval(a, df)).collect::<Result<Vec<_>, _>>()?;
            call_func(*func, &vals, df.len()).map(Cow::Owned)
        }
    }
}

fn broadcast(s: &Scalar, n: usize) -> Column {
    match s {
        Scalar::F64(v) => Column::F64(vec![*v; n]),
        Scalar::I64(v) => Column::I64(vec![*v; n]),
        Scalar::Bool(v) => Column::Bool(vec![*v; n]),
    }
}

fn zip_map<A: Copy, B: Copy, O>(a: &[A], b: &[B], f: impl Fn(A, B) -> O) -> Vec<O> {
    a.iter().zip(b).map(|(&x, &y)| f(x, y)).collect()
}

fn as_f64(col: &Column) -> Result<Cow<'_, [f64]>, LazyError> {
    match col {
        Column::F64(v) => Ok(Cow::Borrowed(v.as_slice())),
        Column::I64(v) => Ok(Cow::Owned(v.iter().map(|&x| x as f64).collect())),
        Column::Bool(_) => Err(LazyError::Type("expected a numeric column, found bool".into())),
    }
}

fn unary(op: UnOp, a: &Column) -> Result<Column, LazyError> {
    match (op, a) {
        (UnOp::Neg, Column::F64(v)) => Ok(Column::F64(v.iter().map(|x| -x).collect())),
        (UnOp::Neg, Column::I64(v)) => {
            Ok(Column::I64(v.iter().map(|x| x.wrapping_neg()).collect()))
        }
        (UnOp::Not, Column::Bool(v)) => Ok(Column::Bool(v.iter().map(|x| !x).collect())),
        _ => Err(LazyError::Type(format!(
            "operator {} does not apply to {}",
            op.symbol(),
            a.dtype()
        ))),
    }
}

fn binary(op: BinOp, l: &Column, r: &Column) -> Result<Column, LazyError> {
    if l.len() != r.len() {
        return Err(LazyError::Invalid("operand lengths differ".into()));
    }
    match op {
        BinOp::And | BinOp::Or => match (l, r) {
            (Column::Bool(a), Column::Bool(b)) => Ok(Column::Bool(if op == BinOp::And {
                zip_map(a, b, |x, y| x && y)
            } else {
                zip_map(a, b, |x, y| x || y)
            })),
            _ => Err(LazyError::Type(format!("operator {} needs bool operands", op.symbol()))),
        },
        BinOp::Add | BinOp::Sub | BinOp::Mul | BinOp::Div => arith(op, l, r),
        _ => compare(op, l, r),
    }
}

fn arith(op: BinOp, lhs: &Column, rhs: &Column) -> Result<Column, LazyError> {
    if let (Column::I64(xs), Column::I64(ys)) = (lhs, rhs) {
        if op != BinOp::Div {
            let func: fn(i64, i64) -> i64 = match op {
                BinOp::Add => i64::wrapping_add,
                BinOp::Sub => i64::wrapping_sub,
                _ => i64::wrapping_mul,
            };
            return Ok(Column::I64(zip_map(xs, ys, func)));
        }
    }
    let xs = as_f64(lhs)?;
    let ys = as_f64(rhs)?;
    let func: fn(f64, f64) -> f64 = match op {
        BinOp::Add => |x, y| x + y,
        BinOp::Sub => |x, y| x - y,
        BinOp::Mul => |x, y| x * y,
        _ => |x, y| x / y,
    };
    Ok(Column::F64(zip_map(&xs, &ys, func)))
}

fn compare(op: BinOp, l: &Column, r: &Column) -> Result<Column, LazyError> {
    match (l, r) {
        (Column::I64(a), Column::I64(b)) => Ok(Column::Bool(cmp_slices(op, a, b))),
        (Column::Bool(a), Column::Bool(b)) => Ok(Column::Bool(cmp_slices(op, a, b))),
        (Column::F64(_) | Column::I64(_), Column::F64(_) | Column::I64(_)) => {
            let a = as_f64(l)?;
            let b = as_f64(r)?;
            Ok(Column::Bool(cmp_slices(op, &a, &b)))
        }
        _ => Err(LazyError::Type("cannot compare bool with a number".into())),
    }
}

fn cmp_slices<T: PartialOrd + Copy>(op: BinOp, a: &[T], b: &[T]) -> Vec<bool> {
    let f: fn(T, T) -> bool = match op {
        BinOp::Eq => |x, y| x == y,
        BinOp::Ne => |x, y| x != y,
        BinOp::Lt => |x, y| x < y,
        BinOp::Le => |x, y| x <= y,
        BinOp::Gt => |x, y| x > y,
        _ => |x, y| x >= y,
    };
    zip_map(a, b, f)
}

fn call_func(func: Func, args: &[Cow<'_, Column>], n: usize) -> Result<Column, LazyError> {
    let cols = args.iter().map(|c| as_f64(c)).collect::<Result<Vec<_>, _>>()?;
    let s: Vec<&[f64]> = cols.iter().map(|c| &c[..]).collect();
    let out: Vec<f64> = match func {
        Func::InvariantMass => (0..n)
            .map(|i| {
                invariant_mass_row(
                    [s[0][i], s[1][i], s[2][i], s[3][i]],
                    [s[4][i], s[5][i], s[6][i], s[7][i]],
                )
            })
            .collect(),
        Func::DeltaR => (0..n).map(|i| delta_r_row(s[0][i], s[1][i], s[2][i], s[3][i])).collect(),
        Func::TransverseMomentum => zip_map(s[0], s[1], f64::hypot),
    };
    Ok(Column::F64(out))
}

/// Four-momentum from `(pt, eta, phi, m)`, returned as `[px, py, pz, E]`.
fn four_momentum([pt, eta, phi, m]: [f64; 4]) -> [f64; 4] {
    let px = pt * phi.cos();
    let py = pt * phi.sin();
    let pz = pt * eta.sinh();
    let e = (px * px + py * py + pz * pz + m * m).sqrt();
    [px, py, pz, e]
}

fn invariant_mass_row(p1: [f64; 4], p2: [f64; 4]) -> f64 {
    let a = four_momentum(p1);
    let b = four_momentum(p2);
    let px = a[0] + b[0];
    let py = a[1] + b[1];
    let pz = a[2] + b[2];
    let e = a[3] + b[3];
    (e * e - px * px - py * py - pz * pz).max(0.0).sqrt()
}

fn delta_r_row(eta1: f64, phi1: f64, eta2: f64, phi2: f64) -> f64 {
    let dphi = (phi1 - phi2 + PI).rem_euclid(TAU) - PI;
    let deta = eta1 - eta2;
    deta.hypot(dphi)
}

/// Folds operators whose operands are all literals into a single literal.
///
/// The fold reuses the evaluator's own `unary` and `binary` kernels on one-row
/// columns, so a folded constant is bit-identical to what the runtime would have
/// produced for every row. Anything that would fail at runtime (a type error, for
/// example) is left in place, so it still fails when the plan runs. Column
/// references and physics calls are never folded.
pub(crate) fn fold_constants(expr: Expr) -> Expr {
    match expr {
        Expr::Col(_) | Expr::Lit(_) => expr,
        Expr::Alias { expr, name } => Expr::Alias { expr: Box::new(fold_constants(*expr)), name },
        Expr::Call { func, args } => {
            Expr::Call { func, args: args.into_iter().map(fold_constants).collect() }
        }
        Expr::Unary { op, arg } => {
            let arg = fold_constants(*arg);
            if let Expr::Lit(s) = &arg {
                if let Ok(c) = unary(op, &scalar_column(s)) {
                    if let Some(folded) = column_scalar(&c) {
                        return Expr::Lit(folded);
                    }
                }
            }
            Expr::Unary { op, arg: Box::new(arg) }
        }
        Expr::Binary { op, lhs, rhs } => {
            let lhs = fold_constants(*lhs);
            let rhs = fold_constants(*rhs);
            if let (Expr::Lit(a), Expr::Lit(b)) = (&lhs, &rhs) {
                if let Ok(c) = binary(op, &scalar_column(a), &scalar_column(b)) {
                    if let Some(folded) = column_scalar(&c) {
                        return Expr::Lit(folded);
                    }
                }
            }
            Expr::Binary { op, lhs: Box::new(lhs), rhs: Box::new(rhs) }
        }
    }
}

/// A one-row column holding `s`.
fn scalar_column(s: &Scalar) -> Column {
    broadcast(s, 1)
}

/// The value of a one-row column, if it has exactly one row.
fn column_scalar(c: &Column) -> Option<Scalar> {
    match c {
        Column::F64(v) if v.len() == 1 => Some(Scalar::F64(v[0])),
        Column::I64(v) if v.len() == 1 => Some(Scalar::I64(v[0])),
        Column::Bool(v) if v.len() == 1 => Some(Scalar::Bool(v[0])),
        _ => None,
    }
}
