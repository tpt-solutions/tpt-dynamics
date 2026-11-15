//! Lazy data frames (ADR-002).
//!
//! Build a logical plan with [`LazyFrame`] (scan, filter, select, `with_column`,
//! `group_by` + `agg`), inspect it with [`LazyFrame::explain`], then run it with
//! [`LazyFrame::collect`]. The optimiser applies predicate and projection
//! pushdown; the executor streams the source in chunks of bounded size and
//! parallelises across chunks with rayon. Results equal eager evaluation.
//!
//! ```
//! use tpt_dynamics_data::lazy::{col, AggExpr, Column, DataFrame, LazyFrame};
//!
//! let df = DataFrame::new(vec![
//!     ("pt", Column::F64(vec![10.0, 30.0, 50.0])),
//!     ("g", Column::I64(vec![1, 2, 2])),
//! ])
//! .unwrap();
//! let out = LazyFrame::scan(df)
//!     .filter(col("pt").gt(20.0))
//!     .group_by(["g"])
//!     .agg([AggExpr::sum(col("pt")).alias("pt_sum")])
//!     .collect()
//!     .unwrap();
//! assert_eq!(out.len(), 1);
//! ```

mod exec;
mod expr;
mod frame;
mod optimize;
mod plan;

#[cfg(test)]
mod tests;

pub use exec::CollectOptions;
pub use expr::{
    col, delta_r, invariant_mass, lit, transverse_momentum, BinOp, Expr, Func, Scalar, UnOp,
};
pub use frame::{Column, DType, DataFrame, Field, GroupBy, LazyError, LazyFrame, Schema, Source};
pub use plan::{AggExpr, AggKind, Plan};
