//! Filling a [`Histogram1D`] from a lazy plan.

use super::expr::Expr;
use super::frame::{Column, DataFrame, LazyError, LazyFrame};
use crate::hist1d::Histogram1D;

/// Temporary output name for the fill variable.
const X_NAME: &str = "__fill_x";
/// Temporary output name for the fill weight.
const W_NAME: &str = "__fill_w";

impl LazyFrame {
    /// Runs the plan and fills `hist` with the values of `x`, optionally weighted
    /// by `weight`.
    ///
    /// `x` and `weight` are expressions over the plan's output (for example
    /// `transverse_momentum(col("px"), col("py"))`). Filters, `with_column` and
    /// group-free projections in the plan apply before the fill. Integer values
    /// are converted to `f64`. Boolean expressions are rejected.
    ///
    /// Only the two fill columns are materialised, so memory is proportional to
    /// the number of selected rows, not to the full input.
    ///
    /// # Errors
    ///
    /// Returns [`LazyError::Type`] if `x` or `weight` is boolean, and any error
    /// from executing the plan.
    ///
    /// # Examples
    ///
    /// ```
    /// use tpt_dynamics_data::hist1d::Histogram1D;
    /// use tpt_dynamics_data::lazy::{col, Column, DataFrame, LazyFrame};
    ///
    /// let df = DataFrame::new(vec![
    ///     ("pt", Column::F64(vec![5.0, 15.0, 25.0, 35.0])),
    ///     ("keep", Column::Bool(vec![true, true, false, true])),
    /// ])
    /// .unwrap();
    /// let mut h = Histogram1D::uniform(4, 0.0, 40.0).unwrap();
    /// LazyFrame::scan(df)
    ///     .filter(col("keep"))
    ///     .fill_histogram(&mut h, col("pt"), None)
    ///     .unwrap();
    /// assert_eq!(h.integral(), 3.0);
    /// ```
    pub fn fill_histogram(
        self,
        hist: &mut Histogram1D,
        x: Expr,
        weight: Option<Expr>,
    ) -> Result<(), LazyError> {
        let mut exprs = vec![x.alias(X_NAME)];
        if let Some(w) = weight {
            exprs.push(w.alias(W_NAME));
        }
        let out = self.select(exprs).collect()?;
        fill_from_frame(hist, &out)
    }
}

/// Fills `hist` from the `__fill_x` column and the optional `__fill_w` column of
/// an already-evaluated frame.
fn fill_from_frame(hist: &mut Histogram1D, df: &DataFrame) -> Result<(), LazyError> {
    let xs = numeric_column(df, X_NAME)?;
    let ws = if df.column(W_NAME).is_some() { Some(numeric_column(df, W_NAME)?) } else { None };
    match ws {
        Some(ws) => {
            for (x, w) in xs.iter().zip(ws.iter()) {
                hist.fill_weighted(*x, *w);
            }
        }
        None => {
            for x in &xs {
                hist.fill(*x);
            }
        }
    }
    Ok(())
}

fn numeric_column(df: &DataFrame, name: &str) -> Result<Vec<f64>, LazyError> {
    match df.column(name) {
        Some(Column::F64(v)) => Ok(v.clone()),
        Some(Column::I64(v)) => Ok(v.iter().map(|&i| i as f64).collect()),
        Some(Column::Bool(_)) => Err(LazyError::Type(format!(
            "histogram fill needs a numeric expression, got boolean for `{name}`"
        ))),
        None => Err(LazyError::UnknownColumn(name.to_owned())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lazy::{col, lit};

    fn frame() -> DataFrame {
        DataFrame::new(vec![
            ("pt", Column::F64(vec![5.0, 15.0, 15.0, 25.0, 35.0, 99.0])),
            ("w", Column::F64(vec![1.0, 2.0, 0.5, 1.0, 3.0, 1.0])),
            ("n", Column::I64(vec![0, 1, 1, 2, 3, 4])),
            ("ok", Column::Bool(vec![true, true, true, false, true, true])),
        ])
        .unwrap()
    }

    #[test]
    fn lazy_fill_matches_eager_fill() {
        let df = frame();
        // Eager reference: apply the filter by hand, then fill.
        let mut eager = Histogram1D::uniform(4, 0.0, 40.0).unwrap();
        let (Column::F64(pt), Column::F64(w), Column::Bool(ok)) =
            (df.column("pt").unwrap(), df.column("w").unwrap(), df.column("ok").unwrap())
        else {
            panic!("unexpected column types");
        };
        for i in 0..pt.len() {
            if ok[i] {
                eager.fill_weighted(pt[i], w[i]);
            }
        }

        let mut lazy = Histogram1D::uniform(4, 0.0, 40.0).unwrap();
        LazyFrame::scan(df)
            .filter(col("ok"))
            .fill_histogram(&mut lazy, col("pt"), Some(col("w")))
            .unwrap();

        assert_eq!(lazy.contents(), eager.contents());
        assert_eq!(lazy.sum_weights_squared(), eager.sum_weights_squared());
        assert_eq!(lazy.overflow(), eager.overflow());
        assert_eq!(lazy.entries(), eager.entries());
    }

    #[test]
    fn integer_and_derived_expressions_are_accepted() {
        let mut h = Histogram1D::uniform(5, 0.0, 5.0).unwrap();
        LazyFrame::scan(frame()).fill_histogram(&mut h, col("n") + lit(0.5), None).unwrap();
        assert_eq!(h.entries(), 6);
        assert_eq!(h.integral(), 6.0);
    }

    #[test]
    fn boolean_fill_variable_is_rejected() {
        let mut h = Histogram1D::uniform(4, 0.0, 4.0).unwrap();
        let err = LazyFrame::scan(frame()).fill_histogram(&mut h, col("ok"), None).unwrap_err();
        assert!(matches!(err, LazyError::Type(_)), "{err:?}");
        assert_eq!(h.entries(), 0);
    }

    #[test]
    fn unknown_column_is_reported() {
        let mut h = Histogram1D::uniform(4, 0.0, 4.0).unwrap();
        let err =
            LazyFrame::scan(frame()).fill_histogram(&mut h, col("missing"), None).unwrap_err();
        assert!(matches!(err, LazyError::UnknownColumn(_)), "{err:?}");
    }
}
