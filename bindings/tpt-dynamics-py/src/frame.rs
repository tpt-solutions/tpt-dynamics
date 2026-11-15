//! Lazy data frames: expressions, aggregates and the `LazyFrame` pipeline.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use pyo3::exceptions::{PyOSError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use tpt_dynamics_data::csv::{CsvOptions, CsvReader, CsvRecord};
use tpt_dynamics_data::lazy::{self, AggExpr, Column, DType, DataFrame, Expr, LazyFrame};

use crate::convert::{column_from_any, columns_to_dict};
use crate::errors::value_err;

/// Converts a Python scalar or [`PyExpr`] into an expression.
fn to_expr(obj: &Bound<'_, PyAny>) -> PyResult<Expr> {
    if let Ok(e) = obj.extract::<PyExpr>() {
        return Ok(e.inner);
    }
    if let Ok(b) = obj.extract::<bool>() {
        return Ok(lazy::lit(b));
    }
    if let Ok(i) = obj.extract::<i64>() {
        return Ok(lazy::lit(i));
    }
    if let Ok(f) = obj.extract::<f64>() {
        return Ok(lazy::lit(f));
    }
    Err(PyTypeError::new_err("expected an Expr, int, float or bool"))
}

/// An expression over the columns of a frame (built with `col` and `lit`).
#[pyclass(name = "Expr", module = "tpt_dynamics", from_py_object)]
#[derive(Clone)]
pub struct PyExpr {
    inner: Expr,
}

#[pymethods]
impl PyExpr {
    /// Renames the output of this expression.
    fn alias(&self, name: &str) -> Self {
        Self { inner: self.inner.clone().alias(name) }
    }

    /// Element-wise `self > other`.
    fn gt(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().gt(to_expr(other)?) })
    }

    /// Element-wise `self >= other`.
    fn ge(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().ge(to_expr(other)?) })
    }

    /// Element-wise `self < other`.
    fn lt(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().lt(to_expr(other)?) })
    }

    /// Element-wise `self <= other`.
    fn le(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().le(to_expr(other)?) })
    }

    /// Element-wise equality. A method rather than `==`, so that Python
    /// equality on expression objects keeps its usual meaning.
    fn eq_to(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().eq_to(to_expr(other)?) })
    }

    /// Element-wise inequality.
    fn ne_to(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().ne_to(to_expr(other)?) })
    }

    /// Logical and of two boolean expressions.
    fn __and__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().and(to_expr(other)?) })
    }

    /// Logical or of two boolean expressions.
    fn __or__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone().or(to_expr(other)?) })
    }

    fn __gt__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        self.gt(other)
    }

    fn __ge__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        self.ge(other)
    }

    fn __lt__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        self.lt(other)
    }

    fn __le__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        self.le(other)
    }

    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone() + to_expr(other)? })
    }

    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: to_expr(other)? + self.inner.clone() })
    }

    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone() - to_expr(other)? })
    }

    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: to_expr(other)? - self.inner.clone() })
    }

    fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone() * to_expr(other)? })
    }

    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: to_expr(other)? * self.inner.clone() })
    }

    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: self.inner.clone() / to_expr(other)? })
    }

    fn __rtruediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: to_expr(other)? / self.inner.clone() })
    }

    fn __neg__(&self) -> Self {
        Self { inner: -self.inner.clone() }
    }

    fn __invert__(&self) -> Self {
        Self { inner: !self.inner.clone() }
    }

    fn __repr__(&self) -> String {
        format!("Expr({})", self.inner)
    }
}

/// References a column by name.
#[pyfunction]
pub(crate) fn col(name: &str) -> PyExpr {
    PyExpr { inner: lazy::col(name) }
}

/// A constant broadcast to every row (int, float or bool).
#[pyfunction]
pub(crate) fn lit(value: &Bound<'_, PyAny>) -> PyResult<PyExpr> {
    Ok(PyExpr { inner: to_expr(value)? })
}

/// Invariant mass of two particles, each given as `(pt, eta, phi, m)`.
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub(crate) fn invariant_mass(
    pt1: &Bound<'_, PyAny>,
    eta1: &Bound<'_, PyAny>,
    phi1: &Bound<'_, PyAny>,
    m1: &Bound<'_, PyAny>,
    pt2: &Bound<'_, PyAny>,
    eta2: &Bound<'_, PyAny>,
    phi2: &Bound<'_, PyAny>,
    m2: &Bound<'_, PyAny>,
) -> PyResult<PyExpr> {
    Ok(PyExpr {
        inner: lazy::invariant_mass(
            to_expr(pt1)?,
            to_expr(eta1)?,
            to_expr(phi1)?,
            to_expr(m1)?,
            to_expr(pt2)?,
            to_expr(eta2)?,
            to_expr(phi2)?,
            to_expr(m2)?,
        ),
    })
}

/// Angular distance `sqrt(d_eta^2 + d_phi^2)` with `d_phi` wrapped to [-pi, pi).
#[pyfunction]
pub(crate) fn delta_r(
    eta1: &Bound<'_, PyAny>,
    phi1: &Bound<'_, PyAny>,
    eta2: &Bound<'_, PyAny>,
    phi2: &Bound<'_, PyAny>,
) -> PyResult<PyExpr> {
    Ok(PyExpr {
        inner: lazy::delta_r(to_expr(eta1)?, to_expr(phi1)?, to_expr(eta2)?, to_expr(phi2)?),
    })
}

/// Transverse momentum `sqrt(px^2 + py^2)`.
#[pyfunction]
pub(crate) fn transverse_momentum(
    px: &Bound<'_, PyAny>,
    py: &Bound<'_, PyAny>,
) -> PyResult<PyExpr> {
    Ok(PyExpr { inner: lazy::transverse_momentum(to_expr(px)?, to_expr(py)?) })
}

/// An aggregate for `group_by(...).agg(...)`. Build with the static methods.
#[pyclass(name = "Agg", module = "tpt_dynamics", from_py_object)]
#[derive(Clone)]
pub struct PyAgg {
    inner: AggExpr,
}

#[pymethods]
impl PyAgg {
    /// Number of rows in each group.
    #[staticmethod]
    fn count() -> Self {
        Self { inner: AggExpr::count() }
    }

    /// Sum of an expression within each group.
    #[staticmethod]
    fn sum(expr: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: AggExpr::sum(to_expr(expr)?) })
    }

    /// Mean of an expression within each group (float64).
    #[staticmethod]
    fn mean(expr: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: AggExpr::mean(to_expr(expr)?) })
    }

    /// Minimum of an expression within each group.
    #[staticmethod]
    fn min(expr: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: AggExpr::min(to_expr(expr)?) })
    }

    /// Maximum of an expression within each group.
    #[staticmethod]
    fn max(expr: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self { inner: AggExpr::max(to_expr(expr)?) })
    }

    /// Renames the output column.
    fn alias(&self, name: &str) -> Self {
        Self { inner: self.inner.clone().alias(name) }
    }
}

/// A deferred pipeline over a data source. Nothing runs until `collect()`.
#[pyclass(name = "LazyFrame", module = "tpt_dynamics")]
#[derive(Clone)]
pub struct PyLazyFrame {
    inner: LazyFrame,
}

fn frame_from_dict(columns: &Bound<'_, PyDict>) -> PyResult<DataFrame> {
    let mut cols = Vec::with_capacity(columns.len());
    for (key, value) in columns.iter() {
        let name: String = key.extract()?;
        let column = column_from_any(&name, &value)?;
        cols.push((name, column));
    }
    DataFrame::new(cols).map_err(value_err)
}

fn parse_dtype(name: &str) -> PyResult<DType> {
    match name {
        "f64" | "float64" => Ok(DType::F64),
        "i64" | "int64" => Ok(DType::I64),
        "bool" => Ok(DType::Bool),
        other => Err(PyValueError::new_err(format!(
            "unknown dtype `{other}` (use \"f64\", \"i64\" or \"bool\")"
        ))),
    }
}

fn push_field(
    column: &mut Column,
    dtype: DType,
    record: &CsvRecord,
    i: usize,
) -> Result<(), String> {
    match (column, dtype) {
        (Column::F64(v), DType::F64) => v.push(record.f64(i).map_err(|e| e.to_string())?),
        (Column::I64(v), DType::I64) => v.push(record.i64(i).map_err(|e| e.to_string())?),
        (Column::Bool(v), DType::Bool) => {
            let raw = record.get(i).unwrap_or("").trim();
            if raw.eq_ignore_ascii_case("true") {
                v.push(true);
            } else if raw.eq_ignore_ascii_case("false") {
                v.push(false);
            } else {
                return Err(format!("line {}: `{raw}` is not a bool", record.line()));
            }
        }
        _ => return Err("internal column type mismatch".to_owned()),
    }
    Ok(())
}

fn read_csv_frame(path: &str, dtypes: &[(String, DType)], delimiter: u8) -> PyResult<DataFrame> {
    let file = File::open(Path::new(path))
        .map_err(|e| PyOSError::new_err(format!("cannot open `{path}`: {e}")))?;
    let opts = CsvOptions { delimiter, ..CsvOptions::default() };
    let mut reader = CsvReader::new(BufReader::new(file), opts).map_err(value_err)?;
    let header: Vec<String> = reader
        .header()
        .ok_or_else(|| PyValueError::new_err("CSV file has no header line"))?
        .to_vec();
    for (name, _) in dtypes {
        if !header.contains(name) {
            return Err(PyValueError::new_err(format!("dtype given for unknown column `{name}`")));
        }
    }
    let types: Vec<DType> = header
        .iter()
        .map(|name| dtypes.iter().find(|(n, _)| n == name).map_or(DType::F64, |(_, d)| *d))
        .collect();
    let mut columns: Vec<Column> = types.iter().map(|d| Column::empty(*d)).collect();
    for record in &mut reader {
        let record = record.map_err(value_err)?;
        if record.len() != header.len() {
            return Err(value_err(format!(
                "line {}: expected {} fields, found {}",
                record.line(),
                header.len(),
                record.len()
            )));
        }
        for (i, dtype) in types.iter().enumerate() {
            push_field(&mut columns[i], *dtype, &record, i).map_err(value_err)?;
        }
    }
    DataFrame::new(header.into_iter().zip(columns)).map_err(value_err)
}

#[pymethods]
impl PyLazyFrame {
    /// Scans an in-memory mapping of column name to 1D numpy array
    /// (float64, int64 or bool). All arrays must have the same length.
    #[staticmethod]
    fn scan_dict(columns: &Bound<'_, PyDict>) -> PyResult<Self> {
        Ok(Self { inner: LazyFrame::scan(frame_from_dict(columns)?) })
    }

    /// Scans a CSV file with a header line. The whole file is parsed when the
    /// frame is created (not streamed). Columns not listed in `dtypes` are
    /// float64; supported dtypes are "f64", "i64" and "bool".
    #[staticmethod]
    #[pyo3(signature = (path, dtypes=None, delimiter=","))]
    fn scan_csv(path: &str, dtypes: Option<&Bound<'_, PyDict>>, delimiter: &str) -> PyResult<Self> {
        let mut parsed = Vec::new();
        if let Some(d) = dtypes {
            for (k, v) in d.iter() {
                parsed.push((k.extract::<String>()?, parse_dtype(&v.extract::<String>()?)?));
            }
        }
        let bytes = delimiter.as_bytes();
        if bytes.len() != 1 {
            return Err(PyValueError::new_err("delimiter must be a single ASCII character"));
        }
        let df = read_csv_frame(path, &parsed, bytes[0])?;
        Ok(Self { inner: LazyFrame::scan(df) })
    }

    /// Keeps the rows where `predicate` is true.
    fn filter(&self, predicate: &PyExpr) -> Self {
        Self { inner: self.inner.clone().filter(predicate.inner.clone()) }
    }

    /// Keeps only the listed expressions (each must have an output name).
    fn select(&self, exprs: Vec<PyExpr>) -> Self {
        Self { inner: self.inner.clone().select(exprs.into_iter().map(|e| e.inner)) }
    }

    /// Adds or replaces a column computed from an expression.
    fn with_column(&self, expr: &PyExpr) -> Self {
        Self { inner: self.inner.clone().with_column(expr.inner.clone()) }
    }

    /// Sets the number of rows per processing chunk.
    fn with_chunk_size(&self, chunk_size: usize) -> Self {
        Self { inner: self.inner.clone().with_chunk_size(chunk_size) }
    }

    /// Starts a group-by on the named key columns.
    fn group_by(&self, keys: Vec<String>) -> PyGroupBy {
        PyGroupBy { frame: self.inner.clone(), keys }
    }

    /// Returns the optimised plan as text.
    fn explain(&self) -> String {
        self.inner.explain()
    }

    /// Executes the pipeline and returns a dict of numpy arrays.
    fn collect<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let out = self.inner.clone().collect().map_err(value_err)?;
        let names: Vec<String> = out.names().to_vec();
        let cols = names.iter().filter_map(|n| out.column(n));
        columns_to_dict(py, &names, cols)
    }
}

/// The pending group-by of a `LazyFrame`; call `agg` to get a frame back.
#[pyclass(name = "GroupBy", module = "tpt_dynamics")]
pub struct PyGroupBy {
    frame: LazyFrame,
    keys: Vec<String>,
}

#[pymethods]
impl PyGroupBy {
    /// Applies aggregates per group and returns the result frame.
    fn agg(&self, aggs: &Bound<'_, PyList>) -> PyResult<PyLazyFrame> {
        let mut list = Vec::with_capacity(aggs.len());
        for item in aggs.iter() {
            list.push(item.extract::<PyAgg>()?.inner);
        }
        let inner = self.frame.clone().group_by(self.keys.clone()).agg(list);
        Ok(PyLazyFrame { inner })
    }
}
