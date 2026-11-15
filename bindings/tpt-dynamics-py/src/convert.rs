//! Conversions between numpy arrays and lazy-frame columns.

use numpy::{IntoPyArray, PyArray1, PyArrayMethods, PyReadonlyArray1};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use tpt_dynamics_data::lazy::Column;

use crate::errors::value_err;

/// Copies a 1D numpy array (float64, int64 or bool) into a [`Column`].
pub(crate) fn column_from_any(name: &str, value: &Bound<'_, PyAny>) -> PyResult<Column> {
    if let Ok(a) = value.extract::<PyReadonlyArray1<'_, f64>>() {
        return Ok(Column::F64(a.as_array().iter().copied().collect()));
    }
    if let Ok(a) = value.extract::<PyReadonlyArray1<'_, i64>>() {
        return Ok(Column::I64(a.as_array().iter().copied().collect()));
    }
    if let Ok(a) = value.extract::<PyReadonlyArray1<'_, bool>>() {
        return Ok(Column::Bool(a.as_array().iter().copied().collect()));
    }
    Err(PyTypeError::new_err(format!(
        "column `{name}` must be a 1D numpy array of float64, int64 or bool"
    )))
}

/// Builds a dict of numpy arrays from named columns.
pub(crate) fn columns_to_dict<'py, 'a>(
    py: Python<'py>,
    names: &[String],
    columns: impl Iterator<Item = &'a Column>,
) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);
    for (name, col) in names.iter().zip(columns) {
        let arr: Bound<'py, PyAny> = match col {
            Column::F64(v) => v.clone().into_pyarray(py).into_any(),
            Column::I64(v) => v.clone().into_pyarray(py).into_any(),
            Column::Bool(v) => v.clone().into_pyarray(py).into_any(),
        };
        out.set_item(name, arr)?;
    }
    Ok(out)
}

/// Converts an `f64` vector to a 1D numpy array.
pub(crate) fn f64_array(py: Python<'_>, v: Vec<f64>) -> Bound<'_, PyArray1<f64>> {
    v.into_pyarray(py)
}

/// Converts a flat row-major `f64` vector to a 2D numpy array of shape `(nx, ny)`.
pub(crate) fn f64_matrix(
    py: Python<'_>,
    flat: Vec<f64>,
    nx: usize,
    ny: usize,
) -> PyResult<Bound<'_, numpy::PyArray2<f64>>> {
    f64_array(py, flat).reshape([nx, ny]).map_err(value_err)
}

/// Checks that an optional weights array matches the length of `x`.
pub(crate) fn check_len(name: &str, expected: usize, got: usize) -> PyResult<()> {
    if expected == got {
        Ok(())
    } else {
        Err(value_err(format!("{name} has length {got}, expected {expected}")))
    }
}
