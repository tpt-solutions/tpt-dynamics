//! Physics histograms exposed to Python.

use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use tpt_dynamics_data::{Histogram1D, Histogram2D};

use crate::convert::{check_len, f64_array, f64_matrix};
use crate::errors::{index_err, value_err};

/// One histogram axis, exposed to Python as `Axis`.
///
/// Holds the bin edges. `edges` and `centers` return numpy arrays, and `len()`
/// gives the number of regular bins.
#[pyclass(name = "Axis", module = "tpt_dynamics", from_py_object)]
#[derive(Clone, Debug)]
pub struct PyAxis {
    edges: Vec<f64>,
}

#[pymethods]
impl PyAxis {
    /// Bin edges as a numpy array of length `nbins + 1`.
    #[getter]
    fn edges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.edges.clone())
    }

    /// Bin centres as a numpy array of length `nbins`.
    #[getter]
    fn centers<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.edges.windows(2).map(|w| f64::midpoint(w[0], w[1])).collect())
    }

    /// Number of regular bins.
    fn __len__(&self) -> usize {
        self.edges.len().saturating_sub(1)
    }

    fn __repr__(&self) -> String {
        format!(
            "Axis(nbins={}, lo={}, hi={})",
            self.__len__(),
            self.edges[0],
            self.edges[self.edges.len() - 1]
        )
    }
}

/// `(values, x_edges, y_edges)` returned by `Histogram2D.to_numpy`.
type Bound3<'py> =
    (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>);

/// Copies a 2D numpy array into a flat row-major vector after a shape check.
fn array2_flat(
    a: &PyReadonlyArray2<'_, f64>,
    nx: usize,
    ny: usize,
    name: &str,
) -> PyResult<Vec<f64>> {
    let arr = a.as_array();
    let (r, c) = arr.dim();
    if (r, c) != (nx, ny) {
        return Err(value_err(format!("{name} has shape ({r}, {c}), expected ({nx}, {ny})")));
    }
    Ok(arr.iter().copied().collect())
}

/// Fills `hist` from a numpy array of values and optional weights.
fn fill_1d(
    hist: &mut Histogram1D,
    x: &PyReadonlyArray1<'_, f64>,
    weights: Option<&PyReadonlyArray1<'_, f64>>,
) -> PyResult<()> {
    let x = x.as_array();
    match weights {
        None => x.iter().for_each(|&v| hist.fill(v)),
        Some(w) => {
            let w = w.as_array();
            check_len("weights", x.len(), w.len())?;
            x.iter().zip(w.iter()).for_each(|(&v, &wt)| hist.fill_weighted(v, wt));
        }
    }
    Ok(())
}

/// A 1D histogram with weighted fills, under/overflow and error tracking.
///
/// Construct with uniform bins `Histogram1D(nbins, lo, hi)` or with explicit
/// edges via `Histogram1D.variable(edges)`.
#[pyclass(name = "Histogram1D", module = "tpt_dynamics", from_py_object)]
#[derive(Clone)]
pub struct PyHistogram1D {
    inner: Histogram1D,
}

#[pymethods]
impl PyHistogram1D {
    /// Creates `nbins` uniform bins over `[lo, hi)`.
    #[new]
    fn new(nbins: usize, lo: f64, hi: f64) -> PyResult<Self> {
        Ok(Self { inner: Histogram1D::uniform(nbins, lo, hi).map_err(value_err)? })
    }

    /// Creates a histogram with explicit, strictly increasing bin edges.
    #[staticmethod]
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts Vec<f64> by value.
    fn variable(edges: Vec<f64>) -> PyResult<Self> {
        Ok(Self { inner: Histogram1D::variable(&edges).map_err(value_err)? })
    }

    /// Fills from a 1D float64 numpy array, with optional per-entry weights.
    #[pyo3(signature = (x, weights=None))]
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts numpy views by value.
    fn fill(
        &mut self,
        x: PyReadonlyArray1<'_, f64>,
        weights: Option<PyReadonlyArray1<'_, f64>>,
    ) -> PyResult<()> {
        fill_1d(&mut self.inner, &x, weights.as_ref())
    }

    /// Number of bins (excluding under/overflow).
    fn nbins(&self) -> usize {
        self.inner.nbins()
    }

    /// Bin edges as a numpy array of length `nbins + 1`.
    fn edges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.inner.edges().to_vec())
    }

    /// Hist-style axis list. A 1D histogram has one axis.
    #[getter]
    fn axes(&self) -> (PyAxis,) {
        (PyAxis { edges: self.inner.edges().to_vec() },)
    }

    /// Bin contents (sum of weights) as a numpy array. Same as `contents()`.
    fn values<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.contents(py)
    }

    /// Per-bin variance, `sum of w^2`, as a numpy array.
    fn variances<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.inner.sum_weights_squared().to_vec())
    }

    /// Returns `(values, edges)` as a tuple of numpy arrays.
    fn to_numpy<'py>(
        &self,
        py: Python<'py>,
    ) -> (Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>) {
        (self.values(py), self.edges(py))
    }

    /// Builds a histogram from bin values and edges.
    ///
    /// `variances` defaults to `values`, which is the Poisson variance of
    /// unweighted counts. Pass it explicitly for weighted histograms.
    #[staticmethod]
    #[pyo3(signature = (values, edges, variances=None))]
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts numpy views by value.
    fn from_numpy(
        values: PyReadonlyArray1<'_, f64>,
        edges: PyReadonlyArray1<'_, f64>,
        variances: Option<PyReadonlyArray1<'_, f64>>,
    ) -> PyResult<Self> {
        let v: Vec<f64> = values.as_array().iter().copied().collect();
        let e: Vec<f64> = edges.as_array().iter().copied().collect();
        let var: Vec<f64> = match variances {
            Some(w) => w.as_array().iter().copied().collect(),
            None => v.clone(),
        };
        Ok(Self { inner: Histogram1D::from_bins(&e, &v, &var).map_err(value_err)? })
    }

    /// Bin centres as a numpy array.
    fn bin_centers<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, (0..self.inner.nbins()).map(|i| self.inner.bin_center(i)).collect())
    }

    /// Weighted bin contents as a numpy array.
    fn contents<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.inner.contents().to_vec())
    }

    /// Statistical errors per bin, `sqrt(sum of w^2)`.
    fn errors<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.inner.sum_weights_squared().iter().map(|w| w.sqrt()).collect())
    }

    /// Content of bin `i`.
    fn bin(&self, i: usize) -> PyResult<f64> {
        self.check_index(i)?;
        Ok(self.inner.bin(i))
    }

    /// Error of bin `i`.
    fn bin_error(&self, i: usize) -> PyResult<f64> {
        self.check_index(i)?;
        Ok(self.inner.bin_error(i))
    }

    /// Weighted sum of entries below the first edge.
    #[getter]
    fn underflow(&self) -> f64 {
        self.inner.underflow()
    }

    /// Weighted sum of entries at or above the last edge.
    #[getter]
    fn overflow(&self) -> f64 {
        self.inner.overflow()
    }

    /// Number of fill calls (all entries, including out-of-range).
    #[getter]
    fn entries(&self) -> u64 {
        self.inner.entries()
    }

    /// Number of NaN values that were filled and skipped.
    #[getter]
    fn nan_fills(&self) -> u64 {
        self.inner.nan_fills()
    }

    /// Sum of in-range bin contents.
    fn integral(&self) -> f64 {
        self.inner.integral()
    }

    /// Weighted mean of in-range fills, or None if empty.
    fn mean(&self) -> Option<f64> {
        self.inner.mean()
    }

    /// Weighted standard deviation of in-range fills, or None if undefined.
    fn std_dev(&self) -> Option<f64> {
        self.inner.std_dev()
    }

    /// Approximate quantile `p` in [0, 1], or None if the histogram is empty.
    fn quantile(&self, p: f64) -> Option<f64> {
        self.inner.quantile(p)
    }

    /// Adds another histogram with identical binning into this one.
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts PyRef by value.
    fn merge(&mut self, other: PyRef<'_, Self>) -> PyResult<()> {
        self.inner.merge(&other.inner).map_err(value_err)
    }

    /// Returns an independent copy.
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Exports every field to a dict of numpy arrays and scalars.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("edges", self.edges(py))?;
        d.set_item("contents", self.contents(py))?;
        d.set_item("errors", self.errors(py))?;
        d.set_item("underflow", self.inner.underflow())?;
        d.set_item("overflow", self.inner.overflow())?;
        d.set_item("entries", self.inner.entries())?;
        d.set_item("nan_fills", self.inner.nan_fills())?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "Histogram1D(nbins={}, entries={}, integral={})",
            self.inner.nbins(),
            self.inner.entries(),
            self.inner.integral()
        )
    }
}

impl PyHistogram1D {
    fn check_index(&self, i: usize) -> PyResult<()> {
        let n = self.inner.nbins();
        if i < n {
            Ok(())
        } else {
            Err(index_err(i, n))
        }
    }
}

/// A 2D histogram with weighted fills.
///
/// Construct with uniform bins `Histogram2D(nx, xlo, xhi, ny, ylo, yhi)` or with
/// explicit edges via `Histogram2D.variable(xedges, yedges)`. Contents are
/// indexed `[ix, iy]`.
#[pyclass(name = "Histogram2D", module = "tpt_dynamics", from_py_object)]
#[derive(Clone)]
pub struct PyHistogram2D {
    inner: Histogram2D,
}

#[pymethods]
impl PyHistogram2D {
    /// Creates `nx * ny` uniform bins over `[xlo, xhi) x [ylo, yhi)`.
    #[new]
    fn new(nx: usize, xlo: f64, xhi: f64, ny: usize, ylo: f64, yhi: f64) -> PyResult<Self> {
        Ok(Self { inner: Histogram2D::uniform(nx, xlo, xhi, ny, ylo, yhi).map_err(value_err)? })
    }

    /// Creates a histogram with explicit edges on each axis.
    #[staticmethod]
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts Vec<f64> by value.
    fn variable(xedges: Vec<f64>, yedges: Vec<f64>) -> PyResult<Self> {
        Ok(Self { inner: Histogram2D::variable(&xedges, &yedges).map_err(value_err)? })
    }

    /// Fills from two 1D float64 numpy arrays of equal length, with optional weights.
    #[pyo3(signature = (x, y, weights=None))]
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts numpy views by value.
    fn fill(
        &mut self,
        x: PyReadonlyArray1<'_, f64>,
        y: PyReadonlyArray1<'_, f64>,
        weights: Option<PyReadonlyArray1<'_, f64>>,
    ) -> PyResult<()> {
        let (xa, ya) = (x.as_array(), y.as_array());
        check_len("y", xa.len(), ya.len())?;
        match weights {
            None => {
                for (&xv, &yv) in xa.iter().zip(ya.iter()) {
                    self.inner.fill(xv, yv);
                }
            }
            Some(w) => {
                let wa = w.as_array();
                check_len("weights", xa.len(), wa.len())?;
                for ((&xv, &yv), &wt) in xa.iter().zip(ya.iter()).zip(wa.iter()) {
                    self.inner.fill_weighted(xv, yv, wt);
                }
            }
        }
        Ok(())
    }

    /// Number of bins along x.
    fn nbins_x(&self) -> usize {
        self.inner.nbins_x()
    }

    /// Number of bins along y.
    fn nbins_y(&self) -> usize {
        self.inner.nbins_y()
    }

    /// Weighted content of bin `(ix, iy)`.
    fn bin(&self, ix: usize, iy: usize) -> PyResult<f64> {
        self.check_index(ix, iy)?;
        Ok(self.inner.bin(ix, iy))
    }

    /// Error of bin `(ix, iy)`.
    fn bin_error(&self, ix: usize, iy: usize) -> PyResult<f64> {
        self.check_index(ix, iy)?;
        Ok(self.inner.bin_error(ix, iy))
    }

    /// Bin contents as a 2D numpy array of shape `(nbins_x, nbins_y)`.
    fn contents<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
        let (nx, ny) = (self.inner.nbins_x(), self.inner.nbins_y());
        let flat = (0..nx)
            .flat_map(|ix| (0..ny).map(move |iy| (ix, iy)))
            .map(|(ix, iy)| self.inner.bin(ix, iy))
            .collect();
        f64_matrix(py, flat, nx, ny)
    }

    /// Bin errors as a 2D numpy array of shape `(nbins_x, nbins_y)`.
    fn errors<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
        let (nx, ny) = (self.inner.nbins_x(), self.inner.nbins_y());
        let flat = (0..nx)
            .flat_map(|ix| (0..ny).map(move |iy| (ix, iy)))
            .map(|(ix, iy)| self.inner.bin_error(ix, iy))
            .collect();
        f64_matrix(py, flat, nx, ny)
    }

    /// Hist-style axis list: `(x_axis, y_axis)`.
    #[getter]
    fn axes(&self) -> (PyAxis, PyAxis) {
        (
            PyAxis { edges: self.inner.x_edges().to_vec() },
            PyAxis { edges: self.inner.y_edges().to_vec() },
        )
    }

    /// Bin contents (sum of weights), shape `(nbins_x, nbins_y)`. Same as `contents()`.
    fn values<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        self.contents(py)
    }

    /// Per-bin variance, `sum of w^2`, shape `(nbins_x, nbins_y)`.
    fn variances<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let (nx, ny) = (self.inner.nbins_x(), self.inner.nbins_y());
        let flat = (0..nx)
            .flat_map(|ix| (0..ny).map(move |iy| (ix, iy)))
            .map(|(ix, iy)| self.inner.bin_variance(ix, iy))
            .collect();
        f64_matrix(py, flat, nx, ny)
    }

    /// Returns `(values, x_edges, y_edges)` as numpy arrays.
    fn to_numpy<'py>(&self, py: Python<'py>) -> PyResult<Bound3<'py>> {
        Ok((self.values(py)?, self.x_edges(py), self.y_edges(py)))
    }

    /// Builds a histogram from a `(nbins_x, nbins_y)` array and the two edge arrays.
    ///
    /// `variances` defaults to `values` (Poisson variance of unweighted counts).
    #[staticmethod]
    #[pyo3(signature = (values, xedges, yedges, variances=None))]
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts numpy views by value.
    fn from_numpy(
        values: PyReadonlyArray2<'_, f64>,
        xedges: PyReadonlyArray1<'_, f64>,
        yedges: PyReadonlyArray1<'_, f64>,
        variances: Option<PyReadonlyArray2<'_, f64>>,
    ) -> PyResult<Self> {
        let xe: Vec<f64> = xedges.as_array().iter().copied().collect();
        let ye: Vec<f64> = yedges.as_array().iter().copied().collect();
        let (nx, ny) = (xe.len().saturating_sub(1), ye.len().saturating_sub(1));
        let v = array2_flat(&values, nx, ny, "values")?;
        let var = match &variances {
            Some(w) => array2_flat(w, nx, ny, "variances")?,
            None => v.clone(),
        };
        Ok(Self { inner: Histogram2D::from_bins(&xe, &ye, &v, &var).map_err(value_err)? })
    }

    /// Bin edges along x.
    fn x_edges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.inner.x_edges().to_vec())
    }

    /// Bin edges along y.
    fn y_edges<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        f64_array(py, self.inner.y_edges().to_vec())
    }

    /// Sum of in-range bin contents.
    fn integral(&self) -> f64 {
        self.inner.integral()
    }

    /// Number of fill calls.
    #[getter]
    fn entries(&self) -> u64 {
        self.inner.entries()
    }

    /// Number of fills skipped because a coordinate was NaN.
    #[getter]
    fn nan_fills(&self) -> u64 {
        self.inner.nan_fills()
    }

    /// Projection onto the x axis (sums over y).
    fn project_x(&self) -> PyHistogram1D {
        PyHistogram1D { inner: self.inner.project_x() }
    }

    /// Projection onto the y axis (sums over x).
    fn project_y(&self) -> PyHistogram1D {
        PyHistogram1D { inner: self.inner.project_y() }
    }

    /// Adds another histogram with identical binning into this one.
    #[allow(clippy::needless_pass_by_value)] // PyO3 extracts PyRef by value.
    fn merge(&mut self, other: PyRef<'_, Self>) -> PyResult<()> {
        self.inner.merge(&other.inner).map_err(value_err)
    }

    /// Returns an independent copy.
    fn copy(&self) -> Self {
        self.clone()
    }

    /// Exports the histogram to a dict of numpy arrays and scalars.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("x_edges", self.x_edges(py))?;
        d.set_item("y_edges", self.y_edges(py))?;
        d.set_item("contents", self.contents(py)?)?;
        d.set_item("errors", self.errors(py)?)?;
        d.set_item("entries", self.inner.entries())?;
        d.set_item("nan_fills", self.inner.nan_fills())?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "Histogram2D(nx={}, ny={}, entries={}, integral={})",
            self.inner.nbins_x(),
            self.inner.nbins_y(),
            self.inner.entries(),
            self.inner.integral()
        )
    }
}

impl PyHistogram2D {
    fn check_index(&self, ix: usize, iy: usize) -> PyResult<()> {
        if ix >= self.inner.nbins_x() {
            return Err(index_err(ix, self.inner.nbins_x()));
        }
        if iy >= self.inner.nbins_y() {
            return Err(index_err(iy, self.inner.nbins_y()));
        }
        Ok(())
    }
}
