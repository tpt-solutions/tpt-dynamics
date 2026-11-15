//! Python bindings for tpt-dynamics (`import tpt_dynamics`).
//!
//! Exposes three areas of the Rust workspace:
//!
//! * [`frame`]: lazy data frames (scan from numpy arrays or CSV, filter,
//!   select, `with_column`, `group_by` + `agg`, `collect` to numpy arrays).
//! * [`hist`]: 1D and 2D physics histograms with numpy fill, merge and dict export.
//! * [`sim`]: molecular dynamics runs configured from a TOML string.
//!
//! The crate is a thin layer: all physics and data logic lives in
//! `tpt-dynamics-data` and `tpt-dynamics-sim`. Hand-written code here contains
//! no `unsafe`. The generated `PyO3` and rust-numpy glue is covered by ADR-007.
#![deny(unsafe_code)]

mod convert;
mod errors;
mod frame;
mod hist;
mod sim;

use pyo3::prelude::*;

/// Returns the crate version.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
#[pyo3(name = "_tpt_dynamics")]
fn tpt_dynamics_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_class::<frame::PyExpr>()?;
    m.add_class::<frame::PyAgg>()?;
    m.add_class::<frame::PyLazyFrame>()?;
    m.add_class::<frame::PyGroupBy>()?;
    m.add_function(wrap_pyfunction!(frame::col, m)?)?;
    m.add_function(wrap_pyfunction!(frame::lit, m)?)?;
    m.add_function(wrap_pyfunction!(frame::invariant_mass, m)?)?;
    m.add_function(wrap_pyfunction!(frame::delta_r, m)?)?;
    m.add_function(wrap_pyfunction!(frame::transverse_momentum, m)?)?;
    m.add_class::<hist::PyAxis>()?;
    m.add_class::<hist::PyHistogram1D>()?;
    m.add_class::<hist::PyHistogram2D>()?;
    m.add_function(wrap_pyfunction!(sim::run_simulation, m)?)?;
    m.add_function(wrap_pyfunction!(sim::validate_config, m)?)?;
    Ok(())
}
