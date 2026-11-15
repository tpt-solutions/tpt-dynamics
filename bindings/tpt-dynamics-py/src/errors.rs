//! Mapping from workspace error types to Python exceptions.

use pyo3::exceptions::{PyIndexError, PyOSError, PyValueError};
use pyo3::PyErr;
use std::fmt::Display;

/// Maps any displayable workspace error to `ValueError`.
pub(crate) fn value_err<E: Display>(e: E) -> PyErr {
    PyValueError::new_err(e.to_string())
}

/// Maps a simulation error: I/O failures become `OSError`, the rest `ValueError`.
pub(crate) fn sim_err(e: tpt_dynamics_sim::SimError) -> PyErr {
    match e {
        tpt_dynamics_sim::SimError::Io(io) => PyOSError::new_err(io.to_string()),
        other => PyValueError::new_err(other.to_string()),
    }
}

/// Maps an out-of-range bin index to `IndexError`.
pub(crate) fn index_err(index: usize, len: usize) -> PyErr {
    PyIndexError::new_err(format!("bin index {index} out of range for {len} bins"))
}
