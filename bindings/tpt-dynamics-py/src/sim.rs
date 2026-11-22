//! Molecular dynamics runs configured from a TOML string.

use std::io::Cursor;

use numpy::IntoPyArray;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use tpt_dynamics_sim::config::SimConfig;
use tpt_dynamics_sim::output::Energies;
use tpt_dynamics_sim::{gro, pdb};

use crate::convert::f64_array;
use crate::errors::{sim_err, value_err};

/// Parses and validates a TOML run configuration.
///
/// Raises `ValueError` if the configuration is malformed or physically invalid.
#[pyfunction]
pub(crate) fn validate_config(config: &str) -> PyResult<()> {
    SimConfig::from_toml_str(config).map(|_| ()).map_err(sim_err)
}

/// Runs a simulation and returns the sampled energies as numpy arrays.
///
/// `config` is the TOML run configuration (see the `tpt_dynamics` docs for keys).
/// `structure` is the text of a GRO (`format="gro"`) or PDB (`format="pdb"`)
/// file that gives the initial positions, box and masses. Output-file settings
/// in `config` are ignored; the returned arrays replace them.
///
/// Energies are sampled at step 0 and then every `sample_every` steps. The
/// returned dict has the keys `step`, `time`, `kinetic`, `potential`, `total`
/// and `temperature`, each a 1D numpy array (energies in kJ/mol, time in ps).
///
/// The optional keyword arguments override the matching TOML settings for this
/// call: `minimize` sets `[minimize] enabled`, `checkpoint_restart` sets
/// `[checkpoint] restart` (the run then starts from that state, and the first
/// sampled step is the checkpoint's step), and `checkpoint_write` sets
/// `[checkpoint] write`, so the final state is written to that file.
#[pyfunction]
#[pyo3(signature = (
    config,
    structure,
    *,
    format="gro",
    sample_every=1,
    minimize=None,
    checkpoint_restart=None,
    checkpoint_write=None
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_simulation<'py>(
    py: Python<'py>,
    config: &str,
    structure: &str,
    format: &str,
    sample_every: u64,
    minimize: Option<bool>,
    checkpoint_restart: Option<String>,
    checkpoint_write: Option<String>,
) -> PyResult<Bound<'py, PyDict>> {
    if sample_every == 0 {
        return Err(PyValueError::new_err("sample_every must be at least 1"));
    }
    let mut cfg = SimConfig::from_toml_str(config).map_err(sim_err)?;
    if let Some(on) = minimize {
        cfg.minimize.enabled = on;
    }
    if let Some(path) = checkpoint_restart {
        cfg.checkpoint.restart = Some(path.into());
    }
    if let Some(path) = checkpoint_write {
        cfg.checkpoint.write = Some(path.into());
    }
    cfg.validate().map_err(sim_err)?;
    let reader = Cursor::new(structure.as_bytes());
    let system = match format {
        "gro" => gro::read_gro(reader).map_err(sim_err)?,
        "pdb" => pdb::read_pdb(reader).map_err(sim_err)?,
        other => {
            return Err(value_err(format!(
                "unknown structure format `{other}` (use \"gro\" or \"pdb\")"
            )))
        }
    };
    let mut sim = cfg.build_simulation(system).map_err(sim_err)?;

    let mut samples: Vec<Energies> = vec![sim.energies()];
    for _ in 0..cfg.nsteps {
        sim.advance().map_err(sim_err)?;
        if sim.step % sample_every == 0 {
            samples.push(sim.energies());
        }
    }

    cfg.write_checkpoint(&mut sim).map_err(sim_err)?;

    let out = PyDict::new(py);
    let steps: Vec<u64> = samples.iter().map(|e| e.step).collect();
    out.set_item("step", steps.into_pyarray(py))?;
    out.set_item("time", f64_array(py, samples.iter().map(|e| e.time).collect()))?;
    out.set_item("kinetic", f64_array(py, samples.iter().map(|e| e.kinetic).collect()))?;
    out.set_item("potential", f64_array(py, samples.iter().map(|e| e.potential).collect()))?;
    out.set_item("total", f64_array(py, samples.iter().map(|e| e.total).collect()))?;
    out.set_item("temperature", f64_array(py, samples.iter().map(|e| e.temperature).collect()))?;
    Ok(out)
}
