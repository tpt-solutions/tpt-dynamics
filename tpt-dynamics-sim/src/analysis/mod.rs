//! Structural analysis of trajectories: RMSD after optimal superposition,
//! radius of gyration, geometric hydrogen-bond counts and the radial
//! distribution function.
//!
//! Every analysis works on [`Snapshot`]s, one per frame. A trajectory is any
//! iterator of snapshots, so frames can be streamed from a file one at a time
//! and only per-frame results are kept. The series functions take
//! `IntoIterator<Item = Snapshot>` for that reason.
//!
//! - [`fn@rmsd`] / [`rmsd_series`]: root-mean-square deviation after least-squares
//!   superposition (Horn's quaternion method, eigenvalues by Jacobi rotations).
//! - [`radius_of_gyration`] / [`rg_series`]: optionally mass-weighted.
//! - [`hbond`]: geometric hydrogen-bond counts from donor–hydrogen and acceptor
//!   index lists.
//! - [`fn@rdf`]: radial distribution function with minimum-image distances in an
//!   orthorhombic box, normalised by the ideal-gas pair count.
//!
//! Positions are used exactly as given. Molecules that straddle the periodic
//! boundary must be made whole before calling [`radius_of_gyration`] or
//! [`fn@rmsd`]; only the distance-based analyses apply the minimum-image convention.

use std::fmt;

use tpt_dynamics_core::Vec3;

use crate::pbc::SimBox;

pub mod hbond;
pub mod rdf;
pub mod rg;
pub mod rmsd;

pub use hbond::{count_hbonds, hbond_series, Donor, HBondCriteria};
pub use rdf::{rdf, Rdf};
pub use rg::{radius_of_gyration, rg_series};
pub use rmsd::{rmsd, rmsd_series};

/// Positions of one frame together with the box they were recorded in.
#[derive(Clone, Debug, PartialEq)]
pub struct Snapshot {
    /// Positions in nm, one per atom.
    pub positions: Vec<Vec3>,
    /// Periodic box of this frame.
    pub sim_box: SimBox,
}

/// Errors from the analysis functions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AnalysisError {
    /// Two inputs that must have the same number of atoms do not.
    AtomCountMismatch {
        /// Expected atom count.
        expected: usize,
        /// Atom count found.
        found: usize,
    },
    /// An input that must be non-empty is empty, or a weight sum is zero.
    EmptyInput(&'static str),
    /// A parameter is out of range or not finite.
    InvalidParameter(&'static str),
    /// An atom index is not below the number of atoms.
    IndexOutOfRange {
        /// Offending index.
        index: usize,
        /// Atoms in the frame.
        n_atoms: usize,
    },
}

impl fmt::Display for AnalysisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::AtomCountMismatch { expected, found } => {
                write!(f, "atom count mismatch: expected {expected}, found {found}")
            }
            Self::EmptyInput(what) => write!(f, "empty input: {what}"),
            Self::InvalidParameter(what) => write!(f, "invalid parameter: {what}"),
            Self::IndexOutOfRange { index, n_atoms } => {
                write!(f, "atom index {index} out of range for {n_atoms} atoms")
            }
        }
    }
}

impl std::error::Error for AnalysisError {}

/// Checks that an optional per-atom weight slice matches `n` atoms.
fn check_weights(weights: Option<&[f64]>, n: usize) -> Result<(), AnalysisError> {
    if let Some(w) = weights {
        if w.len() != n {
            return Err(AnalysisError::AtomCountMismatch { expected: n, found: w.len() });
        }
        if w.iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err(AnalysisError::InvalidParameter("weights must be finite and non-negative"));
        }
    }
    Ok(())
}
