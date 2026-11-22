//! Classical molecular dynamics engine for tpt-dynamics (the GROMACS replacement).
//!
//! Contents (Phases 2 and 5 of `todo.md`):
//!
//! - [`system`]: structure-of-arrays particle state, masses, types, topology.
//! - [`pbc`]: orthorhombic and triclinic periodic boxes (GROMACS lower-triangular form),
//!   minimum-image convention. GRO and PDB files with triclinic boxes are read and written.
//! - [`pdb`] and [`gro`]: fixed-width structure readers and writers.
//! - [`config`]: TOML run configuration: integrator, thermostat, barostat, Lennard-Jones or
//!   bonded (ffparams) force field with optional PME, SHAKE/LINCS/SETTLE constraints,
//!   minimisation before the run, and checkpoint write/restart. Unknown keys are rejected.
//! - [`topology`] and [`ffparams`]: bond graph to angles, dihedrals, impropers,
//!   1-4 pairs and exclusions; TOML force-field parameter files. The shipped
//!   parameter file is synthetic (see `data/synthetic-test.toml`).
//! - [`bonded`]: harmonic bonds and angles (with Urey-Bradley), periodic dihedrals, impropers.
//! - [`nonbonded`]: Lennard-Jones with exclusions and scaled 1-4 pairs, PME Coulomb.
//! - [`neighbour`]: cell lists, Verlet lists with skin, parallel (rayon) builds.
//! - [`forcefield`]: the [`forcefield::ForceField`] trait and Lennard-Jones (cutoff, shift, switch).
//! - [`integrator`]: [`integrator::Integrator`] trait with velocity Verlet and leap-frog.
//! - [`constraints`]: SHAKE, LINCS-style matrix expansion, SETTLE (rigid water).
//! - [`thermostat`]: Berendsen, stochastic v-rescale, Nosé–Hoover (chain length 1).
//! - [`barostat`]: isotropic Berendsen, C-rescale and Parrinello–Rahman pressure coupling.
//! - [`minimize`]: adaptive steepest descent and L-BFGS.
//! - [`solvate`]: placeholder water grid and ion placement.
//! - [`checkpoint`]: versioned, CRC-32-protected binary checkpoints (bitwise restart).
//! - [`crc32`] and [`rng`]: checksum and deterministic PRNG used by the above.
//! - [`output`]: observer hook, binary trajectory format, CSV energy log.
//! - [`simulation`]: the run loop tying everything together.
//!
//! Units are GROMACS-style (nm, ps, amu, kJ/mol, K). Reduced Lennard-Jones
//! units (σ = ε = m = 1) work equally well.
//!
//! # Example
//!
//! ```
//! use tpt_dynamics_sim::config::SimConfig;
//! use tpt_dynamics_sim::system::System;
//! use tpt_dynamics_sim::pbc::SimBox;
//!
//! let cfg = SimConfig::from_toml_str("nsteps = 10\ntemperature = 0.7\n").unwrap();
//! let sys = System::new(SimBox::new(3.0, 3.0, 3.0).unwrap());
//! // An empty system cannot be run; a real system is built from `pdb` or `gro`
//! // input or by pushing atoms with `System::push_atom`.
//! assert!(cfg.build_simulation(sys).is_err());
//! ```
#![forbid(unsafe_code)]
// Physics code uses short single-letter names (x, y, z, r, v, m) by convention.
#![allow(clippy::many_single_char_names)]
// Physics code: symbols such as r_ij, r_kj, d_oh and d_hh are deliberately similar,
// and exact float comparisons are intended (zero-tilt detection, bitwise checks).
#![allow(clippy::similar_names, clippy::float_cmp)]

pub mod analysis;
pub mod barostat;
pub mod bonded;
pub mod checkpoint;
pub mod config;
pub mod constraints;
/// CRC-32, re-exported from `tpt-dynamics-data`, which holds the only implementation.
pub mod crc32 {
    pub use tpt_dynamics_data::crc32::{checksum, Crc32};

    /// CRC-32 of `data`. Same as [`checksum`]; kept for source compatibility.
    pub fn crc32(data: &[u8]) -> u32 {
        checksum(data)
    }
}
pub mod domain;
pub mod error;
pub mod ffparams;
pub mod forcefield;
pub mod gro;
pub mod integrator;
pub mod minimize;
pub mod neighbour;
pub mod nonbonded;
pub mod output;
pub mod pbc;
pub mod pdb;
pub mod rng;
pub mod simulation;
pub mod solvate;
pub mod system;
pub mod textio;
pub mod thermostat;
pub mod topology;
pub mod water;

pub use error::SimError;

/// Returns the crate version.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
