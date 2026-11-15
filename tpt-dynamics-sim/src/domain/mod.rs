//! Spatial domain decomposition for shared-memory parallel MD.
//!
//! - [`grid`]: the separable sub-domain grid, owner lookup, quantile rebalancing.
//! - [`engine`]: per-domain workers with halo (ghost) exchange and atom
//!   migration over `std::sync::mpsc` channels, and the NVE driver.
//!
//! The decomposed engine supports orthorhombic boxes and Lennard-Jones only.
//! See [`engine`] for the execution model and determinism guarantees.

pub mod engine;
pub mod grid;

pub use engine::{decomposed_forces, run_decomposed_nve, DecompOptions, DecompSummary};
pub use grid::{auto_dims, imbalance, DomainGrid};
