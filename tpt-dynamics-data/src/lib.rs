//! Data analysis for tpt-dynamics (the ROOT replacement).
//!
//! * [`hist1d`], [`hist2d`], [`profile`]: 1D, 2D and profile histograms with
//!   fixed or variable binning, under/overflow, weights, and associative,
//!   commutative `merge` for distributed use.
//! * [`serialise`]: binary and JSON round-trips for all histogram types. The
//!   binary form is checksummed.
//! * [`csv`]: a streaming CSV reader with bounded memory and line-numbered errors.
//! * [`columnar`]: a versioned, checksummed binary columnar format for detector
//!   data.
//! * [`physics`]: invariant mass, transverse momentum and ΔR helpers.
//! * [`json`], [`crc32`], [`rng`]: dependency-free building blocks used by the
//!   modules above.
//!
//! * [`parquet`] (feature `parquet`): read-only Parquet import, clean-room, with no
//!   Arrow or Parquet crates (ADR-008). Writing Parquet is not planned. See `todo.md` for the plan.
#![forbid(unsafe_code)]

pub mod arrow_ipc;
pub mod columnar;
pub mod crc32;
pub mod csv;
pub mod hist1d;
pub mod hist2d;
pub mod jagged;
pub mod json;
pub mod lazy;
#[cfg(feature = "parquet")]
pub mod parquet;
pub mod physics;
pub mod profile;
pub mod rng;
pub mod serialise;
pub mod trajectory;

pub use hist1d::{Histogram1D, HistogramError, Slot};
pub use hist2d::Histogram2D;
pub use profile::{Profile1D, Profile2D, ProfileCell};
