//! SIMD and mesh-based acceleration kernels for tpt-dynamics.
//!
//! This crate provides force kernels that operate on plain structure-of-arrays
//! slices, so `tpt-dynamics-sim` can call them without this crate depending on
//! simulation types.
//!
//! # Contents
//!
//! * [`pairwise`]: scalar Lennard-Jones, real-space Ewald Coulomb (`erfc(alpha r)/r`)
//!   and excluded-pair formulas. These define the reference behaviour.
//! * [`scalar`]: scalar pair-list kernel over `(i, j)` lists (reference path).
//! * [`simd_kernel`]: the same kernel vectorised 4-wide with `tpt-simd` (ADR-003).
//! * [`listbuild`]: SIMD cell-list pair search (plain and Verlet lists).
//! * [`cluster`]: padded 4-atom cluster layout, a bounding-sphere cluster-pair
//!   builder and a SIMD cluster-pair kernel.
//! * [`fft`]: pure-Rust radix-2 complex FFT (1D and 3D, rayon-parallel).
//! * [`pme`]: smooth particle-mesh Ewald reciprocal solver and parameter selection.
//! * [`reference`]: brute-force Ewald sums used to validate PME.
//!
//! Units: the Coulomb constant is 1 and the caller chooses length/energy units.
//! Forces are *added* into [`ForceBuffers`], so callers zero them first.
//!
//! The crate contains no `unsafe` code.

#![forbid(unsafe_code)]
// Numerical kernels index several arrays with one loop variable and use short
// physics symbols (r, k, p, u). These pedantic lints would obscure the maths.
#![allow(
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::needless_range_loop,
    clippy::manual_midpoint,
    clippy::cast_possible_wrap
)]

pub mod cluster;
pub mod fft;
pub mod listbuild;
pub mod pairwise;
pub mod pme;
pub mod reference;
pub mod scalar;
pub mod simd_kernel;
pub mod types;

pub use types::{AtomParams, Energies, ForceBuffers, NonbondedConfig, Positions};

/// Returns the crate version.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
