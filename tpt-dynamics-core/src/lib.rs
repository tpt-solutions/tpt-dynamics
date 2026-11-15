//! Foundation crate for tpt-dynamics: physical constants, unit conversions,
//! typed quantities and a small 3-vector type.
//!
//! This crate is `no_std` (with `alloc` not required). It has no FFI and no
//! `unsafe` code, and depends only on [`libm`] for transcendental functions.
//!
//! # Example
//!
//! ```
//! use tpt_dynamics_core::units::{ev_to_kj_per_mol, kcal_per_mol_to_kj_per_mol};
//!
//! // 1 eV per particle is about 96.485 kJ/mol.
//! let kj = ev_to_kj_per_mol(1.0);
//! assert!((kj - 96.485_332_12).abs() < 1e-6);
//!
//! assert!((kcal_per_mol_to_kj_per_mol(1.0) - 4.184).abs() < 1e-12);
//! ```
#![no_std]
#![forbid(unsafe_code)]
#![cfg_attr(docsrs, feature(doc_cfg))]

#[cfg(feature = "std")]
extern crate std;

pub mod constants;
pub mod quantity;
pub mod units;
pub mod vec3;

pub use vec3::Vec3;
