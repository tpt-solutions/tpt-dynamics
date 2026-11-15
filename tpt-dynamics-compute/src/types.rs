//! Shared structure-of-arrays (`SoA`) views used by every kernel in this crate.
//!
//! All kernels take plain slices so that `tpt-dynamics-sim` can call them
//! without this crate depending on the simulation data structures.

/// Positions of `n` atoms in separate x, y and z slices of equal length.
#[derive(Clone, Copy, Debug)]
pub struct Positions<'a> {
    /// x coordinates.
    pub x: &'a [f64],
    /// y coordinates.
    pub y: &'a [f64],
    /// z coordinates.
    pub z: &'a [f64],
}

impl<'a> Positions<'a> {
    /// Creates a view. Panics if the three slices have different lengths.
    pub fn new(x: &'a [f64], y: &'a [f64], z: &'a [f64]) -> Self {
        assert!(x.len() == y.len() && y.len() == z.len(), "coordinate slices differ in length");
        Self { x, y, z }
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.x.len()
    }

    /// True when there are no atoms.
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

/// Per-atom force-field parameters in `SoA` form.
#[derive(Clone, Copy, Debug)]
pub struct AtomParams<'a> {
    /// Lennard-Jones sigma per atom (combined with Lorentz-Berthelot rules).
    pub sigma: &'a [f64],
    /// Lennard-Jones epsilon per atom (combined as `sqrt(e_i * e_j)`).
    pub epsilon: &'a [f64],
    /// Partial charge per atom, in units where the Coulomb constant is 1.
    pub charge: &'a [f64],
}

/// Mutable force accumulators in `SoA` form. Kernels add into these buffers.
#[derive(Debug)]
pub struct ForceBuffers<'a> {
    /// x force components.
    pub fx: &'a mut [f64],
    /// y force components.
    pub fy: &'a mut [f64],
    /// z force components.
    pub fz: &'a mut [f64],
}

impl<'a> ForceBuffers<'a> {
    /// Creates a view. Panics if the slices differ in length.
    pub fn new(fx: &'a mut [f64], fy: &'a mut [f64], fz: &'a mut [f64]) -> Self {
        assert!(fx.len() == fy.len() && fy.len() == fz.len(), "force slices differ in length");
        Self { fx, fy, fz }
    }
}

/// Energies returned by the pairwise kernels, split by term.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Energies {
    /// Lennard-Jones energy.
    pub lj: f64,
    /// Coulomb energy, real-space Ewald term (`qq * erfc(alpha r) / r`).
    pub coulomb: f64,
}

impl Energies {
    /// Sum of all terms.
    pub fn total(&self) -> f64 {
        self.lj + self.coulomb
    }
}

/// Settings shared by the real-space pairwise kernels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonbondedConfig {
    /// Interaction cutoff `rc` (pairs with `r >= rc` contribute nothing). No energy shift.
    pub cutoff: f64,
    /// Ewald splitting parameter `alpha`. Use `0.0` for plain (unscreened) Coulomb.
    pub alpha: f64,
    /// Orthorhombic periodic box edge lengths. `None` disables minimum-image wrapping.
    pub box_len: Option<[f64; 3]>,
}

/// Applies the minimum-image convention to one separation component.
///
/// Valid for orthorhombic boxes; correct for any separation when `cutoff <= box/2`.
#[inline]
pub fn min_image(d: f64, box_edge: f64) -> f64 {
    d - box_edge * (d / box_edge).round()
}
