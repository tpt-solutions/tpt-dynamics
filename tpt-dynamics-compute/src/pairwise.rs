//! Scalar pair interaction formulas (the reference definitions).
//!
//! Each function takes the squared separation `r2` (or the charge product) and
//! returns `(energy, force_over_r)` where the force vector on atom `i` from
//! pair `(i, j)` is `force_over_r * (r_i - r_j)`.

use libm::{erf, erfc, exp};

/// Inverse square root of pi.
const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;

/// Lennard-Jones 12-6 pair term, `4 eps ((s/r)^12 - (s/r)^6)`.
///
/// Returns `(energy, f_over_r)`. `r2` must be positive.
#[inline]
pub fn lj_pair(r2: f64, sigma: f64, epsilon: f64) -> (f64, f64) {
    let sr2 = sigma * sigma / r2;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    let energy = 4.0 * epsilon * (sr12 - sr6);
    let f_over_r = 24.0 * epsilon * (2.0 * sr12 - sr6) / r2;
    (energy, f_over_r)
}

/// Real-space Ewald Coulomb term, `qq * erfc(alpha r) / r`.
///
/// `qq` is the product of the two charges. With `alpha = 0` this is plain Coulomb.
/// Returns `(energy, f_over_r)`. `r2` must be positive.
#[inline]
pub fn coulomb_ewald_pair(r2: f64, qq: f64, alpha: f64) -> (f64, f64) {
    let r = r2.sqrt();
    let ar = alpha * r;
    let ec = erfc(ar);
    let energy = qq * ec / r;
    let f_over_r = qq * (ec + 2.0 * ar * FRAC_1_SQRT_PI * exp(-ar * ar)) / (r2 * r);
    (energy, f_over_r)
}

/// Excluded-pair correction for PME: `-qq * erf(alpha r) / r`.
///
/// Removes the reciprocal-space contribution of a pair that is excluded from the
/// real-space list (for example a bonded pair). Returns `(energy, f_over_r)`.
#[inline]
pub fn excluded_ewald_pair(r2: f64, qq: f64, alpha: f64) -> (f64, f64) {
    let r = r2.sqrt();
    let ar = alpha * r;
    let energy = -qq * erf(ar) / r;
    // From dE/dr = -qq (c - erf) / r^2 with c = 2 alpha r exp(-a^2 r^2) / sqrt(pi),
    // the radial force is qq (c - erf) / r^2, so F/r = qq (c - erf) / r^3.
    let c = 2.0 * ar * FRAC_1_SQRT_PI * exp(-ar * ar);
    let f_over_r = qq * (c - erf(ar)) / (r2 * r);
    (energy, f_over_r)
}

/// Lorentz-Berthelot combination: arithmetic mean of sigmas, geometric mean of epsilons.
#[inline]
pub fn lorentz_berthelot(sigma_i: f64, eps_i: f64, sigma_j: f64, eps_j: f64) -> (f64, f64) {
    (0.5 * (sigma_i + sigma_j), (eps_i * eps_j).sqrt())
}
