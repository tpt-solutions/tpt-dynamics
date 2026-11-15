//! Particle-mesh Ewald (PME) reciprocal-space solver for orthorhombic boxes.
//!
//! Implements the smooth-PME scheme of Essmann et al. (1995): charges are spread
//! onto a grid with cardinal B-splines, the grid is transformed with the pure-Rust
//! FFT in [`crate::fft`], the reciprocal energy is evaluated with the Ewald kernel
//! and the B-spline moduli, and forces are obtained by differentiating the spline
//! weights. Spread, gather and FFT lines run in parallel with rayon.
//!
//! Conventions: energies are in units where the Coulomb constant is 1, so a
//! unit-charge pair at distance `r` has energy `1/r`. The short-range part (real
//! space, `erfc`) and the excluded-pair correction live in [`crate::scalar`] and
//! [`crate::simd_kernel`]; this module provides the reciprocal sum, the self-energy
//! term and the parameter selection.
//!
//! Total Ewald energy for a periodic system is
//! `E = E_real + E_reciprocal + E_self (+ E_excluded correction)`.

use rayon::prelude::*;

use crate::fft::{Complex, Fft3d};
use crate::types::{ForceBuffers, Positions};

/// Largest supported B-spline order.
pub const MAX_ORDER: usize = 12;

/// Target maximum grid spacing used by [`PmeParams::from_accuracy`] (in the same
/// length units as the box, e.g. nm).
pub const DEFAULT_MAX_SPACING: f64 = 0.12;

/// Solves `erfc(alpha * rc) = tol` for `alpha` by bisection.
///
/// This is the standard heuristic: the real-space term at the cutoff is then
/// `tol` times the bare Coulomb value.
pub fn alpha_for_tolerance(rc: f64, tol: f64) -> f64 {
    assert!(rc > 0.0 && tol > 0.0 && tol < 1.0, "need rc > 0 and 0 < tol < 1");
    let (mut lo, mut hi) = (0.0_f64, 20.0 / rc);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        // erfc is decreasing, so too-small alpha gives a value above tol.
        if libm::erfc(mid * rc) > tol {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// PME parameters: Ewald splitting, grid dimensions and spline order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PmeParams {
    /// Ewald splitting parameter `alpha`.
    pub alpha: f64,
    /// FFT grid `[nx, ny, nz]`. Each entry is a power of two.
    pub grid: [usize; 3],
    /// Cardinal B-spline order (number of grid points per atom per axis), 2..=12.
    pub order: usize,
}

impl PmeParams {
    /// Chooses `alpha` from the accuracy target `erfc(alpha rc) = tol` and grid
    /// dimensions so that the spacing is at most [`DEFAULT_MAX_SPACING`], rounded up
    /// to a power of two (and at least `order`).
    pub fn from_accuracy(box_len: [f64; 3], rc: f64, tol: f64, order: usize) -> Self {
        let alpha = alpha_for_tolerance(rc, tol);
        let grid = box_len.map(|l| {
            let n = (l / DEFAULT_MAX_SPACING).ceil() as usize;
            n.max(order).next_power_of_two()
        });
        Self { alpha, grid, order }
    }
}

/// Cardinal B-spline weights and derivatives for one atom along one axis.
#[derive(Clone, Copy, Debug)]
struct Spline {
    /// Wrapped grid indices, `idx[t]` pairs with `w[t]`.
    idx: [usize; MAX_ORDER],
    /// Spline weights `M_p(w + t)`.
    w: [f64; MAX_ORDER],
    /// Derivatives `M_p'(w + t)`.
    dw: [f64; MAX_ORDER],
}

/// Computes the spline of order `p` for scaled coordinate `u` on a grid of `k` points.
fn spline(u: f64, p: usize, k: usize) -> Spline {
    let fl = u.floor();
    let w = u - fl;
    // Order 2: M_2(w) = w, M_2(w + 1) = 1 - w.
    let mut th = [0.0; MAX_ORDER];
    th[0] = w;
    th[1] = 1.0 - w;
    let mut prev = [0.0; MAX_ORDER];
    // Build up to order p with M_n(x) = x/(n-1) M_{n-1}(x) + (n-x)/(n-1) M_{n-1}(x-1).
    for n in 3..=p {
        prev = th;
        let nf = n as f64;
        for j in (0..n).rev() {
            let jf = j as f64;
            let a = (w + jf) * th[j];
            let b = if j >= 1 { (nf - w - jf) * th[j - 1] } else { 0.0 };
            th[j] = (a + b) / (nf - 1.0);
        }
    }
    let mut dw = [0.0; MAX_ORDER];
    if p == 2 {
        dw[0] = 1.0;
        dw[1] = -1.0;
    } else {
        // M_p'(x) = M_{p-1}(x) - M_{p-1}(x - 1); `prev` holds order p-1.
        for t in 0..p {
            dw[t] = prev[t] - if t >= 1 { prev[t - 1] } else { 0.0 };
        }
    }
    let mut idx = [0usize; MAX_ORDER];
    for t in 0..p {
        idx[t] = (fl as i64 - t as i64).rem_euclid(k as i64) as usize;
    }
    Spline { idx, w: th, dw }
}

/// Precomputed PME state for one box and parameter set.
#[derive(Clone, Debug)]
pub struct Pme {
    box_len: [f64; 3],
    params: PmeParams,
    fft: Fft3d,
    /// Combined Ewald kernel and B-spline moduli `C(m) B(m)`, flat, zero at `m = 0`.
    influence: Vec<f64>,
}

impl Pme {
    /// Builds the solver. Panics if the grid is not a power of two or the order is
    /// outside `2..=MAX_ORDER`.
    pub fn new(box_len: [f64; 3], params: PmeParams) -> Self {
        assert!((2..=MAX_ORDER).contains(&params.order), "PME order must be in 2..={MAX_ORDER}");
        let [kx, ky, kz] = params.grid;
        let volume = box_len[0] * box_len[1] * box_len[2];
        let fft = Fft3d::new(params.grid);

        // |b(k)|^2 per axis, with b(k) = e^{2 pi i (p-1) k / K} / sum_t M_p(t+1) e^{2 pi i k t / K}.
        let bmod: Vec<Vec<f64>> =
            (0..3).map(|axis| bspline_modulus(params.order, params.grid[axis])).collect();

        let alpha = params.alpha;
        let pi = std::f64::consts::PI;
        let mut influence = vec![0.0; kx * ky * kz];
        for ix in 0..kx {
            for iy in 0..ky {
                for iz in 0..kz {
                    let flat = (ix * ky + iy) * kz + iz;
                    if ix == 0 && iy == 0 && iz == 0 {
                        continue;
                    }
                    let idx = [ix, iy, iz];
                    let mut m2 = 0.0;
                    for axis in 0..3 {
                        let k = params.grid[axis];
                        let s = if idx[axis] > k / 2 {
                            idx[axis] as f64 - k as f64
                        } else {
                            idx[axis] as f64
                        };
                        let m = s / box_len[axis];
                        m2 += m * m;
                    }
                    let b = 1.0 / (bmod[0][ix] * bmod[1][iy] * bmod[2][iz]);
                    let c = (-pi * pi * m2 / (alpha * alpha)).exp() / (pi * volume * m2);
                    influence[flat] = c * b;
                }
            }
        }
        Self { box_len, params, fft, influence }
    }

    /// Parameters this solver was built with.
    pub fn params(&self) -> PmeParams {
        self.params
    }

    /// Self-energy correction `-alpha/sqrt(pi) * sum_i q_i^2`.
    pub fn self_energy(&self, charges: &[f64]) -> f64 {
        let q2: f64 = charges.iter().map(|q| q * q).sum();
        -self.params.alpha / std::f64::consts::PI.sqrt() * q2
    }

    /// Reciprocal-space energy, with forces added to `forces`.
    pub fn reciprocal(
        &self,
        pos: Positions<'_>,
        charges: &[f64],
        forces: &mut ForceBuffers<'_>,
    ) -> f64 {
        assert_eq!(pos.len(), charges.len(), "charges length does not match positions");
        let mut grid = self.spread(pos, charges);
        self.fft.transform(&mut grid, false);

        // E = 1/2 sum_m G(m) |F(Q)(m)|^2
        let energy = 0.5
            * grid
                .par_iter()
                .zip(self.influence.par_iter())
                .map(|(g, &gi)| gi * g.norm_sqr())
                .sum::<f64>();

        grid.par_iter_mut().zip(self.influence.par_iter()).for_each(|(g, &gi)| *g = g.scale(gi));
        self.fft.transform(&mut grid, true);

        let f = self.gather(pos, charges, &grid);
        for (a, fa) in f.iter().enumerate() {
            forces.fx[a] += fa[0];
            forces.fy[a] += fa[1];
            forces.fz[a] += fa[2];
        }
        energy
    }

    /// Scaled spline coordinates of atom `a` along each axis.
    fn splines(&self, pos: &Positions<'_>, a: usize) -> [Spline; 3] {
        let p = self.params.order;
        let k = self.params.grid;
        let l = self.box_len;
        let coords = [pos.x[a], pos.y[a], pos.z[a]];
        std::array::from_fn(|axis| {
            let u = coords[axis].rem_euclid(l[axis]) / l[axis] * k[axis] as f64;
            spline(u, p, k[axis])
        })
    }

    /// Spreads charges onto the grid (parallel over atoms with per-thread grids).
    fn spread(&self, pos: Positions<'_>, charges: &[f64]) -> Vec<Complex> {
        let [kx, ky, kz] = self.params.grid;
        let p = self.params.order;
        let len = kx * ky * kz;
        let n = pos.len();
        // One private grid per chunk of atoms (about one per thread), summed at the end.
        // A grid per rayon split would allocate far more memory than the work needs.
        let chunks = rayon::current_num_threads().max(1);
        let per_chunk = n.div_ceil(chunks).max(1);
        (0..n.div_ceil(per_chunk))
            .into_par_iter()
            .map(|c| {
                let mut grid = vec![Complex::default(); len];
                for a in c * per_chunk..((c + 1) * per_chunk).min(n) {
                    let s = self.splines(&pos, a);
                    let q = charges[a];
                    for tx in 0..p {
                        for ty in 0..p {
                            let wxy = q * s[0].w[tx] * s[1].w[ty];
                            for tz in 0..p {
                                let flat = (s[0].idx[tx] * ky + s[1].idx[ty]) * kz + s[2].idx[tz];
                                grid[flat].re += wxy * s[2].w[tz];
                            }
                        }
                    }
                }
                grid
            })
            .reduce_with(|mut acc, part| {
                for (x, y) in acc.iter_mut().zip(part) {
                    *x += y;
                }
                acc
            })
            .unwrap_or_else(|| vec![Complex::default(); len])
    }

    /// Interpolates the potential grid back to atoms: `F_i = -q_i (K/L) grad_u phi`.
    fn gather(&self, pos: Positions<'_>, charges: &[f64], phi: &[Complex]) -> Vec<[f64; 3]> {
        let [_, ky, kz] = self.params.grid;
        let p = self.params.order;
        let scale: [f64; 3] =
            std::array::from_fn(|axis| self.params.grid[axis] as f64 / self.box_len[axis]);
        (0..pos.len())
            .into_par_iter()
            .map(|a| {
                let s = self.splines(&pos, a);
                let mut g = [0.0; 3];
                for tx in 0..p {
                    for ty in 0..p {
                        for tz in 0..p {
                            let v = phi[(s[0].idx[tx] * ky + s[1].idx[ty]) * kz + s[2].idx[tz]].re;
                            let (wx, wy, wz) = (s[0].w[tx], s[1].w[ty], s[2].w[tz]);
                            let (dx, dy, dz) = (s[0].dw[tx], s[1].dw[ty], s[2].dw[tz]);
                            g[0] += v * dx * wy * wz;
                            g[1] += v * wx * dy * wz;
                            g[2] += v * wx * wy * dz;
                        }
                    }
                }
                let qa = -charges[a];
                [qa * scale[0] * g[0], qa * scale[1] * g[1], qa * scale[2] * g[2]]
            })
            .collect()
    }
}

/// `|b(k)|^2` along one axis of a `k_len` grid for B-spline order `p`.
fn bspline_modulus(p: usize, k_len: usize) -> Vec<f64> {
    let theta = spline(0.0, p, k_len);
    // M_p(t + 1) for t = 0..=p-2 are the spline values at integers 1..p-1.
    let m: Vec<f64> = (0..p - 1).map(|t| theta.w[t + 1]).collect();
    (0..k_len)
        .map(|k| {
            let (mut re, mut im) = (0.0, 0.0);
            for (t, &mt) in m.iter().enumerate() {
                let ang = 2.0 * std::f64::consts::PI * (k * t) as f64 / k_len as f64;
                re += mt * libm::cos(ang);
                im += mt * libm::sin(ang);
            }
            let denom = re * re + im * im;
            if denom > 1e-300 {
                denom
            } else {
                f64::INFINITY
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spline_is_partition_of_unity_and_derivative_is_consistent() {
        for p in 2..=MAX_ORDER {
            for &u in &[0.01, 0.3, 0.77, 5.5] {
                let s = spline(u, p, 16);
                let sum: f64 = s.w[..p].iter().sum();
                assert!((sum - 1.0).abs() < 1e-12, "order {p} u {u}: sum {sum}");
                let dsum: f64 = s.dw[..p].iter().sum();
                assert!(dsum.abs() < 1e-12, "derivatives must sum to zero");
                // Derivative vs central difference of the weight at the same grid index.
                let h = 1e-6;
                let sp = spline(u + h, p, 16);
                let sm = spline(u - h, p, 16);
                for t in 0..p {
                    let num = (sp.w[t] - sm.w[t]) / (2.0 * h);
                    assert!((num - s.dw[t]).abs() < 1e-6, "order {p} t {t}: {num} vs {}", s.dw[t]);
                }
            }
        }
    }

    #[test]
    fn spline_indices_wrap_periodically() {
        let s = spline(0.5, 4, 8);
        // Grid points floor(u)-t: 0, -1, -2, -3 wrapped to 0, 7, 6, 5.
        assert_eq!(&s.idx[..4], &[0, 7, 6, 5]);
    }
}
