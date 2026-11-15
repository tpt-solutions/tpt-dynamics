//! Pure-Rust radix-2 complex FFT, 1D and 3D.
//!
//! Only power-of-two lengths are supported. PME grid sizes are chosen as powers of
//! two (see [`crate::pme::PmeParams::from_accuracy`]), which keeps this
//! implementation small and dependency-free. Transforms are unnormalised in both
//! directions: `inverse(forward(x)) == n * x`.
//!
//! The 3D transform runs 1D transforms along the contiguous last axis in parallel
//! (rayon), then rotates the axes so that each of the three dimensions takes a turn
//! at being last.

use rayon::prelude::*;

use std::ops::{Add, AddAssign, Mul, Sub};

/// Complex number with `f64` parts.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex {
    /// Creates `re + i im`.
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// Complex conjugate.
    #[must_use]
    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    /// Squared modulus `re^2 + im^2`.
    pub fn norm_sqr(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// Scales both parts by a real factor.
    #[must_use]
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.re * s, self.im * s)
    }
}

impl Add for Complex {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, o: Self) {
        self.re += o.re;
        self.im += o.im;
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
}

impl Mul for Complex {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
}

/// Precomputed plan for a radix-2 FFT of one length.
#[derive(Clone, Debug)]
pub struct Fft1d {
    n: usize,
    /// Bit-reversal permutation.
    rev: Vec<usize>,
    /// `exp(-2 pi i k / n)` for `k < n/2`.
    twiddle: Vec<Complex>,
}

impl Fft1d {
    /// Plans a transform of length `n`. Panics unless `n` is a power of two.
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two(), "FFT length {n} is not a power of two");
        let bits = n.trailing_zeros();
        let rev = (0..n)
            .map(|i| if bits == 0 { 0 } else { i.reverse_bits() >> (usize::BITS - bits) })
            .collect();
        let twiddle = (0..n / 2)
            .map(|k| {
                let ang = -2.0 * std::f64::consts::PI * k as f64 / n as f64;
                Complex::new(libm::cos(ang), libm::sin(ang))
            })
            .collect();
        Self { n, rev, twiddle }
    }

    /// Transform length.
    pub fn len(&self) -> usize {
        self.n
    }

    /// True for the degenerate length-zero plan (never constructed by [`Fft1d::new`]).
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// In-place transform of `buf` (length `n`). `inverse` selects `exp(+i ...)`.
    /// Neither direction is normalised.
    pub fn transform(&self, buf: &mut [Complex], inverse: bool) {
        let n = self.n;
        assert_eq!(buf.len(), n, "buffer length does not match plan");
        for i in 0..n {
            let j = self.rev[i];
            if j > i {
                buf.swap(i, j);
            }
        }
        let mut len = 2;
        while len <= n {
            let half = len / 2;
            let step = n / len;
            for start in (0..n).step_by(len) {
                for k in 0..half {
                    let mut w = self.twiddle[k * step];
                    if inverse {
                        w = w.conj();
                    }
                    let t = w * buf[start + k + half];
                    let u = buf[start + k];
                    buf[start + k] = u + t;
                    buf[start + k + half] = u - t;
                }
            }
            len <<= 1;
        }
    }
}

/// Plan for a 3D transform of a row-major `[nx, ny, nz]` grid (z fastest).
#[derive(Clone, Debug)]
pub struct Fft3d {
    dims: [usize; 3],
    /// Plans in the order they are needed: last axis of the original, then the
    /// rotated layouts (dims[2], dims[0], dims[1]).
    plans: [Fft1d; 3],
}

impl Fft3d {
    /// Plans a 3D transform. Every dimension must be a power of two.
    pub fn new(dims: [usize; 3]) -> Self {
        let plans = [Fft1d::new(dims[2]), Fft1d::new(dims[0]), Fft1d::new(dims[1])];
        Self { dims, plans }
    }

    /// Grid dimensions `[nx, ny, nz]`.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// In-place 3D transform. `inverse` selects the `exp(+i ...)` kernel.
    /// Unnormalised: a forward plus an inverse pass multiplies the data by `nx*ny*nz`.
    pub fn transform(&self, data: &mut [Complex], inverse: bool) {
        let [a, b, c] = self.dims;
        assert_eq!(data.len(), a * b * c, "grid size does not match plan");
        let mut dims = [a, b, c];
        let mut tmp = vec![Complex::default(); data.len()];
        for plan in &self.plans {
            let n = dims[2];
            data.par_chunks_mut(n).for_each(|line| plan.transform(line, inverse));
            rotate(data, &mut tmp, dims);
            data.copy_from_slice(&tmp);
            dims = [dims[1], dims[2], dims[0]];
        }
    }
}

/// Rotates axes `[a, b, c] -> [b, c, a]`: `out[j, k, i] = in[i, j, k]`.
fn rotate(src: &[Complex], dst: &mut [Complex], dims: [usize; 3]) {
    let [a, b, c] = dims;
    dst.par_chunks_mut(c * a).enumerate().for_each(|(j, out)| {
        for k in 0..c {
            for i in 0..a {
                out[k * a + i] = src[(i * b + j) * c + k];
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_naive_dft_and_round_trips() {
        let n = 16;
        let x: Vec<Complex> = (0..n)
            .map(|i| Complex::new((i as f64 * 0.37).sin(), (i as f64 * 0.11).cos()))
            .collect();
        let mut y = x.clone();
        Fft1d::new(n).transform(&mut y, false);
        for k in 0..n {
            let mut acc = Complex::default();
            for (j, xj) in x.iter().enumerate() {
                let ang = -2.0 * std::f64::consts::PI * (j * k) as f64 / n as f64;
                acc += *xj * Complex::new(libm::cos(ang), libm::sin(ang));
            }
            assert!((acc.re - y[k].re).abs() < 1e-12 && (acc.im - y[k].im).abs() < 1e-12);
        }
        Fft1d::new(n).transform(&mut y, true);
        for (a, b) in y.iter().zip(&x) {
            assert!((a.re / n as f64 - b.re).abs() < 1e-12);
            assert!((a.im / n as f64 - b.im).abs() < 1e-12);
        }
    }

    #[test]
    fn three_d_round_trip_and_separable_check() {
        let dims = [4, 8, 2];
        let plan = Fft3d::new(dims);
        let total: usize = 4 * 8 * 2;
        let x: Vec<Complex> = (0..total).map(|i| Complex::new(i as f64 % 7.0, 0.0)).collect();
        let mut y = x.clone();
        plan.transform(&mut y, false);
        // DC term equals the sum of all entries.
        let sum: f64 = x.iter().map(|c| c.re).sum();
        assert!((y[0].re - sum).abs() < 1e-9 && y[0].im.abs() < 1e-9);
        plan.transform(&mut y, true);
        for (a, b) in y.iter().zip(&x) {
            assert!((a.re / total as f64 - b.re).abs() < 1e-10);
        }
    }
}
