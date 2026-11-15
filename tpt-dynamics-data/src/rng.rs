//! Small deterministic pseudo-random number generator (`SplitMix64`).
//!
//! Used by the synthetic data generator and by the fuzz-style tests. It is not
//! cryptographic. The same seed always yields the same sequence on every platform.

/// `SplitMix64` generator with a 64-bit state.
///
/// ```
/// use tpt_dynamics_data::rng::SplitMix64;
/// let mut a = SplitMix64::new(42);
/// let mut b = SplitMix64::new(42);
/// assert_eq!(a.next_u64(), b.next_u64());
/// ```
#[derive(Clone, Debug)]
pub struct SplitMix64 {
    state: u64,
    // Second Box-Muller deviate, cached between calls to `gaussian`.
    spare: Option<f64>,
}

impl SplitMix64 {
    /// Create a generator from `seed`. Every seed, including 0, is valid.
    pub fn new(seed: u64) -> Self {
        Self { state: seed, spare: None }
    }

    /// Next 64 random bits.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform value in `[0, 1)` with 53 bits of precision.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform value in `[lo, hi)`.
    #[inline]
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    /// Uniform integer in `0..n`. Returns 0 when `n == 0`.
    ///
    /// Uses rejection sampling, so the result is unbiased.
    pub fn below(&mut self, n: u64) -> u64 {
        if n == 0 {
            return 0;
        }
        let zone = u64::MAX - (u64::MAX % n);
        loop {
            let v = self.next_u64();
            if v < zone {
                return v % n;
            }
        }
    }

    /// Standard normal deviate (Box-Muller transform).
    pub fn gaussian(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        // Map to (0, 1] so that ln() is finite.
        let u1 = 1.0 - self.next_f64();
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = core::f64::consts::TAU * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }

    /// Exponential deviate with the given slope `tau` (mean `tau`).
    pub fn exponential(&mut self, tau: f64) -> f64 {
        -tau * (1.0 - self.next_f64()).ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_and_seed_sensitive() {
        let a: Vec<u64> = {
            let mut r = SplitMix64::new(7);
            (0..8).map(|_| r.next_u64()).collect()
        };
        let b: Vec<u64> = {
            let mut r = SplitMix64::new(7);
            (0..8).map(|_| r.next_u64()).collect()
        };
        let c = SplitMix64::new(8).next_u64();
        assert_eq!(a, b);
        assert_ne!(a[0], c);
    }

    #[test]
    fn reference_first_output() {
        // First SplitMix64 output for seed 0 (published reference value).
        assert_eq!(SplitMix64::new(0).next_u64(), 0xE220_A839_7B1D_CDAF);
    }

    #[test]
    fn uniform_range_and_mean() {
        let mut r = SplitMix64::new(1);
        let n = 100_000;
        let mut sum = 0.0;
        for _ in 0..n {
            let x = r.next_f64();
            assert!((0.0..1.0).contains(&x));
            sum += x;
        }
        assert!((sum / f64::from(n) - 0.5).abs() < 0.01);
    }

    #[test]
    fn below_is_in_range() {
        let mut r = SplitMix64::new(3);
        assert_eq!(r.below(0), 0);
        for _ in 0..1000 {
            assert!(r.below(7) < 7);
        }
    }

    #[test]
    fn gaussian_moments() {
        let mut r = SplitMix64::new(11);
        let n = 200_000;
        let (mut s, mut s2) = (0.0, 0.0);
        for _ in 0..n {
            let z = r.gaussian();
            s += z;
            s2 += z * z;
        }
        let mean = s / f64::from(n);
        let var = s2 / f64::from(n) - mean * mean;
        assert!(mean.abs() < 0.02);
        assert!((var - 1.0).abs() < 0.02);
    }
}
