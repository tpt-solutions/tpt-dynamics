//! Small deterministic pseudo-random number generator (xoshiro256** seeded by
//! `SplitMix64`). Written in-crate to avoid a dependency; the stream is stable
//! across platforms, so seeded runs are reproducible.

/// Deterministic random number generator.
#[derive(Clone, Debug)]
pub struct Rng {
    state: [u64; 4],
}

impl Rng {
    /// Creates a generator from a 64-bit seed.
    pub fn seed_from_u64(seed: u64) -> Self {
        let mut sm = seed;
        let mut next = || {
            sm = sm.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = sm;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        let state = [next(), next(), next(), next()];
        Self { state }
    }

    /// The full internal state, for checkpointing.
    pub fn state(&self) -> [u64; 4] {
        self.state
    }

    /// Restores a generator from a state saved with [`Rng::state`]. An all-zero
    /// state (which xoshiro cannot leave) is replaced by a fixed non-zero state.
    pub fn from_state(state: [u64; 4]) -> Self {
        if state == [0; 4] {
            Self::seed_from_u64(0)
        } else {
            Self { state }
        }
    }

    /// Returns the next raw 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        let s = &mut self.state;
        let result = s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    /// Returns a uniform sample in `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Returns a standard normal sample (Box-Muller transform).
    pub fn gaussian(&mut self) -> f64 {
        loop {
            let u1 = self.uniform();
            if u1 > 0.0 {
                let u2 = self.uniform();
                return (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Rng;

    #[test]
    fn same_seed_gives_same_stream() {
        let mut a = Rng::seed_from_u64(7);
        let mut b = Rng::seed_from_u64(7);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn uniform_in_unit_interval() {
        let mut r = Rng::seed_from_u64(1);
        for _ in 0..10_000 {
            let u = r.uniform();
            assert!((0.0..1.0).contains(&u));
        }
    }

    #[test]
    fn gaussian_moments() {
        let mut r = Rng::seed_from_u64(99);
        let n: usize = 200_000;
        let (mut s1, mut s2) = (0.0, 0.0);
        for _ in 0..n {
            let g = r.gaussian();
            s1 += g;
            s2 += g * g;
        }
        let mean = s1 / n as f64;
        let var = s2 / n as f64 - mean * mean;
        assert!(mean.abs() < 0.01, "mean {mean}");
        assert!((var - 1.0).abs() < 0.02, "variance {var}");
    }
}
