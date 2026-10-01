//! Small, deterministic, dependency-free pseudo-random number generation.
//!
//! Used by the Langevin thermostat (`tpt-chem-md`) and the Gillespie
//! stochastic simulation algorithm (`tpt-chem-kinetics`). Keeping the RNG
//! in-house means simulations are bit-for-bit reproducible across platforms
//! with a given seed, and the dependency graph stays license-clean.
//!
//! The generator is xoshiro256\*\* seeded through SplitMix64
//! (Steele, Lea & Flood 2021). It is excellent for simulation purposes but
//! **not cryptographically secure**.

use crate::num;

/// xoshiro256\*\* generator with SplitMix64 seeding.
#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// Create a generator from a 64-bit seed (any seed is valid).
    pub fn new(seed: u64) -> Self {
        let mut sm = SplitMix64 { state: seed };
        Rng {
            s: [sm.next(), sm.next(), sm.next(), sm.next()],
        }
    }

    /// Next raw 64-bit value (xoshiro256\*\* scrambler).
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;

        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);

        result
    }

    /// Uniform sample in `[0, 1)`.
    pub fn uniform(&mut self) -> f64 {
        // 53 significant bits — highest precision an f64 fraction can hold.
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Uniform sample in `[lo, hi)`. Panics if `lo >= hi`.
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        assert!(lo < hi, "range: lo must be < hi");
        lo + (hi - lo) * self.uniform()
    }

    /// Standard normal sample via Box–Muller (one normal per call; the
    /// polar mate is discarded for simplicity and reproducibility).
    pub fn normal(&mut self, mean: f64, stddev: f64) -> f64 {
        let u1 = self.uniform().max(f64::MIN_POSITIVE);
        let u2 = self.uniform();
        let r = num::sqrt(-2.0 * num::ln(u1));
        let theta = core::f64::consts::TAU * u2;
        mean + stddev * r * num::cos(theta)
    }

    /// Uniform integer in `[0, n)`. Panics if `n == 0`.
    pub fn usize_below(&mut self, n: usize) -> usize {
        assert!(n > 0, "usize_below: n must be nonzero");
        // Lemire's nearly-divisionless method.
        let n64 = n as u64;
        let mut x = self.next_u64();
        let mut m = u128::from(x) * u128::from(n64);
        let mut l = m as u64;
        if l < n64 {
            let threshold = n64.wrapping_neg() % n64;
            while l < threshold {
                x = self.next_u64();
                m = u128::from(x) * u128::from(n64);
                l = m as u64;
            }
        }
        (m >> 64) as usize
    }
}

/// SplitMix64 — used only to expand a single seed into the xoshiro state.
struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_with_seed() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn uniform_in_range_and_reasonably_spread() {
        let mut rng = Rng::new(7);
        let mut sum = 0.0;
        for _ in 0..10_000 {
            let u = rng.uniform();
            assert!((0.0..1.0).contains(&u));
            sum += u;
        }
        let mean = sum / 10_000.0;
        assert!((mean - 0.5).abs() < 0.02, "mean was {mean}");
    }

    #[test]
    fn normal_has_unit_mean_stddev() {
        let mut rng = Rng::new(123);
        let (mut s1, mut s2) = (0.0, 0.0);
        let n = 50_000;
        for _ in 0..n {
            let x = rng.normal(0.0, 1.0);
            s1 += x;
            s2 += x * x;
        }
        let mean = s1 / f64::from(n);
        let var = s2 / f64::from(n) - mean * mean;
        assert!(mean.abs() < 0.02, "mean {mean}");
        assert!((var - 1.0).abs() < 0.05, "var {var}");
    }

    #[test]
    fn usize_below_covers_all_buckets() {
        let mut rng = Rng::new(99);
        let mut seen = [false; 5];
        for _ in 0..1000 {
            let i = rng.usize_below(5);
            seen[i] = true;
        }
        assert!(seen.iter().all(|&s| s));
    }
}
