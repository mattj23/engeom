//! Deterministic pseudo-random generation for tests.
//!
//! This is a simple alternative to using `rand` which I think is fine for now. Nothing
//! we're testing depends on statistical quality, we just need reproducible scattered
//! values across the range, so Xorshift64 should be fine.  If it's not, we'll switch
//! to `rand.`

/// A seeded xorshift64 generator.
pub struct Rng(u64);

impl Rng {
    /// Seed the generator. Any seed is accepted; zero is nudged off the fixed point.
    pub fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// A value in `[0, 1)`, using the 53 bits an `f64` can hold exactly.
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A value in `[lo, hi)`.
    pub fn in_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    /// A point whose every coordinate is drawn from `[lo, hi)`.
    pub fn point<const N: usize>(&mut self, lo: f64, hi: f64) -> [f64; N] {
        std::array::from_fn(|_| self.in_range(lo, hi))
    }

    /// `count` points whose every coordinate is drawn from `[lo, hi)`.
    pub fn points<const N: usize>(&mut self, count: usize, lo: f64, hi: f64) -> Vec<[f64; N]> {
        (0..count).map(|_| self.point(lo, hi)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_reproducible_from_a_seed() {
        let a: Vec<u64> = (0..8).map(|_| Rng::new(42).next_u64()).collect();
        assert!(a.windows(2).all(|w| w[0] == w[1]), "same seed must replay");

        let mut x = Rng::new(1);
        let mut y = Rng::new(2);
        assert_ne!(x.next_u64(), y.next_u64());
    }

    #[test]
    fn stays_in_range() {
        let mut rng = Rng::new(7);
        for _ in 0..10_000 {
            let v = rng.next_f64();
            assert!((0.0..1.0).contains(&v));

            let w = rng.in_range(-3.0, 5.0);
            assert!((-3.0..5.0).contains(&w));
        }
    }

    /// A generator stuck at zero would silently make every "random" test a constant test.
    #[test]
    fn does_not_collapse_to_a_fixed_point() {
        let mut rng = Rng::new(0);
        let first = rng.next_u64();
        assert_ne!(first, 0);
        assert_ne!(rng.next_u64(), first);
    }
}
