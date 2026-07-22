//! A dimension-generic helper for generating random geometric entities in tests.
//!
//! [`RandomGeometry<D>`] holds the shared, dimension-independent random-generation machinery (a
//! seeded-or-thread RNG source plus scalar, angle, point, and vector sampling). The 2D and 3D
//! specializations [`RandomGeometry2`] and [`RandomGeometry3`] are type aliases that add the few
//! genuinely dimension-specific constructors (isometries, circles) via additional inherent `impl`
//! blocks.

// This entire module is test-only
#![cfg(test)]

use crate::na::{Point, Quaternion, SVector, Translation3, Unit, UnitQuaternion};
use crate::{Circle2, Iso2, Iso3};
use rand::distr::{Distribution, Uniform};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use rand_distr::Normal;
use std::f64::consts::PI;

enum RngSource {
    Seeded(StdRng),
    Thread,
}

/// A helper for generating random geometric entities in tests. Can be constructed with either a
/// fixed seed (for reproducible tests) or using the default thread RNG.
pub struct RandomGeometry<const D: usize> {
    rng: RngSource,
}

impl<const D: usize> Default for RandomGeometry<D> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const D: usize> RandomGeometry<D> {
    /// Create a `RandomGeometry` that uses the default thread RNG.
    pub fn new() -> Self {
        Self {
            rng: RngSource::Thread,
        }
    }

    /// Create a `RandomGeometry` seeded with a fixed value for reproducible tests.
    pub fn from_seed(seed: u64) -> Self {
        Self {
            rng: RngSource::Seeded(StdRng::seed_from_u64(seed)),
        }
    }

    /// Returns a uniformly-distributed `f64` in `[lo, hi)`.
    pub fn f64(&mut self, lo: f64, hi: f64) -> f64 {
        let u = Uniform::new(lo, hi).unwrap();
        match &mut self.rng {
            RngSource::Seeded(r) => u.sample(r),
            RngSource::Thread => u.sample(&mut rand::rng()),
        }
    }

    /// Returns a uniformly-distributed `f64` in `[-hi, hi)`.
    pub fn f64_sym(&mut self, hi: f64) -> f64 {
        self.f64(-hi, hi)
    }

    /// Returns a uniformly-distributed `f64` in `[0, limit)`.
    pub fn positive(&mut self, limit: f64) -> f64 {
        self.f64(0.0, limit)
    }

    /// Returns a random boolean with equal probability of `true` and `false`.
    pub fn bool(&mut self) -> bool {
        match &mut self.rng {
            RngSource::Seeded(r) => r.random_bool(0.5),
            RngSource::Thread => rand::rng().random_bool(0.5),
        }
    }

    /// Returns a normally-distributed `f64` with the given mean and standard deviation.
    pub fn gaussian_f64(&mut self, mean: f64, std_dev: f64) -> f64 {
        let n = Normal::new(mean, std_dev).unwrap();
        match &mut self.rng {
            RngSource::Seeded(r) => n.sample(r),
            RngSource::Thread => n.sample(&mut rand::rng()),
        }
    }

    /// Returns a uniformly-distributed angle in `[-π, π)`.
    pub fn angle_sym_pi(&mut self) -> f64 {
        self.f64(-PI, PI)
    }

    /// Returns a uniformly-distributed angle in `[-2π, 2π)`.
    pub fn angle_sym_2pi(&mut self) -> f64 {
        self.f64(-2.0 * PI, 2.0 * PI)
    }

    /// Returns a uniformly-distributed angle in `[0, 2π)`.
    pub fn angle_pos_2pi(&mut self) -> f64 {
        self.f64(0.0, 2.0 * PI)
    }

    /// Returns a random `Point<f64, D>` with each component in `[-limit, limit)`.
    pub fn point(&mut self, limit: f64) -> Point<f64, D> {
        let mut p = Point::origin();
        for i in 0..D {
            p[i] = self.f64(-limit, limit);
        }
        p
    }

    /// Returns a random `SVector<f64, N>` with each component in `[-limit, limit)`. The vector
    /// dimension `N` is independent of the geometry dimension `D`, which is convenient for sampling
    /// parameter vectors (e.g. a 6-DOF transform storage) from a `RandomGeometry<3>`.
    pub fn vector<const N: usize>(&mut self, limit: f64) -> SVector<f64, N> {
        let mut v = SVector::zeros();
        for i in 0..N {
            v[i] = self.f64(-limit, limit);
        }
        v
    }

    /// Returns an `SVector<f64, N>` whose components are each drawn independently from a zero-mean
    /// normal distribution with the given standard deviation. Useful for adding isotropic Gaussian
    /// noise to points.
    pub fn gaussian_vector<const N: usize>(&mut self, std_dev: f64) -> SVector<f64, N> {
        let mut v = SVector::zeros();
        for i in 0..N {
            v[i] = self.gaussian_f64(0.0, std_dev);
        }
        v
    }

    /// Returns a random unit vector whose direction is uniformly distributed over the unit sphere
    /// in `D` dimensions. This works by drawing each component from a standard normal distribution
    /// and normalizing; a degenerate (near-zero) draw is resampled.
    pub fn unit_vec(&mut self) -> Unit<SVector<f64, D>> {
        loop {
            let v: SVector<f64, D> = self.gaussian_vector(1.0);
            if v.norm() > 1e-12 {
                return Unit::new_normalize(v);
            }
        }
    }
}

/// The 2D specialization of [`RandomGeometry`]. The dimension-independent sampling lives on
/// `RandomGeometry<D>`; the 2D-specific constructors are added below.
pub type RandomGeometry2 = RandomGeometry<2>;

/// The 3D specialization of [`RandomGeometry`]. The dimension-independent sampling lives on
/// `RandomGeometry<D>`; the 3D-specific constructors are added below.
pub type RandomGeometry3 = RandomGeometry<3>;

impl RandomGeometry2 {
    /// Returns a random `Iso2` with translation components in `[-t, t]` and arbitrary rotation.
    pub fn iso2(&mut self, t: f64) -> Iso2 {
        Iso2::new(self.vector(t), self.angle_sym_pi())
    }

    /// Returns a random `Circle2` with its center in `[-center_limit, center_limit]` on each axis
    /// and a radius in `[r_min, r_max)`.
    pub fn circle2(&mut self, center_limit: f64, r_min: f64, r_max: f64) -> Circle2 {
        let c = self.point(center_limit);
        let r = self.f64(r_min, r_max);
        Circle2::new(c.x, c.y, r)
    }
}

impl RandomGeometry3 {
    /// Returns a random `Iso3` with translation components in `[-t, t]` and a rotation drawn
    /// uniformly from `SO(3)` (Haar measure).
    ///
    /// The rotation is sampled by normalizing a 4D standard-normal vector into a unit quaternion:
    /// an isotropic Gaussian projects to a uniform point on the unit 3-sphere, whose double cover
    /// of `SO(3)` yields a uniformly-distributed rotation. This is unbiased, unlike independently
    /// sampled Euler angles, which concentrate near the poles.
    pub fn iso3(&mut self, t: f64) -> Iso3 {
        let translation = Translation3::from(self.vector::<3>(t));
        let g = self.gaussian_vector::<4>(1.0);
        let rotation = UnitQuaternion::new_normalize(Quaternion::new(g[0], g[1], g[2], g[3]));
        Iso3::from_parts(translation, rotation)
    }
}
