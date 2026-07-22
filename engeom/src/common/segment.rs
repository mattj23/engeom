use crate::Result;
use crate::common::Line;
use crate::common::PCoords;
use crate::common::consensus::Magsac;
use crate::common::points::dist;
use crate::na::{AbstractRotation, Isometry, Point, SVector};
use serde::{Deserialize, Serialize};
use std::ops;

/// A line segment in D-dimensional space, defined by two endpoints.
///
/// `Segment<D>` is the base for `Segment2` and `Segment3`, which are two of `engeom`'s geometric
/// primitives.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Segment<const D: usize> {
    pub a: Point<f64, D>,
    pub b: Point<f64, D>,
}

impl<const D: usize> Segment<D> {
    pub fn new_unchecked(a: Point<f64, D>, b: Point<f64, D>) -> Self {
        Self { a, b }
    }

    /// Create a new segment from two points, returning an error if the points are coincident
    /// (within a tolerance of `1e-12`).
    pub fn new(a: &impl PCoords<D>, b: &impl PCoords<D>) -> Result<Self> {
        if dist(a, b) < 1e-12 {
            Err("The two points are too close to each other".into())
        } else {
            Ok(Self::new_unchecked(
                Point::from(a.coords()),
                Point::from(b.coords()),
            ))
        }
    }

    /// Fit a segment to a set of points by ordinary least squares. An infinite line is fit to the
    /// points via [`Line::from_fit`] (a weighted SVD fit), and the segment's endpoints are then set
    /// to the extreme projections of the points onto that line, so the segment spans exactly the
    /// range covered by the input.
    ///
    /// This is not robust to gross outliers, which corrupt both the fitted line and the extent; for
    /// that, use [`Segment::from_consensus`].
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of at least two distinct coordinates to fit the segment to
    /// * `weights`: if `Some`, this must be a slice of floating points the same length as `points`,
    ///   with the weight value to multiply each point residual by. Weights bias the fitted line
    ///   only; the endpoints are still the extreme projections of every point.
    ///
    /// returns: Result<Segment<{ D }>, Box<dyn Error, Global>>
    pub fn from_fit(points: &[impl PCoords<D>], weights: Option<&[f64]>) -> Result<Self> {
        let line = Line::from_fit(points, weights)?;
        let (mut t_min, mut t_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for p in points {
            let t = line.scalar_project(p);
            t_min = t_min.min(t);
            t_max = t_max.max(t);
        }
        if !t_min.is_finite() || !t_max.is_finite() {
            return Err("Failed to determine segment endpoints from the fitted line".into());
        }
        Self::new(&line.at(t_min), &line.at(t_max))
    }

    /// Fit a segment to a set of points using MAGSAC++ robust consensus estimation.
    ///
    /// A robust infinite line is estimated with the same MAGSAC++ consensus fit as
    /// [`Line::from_consensus`], rejecting gross outliers. The segment's endpoints are then set to
    /// the extreme projections of the *inlier* points onto that line, so outliers influence neither
    /// the line nor the segment's extent.
    ///
    /// # Arguments
    ///
    /// * `points`: the points to fit the segment to
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `options`: an optional [`Magsac`] configuration to override the iteration count, refinement
    ///   steps, confidence, or RNG seed. Its `sigma_max` field is overridden by the `sigma_max`
    ///   argument.
    ///
    /// returns: Result<Segment<{ D }>, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[Point<f64, D>],
        sigma_max: f64,
        options: Option<Magsac>,
    ) -> Result<Self> {
        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit = magsac.fit::<D, Line<D>>(points)?;
        let line = fit.model;

        let (mut t_min, mut t_max) = (f64::INFINITY, f64::NEG_INFINITY);
        for &i in &fit.inliers {
            let t = line.scalar_project(&points[i]);
            t_min = t_min.min(t);
            t_max = t_max.max(t);
        }
        if !t_min.is_finite() || !t_max.is_finite() {
            return Err("Consensus fit produced no inliers to bound the segment".into());
        }
        Self::new(&line.at(t_min), &line.at(t_max))
    }

    pub fn dir(&self) -> SVector<f64, D> {
        self.b - self.a
    }

    pub fn at(&self, t: f64) -> Point<f64, D> {
        self.a + t * self.dir()
    }

    pub fn length(&self) -> f64 {
        self.dir().norm()
    }

    /// Returns a new segment with the endpoints reversed.
    pub fn reversed(&self) -> Self {
        Self::new_unchecked(self.b, self.a)
    }

    /// Returns the infinite line passing through this segment's endpoints, in the direction from
    /// `a` to `b`.
    pub fn to_line(&self) -> Line<D> {
        Line::from_points(&self.a, &self.b)
    }

    /// Calculate the scalar projection of a set of coordinates onto the line segment, in which
    /// 0.0 represents a point at the segment's starting point `a` and 1.0 represents a point at
    /// the segment's end point `b`.  The result can be any finite value, including negative ones
    /// or ones greater than zero.
    ///
    /// # Arguments
    ///
    /// * `other`: an entity with coordinates to project onto the segment
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn scalar_projection(&self, other: &impl PCoords<D>) -> f64 {
        let dir = self.dir();
        let test = other.coords() - self.a.coords();
        dir.dot(&test) / dir.norm_squared()
    }

    pub fn closest_point(&self, other: &impl PCoords<D>) -> Point<f64, D> {
        let t = self.scalar_projection(other).clamp(0.0, 1.0);
        self.at(t)
    }

    /// Returns a new segment with both endpoints transformed by the given isometry.
    pub fn transformed_by<R: AbstractRotation<f64, D>>(&self, iso: &Isometry<f64, R, D>) -> Self {
        Self {
            a: iso * self.a,
            b: iso * self.b,
        }
    }
}

impl<const D: usize, R: AbstractRotation<f64, D>> ops::Mul<Segment<D>> for Isometry<f64, R, D> {
    type Output = Segment<D>;

    fn mul(self, rhs: Segment<D>) -> Self::Output {
        rhs.transformed_by(&self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::consensus::Magsac;
    use crate::common::random_geometry::{RandomGeometry, RandomGeometry2};
    use crate::{Point2, Point3};
    use approx::assert_relative_eq;

    /// Assert that `seg` has the same pair of endpoints as `(a, b)`, allowing either ordering (the
    /// fitted direction may run either way along the line).
    fn assert_endpoints<const D: usize>(
        seg: &Segment<D>,
        a: &Point<f64, D>,
        b: &Point<f64, D>,
        eps: f64,
    ) {
        let forward = (seg.a - a).norm() < eps && (seg.b - b).norm() < eps;
        let reversed = (seg.a - b).norm() < eps && (seg.b - a).norm() < eps;
        assert!(
            forward || reversed,
            "segment endpoints {:?}..{:?} do not match {:?}..{:?}",
            seg.a,
            seg.b,
            a,
            b
        );
    }

    // from_fit tests ────────────────────────────────────────────────────────

    #[test]
    fn from_fit_recovers_endpoints_of_clean_span() {
        let a = Point3::new(1.0, 2.0, 3.0);
        let b = Point3::new(4.0, 6.0, 3.0);
        let seg = Segment::<3>::new(&a, &b).unwrap();

        // Evenly spaced points covering the whole span, including both endpoints.
        let points: Vec<Point3> = (0..=10).map(|i| seg.at(i as f64 / 10.0)).collect();

        let fit = Segment::<3>::from_fit(&points, None).unwrap();
        assert_endpoints(&fit, &a, &b, 1e-9);
    }

    #[test]
    fn from_fit_with_noise_is_close() {
        let mut rand = RandomGeometry::<3>::new();
        let a = Point3::new(-2.0, 0.0, 1.0);
        let b = Point3::new(5.0, 3.0, -1.0);
        let seg = Segment::<3>::new(&a, &b).unwrap();

        let points: Vec<Point3> = (0..=40)
            .map(|i| seg.at(i as f64 / 40.0) + rand.gaussian_vector(0.01))
            .collect();

        let fit = Segment::<3>::from_fit(&points, None).unwrap();
        // Endpoints won't be exact under noise, but should be within a few multiples of sigma.
        assert_endpoints(&fit, &a, &b, 0.1);
    }

    #[test]
    fn from_fit_uniform_weights_match_unweighted() {
        let a = Point2::new(0.0, 0.0);
        let b = Point2::new(6.0, 2.0);
        let seg = Segment::<2>::new(&a, &b).unwrap();
        let points: Vec<Point2> = (0..=12).map(|i| seg.at(i as f64 / 12.0)).collect();

        let unweighted = Segment::<2>::from_fit(&points, None).unwrap();
        let weights = vec![1.0; points.len()];
        let weighted = Segment::<2>::from_fit(&points, Some(&weights)).unwrap();
        assert_relative_eq!(weighted.a, unweighted.a, epsilon = 1e-10);
        assert_relative_eq!(weighted.b, unweighted.b, epsilon = 1e-10);
    }

    #[test]
    fn from_fit_single_point_is_error() {
        let points = vec![Point3::new(1.0, 2.0, 3.0)];
        assert!(Segment::<3>::from_fit(&points, None).is_err());
    }

    // from_consensus tests ────────────────────────────────────────────────────

    #[test]
    fn from_consensus_endpoints_use_only_inliers() {
        let mut rand = RandomGeometry2::new();
        let a = Point2::new(0.0, 0.5);
        let b = Point2::new(8.0, 0.5);
        let seg = Segment::<2>::new(&a, &b).unwrap();

        // Inliers spanning the segment with small perpendicular noise.
        let mut points: Vec<Point2> = (0..=80)
            .map(|i| seg.at(i as f64 / 80.0) + rand.gaussian_vector(0.01))
            .collect();

        // A dense cluster of gross outliers well off the line and far beyond endpoint `b`, which
        // would badly stretch the segment if they were allowed to influence the extent.
        let outlier_center = Point2::new(20.0, 6.0);
        for _ in 0..40 {
            points.push(outlier_center + rand.gaussian_vector(0.5));
        }

        let magsac = Magsac {
            sigma_max: 0.03,
            max_iterations: Some(400),
            refinement_steps: 4,
            confidence: 0.99,
            seed: Some(42),
        };
        let fit = Segment::<2>::from_consensus(&points, 0.03, Some(magsac)).unwrap();

        // The endpoints should track the inlier span, not be dragged toward the outlier cluster.
        assert_endpoints(&fit, &a, &b, 0.1);
    }

    #[test]
    fn from_consensus_recovers_clean_segment() {
        let a = Point3::new(1.0, -2.0, 3.0);
        let b = Point3::new(4.0, 2.0, -1.0);
        let seg = Segment::<3>::new(&a, &b).unwrap();
        let points: Vec<Point3> = (0..=60).map(|i| seg.at(i as f64 / 60.0)).collect();

        let fit = Segment::<3>::from_consensus(&points, 0.01, None).unwrap();
        assert_endpoints(&fit, &a, &b, 1e-3);
    }
}
