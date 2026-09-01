//! This module contains a dimension-generic parameterized line, `Line<D>`, which serves as the
//! shared foundation for the 2D and 3D line types. It is structured similarly to
//! [`SurfacePoint`](crate::common::SurfacePoint): the generic features live here, while the
//! dimension-specific behavior (normals and intersections in 2D, plane/sphere/circle
//! intersections in 3D, axis constructors, and isometry multiplication operators) lives in the
//! `geom2` and `geom3` modules.

use crate::Result;
use crate::common::PCoords;
use crate::common::consensus::{ConsensusModel, Magsac};
use crate::common::svd_basis::SvdBasis;
use crate::na::{AbstractRotation, Isometry, Point, SVector};
use serde::{Deserialize, Serialize};

/// A parameterized line in D-dimensional space: `P(t) = origin + t * direction`.
///
/// The direction is not required to be normalized; use [`new_normalize`](Line::new_normalize)
/// for a unit-speed parameterization where `t` equals arc length from the origin.
///
/// `Line<D>` is the base for `Line2` and `Line3`, which are two of `engeom`'s geometric primitives.
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct Line<const D: usize> {
    pub origin: Point<f64, D>,
    pub direction: SVector<f64, D>,
}

impl<const D: usize> Line<D> {
    /// Create a line from an origin point and a direction vector (stored as-is, not normalized).
    pub fn new(origin: Point<f64, D>, direction: SVector<f64, D>) -> Self {
        Self { origin, direction }
    }

    /// Create a line from an origin point and a direction vector, normalizing the direction so
    /// that the parameter `t` equals arc length from the origin.
    pub fn new_normalize(origin: Point<f64, D>, direction: SVector<f64, D>) -> Self {
        Self::new(origin, direction.normalize())
    }

    /// Create a line through two points. The direction is `p2 - p1` (not normalized).
    pub fn from_points(p1: &impl PCoords<D>, p2: &impl PCoords<D>) -> Self {
        let origin = Point::from(p1.coords());
        Self::new(origin, p2.coords() - p1.coords())
    }

    /// Fit a line to a set of points using singular value decomposition, resulting in a
    /// least-squares fitting. The resulting parameterized line will have its t=0 sitting at the
    /// center of the SVD result.
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of coordinates to fit the line to
    /// * `weights`: if `Some`, this must be a slice of floating points the same length as `points`,
    ///   with the weight value to multiply each point residual by.
    ///
    /// returns: Result<Line<{ D }>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_fit(points: &[impl PCoords<D>], weights: Option<&[f64]>) -> Result<Self> {
        let basis = SvdBasis::from_points(points, weights)
            .ok_or("Failed to fit line with singular value decomposition")?;
        Ok(Line::new(basis.center, basis.largest().into_inner()))
    }

    /// Fit a line to a set of points using MAGSAC++ robust consensus estimation.
    ///
    /// Unlike an ordinary least-squares fit ([`Line::from_fit`]), this rejects gross outliers by
    /// taking an upper bound on the inlier noise (`sigma_max`) rather than a hard inlier/outlier
    /// threshold, and refines each candidate with noise-marginalized iteratively reweighted least
    /// squares. It is substantially less sensitive to `sigma_max` than RANSAC is to its threshold,
    /// as long as `sigma_max` is not chosen smaller than the actual noise.
    ///
    /// The resulting line has its origin at the centroid of the inlier set and a unit-length
    /// direction vector.
    ///
    /// # Arguments
    ///
    /// * `points`: the points to fit the line to
    /// * `sigma_max`: the upper bound on the expected inlier noise, in the same units as the points
    /// * `options`: an optional [`Magsac`] configuration to override the iteration count, refinement
    ///   steps, confidence, or RNG seed. Its `sigma_max` field is overridden by the `sigma_max`
    ///   argument.
    ///
    /// returns: Result<Line<{ D }>, Box<dyn Error, Global>>
    pub fn from_consensus(
        points: &[Point<f64, D>],
        sigma_max: f64,
        options: Option<Magsac>,
    ) -> Result<Self> {
        let mut magsac = options.unwrap_or_else(|| Magsac::new(sigma_max));
        magsac.sigma_max = sigma_max;

        let fit = magsac.fit::<D, Line<D>>(points)?;
        Ok(fit.model)
    }

    /// Returns a new line with the same origin, but with the direction inverted.
    pub fn reversed(&self) -> Self {
        Self::new(self.origin, -self.direction)
    }

    /// Returns a new line with the same origin but a normalized direction, so that `t` equals arc
    /// length from the origin.
    pub fn normalized(&self) -> Self {
        Self::new(self.origin, self.direction.normalize())
    }

    /// Returns the point at parameter `t`: `P(t) = origin + t * direction`.
    pub fn at(&self, t: f64) -> Point<f64, D> {
        self.origin + self.direction * t
    }

    /// Returns a new line with the origin shifted by a given amount along the direction of the
    /// line. The direction of the new line is the same as the original line. The original is left
    /// unchanged.
    ///
    /// If the direction is not of unit length, keep in mind this shift will be proportional to
    /// the length of the direction vector.
    pub fn shifted_origin(&self, delta_t: f64) -> Self {
        Self::new(self.origin + self.direction * delta_t, self.direction)
    }

    /// Returns the parameter `t` such that `P(t)` is the closest point on the line to `point`.
    /// Equivalent to the scalar projection of `(point - origin)` onto `direction`, divided by
    /// `|direction|²`.
    pub fn scalar_project(&self, point: &impl PCoords<D>) -> f64 {
        (point.coords() - self.origin.coords).dot(&self.direction) / self.direction.norm_squared()
    }

    /// Returns the closest point on the line to `point`.
    pub fn closest_point(&self, point: &impl PCoords<D>) -> Point<f64, D> {
        self.at(self.scalar_project(point))
    }

    /// Returns the perpendicular distance from `point` to the line.
    pub fn distance_to(&self, point: &impl PCoords<D>) -> f64 {
        let pt = Point::from(point.coords());
        (pt - self.closest_point(&pt)).norm()
    }

    /// Returns the parameters on this line and on `other` at which the two come closest to each
    /// other, or `None` if they are parallel.
    ///
    /// The points [`at`](Line::at) these two parameters are the feet of the common perpendicular:
    /// the segment joining them is orthogonal to both directions, and its length is the distance
    /// between the lines. Where the lines actually meet, which is everywhere they are not parallel
    /// in 2D, both parameters land on the single intersection point and this is equivalent to
    /// [`intersection_param`](crate::geom2::intersection_param).
    ///
    /// Each parameter is in its own line's parameterization, so it is an arc length only if that
    /// line's direction is unit length. A line built by [`to_line`](crate::common::Segment::to_line)
    /// has direction `b - a`, so its parameter is the segment's own normalized parameter; that is
    /// what makes "where, if anywhere, do these two segments cross" this call plus a range check on
    /// both results.
    ///
    /// Parallelism is judged relative to the lengths of the two direction vectors rather than
    /// against a fixed distance, so the answer does not change when the same geometry is expressed
    /// at a different scale. A line with a zero-length direction is treated as parallel to
    /// everything, since it has no direction to be perpendicular to.
    ///
    /// This is the closest approach of two *infinite* lines. For the closest approach of the two
    /// finite spans, see [`Segment::closest_approach`](crate::common::Segment::closest_approach),
    /// which is not simply this result clamped.
    ///
    /// # Arguments
    ///
    /// * `other`: the line to find the closest approach to
    ///
    /// returns: Option<(f64, f64)>, the parameter on `self` followed by the parameter on `other`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Line3, Point3, Vector3};
    /// use approx::assert_relative_eq;
    ///
    /// // The x axis, and a line running parallel to the y axis which passes over it at x = 2,
    /// // five units up in z.
    /// let a = Line3::new(Point3::origin(), Vector3::x());
    /// let b = Line3::new(Point3::new(2.0, 0.0, 5.0), Vector3::y());
    ///
    /// let (ta, tb) = a.closest_approach(&b).unwrap();
    /// assert_relative_eq!(ta, 2.0, epsilon = 1e-12);
    /// assert_relative_eq!(tb, 0.0, epsilon = 1e-12);
    /// assert_relative_eq!((a.at(ta) - b.at(tb)).norm(), 5.0, epsilon = 1e-12);
    ///
    /// // Parallel lines have no unique closest approach.
    /// let c = Line3::new(Point3::new(0.0, 1.0, 0.0), Vector3::x());
    /// assert!(a.closest_approach(&c).is_none());
    /// ```
    pub fn closest_approach(&self, other: &Line<D>) -> Option<(f64, f64)> {
        let w = self.origin - other.origin;
        let a = self.direction.norm_squared();
        let b = self.direction.dot(&other.direction);
        let c = other.direction.norm_squared();
        let d = self.direction.dot(&w);
        let e = other.direction.dot(&w);

        // `a * c - b^2` is `|d1|^2 |d2|^2 sin^2(theta)`, so comparing it against `a * c` tests the
        // angle between the lines and nothing else. A zero-length direction makes `a * c` zero and
        // fails the test, which is the intended answer.
        let denom = a * c - b * b;
        if denom <= 1.0e-14 * a * c {
            return None;
        }

        Some(((b * e - c * d) / denom, (a * e - b * d) / denom))
    }

    /// Return a new line whose direction is spherically interpolated between this instance and
    /// `other`, and whose origin is _linearly interpolated_ between this instance and `other`.
    ///
    /// A value of `t=0` will return this line, and a value of `t=1` will return `other`.
    ///
    /// # Arguments
    ///
    /// * `other`: the other line to interpolate to
    /// * `t`: the interpolation parameter, does _not_ need to be bounded between [0, 1]
    ///
    /// returns: Line<{ D }>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn slerp(&self, other: &Line<D>, t: f64) -> Self {
        let new_direction = self.direction.slerp(&other.direction, t);
        let shift = other.origin - self.origin;
        Self::new(self.origin + shift * t, new_direction)
    }

    /// Returns a new line with both origin and direction transformed by the given isometry.
    pub fn transformed_by<R>(&self, iso: &Isometry<f64, R, D>) -> Self
    where
        R: AbstractRotation<f64, D>,
    {
        Self::new(
            iso * self.origin,
            iso.rotation.transform_vector(&self.direction),
        )
    }
}

impl<const D: usize> ConsensusModel<D> for Line<D> {
    type Point = Point<f64, D>;
    const SAMPLE_SIZE: usize = 2;

    fn from_sample(sample: &[Point<f64, D>]) -> Option<Self> {
        let line = Line::from_points(&sample[0], &sample[1]);
        // Reject degenerate (coincident) samples, whose direction vector is effectively zero.
        (line.direction.norm() > 1e-12).then_some(line)
    }

    fn residual(&self, point: &Point<f64, D>) -> f64 {
        self.distance_to(point)
    }

    fn refine_weighted(points: &[Point<f64, D>], weights: &[f64], _initial: &Self) -> Option<Self> {
        // The MAGSAC++ refinement step is a single weighted least-squares fit, which for a line is
        // exactly the weighted SVD fit provided by `from_fit`.
        Line::from_fit(points, Some(weights)).ok()
    }
}

impl<const D: usize> PCoords<D> for Line<D> {
    fn coords(&self) -> SVector<f64, D> {
        self.origin.coords
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::{RandomGeometry, RandomGeometry2};
    use crate::{Point2, Point3, Vector2, Vector3};
    use approx::assert_relative_eq;
    use rand::RngExt;

    #[test]
    fn at_endpoints() {
        let line = Line::<2>::new(Point2::new(1.0, 2.0), Vector2::new(0.0, 1.0));
        assert_relative_eq!(line.at(0.0), Point2::new(1.0, 2.0), epsilon = 1e-12);
        assert_relative_eq!(line.at(1.0), Point2::new(1.0, 3.0), epsilon = 1e-12);
    }

    #[test]
    fn new_normalize_gives_unit_direction() {
        let line = Line::<3>::new_normalize(Point3::origin(), Vector3::new(3.0, 0.0, 0.0));
        assert_relative_eq!(line.direction.norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn from_points_direction_is_difference() {
        let line = Line::<3>::from_points(&Point3::new(1.0, 0.0, 0.0), &Point3::new(4.0, 0.0, 0.0));
        assert_relative_eq!(line.direction, Vector3::new(3.0, 0.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn closest_point_perpendicular_drop() {
        let line = Line::<3>::new(Point3::origin(), Vector3::x());
        assert_relative_eq!(
            line.closest_point(&Point3::new(4.0, 3.0, 0.0)),
            Point3::new(4.0, 0.0, 0.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn distance_to_known_value() {
        let line = Line::<3>::new(Point3::origin(), Vector3::x());
        assert_relative_eq!(
            line.distance_to(&Point3::new(0.0, 3.0, 4.0)),
            5.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn shifted_origin_moves_along_direction() {
        let line = Line::<2>::new(Point2::new(0.0, 0.0), Vector2::new(2.0, 0.0));
        let shifted = line.shifted_origin(1.5);
        assert_relative_eq!(shifted.origin, Point2::new(3.0, 0.0), epsilon = 1e-12);
        let shifted_again = shifted.shifted_origin(-1.0);
        assert_relative_eq!(shifted_again.origin, Point2::new(1.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn from_fit_recovers_axis_aligned_line() {
        let points = vec![
            Point3::new(0.0, 1.0, 2.0),
            Point3::new(1.0, 1.0, 2.0),
            Point3::new(2.0, 1.0, 2.0),
            Point3::new(3.0, 1.0, 2.0),
        ];
        let line = Line::<3>::from_fit(&points, None).unwrap();
        assert_relative_eq!(line.origin, Point3::new(1.5, 1.0, 2.0), epsilon = 1e-10);
        let dir = line.direction.normalize();
        assert_relative_eq!(dir.x.abs(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(dir.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(dir.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn from_fit_origin_is_at_centroid() {
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(6.0, 0.0),
        ];
        let line = Line::<2>::from_fit(&points, None).unwrap();
        assert_relative_eq!(line.origin, Point2::new(3.0, 0.0), epsilon = 1e-10);
    }

    #[test]
    fn from_fit_minimizes_perpendicular_residuals() {
        // A noisy but clearly linear point set: every point should end up close to the fitted
        // line, much closer than to a line built from just two of the points.
        let points = vec![
            Point2::new(0.0, 0.05),
            Point2::new(1.0, 0.95),
            Point2::new(2.0, 2.05),
            Point2::new(3.0, 2.95),
            Point2::new(4.0, 4.05),
        ];
        let line = Line::<2>::from_fit(&points, None).unwrap();
        for p in &points {
            assert!(
                line.distance_to(p) < 0.1,
                "point {p:?} too far from fitted line"
            );
        }
    }

    #[test]
    fn from_fit_uniform_weights_match_unweighted() {
        let points = vec![
            Point3::new(0.0, 0.0, 5.0),
            Point3::new(1.0, 0.3, 5.2),
            Point3::new(0.2, 1.0, 4.8),
            Point3::new(1.0, 1.0, 5.1),
        ];
        let unweighted = Line::<3>::from_fit(&points, None).unwrap();
        let weights = vec![1.0; points.len()];
        let weighted = Line::<3>::from_fit(&points, Some(&weights)).unwrap();
        assert_relative_eq!(weighted.origin, unweighted.origin, epsilon = 1e-10);
        assert_relative_eq!(weighted.direction, unweighted.direction, epsilon = 1e-10);
    }

    #[test]
    fn from_fit_heavily_weighted_points_pull_fit_toward_them() {
        // Two clusters of points on either side of the origin; heavily weighting one cluster
        // should pull the fitted line's origin towards it.
        let points = vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 10.0),
            Point2::new(1.0, 10.0),
        ];
        let weights = vec![1.0, 1.0, 100.0, 100.0];
        let line = Line::<2>::from_fit(&points, Some(&weights)).unwrap();
        assert!(line.origin.y > 5.0);
    }

    #[test]
    fn from_fit_empty_points_is_error() {
        let points: Vec<Point3> = vec![];
        assert!(Line::<3>::from_fit(&points, None).is_err());
    }

    #[test]
    fn from_fit_single_point_is_error() {
        let points = vec![Point3::new(1.0, 2.0, 3.0)];
        assert!(Line::<3>::from_fit(&points, None).is_err());
    }

    #[test]
    fn from_fit_coincident_points_is_error() {
        let points = vec![Point3::new(1.0, 1.0, 1.0), Point3::new(1.0, 1.0, 1.0)];
        assert!(Line::<3>::from_fit(&points, None).is_err());
    }

    fn line_noise<const D: usize>(line: &Line<D>, n: usize, sigma: f64) -> Vec<Point<f64, D>> {
        let mut rand = RandomGeometry::<D>::new();
        (0..n)
            .map(|_| line.at(rand.f64_sym(10.0)) + rand.gaussian_vector(sigma))
            .collect()
    }

    #[test]
    fn stress_from_fit_recovers_known_line() {
        let mut rng = rand::rng();
        for _ in 0..200 {
            let origin = Point3::new(
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
                rng.random_range(-10.0..10.0),
            );
            let direction = Vector3::new(
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
                rng.random_range(-1.0..1.0),
            )
            .normalize();
            let true_line = Line::<3>::new(origin, direction);

            let points: Vec<Point3> = (0..20).map(|i| true_line.at(i as f64 - 10.0)).collect();
            let fit = Line::<3>::from_fit(&points, None).unwrap();

            // The fitted line should pass through (or very near) every input point, since they
            // are exactly collinear.
            for p in &points {
                assert_relative_eq!(fit.distance_to(p), 0.0, epsilon = 1e-8);
            }

            // The fitted direction should be parallel (or anti-parallel) to the true direction.
            let fit_dir = fit.direction.normalize();
            assert_relative_eq!(fit_dir.dot(&direction).abs(), 1.0, epsilon = 1e-8);
        }
    }

    // from_consensus tests ────────────────────────────────────────────────────────

    #[test]
    fn from_consensus_rejects_outliers() {
        let mut rand = RandomGeometry2::new();

        // Inliers lie along the line y = 0.5 with a small deterministic perturbation.
        let true_line = Line::<2>::new(Point2::new(0.0, 0.5), Vector2::new(1.0, 0.0));
        let decoy_line = Line::<2>::new(Point2::new(0.0, 0.0), Vector2::new(1.0, 1.0));
        let inliers = line_noise(&true_line, 150, 0.01);
        let mut points = [inliers.clone(), line_noise(&decoy_line, 50, 0.05)].concat();

        // A dense cluster of gross outliers well off the line.
        let center = Point2::new(-2.0, 3.0);
        for _i in 0..50 {
            points.push(center + rand.gaussian_vector(1.0));
        }

        let magsac = Magsac {
            sigma_max: 0.02,
            max_iterations: Some(400),
            refinement_steps: 4,
            confidence: 0.99,
            seed: Some(42),
        };
        let fit = magsac
            .fit_filtered::<2, Line<2>, _>(&points, |_| true)
            .unwrap();

        // Every inlier should lie very close to the recovered line
        for i in inliers.iter() {
            assert!(fit.model.distance_to(i) < 0.01 * 6.0);
        }

        // The direction should be parallel (or anti-parallel) to the true horizontal line.
        let dir = fit.model.direction.normalize();
        assert_relative_eq!(
            dir.dot(&true_line.direction.normalize()).abs(),
            1.0,
            epsilon = 1e-2
        );
    }

    #[test]
    fn from_consensus_convenience_recovers_line_3d() {
        // A clean set of collinear 3D points; the convenience method should recover the line.
        let true_line =
            Line::<3>::new_normalize(Point3::new(1.0, -2.0, 3.0), Vector3::new(0.5, 1.0, -0.25));
        let points: Vec<Point3> = (0..60)
            .map(|i| true_line.at(i as f64 * 0.1 - 3.0))
            .collect();

        let line = Line::<3>::from_consensus(&points, 0.01, None).unwrap();

        for p in &points {
            assert!(line.distance_to(p) < 1e-3);
        }
        let dir = line.direction.normalize();
        assert_relative_eq!(
            dir.dot(&true_line.direction.normalize()).abs(),
            1.0,
            epsilon = 1e-6
        );
    }

    // closest_approach tests ────────────────────────────────────────────────

    #[test]
    fn closest_approach_of_skew_axes() {
        // The x axis and a line along y passing over it at x = 2, three units up in z. The common
        // perpendicular is the z segment between them, so both feet are at the crossing in plan.
        let a = Line::<3>::new(Point3::origin(), Vector3::x());
        let b = Line::<3>::new(Point3::new(2.0, 0.0, 3.0), Vector3::y());

        let (ta, tb) = a.closest_approach(&b).unwrap();
        assert_relative_eq!(ta, 2.0, epsilon = 1e-12);
        assert_relative_eq!(tb, 0.0, epsilon = 1e-12);
        assert_relative_eq!((a.at(ta) - b.at(tb)).norm(), 3.0, epsilon = 1e-12);
    }

    #[test]
    fn closest_approach_of_intersecting_lines_is_the_intersection() {
        // Two lines which genuinely meet, at (3, 4, 0). Both feet must land on that single point.
        let a = Line::<3>::new(Point3::origin(), Vector3::new(3.0, 4.0, 0.0));
        let b = Line::<3>::new(Point3::new(3.0, 0.0, 0.0), Vector3::new(0.0, 4.0, 0.0));

        let (ta, tb) = a.closest_approach(&b).unwrap();
        assert_relative_eq!(a.at(ta), Point3::new(3.0, 4.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(b.at(tb), Point3::new(3.0, 4.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn closest_approach_in_2d_matches_the_intersection_parameters() {
        // In 2D non-parallel lines always meet, so the closest approach and the intersection are
        // the same question asked twice. This pins the two answers together.
        let a = Line::<2>::new(Point2::new(1.0, 1.0), Vector2::new(2.0, 1.0));
        let b = Line::<2>::new(Point2::new(4.0, 0.0), Vector2::new(-1.0, 3.0));

        let (ta, tb) = a.closest_approach(&b).unwrap();
        let (ia, ib) =
            crate::geom2::intersection_param(&a.origin, &a.direction, &b.origin, &b.direction)
                .unwrap();

        assert_relative_eq!(ta, ia, epsilon = 1e-12);
        assert_relative_eq!(tb, ib, epsilon = 1e-12);
        assert_relative_eq!(a.at(ta), b.at(tb), epsilon = 1e-12);
    }

    #[test]
    fn closest_approach_refuses_parallel_and_degenerate_lines() {
        let x = Line::<3>::new(Point3::origin(), Vector3::x());

        // Parallel, same direction.
        let parallel = Line::<3>::new(Point3::new(0.0, 1.0, 0.0), Vector3::x());
        assert!(x.closest_approach(&parallel).is_none());

        // Antiparallel, and scaled, which must not change the verdict.
        let anti = Line::<3>::new(Point3::new(0.0, 1.0, 0.0), Vector3::new(-7.0, 0.0, 0.0));
        assert!(x.closest_approach(&anti).is_none());

        // Collinear with itself.
        assert!(x.closest_approach(&x).is_none());

        // A zero-length direction has no direction to be perpendicular to.
        let point = Line::<3>::new(Point3::new(0.0, 1.0, 0.0), Vector3::zeros());
        assert!(x.closest_approach(&point).is_none());
        assert!(point.closest_approach(&x).is_none());
    }

    #[test]
    fn closest_approach_is_scale_invariant() {
        // The same two lines, one of them with a direction vector a million times longer. The feet
        // must not move, so the parameter must shrink by exactly that factor.
        let a = Line::<3>::new(Point3::origin(), Vector3::x());
        let b = Line::<3>::new(Point3::new(2.0, 0.0, 3.0), Vector3::y());
        let b_long = Line::<3>::new(b.origin, b.direction * 1.0e6);

        let (ta, tb) = a.closest_approach(&b).unwrap();
        let (ta_long, tb_long) = a.closest_approach(&b_long).unwrap();

        assert_relative_eq!(ta, ta_long, epsilon = 1e-12);
        assert_relative_eq!(tb, tb_long * 1.0e6, epsilon = 1e-12);
    }

    #[test]
    fn stress_closest_approach_is_perpendicular_to_both() {
        // The defining property: the segment joining the two feet is orthogonal to both directions.
        // Nothing here is compared against another implementation, only against the definition.
        let mut rand = RandomGeometry::<3>::from_seed(0x5eed_1abe);

        for _ in 0..2000 {
            let a = Line::<3>::new(rand.point(10.0), rand.unit_vec().into_inner());
            let b = Line::<3>::new(rand.point(10.0), rand.unit_vec().into_inner());

            // Skip the near-parallel draws, where the feet are ill-conditioned by nature rather
            // than by any fault of the arithmetic.
            let sin_sq = 1.0 - a.direction.dot(&b.direction).powi(2);
            if sin_sq < 1.0e-6 {
                continue;
            }

            let (ta, tb) = a.closest_approach(&b).unwrap();
            let joining = b.at(tb) - a.at(ta);

            assert_relative_eq!(joining.dot(&a.direction), 0.0, epsilon = 1e-9);
            assert_relative_eq!(joining.dot(&b.direction), 0.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn stress_closest_approach_beats_dense_sampling() {
        // The independent check: no pair of points sampled along the two lines may be closer than
        // the analytic answer claims. Measured rather than asserted, in the manner of the rest of
        // the crate.
        let mut rand = RandomGeometry::<3>::from_seed(0xc10_5e57);

        for _ in 0..200 {
            let a = Line::<3>::new(rand.point(5.0), rand.unit_vec().into_inner());
            let b = Line::<3>::new(rand.point(5.0), rand.unit_vec().into_inner());

            let sin_sq = 1.0 - a.direction.dot(&b.direction).powi(2);
            if sin_sq < 1.0e-4 {
                continue;
            }

            let (ta, tb) = a.closest_approach(&b).unwrap();
            let best = (a.at(ta) - b.at(tb)).norm();

            // Sample a window around the analytic feet; nothing in it may be closer.
            let n = 60;
            for i in 0..=n {
                let sa = ta + 4.0 * (i as f64 / n as f64 - 0.5);
                for j in 0..=n {
                    let sb = tb + 4.0 * (j as f64 / n as f64 - 0.5);
                    assert!(
                        (a.at(sa) - b.at(sb)).norm() >= best - 1.0e-9,
                        "sampled pair at ({sa}, {sb}) beat the analytic closest approach {best}"
                    );
                }
            }
        }
    }
}
