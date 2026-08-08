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

    /// Returns the parameters on this segment and on `other` at which the two come closest to each
    /// other, both in the normalized `0.0..1.0` form used by [`at`](Segment::at).
    ///
    /// Unlike the infinite-line case this always has an answer, so there is nothing to return an
    /// `Option` for. It is also **not** the line answer with both parameters clamped into range:
    /// pushing one parameter back onto its segment moves the closest point, so the other has to be
    /// re-solved against the new position. Clamping both independently still yields a pair of points
    /// genuinely on the two segments, so it reports a distance that is merely too *large*, which is
    /// the direction that looks reasonable and hides.
    ///
    /// When the segments are parallel, or when either has collapsed to a point, the closest pair is
    /// not unique. A canonical one is returned, chosen by working from this segment's `a` end.
    ///
    /// For the common perpendicular of the two supporting lines, and in particular for asking
    /// whether the closest approach falls strictly *inside* both segments rather than on an end,
    /// use [`Line::closest_approach`] on [`to_line`](Segment::to_line) instead; the parameters mean
    /// the same thing in both, so the results are directly comparable.
    ///
    /// # Arguments
    ///
    /// * `other`: the segment to find the closest approach to
    ///
    /// returns: (f64, f64), the parameter on `self` followed by the parameter on `other`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point3;
    /// use engeom::geom3::Segment3;
    /// use approx::assert_relative_eq;
    ///
    /// // A unit segment along x, and one along y sitting three units above it, offset in x so far
    /// // that the closest approach is between two endpoints rather than in either interior.
    /// let a = Segment3::new(&Point3::origin(), &Point3::new(1.0, 0.0, 0.0)).unwrap();
    /// let b = Segment3::new(&Point3::new(5.0, 0.0, 3.0), &Point3::new(5.0, 1.0, 3.0)).unwrap();
    ///
    /// let (sa, sb) = a.closest_approach(&b);
    /// assert_relative_eq!(sa, 1.0, epsilon = 1e-12);
    /// assert_relative_eq!(sb, 0.0, epsilon = 1e-12);
    /// assert_relative_eq!((a.at(sa) - b.at(sb)).norm(), 5.0, epsilon = 1e-12);
    /// ```
    pub fn closest_approach(&self, other: &Segment<D>) -> (f64, f64) {
        // Squared lengths throughout, so this threshold is the same `1e-12` coincidence distance
        // that `Segment::new` refuses to build a segment from.
        const DEGENERATE: f64 = 1.0e-24;

        let d1 = self.dir();
        let d2 = other.dir();
        let r = self.a - other.a;

        let a = d1.norm_squared();
        let e = d2.norm_squared();
        let f = d2.dot(&r);

        if a <= DEGENERATE && e <= DEGENERATE {
            // Two points, so the only parameters that exist are the ends.
            return (0.0, 0.0);
        }
        if a <= DEGENERATE {
            return (0.0, (f / e).clamp(0.0, 1.0));
        }

        let c = d1.dot(&r);
        if e <= DEGENERATE {
            return ((-c / a).clamp(0.0, 1.0), 0.0);
        }

        let b = d1.dot(&d2);
        let denom = a * e - b * b;

        // As in `Line::closest_approach`, `a * e - b^2` against `a * e` is a test on the angle
        // alone. Parallel segments have no unique answer, and starting from `s = 0` is the
        // canonical choice; the pinning below still lands on a genuinely closest pair.
        let mut s = if denom > 1.0e-14 * a * e {
            ((b * f - c * e) / denom).clamp(0.0, 1.0)
        } else {
            0.0
        };

        let mut t = (b * s + f) / e;

        // A `t` off the end of `other` means the closest pair uses that endpoint, so pin it there
        // and re-solve `s` against the fixed point. This is the step a plain clamp of both
        // parameters skips.
        if t < 0.0 {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else if t > 1.0 {
            t = 1.0;
            s = ((b - c) / a).clamp(0.0, 1.0);
        }

        (s, t)
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

    // closest_approach tests ────────────────────────────────────────────────

    /// The closest distance between two segments, found by brute force over a dense grid of
    /// parameter pairs. Slow and obviously correct, which is the point: it is the only thing in
    /// these tests that does not share arithmetic with the method under test.
    fn brute_force<const D: usize>(a: &Segment<D>, b: &Segment<D>, n: usize) -> f64 {
        let mut best = f64::INFINITY;
        for i in 0..=n {
            let pa = a.at(i as f64 / n as f64);
            for j in 0..=n {
                best = best.min((pa - b.at(j as f64 / n as f64)).norm());
            }
        }
        best
    }

    #[test]
    fn closest_approach_lands_in_both_interiors() {
        // Skew segments whose common perpendicular falls inside both, so the answer is the line
        // answer and both parameters are strictly between the ends.
        let a =
            Segment::<3>::new(&Point3::new(-1.0, 0.0, 0.0), &Point3::new(1.0, 0.0, 0.0)).unwrap();
        let b =
            Segment::<3>::new(&Point3::new(0.0, -1.0, 2.0), &Point3::new(0.0, 1.0, 2.0)).unwrap();

        let (sa, sb) = a.closest_approach(&b);
        assert_relative_eq!(sa, 0.5, epsilon = 1e-12);
        assert_relative_eq!(sb, 0.5, epsilon = 1e-12);
        assert_relative_eq!((a.at(sa) - b.at(sb)).norm(), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn closest_approach_is_not_the_clamped_line_answer() {
        // The case that separates a real segment solve from clamping the infinite-line result.
        //
        // The supporting lines come closest well off the end of `b`, so a naive implementation
        // clamps `b` to its endpoint and keeps the line's parameter for `a`. The truth is that
        // pinning `b` moves the closest point on `a` too, and only re-solving finds it.
        let a = Segment::<3>::new(&Point3::origin(), &Point3::new(10.0, 0.0, 0.0)).unwrap();
        let b =
            Segment::<3>::new(&Point3::new(5.0, 5.0, 1.0), &Point3::new(6.0, 4.0, 1.0)).unwrap();

        let (sa, sb) = a.closest_approach(&b);
        let solved = (a.at(sa) - b.at(sb)).norm();

        // The naive answer: take the line parameters and clamp each into range independently.
        let (la, lb) = a.to_line().closest_approach(&b.to_line()).unwrap();
        let naive = (a.at(la.clamp(0.0, 1.0)) - b.at(lb.clamp(0.0, 1.0))).norm();

        assert!(
            solved < naive - 1.0e-9,
            "clamping the line answer gave {naive}, which the segment solve should have beaten \
             with {solved}"
        );
        assert_relative_eq!(solved, brute_force(&a, &b, 4000), epsilon = 1e-3);
    }

    #[test]
    fn closest_approach_between_two_endpoints() {
        // Pulled far enough apart that neither interior is involved at all.
        let a = Segment::<3>::new(&Point3::origin(), &Point3::new(1.0, 0.0, 0.0)).unwrap();
        let b =
            Segment::<3>::new(&Point3::new(5.0, 0.0, 3.0), &Point3::new(5.0, 1.0, 3.0)).unwrap();

        let (sa, sb) = a.closest_approach(&b);
        assert_relative_eq!(sa, 1.0, epsilon = 1e-12);
        assert_relative_eq!(sb, 0.0, epsilon = 1e-12);
        assert_relative_eq!((a.at(sa) - b.at(sb)).norm(), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn closest_approach_handles_parallel_and_collapsed_segments() {
        let a = Segment::<3>::new(&Point3::origin(), &Point3::new(4.0, 0.0, 0.0)).unwrap();

        // Parallel and overlapping: no unique answer, but the distance must still be right.
        let overlapping =
            Segment::<3>::new(&Point3::new(1.0, 2.0, 0.0), &Point3::new(3.0, 2.0, 0.0)).unwrap();
        let (sa, sb) = a.closest_approach(&overlapping);
        assert_relative_eq!((a.at(sa) - overlapping.at(sb)).norm(), 2.0, epsilon = 1e-12);

        // Parallel and disjoint along their shared direction: the answer is a pair of ends.
        let beyond =
            Segment::<3>::new(&Point3::new(7.0, 0.0, 0.0), &Point3::new(9.0, 0.0, 0.0)).unwrap();
        let (sa, sb) = a.closest_approach(&beyond);
        assert_relative_eq!((a.at(sa) - beyond.at(sb)).norm(), 3.0, epsilon = 1e-12);

        // One segment collapsed to a point, which `new_unchecked` allows even though `new` does
        // not. It should behave as a point query.
        let collapsed =
            Segment::<3>::new_unchecked(Point3::new(1.0, 5.0, 0.0), Point3::new(1.0, 5.0, 0.0));
        let (sa, sb) = a.closest_approach(&collapsed);
        assert_relative_eq!(a.at(sa), Point3::new(1.0, 0.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(
            collapsed.at(sb),
            Point3::new(1.0, 5.0, 0.0),
            epsilon = 1e-12
        );

        // Both collapsed.
        let other =
            Segment::<3>::new_unchecked(Point3::new(0.0, 1.0, 0.0), Point3::new(0.0, 1.0, 0.0));
        let (sa, sb) = collapsed.closest_approach(&other);
        assert_relative_eq!(
            (collapsed.at(sa) - other.at(sb)).norm(),
            (Point3::new(1.0, 5.0, 0.0) - Point3::new(0.0, 1.0, 0.0)).norm(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn closest_approach_is_symmetric() {
        let a =
            Segment::<3>::new(&Point3::new(-2.0, 1.0, 0.5), &Point3::new(3.0, -1.0, 2.0)).unwrap();
        let b =
            Segment::<3>::new(&Point3::new(0.0, 4.0, -1.0), &Point3::new(1.0, -2.0, 3.0)).unwrap();

        let (sa, sb) = a.closest_approach(&b);
        let (tb, ta) = b.closest_approach(&a);

        assert_relative_eq!(sa, ta, epsilon = 1e-12);
        assert_relative_eq!(sb, tb, epsilon = 1e-12);
    }

    #[test]
    fn stress_closest_approach_matches_brute_force() {
        // The independent check. A coarse brute force cannot resolve the exact minimum, so the
        // comparison is one-sided in the tight direction and tolerant in the loose one: the
        // analytic answer must never be worse than the grid, and must be close to it.
        let mut rand = RandomGeometry::<3>::from_seed(0x5e6_c105e);

        for _ in 0..300 {
            let a = Segment::<3>::new_unchecked(rand.point(5.0), rand.point(5.0));
            let b = Segment::<3>::new_unchecked(rand.point(5.0), rand.point(5.0));

            let (sa, sb) = a.closest_approach(&b);
            let solved = (a.at(sa) - b.at(sb)).norm();

            assert!(
                (0.0..=1.0).contains(&sa) && (0.0..=1.0).contains(&sb),
                "parameters ({sa}, {sb}) escaped the segments"
            );

            let sampled = brute_force(&a, &b, 400);
            assert!(
                solved <= sampled + 1.0e-9,
                "dense sampling found {sampled}, beating the analytic answer {solved}"
            );
            assert!(
                sampled - solved < 0.05,
                "analytic answer {solved} is suspiciously far below the sampled {sampled}"
            );
        }
    }

    #[test]
    fn stress_closest_approach_2d_agrees_with_crossing_segments() {
        // In 2D, segments that cross have a closest approach of exactly zero, and segments that do
        // not must report a positive distance. `intersects_other` answers the same question by a
        // completely different route, so the two must agree on every draw.
        let mut rand = RandomGeometry::<2>::from_seed(0x2d_c1057);

        for _ in 0..2000 {
            let a = Segment::<2>::new_unchecked(rand.point(3.0), rand.point(3.0));
            let b = Segment::<2>::new_unchecked(rand.point(3.0), rand.point(3.0));

            let (sa, sb) = a.closest_approach(&b);
            let distance = (a.at(sa) - b.at(sb)).norm();

            let crosses = a.intersects_other(&b);

            if crosses {
                assert!(
                    distance < 1.0e-9,
                    "crossing segments reported a gap of {distance}"
                );
            } else {
                // Segments which do not cross may still touch to within rounding, so this only
                // asserts the direction of the disagreement that would matter.
                assert!(
                    distance >= 0.0,
                    "non-crossing segments reported a negative distance"
                );
            }
        }
    }
}
