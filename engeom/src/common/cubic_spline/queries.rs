//! This module is to group spatial query related actions on the cubic Bezier spline struct.
//!
//! The main objective of this module is to get an accurate, reliable, and reasonably fast
//! implementation of the closest point query algorithm.  The common approaches of lookup table
//! pre-search can produce incorrect results when a query is done near spots where the curve
//! loops around or crosses itself, and you end up pruning the interval with the true minimum value
//! because an endpoint happened to be closer.
//!
//! Besides the incorrect pruning problem, the other danger is not handling searches correctly in
//! areas of the curve with extreme changes that aren't represented well by the discrete points at
//! interval boundaries.
//!
//! Initial thoughts:
//! - I want to avoid having to supply a tolerance for an acceleration structure
//! - Generalizing over two and three dimensions may be making this more complicated
//! - If there's some way I could break the spline up into simpler segments with some known maximum
//!   complexity, that would probably let me accomplish what I want
//! - It would be nice to have a single acceleration structure that could also be used as a
//!   foundation for intersection queries
//! - Bounding boxes feel heavy, but may be necessary
//!
//! Some stuff I think I know about cubic Béziers:
//! - Derivative roots for the individual components can be found via quadratic formula
//! - In two or three dimensions, it's only possible to have one cusp, and that cusp can be found
//!   (I think?) where the derivatives of all dimensions go to zero at the same place.
//! - In 2D, splitting at a cusp leaves you with two extremely simple arc-like curves
//! - What if we split so that we get segments that are convex and non-intersecting?
//!     - If there's a cusp, split there
//!     - There should be at most two inflection points on both 2D or 3D(?) splines
//!     - Curves that are loops can have no inflection points, but they aren't convex
//! - The goal should be to split into segments where extreme points are easily found in arbitrary
//!   directions?
//!
//! What types of queries and so what types of extremes?
//! - Closest point query: minimal distance to coordinates
//! - Closest point to halfspace: minimum distance to halfspace
//! - Halfspace intersection: minimum and maximum distance to halfspace
//! - Closest point to circle/sphere: minimum and maximum distance to coordinates
//! - Circle/sphere intersection: minimum and maximum distance to coordinates
//! - Intersections with other spline?
//! - Closest distance to other spline?
//!
//! What condition would create more than one local minima? How would you find both?
//! - Under any(?) circumstances, a local minima _or_ maxima will be at one of these:
//!     - A point where the derivative is orthogonal to the query direction
//!     - One of the interval endpoints
//!
//! Ok, how about this. We create bounding capsules, recursively breaking up the spline until we
//! reach a fixed number of divisions. We know that a cubic spline can only change directions
//! twice, so the priority for splitting should be (I think):
//! - Split on a cusp if one exists
//! - Find any inflection points (there will be 0, 1, or 2) and split between them
//! - Split at the maximum distance from the base leg
//!
//! After the first split, there should be no cusps or inflection points, so from there it should
//! just be splitting at max distance

use super::*;
use crate::common::{Line, PCoords, Segment, dist, linear_space};

const N_INTR: usize = 6;

#[derive(Debug, Clone)]
struct IntervalData<const D: usize> {
    /// Parameter at the beginning of the interval
    t0: f64,

    /// Point and direction at the beginning of the interval,
    line0: Line<D>,

    /// A guaranteed upper bound on the deviation of the spline from the base leg over the interval.
    /// This is the radius of the capsule bounding volume, so it must never under-estimate the true
    /// deviation or the branch-and-bound pruning could discard the interval that actually holds the
    /// closest point. See [`interval_error_bound`].
    e_max: f64,
}

/// This function finds the closest and farthest possible distances between a test point and the
/// capsule shape formed by the line segment of an interval and the radius of its maximum known
/// error. Returns the closest distance as the first element in the tuple, and the farthest distance
/// as the second.
fn dist_min_max<const D: usize>(
    test: &impl PCoords<D>,
    data: &IntervalData<D>,
    end: &Point<f64, D>,
) -> (f64, f64) {
    let dir = end - data.line0.origin;

    // First, we find the shortest distance to the line segment, which is the distance to the
    // projection of the test point onto the segment clamped to its bounds.
    let t =
        ((test.coords() - data.line0.origin.coords).dot(&dir) / dir.norm_squared()).clamp(0.0, 1.0);
    let shortest = (test.coords() - (data.line0.origin + dir * t).coords).norm();

    // Now we find the largest distance to the line segment, which is the larger of the distances
    // to the end points.
    let da = test.coords() - data.line0.origin.coords;
    let db = test.coords() - end.coords;
    let longest = da.norm_squared().max(db.norm_squared()).sqrt();

    // Finally, we return the closest possible distance from the test point to the capsule and
    // the farthest possible distance from the test point to the capsule
    ((shortest - data.e_max).max(0.0), longest + data.e_max)
}

/// This function finds the closest distance between a test point and the region of the spline
/// bounded by `t0` and `t1`. This should only be used on sections of the spline that have no
/// cusps or inflections, and are convex and flat*
fn closest_to_point<const D: usize>(
    test: &impl PCoords<D>,
    spline: &CubicSpline<D>,
    t0: f64,
    t1: f64,
    l0: &Line<D>,
    l1: &Line<D>,
) -> SplineValue<f64> {
    // First we can check the point to figure out what Voronoi region it's in using the
    // direction of the interval ends
    let before_front = l0.scalar_project(test) < 0.0;
    let after_back = l1.scalar_project(test) > 0.0;
    match (before_front, after_back) {
        // The test point lies unambiguously before the front of the interval, so the closest point
        // is the very front
        (true, false) => SplineValue::new(t0, dist(test, &l0.origin)),

        // The test point lies unambiguously after the end of the interval, so the closest point
        // is the very back
        (false, true) => SplineValue::new(t1, dist(test, &l1.origin)),

        // The test point is both before the front _and_ after the back, which is possible when
        // it is on the concave side beyond the focus of the concavity. The closest point is either
        // the front or back.
        (true, true) => {
            let q0 = SplineValue::new(t0, dist(test, &l0.origin));
            let q1 = SplineValue::new(t1, dist(test, &l1.origin));
            if q0.value < q1.value { q0 } else { q1 }
        }

        // The test point is within the Voronoi region between the endpoints. If we're on the
        // concave side of the curvature, there may be up to two local distance minima, so we need
        // to identify both and then find the one which is closer to the test point.
        //
        // Seeding the iteration from each end biases it toward the minimum in that end's basin, so
        // when two minima are present (one near each endpoint) both get found and the closer is
        // returned. When there is only one, both seeds converge to it.
        (false, false) => {
            let q0 = iterate_to_closest(test, spline, t0, t0, t1);
            let q1 = iterate_to_closest(test, spline, t1, t0, t1);
            if q0.value < q1.value { q0 } else { q1 }
        }
    }
}

/// Iterates to a local distance minimum, starting from `t_start` and confined to `[lo, hi]`.
///
/// The minimum of the distance to the spline is a root of the footpoint function
/// `g(t) = (p(t) - test) · p'(t)`, whose own derivative is
/// `g'(t) = |p'(t)|² + (p(t) - test) · p''(t)`. Both are available analytically from the spline,
/// so the root is refined with Newton's method. Each step is clamped to `[lo, hi]`; if a step
/// would leave the bracket or the curvature `g'` is too small to trust, the step is replaced by a
/// bisection of the current sign bracket so the iteration cannot diverge.
fn iterate_to_closest<const D: usize>(
    test: &impl PCoords<D>,
    spline: &CubicSpline<D>,
    t_start: f64,
    lo: f64,
    hi: f64,
) -> SplineValue<f64> {
    const MAX_ITER: usize = 32;

    let tc = test.coords();
    let g = |t: f64| -> f64 { (spline.position(t).coords - tc).dot(&spline.derivative(t)) };

    // Maintain a sign bracket `[a, b]` with `g(a) <= 0 <= g(b)` as a fallback for Newton steps
    // that misbehave. It is only usable once we actually observe opposite signs at the ends.
    let mut a = lo;
    let mut b = hi;
    let have_bracket = g(lo) <= 0.0 && g(hi) >= 0.0;

    let mut t = t_start.clamp(lo, hi);
    for _ in 0..MAX_ITER {
        let r = spline.position(t).coords - tc;
        let d1 = spline.derivative(t);
        let gt = r.dot(&d1);
        if gt.abs() < 1e-12 {
            break;
        }

        // Tighten the fallback bracket using the sign of g at the current estimate, preserving
        // the `g(a) <= 0 <= g(b)` invariant.
        if have_bracket {
            if gt < 0.0 {
                a = t;
            } else {
                b = t;
            }
        }

        let gpt = d1.norm_squared() + r.dot(&spline.second_derivative(t));
        let newton = t - gt / gpt;
        let next = if gpt.abs() > 1e-14 && newton > lo && newton < hi {
            newton
        } else if have_bracket {
            0.5 * (a + b)
        } else {
            // No usable curvature and no bracket to fall back on; the estimate is as good as it
            // gets without more structure, so stop here.
            break;
        };

        let converged = (next - t).abs() < 1e-12;
        t = next;
        if converged {
            break;
        }
    }

    SplineValue::new(t, dist(test, &spline.position(t)))
}

/// A prebuilt acceleration structure for running repeated spatial queries against a
/// [`CubicSpline`].
///
/// Construction partitions the curve's parameter domain into a fixed number of intervals, each
/// wrapped in a capsule-shaped bounding volume (the chord between the interval's endpoints, plus a
/// guaranteed-conservative radius bounding the curve's deviation from that chord). Splits are
/// placed at structural features first (a cusp, then curvature zeros), so every interval is a
/// simple convex, inflection-free arc; the remaining budget is spent halving whichever interval
/// deviates most from its chord. Queries then use branch-and-bound over the capsules to discard
/// intervals that cannot contain the answer before solving the survivors exactly.
///
/// Build the structure once and reuse it across many queries. Construct it with
/// [`CubicSpline::into_query`] or the [`From`] impls (`CubicSplineQueries::from(spline)`).
#[derive(Debug, Clone)]
pub struct CubicSplineQueries<const D: usize> {
    /// The underlying spline that gets queried
    spline: CubicSpline<D>,

    /// The point and direction at the end of the spline
    line1: Line<D>,

    /// The intervals
    intervals: [IntervalData<D>; N_INTR],
}

impl<const D: usize> CubicSplineQueries<D> {
    /// Returns the closest point on the curve to `point` as a [`SplineValue<f64>`], holding the
    /// parameter `t` in `[0, 1]` and, as its value, the distance from `point` to that location.
    /// Recover the point itself with [`CubicSpline::position`] on the spline the structure was
    /// built from.
    ///
    /// This is a true global closest point over the whole curve: branch-and-bound over the interval
    /// bounding volumes guarantees no interval that could contain the minimum is discarded, so the
    /// result is correct even where the curve loops or crosses itself and a naive nearest-sample
    /// search would be fooled. When two points on the curve are equidistant from `point`, which of
    /// the tied projections is returned is unspecified.
    pub fn project_point(&self, point: &impl PCoords<D>) -> SplineValue<f64> {
        // To do the initial pruning, we're going to use the closest/farthest method. Each interval
        // has a capsule shaped bounding volume formed by its two endpoints and the known maximum
        // error value of the curve to the segment between the endpoints. For any test point, we
        // can compute the distance to the nearest point in the volume (minimum of 0.0) and the
        // distance to the farthest point in the volume.
        //
        // Of the bounding volumes, one of them will have the smallest value of the farthest
        // possible distance.  Any volume who does not have a _closest_ distance less than that
        // value cannot possibly contain the projection of the test point.
        let mut prune = [(f64::INFINITY, f64::INFINITY); N_INTR];
        for (i, slot) in prune.iter_mut().enumerate() {
            let end = if i < N_INTR - 1 {
                self.intervals[i + 1].line0.origin
            } else {
                self.line1.origin
            };
            *slot = dist_min_max(point, &self.intervals[i], &end);
        }
        let min_farthest = prune
            .iter()
            .map(|(_, far)| *far)
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        // Once the bounds have been found, we'll find the closest projection for each interval
        // that has a closest distance less than the minimum farthest
        let mut closest = [SplineValue::new(f64::NAN, f64::INFINITY); N_INTR];
        for i in 0..N_INTR {
            // `<=` (not `<`) so the interval that attains `min_farthest` is never pruned: its
            // closest bound can equal `min_farthest` in degenerate cases, and pruning every
            // interval would leave nothing to return.
            if prune[i].0 <= min_farthest {
                let (t1, line1) = if i < N_INTR - 1 {
                    (self.intervals[i + 1].t0, self.intervals[i + 1].line0)
                } else {
                    (1.0, self.line1)
                };
                closest[i] = closest_to_point(
                    point,
                    &self.spline,
                    self.intervals[i].t0,
                    t1,
                    &self.intervals[i].line0,
                    &line1,
                );
            }
        }

        // Finally, we return the projection with the smallest distance
        *closest
            .iter()
            .min_by(|a, b| a.value.partial_cmp(&b.value).unwrap())
            .unwrap()
    }

    /// Builds the fixed acceleration structure for `spline` by partitioning its domain into
    /// `N_INTR` intervals.
    ///
    /// The first splits are chosen in priority order: a cusp (if one exists), then any
    /// inflection/curvature-zero points (there are at most two), and then the parameter of
    /// maximum deviation from the base leg. Once those structural splits are exhausted, the
    /// interval array is filled by repeatedly splitting whichever interval currently has the
    /// largest deviation from its own local base leg, at the parameter where that deviation is
    /// greatest, until all `N_INTR` slots are used.
    pub(crate) fn new(spline: CubicSpline<D>) -> Self {
        let mut working = vec![WorkingInterval::new(&spline, 0.0, 1.0)];
        let line1 = spline.line_at(1.0);

        if let Some(tc) = spline.find_cusp() {
            Self::try_split_at(&spline, &mut working, tc);
        }

        for &t in spline.find_curvature_zeros().iter() {
            if t.is_finite() {
                Self::try_split_at(&spline, &mut working, t);
            }
        }

        while working.len() < N_INTR {
            let worst = working
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.e_max.partial_cmp(&b.e_max).unwrap())
                .map(|(i, _)| i)
                .expect("working always has at least one interval");
            let t_split = working[worst].t_max;
            Self::split_at_index(&spline, &mut working, worst, t_split);
        }

        let intervals: [IntervalData<D>; N_INTR] = working
            .into_iter()
            .map(|w| IntervalData {
                t0: w.t0,
                line0: spline.line_at(w.t0),
                e_max: interval_error_bound(&spline, w.t0, w.t1),
            })
            .collect::<Vec<_>>()
            .try_into()
            .unwrap_or_else(|_| panic!("working did not have exactly N_INTR intervals"));

        Self {
            spline,
            line1,
            intervals,
        }
    }

    /// If `t` falls strictly inside one of the `working` intervals (and there is still room left
    /// in the fixed-size array), splits that interval into two at `t`. Otherwise does nothing.
    fn try_split_at(spline: &CubicSpline<D>, working: &mut Vec<WorkingInterval>, t: f64) {
        if working.len() >= N_INTR {
            return;
        }

        let index = working
            .iter()
            .position(|w| t > w.t0 + 1e-12 && t < w.t1 - 1e-12);
        if let Some(index) = index {
            Self::split_at_index(spline, working, index, t);
        }
    }

    /// Splits the working interval at `index` into two at parameter `t`, replacing it in place
    /// with its left and right halves.
    fn split_at_index(
        spline: &CubicSpline<D>,
        working: &mut Vec<WorkingInterval>,
        index: usize,
        t: f64,
    ) {
        let old = working.remove(index);
        let left = WorkingInterval::new(spline, old.t0, t);
        let right = WorkingInterval::new(spline, t, old.t1);
        working.insert(index, right);
        working.insert(index, left);
    }

    pub fn into_spline(self) -> CubicSpline<D> {
        self.spline
    }
}

/// Builds the acceleration structure from a borrowed spline, cloning it internally. Use the
/// owning [`From<CubicSpline>`] impl to avoid the clone when the spline is no longer needed.
impl<const D: usize> From<&CubicSpline<D>> for CubicSplineQueries<D> {
    fn from(value: &CubicSpline<D>) -> Self {
        CubicSplineQueries::new(value.clone())
    }
}

/// Builds the acceleration structure from an owned spline, taking ownership of it. Equivalent to
/// [`CubicSpline::into_query`].
impl<const D: usize> From<CubicSpline<D>> for CubicSplineQueries<D> {
    fn from(value: CubicSpline<D>) -> Self {
        CubicSplineQueries::new(value)
    }
}

/// Bookkeeping for an in-progress interval during the construction of [`CubicSplineQueries`].
/// Unlike the final [`IntervalData`], this also tracks the interval's end parameter and the
/// location of its maximum error, since both are needed to decide where to split next.
#[derive(Debug, Clone, Copy)]
struct WorkingInterval {
    /// Parameter at the beginning of the interval
    t0: f64,

    /// Parameter at the end of the interval
    t1: f64,

    /// Parameter within `[t0, t1]` where the deviation from the interval's own base leg (the
    /// chord between the curve's position at `t0` and at `t1`) is greatest
    t_max: f64,

    /// The deviation from the base leg at `t_max`
    e_max: f64,
}

impl WorkingInterval {
    fn new<const D: usize>(spline: &CubicSpline<D>, t0: f64, t1: f64) -> Self {
        let (t_max, e_max) = interval_max_error(spline, t0, t1);
        Self {
            t0,
            t1,
            t_max,
            e_max,
        }
    }
}

/// Returns the sub-curve of `spline` corresponding to the parameter range `[t0, t1]`, with its
/// own parameterization re-based to `[0, 1]`, via two applications of de Casteljau splitting.
fn sub_curve<const D: usize>(spline: &CubicSpline<D>, t0: f64, t1: f64) -> CubicSpline<D> {
    let right = if t0 <= 0.0 {
        spline.clone()
    } else {
        spline.split(t0).1
    };

    if t1 >= 1.0 {
        right
    } else {
        let local_t1 = (t1 - t0) / (1.0 - t0);
        right.split(local_t1).0
    }
}

/// Finds the parameter (in the original spline's domain) and magnitude of the maximum deviation
/// of `spline` from the straight leg between its positions at `t0` and `t1`, restricted to the
/// sub-range `[t0, t1]`.
fn interval_max_error<const D: usize>(spline: &CubicSpline<D>, t0: f64, t1: f64) -> (f64, f64) {
    if t1 - t0 < 1e-12 {
        return (t0, 0.0);
    }

    let sub = sub_curve(spline, t0, t1);
    let helper = CubicSplineBaseDist::new(&sub);
    let local_t = helper.find_max_t(0.0, 1.0);
    let e_max = helper.e(local_t);
    (t0 + local_t * (t1 - t0), e_max)
}

/// Returns a guaranteed upper bound on the deviation of the sub-curve over `[t0, t1]` from its base
/// leg (the chord between the sub-curve's endpoints).
///
/// Unlike [`interval_max_error`], which uses a sampling-plus-Newton search that can under-estimate
/// the true peak, this is provably conservative and so is what the capsule bounding volume must use
/// as its radius. A Bézier curve lies inside the convex hull of its control points, and the
/// distance to the (convex) base segment is itself a convex function of position, so it attains its
/// maximum over the hull at a control point. The endpoints `P0` and `P3` lie on the base leg, so
/// only the inner control points `P1` and `P2` can contribute, and their larger distance to the leg
/// is an upper bound on the deviation of every point of the sub-curve.
fn interval_error_bound<const D: usize>(spline: &CubicSpline<D>, t0: f64, t1: f64) -> f64 {
    if t1 - t0 < 1e-12 {
        return 0.0;
    }

    let sub = sub_curve(spline, t0, t1);
    if dist(&sub.p0, &sub.p3) < 1e-12 {
        // The base leg has collapsed to a point, so measure the inner control points from it.
        return dist(&sub.p1, &sub.p0).max(dist(&sub.p2, &sub.p0));
    }

    let seg = Segment::new_unchecked(sub.p0, sub.p3);
    let d1 = dist(&seg.closest_point(&sub.p1), &sub.p1);
    let d2 = dist(&seg.closest_point(&sub.p2), &sub.p2);
    d1.max(d2)
}

/// This struct is a helper to find the distance between a cubic spline of arbitrary dimension D
/// and the straight line segment running between its endpoints.  This is used for finding the
/// maximum distances for partitioning.
pub(crate) struct CubicSplineBaseDist<'a, const D: usize> {
    spline: &'a CubicSpline<D>,
    seg: Option<Segment<D>>,
}

impl<'a, const D: usize> CubicSplineBaseDist<'a, D> {
    pub fn new(spline: &'a CubicSpline<D>) -> Self {
        Self {
            spline,
            seg: if dist(&spline.p0, &spline.p3) < 1e-12 {
                None
            } else {
                Some(Segment::new_unchecked(spline.p0, spline.p3))
            },
        }
    }

    /// Get the error distance from the spline at parameter `t` to the nearest point on the base
    /// segment running between `p0` and `p3`
    pub fn e(&self, t: f64) -> f64 {
        let p = self.spline.position(t);
        let cp = if let Some(seg) = self.seg {
            seg.closest_point(&p)
        } else {
            self.spline.p0
        };

        dist(&cp, &p)
    }

    /// Get the first derivative of the error distance `e` at the parameter `t`
    pub fn dedt(&self, t: f64) -> f64 {
        let p = self.spline.position(t);
        let cp = if let Some(seg) = self.seg {
            seg.closest_point(&p)
        } else {
            self.spline.p0
        };

        let n = (p - cp).normalize();
        n.dot(&self.spline.derivative(t))
    }

    /// Finds the parameter `t` in `[t0, t1]` (inclusive) at which the error distance `e` between
    /// the spline and its base segment is maximized.
    ///
    /// The search begins by evaluating `e` at five evenly spaced parameters spanning the
    /// interval (including both endpoints). The largest of those samples is then refined with
    /// Newton's method applied to the root of `dedt`, the slope of which is estimated with a
    /// central difference. If the refinement fails to improve on the best sample (for example
    /// because it diverged or converged to a minimum instead of a maximum), the best sample is
    /// returned instead.
    ///
    /// When `e` is (numerically) flat across the interval, e.g. because the spline happens to run
    /// exactly along its own base leg, every sample ties and the middle sample is kept as the
    /// default. This guarantees a midpoint split rather than a degenerate, zero-width one when
    /// this is used to choose where to subdivide an interval.
    pub(crate) fn find_max_t(&self, t0: f64, t1: f64) -> f64 {
        const N_SAMPLES: usize = 5;
        const MAX_ITER: usize = 20;

        let samples = linear_space(t0, t1, N_SAMPLES);
        let mut best_t = samples[N_SAMPLES / 2];
        let mut best_e = self.e(best_t);
        for &t in samples.values() {
            let e = self.e(t);
            if e > best_e {
                best_e = e;
                best_t = t;
            }
        }

        let h = ((t1 - t0).abs() * 1e-4).max(1e-9);
        let mut t = best_t;
        for _ in 0..MAX_ITER {
            let d1 = self.dedt(t);
            if d1.abs() < 1e-12 {
                break;
            }

            let t_minus = (t - h).max(t0);
            let t_plus = (t + h).min(t1);
            if t_plus - t_minus < 1e-15 {
                break;
            }
            let d2 = (self.dedt(t_plus) - self.dedt(t_minus)) / (t_plus - t_minus);
            if d2.abs() < 1e-12 {
                break;
            }

            let next_t = (t - d1 / d2).clamp(t0, t1);
            let converged = (next_t - t).abs() < 1e-12;
            t = next_t;
            if converged {
                break;
            }
        }

        let e_t = self.e(t);
        if e_t > best_e { t } else { best_t }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point2;
    use approx::assert_relative_eq;

    fn sample_2d() -> CubicSpline<2> {
        CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 1.0),
            Point2::new(3.0, 0.0),
        )
    }

    #[test]
    fn find_max_t_locates_symmetric_peak() {
        // sample_2d bulges symmetrically above its base chord, so the maximum error should sit
        // at the midpoint of the parameter range.
        let c = sample_2d();
        let helper = CubicSplineBaseDist::new(&c);
        let t = helper.find_max_t(0.0, 1.0);
        assert_relative_eq!(t, 0.5, epsilon = 1e-6);
    }

    #[test]
    fn find_max_t_matches_dense_sampling() {
        // Cross-check against a dense brute-force scan over the interval.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 2.0),
            Point2::new(2.0, -3.0),
            Point2::new(3.0, 0.0),
        );
        let helper = CubicSplineBaseDist::new(&c);
        let t = helper.find_max_t(0.0, 1.0);

        let mut best_t = 0.0;
        let mut best_e = f64::NEG_INFINITY;
        for i in 0..=10_000 {
            let s = i as f64 / 10_000.0;
            let e = helper.e(s);
            if e > best_e {
                best_e = e;
                best_t = s;
            }
        }

        assert_relative_eq!(t, best_t, epsilon = 1e-3);
        assert_relative_eq!(helper.e(t), best_e, epsilon = 1e-6);
    }

    #[test]
    fn find_max_t_restricted_to_sub_interval() {
        // Restricting the search window should confine the result to that window even when the
        // global maximum lies outside it.
        let c = sample_2d();
        let helper = CubicSplineBaseDist::new(&c);
        let t = helper.find_max_t(0.6, 1.0);
        assert!((0.6 - 1e-9..=1.0 + 1e-9).contains(&t));
        // Error is monotonically decreasing from the t=0.5 peak out to t=1, so the restricted
        // maximum should sit at the left edge of the window.
        assert_relative_eq!(t, 0.6, epsilon = 1e-6);
    }

    #[test]
    fn find_max_t_degenerate_interval_returns_endpoint() {
        let c = sample_2d();
        let helper = CubicSplineBaseDist::new(&c);
        let t = helper.find_max_t(0.3, 0.3);
        assert_relative_eq!(t, 0.3, epsilon = 1e-12);
    }

    /// Checks that the `t0` values of the built intervals are strictly ascending, start at 0.0,
    /// and never reach 1.0 (since the last interval implicitly runs to the end of the spline).
    fn assert_well_formed<const D: usize>(q: &CubicSplineQueries<D>) {
        assert_eq!(q.intervals.len(), N_INTR);
        assert_relative_eq!(q.intervals[0].t0, 0.0, epsilon = 1e-12);
        for w in q.intervals.windows(2) {
            assert!(
                w[1].t0 > w[0].t0,
                "t0 values must be strictly ascending: {} then {}",
                w[0].t0,
                w[1].t0
            );
        }
        assert!(q.intervals[N_INTR - 1].t0 < 1.0);
    }

    #[test]
    fn new_builds_n_intr_well_formed_intervals() {
        let q = CubicSplineQueries::new(sample_2d());
        assert_well_formed(&q);
    }

    #[test]
    fn new_splits_at_a_cusp() {
        // p1 == p2 with mirrored arms produces a cusp at exactly t = 0.5.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 0.0),
        );
        let cusp = c.find_cusp().expect("expected a cusp");
        let q = CubicSplineQueries::new(c);
        assert_well_formed(&q);
        assert!(
            q.intervals.iter().any(|iv| (iv.t0 - cusp).abs() < 1e-9),
            "expected an interval boundary at the cusp t = {}",
            cusp
        );
    }

    #[test]
    fn new_splits_at_an_inflection() {
        // A symmetric S-curve with a single inflection at t = 0.5.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, -1.0),
            Point2::new(3.0, 0.0),
        );
        let infl = c.find_inflections()[0];
        assert!(infl.is_finite());
        let q = CubicSplineQueries::new(c);
        assert_well_formed(&q);
        assert!(
            q.intervals.iter().any(|iv| (iv.t0 - infl).abs() < 1e-9),
            "expected an interval boundary at the inflection t = {}",
            infl
        );
    }

    #[test]
    fn new_intervals_have_bounded_error() {
        // Once the fixed budget of intervals is spent splitting at points of maximum deviation,
        // each interval's local error should be small relative to the deviation of the whole
        // curve from its single base chord.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 3.0),
            Point2::new(2.5, -3.0),
            Point2::new(3.0, 0.0),
        );
        let whole_error = {
            let helper = CubicSplineBaseDist::new(&c);
            let t = helper.find_max_t(0.0, 1.0);
            helper.e(t)
        };

        let q = CubicSplineQueries::new(c);
        assert_well_formed(&q);
        for iv in &q.intervals {
            assert!(
                iv.e_max <= whole_error + 1e-9,
                "interval error {} exceeded whole-curve error {}",
                iv.e_max,
                whole_error
            );
        }
        // The worst remaining interval error should have shrunk substantially from subdividing.
        let worst = q.intervals.iter().map(|iv| iv.e_max).fold(0.0, f64::max);
        assert!(
            worst < whole_error * 0.5,
            "worst interval error {} did not shrink enough from whole-curve error {}",
            worst,
            whole_error
        );
    }

    #[test]
    fn new_works_in_3d() {
        use crate::Point3;
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.5),
            Point3::new(2.0, -1.0, -0.5),
            Point3::new(3.0, 0.0, 0.0),
        );
        let q = CubicSplineQueries::new(c);
        assert_well_formed(&q);
    }

    /// Brute-force closest parameter over a dense scan, used to cross-check `project_point`.
    fn brute_closest<const D: usize>(spline: &CubicSpline<D>, p: &Point<f64, D>) -> f64 {
        let mut best_t = 0.0;
        let mut best_d = f64::INFINITY;
        for i in 0..=20_000 {
            let t = i as f64 / 20_000.0;
            let d = dist(p, &spline.position(t));
            if d < best_d {
                best_d = d;
                best_t = t;
            }
        }
        best_t
    }

    fn assert_projects_like_brute(spline: &CubicSpline<2>, p: Point<f64, 2>) {
        let q = CubicSplineQueries::new(spline.clone());
        let proj = q.project_point(&p);
        let t_ref = brute_closest(spline, &p);
        // Compare on distance rather than parameter, since near-flat regions can have very
        // different parameters at nearly identical distances.
        let d = dist(&p, &spline.position(proj.t));
        let d_ref = dist(&p, &spline.position(t_ref));
        // The returned distance must match the distance to the returned parameter.
        assert_relative_eq!(proj.value, d, epsilon = 1e-9);
        assert!(
            d <= d_ref + 1e-6,
            "project_point distance {} worse than brute-force {} (t={}, t_ref={})",
            d,
            d_ref,
            proj.t,
            t_ref
        );
    }

    #[test]
    fn interval_e_max_is_a_true_upper_bound() {
        // The capsule radius must never under-estimate the curve's deviation from its base leg,
        // or pruning could discard the interval holding the closest point. Check it against a
        // dense scan of the actual deviation for several shapes, including a high-deflection curve.
        let curves = [
            sample_2d(),
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(0.5, 3.0),
                Point2::new(2.5, -3.0),
                Point2::new(3.0, 0.0),
            ),
            CubicSpline::new(
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 2.0),
                Point2::new(2.0, -2.0),
                Point2::new(3.0, 0.0),
            ),
        ];

        for c in &curves {
            let q = CubicSplineQueries::new(c.clone());
            for (i, iv) in q.intervals.iter().enumerate() {
                let t0 = iv.t0;
                let t1 = if i < N_INTR - 1 {
                    q.intervals[i + 1].t0
                } else {
                    1.0
                };
                let p0 = c.position(t0);
                let p1 = c.position(t1);
                let leg = Segment::new_unchecked(p0, p1);
                let mut worst = 0.0_f64;
                for k in 0..=2000 {
                    let t = t0 + (t1 - t0) * (k as f64 / 2000.0);
                    let p = c.position(t);
                    worst = worst.max(dist(&leg.closest_point(&p), &p));
                }
                assert!(
                    iv.e_max >= worst - 1e-9,
                    "interval {} radius {} under-bounds observed deviation {}",
                    i,
                    iv.e_max,
                    worst
                );
            }
        }
    }

    #[test]
    fn project_point_matches_brute_force_arch() {
        let c = sample_2d();
        // A spread of query points: above, below, off each end, and near the curve.
        for &p in &[
            Point2::new(1.5, 5.0),
            Point2::new(1.5, -3.0),
            Point2::new(-2.0, 0.0),
            Point2::new(5.0, 0.0),
            Point2::new(1.5, 0.75),
            Point2::new(0.5, 0.5),
        ] {
            assert_projects_like_brute(&c, p);
        }
    }

    #[test]
    fn project_point_matches_brute_force_s_curve() {
        // An S-curve exercises the inflection split and concave-side queries.
        let c = CubicSpline::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(2.0, -2.0),
            Point2::new(3.0, 0.0),
        );
        for &p in &[
            Point2::new(0.5, 2.0),
            Point2::new(2.5, -2.0),
            Point2::new(1.5, 0.0),
            Point2::new(1.5, 3.0),
            Point2::new(1.5, -3.0),
            Point2::new(-1.0, 1.0),
            Point2::new(4.0, -1.0),
        ] {
            assert_projects_like_brute(&c, p);
        }
    }
}
