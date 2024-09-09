//! Common operations on f64 points in D-dimensional space.

use parry3d_f64::na::{Point, SVector, Isometry, AbstractRotation};
use crate::common::surface_point::SurfacePoint;


/// Returns the distance between two points in D-dimensional space.
///
/// # Arguments
///
/// * `a`: the first point
/// * `b`: the second point
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::points::dist;
/// use engeom::Point2;
/// let a = Point2::new(1.0, 2.0);
/// let b = Point2::new(3.0, 2.0);
/// let d = dist(&a, &b);
/// assert_eq!(d, 2.0);
/// ```
pub fn dist<const D: usize>(a: &Point<f64, D>, b: &Point<f64, D>) -> f64 {
    (a - b).norm()
}

/// Returns the midpoint between two points in D-dimensional space.
///
/// # Arguments
///
/// * `a`: the first point
/// * `b`: the second point
///
/// returns: OPoint<f64, Const<{ D }>>
///
/// # Examples
///
/// ```
/// use engeom::common::points::mid_point;
/// use engeom::Point2;
/// let a = Point2::new(1.0, 2.0);
/// let b = Point2::new(3.0, 4.0);
/// let mid = mid_point(&a, &b);
/// assert_eq!(mid, Point2::new(2.0, 3.0));
/// ```
pub fn mid_point<const D: usize>(a: &Point<f64, D>, b: &Point<f64, D>) -> Point<f64, D> {
    b + (a - b) * 0.5
}

/// Returns the mean point of a set of points in D-dimensional space.  The mean is found by summing
/// the coordinates of all points and dividing by the number of points.  There is no weighting of
/// the points in this calculation.
///
/// # Arguments
///
/// * `points`: a slice of points to compute the mean of
///
/// returns: OPoint<f64, Const<{ D }>>
///
/// # Examples
///
/// ```
/// use engeom::common::points::mean_point;
/// use engeom::Point2;
/// let points = vec![Point2::new(1.0, 2.0), Point2::new(3.0, 4.0), Point2::new(5.0, 6.0)];
/// let mean = mean_point(&points);
/// assert_eq!(mean, Point2::new(3.0, 4.0));
/// ```
pub fn mean_point<const D: usize>(points: &[Point<f64, D>]) -> Point<f64, D> {
    let mut sum = SVector::<f64, D>::zeros();
    for p in points {
        sum += p.coords;
    }
    Point::<f64, D>::from(sum / points.len() as f64)
}

/// Computes the weighted mean point of a set of points in D-dimensional space
///
/// # Arguments
///
/// * `points`: the points to compute the weighted mean of
/// * `weights`: a slice of weights for each point. The length of this slice must be the same as
/// the length of the `points` slice
///
/// returns: OPoint<f64, Const<{ D }>>
pub fn mean_point_weighted<const D: usize>(
    points: &[Point<f64, D>],
    weights: &[f64],
) -> Point<f64, D> {
    let mut sum = SVector::<f64, D>::zeros();
    let mut total_weight = 0.0;
    for (p, w) in points.iter().zip(weights) {
        sum += p.coords * *w;
        total_weight += *w;
    }
    Point::<f64, D>::from(sum / total_weight)
}

/// Perform the Ramer-Douglas-Peucker algorithm on a set of points in D-dimensional space.  The
/// algorithm simplifies a curve by reducing the number of points while preserving the shape of the
/// curve.  The `tol` parameter is the maximum distance from the simplified curve to the original
/// curve.
///
/// # Arguments
///
/// * `points`: The ordered set of points to simplify
/// * `tol`: The maximum distance from the simplified curve to the original curve
///
/// returns: Vec<OPoint<f64, Const<{ D }>>, Global>
///
/// # Examples
///
/// ```
///
/// ```
pub fn ramer_douglas_peucker<const D: usize>(
    points: &[Point<f64, D>],
    tol: f64,
) -> Vec<Point<f64, D>> {
    let mut rdp = Rdp::new(points, tol);
    rdp.simplify(0, points.len() - 1);
    rdp.generate_points()
}


/// A struct for handling the state of the Ramer-Douglas-Peucker algorithm.
struct Rdp<'a, const D: usize> {
    points: &'a [Point<f64, D>],
    keep: Vec<bool>,
    tol: f64,
}

impl<'a, const D: usize> Rdp<'a, D> {
    fn new(points: &'a [Point<f64, D>], tol: f64) -> Self {
        let keep = points.iter().map(|_| false).collect();
        Rdp { points, keep, tol }
    }

    fn simplify(&mut self, i0: usize, i1: usize) {
        self.keep[i0] = true;
        self.keep[i1] = true;
        if i1 - i0 < 2 {
            return;
        }

        let sp = SurfacePoint::new_normalize(self.points[i0], self.points[i1] - self.points[i0]);
        let mut max_dist = 0.0;
        let mut max_i = 0;

        for i in i0 + 1..i1 {
            let dist = (sp.projection(&self.points[i]) - self.points[i]).norm();
            if dist > max_dist {
                max_dist = dist;
                max_i = i;
            }
        }

        if max_dist > self.tol {
            self.simplify(i0, max_i);
            self.simplify(max_i, i1);
        }
    }

    fn generate_points(&self) -> Vec<Point<f64, D>> {
        let mut result = Vec::new();
        for i in 0..self.points.len() {
            if self.keep[i] {
                result.push(self.points[i]);
            }
        }
        result
    }
}

/// Generic 2 or 3 dimensional transformation of a slice of `Point` entities by an `Isometry`,
/// resulting in an owned `Vec` of new point entities being created and returned.
///
/// TODO: is parry's `transformed(...)` method acceptable to use here instead?
///
/// # Arguments
///
/// * `points`: a slice of `Point` entities to transform
/// * `transform`: the `Isometry` to apply to each point
///
/// returns: Vec<OPoint<f64, Const<{ D }>>, Global>
///
/// # Examples
///
/// ```
/// use engeom::{Point2, Iso2};
/// use engeom::common::points::transform_points;
///
/// let points = vec![Point2::new(1.0, 2.0), Point2::new(3.0, 4.0)];
/// let transform = Iso2::translation(1.0, 2.0);
/// let transformed_points = transform_points(&points, &transform);
/// assert_eq!(transformed_points[0], Point2::new(2.0, 4.0));
/// assert_eq!(transformed_points[1], Point2::new(4.0, 6.0));
/// ```
pub fn transform_points<R, const D: usize>(
    points: &[Point<f64, D>],
    transform: &Isometry<f64, R, D>,
) -> Vec<Point<f64, D>>
where
    R: AbstractRotation<f64, D>,
{
    points.iter().map(|p| transform * p).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{Curve2, Point2};
    use approx::assert_relative_eq;

    #[test]
    fn test_simple_reduce() {
        let mut c = curve_from(&[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]);
        c = insert(&mut c, 0.25);
        c = insert(&mut c, 0.75);
        let points = c.clone_points();
        assert_eq!(points.len(), 5);
        let reduced = ramer_douglas_peucker(&points, 0.001);
        assert_eq!(reduced.len(), 3);
        assert_relative_eq!(reduced[0], Point2::new(0.0, 0.0));
        assert_relative_eq!(reduced[1], Point2::new(1.0, 0.0));
        assert_relative_eq!(reduced[2], Point2::new(1.0, 1.0));
    }

    fn curve_from(points: &[(f64, f64)]) -> Curve2 {
        let points = points
            .iter()
            .map(|(x, y)| Point2::new(*x, *y))
            .collect::<Vec<_>>();
        Curve2::from_points(&points, 0.001, false).unwrap()
    }

    fn insert(c: &mut Curve2, l: f64) -> Curve2 {
        let points = c.clone_points();
        let lengths = c.lengths();
        let mut working = points
            .iter()
            .zip(lengths.iter())
            .map(|(p, l)| (*l, *p))
            .collect::<Vec<_>>();
        working.push((l, c.at_length(l).unwrap().point()));
        working.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let points = working.iter().map(|(_, p)| *p).collect::<Vec<_>>();
        Curve2::from_points(&points, 0.001, false).unwrap()
    }
}
