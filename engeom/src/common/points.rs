//! Common operations on f64 points in D-dimensional space.

use crate::common::PCoords;
use crate::common::surface_point::SurfacePoint;
use parry3d_f64::na::{AbstractRotation, Isometry, Point, SVector};
use crate::common::kd_tree::{KdTree, KdTreeSearch};

/// Calculate the barycentric coordinates of a point `p` with respect to the triangle defined by
/// points `a`, `b`, and `c` in D-dimensional space.  The barycentric coordinates are a way of
/// expressing the position of a point within a triangle as a combination of the triangle's
/// vertices.
///
/// This function was taken from the book "Real-Time Collision Detection" by Christer Ericson.
///
/// # Arguments
///
/// * `a`: the first vertex of the triangle
/// * `b`: the second vertex of the triangle
/// * `c`: the third vertex of the triangle
/// * `p`: the point for which to calculate the barycentric coordinates
///
/// returns: [f64; 3]
pub fn barycentric<const D: usize>(
    a: &impl PCoords<D>,
    b: &impl PCoords<D>,
    c: &impl PCoords<D>,
    p: &impl PCoords<D>,
) -> [f64; 3] {
    let v0 = b.coords() - a.coords();
    let v1 = c.coords() - a.coords();
    let v2 = p.coords() - a.coords();

    let d00 = v0.dot(&v0);
    let d01 = v0.dot(&v1);
    let d11 = v1.dot(&v1);
    let d20 = v2.dot(&v0);

    let d21 = v2.dot(&v1);
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() < f64::EPSILON {
        [0.0, 0.0, 0.0]
    } else {
        let v = (d11 * d20 - d01 * d21) / denom;
        let w = (d00 * d21 - d01 * d20) / denom;
        let u = 1.0 - v - w;
        [u, v, w]
    }
}

/// Measure the angle between `pc`->`p1` and `pc`->`p2` in D-dimensional space.
///
/// # Arguments
///
/// * `pc`: The common point from which the angle is measured.
/// * `p1`: The first point to measure the angle to.
/// * `p2`: The second point to measure the angle to.
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::Point2;
/// use engeom::common::points::three_point_angle;
/// let pc = Point2::new(0.0, 0.0);
/// let p1 = Point2::new(1.0, 0.0);
/// let p2 = Point2::new(0.0, 1.0);
/// let angle = three_point_angle(&pc, &p1, &p2);
///
/// assert_relative_eq!(angle, std::f64::consts::FRAC_PI_2, epsilon = 1e-6);
/// ```
pub fn three_point_angle<const D: usize>(
    pc: &impl PCoords<D>,
    p1: &impl PCoords<D>,
    p2: &impl PCoords<D>,
) -> f64 {
    let v1 = p1.coords() - pc.coords();
    let v2 = p2.coords() - pc.coords();
    v1.angle(&v2)
}

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
pub fn dist<const D: usize>(a: &impl PCoords<D>, b: &impl PCoords<D>) -> f64 {
    (a.coords() - b.coords()).norm()
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
pub fn mid_point<const D: usize>(a: &impl PCoords<D>, b: &impl PCoords<D>) -> Point<f64, D> {
    Point::from(b.coords() + (a.coords() - b.coords()) * 0.5)
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
pub fn mean_point<const D: usize>(points: &[impl PCoords<D>]) -> Point<f64, D> {
    let mut sum = SVector::<f64, D>::zeros();
    for p in points {
        sum += p.coords();
    }
    Point::from(sum / points.len() as f64)
}

/// Computes the weighted mean point of a set of points in D-dimensional space
///
/// # Arguments
///
/// * `points`: the points to compute the weighted mean of
/// * `weights`: a slice of weights for each point. The length of this slice must be the same as
///   the length of the `points` slice
///
/// returns: OPoint<f64, Const<{ D }>>
pub fn mean_point_weighted<const D: usize>(
    points: &[impl PCoords<D>],
    weights: &[f64],
) -> Point<f64, D> {
    let mut sum = SVector::<f64, D>::zeros();
    let mut total_weight = 0.0;
    for (p, w) in points.iter().zip(weights) {
        sum += p.coords() * *w;
        total_weight += *w;
    }
    Point::from(sum / total_weight)
}

/// Produces a new set of points by evenly spacing points between `start` and `end` in
/// D-dimensional space.  The number of points to generate is specified by `num_points`.  The
/// start and end points are included in the result.
///
/// # Arguments
///
/// * `start`: the starting point (will be included in the result)
/// * `end`: the ending point (will be included in the result)
/// * `num_points`: the total number of points to generate
///
/// returns: Vec<OPoint<f64, Const<{ D }>>, Global>
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::common::points::evenly_spaced_points;
/// use engeom::Point2;
/// let start = Point2::new(0.0, 0.0);
/// let end = Point2::new(2.0, 0.0);
/// let points = evenly_spaced_points(&start, &end, 3);
///
/// assert_eq!(points.len(), 3);
/// assert_relative_eq!(points[0], Point2::new(0.0, 0.0));
/// assert_relative_eq!(points[1], Point2::new(1.0, 0.0));
/// assert_relative_eq!(points[2], Point2::new(2.0, 0.0));
/// ```
pub fn evenly_spaced_points<const D: usize>(
    start: &impl PCoords<D>,
    end: &impl PCoords<D>,
    num_points: usize,
) -> Vec<Point<f64, D>> {
    let mut result = Vec::new();
    let step = (end.coords() - start.coords()) / (num_points - 1) as f64;
    for i in 0..num_points {
        result.push((start.coords() + step * i as f64).into());
    }
    result
}

/// Generate a new set of points by evenly spacing points between `start` and `end` in D-dimensional
/// space.  The number of points to generate is specified by `num_points`.  The start and end points
/// are **not** included in the result.
///
/// # Arguments
///
/// * `start`: the starting point (will not be included in the result)
/// * `end`: the ending point (will not be included in the result)
/// * `num_points`: the total number of points to generate
///
/// returns: Vec<OPoint<f64, Const<{ D }>>, Global>
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::common::points::evenly_spaced_points_between;
/// use engeom::Point2;
/// let start = Point2::new(0.0, 0.0);
/// let end = Point2::new(2.0, 0.0);
/// let points = evenly_spaced_points_between(&start, &end, 3);
///
/// assert_eq!(points.len(), 3);
/// assert_relative_eq!(points[0], Point2::new(0.5, 0.0));
/// assert_relative_eq!(points[1], Point2::new(1.0, 0.0));
/// assert_relative_eq!(points[2], Point2::new(1.5, 0.0));
/// ```
pub fn evenly_spaced_points_between<const D: usize>(
    start: &impl PCoords<D>,
    end: &impl PCoords<D>,
    num_points: usize,
) -> Vec<Point<f64, D>> {
    let mut result = Vec::new();
    let step = (end.coords() - start.coords()) / (num_points + 1) as f64;
    for i in 1..num_points + 1 {
        result.push((start.coords() + step * i as f64).into());
    }

    result
}

/// Creates a new vec of points from a slice of points by filling in gaps between points where the
/// distance between them is greater than `max_dist`.  An evenly spaced set of points is inserted
/// between the two points to satisfy the `max_dist` threshold.
///
/// This is the opposite of a simplification algorithm, like Ramer-Douglas-Peucker, and can be used
/// to precondition a set of points before using an algorithm which relies on proximity.
///
/// # Arguments
///
/// * `original`: the original slice of points to fill gaps in.
/// * `max_dist`: the maximum distance between points in the new set of points.
///
/// returns: Vec<OPoint<f64, Const<{ D }>>, Global>
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::common::points::fill_gaps;
/// use engeom::Point2;
/// let points = vec![Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)];
/// let filled = fill_gaps(&points, 1.5);
///
/// assert_eq!(filled.len(), 3);
/// assert_relative_eq!(filled[0], Point2::new(0.0, 0.0));
/// assert_relative_eq!(filled[1], Point2::new(1.0, 0.0));
/// assert_relative_eq!(filled[2], Point2::new(2.0, 0.0));
/// ```
pub fn fill_gaps<const D: usize>(original: &[Point<f64, D>], max_dist: f64) -> Vec<Point<f64, D>> {
    // We fill gaps by iterating through the points, and if the distance between the current point
    // and its predecessor is greater than some threshold, we compute the number of points that
    // should be inserted (evenly spaced) between the two points to satisfy the `max_dist`
    // threshold.
    if original.len() < 2 {
        return original.to_vec();
    }

    let mut result = vec![original[0]];

    for p in original.iter().skip(1) {
        let d = dist(p, result.last().unwrap());
        if d > max_dist {
            let mut n = 1;
            while d / (n + 1) as f64 > max_dist {
                n += 1;
            }
            for x in evenly_spaced_points_between(result.last().unwrap(), p, n) {
                result.push(x);
            }
        }
        result.push(*p);
    }

    result
}

/// Returns the point in a set of points which has the largest value when projected onto a vector.
/// This is effectively the same as finding the point that is the furthest in the direction of the
/// vector.
///
/// There is nothing assumed about the list of points, so it is an O(n) operation where n is the
/// number of points in the set.
///
/// # Arguments
///
/// * `points`: the set of points to search
/// * `vector`: the vector to project the points onto
///
/// returns: Option<(usize, OPoint<f64, Const<{ D }>>)>
///
/// # Examples
///
/// ```
/// use engeom::common::points::max_point_in_direction;
/// use engeom::{Point2, Vector2};
/// use approx::assert_relative_eq;
///
/// let dir = Vector2::new(-10.0, -10.0);
///
/// let points = vec![
///    Point2::new(10.0, 0.0),
///    Point2::new(11.0, 1.0),
///    Point2::new(12.0, 2.0),
///    Point2::new(13.0, 3.0),
/// ];
///
/// let (_, max) = max_point_in_direction(&points, &dir).unwrap();
/// assert_relative_eq!(max, Point2::new(10.0, 0.0));
/// ```
pub fn max_point_in_direction<const D: usize>(
    points: &[impl PCoords<D>],
    vector: &SVector<f64, D>,
) -> Option<(usize, Point<f64, D>)> {
    let mut max_dist = f64::MIN;
    let mut max_i = None;

    for (i, p) in points.iter().enumerate() {
        let dist = p.coords().dot(vector);
        if dist > max_dist {
            max_dist = dist;
            max_i = Some(i);
        }
    }

    max_i.map(|max_i| (max_i, points[max_i].coords().into()))
}

/// Compute the error of a linear interpolation between two points `p0` and `p1` with respect to a
/// test point `p_test`.  The error is the distance between the test point and the projection of the
/// test point onto the line defined by `p0` and `p1`.  This can be thought of as the error that
/// would exist if the test point did not exist.
///
/// # Arguments
///
/// * `p0`: the first point of the linear interpolation
/// * `p1`: the second point of the linear interpolation
/// * `p_test`: a test point to project onto the line defined by `p0` and `p1`
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::Point2;
/// use engeom::common::points::linear_interpolation_error;
///
/// let p0 = Point2::new(0.0, 0.0);
/// let p1 = Point2::new(2.0, 0.0);
/// let p_test = Point2::new(1.0, 1.0);
///
/// let error = linear_interpolation_error(&p0, &p1, &p_test);
/// assert_relative_eq!(error, 1.0);
/// ```
pub fn linear_interpolation_error<const D: usize>(
    p0: &impl PCoords<D>,
    p1: &impl PCoords<D>,
    p_test: &Point<f64, D>,
) -> f64 {
    let sp = SurfacePoint::new_normalize(p0.coords().into(), p1.coords() - p0.coords());
    let proj = sp.projection(p_test);

    dist(p_test, &proj)
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


/// Clusters a set of points by distance tolerance, where all points reachable by a series of hops
/// less than or equal to the specified tolerance are grouped together into a single cluster.
///
/// This implementation does not make use of any spatial data structures, so it is O(n^2) in the
/// number of points.  Use it for relatively small numbers of points.
///
/// # Arguments
///
/// * `points`:
/// * `tol`:
///
/// returns: Vec<Vec<OPoint<f64, Const<{ D }>>, Global>, Global>
pub fn cluster_points_by_tol<const D: usize>(
    points: &[Point<f64, D>],
    tol: f64,
) -> Vec<Vec<Point<f64, D>>> {
    let mut clusters = Vec::new();
    let mut working = points.to_vec();
    let tol2 = tol * tol;

    while !working.is_empty() {
        let mut group = vec![working.pop().unwrap()];
        loop {
            let tree = KdTree::new(&group);
            let mut to_add = Vec::new();
            for (i, p) in working.iter().enumerate() {
                // let (d, _) = tree.nearest_one(&to_slice(p), &squared_euclidean);
                let (_, d) = tree.nearest_one(p);
                if d <= tol2 {
                    to_add.push(i);
                }
            }

            if to_add.is_empty() {
                break;
            } else {
                for i in to_add.iter().rev() {
                    group.push(working.remove(*i));
                }
            }
        }
        clusters.push(group);
    }

    clusters
}

pub fn merge_points_by_tol<const D: usize>(
    points: &[Point<f64, D>],
    tol: f64,
) -> Vec<Point<f64, D>> {
    let clusters = cluster_points_by_tol(points, tol);
    clusters.iter().map(|c| mean_point(c)).collect()
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector2;
    use crate::geom2::{Curve2, Point2};
    use approx::assert_relative_eq;
    use rand::Rng;

    #[test]
    fn distance_calc() {
        let p0 = Point2::new(0.0, 0.0);
        let p1 = Point2::new(3.0, 4.0);
        let d = dist(&p0, &p1);
        assert_eq!(d, 5.0); // 3-4-5 triangle
    }

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

    #[test]
    fn stress_barycentric_round_trip() {
        let a = Point2::new(-1.2, -3.0);
        let b = Point2::new(4.7, 2.3);
        let c = Point2::new(-4.2, 5.0);
        let mut rng = rand::rng();

        for _ in 0..10000 {
            let u = rng.random_range(0.0..1.0);
            let v = rng.random_range(0.0..u);
            let w = 1.0 - u - v;

            let p = Point2::from(a.coords * u + b.coords * v + c.coords * w);

            let bary = barycentric(&a, &b, &c, &p);
            assert_relative_eq!(bary[0], u, epsilon = 1e-6);
            assert_relative_eq!(bary[1], v, epsilon = 1e-6);
            assert_relative_eq!(bary[2], w, epsilon = 1e-6);
        }
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

    #[test]
    fn single_point_between() {
        let p0 = Point2::new(0.0, 0.0);
        let p1 = Point2::new(2.0, 0.0);
        let points = evenly_spaced_points_between(&p0, &p1, 1);

        assert_eq!(points.len(), 1);
        assert_relative_eq!(points[0], Point2::new(1.0, 0.0));
    }

    #[test]
    fn test_fill_gaps() {
        let points = vec![Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)];
        let filled = fill_gaps(&points, 1.5);

        assert_eq!(filled.len(), 3);
        assert_relative_eq!(filled[0], Point2::new(0.0, 0.0));
        assert_relative_eq!(filled[1], Point2::new(1.0, 0.0));
        assert_relative_eq!(filled[2], Point2::new(2.0, 0.0));
    }

    #[test]
    fn max_point_in_dir() {
        let dir = Vector2::new(-10.0, -10.0);

        let points = vec![
            Point2::new(10.0, 0.0),
            Point2::new(11.0, 1.0),
            Point2::new(12.0, 2.0),
            Point2::new(13.0, 3.0),
            Point2::new(14.0, 4.0),
        ];

        let (_, max) = max_point_in_direction(&points, &dir).unwrap();
        assert_relative_eq!(max, Point2::new(10.0, 0.0));
    }
}
