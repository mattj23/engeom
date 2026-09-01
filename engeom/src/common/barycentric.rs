//! Barycentric coordinates on a triangle, in any dimension.
//!
//! Barycentric coordinates express a point as a weighted blend of a triangle's three vertices, with
//! the weights summing to one. They are the natural currency whenever something has to be carried
//! *across* a triangle rather than merely located on it: a position, a normal, a scalar field
//! sampled at the corners, or an error radius. Getting the weights back rather than the point is the
//! whole reason several of these functions exist in this form.
//!
//! Weights are always ordered to match the vertices they were computed against, so `w[0]` belongs to
//! `a`, `w[1]` to `b` and `w[2]` to `c`.
//!
//! # Which one to use
//!
//! - [`barycentric`] locates a point relative to a triangle. Weights may be negative, which is how
//!   it says the point lies outside.
//! - [`closest_barycentric`] locates the nearest point *on* the triangle. Weights are non-negative
//!   and sum to one, so the result is inside by construction however far outside the query is.
//! - [`barycentric_within2`] is the planar case with containment folded in, returning `None`
//!   instead of negative weights.
//! - [`barycentric_point`] goes the other way, turning weights back into a position.
//! - [`barycentric_grid`] covers a triangle with weights at a requested spacing.

use crate::Point3;
use crate::common::{PCoords, linear_space};
use parry3d_f64::na::{Point, SVector};

/// Calculate the barycentric coordinates of a point `p` with respect to the triangle defined by
/// points `a`, `b`, and `c` in D-dimensional space.  The barycentric coordinates are a way of
/// expressing the position of a point within a triangle as a combination of the triangle's
/// vertices.
///
/// A weight may be negative, which means `p` lies outside the triangle on the far side of the edge
/// opposite that vertex. Use [`closest_barycentric`] when the answer has to be a point on the
/// triangle, or [`barycentric_within2`] in the plane when "outside" should be an explicit `None`.
///
/// For a point off the plane of the triangle in 3D, the weights describe the projection of `p` onto
/// that plane; the out-of-plane distance is not recoverable from them.
///
/// A degenerate triangle has no barycentric frame, and this returns `[0.0, 0.0, 0.0]` for one rather
/// than dividing by zero. That is not a valid weight triple, since it does not sum to one, so a
/// caller which can be handed degenerate input should check.
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

/// The point at barycentric weights `w` on the triangle `a`, `b`, `c`.
///
/// The inverse of [`barycentric`] for weights summing to one, and the thing to reach for instead of
/// writing the blend out by hand.
///
/// # Arguments
///
/// * `a`: the first vertex of the triangle
/// * `b`: the second vertex of the triangle
/// * `c`: the third vertex of the triangle
/// * `w`: the barycentric weights, matching the vertex order
///
/// returns: Point<f64, D>
pub fn barycentric_point<const D: usize>(
    a: &impl PCoords<D>,
    b: &impl PCoords<D>,
    c: &impl PCoords<D>,
    w: [f64; 3],
) -> Point<f64, D> {
    Point::from(a.coords() * w[0] + b.coords() * w[1] + c.coords() * w[2])
}

/// The barycentric weights of the closest point on a triangle to `p`.
///
/// Ericson's region test rather than a projection followed by a clamp, because the answer is wanted
/// as weights and not as a position: the weights are what carry per-vertex quantities across the
/// triangle, and recovering them from a clamped position is both slower and less accurate. Use
/// [`barycentric_point`] if the position is what is actually needed.
///
/// The three weights are non-negative and sum to one, so the closest point is inside the triangle by
/// construction even when `p` is far outside it. That is the difference from [`barycentric`], which
/// reports where `p` is rather than where the triangle can reach.
///
/// # Arguments
///
/// * `a`: the first vertex of the triangle
/// * `b`: the second vertex of the triangle
/// * `c`: the third vertex of the triangle
/// * `p`: the point to find the closest point on the triangle to
///
/// returns: [f64; 3]
pub fn closest_barycentric<const D: usize>(
    a: &impl PCoords<D>,
    b: &impl PCoords<D>,
    c: &impl PCoords<D>,
    p: &impl PCoords<D>,
) -> [f64; 3] {
    let (a, b, c, x) = (a.coords(), b.coords(), c.coords(), p.coords());
    let (ab, ac) = (b - a, c - a);

    let ap = x - a;
    let d1 = ab.dot(&ap);
    let d2 = ac.dot(&ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return [1.0, 0.0, 0.0];
    }

    let bp = x - b;
    let d3 = ab.dot(&bp);
    let d4 = ac.dot(&bp);
    if d3 >= 0.0 && d4 <= d3 {
        return [0.0, 1.0, 0.0];
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return [1.0 - v, v, 0.0];
    }

    let cp = x - c;
    let d5 = ab.dot(&cp);
    let d6 = ac.dot(&cp);
    if d6 >= 0.0 && d5 <= d6 {
        return [0.0, 0.0, 1.0];
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return [1.0 - w, 0.0, w];
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return [0.0, 1.0 - w, w];
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    [1.0 - v - w, v, w]
}

/// Twice the signed area of a planar triangle, positive when `a`, `b`, `c` wind counterclockwise.
///
/// Twice rather than the area itself because that is the form the barycentric denominators want, and
/// halving it only to multiply it back is wasted work and a lost bit of precision. Use
/// [`triangle_area`](crate::common::triangle_area) when the actual area is what is wanted.
///
/// The sign is the useful part: it says which way a triangle is wound, so comparing signs across a
/// set of triangles detects a projection that folds.
pub fn signed_area2(a: &impl PCoords<2>, b: &impl PCoords<2>, c: &impl PCoords<2>) -> f64 {
    let (a, b, c) = (a.coords(), b.coords(), c.coords());
    (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
}

/// The barycentric weights of `p` in a planar triangle, or `None` if it lies outside.
///
/// The signed-area formulation rather than [`barycentric`]'s dot products, which is both cheaper and
/// better conditioned in the plane, plus a containment test folded in.
///
/// A small negative tolerance is allowed, so a point sitting exactly on an edge is claimed by one of
/// the two triangles meeting there rather than by neither. Weights are clamped to be non-negative on
/// the way out, so a returned triple is always usable as a blend.
///
/// `None` also covers a degenerate triangle, which has no barycentric frame to answer with.
///
/// # Arguments
///
/// * `a`: the first vertex of the triangle
/// * `b`: the second vertex of the triangle
/// * `c`: the third vertex of the triangle
/// * `p`: the point to locate
///
/// returns: Option<[f64; 3]>
pub fn barycentric_within2(
    a: &impl PCoords<2>,
    b: &impl PCoords<2>,
    c: &impl PCoords<2>,
    p: &impl PCoords<2>,
) -> Option<[f64; 3]> {
    let total = signed_area2(a, b, c);
    if total.abs() < 1.0e-14 {
        return None;
    }

    let (a, b, c, q) = (a.coords(), b.coords(), c.coords(), p.coords());
    let sub = |u: &SVector<f64, 2>, v: &SVector<f64, 2>| {
        (v[0] - u[0]) * (q[1] - u[1]) - (q[0] - u[0]) * (v[1] - u[1])
    };

    let w0 = sub(&b, &c) / total;
    let w1 = sub(&c, &a) / total;
    let w2 = sub(&a, &b) / total;

    let eps = -1.0e-9;
    if w0 < eps || w1 < eps || w2 < eps {
        return None;
    }

    Some([w0.max(0.0), w1.max(0.0), w2.max(0.0)])
}

/// Cover a triangle with barycentric weights spaced no further apart than `max_spacing`.
///
/// The grid is laid out along whichever median is longest, so a sliver is sampled across its short
/// direction rather than being given a grid sized for its long one. Points are offset half a step in
/// from the edges, which keeps the samples of two triangles sharing an edge from landing on top of
/// each other.
///
/// # Arguments
///
/// * `a`: the first vertex of the triangle
/// * `b`: the second vertex of the triangle
/// * `c`: the third vertex of the triangle
/// * `max_spacing`: the largest gap allowed between neighbouring samples
///
/// returns: Vec<[f64; 3]>
pub fn barycentric_grid(a: &Point3, b: &Point3, c: &Point3, max_spacing: f64) -> Vec<[f64; 3]> {
    let mut result = Vec::new();
    let va = a - barycentric_point(a, b, c, [0.0, 0.5, 0.5]);
    let vb = b - barycentric_point(a, b, c, [0.5, 0.0, 0.5]);
    let vc = c - barycentric_point(a, b, c, [0.5, 0.5, 0.0]);

    let na = (va.norm() / max_spacing).ceil() as usize + 3;
    let nb = (vb.norm() / max_spacing).ceil() as usize + 3;
    let nc = (vc.norm() / max_spacing).ceil() as usize + 3;

    if na >= nb && na >= nc {
        let op_edge = (b - c).norm();
        for (bca, bcb, bcc) in bc_order(na, op_edge, max_spacing) {
            result.push([bca, bcb, bcc]);
        }
    } else if nb >= na && nb >= nc {
        let op_edge = (a - c).norm();
        for (bcb, bcc, bca) in bc_order(nb, op_edge, max_spacing) {
            result.push([bca, bcb, bcc]);
        }
    } else {
        let op_edge = (a - b).norm();
        for (bcc, bcb, bca) in bc_order(nc, op_edge, max_spacing) {
            result.push([bca, bcb, bcc]);
        }
    }

    result
}

fn bc_order(n0: usize, op_edge: f64, max_spacing: f64) -> Vec<(f64, f64, f64)> {
    let mut result = Vec::new();
    let spacing = 1.0 / n0 as f64;
    for bc0 in linear_space(spacing * 0.5, 1.0 - spacing * 0.5, n0).iter() {
        let leftover = 1.0 - bc0;
        let width = (1.0 - bc0) * op_edge;
        let nw = (width / max_spacing).ceil() as usize + 3;
        let sw = 1.0 / nw as f64;
        for bc1 in linear_space(sw * 0.5, 1.0 - sw * 0.5, nw).iter() {
            let bc2 = (1.0 - bc1) * leftover;
            let bc1 = bc1 * leftover;
            result.push((*bc0, bc1, bc2));
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::kd_tree::{KdTree, KdTreeSearch};
    use crate::common::points::evenly_spaced_points_between;
    use crate::{Point2, Point3};
    use approx::assert_relative_eq;
    use rand::RngExt;

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

    #[test]
    fn barycentric_point_inverts_barycentric() {
        let (a, b, c) = (
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 1.0),
            Point3::new(0.5, 3.0, -1.0),
        );
        let p = Point3::from(a.coords * 0.2 + b.coords * 0.5 + c.coords * 0.3);

        let w = barycentric(&a, &b, &c, &p);
        assert_relative_eq!(barycentric_point(&a, &b, &c, w), p, epsilon = 1e-12);
    }

    #[test]
    fn barycentric_is_negative_outside_the_triangle() {
        let (a, b, c) = (
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        );

        // Beyond the edge opposite `a`, so `a`'s weight is the one that goes negative.
        let w = barycentric(&a, &b, &c, &Point2::new(1.0, 1.0));
        assert!(w[0] < 0.0, "expected a negative first weight, got {w:?}");
        assert_relative_eq!(w[0] + w[1] + w[2], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn barycentric_of_a_degenerate_triangle_is_not_a_blend() {
        let (a, b, c) = (
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(2.0, 2.0),
        );
        assert_eq!(barycentric(&a, &b, &c, &Point2::new(0.5, 0.5)), [0.0; 3]);
    }

    #[test]
    fn closest_barycentric_stays_inside_the_triangle() {
        let (a, b, c) = (
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        );

        // Every region of Ericson's test: the three vertices, the three edges, the interior, and a
        // point off the plane.
        for q in [
            Point3::new(-5.0, -5.0, 0.0),
            Point3::new(9.0, -3.0, 0.0),
            Point3::new(-3.0, 9.0, 0.0),
            Point3::new(0.5, -4.0, 0.0),
            Point3::new(-4.0, 0.5, 0.0),
            Point3::new(4.0, 4.0, 0.0),
            Point3::new(0.25, 0.25, 0.0),
            Point3::new(0.25, 0.25, 7.0),
        ] {
            let w = closest_barycentric(&a, &b, &c, &q);
            assert!(w.iter().all(|x| *x >= 0.0), "{q:?} gave {w:?}");
            assert_relative_eq!(w[0] + w[1] + w[2], 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn stress_closest_barycentric_beats_dense_sampling() {
        // The independent check: no point sampled across the triangle may be closer to the query
        // than the analytic answer claims.
        let (a, b, c) = (
            Point3::new(-1.0, 0.0, 0.3),
            Point3::new(2.0, 0.5, -0.4),
            Point3::new(0.2, 2.5, 1.1),
        );
        let mut rng = rand::rng();

        for _ in 0..200 {
            let q = Point3::new(
                rng.random_range(-4.0..4.0),
                rng.random_range(-4.0..4.0),
                rng.random_range(-4.0..4.0),
            );
            let best =
                (barycentric_point(&a, &b, &c, closest_barycentric(&a, &b, &c, &q)) - q).norm();

            let n = 60;
            for i in 0..=n {
                for j in 0..=(n - i) {
                    let (u, v) = (i as f64 / n as f64, j as f64 / n as f64);
                    let p = barycentric_point(&a, &b, &c, [1.0 - u - v, u, v]);
                    assert!(
                        (p - q).norm() >= best - 1.0e-9,
                        "a sampled point beat the analytic closest point"
                    );
                }
            }
        }
    }

    #[test]
    fn barycentric_within2_rejects_points_outside() {
        let (a, b, c) = (
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        );

        assert!(barycentric_within2(&a, &b, &c, &Point2::new(0.25, 0.25)).is_some());
        assert!(barycentric_within2(&a, &b, &c, &Point2::new(1.0, 1.0)).is_none());
        assert!(barycentric_within2(&a, &b, &c, &Point2::new(-0.5, 0.5)).is_none());

        // A point exactly on an edge is claimed rather than dropped.
        let w = barycentric_within2(&a, &b, &c, &Point2::new(0.5, 0.0)).unwrap();
        assert_relative_eq!(w[2], 0.0, epsilon = 1e-9);

        // A degenerate triangle has no frame to answer with.
        let d = Point2::new(2.0, 0.0);
        assert!(barycentric_within2(&a, &b, &d, &Point2::new(0.5, 0.0)).is_none());
    }

    #[test]
    fn barycentric_within2_agrees_with_barycentric_inside() {
        let (a, b, c) = (
            Point2::new(-1.0, -0.5),
            Point2::new(2.5, 0.25),
            Point2::new(0.0, 3.0),
        );
        let p = Point2::from(a.coords * 0.3 + b.coords * 0.45 + c.coords * 0.25);

        let want = barycentric(&a, &b, &c, &p);
        let got = barycentric_within2(&a, &b, &c, &p).unwrap();
        for k in 0..3 {
            assert_relative_eq!(got[k], want[k], epsilon = 1e-12);
        }
    }

    #[test]
    fn signed_area2_changes_sign_with_winding() {
        let (a, b, c) = (
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(0.0, 3.0),
        );
        assert_relative_eq!(signed_area2(&a, &b, &c), 6.0, epsilon = 1e-12);
        assert_relative_eq!(signed_area2(&a, &c, &b), -6.0, epsilon = 1e-12);
    }

    #[test]
    fn barycentric_grid_spacing() {
        // The following conditions should be true:
        // 1. No point on the edge of the triangle is more than max_spacing/2 from the nearest grid
        //    point.
        // 2. No point in the grid is more than max_spacing from another point in the grid.
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);

        let max_spacing = 0.1;
        let grid = barycentric_grid(&a, &b, &c, max_spacing);

        // Check for NAN
        for bc in &grid {
            assert!(
                !bc.iter().any(|&x| x.is_nan()),
                "Barycentric coordinate contains NaN: {:?}",
                bc
            );
        }

        let grid_points = grid
            .iter()
            .map(|bc| barycentric_point(&a, &b, &c, *bc))
            .collect::<Vec<_>>();

        // Check for NAN
        for point in &grid_points {
            assert!(
                !point.coords.iter().any(|&x| x.is_nan()),
                "Point contains NaN: {:?}",
                point
            );
        }

        let tree = KdTree::try_new(&grid_points).expect("Tree construction failed");

        // Check that no point in the grid is more than max_spacing from another point in the grid
        for point in &grid_points {
            let neighbors = tree.within(point, max_spacing);
            assert!(neighbors.len() > 1, "Point {:?} has no neighbors", point);
        }

        // Check that no point on the edge of the triangle is more than max_spacing/2 from the
        // nearest grid point
        let mut edge_points = evenly_spaced_points_between(&a, &b, 100);
        edge_points.extend(evenly_spaced_points_between(&b, &c, 100));
        edge_points.extend(evenly_spaced_points_between(&c, &a, 100));

        for edge_point in edge_points {
            let (_, d) = tree.nearest_one(&edge_point);
            assert!(
                d <= max_spacing * 0.7,
                "Edge point {:?} is too far from nearest grid point: {}",
                edge_point,
                d
            );
        }
    }
}
