//! Dimension-generic minimum enclosing balls using Welzl's algorithm. This module implements
//! `Circle2::from_min_enclosing` and `Sphere3::from_min_enclosing`; their public entry points live
//! in `geom2` and `geom3`.
//!
//! The move-to-front formulation has expected O(n) time at fixed dimension and recursion depth
//! bounded by `D + 2`. Support balls are found from a Gram system. A final pass adjusts the radius
//! for floating-point error so it contains every input point.

use crate::Result;
use crate::common::PCoords;
use crate::na::{Point, SVector};
use rand::SeedableRng;
use rand::prelude::{SliceRandom, StdRng};

/// Relative squared-radius tolerance for candidate containment tests. It affects candidate
/// updates only; the final pass guarantees containment.
const REL_TOL: f64 = 1e-10;

/// Relative tolerance for rejecting circumballs that are not equidistant from their support.
const EQ_TOL: f64 = 1e-9;

/// Relative pivot threshold for treating a support's Gram system as singular.
const PIVOT_TOL: f64 = 1e-12;

/// Fixed seed for a deterministic shuffle that avoids adversarial input ordering.
const SHUFFLE_SEED: u64 = 0x5EED_BA11;

#[derive(Copy, Clone, Debug)]
struct Ball<const D: usize> {
    center: SVector<f64, D>,
    r2: f64,
}

impl<const D: usize> Ball<D> {
    fn contains(&self, p: &SVector<f64, D>) -> bool {
        (p - self.center).norm_squared() <= self.r2 * (1.0 + REL_TOL)
    }
}

fn contains<const D: usize>(ball: &Option<Ball<D>>, p: &SVector<f64, D>) -> bool {
    ball.as_ref().is_some_and(|b| b.contains(p))
}

/// Computes the smallest ball through `support`, centered in its affine hull.
///
/// Returns `None` for empty or nearly affinely dependent support.
fn circumball<const D: usize>(support: &[SVector<f64, D>]) -> Option<Ball<D>> {
    match support.len() {
        0 => None,
        1 => Some(Ball {
            center: support[0],
            r2: 0.0,
        }),
        k => {
            // Work relative to the first support point: with v_j = p_j - p0 and the center at
            // p0 + sum(lambda_j * v_j), equidistance from p0 and each p_j gives the Gram system
            // G * lambda = b with G_ij = v_i . v_j and b_i = (v_i . v_i) / 2. The m x m system
            // (m = k - 1 <= D) is small enough that a partial-pivot Gaussian elimination on
            // stack arrays beats a general decomposition and avoids the allocator/trait bounds
            // that nalgebra's factorizations require under const generics.
            let vs: Vec<SVector<f64, D>> = support[1..].iter().map(|p| p - support[0]).collect();
            let m = k - 1;
            let mut a = [[0.0_f64; D]; D];
            let mut b = [0.0_f64; D];
            for i in 0..m {
                for j in 0..m {
                    a[i][j] = vs[i].dot(&vs[j]);
                }
                b[i] = 0.5 * vs[i].norm_squared();
            }

            let mut scale = 0.0_f64;
            for row in a.iter().take(m) {
                for value in row.iter().take(m) {
                    scale = scale.max(value.abs());
                }
            }

            for col in 0..m {
                let mut piv = col;
                for row in (col + 1)..m {
                    if a[row][col].abs() > a[piv][col].abs() {
                        piv = row;
                    }
                }
                if a[piv][col].abs() <= scale * PIVOT_TOL {
                    return None;
                }
                a.swap(col, piv);
                b.swap(col, piv);
                let pivot_row = a[col];
                let pivot_b = b[col];
                for row in (col + 1)..m {
                    let f = a[row][col] / pivot_row[col];
                    for (x, p) in a[row][col..m].iter_mut().zip(&pivot_row[col..m]) {
                        *x -= f * p;
                    }
                    b[row] -= f * pivot_b;
                }
            }

            let mut lambda = [0.0_f64; D];
            for col in (0..m).rev() {
                let mut s = b[col];
                for j in (col + 1)..m {
                    s -= a[col][j] * lambda[j];
                }
                lambda[col] = s / a[col][col];
            }

            let mut offset = SVector::<f64, D>::zeros();
            for i in 0..m {
                offset += lambda[i] * vs[i];
            }
            let r2 = offset.norm_squared();

            // An ill-conditioned Gram solve can return a center which is not actually equidistant
            // from the support points; validate before accepting it.
            let tol = r2 * EQ_TOL;
            for v in &vs {
                if ((v - offset).norm_squared() - r2).abs() > tol {
                    return None;
                }
            }

            Some(Ball {
                center: support[0] + offset,
                r2,
            })
        }
    }
}

/// Computes the minimum ball enclosing `pts` and passing through `support`.
///
/// Each recursion level adds a support point, bounding depth by `D + 2`. Moving outside points to
/// the front limits rescanning in enclosing calls.
fn welzl<const D: usize>(
    pts: &mut [SVector<f64, D>],
    support: &mut Vec<SVector<f64, D>>,
) -> Option<Ball<D>> {
    let mut ball = circumball(support);
    if support.len() == D + 1 {
        return ball;
    }

    for i in 0..pts.len() {
        if !contains(&ball, &pts[i]) {
            support.push(pts[i]);
            let sub = welzl(&mut pts[..i], support);
            support.pop();
            if let Some(b) = sub {
                ball = Some(b);
                pts[..=i].rotate_right(1);
            }
            // A degenerate (affinely dependent) support yields no ball; keep the previous
            // candidate and skip the point, relying on the final inflation pass for containment.
        }
    }

    ball
}

/// Computes the minimum enclosing ball using Welzl's algorithm in expected O(n) time at a fixed
/// dimension.
///
/// The result is deterministic and contains every input point. A single point or identical points
/// produce a zero-radius ball; collinear points produce the diametral ball of the extreme pair.
///
/// # Arguments
///
/// * `points`: a slice of coordinates to enclose; must not be empty
///
/// returns: Result<(Point<f64, { D }>, f64), Box<dyn Error, Global>>
pub fn compute_min_ball<const D: usize>(
    points: &[impl PCoords<D>],
) -> Result<(Point<f64, D>, f64)> {
    if points.is_empty() {
        return Err("Cannot compute the minimum enclosing ball of an empty point set".into());
    }

    let mut pts: Vec<SVector<f64, D>> = points.iter().map(|p| p.coords()).collect();
    let mut rng = StdRng::seed_from_u64(SHUFFLE_SEED);
    pts.shuffle(&mut rng);

    let ball = welzl(&mut pts, &mut Vec::with_capacity(D + 1))
        .ok_or("Failed to compute the minimum enclosing ball")?;

    // Inflate the radius to the farthest point so that containment holds without a tolerance,
    // regardless of floating point jitter in the candidate updates.
    let mut r2 = ball.r2;
    for p in &pts {
        r2 = r2.max((p - ball.center).norm_squared());
    }

    Ok((Point::from(ball.center), r2.sqrt()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::dist;
    use crate::common::random_geometry::{RandomGeometry2, RandomGeometry3};
    use crate::{Point2, Point3};
    use approx::assert_relative_eq;

    /// Brute-force minimum enclosing ball: the smallest circumball of any support subset of size
    /// 1..=D+1 which contains every point. Valid only for small n.
    fn brute_force<const D: usize>(points: &[SVector<f64, D>]) -> (SVector<f64, D>, f64) {
        let n = points.len();
        let mut best: Option<Ball<D>> = None;

        let mut subset = Vec::new();

        // Enumerate subsets of size 1..=D+1 by index combinations
        fn combos(n: usize, k: usize) -> Vec<Vec<usize>> {
            if k == 0 {
                return vec![vec![]];
            }
            let mut out = Vec::new();
            for first in 0..n {
                for mut rest in combos(n - first - 1, k - 1) {
                    for r in rest.iter_mut() {
                        *r += first + 1;
                    }
                    let mut c = vec![first];
                    c.extend(rest);
                    out.push(c);
                }
            }
            out
        }

        for k in 1..=(D + 1).min(n) {
            for combo in combos(n, k) {
                subset.clear();
                subset.extend(combo.iter().map(|&i| points[i]));
                if let Some(ball) = circumball(&subset) {
                    let r_tol = ball.r2 * (1.0 + 1e-9) + 1e-12;
                    if points
                        .iter()
                        .all(|p| (p - ball.center).norm_squared() <= r_tol)
                        && best.map(|b| ball.r2 < b.r2).unwrap_or(true)
                    {
                        best = Some(ball);
                    }
                }
            }
        }

        let ball = best.expect("brute force found no enclosing ball");
        (ball.center, ball.r2.sqrt())
    }

    #[test]
    fn two_points_diametral_2d() {
        let points = [Point2::new(-1.0, 2.0), Point2::new(3.0, 2.0)];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point2::new(1.0, 2.0), epsilon = 1e-12);
        assert_relative_eq!(r, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn two_points_diametral_3d() {
        let points = [Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 2.0, 2.0)];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point3::new(1.0, 1.0, 1.0), epsilon = 1e-12);
        assert_relative_eq!(r, 3.0_f64.sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn equilateral_triangle_circumcircle() {
        let s3 = 3.0_f64.sqrt();
        let points = [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.5, s3 / 2.0),
        ];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point2::new(0.5, s3 / 6.0), epsilon = 1e-12);
        assert_relative_eq!(r, 1.0 / s3, epsilon = 1e-12);
    }

    #[test]
    fn regular_tetrahedron_circumsphere() {
        let points = [
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(-1.0, 1.0, -1.0),
            Point3::new(-1.0, -1.0, 1.0),
        ];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point3::origin(), epsilon = 1e-12);
        assert_relative_eq!(r, 3.0_f64.sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn obtuse_triangle_uses_diametral_ball() {
        // The circumcircle of an obtuse triangle is larger than the diametral ball of its longest
        // edge; the minimum enclosing ball must be the latter, with only two support points.
        let points = [
            Point2::new(0.0, 0.0),
            Point2::new(4.0, 0.0),
            Point2::new(1.0, 1.0),
        ];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point2::new(2.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(r, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn recovers_known_circle_2d() {
        let mut rg = RandomGeometry2::from_seed(1);
        let center = Point2::new(-3.0, 7.5);
        let radius = 2.25;
        let mut points = Vec::new();
        for _ in 0..50 {
            let angle = rg.angle_pos_2pi();
            points.push(center + radius * SVector::<f64, 2>::new(angle.cos(), angle.sin()));
        }
        for _ in 0..200 {
            let angle = rg.angle_pos_2pi();
            let r = rg.f64(0.0, radius * 0.99);
            points.push(center + r * SVector::<f64, 2>::new(angle.cos(), angle.sin()));
        }
        let (c, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(c, center, epsilon = 1e-9);
        assert_relative_eq!(r, radius, epsilon = 1e-9);
    }

    #[test]
    fn recovers_known_sphere_3d() {
        let mut rg = RandomGeometry3::from_seed(2);
        let center = Point3::new(4.0, -1.0, 10.0);
        let radius = 5.0;
        let mut points = Vec::new();
        for _ in 0..100 {
            points.push(center + radius * rg.unit_vec().into_inner());
        }
        for _ in 0..300 {
            points.push(center + rg.f64(0.0, radius * 0.99) * rg.unit_vec().into_inner());
        }
        let (c, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(c, center, epsilon = 1e-9);
        assert_relative_eq!(r, radius, epsilon = 1e-9);
    }

    #[test]
    fn collinear_points_2d() {
        let points = [
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(1.5, 0.0),
            Point2::new(2.0, 0.0),
        ];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point2::new(1.5, 0.0), epsilon = 1e-12);
        assert_relative_eq!(r, 1.5, epsilon = 1e-12);
    }

    #[test]
    fn coplanar_points_3d() {
        // Square corners in the z = 0 plane; the minimum sphere is the circumscribed circle's
        // sphere, not anything larger.
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point3::new(1.0, 1.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(r, 2.0_f64.sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn collinear_points_3d() {
        let points = [
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(2.0, 2.0, 2.0),
            Point3::new(4.0, 4.0, 4.0),
            Point3::new(3.0, 3.0, 3.0),
        ];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point3::new(2.5, 2.5, 2.5), epsilon = 1e-12);
        assert_relative_eq!(r, 1.5 * 3.0_f64.sqrt(), epsilon = 1e-12);
    }

    #[test]
    fn identical_points_zero_radius() {
        let points = [Point2::new(3.0, 4.0); 10];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point2::new(3.0, 4.0), epsilon = 1e-12);
        assert_relative_eq!(r, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn single_point_zero_radius() {
        let points = [Point3::new(1.0, 2.0, 3.0)];
        let (center, r) = compute_min_ball(&points).unwrap();
        assert_relative_eq!(center, Point3::new(1.0, 2.0, 3.0), epsilon = 1e-12);
        assert_relative_eq!(r, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn empty_input_is_error() {
        let points: Vec<Point2> = Vec::new();
        assert!(compute_min_ball(&points).is_err());
        let points: Vec<Point3> = Vec::new();
        assert!(compute_min_ball(&points).is_err());
    }

    #[test]
    fn stress_containment_2d() {
        for seed in 0..5 {
            let mut rg = RandomGeometry2::from_seed(seed);
            let points: Vec<Point2> = (0..2000).map(|_| rg.point(100.0)).collect();
            let (center, r) = compute_min_ball(&points).unwrap();
            for p in &points {
                assert!(dist(p, &center) <= r * (1.0 + 1e-9));
            }
        }
    }

    #[test]
    fn stress_containment_3d() {
        for seed in 0..5 {
            let mut rg = RandomGeometry3::from_seed(seed);
            let points: Vec<Point3> = (0..2000).map(|_| rg.point(100.0)).collect();
            let (center, r) = compute_min_ball(&points).unwrap();
            for p in &points {
                assert!(dist(p, &center) <= r * (1.0 + 1e-9));
            }
        }
    }

    #[test]
    fn stress_minimality_small_sets_2d() {
        let mut rg = RandomGeometry2::from_seed(10);
        for i in 0..200 {
            let n = 2 + (i % 7);
            let points: Vec<SVector<f64, 2>> = (0..n).map(|_| rg.point(10.0).coords).collect();
            let (center, r) = compute_min_ball(&points).unwrap();
            let (bc, br) = brute_force(&points);
            assert_relative_eq!(center.coords, bc, epsilon = 1e-9 * (1.0 + br));
            assert_relative_eq!(r, br, epsilon = 1e-9 * (1.0 + br));
        }
    }

    #[test]
    fn stress_minimality_small_sets_3d() {
        let mut rg = RandomGeometry3::from_seed(11);
        for i in 0..200 {
            let n = 2 + (i % 7);
            let points: Vec<SVector<f64, 3>> = (0..n).map(|_| rg.point(10.0).coords).collect();
            let (center, r) = compute_min_ball(&points).unwrap();
            let (bc, br) = brute_force(&points);
            assert_relative_eq!(center.coords, bc, epsilon = 1e-9 * (1.0 + br));
            assert_relative_eq!(r, br, epsilon = 1e-9 * (1.0 + br));
        }
    }

    #[test]
    fn shuffle_invariance() {
        let mut rg = RandomGeometry3::from_seed(42);
        let points: Vec<Point3> = (0..500).map(|_| rg.point(50.0)).collect();
        let (center, r) = compute_min_ball(&points).unwrap();

        let mut reversed = points.clone();
        reversed.reverse();
        let (c2, r2) = compute_min_ball(&reversed).unwrap();
        assert_relative_eq!(center, c2, epsilon = 1e-9);
        assert_relative_eq!(r, r2, epsilon = 1e-9);

        let mut rotated = points.clone();
        rotated.rotate_left(137);
        let (c3, r3) = compute_min_ball(&rotated).unwrap();
        assert_relative_eq!(center, c3, epsilon = 1e-9);
        assert_relative_eq!(r, r3, epsilon = 1e-9);
    }
}
