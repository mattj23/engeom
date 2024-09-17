use super::signed_angle;
use crate::common::points::mid_point;
use parry2d_f64::bounding_volume::SimdAabb;
use parry2d_f64::math::{SimdBool, SimdReal, DIM, SIMD_WIDTH};
use parry2d_f64::na::{Isometry2, Point2, SimdPartialOrd, SimdValue, Vector2};
use parry2d_f64::partitioning::{SimdVisitStatus, SimdVisitor};
use parry2d_f64::query::{Ray, SimdRay};
use parry2d_f64::shape::{Polyline, SimdCompositeShape};

use crate::geom2::line2::{intersect_rays, Line2};
use serde::{Deserialize, Serialize};

/// A `SpanningRay` is a special case of ray which spans two points in a polyline, typically when
/// there is a closed polyline and a ray that crosses from one side to the other.  It is a wrapper
/// around a ray where both the ray origin and ray origin + ray direction are points on the
/// polyline.
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct SpanningRay {
    ray: Ray,
}

impl SpanningRay {
    pub fn ray(&self) -> Ray {
        self.ray
    }

    pub fn new(p0: Point2<f64>, p1: Point2<f64>) -> Self {
        Self {
            ray: Ray::new(p0, p1 - p0),
        }
    }

    /// Computes and returns a symmetry ray between this and another spanning ray
    pub fn symmetry(&self, other: &SpanningRay) -> Ray {
        let angle = signed_angle(&self.ray.dir, &other.ray.dir) * 0.5;
        Ray::new(
            mid_point(&self.ray.origin, &other.ray.origin),
            Isometry2::rotation(angle) * self.ray.dir,
        )
    }

    pub fn reversed(&self) -> Self {
        Self::new(self.ray.point_at(1.0), self.ray.origin.clone())
    }
}

pub fn ray_intersect_with_edge(line: &Polyline, ray: &Ray, edge_index: usize) -> Option<f64> {
    let v0 = line.vertices()[edge_index].clone();
    let v1 = line.vertices()[edge_index + 1].clone();
    let dir = v1 - &v0;
    let edge_ray = Ray::new(v0, dir);
    if let Some((t0, t1)) = intersect_rays(ray, &edge_ray) {
        if (0.0..=1.0).contains(&t1) {
            Some(t0)
        } else {
            None
        }
    } else {
        None
    }
}

impl Line2 for SpanningRay {
    fn origin(&self) -> Point2<f64> {
        self.ray.origin.clone()
    }

    fn dir(&self) -> Vector2<f64> {
        self.ray.dir
    }

    fn at(&self, t: f64) -> Point2<f64> {
        self.ray.point_at(t)
    }
}

pub fn max_intersection(line: &Polyline, ray: &Ray) -> Option<f64> {
    let ts: Vec<f64> = polyline_intersections(line, ray)
        .iter()
        .map(|(t, _)| *t)
        .collect();
    ts.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).cloned()
}

/// Finds the projected distance of the farthest point in the polyline from a ray origin in the
/// ray direction
pub fn farthest_point_direction_distance(line: &Polyline, ray: &Ray) -> f64 {
    let mut farthest = f64::MIN;
    let n = ray.dir.normalize();
    for v in line.vertices().iter() {
        farthest = farthest.max(n.dot(&(v - ray.origin.clone())));
    }

    farthest
}

/// Attempts to create a "spanning ray" through the polyline along the parameterized line
/// represented by the ray argument. A "spanning ray" is a ray that starts on the surface of
/// the polyline and passes through it ending at the other side, such that t=0 is an
/// intersection with the polyline on one side, t=1.0 is an intersection on the other side, and
/// there are no additional intersections between them. The spanning ray will have the same
/// direction as the original intersection ray.
pub fn spanning_ray(line: &Polyline, ray: &Ray) -> Option<SpanningRay> {
    let mut results = polyline_intersections(line, ray);
    results.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    if results.len() == 2 {
        Some(SpanningRay::new(
            ray.point_at(results[0].0),
            ray.point_at(results[1].0),
        ))
    } else {
        None
    }
}

pub fn polyline_intersections(polyline: &Polyline, ray: &Ray) -> Vec<(f64, usize)> {
    let mut results = Vec::new();
    let mut visitor = RayVisitor::new(ray);
    polyline.qbvh().traverse_depth_first(&mut visitor);

    for i in visitor.collector.iter() {
        if let Some(t) = ray_intersect_with_edge(polyline, ray, *i as usize) {
            results.push((t, *i as usize));
        }
    }

    results.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    results.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-8);

    results
}

struct RayVisitor {
    ray: SimdRay,
    collector: Vec<u32>,
}

impl RayVisitor {
    pub fn new(ray: &Ray) -> Self {
        Self {
            ray: SimdRay::splat(*ray),
            collector: Vec::new(),
        }
    }
}

impl SimdVisitor<u32, SimdAabb> for RayVisitor {
    fn visit(
        &mut self,
        bv: &SimdAabb,
        data: Option<[Option<&u32>; SIMD_WIDTH]>,
    ) -> SimdVisitStatus {
        let mask = cast_ray(bv, &self.ray).0;
        if let Some(data) = data {
            for (i, &d_opt) in data.iter().enumerate().take(SIMD_WIDTH) {
                if mask.extract(i) {
                    if let Some(d) = d_opt {
                        self.collector.push(*d);
                    }
                }
            }
        }
        SimdVisitStatus::MaybeContinue(mask)
    }
}

/// Performs a ray cast except that it does not fail with negative t values.
fn cast_ray(bv: &SimdAabb, ray: &SimdRay) -> (SimdBool, SimdReal) {
    let zero = SimdReal::splat(0.0);
    let one = SimdReal::splat(1.0);
    let infinity = SimdReal::splat(f64::MAX);

    let mut hit = SimdBool::splat(true);
    let mut tmin = SimdReal::splat(f64::MIN);
    let mut tmax = SimdReal::splat(f64::MAX);

    // TODO: could this be optimized more considering we really just need a boolean answer?
    for i in 0usize..DIM {
        let is_not_zero = ray.dir[i].simd_ne(zero);
        let is_zero_test = ray.origin[i].simd_ge(bv.mins[i]) & ray.origin[i].simd_le(bv.maxs[i]);
        let is_not_zero_test = {
            let denom = one / ray.dir[i];
            let mut inter_with_near_plane =
                ((bv.mins[i] - ray.origin[i]) * denom).select(is_not_zero, -infinity);
            let mut inter_with_far_plane =
                ((bv.maxs[i] - ray.origin[i]) * denom).select(is_not_zero, infinity);

            let gt = inter_with_near_plane.simd_gt(inter_with_far_plane);
            simd_swap(gt, &mut inter_with_near_plane, &mut inter_with_far_plane);

            tmin = tmin.simd_max(inter_with_near_plane);
            tmax = tmax.simd_min(inter_with_far_plane);

            tmin.simd_le(tmax)
        };

        hit = hit & is_not_zero_test.select(is_not_zero, is_zero_test);
    }

    (hit, tmin)
}

fn simd_swap(do_swap: SimdBool, a: &mut SimdReal, b: &mut SimdReal) {
    let _a = *a;
    *a = b.select(do_swap, *a);
    *b = _a.select(do_swap, *b);
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use test_case::test_case;

    #[test]
    fn test_symmetry_ray() {
        let s0 = SpanningRay::new(Point2::new(1.0, 0.0), Point2::new(2.0, 0.0));
        let s1 = SpanningRay::new(Point2::new(0.0, 1.0), Point2::new(0.0, 2.0));
        let r = s0.symmetry(&s1);
        let en = Vector2::new(1.0, 1.0).normalize();

        assert_relative_eq!(0.5, r.origin.x, epsilon = 1e-8);
        assert_relative_eq!(0.5, r.origin.y, epsilon = 1e-8);

        assert_relative_eq!(en.x, r.dir.x, epsilon = 1e-8);
        assert_relative_eq!(en.y, r.dir.y, epsilon = 1e-8);
    }

    #[test_case((3.1, 4.7, 0.9, 3.5), 1.297781)]
    #[test_case((-3.6, -0.6, 2.4, -4.5), 7.517647)]
    #[test_case((1.5, -1.9, -2.6, 0.1), 9.285442)]
    #[test_case((-0.8, 4.5, 3.5, 4.0), 3.076157)]
    #[test_case((3.7, 3.9, 0.4, 3.4), 2.249194)]
    #[test_case((1.4, 4.7, 4.9, -2.0), 6.112484)]
    #[test_case((-0.3, -2.3, 2.3, -0.9), 5.061102)]
    #[test_case((3.6, 2.9, -3.3, -3.3), 12.657211)]
    #[test_case((-1.8, -2.1, -2.7, 0.7), 6.134232)]
    #[test_case((-3.9, -4.5, 1.9, 1.7), 11.441415)]
    fn test_farthest_dist_direction(a: (f64, f64, f64, f64), d: f64) {
        let r = Ray::new(Point2::new(a.0, a.1), Vector2::new(a.2, a.3));
        let result = farthest_point_direction_distance(&sample_polyline(), &r);
        assert_relative_eq!(d, result, epsilon = 1e-5);
    }

    fn naive_ray_intersections(line: &Polyline, ray: &Ray) -> Vec<f64> {
        let mut results = Vec::new();
        for i in 0..line.vertices().len() - 1 {
            if let Some(point) = ray_intersect_with_edge(line, ray, i) {
                results.push(point);
            }
        }

        results
    }

    #[test]
    fn test_intersections_against_naive() {
        use std::f64::consts::PI;

        let line = sample_polyline();

        for i in 1..360 {
            let ai = Isometry2::rotation(i as f64 / 180.0 * PI) * Point2::new(10.0, 0.0);
            for j in 1..360 {
                let aj = Isometry2::rotation(j as f64 / 180.0 * PI) * Vector2::new(1.0, 0.0);
                let ray = Ray::new(ai, aj);

                let mut naive = naive_ray_intersections(&line, &ray);
                let mut fast: Vec<f64> = polyline_intersections(&line, &ray)
                    .iter()
                    .map(|(t, _)| *t)
                    .collect();
                naive.sort_by(|a, b| a.partial_cmp(b).unwrap());
                naive.dedup_by(|a, b| (*a - *b).abs() < 1e-5);
                fast.sort_by(|a, b| a.partial_cmp(b).unwrap());

                assert_eq!(naive, fast);
            }
        }
    }

    fn sample_polyline() -> Polyline {
        let points: Vec<Point2<f64>> = vec![
            Point2::new(5.0, 0.0),
            Point2::new(3.5, 0.9),
            Point2::new(4.0, 2.3),
            Point2::new(3.3, 3.3),
            Point2::new(2.7, 4.7),
            Point2::new(1.7, 6.4),
            Point2::new(0.0, 5.9),
            Point2::new(-1.5, 5.7),
            Point2::new(-3.7, 6.4),
            Point2::new(-5.3, 5.3),
            Point2::new(-6.4, 3.7),
            Point2::new(-7.1, 1.9),
            Point2::new(-7.3, 0.0),
            Point2::new(-7.8, -2.1),
            Point2::new(-6.3, -3.7),
            Point2::new(-5.7, -5.7),
            Point2::new(-3.7, -6.3),
            Point2::new(-1.7, -6.2),
            Point2::new(-0.0, -7.2),
            Point2::new(1.5, -5.6),
            Point2::new(2.4, -4.2),
            Point2::new(3.9, -3.9),
            Point2::new(4.9, -2.9),
            Point2::new(4.9, -1.3),
            Point2::new(5.0, 0.0),
        ];
        Polyline::new(points, None)
    }
}
