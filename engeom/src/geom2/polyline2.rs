use super::{Aabb2, Iso2, Point2, Vector2, signed_angle};
use crate::common::PCoords;
use crate::common::points::mid_point;
use crate::geom2::line2::{Line2, intersect_rays};
use parry2d_f64::partitioning::TraversalAction;
// use parry2d_f64::na::{, SimdPartialOrd, SimdValue};
use parry2d_f64::query::{Ray, RayCast};
use parry2d_f64::shape::Polyline;
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

    pub fn new(p0: &impl PCoords<2>, p1: &impl PCoords<2>) -> Self {
        let p0 = Point2::from(p0.coords());
        let p1 = Point2::from(p1.coords());
        Self {
            ray: Ray::new(p0, p1 - p0),
        }
    }

    /// Computes and returns a symmetry ray between this and another spanning ray
    pub fn symmetry(&self, other: &SpanningRay) -> Ray {
        let angle = signed_angle(&self.ray.dir, &other.ray.dir) * 0.5;
        Ray::new(
            mid_point(&self.ray.origin, &other.ray.origin),
            Iso2::rotation(angle) * self.ray.dir,
        )
    }

    pub fn reversed(&self) -> Self {
        Self::new(&self.ray.point_at(1.0), &self.ray.origin)
    }
}

pub fn ray_intersect_with_edge(line: &Polyline, ray: &Ray, edge_index: usize) -> Option<f64> {
    let v0 = line.vertices()[edge_index];
    let v1 = line.vertices()[edge_index + 1];
    let dir = v1 - v0;
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
    fn origin(&self) -> Point2 {
        self.ray.origin
    }

    fn dir(&self) -> Vector2 {
        self.ray.dir
    }

    fn at(&self, t: f64) -> Point2 {
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
        farthest = farthest.max(n.dot(&(v - ray.origin)));
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
            &ray.point_at(results[0].0),
            &ray.point_at(results[1].0),
        ))
    } else {
        None
    }
}

pub fn polyline_intersections(polyline: &Polyline, ray: &Ray) -> Vec<(f64, usize)> {
    let mut candidates = Vec::new();
    let r_inv = Vector2::new(1.0 / ray.dir.x, 1.0 / ray.dir.y);

    polyline.bvh().traverse(|node| {
        if let Some(data) = node.leaf_data() {
            candidates.push(data)
        };

        if !slab_method(&node.aabb(), ray, &r_inv) {
            TraversalAction::Prune
        } else {
            TraversalAction::Continue
        }
    });

    let mut results = Vec::new();
    for i in candidates.iter() {
        if let Some(t) = ray_intersect_with_edge(polyline, ray, *i as usize) {
            results.push((t, *i as usize));
        }
    }

    results.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    results.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-8);

    results
}

fn slab_method(bv: &Aabb2, ray: &Ray, n_inv: &Vector2) -> bool {
    let mut t1 = (bv.mins.x - ray.origin.x) * n_inv.x;
    let mut t2 = (bv.maxs.x - ray.origin.x) * n_inv.x;

    let tmin = t1.min(t2);
    let tmax = t1.max(t2);

    t1 = (bv.mins.y - ray.origin.y) * n_inv.y;
    t2 = (bv.maxs.y - ray.origin.y) * n_inv.y;

    let tmin = tmin.max(t1.min(t2).min(tmax));
    let tmax = tmax.min(t1.max(t2).max(tmin));

    tmax >= tmin
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use test_case::test_case;

    #[test]
    fn test_symmetry_ray() {
        let s0 = SpanningRay::new(&Point2::new(1.0, 0.0), &Point2::new(2.0, 0.0));
        let s1 = SpanningRay::new(&Point2::new(0.0, 1.0), &Point2::new(0.0, 2.0));
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
            let ai = Iso2::rotation(i as f64 / 180.0 * PI) * Point2::new(10.0, 0.0);
            for j in 1..360 {
                let aj = Iso2::rotation(j as f64 / 180.0 * PI) * Vector2::new(1.0, 0.0);
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
        let points: Vec<Point2> = vec![
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
