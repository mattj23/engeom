use super::Mesh;
use crate::common::SurfacePointCollection;
use crate::common::indices::index_vec;
use crate::common::kd_tree::KdTreeSearch;
use crate::common::points::{dist, mean_point};
use crate::common::poisson_disk::sample_poisson_disk;
use crate::{Iso3, KdTree3, Point3, SurfacePoint3, SvdBasis3, To2D, TransformBy};
use parry2d_f64::transformation::convex_hull;
use rand::prelude::SliceRandom;
use std::f64::consts::PI;
use std::num::NonZero;

impl Mesh {
    pub fn sample_uniform(&self, n: usize) -> Vec<SurfacePoint3> {
        let mut cumulative_areas = Vec::new();
        let mut total_area = 0.0;
        for tri in self.shape.triangles() {
            total_area += tri.area();
            cumulative_areas.push(total_area);
        }

        let mut result = Vec::new();
        for _ in 0..n {
            let r = rand::random::<f64>() * total_area;
            let tri_id = cumulative_areas
                .binary_search_by(|a| a.partial_cmp(&r).unwrap())
                .unwrap_or_else(|i| i);
            let tri = self.shape.triangle(tri_id as u32);
            let r1 = rand::random::<f64>();
            let r2 = rand::random::<f64>();
            let a = 1.0 - r1.sqrt();
            let b = r1.sqrt() * (1.0 - r2);
            let c = r1.sqrt() * r2;
            let v = tri.a.coords * a + tri.b.coords * b + tri.c.coords * c;
            result.push(SurfacePoint3::new(Point3::from(v), tri.normal().unwrap()));
        }

        result
    }

    pub fn sample_poisson(&self, radius: f64) -> Vec<SurfacePoint3> {
        let starting = self.sample_dense(radius * 0.5);
        // TODO: this can be more efficient without all the copying
        let points = (&starting).clone_points();
        let mut rng = rand::rng();
        let mut indices = index_vec(None, starting.len());
        indices.shuffle(&mut rng);

        // Deal with kiddo bug
        let to_take = sample_poisson_disk(&points, &indices, radius);
        let to_take = sample_poisson_disk(&points, &to_take, radius);
        to_take.into_iter().map(|i| starting[i]).collect()
    }

    pub fn sample_dense(&self, max_spacing: f64) -> Vec<SurfacePoint3> {
        let mut sampled = Vec::new();
        for face in self.shape.triangles() {
            // If the triangle is too small, just add the center point.
            let center = mean_point(&[face.a, face.b, face.c]);
            if dist(&face.a, &center) < max_spacing
                && dist(&face.b, &center) < max_spacing
                && dist(&face.c, &center) < max_spacing
            {
                sampled.push(SurfacePoint3::new(center, face.normal().unwrap()));
                continue;
            }

            // Find the angle closest to 90 degrees
            let ua = face.b - face.a;
            let va = face.c - face.a;

            let ub = face.a - face.b;
            let vb = face.c - face.b;

            let uc = face.a - face.c;
            let vc = face.b - face.c;

            let aa = ua.angle(&va).abs() - PI / 2.0;
            let ab = ub.angle(&vb).abs() - PI / 2.0;
            let ac = uc.angle(&vc).abs() - PI / 2.0;

            let (u, v, p) = if aa < ab && aa < ac {
                (ua, va, face.a)
            } else if ab < aa && ab < ac {
                (ub, vb, face.b)
            } else {
                (uc, vc, face.c)
            };

            let nu = u.norm() / max_spacing;
            let nv = v.norm() / max_spacing;

            for ui in 0..nu as usize {
                for vi in 0..nv as usize {
                    let uf = ui as f64 / nu;
                    let vf = vi as f64 / nv;
                    if uf + vf <= 1.0 {
                        let p = p + u * uf + v * vf;
                        let sp = SurfacePoint3::new(p, face.normal().unwrap());
                        sampled.push(sp);
                    }
                }
            }
        }

        sampled
    }

    pub fn sample_alignment_candidates(&self, max_spacing: f64) -> Vec<ACPoint> {
        let surf_points = self.sample_poisson(max_spacing);
        let points = surf_points.iter().map(|sp| sp.point).collect::<Vec<_>>();
        let tree = KdTree3::new(&points);
        let mut results = Vec::new();
        for (i, sp) in surf_points.iter().enumerate() {
            let n = tree.nearest(&sp.point, NonZero::new(7).unwrap());
            let indices = n
                .iter()
                .filter_map(|(j, _)| if *j != i { Some(*j) } else { None });
            let sps = indices
                .into_iter()
                .map(|j| surf_points[j])
                .collect::<Vec<_>>();

            if sac_check(sp, &sps, max_spacing) {
                results.push(ACPoint {
                    sp: *sp,
                    neighbors: sps,
                });
            }
        }

        results
    }

    pub fn sample_alignment_points(
        &self,
        max_spacing: f64,
        reference: &Mesh,
        iso: &Iso3,
    ) -> Vec<Point3> {
        let surf_points = self.sample_poisson(max_spacing);
        let points = surf_points.iter().map(|sp| sp.point).collect::<Vec<_>>();
        let tree = KdTree3::new(&points);
        let mut candidates: Vec<Point3> = Vec::new();
        for (i, sp) in surf_points.iter().enumerate() {
            let n = tree.nearest(&sp.point, NonZero::new(7).unwrap());
            let indices = n
                .iter()
                .filter_map(|(j, _)| if *j != i { Some(*j) } else { None });
            let sps = indices
                .into_iter()
                .map(|j| surf_points[j])
                .collect::<Vec<_>>();
            if smpl_check(sp, &sps, max_spacing, reference, iso) {
                candidates.push(sp.point);
            }
        }

        // Get the distances so that we can filter all points more than 3 standard deviations
        // away from the mean.
        let distances = candidates
            .iter()
            .map(|p| dist(p, &reference.point_closest_to(&(iso * p))))
            .collect::<Vec<_>>();

        let mean_distance = distances.iter().sum::<f64>() / distances.len() as f64;
        let std_dev = (distances
            .iter()
            .map(|d| (d - mean_distance).powi(2))
            .sum::<f64>()
            / distances.len() as f64)
            .sqrt();

        candidates
            .iter()
            .zip(distances.iter())
            .filter_map(|(c, &d)| {
                if d < mean_distance + 3.0 * std_dev {
                    Some(*c)
                } else {
                    None
                }
            })
            .collect()
    }
}

/// A candidate for an alignment point on the surface of a Poisson disk sampled mesh, along with
/// its nearest neighbors.
pub struct ACPoint {
    /// The surface point at the location of the candidate.
    pub sp: SurfacePoint3,

    /// The nearest neighbors of the candidate surface point.
    pub neighbors: Vec<SurfacePoint3>,
}

fn smpl_check(
    check: &SurfacePoint3,
    neighbors: &[SurfacePoint3],
    max_spacing: f64,
    reference: &Mesh,
    iso: &Iso3,
) -> bool {
    // Actual points check
    if !sac_check(check, neighbors, max_spacing) {
        return false;
    }

    // If the points on the test mesh pass, we project the points to the reference mesh and
    // run the same check.
    let check_ref = reference.surf_closest_to(&(iso * check.point));

    // Normals must be facing the same direction
    if check_ref.normal.dot(&check.normal) < 0.0 {
        return false;
    }

    let neighbors_ref = neighbors
        .iter()
        .map(|sp| reference.surf_closest_to(&(iso * sp.point)))
        .collect::<Vec<_>>();

    // The minimum spacing to the check_ref point should be max_spacing
    for sp in &neighbors_ref {
        if dist(&sp.point, &check_ref.point) < max_spacing {
            return false;
        }
    }

    sac_check(&check_ref, &neighbors_ref, max_spacing)
}

pub fn sac_ref_check(ac: &ACPoint, max_spacing: f64, reference: &Mesh, iso: &Iso3) -> bool {
    // If the points on the test mesh pass, we project the points to the reference mesh and
    // run the same check.
    let check_ref = reference.surf_closest_to(&(iso * ac.sp.point));
    let neighbors_ref = ac
        .neighbors
        .iter()
        .map(|sp| reference.surf_closest_to(&(iso * sp.point)))
        .collect::<Vec<_>>();

    // The minimum spacing to the check_ref point should be max_spacing
    for sp in &neighbors_ref {
        if dist(&sp.point, &check_ref.point) < max_spacing {
            return false;
        }
    }

    sac_check(&check_ref, &neighbors_ref, max_spacing)
}

/// Perform a sample alignment candidate check on a single surface point and its neighbors. This
/// check looks for a minimum number of neighbors, ensures that the neighbors are within a
/// certain distance from the check point, that the normals of the neighbors are within PI/3 of the
/// test point, that the curvature is low, and that the check point is near the centroid of the
/// convex hull of the neighbors.
///
/// # Arguments
///
/// * `check_point`: the point being checked
/// * `surface_points`: a collection of surface points that are neighbors to the check point
/// * `max_spacing`: the Poisson disk sampling spacing used to create the full set of surface
///   points from which the neighbors were selected.
///
/// returns: bool
pub fn sac_check(
    check_point: &SurfacePoint3,
    surface_points: &[SurfacePoint3],
    max_spacing: f64,
) -> bool {
    if surface_points.len() < 5 {
        return false;
    }

    for sp in surface_points {
        if dist(&sp.point, &check_point.point) > max_spacing * 2.0
            || sp.normal.angle(&check_point.normal) > PI / 3.0
            || check_point.scalar_projection(&sp.point).abs() > max_spacing / 2.0
        {
            return false;
        }
    }

    let mut points = surface_points.clone_points();
    points.push(check_point.point);

    let basis = SvdBasis3::from_points(&points, None);
    let iso = Iso3::from(&basis);
    let points = (&points).transform_by(&iso);

    if points
        .iter()
        .map(|p| p.z.abs())
        .any(|z| z > max_spacing / 20.0)
    {
        return false;
    }

    let check2 = (iso * check_point.point).to_2d();
    let points2 = points.to_2d();
    let centroid = mean_point(&convex_hull(&points2));
    if dist(&centroid, &check2) > max_spacing {
        return false;
    }

    true
}
