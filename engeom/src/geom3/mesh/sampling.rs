use super::Mesh;
use crate::common::SurfacePointCollection;
use crate::common::indices::index_vec;
use crate::common::kd_tree::KdTreeSearch;
use crate::common::points::{dist, mean_point};
use crate::common::poisson_disk::sample_poisson_disk;
use crate::{Iso3, KdTree3, Point3, SurfacePoint3, SvdBasis3, To2D, TransformBy};
use rand::prelude::SliceRandom;
use std::f64::consts::PI;
use std::num::NonZero;
use parry2d_f64::transformation::convex_hull;
use crate::geom2::hull::convex_hull_2d;

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
        let points = starting.clone_points();
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

    pub fn sample_alignment_candidates(
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
            if smpl_check(i, &surf_points, max_spacing, reference, iso, &n) {
                candidates.push(sp.point);
            }
        }

        candidates
    }
}

fn smpl_check(
    i: usize,
    all_sps: &[SurfacePoint3],
    max_spacing: f64,
    reference: &Mesh,
    iso: &Iso3,
    n: &[(usize, f64)],
) -> bool {
    // Less than 6 points means we can't have a good alignment
    if n.len() < 6 {
        return false;
    }

    // Points must be within a certain distance and angle of the test point
    let this_sp = &all_sps[i];
    for (i, d) in n.iter() {
        if *d > max_spacing * 2.0
            || all_sps[*i].normal.angle(&this_sp.normal) > PI / 3.0
            || this_sp.scalar_projection(&all_sps[*i].point).abs() > max_spacing / 2.0
        {
            return false;
        }
    }

    let (mapped_i, points) = map_indices(i, n, all_sps);
    let basis = SvdBasis3::from_points(&points, None);
    let iso = Iso3::from(&basis);
    let points = (&points).transform_by(&iso);

    if points.iter().map(|p| p.z.abs()).any(|z| z > max_spacing / 10.0) {
        return false;
    }

    let points2 = (&points).to_2d();
    let hull = convex_hull_2d(&points2);
    if hull.contains(&mapped_i) {
        return false;
    }

    true
}

fn map_indices(i: usize, n: &[(usize, f64)], all_sps: &[SurfacePoint3]) -> (usize, Vec<Point3>) {
    let mut result = Vec::new();
    let mut mapped_index = None;
    for (j, (k, _)) in n.iter().enumerate() {
        if *k == i {
            mapped_index = Some(j);
        }
        result.push(all_sps[*k].point);
    }

    (mapped_index.unwrap(), result)
}
