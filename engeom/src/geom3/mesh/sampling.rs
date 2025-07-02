use super::Mesh;
use crate::common::SurfacePointCollection;
use crate::common::indices::index_vec;
use crate::common::points::{dist, mean_point};
use crate::common::poisson_disk::sample_poisson_disk;
use crate::{Point3, SurfacePoint3};
use rand::prelude::SliceRandom;
use std::f64::consts::PI;

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

        let to_take = sample_poisson_disk(&points, &indices, radius);
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
}
