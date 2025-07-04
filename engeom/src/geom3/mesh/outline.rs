//! This module exists to help generate a visual outline of a mesh

use super::Mesh;
use crate::common::points::{fill_gaps, mid_point};
use crate::geom3::mesh::edges::{edge_key, naive_edges, unique_edges};
use crate::{Point3, UnitVec3};
use parry3d_f64::query::{Ray, RayCast};
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
type CEdgeTypes = (HashSet<[u32; 2]>, HashMap<[u32; 2], [u32; 2]>);

impl Mesh {
    pub fn visual_outline(
        &self,
        facing: UnitVec3,
        max_edge_length: f64,
        corner_angle: Option<f64>,
    ) -> Vec<(Point3, Point3, u8)> {
        let corner_angle = corner_angle.unwrap_or(PI / 4.0 - 1e-2);
        let (boundaries, mut corners) = self.classified_edge_types();
        // let mut working = KeyChainer::new();
        let mut working = Vec::new();

        for (i, indices) in self.shape.indices().iter().enumerate() {
            for (i0, i1) in [(0, 1), (1, 2), (2, 0)] {
                let k = edge_key(&[indices[i0], indices[i1]]);

                if boundaries.contains(&k) {
                    working.push(k);
                } else if let Some(corner) = corners.get_mut(&k) {
                    if corner[0] == u32::MAX {
                        corner[0] = i as u32;
                    } else {
                        corner[1] = i as u32;
                    }
                }
            }
        }

        // At this point, working contains boundary edges and corners contains corner face pairs
        // Now we need to process the corners
        for (key, corner) in corners.iter() {
            if corner[0] == u32::MAX || corner[1] == u32::MAX {
                continue;
            }

            let n0u = self.shape.triangle(corner[0]).normal();
            let n1u = self.shape.triangle(corner[1]).normal();

            if let (Some(n0), Some(n1)) = (n0u, n1u) {
                if n0.angle(&n1) > corner_angle {
                    // Is this a corner?
                    working.push(*key);
                } else {
                    let f0 = facing.dot(&n0);
                    let f1 = facing.dot(&n1);
                    let f_max = f0.max(f1);
                    let f_min = f0.min(f1);

                    if f_max >= 0.0 && f_min < 0.0 {
                        // Is this a silhouette?
                        working.push(*key);
                    }
                }
            }
        }

        let vert_normals = self.get_vertex_normals();
        let mut edges = Vec::new();
        for k in working {
            let k0 = k[0];
            let k1 = k[1];

            let p0: Point3 = self.shape.vertices()[k0 as usize] + vert_normals[k0 as usize] * 1e-2;
            let p1: Point3 = self.shape.vertices()[k1 as usize] + vert_normals[k1 as usize] * 1e-2;

            let points = fill_gaps(&[p0, p1], max_edge_length);

            for (p0, p1) in points.iter().zip(points.iter().skip(1)) {
                let p = mid_point(p0, p1) + facing.into_inner() * 1e-2;

                let ray = Ray::new(p, facing.into_inner());

                if self.shape.intersects_local_ray(&ray, f64::MAX) {
                    edges.push((*p0, *p1, 1))
                } else {
                    edges.push((*p0, *p1, 0))
                }
            }
        }

        edges
    }

    fn classified_edge_types(&self) -> CEdgeTypes {
        let naive = naive_edges(self.shape.indices());
        let unique = unique_edges(&naive);

        let mut boundaries = HashSet::new();
        let mut corners = HashMap::new();

        for (key, count) in unique {
            if count == 1 {
                boundaries.insert(key);
            } else if count == 2 {
                corners.insert(key, [u32::MAX, u32::MAX]);
            }
        }

        (boundaries, corners)
    }
}
