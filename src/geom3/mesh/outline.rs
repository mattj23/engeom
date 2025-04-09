//! This module exists to help generate a visual outline of a mesh

use super::Mesh;
use crate::common::points::{fill_gaps, mid_point};
use crate::geom3::mesh::edges::{edge_key, naive_edges, unique_edges};
use crate::{Point3, UnitVec3};
use parry3d_f64::query::{Ray, RayCast};
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;

// struct KeyChainer {
//     forward: HashMap<u32, HashSet<u32>>,
//     reverse: HashMap<u32, HashSet<u32>>,
// }
//
// fn first_single(map: &mut HashMap<u32, HashSet<u32>>) -> Option<(u32, u32)> {
//     for (k, v) in map.iter_mut() {
//         if v.len() == 1 {
//             let k0 = *k;
//             let k1 = *v.iter().next().unwrap();
//             return Some((k0, k1));
//         }
//     }
//     None
// }
//
// impl KeyChainer {
//     pub fn new() -> Self {
//         KeyChainer {
//             forward: HashMap::new(),
//             reverse: HashMap::new(),
//         }
//     }
//
//     pub fn push(&mut self, key: &[u32; 2]) {
//         let k_min = key[0].min(key[1]);
//         let k_max = key[0].max(key[1]);
//         self.forward.entry(k_min).or_default().insert(k_max);
//         self.reverse.entry(k_max).or_default().insert(k_min);
//     }
//
//     pub fn remove_pair(&mut self, k0: u32, k1: u32) {
//         let k_min = k0.min(k1);
//         let k_max = k0.max(k1);
//
//         if let Some(v) = self.forward.get_mut(&k_min) {
//             v.remove(&k_max);
//             if v.is_empty() {
//                 self.forward.remove(&k_min);
//             }
//         }
//
//         if let Some(v) = self.reverse.get_mut(&k_max) {
//             v.remove(&k_min);
//             if v.is_empty() {
//                 self.reverse.remove(&k_max);
//             }
//         }
//     }
//
//     pub fn pop_new(&mut self) -> Option<(u32, u32)> {
//         if let Some((k0, k1)) = first_single(&mut self.forward) {
//             self.remove_pair(k0, k1);
//             Some((k0, k1))
//         } else if let Some((k0, k1)) = first_single(&mut self.reverse) {
//             self.remove_pair(k0, k1);
//             Some((k0, k1))
//         } else {
//             // Find the lowest count from forward
//             if let Some(n) = self
//                 .forward
//                 .iter()
//                 .min_by(|a, b| a.1.len().cmp(&b.1.len()))
//                 .clone()
//             {
//                 let k0 = *n.0;
//                 let k1 = *n.1.iter().next().unwrap();
//                 self.remove_pair(k0, k1);
//                 Some((k0, k1))
//             } else {
//                 None
//             }
//         }
//     }
//
//     pub fn find_next(&mut self, k0: u32) -> Option<u32> {
//         if let Some(v) = self.forward.get_mut(&k0) {
//             if v.len() == 1 {
//                 let k1 = *v.iter().next().unwrap();
//                 self.remove_pair(k0, k1);
//                 return Some(k1);
//             }
//         }
//         if let Some(v) = self.reverse.get_mut(&k0) {
//             if v.len() == 1 {
//                 let k1 = *v.iter().next().unwrap();
//                 self.remove_pair(k0, k1);
//                 return Some(k1);
//             }
//         }
//         None
//     }
//
//     pub fn sequences(&mut self) -> Vec<Vec<u32>> {
//         let mut chains = Vec::new();
//         let mut working = Vec::new();
//
//         // Find first
//         while !self.forward.is_empty() {
//             if working.is_empty() {
//                 if let Some((k0, k1)) = self.pop_new() {
//                     working.push(k0);
//                     working.push(k1);
//                 } else {
//                     for (a, b) in self.forward.iter() {
//                         println!("{} -> {:?}", a, b);
//                     }
//                     for (a, b) in self.reverse.iter() {
//                         println!("{} <- {:?}", a, b);
//                     }
//                     panic!("This should not happen");
//                 }
//             }
//
//             // Now we have a working chain
//             let last = *working.last().unwrap();
//             if let Some(next) = self.find_next(last) {
//                 self.remove_pair(last, next);
//                 working.push(next);
//             } else {
//                 // We are done with this chain
//                 chains.push(working.clone());
//                 working.clear();
//             }
//         }
//
//         chains
//     }
// }

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

            let p0: Point3 =
                self.shape.vertices()[k0 as usize].clone() + vert_normals[k0 as usize] * 1e-2;
            let p1: Point3 =
                self.shape.vertices()[k1 as usize].clone() + vert_normals[k1 as usize] * 1e-2;

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

        // // We want to chain vertices together
        // let sequences = working
        //     .sequences()
        //     .into_iter()
        //     .map(|indices| {
        //         indices
        //             .iter()
        //             .map(|i| self.shape.vertices()[*i as usize])
        //             .collect::<Vec<_>>()
        //     })
        //     .collect::<Vec<_>>();
        //
        // let mut edges = Vec::new();
        // let shift = -facing.into_inner();
        // for c in sequences.iter() {
        //     if let Ok(chain) = Curve3::from_points(c, 1e-3) {
        //         if chain.count() < 3 || chain.length() < max_edge_length * 2.0  {
        //             continue;
        //         }
        //
        //         // println!("before l={} n={}", chain.length(), chain.count());
        //         let chain = chain.resample(Resample::ByMaxSpacing(max_edge_length));
        //         // println!("after");
        //
        //         for (p0, p1) in chain.points().iter().zip(chain.points().iter().skip(1)) {
        //             let p = mid_point(p0, p1) + shift * 1e-2;
        //
        //             let ray = Ray::new(p, shift);
        //
        //             if self.shape.intersects_local_ray(&ray, f64::MAX) {
        //                 edges.push((*p0, *p1, 1))
        //             } else {
        //                 edges.push((*p0, *p1, 0))
        //             }
        //         }
        //     }
        // }

        edges
    }

    fn classified_edge_types(&self) -> (HashSet<[u32; 2]>, HashMap<[u32; 2], [u32; 2]>) {
        let naive = naive_edges(&self.shape.indices());
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
