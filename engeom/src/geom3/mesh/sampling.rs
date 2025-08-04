use super::{Mesh, MeshSurfPoint};
use crate::common::linear_space;
use crate::common::points::dist;
use crate::common::poisson_disk::sample_poisson_disk_all;
use crate::{Point3, SurfacePoint3};

fn bc_to_point(bc: [f64; 3], a: &Point3, b: &Point3, c: &Point3) -> Point3 {
    (a.coords * bc[0] + b.coords * bc[1] + c.coords * bc[2]).into()
}

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

    pub fn sample_poisson(&self, radius: f64) -> Vec<MeshSurfPoint> {
        let starting = self.sample_surface_dense(radius * 0.5);
        let mask = sample_poisson_disk_all(&starting, radius);

        let mut result = Vec::new();
        for i in mask.iter_true() {
            result.push(starting[i]);
        }

        result
    }

    pub fn sample_surface_dense(&self, max_spacing: f64) -> Vec<MeshSurfPoint> {
        let mut sampled = Vec::with_capacity(self.faces().len());
        for (face_i, vert) in self.faces().iter().enumerate() {
            let a = self.vertices()[vert[0] as usize];
            let b = self.vertices()[vert[1] as usize];
            let c = self.vertices()[vert[2] as usize];
            let face_index = face_i as u32;

            if dist(&a, &b) < max_spacing
                && dist(&a, &c) < max_spacing
                && dist(&b, &c) < max_spacing
            {
                // If all distances between vertices are less than the max spacing, we'll just
                // sample the centroid of the face, as an equally sized neighbor should have its
                // centroid within the max spacing distance of this triangle's centroid.
                let bc = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0];
                let sp = self.at_barycentric(face_index, bc).unwrap().sp;
                sampled.push(MeshSurfPoint { face_index, bc, sp });
            } else {
                let grid = barycentric_grid(&a, &b, &c, max_spacing);
                for bc in grid {
                    let sp = self.at_barycentric(face_index, bc).unwrap().sp;
                    sampled.push(MeshSurfPoint { face_index, bc, sp });
                }
            }
        }

        sampled
    }

    pub fn sample_dense(&self, max_spacing: f64) -> Vec<SurfacePoint3> {
        self.sample_surface_dense(max_spacing)
            .into_iter()
            .map(|msp| msp.sp)
            .collect()
    }

    // pub fn sample_alignment_candidates(&self, max_spacing: f64) -> Vec<ACPoint> {
    //     let surf_points = self.sample_poisson(max_spacing);
    //     let points = surf_points.iter().map(|sp| sp.point).collect::<Vec<_>>();
    //     let tree = KdTree3::new(&points);
    //     let mut results = Vec::new();
    //     for (i, sp) in surf_points.iter().enumerate() {
    //         let n = tree.nearest(&sp.point, NonZero::new(7).unwrap());
    //         let indices = n
    //             .iter()
    //             .filter_map(|(j, _)| if *j != i { Some(*j) } else { None });
    //         let sps = indices
    //             .into_iter()
    //             .map(|j| surf_points[j])
    //             .collect::<Vec<_>>();
    //
    //         if sac_check(sp, &sps, max_spacing) {
    //             results.push(ACPoint {
    //                 sp: *sp,
    //                 neighbors: sps,
    //             });
    //         }
    //     }
    //
    //     results
    // }
    //
    // pub fn sample_alignment_points(
    //     &self,
    //     max_spacing: f64,
    //     reference: &Mesh,
    //     iso: &Iso3,
    // ) -> Vec<SurfacePoint3> {
    //     let surf_points = self.sample_poisson(max_spacing);
    //
    //     let points = surf_points.iter().map(|sp| sp.point).collect::<Vec<_>>();
    //     let tree = KdTree3::new(&points);
    //
    //     let mut candidates: Vec<SurfacePoint3> = Vec::new();
    //     for (i, sp) in surf_points.iter().enumerate() {
    //         let n = tree.nearest(&sp.point, NonZero::new(7).unwrap());
    //         let indices = n
    //             .iter()
    //             .filter_map(|(j, _)| if *j != i { Some(*j) } else { None });
    //         let sps = indices
    //             .into_iter()
    //             .map(|j| surf_points[j])
    //             .collect::<Vec<_>>();
    //         if smpl_check(sp, &sps, max_spacing, reference, iso) {
    //             candidates.push(*sp);
    //         }
    //     }
    //
    //     // Get the distances so that we can filter all points more than 3 standard deviations
    //     // away from the mean.
    //     let distances = candidates
    //         .iter()
    //         .map(|p| dist(&p.point, &reference.point_closest_to(&(iso * p.point))))
    //         .collect::<Vec<_>>();
    //
    //     let mean_distance = distances.iter().sum::<f64>() / distances.len() as f64;
    //     let std_dev = (distances
    //         .iter()
    //         .map(|d| (d - mean_distance).powi(2))
    //         .sum::<f64>()
    //         / distances.len() as f64)
    //         .sqrt();
    //
    //     candidates
    //         .iter()
    //         .zip(distances.iter())
    //         .filter_map(|(c, &d)| {
    //             if d < mean_distance + 3.0 * std_dev {
    //                 Some(*c)
    //             } else {
    //                 None
    //             }
    //         })
    //         .collect()
    // }
}

// /// A candidate for an alignment point on the surface of a Poisson disk sampled mesh, along with
// /// its nearest neighbors.
// pub struct ACPoint {
//     /// The surface point at the location of the candidate.
//     pub sp: SurfacePoint3,
//
//     /// The nearest neighbors of the candidate surface point.
//     pub neighbors: Vec<SurfacePoint3>,
// }

fn barycentric_grid(a: &Point3, b: &Point3, c: &Point3, max_spacing: f64) -> Vec<[f64; 3]> {
    let mut result = Vec::new();
    let va = a - bc_to_point([0.0, 0.5, 0.5], a, b, c);
    let vb = b - bc_to_point([0.5, 0.0, 0.5], a, b, c);
    let vc = c - bc_to_point([0.5, 0.5, 0.0], a, b, c);

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
    use crate::KdTree3;
    use crate::common::kd_tree::*;
    use crate::common::points::evenly_spaced_points_between;

    #[test]
    fn check_kiddo_bug() {
        let mesh = Mesh::create_sphere(100.0, 300, 300);
        let r = 5.0;
        let sampled = mesh.sample_poisson(r);

        let points = sampled.iter().map(|mp| mp.sp.point).collect::<Vec<_>>();

        let tree = KdTree3::new(&points);
        for mp in &sampled {
            let neighbors = tree.within(&mp.sp.point, r);
            assert_eq!(neighbors.len(), 1, "Missed duplicate");
        }
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
            .map(|bc| bc_to_point(*bc, &a, &b, &c))
            .collect::<Vec<_>>();

        // Check for NAN
        for point in &grid_points {
            assert!(
                !point.coords.iter().any(|&x| x.is_nan()),
                "Point contains NaN: {:?}",
                point
            );
        }

        let tree = KdTree::new(&grid_points);

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
