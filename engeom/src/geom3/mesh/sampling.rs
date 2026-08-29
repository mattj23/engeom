use super::{Mesh3, MeshSurfPoint};
use crate::common::IndexMask;
use crate::common::barycentric::barycentric_grid;
use crate::common::points::dist;
use crate::common::poisson_disk::sample_poisson_disk_all;
use crate::{Point3, SurfacePoint3};

impl Mesh3 {
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

    pub fn sample_poisson(
        &self,
        radius: f64,
        index_mask: Option<&IndexMask>,
    ) -> Vec<MeshSurfPoint> {
        let starting = self.sample_surface_dense(radius * 0.5, index_mask);
        let mask = sample_poisson_disk_all(&starting, radius);

        let mut result = Vec::new();
        for i in mask.iter_true() {
            result.push(starting[i]);
        }

        result
    }

    /// Generate a dense sampling of the surface of a mesh, where points are guaranteed to be no
    /// more than `max_spacing` units apart but have no guarantees about how close together they
    /// may be.
    ///
    /// TODO: test these guarantees
    ///
    /// # Arguments
    ///
    /// * `max_spacing`: the maximum distance between any two points in the sampling
    /// * `index_mask`: Optional mask to apply to the sampling, only sampling faces that are true
    ///   in the mask
    ///
    /// returns: Vec<MeshSurfPoint, Global>
    pub fn sample_surface_dense(
        &self,
        max_spacing: f64,
        index_mask: Option<&IndexMask>,
    ) -> Vec<MeshSurfPoint> {
        let mut sampled = Vec::with_capacity(self.faces().len());

        for (face_i, vert) in self.faces().iter().enumerate() {
            if let Some(index_mask) = index_mask
                && !index_mask.get(face_i)
            {
                continue;
            }

            let a = self.points()[vert[0] as usize];
            let b = self.points()[vert[1] as usize];
            let c = self.points()[vert[2] as usize];
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

    pub fn sample_dense(
        &self,
        max_spacing: f64,
        index_mask: Option<&IndexMask>,
    ) -> Vec<SurfacePoint3> {
        self.sample_surface_dense(max_spacing, index_mask)
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
    //     reference: &Mesh3,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::KdTree3;
    use crate::common::kd_tree::*;

    #[test]
    fn check_kiddo_bug() {
        let mesh = Mesh3::create_sphere(100.0, 0.011).unwrap();
        let r = 5.0;
        let sampled = mesh.sample_poisson(r, None);

        let points = sampled.iter().map(|mp| mp.sp.point).collect::<Vec<_>>();

        let tree = KdTree3::try_new(&points).expect("Tree construction failed");
        for mp in &sampled {
            let neighbors = tree.within(&mp.sp.point, r);
            assert_eq!(neighbors.len(), 1, "Missed duplicate");
        }
    }

    #[test]
    fn sample_poisson_index_mask_restricts_to_masked_faces() {
        let mesh = Mesh3::create_sphere(100.0, 1.1).unwrap();
        let n_faces = mesh.faces().len();

        // Mask only the first half of the faces
        let masked_indices: Vec<usize> = (0..n_faces / 2).collect();
        let mask = IndexMask::try_from_indices(&masked_indices, n_faces).unwrap();

        let sampled = mesh.sample_poisson(5.0, Some(&mask));

        assert!(!sampled.is_empty(), "Expected at least one sample");
        for mp in &sampled {
            assert!(
                mask.get(mp.face_index as usize),
                "face_index {} is not in the mask",
                mp.face_index
            );
        }
    }
}
