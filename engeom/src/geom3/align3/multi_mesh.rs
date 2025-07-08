//! This module contains an optimization for multiple meshes to each other in one combined
//! Levenberg-Marquardt minimization.  This is different from a pose graph optimization because it
//! contains only a single transformation for each mesh, minus the first.
//!
//! In general this technique will be stable and produce extremely accurate results for high
//! quality, low-noise meshes which have already been pre-aligned to be close to each other
//! with a relatively large amount of overlap between meshes.  This code was implemented to perform
//! bundle adjustment between metrology quality scans of objects with unambiguous morphology.

use crate::common::points::dist;
use crate::Result;
use crate::geom3::Align3;
use crate::geom3::align3::jacobian::{point_plane_jacobian, point_plane_jacobian_rev};
use crate::geom3::align3::multi_param::ParamHandler;
use crate::geom3::align3::{distance_weight, normal_weight};
use crate::na::{Dyn, Matrix, Owned, U1, Vector, DMatrix};
use crate::{Iso3, Mesh, Point3, SurfacePoint3};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use crate::geom3::mesh::sampling::sac_ref_check;

/// Options for the multi-mesh simultaneous alignment algorithm.
#[derive(Debug, Clone, Copy)]
pub struct MMOpts {
    pub search_radius: f64,
    pub sample_radius: f64,
    pub respect_normals: bool,
}

impl MMOpts {
    pub fn new(search_radius: f64, sample_radius: f64, respect_normals: bool) -> Self {
        Self {
            search_radius,
            sample_radius,
            respect_normals,
        }
    }
}

pub fn multi_mesh_adjustment(
    meshes: &[Mesh],
    opts: MMOpts,
    initial: Option<&[Iso3]>,
) -> Result<Vec<Align3>> {
    // Produce the sample candidate points
    let sample_candidates = meshes
        .iter()
        .map(|m| m.sample_alignment_candidates(opts.sample_radius))
        .collect::<Vec<_>>();

    let transforms = match initial {
        Some(transforms) => transforms.to_vec(),
        None => vec![Iso3::identity(); meshes.len()],
    };

    // We want to build a correspondence matrix which will help us determine which mesh will be the
    // static reference mesh.  In the matrix, each i, j entry will be the number of sample points
    // in mesh j which are a good match for mesh i.  The row with the highest sum of its columns
    // has the most points which reference it, however we don't simply want to find the highest
    // count because two very overlapping meshes will inflate the numbers without being a good
    // candidate for the static reference.
    //
    // Instead, we want to preference meshes which have a higher number of other meshes that
    // reference them, so we will divide each cell by the maximum value in the matrix, scaling
    // everything from 0 to 1, and then take the square root of each cell, essentially granting
    // diminishing returns to the number of points from a single mesh and boosting meshes that
    // have a moderate number of points from a large number of other meshes.
    let mut matrix = DMatrix::<f64>::zeros(meshes.len(), meshes.len());

    let mut handles = Vec::new();

    for i in 0..meshes.len() {
        for j in 0..meshes.len() {
            if i == j {
                continue;
            }

            // Get the relative transform between the meshes, allowing a point in the test mesh
            // to be transformed to the reference mesh with the same distances between points as
            // if both clouds were transformed by their respective transforms.
            let t = transforms[i].inv_mul(&transforms[j]);

            for (si, c) in sample_candidates[j].iter().enumerate() {
                if sac_ref_check(c, opts.sample_radius, &meshes[i], &t) {
                    // If the candidate point is a good match for the reference mesh, then we
                    // increment the count in the matrix
                    matrix[(i, j)] += 1.0;

                    // We'll also create a handle for this point, which will be used later assuming
                    // that mesh j does not end up being the static reference mesh.
                    handles.push(TestPoint {
                        mesh_i: j,
                        point_i: si,
                        ref_i: i,
                    })
                }
            }
        }
    }


    let mut corr = &matrix / matrix.max();
    corr.apply(|x| *x = x.sqrt());

    // Now we want to sort the column sums in reverse order to get the reference priority. The
    // first element in the list will be the static reference cloud.
    let mut corr_pairs = corr
        .column_sum()
        .iter()
        .enumerate()
        .map(|(i, x)| (i, *x))
        .collect::<Vec<_>>();
    corr_pairs.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let reference_order = corr_pairs.iter().map(|(i, _)| *i).collect::<Vec<_>>();
    let static_i = reference_order[0];
    println!("correspondence matrix: {}", matrix);
    println!("static_i: {}", static_i);
    println!("corr: {:?}", corr_pairs);

    todo!("i'll take it from here!")
    /*
    // Now we want to generate the test points. Each test point is a point in a point cloud which
    // is being matched to another point cloud.  We want to generate these such that for each
    // unique pair of clouds there are only test points which go from one cloud to the other, and
    // none which go in reverse.
    // let start = Instant::now();
    // let mut test_points = Vec::new();
    // let mut clouds_to_test = (0..clouds.len()).collect::<Vec<_>>();
    // for i in reference_order {
    //     // Remove the current reference cloud from the list of clouds to test
    //     clouds_to_test.retain(|j| *j != i);
    //
    //     // Get all the clouds which reference the current working reference cloud and create
    //     // test points for them
    //     for test_i in &clouds_to_test {
    //         for j in overlap.indices[&(i, *test_i)].iter() {
    //             test_points.push(TestPoint::new(*test_i, *j, i));
    //         }
    //     }
    // }

    // println!("test_points: {:?}", start.elapsed());
    // println!("test_points: {:?}", test_points.len());

    // Now we want to create the problem and solve it
    // let start = Instant::now();
    let smpl = sample_candidates
        .into_iter()
        .map(|row| row.iter().map(|c| c.sp).collect::<Vec<_>>())
        .collect::<Vec<_>>();

    let problem = MultiMeshProblem::new(meshes, smpl, handles, static_i, opts, initial);
    let (result, report) = LevenbergMarquardt::new().minimize(problem);
    // println!("minimize: {:?}", start.elapsed());
    if report.termination.was_successful() {
        let alignments = (0..meshes.len())
            .map(|i| result.params.get_transform(i))
            .collect::<Vec<_>>();

        let mut grouped = (0..meshes.len()).map(|_| Vec::new()).collect::<Vec<_>>();
        let residuals = result.residuals().unwrap();

        for (i, p) in result.point_handles.iter().enumerate() {
            grouped[p.mesh_i].push(residuals[i]);
        }

        Ok(alignments
            .iter()
            .zip(grouped.into_iter())
            .map(|(a, g)| Align3::new(*a, g))
            .collect())
    } else {
        println!("{:?}", report.termination);
        Err("Failed to converge".into())
    }
     */
}

struct TestPoint {
    /// The index of the mesh this point belongs to
    mesh_i: usize,

    /// The index of the point in the alignment points denoted by `mesh_i`
    point_i: usize,

    /// The index of the mesh this point is being matched to
    ref_i: usize,
}

impl TestPoint {
    fn new(mesh_i: usize, point_i: usize, ref_i: usize) -> Self {
        Self {
            mesh_i,
            point_i,
            ref_i,
        }
    }
}

struct MultiMeshProblem<'a> {
    /// The meshes that are being aligned
    meshes: &'a [Mesh],

    /// The collections of sample points (one for each mesh) that are used to align the meshes to
    /// each other. The i-th collection contains the points from the i-th mesh.
    sample_points: Vec<Vec<SurfacePoint3>>,

    /// The collection of alignment point handles, each specifying which mesh/index it belongs to,
    /// and which mesh it is being matched to.  This allows the same point to be used against
    /// multiple targets.
    point_handles: Vec<TestPoint>,

    /// The collection of sample points after they've been moved by the optimizer. These correspond
    /// to the point handles, and are used to compute the residuals.
    moved: Vec<SurfacePoint3>,

    /// A collection of the closest points on the mesh surfaces which correspond with the point
    /// handles. The i-th entry corresponds to the i-th point handle.
    closest: Vec<SurfacePoint3>,

    /// A collection of weights for each point handle, which is used to scale the residuals. The
    /// i-th entry corresponds to the i-th point handle.
    weight: Vec<f64>,

    /// The internal parameters for the optimization, which is an encoding of the relative
    /// transformations between the meshes.
    params: ParamHandler,

    /// The options for the multi-mesh optimization, such as the search radius and sample radius.
    options: MMOpts,
}

impl<'a> MultiMeshProblem<'a> {
    fn new(
        meshes: &'a [Mesh],
        sample_points: Vec<Vec<SurfacePoint3>>,
        point_handles: Vec<TestPoint>,
        static_i: usize,
        options: MMOpts,
        initial: Option<&[Iso3]>,
    ) -> Self {
        let mean_points = meshes.iter().map(|m| m.aabb().center()).collect::<Vec<_>>();
        let params = ParamHandler::new(static_i, mean_points, initial);

        let mut item = Self {
            meshes,
            sample_points,
            point_handles,
            moved: Vec::new(),
            closest: Vec::new(),
            weight: Vec::new(),
            params,
            options,
        };

        item.move_points();
        item
    }

    fn move_points(&mut self) {
        self.moved.clear();
        self.closest.clear();
        self.weight.clear();

        for h in self.point_handles.iter() {
            // let mesh = &self.meshes[h.mesh_i];
            let point = &self.sample_points[h.mesh_i][h.point_i];
            // let normal = cloud.normals().unwrap()[p.point_i];
            let ref_cloud = &self.meshes[h.ref_i];
            let t = self.params.relative_transform(h.mesh_i, h.ref_i);

            let moved = &t * point;
            let closest = ref_cloud.surf_closest_to(&moved.point);

            let w0 = distance_weight(
                dist(&moved.point, &closest.point),
                self.options.search_radius,
            );

            let w1 = normal_weight(&moved.normal.into_inner(), &closest.normal.into_inner());

            self.closest.push(closest);
            self.moved.push(moved);
            self.weight.push(w0 * w1);
        }
    }
}

impl<'a> LeastSquaresProblem<f64, Dyn, Dyn> for MultiMeshProblem<'a> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    type ParameterStorage = Owned<f64, Dyn>;

    fn set_params(&mut self, x: &Vector<f64, Dyn, Self::ParameterStorage>) {
        self.params.set_param(x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, Dyn, Self::ParameterStorage> {
        (*self.params.params()).clone()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res =
            Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.point_handles.len());

        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            // Within this code block:
            //  - i is the index of the handle and the residual we're working on
            //  - p is the moved sample surface point
            //  - c is the closest point
            let v = p.point - c.point;
            let d = v.dot(&c.normal);
            res[i] = self.weight[i] * d.abs();
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, Dyn, Self::JacobianStorage>> {
        // rows are each of the test samples
        // columns are the parameters
        let mut jac = Matrix::<f64, Dyn, Dyn, Self::JacobianStorage>::zeros(
            self.point_handles.len(),
            self.params.params().len(),
        );

        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            // p is the moved surface point
            // c is the closest surface point on the reference
            let handle = &self.point_handles[i];
            let test_i = handle.mesh_i;
            let ref_i = handle.ref_i;

            // values0 contains the derivatives of the residual with respect to the transform
            // of the test cloud
            let values0 = point_plane_jacobian(&p.point, &c, &self.params.params[test_i]);
            self.params
                .set_jacobian(&mut jac, i, test_i, &(values0 * self.weight[i]));

            // values1 contains the derivatives of the residual with respect to the transform
            // of the reference cloud
            let values1 = point_plane_jacobian_rev(&p.point, &c, &self.params.params[ref_i]);
            self.params
                .set_jacobian(&mut jac, i, ref_i, &(values1 * self.weight[i]));
        }

        Some(jac)
    }
}

// /// Calculates a matrix containing "correspondences" between point clouds, where the number of
// /// points in cloud j which have a distance less than `search_radius` from a point in cloud i is
// /// stored in the (i, j) entry of the matrix
// pub fn correspondences(clouds: &[PointCloudWithTree], search_radius: f64) -> Result<DMatrix<f64>> {
//     // First, we need to find the point cloud that has the most points referencing it. That one
//     // will be the static reference
//     let mut references = DMatrix::zeros(clouds.len(), clouds.len());
//     for (i, c_ref) in clouds.iter().enumerate() {
//         if c_ref.points().is_empty() {
//             return Err("Point cloud has no points".into());
//         }
//
//         for (j, c_test) in clouds.iter().enumerate() {
//             if i == j {
//                 continue;
//             }
//             references[(i, j)] = c_ref.matching_point_count(c_test.points(), search_radius) as f64;
//         }
//     }
//
//     Ok(references)
// }
//
