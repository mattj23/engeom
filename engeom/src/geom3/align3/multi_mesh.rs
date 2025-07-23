//! This module contains an optimization for multiple meshes to each other in one combined
//! Levenberg-Marquardt minimization.  This is different from a pose graph optimization because it
//! contains only a single transformation for each mesh, minus the first.
//!
//! In general this technique will be stable and produce extremely accurate results for high
//! quality, low-noise meshes which have already been pre-aligned to be close to each other
//! with a relatively large amount of overlap between meshes.  This code was implemented to perform
//! bundle adjustment between metrology quality scans of objects with unambiguous morphology.

use crate::Result;
use crate::common::points::dist;
use crate::geom3::Align3;
use crate::geom3::align3::jacobian::{point_plane_jacobian, point_plane_jacobian_rev};
use crate::geom3::align3::multi_param::ParamHandler;
use crate::geom3::align3::{GAPParams, distance_weight, normal_weight};
use crate::na::{DMatrix, Dyn, Matrix, Owned, U1, Vector};
use faer::prelude::default;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use rayon::prelude::*;
use std::time::Instant;

use crate::geom3::align3::mesh::{AlignmentMesh, generate_alignment_points};
use crate::geom3::mesh::MeshSurfPoint;

/// Options for the multi-mesh simultaneous alignment algorithm.
#[derive(Debug, Clone, Copy)]
pub struct MMOpts {
    pub search_radius: f64,
    pub respect_normals: bool,
    pub sample: GAPParams,
}

impl MMOpts {
    pub fn new(search_radius: f64, respect_normals: bool, sample: GAPParams) -> Self {
        Self {
            search_radius,
            respect_normals,
            sample,
        }
    }
}

pub fn multi_mesh_adjustment(meshes: &[AlignmentMesh], opts: MMOpts) -> Result<Vec<Align3>> {
    let matrix = correspondence_matrix(meshes, &opts.sample);
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
    // println!("correspondence matrix: {}", matrix);
    // println!("static_i: {}", static_i);
    // println!("corr: {:?}", corr_pairs);

    // Now we want to generate the test points. Each test point is a point in a point cloud which
    // is being matched to another point cloud.  We want to generate these such that for each
    // unique pair of clouds there are only test points which go from one cloud to the other, and
    // none which go in reverse.
    let mut work_list = Vec::new();
    let mut meshes_to_test = (0..meshes.len()).collect::<Vec<_>>();
    for ref_i in reference_order {
        // Remove the current reference cloud from the list of clouds to test
        meshes_to_test.retain(|j| *j != ref_i);
        // Get all the clouds which reference the current working reference cloud and create
        // test points for them
        for &mesh_i in meshes_to_test.iter() {
            work_list.push((mesh_i, ref_i));
        }
    }

    let start = Instant::now();
    let handles = work_list
        .par_iter()
        .map(|(mesh_i, ref_i)| {
            let t = meshes[*ref_i]
                .transform()
                .inv_mul(&meshes[*mesh_i].transform());
            let mut to_test = Vec::new();
            let samples = generate_alignment_points(
                meshes[*mesh_i].mesh,
                meshes[*ref_i].mesh,
                &t,
                &opts.sample,
            );
            for mp in samples {
                // Check if there is weight associated with this point
                let weight = if let Some(providers) = &meshes[*mesh_i].weights {
                    let mut w = 1.0;
                    for item in providers.iter() {
                        w *= item.weight(&mp);
                    }
                    w
                } else {
                    1.0
                };

                to_test.push(TestPoint::new(*mesh_i, mp, *ref_i, weight, 0.0));
            }
            to_test
        })
        .flatten()
        .collect::<Vec<_>>();

    // println!("test_points: {:?}", start.elapsed());
    // println!("handles: {:?}", handles.len());
    let weighted_count = handles.iter().filter(|h| h.weight > 1.01).count();
    // println!("weighted n={weighted_count}");

    // Now we want to create the problem and solve it
    let start = Instant::now();

    let problem = MultiMeshProblem::new(meshes, handles, static_i, opts);
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
            .zip(grouped)
            .map(|(a, g)| Align3::new(*a, g))
            .collect())
    } else {
        // println!("{:?}", report.termination);
        Err("Failed to converge".into())
    }
}

struct TestPoint {
    /// The index of the mesh this point belongs to
    pub mesh_i: usize,

    pub mp: MeshSurfPoint,

    /// The index of the mesh this point is being matched to
    pub ref_i: usize,

    /// The base weight for this point, which is used to scale the residuals
    pub weight: f64,

    #[allow(dead_code)]
    pub uncert: f64,
}

impl TestPoint {
    fn new(mesh_i: usize, mp: MeshSurfPoint, ref_i: usize, weight: f64, uncert: f64) -> Self {
        Self {
            mesh_i,
            mp,
            ref_i,
            weight,
            uncert,
        }
    }
}

struct MultiMeshProblem<'a> {
    /// The meshes that are being aligned
    meshes: &'a [AlignmentMesh<'a>],

    /// The collection of alignment point handles, each specifying which mesh/index it belongs to,
    /// and which mesh it is being matched to.  This allows the same point to be used against
    /// multiple targets.
    point_handles: Vec<TestPoint>,

    /// The collection of sample points after they've been moved by the optimizer. These correspond
    /// to the point handles, and are used to compute the residuals.
    moved: Vec<MeshSurfPoint>,

    /// A collection of the closest points on the mesh surfaces which correspond with the point
    /// handles. The i-th entry corresponds to the i-th point handle.
    closest: Vec<MeshSurfPoint>,

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
        meshes: &'a [AlignmentMesh],
        point_handles: Vec<TestPoint>,
        static_i: usize,
        options: MMOpts,
    ) -> Self {
        let mean_points = meshes
            .iter()
            .map(|m| m.mesh.aabb().center())
            .collect::<Vec<_>>();
        let initial = meshes.iter().map(|m| m.transform()).collect::<Vec<_>>();

        let params = ParamHandler::new(static_i, mean_points, Some(&initial));
        let count: usize = point_handles.len();

        // TODO: we can pre-compute the uncertainties for point handles
        // If uncertainties are provided, we need to measure the uncertainties of each test point
        // let handle_uncertainties = if let Some(u) = uncertainties {
        //     let mut working = vec![0.0; count];
        //     for (i, h) in point_handles.iter().enumerate() {
        //         let m = &meshes[h.mesh_i];
        //         working[i] = un_cert(&h.sp.point, m, &u[h.mesh_i]).0;
        //     }
        //     working
        // } else {
        //     Vec::new()
        // };

        let mut item = Self {
            meshes,
            point_handles,
            moved: vec![default(); count],
            closest: vec![default(); count],
            weight: vec![0.0; count],
            params,
            options,
        };

        item.move_points();
        item
    }

    fn move_points(&mut self) {
        let indices = (0..self.point_handles.len()).collect::<Vec<_>>();
        let collected = indices
            .par_iter()
            .map(|i| {
                let h = &self.point_handles[*i];
                let ref_mesh = &self.meshes[h.ref_i];
                let t = self.params.relative_transform(h.mesh_i, h.ref_i);

                let moved = h.mp.transformed_by(&t);
                let closest = ref_mesh.mesh.surf_closest_to(&moved.point());

                // TODO: Uncertainties here for closest

                // let (wu, closest) = match self.uncertainties {
                //     Some(u) => {
                //         // Get the standard deviation (uncertainty) for the point along with the
                //         // actual closest surface point/normal on the reference mesh.
                //         let (sd, sp) = un_cert(&moved.point, ref_mesh, &u[h.ref_i]);
                //         let hsd = self.handle_uncertainties[*i];
                //
                //         // Now we combine the uncertainties and get the peak height of a gaussian
                //         let peak = 1.0 / (2.0 * PI * (sd + hsd)).sqrt();
                //         (peak, sp)
                //     }
                //     None => (1.0, ref_mesh.surf_closest_to(&moved.point).sp),
                // };

                // Calculate the different weights
                let weight_d = distance_weight(dist(&moved, &closest), self.options.search_radius);
                let weight_n = if self.options.respect_normals {
                    normal_weight(&moved.normal(), &closest.normal())
                } else {
                    1.0
                };

                let weight = h.weight * weight_d * weight_n;

                (*i, moved, closest, weight)
            })
            .collect::<Vec<_>>();

        for (i, moved, closest, w) in collected {
            self.moved[i] = moved;
            self.closest[i] = closest;
            self.weight[i] = w;
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
            let v = p.point() - c.point();
            let d = v.dot(&c.normal());
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
            let values0 = point_plane_jacobian(&p.point(), &c.sp, &self.params.params[test_i]);
            self.params
                .set_jacobian(&mut jac, i, test_i, &(values0 * self.weight[i]));

            // values1 contains the derivatives of the residual with respect to the transform
            // of the reference cloud
            let values1 = point_plane_jacobian_rev(&p.point(), &c.sp, &self.params.params[ref_i]);
            self.params
                .set_jacobian(&mut jac, i, ref_i, &(values1 * self.weight[i]));
        }

        Some(jac)
    }
}

// fn un_cert(point: &Point3, target: &Mesh, by_vertex: &[f64]) -> (f64, SurfacePoint3) {
//     let (prj, (fi, tpl)) = target
//         .tri_mesh()
//         .project_local_point_and_get_location(point, false);
//
//     let face = target.faces()[fi as usize];
//     let normal = target.tri_mesh().triangle(fi).normal().unwrap();
//
//     let u0 = by_vertex[face[0] as usize];
//     let u1 = by_vertex[face[1] as usize];
//     let u2 = by_vertex[face[2] as usize];
//
//     let st_dev = if let Some(bc) = tpl.barycentric_coordinates() {
//         u0 * bc[0] + u1 * bc[1] + u2 * bc[2]
//     } else {
//         u0.max(u1).max(u2)
//     };
//
//     (st_dev.powi(2), SurfacePoint3::new(prj.point, normal))
// }

fn correspondence_matrix(meshes: &[AlignmentMesh], params: &GAPParams) -> DMatrix<f64> {
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

    let mut work_list = Vec::new();
    for i in 0..meshes.len() {
        for j in i..meshes.len() {
            if i != j {
                work_list.push((i, j));
            }
        }
    }

    let collected = work_list
        .par_iter()
        .map(|&(i, j)| {
            let t = meshes[i].transform().inv_mul(&meshes[j].transform());
            let samples = generate_alignment_points(meshes[j].mesh, meshes[i].mesh, &t, params);

            (i, j, samples.len() as f64)
        })
        .collect::<Vec<_>>();

    for (i, j, count) in collected {
        matrix[(i, j)] = count;
        matrix[(j, i)] = count; // Symmetric matrix
    }

    matrix
}
