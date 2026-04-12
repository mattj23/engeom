use super::*;
use crate::common::DistMode;
use crate::geom3::align3::jacobian::{
    copy_jacobian, point_plane_jacobian_full, point_point_jacobian_full, point_surf_jacobian,
};
use crate::geom3::mesh::Mesh;
use crate::geom3::{Align3, Point3, SurfacePoint3};

use crate::Result;
use crate::common::kd_tree::{KdTree, KdTreeSearch};
use crate::common::points::{dist, mean_point};
use crate::geom3::align3::params::AlignParams3;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use num_traits::real::Real;
use parry3d_f64::na::{Dyn, Matrix, Owned, U1, U6, Vector};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};
use rayon::prelude::*;

/// Performs a limited form of RANSAC on a set of points to find a rough alignment to the mesh.
///
/// You provide a number of sample points. It will work best if they are roughly evenly spaced, so
/// consider generating or filtering them with Poisson disk sampling. The RANSAC algorithm will
/// randomly pick three of these points, then find the three nearest neighbors of each, for a total
/// of up to nine points (duplicates are removed) in the subsample. The `points_to_mesh` algorithm
/// is then performed on these subsample points, and if the alignment converges, it iterates
/// through the entire set of original points, transforming them by the result, and counting the
/// number of inliers.  Inliers are points that are less than `inlier_threshold` from the target
/// mesh after transformation.
///
/// The transformation with the highest number of inliers is returned.
///
/// # Arguments
///
/// * `points`: the points to be aligned, in their own local coordinate system
/// * `mesh`: the target mesh entity which the points will be aligned to
/// * `working_iso`: a transform from the points' local coordinate system to a working space.  The
///   points will be aligned as if they were starting at the working position, and the result will
///   be such that `result.transform() * working_iso * points[...]` brings the raw points to the
///   target
/// * `center`: the rotation center for the alignment, if specified, otherwise the mean point of
///   the points will be used
/// * `dof`: a 6-DOF constraint to be applied to the alignment
/// * `inlier_threshold`: the maximum distance between a point and the target mesh that is still
///   considered an inlier
/// * `iterations`: the number of iterations of RANSAC to perform
///
/// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
pub fn ransac_points_to_mesh(
    points: &[Point3],
    mesh: &Mesh,
    working_iso: Iso3,
    mode: DistMode,
    center: Option<Point3>,
    dof: Dof6,
    inlier_threshold: f64,
    iterations: usize,
) -> Result<Iso3> {
    // let tree = KdTree::try_new(points)?;
    // let results = ransac_indices(iterations, points.len())
    //     .par_iter()
    //     .map(|(a, b, c)| {
    //         let mut local_indices = vec![*a, *b, *c];
    //         for (i, _) in tree.nearest(&points[*a], 3) {
    //             local_indices.push(i);
    //         }
    //         for (i, _) in tree.nearest(&points[*b], 3) {
    //             local_indices.push(i);
    //         }
    //         for (i, _) in tree.nearest(&points[*c], 3) {
    //             local_indices.push(i);
    //         }
    //         // Sort and deduplicate
    //         local_indices.sort();
    //         local_indices.dedup();
    //
    //         let local_points = local_indices.iter().map(|&i| points[i]).collect::<Vec<_>>();
    //
    //         if let Ok(align) = points_to_mesh(&local_points, mesh, working_iso, mode, center, dof) {
    //             let t = align.transform() * working_iso;
    //             let mut inliers = 0;
    //             for test_point in points.iter() {
    //                 if mesh.distance_closest_to(&(t * test_point)) < inlier_threshold {
    //                     inliers += 1;
    //                 }
    //             }
    //
    //             Some((*align.transform(), inliers))
    //         } else {
    //             None
    //         }
    //     })
    //     .filter_map(|x| x)
    //     .collect::<Vec<_>>();
    // let (best_transform, _) = results.iter().max_by_key(|(_, inliers)| *inliers).unwrap();
    // Ok(*best_transform)

    todo!()
}

fn ransac_indices(iterations: usize, point_len: usize) -> Vec<(usize, usize, usize)> {
    let mut rng = StdRng::seed_from_u64(0);
    (0..iterations)
        .map(|_| {
            let a = rng.random_range(0..point_len);
            let mut b = rng.random_range(0..point_len);
            while b == a {
                b = rng.random_range(0..point_len);
            }

            let mut c = rng.random_range(0..point_len);
            while c == a || c == b {
                c = rng.random_range(0..point_len);
            }

            (a, b, c)
        })
        .collect()
}

pub fn points_to_mesh(points: &[Point3], mesh: &Mesh, params: AlignParams3) -> Result<Align3> {
    let problem = PointsToMesh::new(points, mesh, params);
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        let residuals = result.residuals().unwrap().as_slice().to_vec();
        Ok(Align3::new(result.params.final_result(), residuals))
    } else {
        Err("Failed to align points to mesh".into())
    }
}

struct PointsToMesh<'a> {
    points: &'a [Point3],
    mesh: &'a Mesh,
    params: AlignParams3,
    moved: Vec<Point3>,
    closest: Vec<SurfacePoint3>,
    residuals: Vec<f64>,
    weights: Vec<f64>,
    parallel: bool,
}

impl<'a> PointsToMesh<'a> {
    fn new(points: &'a [Point3], mesh: &'a Mesh, params: AlignParams3) -> Self {
        let mut x = Self {
            points,
            mesh,
            params,
            moved: vec![Point3::default(); points.len()],
            closest: vec![SurfacePoint3::default(); points.len()],
            residuals: vec![0.0; points.len()],
            weights: vec![1.0; points.len()],
            parallel: false,
        };

        x.move_points();
        x
    }

    /// Internally, this moves the points and computes the closest surface point on the mesh to
    /// each.
    fn move_points(&mut self) {
        let current = self.params.current_values();
        let indices = (0..self.points.len()).collect::<Vec<_>>();

        if self.parallel {
            let collected = indices
                .par_iter()
                .map(|&i| {
                    let m = current.transform * self.points[i];
                    let c = self.mesh.surf_closest_to(&m);
                    (i, m, c)
                })
                .collect::<Vec<_>>();
            for (i, m, c) in collected {
                self.moved[i] = m;
                self.closest[i] = c.sp;
            }
        } else {
            for (i, &j) in indices.iter().enumerate() {
                let m = current.transform * self.points[j];
                let c = self.mesh.surf_closest_to(&m);
                self.moved[i] = m;
                self.closest[i] = c.sp;
            }
        }

        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            // The residual is the distance between the test point and the closest point on the
            // mesh surface, adjusted for the direction of the scalar projection.
            self.residuals[i] = dist(p, &c.point) * c.scalar_projection(p).signum();
        }

        // Weights and thresholding are disabled for now, but this is where they would go
        // let (mean, stdev) = mean_and_stdev(&self.residuals).unwrap();
        // let limit = (stdev * 3.0).max(1e-6);
        // self.weights = (0..self.points.len())
        //     .map(|i| {
        //         if (self.residuals[i] - mean).abs() > limit {
        //             0.0
        //         } else {
        //             1.0
        //         }
        //     })
        //     .collect::<Vec<_>>();
    }
}

impl LeastSquaresProblem<f64, Dyn, U6> for PointsToMesh<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.params.set_storage(*x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.params.get_storage()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for i in 0..self.points.len() {
            res[i] = self.residuals[i] * self.weights[i];
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let current = self.params.current_values();
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = point_surf_jacobian(p, c, &current) * self.weights[i];
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::{clone_points, transform_points};
    use crate::tests::engine_blade;
    use crate::{SelectOp, Selection};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn simple_box_disturbed() -> Result<()> {
        // This test is to verify that a simple test against a box that doesn't have large rotations
        // produces a result that is roughly the inverse of the disturbance
        let mesh = Mesh::create_box(10.0, 5.0, 2.0, false);
        let points = clone_points(&mesh.sample_poisson(0.1, None));
        let disturb = Iso3::from_parts(
            Translation3::new(3.0, 2.0, 1.0),
            UnitQuaternion::from_euler_angles(PI / 8.0, PI / 12.0, PI / 16.0),
        );

        let params = AlignParams3::new_at_origin(None);
        let to_align = transform_points(&points, &disturb);
        let result = points_to_mesh(&to_align, &mesh, params)?;

        assert_relative_eq!(disturb.inverse(), result.transform(), epsilon = 1e-8);
        Ok(())
    }

    #[test]
    fn blade_example() -> Result<()> {
        let mesh = engine_blade();
        let mask = mesh
            .face_select(Selection::None)
            .facing(&Vector3::y(), PI / 4.0, SelectOp::Add)
            .take_mask();
        let expected_points = clone_points(&mesh.sample_poisson(2.0, Some(&mask)));

        // TODO: figure out why this thing lost robustness
        let disturb = Iso3::from_parts(
            Translation3::new(-100.0, 150.0, 0.0),
            UnitQuaternion::new(Vector3::new(1.0, 1.0, 1.0).normalize() * PI / 6.0),
        );

        let to_align = transform_points(&expected_points, &disturb);

        let params = AlignParams3::new_at_center(mean_point(&to_align), None);
        let result = points_to_mesh(&to_align, &mesh, params)?;

        let aligned = transform_points(&to_align, result.transform());

        let max_deviation = aligned
            .iter()
            .zip(expected_points.iter())
            .map(|(a, e)| (a - e).norm())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        assert!(
            max_deviation < 1e-6,
            "Max deviation is too high: {}",
            max_deviation
        );

        Ok(())
    }
}
