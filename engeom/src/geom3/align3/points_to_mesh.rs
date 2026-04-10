use super::*;
use crate::common::DistMode;
use crate::geom3::align3::jacobian::{
    copy_jacobian, point_plane_jacobian, point_plane_jacobian_full, point_point_jacobian,
    point_point_jacobian_full,
};
use crate::geom3::mesh::Mesh;
use crate::geom3::{Align3, Point3, SurfacePoint3};

use crate::Result;
use crate::common::points::{dist, mean_point};
use crate::geom3::align3::params::AlignParams3;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Dyn, Matrix, Owned, U1, U6, Vector};
use rayon::prelude::*;
use crate::common::vec_f64::mean_and_stdev;

/// Attempts to compute the alignment of a set of points to the surface of a mesh using a
/// Levenberg-Marquardt solver.  The points are projected onto their closest matching surface point
/// on the mesh, and the residuals being minimized are the distance between the projected point and
/// the surface point.
///
/// The `mode` parameter determines whether the residuals being minimized are the entire Euclidean
/// distance between the points and their closest corresponding point on the surface of the mesh or
/// just the component of that distance orthogonal to the surface normal:
///
/// `DistMode::ToPoint` is often useful when you know that the points will match well with the
/// surface, or if you are less sure about how close the initial guess is to the correct alignment.
///
/// `DistMode::ToPlane` is typically faster and will not penalize points which slip off the edge of
/// a triangle at the boundary of the mesh so long as they are still close to the plane of the
/// triangle.
///
/// # Arguments
///
/// * `points`: the points to be aligned, in their own local coordinate system
/// * `mesh`: the target mesh entity which the points will be aligned to
/// * `working_iso`: a transform from the points' local coordinate system to a working space.  The
///   points will be aligned as if they were starting at the working position, and the result will
///   be such that `result.transform() * working_iso * points[...]` brings the raw points to the
///   target
/// * `mode`: either `DistMode::ToPoint` or `DistMode::ToPlane`
/// * `center`: the rotation center for the alignment, if specified, otherwise the mean point of
///   the points will be used
/// * `dof`: a 6-DOF constraint to be applied to the alignment
///
/// returns: Result<Alignment<Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
pub fn points_to_mesh(
    points: &[Point3],
    mesh: &Mesh,
    working_iso: Iso3,
    mode: DistMode,
    center: Option<Point3>,
    dof: Dof6,
) -> Result<Align3> {
    let problem = if let Some(c) = center {
        PointsToMesh::new_with_center(points, mesh, working_iso, mode, c, dof)
    } else {
        PointsToMesh::new_from_mean(points, mesh, working_iso, mode, dof)
    };

    let (result, report) = LevenbergMarquardt::new().minimize(problem);
    if report.termination.was_successful() {
        let residuals = result.residuals().unwrap().as_slice().to_vec();
        Ok(Align3::new(result.params.final_result(), residuals))
    } else {
        Err("Failed to align points to mesh".into())
    }
}

/// This is the internal implementation of a point-to-mesh Levenberg-Marquardt alignment problem.
/// It handles a transform to a local working space, an arbitrary center of rotation, the ability
/// to apply a general 6-DOF constraint, and two different modes of computing distance.
///
/// This is the simplest implementation of a general R^3 alignment problem and is a reference for
/// more complicated alignments.
struct PointsToMesh<'a> {
    points: &'a [Point3],
    mesh: &'a Mesh,
    params: AlignParams3,
    moved: Vec<Point3>,
    closest: Vec<SurfacePoint3>,
    mode: DistMode,

    residuals: Vec<f64>,
    weights: Vec<f64>,
}

impl<'a> PointsToMesh<'a> {
    /// Create a new alignment problem with a given center of rotation and a 6-DOF constraint. The
    /// `working_iso` argument is used to transform the points into a working space, and the
    /// result of the alignment will be a transform from `working_iso` to the target.
    ///
    /// # Arguments
    ///
    /// * `points`: the points to be aligned, in their original coordinate system
    /// * `mesh`: the target mesh entity which the points will be aligned to
    /// * `working_iso`: a transform from the points' local coordinate system to a working space.
    ///   The points will be aligned as if they were starting at the working position.
    /// * `mode`: using `DistMode::ToPoint` will minimize the distance between the points and their
    ///   closest corresponding point on the mesh, while using `DistMode::ToPlane` will do the same
    ///   except that it will ignore the component of the distance orthogonal to the surface normal.
    /// * `center`: the center of rotation for the alignment
    /// * `constraint`: a 6-DOF constraint to be applied to the alignment
    ///
    /// returns: PointsToMesh
    fn new_with_center(
        points: &'a [Point3],
        mesh: &'a Mesh,
        working_iso: Iso3,
        mode: DistMode,
        center: Point3,
        constraint: Dof6,
    ) -> Self {
        let params = AlignParams3::new(center, working_iso, constraint);
        let count = points.len();

        let mut item = Self {
            points,
            mesh,
            params,
            moved: vec![Point3::origin(); count],
            closest: vec![Default::default(); count],
            mode,
            residuals: vec![0.0; count],
            weights: vec![1.0; count],
        };

        item.move_points();
        item
    }

    /// Create a new alignment problem with a given center of rotation and a 6-DOF constraint. The
    /// `working_iso` argument is used to transform the points into a working space, and the
    /// result of the alignment will be a transform from `working_iso` to the target.  The center of
    /// rotation is computed from the mean position of the points.
    ///
    /// # Arguments
    ///
    /// * `points`: slice containing the points to be aligned
    /// * `mesh`: the target mesh entity which the points will be aligned to
    /// * `working_iso`: a transform from the points' local coordinate system to a working space.
    ///   The points will be aligned as if they were starting at the working position.
    /// * `mode`: using `DistMode::ToPoint` will minimize the distance between the points and their
    ///   closest corresponding point on the mesh, while using `DistMode::ToPlane` will do the same
    ///   except that it will ignore the component of the distance orthogonal to the surface normal.
    /// * `constraint`: a 6-DOF constraint to be applied to the alignment
    ///
    /// returns: PointsToMesh
    fn new_from_mean(
        points: &'a [Point3],
        mesh: &'a Mesh,
        working_iso: Iso3,
        mode: DistMode,
        constraint: Dof6,
    ) -> Self {
        let mean_point = mean_point(points);
        Self::new_with_center(points, mesh, working_iso, mode, mean_point, constraint)
    }

    /// Internally, this moves
    fn move_points(&mut self) {
        let current = self.params.current_values();
        let indices = (0..self.points.len()).collect::<Vec<_>>();
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

        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            self.residuals[i] = match self.mode {
                DistMode::ToPoint => dist(p, &c.point),
                DistMode::ToPlane => c.scalar_projection(p).abs(),
            };
        }

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
            let values = match self.mode {
                DistMode::ToPoint => point_point_jacobian_full(p, &c.point, &current),
                DistMode::ToPlane => point_plane_jacobian_full(p, c, &current),
            } * self.weights[i];
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

    // /// This tests whether the initial transform is correctly re-expressed in the problem's frame
    // /// as rotations around the mean point.
    // #[test]
    // fn initial_round_trip() {
    //     let box_mesh = Mesh::create_box(10.0, 5.0, 2.0, false);
    //     let points = box_mesh
    //         .sample_uniform(1000)
    //         .into_iter()
    //         .map(|p| p.point)
    //         .collect::<Vec<_>>();
    //
    //     let initial = iso3_from_param(&T3Storage::new(1.0, 2.0, 3.0, 0.1, 0.2, 0.3));
    //     let problem = PointsToMesh::new(&points, &box_mesh, &initial, DistMode::ToPlane, None);
    //     let result = problem.current_transform();
    //     assert_relative_eq!(result.to_matrix(), initial.to_matrix(), epsilon = 1e-8);
    // }

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

        let to_align = transform_points(&points, &disturb);
        let result = points_to_mesh(
            &to_align,
            &mesh,
            Iso3::identity(),
            DistMode::ToPoint,
            None,
            Default::default(),
        )?;

        assert_relative_eq!(disturb.inverse(), result.transform(), epsilon = 1e-8);
        Ok(())
    }

    #[test]
    fn simple_box_disturbed_working() -> Result<()> {
        // This test is to verify that when using a working transform, the overall disturbance
        // transform is the inverse of `result * working`

        let mesh = Mesh::create_box(10.0, 5.0, 2.0, false);
        let points = clone_points(&mesh.sample_poisson(0.1, None));
        let disturb = Iso3::from_parts(
            Translation3::new(3.0, 2.0, 1.0),
            UnitQuaternion::from_euler_angles(PI / 8.0, PI / 12.0, PI / 16.0),
        );
        let to_align = transform_points(&points, &disturb);

        let working = Iso3::from_parts(
            Translation3::new(0.5, 0.25, 0.125),
            UnitQuaternion::from_euler_angles(PI / 16.0, PI / 24.0, PI / 32.0),
        );
        let result = points_to_mesh(
            &to_align,
            &mesh,
            working,
            DistMode::ToPoint,
            None,
            Dof6::all(),
        )?;

        assert_relative_eq!(
            disturb.inverse(),
            result.transform() * working,
            epsilon = 1e-8
        );
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
            // Old version that worked
            // Translation3::new(-100.0, 150.0, 0.0),
            // UnitQuaternion::new(Vector3::new(1.0, 1.0, 1.0).normalize() * PI / 6.0),
            Translation3::new(10.0, 10.0, 10.0),
            UnitQuaternion::new(Vector3::new(1.0, 1.0, 1.0).normalize() * -PI / 16.0),
        );

        let to_align = transform_points(&expected_points, &disturb);

        let result = points_to_mesh(
            &to_align,
            &mesh,
            Iso3::identity(),
            DistMode::ToPoint,
            None,
            Dof6::all(),
        )?;

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
