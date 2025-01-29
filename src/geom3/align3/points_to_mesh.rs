use super::*;
use crate::common::DistMode;
use crate::geom3::align3::jacobian::{copy_jacobian, point_plane_jacobian, point_point_jacobian};
use crate::geom3::mesh::Mesh;
use crate::geom3::{Align3, Point3, SurfacePoint3};

use crate::common::points::{dist, mean_point};
use crate::Result;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Dyn, Matrix, Owned, Vector, U1, U6};

/// Attempts to compute the alignment of a set of points to the surface of a mesh using a
/// Levenberg-Marquardt solver.  The points are projected onto their closest matching surface point
/// on the mesh, and the residuals being minimized are the distance between the projected point and
/// the surface point.
///
/// The `mode` parameter determines whether the residuals being minimized are the entire euclidean
/// distance between the points and their closest corresponding point on the surface of the mesh, or
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
/// * `points`: a slice containing the points to be aligned
/// * `mesh`: the reference mesh entity onto which the points are being projected
/// * `initial`: an initial guess for the alignment transform
/// * `mode`: using `DistMode::ToPoint` will minimize the distance between the points and their
/// closest corresponding point on the mesh, while using `DistMode::ToPlane` will do the same
/// except that it will ignore the component of the distance orthogonal to the surface normal.
///
/// returns: Result<Alignment<Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use parry3d_f64::na::{Translation3, UnitQuaternion};
/// use engeom::common::DistMode;
/// use engeom::geom3::{Align3, Point3, Mesh, Iso3, Vector3};
/// use engeom::geom3::align3::points_to_mesh;
/// let mesh = Mesh::create_box(10.0, 5.0, 2.0, false);
///
/// let disturb = Iso3::from_parts(
///     Translation3::new(1.0, 2.0, 3.0),
///     UnitQuaternion::from_euler_angles(0.1, 0.2, 0.3)
/// );
///
/// let points = mesh.sample_uniform(1000)
///     .into_iter()
///     .map(|p| disturb * p.point)
///     .collect::<Vec<_>>();
///
/// let result = points_to_mesh(&points, &mesh, &Iso3::identity(), DistMode::ToPlane).unwrap();
/// let expected = disturb.inverse();
/// assert_relative_eq!(result.transform().to_matrix(), expected.to_matrix(), epsilon = 1e-6);
/// ```
pub fn points_to_mesh(
    points: &[Point3],
    mesh: &Mesh,
    initial: &Iso3,
    mode: DistMode,
) -> Result<Align3> {
    let problem = PointsToMesh::new(points, mesh, initial, mode);
    let (result, report) = LevenbergMarquardt::new().minimize(problem);
    if report.termination.was_successful() {
        let residuals = result.residuals().unwrap().as_slice().to_vec();
        Ok(Align3::new(result.current_transform(), residuals))
    } else {
        Err("Failed to align points to mesh".into())
    }
}

struct PointsToMesh<'a> {
    points: &'a [Point3],
    mesh: &'a Mesh,
    params: RcParams3,
    moved: Vec<Point3>,
    closest: Vec<SurfacePoint3>,
    mode: DistMode,
}

impl<'a> PointsToMesh<'a> {
    fn new(points: &'a [Point3], mesh: &'a Mesh, initial: &Iso3, mode: DistMode) -> Self {
        let mean_point = mean_point(points);
        let params = RcParams3::from_initial(initial, &mean_point);

        let mut item = Self {
            points,
            mesh,
            params,
            moved: Vec::with_capacity(points.len()),
            closest: Vec::with_capacity(points.len()),
            mode,
        };

        item.move_points();
        item
    }

    fn move_points(&mut self) {
        let t = self.current_transform();
        self.moved.clear();
        self.closest.clear();
        for p in self.points {
            let m = t * *p;
            self.closest.push(self.mesh.surf_closest_to(&m));
            self.moved.push(m);
        }
    }

    fn current_transform(&self) -> Iso3 {
        *self.params.transform()
    }
}

impl<'a> LeastSquaresProblem<f64, Dyn, U6> for PointsToMesh<'a> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.params.set(x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.params.x
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            res[i] = match self.mode {
                DistMode::ToPoint => dist(p, &c.point),
                DistMode::ToPlane => c.scalar_projection(p).abs(),
            };
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let _center = self.params.transform() * self.params.rc;
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = match self.mode {
                DistMode::ToPoint => point_point_jacobian(p, &c.point, &self.params),
                DistMode::ToPlane => point_plane_jacobian(p, c, &self.params),
            };
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// This tests whether the initial transform is correctly re-expressed in the problem's frame
    /// as rotations around the mean point.
    #[test]
    fn test_initial_round_trip() {
        let box_mesh = Mesh::create_box(10.0, 5.0, 2.0, false);
        let points = box_mesh
            .sample_uniform(1000)
            .into_iter()
            .map(|p| p.point)
            .collect::<Vec<_>>();

        let initial = iso3_from_param(&T3Storage::new(1.0, 2.0, 3.0, 0.1, 0.2, 0.3));
        let problem = PointsToMesh::new(&points, &box_mesh, &initial, DistMode::ToPlane);
        let result = problem.current_transform();
        assert_relative_eq!(result.to_matrix(), initial.to_matrix(), epsilon = 1e-8);
    }
}
