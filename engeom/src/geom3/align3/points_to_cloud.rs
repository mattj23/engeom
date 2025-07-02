use super::*;
use crate::geom3::align3::jacobian::{copy_jacobian, point_plane_jacobian, point_point_jacobian};

use crate::Result;
use crate::common::kd_tree::KdTreeSearch;
use crate::common::points::{dist, mean_point};
use crate::geom3::point_cloud::{PointCloudFeatures, PointCloudKdTree};
use crate::geom3::{Align3, Point3, SurfacePoint3, UnitVec3, Vector3};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry3d_f64::na::{Dyn, Matrix, Owned, U1, U6, Vector};

struct PointsToCloud<'a> {
    points: &'a [Point3],
    cloud: &'a PointCloudKdTree<'a>,
    search_radius: f64,
    cloud_normals: &'a [UnitVec3],
    params: RcParams3,
    moved: Vec<Point3>,
    closest: Vec<usize>,
    weight: Vec<f64>,
}

impl<'a> PointsToCloud<'a> {
    fn new(
        points: &'a [Point3],
        cloud: &'a PointCloudKdTree,
        search_radius: f64,
        initial: &Iso3,
    ) -> Self {
        let mean_point = mean_point(points);
        let params = RcParams3::from_initial(initial, &mean_point);

        let cloud_normals = if let Some(normals) = cloud.normals() {
            normals
        } else {
            &[]
        };

        let mut item = Self {
            points,
            cloud,
            search_radius,
            cloud_normals,
            params,
            moved: Vec::with_capacity(points.len()),
            closest: Vec::with_capacity(points.len()),
            weight: Vec::with_capacity(points.len()),
        };

        item.move_points();
        item
    }

    fn move_points(&mut self) {
        let t = self.current_transform();
        self.moved.clear();
        self.closest.clear();
        self.weight.clear();
        let max_dist = self.search_radius * self.search_radius;

        for p in self.points {
            let m = t * *p;
            let (i, d) = self.cloud.tree().nearest_one(&m);
            self.weight.push(distance_weight(d, max_dist));
            self.closest.push(i);
            self.moved.push(m);
        }
    }

    fn current_transform(&self) -> Iso3 {
        *self.params.transform()
    }
}

impl<'a> LeastSquaresProblem<f64, Dyn, U6> for PointsToCloud<'a> {
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
            res[i] = if self.cloud_normals.is_empty() {
                self.weight[i] * dist(p, &self.cloud.points()[*c])
            } else {
                let n = self.cloud_normals[*c];
                let v = p - self.cloud.points()[*c];
                let d = v.dot(&n);
                self.weight[i] * d.abs()
            }
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, ci)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = if self.cloud_normals.is_empty() {
                point_point_jacobian(p, &self.cloud.points()[*ci], &self.params)
            } else {
                let sp = SurfacePoint3::new(self.cloud.points()[*ci], self.cloud_normals[*ci]);
                point_plane_jacobian(p, &sp, &self.params)
            };
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

pub fn points_to_cloud(
    points: &[Point3],
    cloud: &PointCloudKdTree,
    search_radius: f64,
    initial: &Iso3,
) -> Result<Align3> {
    let problem = PointsToCloud::new(points, cloud, search_radius, initial);
    let (result, report) = LevenbergMarquardt::new().minimize(problem);
    if report.termination.was_successful() {
        let residuals = result.residuals().unwrap().as_slice().to_vec();
        Ok(Align3::new(result.current_transform(), residuals))
    } else {
        Err("Failed to align points to cloud".into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::{Mesh, PointCloud};
    use approx::assert_relative_eq;

    fn rotated_box_mesh() -> Mesh {
        let mut box_mesh = Mesh::create_box(10.0, 5.0, 5.0, false);
        let adjust =
            Iso3::rotation(Vector3::new(0.7, 0.7, 0.7)) * Iso3::translation(-5.0, -2.5, -2.5);
        box_mesh.transform_by(&adjust);
        box_mesh
    }

    #[test]
    fn test_points_cloud_box() {
        let box_mesh = rotated_box_mesh();

        let base_cloud = PointCloud::from_surface_points(&box_mesh.sample_poisson(0.1));
        let tree = base_cloud.create_matched_tree().unwrap();
        let ref_cloud = PointCloudKdTree::try_new(&base_cloud, &tree).unwrap();

        let disturb = Iso3::from_parts(
            Translation3::new(-1.0, 0.5, 0.5),
            UnitQuaternion::from_euler_angles(0.1, 0.1, 0.05),
        );

        let points = box_mesh
            .sample_uniform(500)
            .into_iter()
            .map(|p| disturb * p.point)
            .collect::<Vec<_>>();

        let result = points_to_cloud(&points, &ref_cloud, 1.0, &Iso3::identity()).unwrap();
        let test = result.transform();
        let expected = disturb.inverse();

        assert_relative_eq!(expected.translation.x, test.translation.x, epsilon = 5e-2);
        assert_relative_eq!(expected.translation.y, test.translation.y, epsilon = 5e-2);
        assert_relative_eq!(expected.translation.z, test.translation.z, epsilon = 5e-2);

        let test = &test.rotation.euler_angles();
        let expected = &expected.rotation.euler_angles();

        assert_relative_eq!(expected.0, test.0, epsilon = 5e-3);
        assert_relative_eq!(expected.1, test.1, epsilon = 5e-3);
        assert_relative_eq!(expected.2, test.2, epsilon = 5e-3);
    }
}
