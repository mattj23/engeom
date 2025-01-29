pub mod align3;
mod curve3;
pub mod kd_tree3;
pub mod mesh;
mod plane3;
mod points;

use crate::common::surface_point::{SurfacePoint, SurfacePointCollection};
use crate::common::svd_basis::SvdBasis;
use std::ops;

use crate::{Result, TransformBy};
pub use curve3::{Curve3, CurveStation3};
pub use mesh::{Mesh, MeshData, UvMapping};
use parry3d_f64::na::{try_convert, Matrix4, UnitQuaternion};
pub use plane3::Plane3;

pub type Point3 = parry3d_f64::na::Point3<f64>;
pub type Vector3 = parry3d_f64::na::Vector3<f64>;
pub type UnitVec3 = parry3d_f64::na::Unit<Vector3>;
pub type SurfacePoint3 = SurfacePoint<3>;
pub type Iso3 = parry3d_f64::na::Isometry3<f64>;

pub type SvdBasis3 = SvdBasis<3>;
pub type Align3 = crate::common::align::Alignment<UnitQuaternion<f64>, 3>;

impl ops::Mul<SurfacePoint3> for &Iso3 {
    type Output = SurfacePoint3;

    fn mul(self, rhs: SurfacePoint3) -> Self::Output {
        rhs.transformed(self)
    }
}

impl ops::Mul<&SurfacePoint3> for &Iso3 {
    type Output = SurfacePoint3;

    fn mul(self, rhs: &SurfacePoint3) -> Self::Output {
        rhs.transformed(self)
    }
}

impl SurfacePointCollection<3> for Vec<SurfacePoint3> {
    fn clone_points(&self) -> Vec<Point3> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec3> {
        self.iter().map(|sp| sp.normal).collect()
    }
}

impl TransformBy<Iso3, Vec<Point3>> for &[Point3] {
    fn transform_by(&self, transform: &Iso3) -> Vec<Point3> {
        self.iter().map(|p| transform * p).collect()
    }
}

impl TransformBy<Iso3, Vec<Point3>> for &Vec<Point3> {
    fn transform_by(&self, transform: &Iso3) -> Vec<Point3> {
        self.iter().map(|p| transform * p).collect()
    }
}

/// Try to convert a 16-element array representing a 4x4 matrix to an Iso3. The array is assumed to
/// be in row-major order, such that the first four elements are the first row of the matrix, the
/// second four elements are the second row, and so on.
///
/// # Arguments
///
/// * `array`: &[f64; 16] - The array to convert
///
/// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
pub fn iso3_try_from_array(array: &[f64; 16]) -> Result<Iso3> {
    try_convert(Matrix4::from_row_slice(array)).ok_or("Could not convert to Iso3".into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn iso3_try_from_array_simple() {
        let array = [
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 1.0,
        ];
        let iso = iso3_try_from_array(&array).unwrap();
        let m = iso.to_matrix();
        let expected = Matrix4::new(
            1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 1.0,
        );
        assert_relative_eq!(m, expected);
    }
}
