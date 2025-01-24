pub mod kd_tree3;
pub mod mesh;
mod plane3;
mod points;
mod curve3;

use crate::common::surface_point::{SurfacePoint, SurfacePointCollection};
use crate::common::svd_basis::SvdBasis;
use std::ops;

use crate::TransformBy;
pub use mesh::{Mesh, MeshData, UvMapping};
pub use plane3::Plane3;
pub use curve3::{Curve3, CurveStation3};

pub type Point3 = parry3d_f64::na::Point3<f64>;
pub type Vector3 = parry3d_f64::na::Vector3<f64>;
pub type UnitVec3 = parry3d_f64::na::Unit<Vector3>;
pub type SurfacePoint3 = SurfacePoint<3>;
pub type Iso3 = parry3d_f64::na::Isometry3<f64>;

pub type SvdBasis3 = SvdBasis<3>;

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
