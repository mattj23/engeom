pub mod align3;
mod curve3;
mod iso3;
pub mod mesh;
mod plane3;
mod point_cloud;
mod xyzwpr;

use parry3d_f64::na::UnitQuaternion;

use crate::TransformBy;
use crate::common::surface_point::{SurfacePoint, SurfacePointCollection};
use crate::common::svd_basis::SvdBasis;
pub use curve3::{Curve3, CurveStation3};
pub use iso3::IsoExtensions3;
pub use mesh::{Mesh, MeshCollisionSet, UvMapping};
use parry3d_f64::query::Ray;
pub use plane3::Plane3;
pub use point_cloud::{PointCloud, PointCloudFeatures, PointCloudKdTree, PointCloudOverlap};
use std::ops;
pub use xyzwpr::XyzWpr;

pub type Point3 = parry3d_f64::na::Point3<f64>;
pub type Vector3 = parry3d_f64::na::Vector3<f64>;
pub type UnitVec3 = parry3d_f64::na::Unit<Vector3>;
pub type SurfacePoint3 = SurfacePoint<3>;
pub type Iso3 = parry3d_f64::na::Isometry3<f64>;
pub type KdTree3 = crate::common::kd_tree::KdTree<3>;

pub type SvdBasis3 = SvdBasis<3>;
pub type Align3 = crate::common::align::Alignment<UnitQuaternion<f64>, 3>;

pub type Aabb3 = parry3d_f64::bounding_volume::Aabb;

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

impl From<&SurfacePoint3> for Ray {
    fn from(value: &SurfacePoint3) -> Self {
        Ray::new(value.point, value.normal.into_inner())
    }
}
