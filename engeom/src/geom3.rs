pub mod align3;
pub mod attributes3;
mod circle3;
mod cone3;
mod cubic_spline3;
mod curve3;
mod curve_group3;
mod cylinder3;
pub mod half_edge3;
mod iso3;
mod line3;
mod manifold;
pub mod mesh;
mod planar_map;
mod plane3;
pub mod point_cloud;
mod segment3;
mod sphere3;
mod swept_boundary3;
mod unroll_transform;
mod xyz_quat;
mod xyzwpr;

use parry3d_f64::na::UnitQuaternion;

use crate::common::surface_point::{SurfacePoint, SurfacePointCollection};
use crate::common::svd_basis::SvdBasis;
pub use attributes3::{Attr3, FaceAttrSet3, PointAttrSet3};
pub use circle3::Circle3;
pub use cone3::Cone3;
pub use cubic_spline3::CubicSpline3;
pub use curve_group3::CurveGroup3;
pub use curve3::{Curve3, CurveStation3};
pub use cylinder3::Cylinder3;
pub use half_edge3::HalfEdgeMesh3;
pub use iso3::IsoExtensions3;
pub use line3::Line3;
pub use manifold::Manifold1Pos3;
pub use mesh::{FlatDomain, Mesh3, MeshCollisionSet, MeshData3, MeshView3, PlanarSection};
use parry3d_f64::query::Ray;
pub use planar_map::{PlanarMap, PlaneFrame};
pub use plane3::Plane3;
pub use point_cloud::{
    CloudIndex3, PointCloud3, PointCloudOverlap, VOXEL_COHERENCE_ATTR, VOXEL_COUNT_ATTR,
};
pub use segment3::Segment3;
pub use sphere3::Sphere3;
use std::ops;
pub use swept_boundary3::*;
pub use unroll_transform::UnrollTransform;
pub use xyz_quat::XyzQuat;
pub use xyzwpr::XyzWpr;

pub type Point3 = parry3d_f64::na::Point3<f64>;
pub type Vector3 = parry3d_f64::na::Vector3<f64>;
pub type UnitVec3 = parry3d_f64::na::Unit<Vector3>;
pub type Iso3 = parry3d_f64::na::Isometry3<f64>;
pub type KdTree3 = crate::common::kd_tree::KdTree<3>;

/// A point in 3D space paired with a unit normal direction, representing a position on a
/// 2-manifold (surface, mesh face, etc...) embedded in 3D space. Mathematically this is identical
/// to a parameterized line or ray, but represents a different conceptual entity and so has
/// different semantics.
///
/// This is one of `engeom`'s 3D geometric primitives.
pub type SurfacePoint3 = SurfacePoint<3>;

pub type SvdBasis3 = SvdBasis<3>;
pub type Align3 = crate::common::align::Alignment<UnitQuaternion<f64>, 3>;
pub type AlignOutcome3 = crate::common::align::AlignOutcome<UnitQuaternion<f64>, 3>;
pub type MultiOutcome3 = crate::common::align::MultiOutcome<UnitQuaternion<f64>, 3>;

pub type Aabb3 = parry3d_f64::bounding_volume::Aabb;

impl ops::Mul<SurfacePoint3> for &Iso3 {
    type Output = SurfacePoint3;

    fn mul(self, rhs: SurfacePoint3) -> Self::Output {
        rhs.transformed_by(self)
    }
}

impl ops::Mul<&SurfacePoint3> for &Iso3 {
    type Output = SurfacePoint3;

    fn mul(self, rhs: &SurfacePoint3) -> Self::Output {
        rhs.transformed_by(self)
    }
}

impl SurfacePointCollection<3> for &Vec<SurfacePoint3> {
    fn clone_points(&self) -> Vec<Point3> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec3> {
        self.iter().map(|sp| sp.normal).collect()
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

impl SurfacePointCollection<3> for &[SurfacePoint3] {
    fn clone_points(&self) -> Vec<Point3> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec3> {
        self.iter().map(|sp| sp.normal).collect()
    }
}

impl From<&SurfacePoint3> for Ray {
    fn from(value: &SurfacePoint3) -> Self {
        Ray::new(value.point, value.normal.into_inner())
    }
}

impl Default for SurfacePoint3 {
    fn default() -> Self {
        SurfacePoint3::new(Point3::origin(), Vector3::x_axis())
    }
}
