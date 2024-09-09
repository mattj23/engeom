use crate::common::surface_point::SurfacePoint;
use crate::common::svd_basis::SvdBasis;

pub type Point3 = parry3d_f64::na::Point3<f64>;
pub type Vector3 = parry3d_f64::na::Vector3<f64>;
pub type UnitVec3 = parry3d_f64::na::Unit<Vector3>;
pub type SurfacePoint3 = SurfacePoint<3>;
pub type Iso3 = parry3d_f64::na::Isometry3<f64>;

pub type SvdBasis3 = SvdBasis<3>;

