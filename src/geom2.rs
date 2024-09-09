mod curve2;
pub mod hull;
pub mod kd_tree2;
mod polyline2;
mod angles2;
mod line2;

use crate::common::svd_basis::SvdBasis;
use crate::common::surface_point::SurfacePoint;

pub type Point2 = parry2d_f64::na::Point2<f64>;
pub type Vector2 = parry2d_f64::na::Vector2<f64>;
pub type UnitVec2 = parry2d_f64::na::Unit<Vector2>;
pub type SurfacePoint2 = SurfacePoint<2>;
pub type Iso2 = parry2d_f64::na::Isometry2<f64>;
pub type SvdBasis2 = SvdBasis<2>;
pub type Aabb2 = parry2d_f64::bounding_volume::Aabb;
pub type Ray2 = parry2d_f64::query::Ray;

pub use self::curve2::{Curve2, CurveStation2};
pub use self::angles2::{rot90, rot270, signed_angle, directed_angle};