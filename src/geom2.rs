pub mod align2;
mod angles2;
mod circle2;
mod curve2;
pub mod hull;
pub mod kd_tree2;
mod line2;
pub mod polyline2;

use crate::common::surface_point::SurfacePoint;
use crate::common::svd_basis::SvdBasis;
use crate::common::SurfacePointCollection;
use parry2d_f64::na::UnitComplex;
use std::ops;

pub type Point2 = parry2d_f64::na::Point2<f64>;
pub type Vector2 = parry2d_f64::na::Vector2<f64>;
pub type UnitVec2 = parry2d_f64::na::Unit<Vector2>;
pub type SurfacePoint2 = SurfacePoint<2>;
pub type Iso2 = parry2d_f64::na::Isometry2<f64>;
pub type SvdBasis2 = SvdBasis<2>;
pub type Aabb2 = parry2d_f64::bounding_volume::Aabb;
pub type Ray2 = parry2d_f64::query::Ray;
pub type Align2 = crate::common::align::Alignment<UnitComplex<f64>, 2>;

pub use self::angles2::{directed_angle, rot270, rot90, signed_angle};
pub use self::circle2::{Arc2, Circle2};
pub use self::curve2::{Curve2, CurveStation2};
pub use self::line2::{intersect_rays, intersection_param, Line2, Segment2};

impl ops::Mul<SurfacePoint2> for &Iso2 {
    type Output = SurfacePoint2;

    fn mul(self, rhs: SurfacePoint2) -> Self::Output {
        rhs.transformed(self)
    }
}

impl ops::Mul<&SurfacePoint2> for &Iso2 {
    type Output = SurfacePoint2;

    fn mul(self, rhs: &SurfacePoint2) -> Self::Output {
        rhs.transformed(self)
    }
}

impl SurfacePointCollection<2> for &[SurfacePoint2] {
    fn clone_points(&self) -> Vec<Point2> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec2> {
        self.iter().map(|sp| sp.normal).collect()
    }
}

impl SurfacePointCollection<2> for &Vec<SurfacePoint2> {
    fn clone_points(&self) -> Vec<Point2> {
        self.iter().map(|sp| sp.point).collect()
    }

    fn clone_normals(&self) -> Vec<UnitVec2> {
        self.iter().map(|sp| sp.normal).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::linear_space;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn iso2_only_rotates_vector() {
        let iso = Iso2::new(Vector2::new(1.0, 2.0), -PI / 2.0);
        let v = Vector2::new(1.0, 1.0);
        let vt = iso * v;
        assert_relative_eq!(vt, Vector2::new(1.0, -1.0), epsilon = 1e-6);
    }
}
