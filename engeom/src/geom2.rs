pub mod aabb2;
pub mod align2;
mod angles2;
mod circle2;
mod curve2;
pub mod hull;
mod line2;
pub mod polyline2;

use crate::AngleDir;
use crate::AngleDir::Cw;
use crate::common::SurfacePointCollection;
use crate::common::surface_point::SurfacePoint;
use crate::common::svd_basis::SvdBasis;
use parry2d_f64::na::UnitComplex;
use std::ops;

pub type Point2 = parry2d_f64::na::Point2<f64>;
pub type Vector2 = parry2d_f64::na::Vector2<f64>;
pub type UnitVec2 = parry2d_f64::na::Unit<Vector2>;
pub type SurfacePoint2 = SurfacePoint<2>;
pub type Iso2 = parry2d_f64::na::Isometry2<f64>;
pub type SvdBasis2 = SvdBasis<2>;
pub type Ray2 = parry2d_f64::query::Ray;
pub type Align2 = crate::common::align::Alignment<UnitComplex<f64>, 2>;
pub type KdTree2 = crate::common::kd_tree::KdTree<2>;

pub use self::aabb2::Aabb2;
pub use self::angles2::{directed_angle, rot90, rot270, signed_angle};
pub use self::circle2::{Arc2, Circle2};
pub use self::curve2::{Curve2, CurveStation2};
pub use self::line2::{Line2, Segment2, intersect_rays, intersection_param};

pub trait HasBounds2 {
    fn aabb(&self) -> &Aabb2;
}

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

impl SurfacePoint2 {
    /// Shift a surface point in the direction orthogonal to its normal (sideways) by the given
    /// distance. The direction of travel is the surface point's normal vector rotated by 90 degrees
    /// in the *clockwise* direction. If the normal is pointing up, a positive distance will move
    /// the point to the right, and a negative distance will move the point to the left.
    ///
    /// # Arguments
    ///
    /// * `distance`: the distance to shift the point
    ///
    /// returns: SurfacePoint<2>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::SurfacePoint2;
    ///
    /// let sp = SurfacePoint2::new_normalize([0.0, 0.0].into(), [0.0, 1.0].into());
    /// let shifted = sp.shift_orthogonal(1.0);
    /// assert_relative_eq!(shifted.point, [1.0, 0.0].into(), epsilon = 1e-6);
    /// ```
    pub fn shift_orthogonal(&self, distance: f64) -> Self {
        let shift = rot90(Cw) * self.normal.into_inner() * distance;
        Self::new(self.point + shift, self.normal)
    }

    /// Rotate the normal vector of a surface point by the given angle. The angle is in radians and
    /// is measured counterclockwise from the positive x-axis. If the normal is pointing up, a
    /// positive angle will rotate it to the left, and a negative angle will rotate it to the right.
    ///
    /// # Arguments
    ///
    /// * `angle`: the angle to rotate the normal vector by
    ///
    /// returns: SurfacePoint<2>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{SurfacePoint2, Vector2, UnitVec2};
    ///
    /// let sp = SurfacePoint2::new_normalize([0.0, 0.0].into(), [0.0, 1.0].into());
    /// let rotated = sp.rot_normal(std::f64::consts::PI / 2.0);
    /// let expected = UnitVec2::new_normalize(Vector2::new(-1.0, 0.0));
    /// assert_relative_eq!(rotated.normal, expected, epsilon = 1e-6);
    /// ```
    pub fn rot_normal(&self, angle: f64) -> Self {
        let n = Iso2::rotation(angle) * self.normal.into_inner();
        Self::new_normalize(self.point, n)
    }

    /// Rotate the normal vector of a surface point by 90 degrees in the given direction. If the
    /// normal is pointing up, rotating it clockwise will make it point to the right, and rotating
    /// it counterclockwise will make it point to the left.
    ///
    /// # Arguments
    ///
    /// * `dir`: the direction to rotate the normal vector
    ///
    /// returns: SurfacePoint<2>
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{SurfacePoint2, AngleDir, Vector2, UnitVec2};
    ///
    /// let sp = SurfacePoint2::new_normalize([0.0, 0.0].into(), [0.0, 1.0].into());
    /// let rotated = sp.rot_normal_90(AngleDir::Cw);
    /// let expected = UnitVec2::new_normalize(Vector2::new(1.0, 0.0));
    /// assert_relative_eq!(rotated.normal, expected, epsilon = 1e-6);
    /// ```
    pub fn rot_normal_90(&self, dir: AngleDir) -> Self {
        Self::new_normalize(self.point, rot90(dir) * self.normal.into_inner())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
