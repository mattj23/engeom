pub mod aabb2;
pub mod align2;
mod angles2;
mod arc2;
mod boundary2;
mod circle2;
mod cubic_spline2;
mod curve2;
pub mod hull;
mod line2;
pub mod polyline2;
mod segment2;

use crate::AngleDir;
use crate::AngleDir::Cw;
use crate::common::surface_point::SurfacePoint;
use crate::common::svd_basis::SvdBasis;
use crate::common::{PCoords, SurfacePointCollection};
use crate::na::SVector;
use parry2d_f64::na::UnitComplex;
use serde::{Deserialize, Serialize};
use std::ops;

pub type Point2 = parry2d_f64::na::Point2<f64>;
pub type Vector2 = parry2d_f64::na::Vector2<f64>;
pub type UnitVec2 = parry2d_f64::na::Unit<Vector2>;
pub type SurfacePoint2 = SurfacePoint<2>;
pub type Iso2 = parry2d_f64::na::Isometry2<f64>;
pub type SvdBasis2 = SvdBasis<2>;
pub type Ray2 = parry2d_f64::query::Ray;
pub type Align2 = crate::common::align::Align<UnitComplex<f64>, 2>;
pub type KdTree2 = crate::common::kd_tree::KdTree<2>;

pub use self::aabb2::Aabb2;
pub use self::angles2::{directed_angle, rot90, rot270, signed_angle};
pub use self::arc2::Arc2;
pub use self::boundary2::*;
pub use self::circle2::Circle2;
pub use self::cubic_spline2::CubicSpline2;
pub use self::curve2::{Curve2, CurvePartitioner2, CurveStation2};
pub use self::line2::{Line2, LineOps2, intersect_lines, intersect_rays, intersection_param};
pub use self::segment2::Segment2;

/// This struct represents a position along a 1-manifold embedded in 2D space. Examples of such
/// manifolds are lines, line segments, circles, arcs, splines, etc...anything consisting of a
/// contiguous set of 2D points along a single linear dimension.
///
/// In 2D space, a 1-manifold that is infinitely long (such as a line) or is closed (such as a
/// circle) can divide space into two separate regions, and a single position on a 1-manifold
/// divides the local space in half. As a result, there is an additional feature of the
/// `Manifold1Pos2` to have a surface normal direction, indicating the direction of this
/// separation.  In conformance with the counterclockwise winding convention, the surface normal
/// direction of a manifold is to the right of its direction.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Manifold1Pos2 {
    /// The position of the point along the manifold length
    pub l: f64,

    /// The position of the point in full 2D cartesian space
    pub point: Point2,

    /// The direction of the positive manifold length in full 2D cartesian space at the current
    /// position. You can think of this vector as pointing in the direction of the point at
    /// `self.l + epsilon`
    pub direction: UnitVec2,

    /// The direction of the manifold surface normal in full 2D cartesian space at the current
    /// position. This should follow the counter-clockwise winding order convention, so for most
    /// manifolds this will be equivalent to `self.direction` rotated clockwise by 90 degrees.
    pub normal: UnitVec2,
}

impl Manifold1Pos2 {
    pub fn new(l: f64, point: Point2, direction: UnitVec2, normal: UnitVec2) -> Self {
        Self {
            l,
            point,
            direction,
            normal,
        }
    }

    /// Returns a surface point representing the position and surface normal of the manifold at
    /// this position. The surface point normal will be pointing to the right (clockwise) of the
    /// manifold direction.
    pub fn surface_point(&self) -> SurfacePoint2 {
        SurfacePoint2::new(self.point, self.direction)
    }

    /// Returns a line representing the tangent direction of the manifold at this position. The
    /// line direction has been normalized to unit length.
    pub fn direction_line(&self) -> Line2 {
        Line2::new(self.point, self.direction.into_inner())
    }

    /// Get the scalar projection of a position in 2D space onto the surface normal of the manifold
    /// at this position. A positive value means that the projection is to the outside of the
    /// manifold, a negative value means that it is to the inside. The surface normal is of unit
    /// length, so the scalar projection corresponds to a physical distance.
    ///
    /// # Arguments
    ///
    /// * `other`: A position to take the scalar projection of
    ///
    /// returns: f64
    pub fn normal_scalar_projection(&self, other: &impl PCoords<2>) -> f64 {
        self.normal.dot(&(other.coords() - self.point.coords))
    }

    /// Get the scalar projection of a position in 2D space onto the direction vector of the
    /// manifold at this position. A positive value means that the projection is in front of the
    /// position, a negative value means that it is behind. The direction vector is of unit length
    /// so the scalar projection corresponds to a physical distance.
    ///
    /// # Arguments
    ///
    /// * `other`: A position to take the scalar projection of
    ///
    /// returns: f64
    pub fn direction_scalar_project(&self, other: &impl PCoords<2>) -> f64 {
        self.direction.dot(&(other.coords() - self.point.coords))
    }
}

impl PCoords<2> for Manifold1Pos2 {
    fn coords(&self) -> SVector<f64, 2> {
        self.point.coords
    }
}

impl ops::Mul<SurfacePoint2> for Iso2 {
    type Output = SurfacePoint2;

    fn mul(self, rhs: SurfacePoint2) -> Self::Output {
        rhs.transformed(&self)
    }
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

impl ops::Mul<&SurfacePoint2> for Iso2 {
    type Output = SurfacePoint2;

    fn mul(self, rhs: &SurfacePoint2) -> Self::Output {
        rhs.transformed(&self)
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

impl LineOps2 for SurfacePoint2 {
    fn origin(&self) -> Point2 {
        self.point
    }

    fn dir(&self) -> Vector2 {
        self.normal.into_inner()
    }

    fn at(&self, t: f64) -> Point2 {
        self.at_distance(t)
    }

    fn orthogonal(&self) -> Vector2 {
        rot90(Cw) * self.normal.into_inner()
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
pub mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::distr::{Distribution, Uniform};
    use rand::rngs::StdRng;
    use rand::{Rng, RngExt, SeedableRng};
    use std::f64::consts::PI;

    enum RngSource {
        Seeded(StdRng),
        Thread,
    }

    pub struct RandomGeometry2 {
        rng: RngSource,
    }

    impl Default for RandomGeometry2 {
        fn default() -> Self {
            Self::new()
        }
    }

    impl RandomGeometry2 {
        pub fn new() -> Self {
            Self {
                rng: RngSource::Thread,
            }
        }

        pub fn from_seed(seed: u64) -> Self {
            Self {
                rng: RngSource::Seeded(StdRng::seed_from_u64(seed)),
            }
        }

        pub fn f64(&mut self, lo: f64, hi: f64) -> f64 {
            let u = Uniform::new(lo, hi).unwrap();
            match &mut self.rng {
                RngSource::Seeded(r) => u.sample(r),
                RngSource::Thread => u.sample(&mut rand::rng()),
            }
        }

        pub fn bool(&mut self) -> bool {
            match &mut self.rng {
                RngSource::Seeded(r) => r.random_bool(0.5),
                RngSource::Thread => rand::rng().random_bool(0.5),
            }
        }

        pub fn f64_sym(&mut self, hi: f64) -> f64 {
            self.f64(-hi, hi)
        }

        pub fn point(&mut self, limit: f64) -> Point2 {
            Point2::new(self.f64(-limit, limit), self.f64(-limit, limit))
        }

        pub fn vector2(&mut self, limit: f64) -> Vector2 {
            Vector2::new(self.f64(-limit, limit), self.f64(-limit, limit))
        }

        pub fn unit_vec2(&mut self) -> UnitVec2 {
            UnitVec2::new_normalize(Vector2::new(
                self.angle_sym_pi().cos(),
                self.angle_sym_pi().sin(),
            ))
        }

        pub fn iso2(&mut self, t: f64) -> Iso2 {
            Iso2::new(self.vector2(t), self.angle_sym_pi())
        }

        pub fn angle_sym_pi(&mut self) -> f64 {
            self.f64(-PI, PI)
        }

        pub fn angle_sym_2pi(&mut self) -> f64 {
            self.f64(-2.0 * PI, 2.0 * PI)
        }

        pub fn angle_pos_2pi(&mut self) -> f64 {
            self.f64(0.0, 2.0 * PI)
        }

        pub fn positive(&mut self, limit: f64) -> f64 {
            self.f64(0.0, limit)
        }

        pub fn circle2(&mut self, center_limit: f64, r_min: f64, r_max: f64) -> Circle2 {
            let c = self.point(center_limit);
            let r = self.f64(r_min, r_max);
            Circle2::new(c.x, c.y, r)
        }
    }

    #[test]
    fn iso2_only_rotates_vector() {
        let iso = Iso2::new(Vector2::new(1.0, 2.0), -PI / 2.0);
        let v = Vector2::new(1.0, 1.0);
        let vt = iso * v;
        assert_relative_eq!(vt, Vector2::new(1.0, -1.0), epsilon = 1e-6);
    }
}
