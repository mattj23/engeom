//! This module contains a dimension-generic parameterized line, `Line<D>`, which serves as the
//! shared foundation for the 2D and 3D line types. It is structured similarly to
//! [`SurfacePoint`](crate::common::SurfacePoint): the generic features live here, while the
//! dimension-specific behavior (normals and intersections in 2D, plane/sphere/circle
//! intersections in 3D, axis constructors, and isometry multiplication operators) lives in the
//! `geom2` and `geom3` modules.

use crate::Result;
use crate::common::PCoords;
use crate::common::svd_basis::SvdBasis;
use crate::na::{AbstractRotation, Isometry, Point, SVector};
use serde::{Deserialize, Serialize};

/// A parameterized line in D-dimensional space: `P(t) = origin + t * direction`.
///
/// The direction is not required to be normalized; use [`new_normalize`](Line::new_normalize)
/// for a unit-speed parameterization where `t` equals arc length from the origin.
///
/// `Line<D>` is the base for `Line2` and `Line3`, which are two of `engeom`'s geometric primitives.
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct Line<const D: usize> {
    pub origin: Point<f64, D>,
    pub direction: SVector<f64, D>,
}

impl<const D: usize> Line<D> {
    /// Create a line from an origin point and a direction vector (stored as-is, not normalized).
    pub fn new(origin: Point<f64, D>, direction: SVector<f64, D>) -> Self {
        Self { origin, direction }
    }

    /// Create a line from an origin point and a direction vector, normalizing the direction so
    /// that the parameter `t` equals arc length from the origin.
    pub fn new_normalize(origin: Point<f64, D>, direction: SVector<f64, D>) -> Self {
        Self::new(origin, direction.normalize())
    }

    /// Create a line through two points. The direction is `p2 - p1` (not normalized).
    pub fn from_points(p1: &impl PCoords<D>, p2: &impl PCoords<D>) -> Self {
        let origin = Point::from(p1.coords());
        Self::new(origin, p2.coords() - p1.coords())
    }

    /// Fit a line to a set of points using singular value decomposition, resulting in a
    /// least-squares fitting. The resulting parameterized line will have its t=0 sitting at the
    /// center of the SVD result.
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of coordinates to fit the line to
    /// * `weights`: if `Some`, this must be a slice of floating points the same length as `points`,
    ///   with the weight value to multiply each point residual by.
    ///
    /// returns: Result<Line<{ D }>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_fit(points: &[impl PCoords<D>], weights: Option<&[f64]>) -> Result<Self> {
        let basis = SvdBasis::from_points(points, weights)
            .ok_or("Failed to fit line with singular value decomposition")?;
        Ok(Line::new(basis.center, basis.largest().into_inner()))
    }

    /// Returns a new line with the same origin, but with the direction inverted.
    pub fn reversed(&self) -> Self {
        Self::new(self.origin, -self.direction)
    }

    /// Returns a new line with the same origin but a normalized direction, so that `t` equals arc
    /// length from the origin.
    pub fn normalized(&self) -> Self {
        Self::new(self.origin, self.direction.normalize())
    }

    /// Returns the point at parameter `t`: `P(t) = origin + t * direction`.
    pub fn at(&self, t: f64) -> Point<f64, D> {
        self.origin + self.direction * t
    }

    /// Returns a new line with the origin shifted by a given amount along the direction of the
    /// line. The direction of the new line is the same as the original line. The original is left
    /// unchanged.
    ///
    /// If the direction is not of unit length, keep in mind this shift will be proportional to
    /// the length of the direction vector.
    pub fn shifted_origin(&self, delta_t: f64) -> Self {
        Self::new(self.origin + self.direction * delta_t, self.direction)
    }

    /// Returns the parameter `t` such that `P(t)` is the closest point on the line to `point`.
    /// Equivalent to the scalar projection of `(point - origin)` onto `direction`, divided by
    /// `|direction|²`.
    pub fn scalar_project(&self, point: &impl PCoords<D>) -> f64 {
        (point.coords() - self.origin.coords).dot(&self.direction) / self.direction.norm_squared()
    }

    /// Returns the closest point on the line to `point`.
    pub fn closest_point(&self, point: &impl PCoords<D>) -> Point<f64, D> {
        self.at(self.scalar_project(point))
    }

    /// Returns the perpendicular distance from `point` to the line.
    pub fn distance_to(&self, point: &impl PCoords<D>) -> f64 {
        let pt = Point::from(point.coords());
        (pt - self.closest_point(&pt)).norm()
    }

    /// Return a new line whose direction is spherically interpolated between this instance and
    /// `other`, and whose origin is _linearly interpolated_ between this instance and `other`.
    ///
    /// A value of `t=0` will return this line, and a value of `t=1` will return `other`.
    ///
    /// # Arguments
    ///
    /// * `other`: the other line to interpolate to
    /// * `t`: the interpolation parameter, does _not_ need to be bounded between [0, 1]
    ///
    /// returns: Line<{ D }>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn slerp(&self, other: &Line<D>, t: f64) -> Self {
        let new_direction = self.direction.slerp(&other.direction, t);
        let shift = other.origin - self.origin;
        Self::new(self.origin + shift * t, new_direction)
    }

    /// Returns a new line with both origin and direction transformed by the given isometry.
    pub fn transformed_by<R>(&self, iso: &Isometry<f64, R, D>) -> Self
    where
        R: AbstractRotation<f64, D>,
    {
        Self::new(
            iso * self.origin,
            iso.rotation.transform_vector(&self.direction),
        )
    }
}

impl<const D: usize> PCoords<D> for Line<D> {
    fn coords(&self) -> SVector<f64, D> {
        self.origin.coords
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Point2, Point3, Vector2, Vector3};
    use approx::assert_relative_eq;

    #[test]
    fn at_endpoints() {
        let line = Line::<2>::new(Point2::new(1.0, 2.0), Vector2::new(0.0, 1.0));
        assert_relative_eq!(line.at(0.0), Point2::new(1.0, 2.0), epsilon = 1e-12);
        assert_relative_eq!(line.at(1.0), Point2::new(1.0, 3.0), epsilon = 1e-12);
    }

    #[test]
    fn new_normalize_gives_unit_direction() {
        let line = Line::<3>::new_normalize(Point3::origin(), Vector3::new(3.0, 0.0, 0.0));
        assert_relative_eq!(line.direction.norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn from_points_direction_is_difference() {
        let line = Line::<3>::from_points(&Point3::new(1.0, 0.0, 0.0), &Point3::new(4.0, 0.0, 0.0));
        assert_relative_eq!(line.direction, Vector3::new(3.0, 0.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn closest_point_perpendicular_drop() {
        let line = Line::<3>::new(Point3::origin(), Vector3::x());
        assert_relative_eq!(
            line.closest_point(&Point3::new(4.0, 3.0, 0.0)),
            Point3::new(4.0, 0.0, 0.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn distance_to_known_value() {
        let line = Line::<3>::new(Point3::origin(), Vector3::x());
        assert_relative_eq!(
            line.distance_to(&Point3::new(0.0, 3.0, 4.0)),
            5.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn shifted_origin_moves_along_direction() {
        let line = Line::<2>::new(Point2::new(0.0, 0.0), Vector2::new(2.0, 0.0));
        let shifted = line.shifted_origin(1.5);
        assert_relative_eq!(shifted.origin, Point2::new(3.0, 0.0), epsilon = 1e-12);
        let shifted_again = shifted.shifted_origin(-1.0);
        assert_relative_eq!(shifted_again.origin, Point2::new(1.0, 0.0), epsilon = 1e-12);
    }
}
