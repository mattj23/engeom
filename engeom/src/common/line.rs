//! This module contains a dimension-generic parameterized line, `Line<D>`, which serves as the
//! shared foundation for the 2D and 3D line types. It is structured similarly to
//! [`SurfacePoint`](crate::common::SurfacePoint): the generic features live here, while the
//! dimension-specific behavior (normals and intersections in 2D, plane/sphere/circle
//! intersections in 3D, axis constructors, and isometry multiplication operators) lives in the
//! `geom2` and `geom3` modules.

use crate::common::PCoords;
use crate::na::{AbstractRotation, Isometry, Point, SVector};
use serde::{Deserialize, Serialize};

/// A parameterized line in D-dimensional space: `P(t) = origin + t * direction`.
///
/// The direction is not required to be normalized; use [`new_normalize`](Line::new_normalize)
/// for a unit-speed parameterization where `t` equals arc length from the origin.
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

    /// Returns a new line with the same origin, but with the direction inverted.
    pub fn new_reversed(&self) -> Self {
        Self::new(self.origin, -self.direction)
    }

    /// Normalizes the direction vector in place so that `t` equals arc length from the origin.
    pub fn normalize(&mut self) {
        self.direction = self.direction.normalize();
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

    /// Moves the origin of the line by a given amount along the direction of the line. A positive
    /// `delta_t` moves the origin forward along the direction of the line, while a negative
    /// `delta_t` moves it backward. The line is modified in place.
    ///
    /// If the direction is not of unit length, keep in mind this shift will be proportional to
    /// the length of the direction vector.
    pub fn shift_origin(&mut self, delta_t: f64) {
        self.origin += self.direction * delta_t;
    }

    /// Returns a new line with the origin shifted by a given amount along the direction of the
    /// line. The direction of the new line is the same as the original line. The original is left
    /// unchanged.
    ///
    /// If the direction is not of unit length, keep in mind this shift will be proportional to
    /// the length of the direction vector.
    pub fn new_shifted_origin(&self, delta_t: f64) -> Self {
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

    /// Returns a new line with both origin and direction transformed by the given isometry.
    pub fn new_transformed_by<R>(&self, iso: &Isometry<f64, R, D>) -> Self
    where
        R: AbstractRotation<f64, D>,
    {
        let mut clone = *self;
        clone.transform_by(iso);
        clone
    }

    /// Transforms this line in place by the given isometry.
    pub fn transform_by<R>(&mut self, iso: &Isometry<f64, R, D>)
    where
        R: AbstractRotation<f64, D>,
    {
        self.origin = iso * self.origin;
        self.direction = iso.rotation.transform_vector(&self.direction);
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
    fn shift_origin_moves_along_direction() {
        let mut line = Line::<2>::new(Point2::new(0.0, 0.0), Vector2::new(2.0, 0.0));
        line.shift_origin(1.5);
        assert_relative_eq!(line.origin, Point2::new(3.0, 0.0), epsilon = 1e-12);
        let shifted = line.new_shifted_origin(-1.0);
        assert_relative_eq!(shifted.origin, Point2::new(1.0, 0.0), epsilon = 1e-12);
    }
}
