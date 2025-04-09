//! This module has representations of different types of dimensions

use crate::common::points::mid_point;
use crate::common::surface_point::SurfacePoint;
use parry3d_f64::na::{Point, SVector, Unit};

pub trait Measurement {
    fn value(&self) -> f64;
}

/// Represents a signed distance between two points in space along a specific direction. The value
/// of the measurement will be the vector from `a` to `b` projected onto the direction vector,
/// meaning that the value will be positive if `b` is in the direction of the vector and negative
/// if `b` is in the opposite direction.
pub struct Distance<const D: usize> {
    /// The starting point of the distance measurement
    pub a: Point<f64, D>,

    /// The ending point of the distance measurement
    pub b: Point<f64, D>,

    /// The direction of the distance measurement
    pub direction: Unit<SVector<f64, D>>,
}

impl<const D: usize> Distance<D> {
    /// Create a new signed distance measurement between two points in space along a specific
    /// direction. The value of the measurement will be the vector from `a` to `b` projected onto
    /// the direction vector, meaning that the value will be positive if `b - a` is in the
    /// direction of the vector and negative if not.
    ///
    /// # Arguments
    ///
    /// * `a`: The start point of the distance measurement
    /// * `b`: The end point of the distance measurement
    /// * `direction`: A unit vector representing the direction of the distance measurement. The
    ///   default value is `None`, in which case the direction will be calculated as the normalized
    ///   vector from `a` to `b`, resulting in a positive value with the full magnitude of the
    ///   distance between the two points.
    ///
    /// returns: Distance<{ D }>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Point2, UnitVec2, Vector2};
    /// use engeom::metrology::{Distance2, Measurement};
    /// use approx::assert_relative_eq;
    ///
    /// let a = Point2::new(0.0, 0.0);
    /// let b = Point2::new(1.0, 1.0);
    ///
    /// let d1 = Distance2::new(a, b, None);
    /// let d2 = Distance2::new(b, a, None);
    ///
    /// assert_relative_eq!(d1.value(), 2_f64.sqrt());
    /// assert_relative_eq!(d2.value(), 2_f64.sqrt());
    ///
    /// let dir = UnitVec2::new_unchecked(Vector2::x());
    /// let d3 = Distance2::new(a, b, Some(dir));
    /// let d4 = Distance2::new(b, a, Some(dir));
    ///
    /// assert_relative_eq!(d3.value(), 1.0);
    /// assert_relative_eq!(d4.value(), -1.0);
    /// ```
    pub fn new(
        a: Point<f64, D>,
        b: Point<f64, D>,
        direction: Option<Unit<SVector<f64, D>>>,
    ) -> Self {
        let direction = direction.unwrap_or(Unit::new_normalize(b - a));
        Self { a, b, direction }
    }

    /// Reverse the distance measurement by swapping start and end points and flipping the direction
    /// of measurement. The value and sign of the measurement will be the same after the operation,
    /// but the points will have been swapped.
    pub fn reversed(&self) -> Self {
        Self {
            a: self.b,
            b: self.a,
            direction: -self.direction,
        }
    }

    /// Compute a SurfacePoint that is located halfway between point `a` and `b` and whose normal
    /// vector is the direction of the distance measurement.
    pub fn center(&self) -> SurfacePoint<D> {
        SurfacePoint::new(mid_point(&self.a, &self.b), self.direction)
    }
}

impl<const D: usize> Measurement for Distance<D> {
    fn value(&self) -> f64 {
        self.direction.dot(&(self.b - self.a))
    }
}
