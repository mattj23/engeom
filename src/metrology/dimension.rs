//! This module has representations of different types of dimensions

use crate::{Point2, SurfacePoint2, UnitVec2};

pub trait Dimension {
    fn value(&self) -> f64;
}

/// Represents a length measurement in two dimensions
pub struct Length2 {
    pub a: Point2,
    pub b: Point2,
    pub direction: UnitVec2,
}

impl Length2 {
    pub fn new(a: Point2, b: Point2, direction: Option<UnitVec2>) -> Self {
        let direction = direction.unwrap_or(UnitVec2::new_normalize(b - a));
        Self { a, b, direction }
    }

    pub fn reversed(&self) -> Self {
        Self {
            a: self.b,
            b: self.a,
            direction: -self.direction,
        }
    }

    pub fn center(&self) -> SurfacePoint2 {
        let v = self.a - self.b;
        SurfacePoint2::new(self.b + v * 0.5, self.direction)
    }
}

impl Dimension for Length2 {
    fn value(&self) -> f64 {
        self.direction.dot(&(self.b - self.a))
    }
}
