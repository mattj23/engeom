//! This module has representations of different types of dimensions

use crate::common::surface_point::SurfacePoint;
use parry3d_f64::na::{Point, SVector, Unit};

pub trait Measurement {
    fn value(&self) -> f64;
}

/// Represents a length measurement in two dimensions
pub struct Length<const D: usize> {
    pub a: Point<f64, D>,
    pub b: Point<f64, D>,
    pub direction: Unit<SVector<f64, D>>,
}

impl<const D: usize> Length<D> {
    pub fn new(
        a: Point<f64, D>,
        b: Point<f64, D>,
        direction: Option<Unit<SVector<f64, D>>>,
    ) -> Self {
        let direction = direction.unwrap_or(Unit::new_normalize(b - a));
        Self { a, b, direction }
    }

    pub fn reversed(&self) -> Self {
        Self {
            a: self.b,
            b: self.a,
            direction: -self.direction,
        }
    }

    pub fn center(&self) -> SurfacePoint<D> {
        let v = self.a - self.b;
        SurfacePoint::new(self.b + v * 0.5, self.direction)
    }
}

impl<const D: usize> Measurement for Length<D> {
    fn value(&self) -> f64 {
        self.direction.dot(&(self.b - self.a))
    }
}
