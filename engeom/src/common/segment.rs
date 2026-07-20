use crate::common::PCoords;
use crate::na::{AbstractRotation, Isometry, Point, SVector};
use serde::{Deserialize, Serialize};
use std::ops;

/// A line segment in D-dimensional space, defined by two endpoints.
///
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Segment<const D: usize> {
    a: Point<f64, D>,
    b: Point<f64, D>,
}

impl<const D: usize> Segment<D> {
    pub fn new_unchecked(a: Point<f64, D>, b: Point<f64, D>) -> Self {
        Self { a, b }
    }

    pub fn a(&self) -> &Point<f64, D> {
        &self.a
    }

    pub fn b(&self) -> &Point<f64, D> {
        &self.b
    }

    pub fn dir(&self) -> SVector<f64, D> {
        self.b - self.a
    }

    pub fn point_at(&self, t: f64) -> Point<f64, D> {
        self.a + t * self.dir()
    }

    pub fn length(&self) -> f64 {
        self.dir().norm()
    }

    /// Calculate the scalar projection of a set of coordinates onto the line segment, in which
    /// 0.0 represents a point at the segment's starting point `a` and 1.0 represents a point at
    /// the segment's end point `b`.  The result can be any finite value, including negative ones
    /// or ones greater than zero.
    ///
    /// # Arguments
    ///
    /// * `other`: an entity with coordinates to project onto the segment
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn scalar_projection(&self, other: &impl PCoords<D>) -> f64 {
        let dir = self.dir();
        let test = other.coords() - self.a.coords();
        dir.dot(&test) / dir.norm_squared()
    }

    pub fn closest_point(&self, other: &impl PCoords<D>) -> Point<f64, D> {
        let t = self.scalar_projection(other).clamp(0.0, 1.0);
        self.point_at(t)
    }

    pub fn new_transformed<R: AbstractRotation<f64, D>>(&self, iso: &Isometry<f64, R, D>) -> Self {
        Self {
            a: iso * self.a,
            b: iso * self.b,
        }
    }
}

impl<const D: usize, R: AbstractRotation<f64, D>> ops::Mul<Segment<D>> for Isometry<f64, R, D> {
    type Output = Segment<D>;

    fn mul(self, rhs: Segment<D>) -> Self::Output {
        rhs.new_transformed(&self)
    }
}
