use crate::Result;
use crate::common::Line;
use crate::common::PCoords;
use crate::common::points::dist;
use crate::na::{AbstractRotation, Isometry, Point, SVector};
use serde::{Deserialize, Serialize};
use std::ops;

/// A line segment in D-dimensional space, defined by two endpoints.
///
/// `Segment<D>` is the base for `Segment2` and `Segment3`, which are two of `engeom`'s geometric
/// primitives.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct Segment<const D: usize> {
    pub a: Point<f64, D>,
    pub b: Point<f64, D>,
}

impl<const D: usize> Segment<D> {
    pub fn new_unchecked(a: Point<f64, D>, b: Point<f64, D>) -> Self {
        Self { a, b }
    }

    /// Create a new segment from two points, returning an error if the points are coincident
    /// (within a tolerance of `1e-12`).
    pub fn new(a: &impl PCoords<D>, b: &impl PCoords<D>) -> Result<Self> {
        if dist(a, b) < 1e-12 {
            Err("The two points are too close to each other".into())
        } else {
            Ok(Self::new_unchecked(
                Point::from(a.coords()),
                Point::from(b.coords()),
            ))
        }
    }

    pub fn dir(&self) -> SVector<f64, D> {
        self.b - self.a
    }

    pub fn at(&self, t: f64) -> Point<f64, D> {
        self.a + t * self.dir()
    }

    pub fn length(&self) -> f64 {
        self.dir().norm()
    }

    /// Returns a new segment with the endpoints reversed.
    pub fn reversed(&self) -> Self {
        Self::new_unchecked(self.b, self.a)
    }

    /// Returns the infinite line passing through this segment's endpoints, in the direction from
    /// `a` to `b`.
    pub fn to_line(&self) -> Line<D> {
        Line::from_points(&self.a, &self.b)
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
        self.at(t)
    }

    /// Returns a new segment with both endpoints transformed by the given isometry.
    pub fn transformed_by<R: AbstractRotation<f64, D>>(&self, iso: &Isometry<f64, R, D>) -> Self {
        Self {
            a: iso * self.a,
            b: iso * self.b,
        }
    }
}

impl<const D: usize, R: AbstractRotation<f64, D>> ops::Mul<Segment<D>> for Isometry<f64, R, D> {
    type Output = Segment<D>;

    fn mul(self, rhs: Segment<D>) -> Self::Output {
        rhs.transformed_by(&self)
    }
}
