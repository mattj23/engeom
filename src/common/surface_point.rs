use parry3d_f64::na::{AbstractRotation, Isometry, Point, SVector, Unit};
use serde::{Deserialize, Serialize};

/// A `SurfacePoint` is a struct which is used to represent a point on a surface (n-1 dimensional
/// manifold) in n-dimensional space. It is defined by a point and a normal vector. Mathematically,
/// a `SurfacePoint` is identical to a parameterized line or a ray with a unit direction. It also
/// uniquely defines half-spaces (so a plane in 3D and a half-space line in 2D).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SurfacePoint<const D: usize> {
    pub point: Point<f64, D>,
    pub normal: Unit<SVector<f64, D>>,
}

impl<const D: usize> SurfacePoint<D> {
    pub fn new(point: Point<f64, D>, normal: Unit<SVector<f64, D>>) -> Self {
        Self { point, normal }
    }

    pub fn new_normalize(point: Point<f64, D>, normal: SVector<f64, D>) -> Self {
        Self::new(point, Unit::new_normalize(normal))
    }

    /// Returns the point offset from the surface point by the given distance along the normal
    pub fn at_distance(&self, distance: f64) -> Point<f64, D> {
        self.point + self.normal.as_ref() * distance
    }

    /// Returns the scalar projection value of another point onto the line defined by the point and
    /// normal. This can be interpreted as the physical distance along the normal line that the
    /// other point is from the surface point.
    pub fn scalar_projection(&self, other: &Point<f64, D>) -> f64 {
        self.normal.dot(&(other - self.point))
    }

    /// Returns the point on the line defined by the point and normal that is closest to the other
    /// point, aka the projection of the other point onto the line defined by this surface point.
    pub fn projection(&self, other: &Point<f64, D>) -> Point<f64, D> {
        self.at_distance(self.scalar_projection(other))
    }

    /// Returns a new surface point with the same point but with the normal reversed
    pub fn reversed(&self) -> Self {
        Self::new(self.point, -self.normal)
    }

    /// Returns a new surface point transformed by the given isometry
    pub fn transformed<R>(&self, t: &Isometry<f64, R, D>) -> Self
    where
        R: AbstractRotation<f64, D>,
    {
        Self::new(t * self.point, t * self.normal)
    }

    /// Returns the distance between a test point and its projection onto the line defined by the
    /// surface point. This is a complement to the `scalar_projection` method, except that it can
    /// only compute the magnitude of the distance, since the number of other dimensions may be
    /// greater than one.
    pub fn planar_distance(&self, other: &Point<f64, D>) -> f64 {
        let projection = self.projection(other);
        (projection - other).norm()
    }
}

pub trait SurfacePointCollection<const D: usize> {
    fn points(&self) -> Vec<Point<f64, D>>;
    fn normals(&self) -> Vec<Unit<SVector<f64, D>>>;
}
