use crate::Result;
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
    ///
    /// # Arguments
    ///
    /// * `other`: the point to project onto the line defined by the surface point
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Point2, SurfacePoint2, Vector2};
    /// use approx::assert_relative_eq;
    ///
    /// let sp = SurfacePoint2::new_normalize(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
    ///
    /// let other = Point2::new(-1.0, -1.0);
    /// let scalar_projection = sp.scalar_projection(&other);
    ///
    /// assert_relative_eq!(scalar_projection, -1.0, epsilon = 1e-6);
    /// ```
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

/// Created a vector of `SurfacePoint` instances from a vector of points and a vector of normals.
/// If the number of points and normals are not the same, an error is returned.
///
/// # Arguments
///
/// * `points`: the vector of points, ordered to match the normals
/// * `normals`: the vector of normals, ordered to match the points
///
/// returns: Result<Vec<SurfacePoint<{ D }>, Global>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
/// use engeom::{Point2, Vector2};
/// use engeom::common::surface_point::surface_point_vector;
/// use engeom::geom2::UnitVec2;
///
/// let points = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 1.0)];
/// let normals = vec![Vector2::new(0.0, 1.0), Vector2::new(1.0, 0.0)];
///
/// let surface_points = surface_point_vector(&points, &normals).unwrap();
///
/// assert_eq!(surface_points[0].point, Point2::new(0.0, 0.0));
/// assert_eq!(surface_points[0].normal, UnitVec2::new_normalize(Vector2::new(0.0, 1.0)));
/// assert_eq!(surface_points[1].point, Point2::new(1.0, 1.0));
/// assert_eq!(surface_points[1].normal, UnitVec2::new_normalize(Vector2::new(1.0, 0.0)));
/// ```
pub fn surface_point_vector<const D: usize>(
    points: &[Point<f64, D>],
    normals: &[SVector<f64, D>],
) -> Result<Vec<SurfacePoint<D>>> {
    // Check that the number of points and normals are the same
    if points.len() != normals.len() {
        return Err("The number of points and normals must be the same".into());
    }

    Ok(points
        .iter()
        .zip(normals.iter())
        .map(|(p, n)| SurfacePoint::new_normalize(*p, *n))
        .collect())
}

pub trait SurfacePointCollection<const D: usize> {
    fn clone_points(&self) -> Vec<Point<f64, D>>;
    fn clone_normals(&self) -> Vec<Unit<SVector<f64, D>>>;
}
