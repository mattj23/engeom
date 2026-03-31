use crate::Result;
use crate::common::PCoords;
use parry3d_f64::na::{AbstractRotation, Isometry, Point, SVector, Unit};
use serde::{Deserialize, Serialize};

/// A `SurfacePoint` is a struct that is used to represent a point on a surface (n-1 dimensional
/// manifold) in n-dimensional space. It is defined by a point and a normal vector. Mathematically,
/// a `SurfacePoint` is identical to a parameterized line or a ray with a unit direction. It also
/// uniquely defines half-spaces (so a plane in 3D and a half-space line in 2D).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct SurfacePoint<const D: usize> {
    pub point: Point<f64, D>,
    pub normal: Unit<SVector<f64, D>>,
}

impl<const D: usize> SurfacePoint<D> {
    /// Creates a new `SurfacePoint` from a point and an already-normalized unit normal.
    pub fn new(point: Point<f64, D>, normal: Unit<SVector<f64, D>>) -> Self {
        Self { point, normal }
    }

    /// Creates a new `SurfacePoint` from a point and an arbitrary normal vector, normalizing the
    /// vector in the process.
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
    pub fn scalar_projection(&self, other: &impl PCoords<D>) -> f64 {
        self.normal.dot(&(other.coords() - self.point.coords))
    }

    /// Returns the point on the line defined by the point and normal that is closest to the other
    /// point, aka the projection of the other point onto the line defined by this surface point.
    pub fn projection(&self, other: &impl PCoords<D>) -> Point<f64, D> {
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
    pub fn planar_distance(&self, other: &impl PCoords<D>) -> f64 {
        let projection = self.projection(other);
        (projection.coords - other.coords()).norm()
    }

    /// Returns a new surface point shifted from the original surface point by the given distance
    /// along the normal. This is useful for creating a new surface point that is a certain distance
    /// away from the original surface point, in the direction of the normal.
    ///
    /// # Arguments
    ///
    /// * `shift`: the distance to offset the surface point along the normal
    ///
    /// returns: SurfacePoint<{ D }>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Point2, SurfacePoint2, Vector2};
    /// use approx::assert_relative_eq;
    ///
    /// let sp = SurfacePoint2::new_normalize(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
    ///
    /// let shifted = sp.new_shifted(2.0);
    /// assert_relative_eq!(shifted.point, Point2::new(0.0, 2.0), epsilon = 1e-6);
    /// assert_relative_eq!(shifted.normal.into_inner(), Vector2::new(0.0, 1.0), epsilon = 1e-6);
    /// ```
    pub fn new_shifted(&self, offset: f64) -> Self {
        let new_point = self.point + self.normal.as_ref() * offset;
        Self::new(new_point, self.normal)
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

/// A trait for collections of `SurfacePoint` instances that can produce owned copies of their
/// points and normals as flat `Vec`s. Implemented for slices and `Vec` references of typed
/// surface point aliases such as `SurfacePoint2` and `SurfacePoint3`.
pub trait SurfacePointCollection<const D: usize> {
    /// Returns a `Vec` of cloned position points from this collection.
    fn clone_points(&self) -> Vec<Point<f64, D>>;

    /// Returns a `Vec` of cloned unit normals from this collection.
    fn clone_normals(&self) -> Vec<Unit<SVector<f64, D>>>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use parry3d_f64::na::{Point2, Point3, Vector2, Vector3};

    fn sp2(px: f64, py: f64, nx: f64, ny: f64) -> SurfacePoint<2> {
        SurfacePoint::new_normalize(Point2::new(px, py), Vector2::new(nx, ny))
    }

    fn sp3(px: f64, py: f64, pz: f64, nx: f64, ny: f64, nz: f64) -> SurfacePoint<3> {
        SurfacePoint::new_normalize(Point3::new(px, py, pz), Vector3::new(nx, ny, nz))
    }

    #[test]
    fn new_stores_point_and_normal() {
        let sp = sp2(1.0, 2.0, 0.0, 1.0);
        assert_relative_eq!(sp.point, Point2::new(1.0, 2.0));
        assert_relative_eq!(sp.normal.into_inner(), Vector2::new(0.0, 1.0));
    }

    #[test]
    fn new_normalize_normalizes_the_vector() {
        let sp = SurfacePoint::new_normalize(Point2::new(0.0, 0.0), Vector2::new(0.0, 3.0));
        assert_relative_eq!(sp.normal.into_inner().norm(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(
            sp.normal.into_inner(),
            Vector2::new(0.0, 1.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn at_distance_positive_moves_along_normal() {
        let sp = sp2(0.0, 0.0, 0.0, 1.0);
        let p = sp.at_distance(3.0);
        assert_relative_eq!(p, Point2::new(0.0, 3.0), epsilon = 1e-12);
    }

    #[test]
    fn at_distance_negative_moves_opposite_normal() {
        let sp = sp2(1.0, 1.0, 1.0, 0.0);
        let p = sp.at_distance(-2.0);
        assert_relative_eq!(p, Point2::new(-1.0, 1.0), epsilon = 1e-12);
    }

    #[test]
    fn at_distance_zero_returns_the_origin_point() {
        let sp = sp3(3.0, 4.0, 5.0, 0.0, 0.0, 1.0);
        assert_relative_eq!(sp.at_distance(0.0), sp.point, epsilon = 1e-12);
    }

    #[test]
    fn scalar_projection_point_ahead_of_plane() {
        let sp = sp2(0.0, 0.0, 0.0, 1.0); // normal = +y
        let other = Point2::new(0.0, 5.0);
        assert_relative_eq!(sp.scalar_projection(&other), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn scalar_projection_point_behind_plane() {
        let sp = sp2(0.0, 0.0, 0.0, 1.0); // normal = +y
        let other = Point2::new(1.0, -3.0);
        assert_relative_eq!(sp.scalar_projection(&other), -3.0, epsilon = 1e-12);
    }

    #[test]
    fn scalar_projection_point_on_plane_is_zero() {
        let sp = sp2(0.0, 0.0, 1.0, 0.0); // normal = +x, point at origin
        let other = Point2::new(0.0, 7.0); // lies on the plane x = 0
        assert_relative_eq!(sp.scalar_projection(&other), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn scalar_projection_offset_origin() {
        // Surface point not at origin: the projection should be relative to sp.point
        let sp = sp2(0.0, 2.0, 0.0, 1.0); // normal = +y, point at (0, 2)
        let other = Point2::new(0.0, 5.0);
        assert_relative_eq!(sp.scalar_projection(&other), 3.0, epsilon = 1e-12);
    }

    #[test]
    fn projection_lands_on_the_line() {
        let sp = sp2(0.0, 0.0, 0.0, 1.0); // vertical line through origin
        let other = Point2::new(4.0, 3.0);
        let proj = sp.projection(&other);
        assert_relative_eq!(proj, Point2::new(0.0, 3.0), epsilon = 1e-12);
    }

    #[test]
    fn projection_of_point_already_on_line_is_itself() {
        let sp = sp2(0.0, 0.0, 0.0, 1.0);
        let on_line = Point2::new(0.0, 7.0);
        assert_relative_eq!(sp.projection(&on_line), on_line, epsilon = 1e-12);
    }

    #[test]
    fn reversed_flips_normal_keeps_point() {
        let sp = sp2(1.0, 2.0, 0.0, 1.0);
        let rev = sp.reversed();
        assert_relative_eq!(rev.point, sp.point);
        assert_relative_eq!(
            rev.normal.into_inner(),
            -sp.normal.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn reversed_twice_is_identity() {
        let sp = sp3(1.0, 2.0, 3.0, 1.0, 1.0, 0.0);
        let twice = sp.reversed().reversed();
        assert_relative_eq!(twice.point, sp.point);
        assert_relative_eq!(
            twice.normal.into_inner(),
            sp.normal.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn planar_distance_measures_lateral_offset() {
        // Normal is +y; lateral distance is along x
        let sp = sp2(0.0, 0.0, 0.0, 1.0);
        let other = Point2::new(3.0, 7.0);
        assert_relative_eq!(sp.planar_distance(&other), 3.0, epsilon = 1e-12);
    }

    #[test]
    fn planar_distance_zero_for_point_on_line() {
        let sp = sp2(0.0, 0.0, 0.0, 1.0);
        let on_line = Point2::new(0.0, 5.0);
        assert_relative_eq!(sp.planar_distance(&on_line), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn planar_distance_3d() {
        // Normal is +z; lateral distance is in xy-plane
        let sp = sp3(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let other = Point3::new(3.0, 4.0, 10.0); // distance from z-axis in xy = 5
        assert_relative_eq!(sp.planar_distance(&other), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn new_shifted_moves_point_keeps_normal() {
        let sp = sp2(0.0, 0.0, 0.0, 1.0); // normal = +y
        let shifted = sp.new_shifted(4.0);
        assert_relative_eq!(shifted.point, Point2::new(0.0, 4.0), epsilon = 1e-12);
        assert_relative_eq!(
            shifted.normal.into_inner(),
            sp.normal.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn new_shifted_negative_offset() {
        let sp = sp2(0.0, 5.0, 0.0, 1.0); // normal = +y, point at (0, 5)
        let shifted = sp.new_shifted(-3.0);
        assert_relative_eq!(shifted.point, Point2::new(0.0, 2.0), epsilon = 1e-12);
    }

    #[test]
    fn surface_point_vector_builds_correctly() {
        let points = vec![Point2::new(0.0, 0.0), Point2::new(1.0, 0.0)];
        let normals = vec![Vector2::new(0.0, 1.0), Vector2::new(0.0, 1.0)];
        let spv = surface_point_vector(&points, &normals).unwrap();
        assert_eq!(spv.len(), 2);
        assert_relative_eq!(spv[0].point, points[0]);
        assert_relative_eq!(spv[1].point, points[1]);
        assert_relative_eq!(spv[0].normal.into_inner(), normals[0]);
        assert_relative_eq!(spv[1].normal.into_inner(), normals[1]);
    }

    #[test]
    fn surface_point_vector_length_mismatch_returns_error() {
        let points = vec![Point2::new(0.0, 0.0)];
        let normals: Vec<Vector2<f64>> = vec![];
        assert!(surface_point_vector(&points, &normals).is_err());
    }
}
