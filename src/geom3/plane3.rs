use crate::geom3::UnitVec3;
use crate::{Iso3, Point3, SurfacePoint3};

#[derive(Debug, Clone)]
pub struct Plane3 {
    pub normal: UnitVec3,
    pub d: f64,
}

impl Plane3 {
    pub fn new(normal: UnitVec3, d: f64) -> Self {
        Self { normal, d }
    }

    /// Create a new plane which is in the same position as the input plane, but with the normal
    /// direction inverted.
    pub fn inverted_normal(&self) -> Self {
        Self::new(-self.normal, -self.d)
    }

    /// Measure and return the signed distance from the plane to a point in 3D space. The sign of
    /// the distance indicates whether the point is above or below the plane according to the
    /// plane's normal vector.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn signed_distance_to_point(&self, point: &Point3) -> f64 {
        self.normal.dot(&point.coords) - self.d
    }


    /// Measure and return the distance from the plane to a point in 3D space. The distance is
    /// always positive, and indicates the shortest distance from the point to the plane. If you
    /// need to know whether the point is above or below the plane, use `signed_distance_to_point`.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn distance_to_point(&self, point: &Point3) -> f64 {
        self.signed_distance_to_point(point).abs()
    }

    /// Project a point onto the plane, returning a point in 3D space which lies on the plane. This
    /// is also the closest point on the plane to the input point.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: OPoint<f64, Const<3>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn project_point(&self, point: &Point3) -> Point3 {
        point - self.normal.into_inner() * self.signed_distance_to_point(point)
    }

    /// Transform the plane by an isometry
    ///
    /// # Arguments
    ///
    /// * `iso`: The isometry to transform the plane by
    ///
    /// returns: Plane3
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn transform_by(&self, iso: &Iso3) -> Self {
        let pos = self.normal.into_inner() * self.d;
        let repr = SurfacePoint3::new(pos.into(), self.normal);

        let new_repr = repr.transformed(&iso);
        Self::from(&new_repr)
    }
}

// TODO: should this be a Result?
impl From<(&Point3, &Point3, &Point3)> for Plane3 {
    /// Create a Plane3 from three points
    ///
    /// # Arguments
    ///
    /// * `(p1, p2, p3)`:
    ///
    /// returns: Plane3
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn from((p1, p2, p3): (&Point3, &Point3, &Point3)) -> Self {
        let normal = UnitVec3::new_normalize((p2 - p1).cross(&(p3 - p1)));
        Self::from((&normal, p1))
    }
}

impl From<(&UnitVec3, &Point3)> for Plane3 {
    fn from((normal, point): (&UnitVec3, &Point3)) -> Self {
        let d = normal.dot(&point.coords);
        Self::new(*normal, d)
    }
}

impl From<&SurfacePoint3> for Plane3 {
    fn from(surface_point: &SurfacePoint3) -> Self {
        Self::from((&surface_point.normal, &surface_point.point))
    }
}
