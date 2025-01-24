use crate::geom3::UnitVec3;
use crate::{Point3, SurfacePoint3};

pub struct Plane3 {
    pub normal: UnitVec3,
    pub d: f64,
}

impl Plane3 {
    pub fn new(normal: UnitVec3, d: f64) -> Self {
        Self { normal, d }
    }

    pub fn signed_distance_to_point(&self, point: &Point3) -> f64 {
        self.normal.dot(&point.coords) - self.d
    }

    pub fn distance_to_point(&self, point: &Point3) -> f64 {
        self.signed_distance_to_point(point).abs()
    }

    pub fn project_point(&self, point: &Point3) -> Point3 {
        point - self.normal.into_inner() * self.signed_distance_to_point(point)
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