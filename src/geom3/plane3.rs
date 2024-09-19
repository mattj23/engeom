use crate::geom3::UnitVec3;
use crate::Point3;

pub struct Plane3 {
    pub normal: UnitVec3,
    pub d: f64,
}

impl Plane3 {
    pub fn new(normal: UnitVec3, d: f64) -> Self {
        Self { normal, d }
    }

    pub fn from_point_and_normal(point: &Point3, normal: &UnitVec3) -> Self {
        let d = normal.dot(&point.coords);
        Self::new(*normal, d)
    }

    pub fn from_points(p1: &Point3, p2: &Point3, p3: &Point3) -> Self {
        let normal = UnitVec3::new_normalize((p2 - p1).cross(&(p3 - p1)));
        Self::from_point_and_normal(p1, &normal)
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