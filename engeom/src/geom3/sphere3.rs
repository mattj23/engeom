use crate::common::PCoords;
use crate::geom3::circle3::Circle3;
use crate::geom3::plane3::Plane3;
use crate::{Iso3, Point3, Result, SurfacePoint3, UnitVec3};
use parry3d_f64::na::{Translation3, UnitQuaternion};
use parry3d_f64::query::{Ray, RayCast};
use parry3d_f64::shape::Ball;
use std::ops;

/// A sphere in 3D space, defined by a center point and a radius.
#[derive(Debug, Clone)]
pub struct Sphere3 {
    center: Point3,
    radius: f64,
}

impl Sphere3 {
    pub fn new(center: Point3, radius: f64) -> Self {
        Self { center, radius }
    }

    pub fn center(&self) -> Point3 {
        self.center
    }

    pub fn r(&self) -> f64 {
        self.radius
    }

    /// Returns a new sphere transformed by the given isometry. Only the translation component
    /// affects the sphere (rotation does not change a sphere).
    ///
    /// # Arguments
    ///
    /// * `iso`: an isometry that will be applied to the sphere's center point
    ///
    /// returns: Sphere3
    pub fn new_transformed_by(&self, iso: &Iso3) -> Self {
        Self::new(iso * self.center, self.radius)
    }

    /// Transforms this sphere in place by the given isometry. Only the translation component
    /// affects the sphere (rotation does not change a sphere).
    ///
    /// # Arguments
    ///
    /// * `iso`: an isometry that will be applied to the sphere's center point
    ///
    /// returns: ()
    pub fn transform_by(&mut self, iso: &Iso3) {
        self.center = iso * self.center;
    }

    /// Returns the closest point on the sphere's surface to `test_point`, along with the outward
    /// surface normal at that point, or `None` if `test_point` is at the exact center of the
    /// sphere.
    ///
    ///
    /// # Arguments
    ///
    /// * `test_point`: a point in world space to test
    ///
    /// returns: Option<SurfacePoint<3>>
    pub fn closest_point(&self, test_point: &impl PCoords<3>) -> Option<SurfacePoint3> {
        let v = Point3::from(test_point.coords()) - self.center;
        let dist = v.norm();
        if dist < 1e-10 {
            return None;
        }
        let normal = UnitVec3::new_normalize(v);
        let point = self.center + normal.into_inner() * self.radius;
        Some(SurfacePoint3::new(point, normal))
    }

    /// Intersects the sphere with a plane, returning the resulting `Circle3`, or `None` if the
    /// plane does not intersect the sphere.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane to intersect with
    ///
    /// returns: Option<Circle3>
    pub fn intersect_plane(&self, plane: &Plane3) -> Option<Circle3> {
        let dist = plane.signed_distance_to_point(&self.center);
        if dist.abs() > self.radius {
            return None;
        }
        let circle_radius = (self.radius * self.radius - dist * dist).sqrt();
        let circle_center = self.center - plane.normal.into_inner() * dist;
        Circle3::from_point_normal(&circle_center, &plane.normal, circle_radius).ok()
    }

    /// Intersects this sphere with another sphere, returning the resulting `Circle3`, or `None` if
    /// the spheres do not intersect (too far apart, one inside the other, or concentric).
    ///
    /// # Arguments
    ///
    /// * `other`: the other sphere to intersect with
    ///
    /// returns: Option<Circle3>
    pub fn intersect_sphere(&self, other: &Sphere3) -> Option<Circle3> {
        let axis = other.center - self.center;
        let d = axis.norm();
        if d < 1e-10 || d > self.radius + other.radius || d < (self.radius - other.radius).abs() {
            return None;
        }
        let h = (d * d + self.radius * self.radius - other.radius * other.radius) / (2.0 * d);
        let circle_radius = (self.radius * self.radius - h * h).max(0.0).sqrt();
        let normal = UnitVec3::new_normalize(axis);
        let circle_center = self.center + normal.into_inner() * h;
        Circle3::from_point_normal(&circle_center, &normal, circle_radius).ok()
    }

    /// Intersects a ray with the sphere. Returns the first intersection as a surface point with
    /// the outward normal, or `None` if the ray does not intersect.
    pub fn ray_intersection(&self, ray: &Ray) -> Option<SurfacePoint3> {
        let iso = Iso3::from_parts(
            Translation3::from(self.center.coords),
            UnitQuaternion::identity(),
        );
        let hit = Ball::new(self.radius).cast_ray_and_get_normal(&iso, ray, f64::MAX, false)?;
        let point = ray.origin + ray.dir * hit.time_of_impact;
        let normal = UnitVec3::new_normalize(hit.normal);
        Some(SurfacePoint3::new(point, normal))
    }
}

impl ops::Mul<Sphere3> for Iso3 {
    type Output = Sphere3;
    fn mul(self, rhs: Sphere3) -> Sphere3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<&Sphere3> for Iso3 {
    type Output = Sphere3;
    fn mul(self, rhs: &Sphere3) -> Sphere3 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<Sphere3> for &Iso3 {
    type Output = Sphere3;
    fn mul(self, rhs: Sphere3) -> Sphere3 {
        rhs.new_transformed_by(self)
    }
}

impl ops::Mul<&Sphere3> for &Iso3 {
    type Output = Sphere3;
    fn mul(self, rhs: &Sphere3) -> Sphere3 {
        rhs.new_transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn unit_sphere() -> Sphere3 {
        Sphere3::new(Point3::origin(), 1.0)
    }

    fn offset_sphere() -> Sphere3 {
        Sphere3::new(Point3::new(1.0, 2.0, 3.0), 2.0)
    }

    #[test]
    fn transform_by_moves_center() {
        let sphere = Sphere3::new(Point3::origin(), 2.0);
        let iso = Iso3::translation(1.0, 2.0, 3.0);
        let moved = sphere.new_transformed_by(&iso);
        assert_relative_eq!(moved.center(), Point3::new(1.0, 2.0, 3.0), epsilon = 1e-12);
        assert_relative_eq!(moved.r(), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn closest_point_outside_is_on_surface() {
        let sphere = unit_sphere();
        let test = Point3::new(3.0, 0.0, 0.0);
        let sp = sphere.closest_point(&test).unwrap();
        assert_relative_eq!(sp.point, Point3::new(1.0, 0.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(sp.normal.into_inner(), crate::Vector3::x(), epsilon = 1e-12);
    }

    #[test]
    fn closest_point_inside_is_on_surface() {
        let sphere = unit_sphere();
        let test = Point3::new(0.5, 0.0, 0.0);
        let sp = sphere.closest_point(&test).unwrap();
        // should still project to the surface in the same direction
        assert_relative_eq!(
            (sp.point - sphere.center()).norm(),
            sphere.r(),
            epsilon = 1e-12
        );
        assert_relative_eq!(sp.point, Point3::new(1.0, 0.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn closest_point_offset_sphere() {
        let sphere = offset_sphere();
        // test point directly above the center along Z
        let test = Point3::new(1.0, 2.0, 10.0);
        let sp = sphere.closest_point(&test).unwrap();
        assert_relative_eq!(sp.point, Point3::new(1.0, 2.0, 5.0), epsilon = 1e-12);
        assert_relative_eq!(sp.normal.into_inner(), crate::Vector3::z(), epsilon = 1e-12);
    }

    #[test]
    fn closest_point_at_center_returns_none() {
        let sphere = unit_sphere();
        assert!(sphere.closest_point(&Point3::origin()).is_none());
    }

    #[test]
    fn plane_through_center_gives_great_circle() {
        let sphere = unit_sphere();
        let plane = Plane3::xy(); // z=0 through origin
        let circle = sphere.intersect_plane(&plane).unwrap();
        assert_relative_eq!(circle.r(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(circle.center(), Point3::origin(), epsilon = 1e-12);
        assert_relative_eq!(
            circle.normal().into_inner(),
            crate::Vector3::z(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn plane_tangent_gives_point_circle() {
        let sphere = unit_sphere();
        // plane z=1 is tangent at the top of the unit sphere
        let plane = Plane3::new(crate::Vector3::z_axis(), 1.0);
        let circle = sphere.intersect_plane(&plane).unwrap();
        assert_relative_eq!(circle.r(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(circle.center(), Point3::new(0.0, 0.0, 1.0), epsilon = 1e-12);
    }

    #[test]
    fn plane_misses_sphere_returns_none() {
        let sphere = unit_sphere();
        let plane = Plane3::new(crate::Vector3::z_axis(), 2.0); // z=2, outside sphere
        let result = sphere.intersect_plane(&plane);
        assert!(result.is_none());
    }

    #[test]
    fn plane_intersection_circle_radius_correct() {
        // Sphere radius 5, plane at distance 3 from center → circle radius = sqrt(25-9) = 4
        let sphere = Sphere3::new(Point3::origin(), 5.0);
        let plane = Plane3::new(crate::Vector3::z_axis(), 3.0);
        let circle = sphere.intersect_plane(&plane).unwrap();
        assert_relative_eq!(circle.r(), 4.0, epsilon = 1e-12);
    }

    #[test]
    fn plane_intersection_circle_center_on_plane() {
        let sphere = offset_sphere(); // center (1,2,3), radius 2
        let plane = Plane3::new(crate::Vector3::z_axis(), 3.0); // z=3 passes through center
        let circle = sphere.intersect_plane(&plane).unwrap();
        // circle center should be the projection of sphere center onto the plane
        assert_relative_eq!(circle.center(), Point3::new(1.0, 2.0, 3.0), epsilon = 1e-12);
        assert_relative_eq!(circle.r(), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn plane_intersection_circle_normal_matches_plane() {
        let sphere = unit_sphere();
        let normal = crate::UnitVec3::new_normalize(crate::Vector3::new(1.0, 1.0, 0.0));
        let plane = Plane3::new(normal, 0.0);
        let circle = sphere.intersect_plane(&plane).unwrap();
        assert_relative_eq!(
            circle.normal().into_inner(),
            normal.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn two_unit_spheres_adjacent_gives_great_circle() {
        // Two unit spheres with centers 2 apart touch at origin; intersection is a point (r=0)
        // But two unit spheres with centers sqrt(2) apart: h = sqrt(2)/2, circle_r = sqrt(1 - 0.5) = sqrt(0.5)
        let s1 = Sphere3::new(Point3::new(-1.0, 0.0, 0.0), 1.0);
        let s2 = Sphere3::new(Point3::new(1.0, 0.0, 0.0), 1.0);
        // d=2 = r1+r2, tangent externally → point circle at origin
        let circle = s1.intersect_sphere(&s2).unwrap();
        assert_relative_eq!(circle.r(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(circle.center(), Point3::origin(), epsilon = 1e-12);
    }

    #[test]
    fn equal_spheres_offset_along_x() {
        // Two equal spheres radius 5, centers 6 apart: h = 3, circle_r = sqrt(25-9) = 4
        let s1 = Sphere3::new(Point3::new(-3.0, 0.0, 0.0), 5.0);
        let s2 = Sphere3::new(Point3::new(3.0, 0.0, 0.0), 5.0);
        let circle = s1.intersect_sphere(&s2).unwrap();
        assert_relative_eq!(circle.r(), 4.0, epsilon = 1e-12);
        // circle center should be at the midpoint
        assert_relative_eq!(circle.center(), Point3::origin(), epsilon = 1e-12);
        // normal should point along X
        assert_relative_eq!(circle.normal().into_inner().x.abs(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn unequal_spheres_circle_center_correct() {
        // s1 radius 3 at origin, s2 radius 4 at (5,0,0)
        // d=5, h=(25+9-16)/10 = 18/10 = 1.8, circle_r = sqrt(9-3.24) = sqrt(5.76) = 2.4
        let s1 = Sphere3::new(Point3::origin(), 3.0);
        let s2 = Sphere3::new(Point3::new(5.0, 0.0, 0.0), 4.0);
        let circle = s1.intersect_sphere(&s2).unwrap();
        assert_relative_eq!(circle.r(), 2.4, epsilon = 1e-10);
        assert_relative_eq!(circle.center(), Point3::new(1.8, 0.0, 0.0), epsilon = 1e-10);
    }

    #[test]
    fn spheres_too_far_apart_returns_none() {
        let s1 = Sphere3::new(Point3::origin(), 1.0);
        let s2 = Sphere3::new(Point3::new(10.0, 0.0, 0.0), 1.0);
        assert!(s1.intersect_sphere(&s2).is_none());
    }

    #[test]
    fn sphere_inside_other_returns_none() {
        let s1 = Sphere3::new(Point3::origin(), 5.0);
        let s2 = Sphere3::new(Point3::new(1.0, 0.0, 0.0), 1.0);
        assert!(s1.intersect_sphere(&s2).is_none());
    }

    #[test]
    fn concentric_spheres_returns_none() {
        let s1 = Sphere3::new(Point3::origin(), 1.0);
        let s2 = Sphere3::new(Point3::origin(), 2.0);
        assert!(s1.intersect_sphere(&s2).is_none());
    }

    #[test]
    fn intersection_circle_lies_on_both_spheres() {
        // Verify several points on the circle lie on both sphere surfaces
        let s1 = Sphere3::new(Point3::new(0.0, 0.0, 0.0), 5.0);
        let s2 = Sphere3::new(Point3::new(4.0, 0.0, 0.0), 4.0);
        let circle = s1.intersect_sphere(&s2).unwrap();
        for angle in [0.0_f64, 1.0, 2.0, 3.0, 4.0, 5.0] {
            let pt = circle.at_angle(angle).point;
            assert_relative_eq!((pt - s1.center()).norm(), s1.r(), epsilon = 1e-10);
            assert_relative_eq!((pt - s2.center()).norm(), s2.r(), epsilon = 1e-10);
        }
    }

    #[test]
    fn ray_hits_front_of_unit_sphere() {
        let sphere = unit_sphere();
        let ray = Ray::new(
            Point3::new(5.0, 0.0, 0.0),
            (-crate::Vector3::x()).normalize(),
        );
        let sp = sphere.ray_intersection(&ray).unwrap();
        assert_relative_eq!(sp.point, Point3::new(1.0, 0.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(sp.normal.into_inner(), crate::Vector3::x(), epsilon = 1e-12);
    }

    #[test]
    fn ray_misses_returns_none() {
        let sphere = unit_sphere();
        let ray = Ray::new(
            Point3::new(5.0, 5.0, 0.0),
            (-crate::Vector3::x()).normalize(),
        );
        assert!(sphere.ray_intersection(&ray).is_none());
    }

    #[test]
    fn ray_intersection_point_is_on_surface() {
        let sphere = offset_sphere();
        let origin = Point3::new(10.0, 2.0, 3.0);
        let dir = (sphere.center() - origin).normalize();
        let ray = Ray::new(origin, dir);
        let sp = sphere.ray_intersection(&ray).unwrap();
        let dist = (sp.point - sphere.center()).norm();
        assert_relative_eq!(dist, sphere.r(), epsilon = 1e-12);
    }

    #[test]
    fn ray_intersection_normal_points_outward() {
        let sphere = unit_sphere();
        let ray = Ray::new(
            Point3::new(0.0, 5.0, 0.0),
            (-crate::Vector3::y()).normalize(),
        );
        let sp = sphere.ray_intersection(&ray).unwrap();
        // normal at the top of the sphere should point in +Y
        assert_relative_eq!(sp.normal.into_inner(), crate::Vector3::y(), epsilon = 1e-12);
    }

    #[test]
    fn ray_from_inside_hits_back_surface() {
        let sphere = unit_sphere();
        // Ray starting inside, heading +X, should hit the back (+X side)
        let ray = Ray::new(Point3::new(0.0, 0.0, 0.0), crate::Vector3::x().normalize());
        let sp = sphere.ray_intersection(&ray).unwrap();
        assert_relative_eq!(sp.point, Point3::new(1.0, 0.0, 0.0), epsilon = 1e-12);
    }
}
