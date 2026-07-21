use crate::common::PCoords;
use crate::geom3::circle3::Circle3;
use crate::geom3::line3::Line3;
use crate::{Iso3, Point3, Result, SurfacePoint3, UnitVec3};
use std::f64::consts::PI;
use std::ops;

// TODO: Can this be replaced with another method or computed adaptively
/// Epsilon for bounds check, it would be good to be able to replace this with something else
const EPSILON: f64 = 1e-10;

/// A cylinder in 3D space, defined by a starting point, a unit direction, a radius, and a length.
/// The cylinder's axis runs from `center` to `center + direction * length`, and the cylindrical
/// surface is every point at distance `radius` from that axis within its extent.
///
/// Because the axis direction is stored as a unit vector rather than being implied by a
/// non-zero length, `length` and `radius` may independently be zero: a zero-length cylinder is a
/// flat disc, a zero-radius cylinder is a line segment, and a cylinder with both zero is a point.
///
/// This is one of `engeom`'s 3D geometric primitives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder3 {
    pub center: Point3,
    pub direction: UnitVec3,
    pub radius: f64,
    pub length: f64,
}

impl Cylinder3 {
    /// Create a cylinder from a starting point, a unit axis direction, a radius, and a length.
    ///
    /// # Arguments
    ///
    /// * `center`: the point at the base of the cylinder's axis, in world space
    /// * `direction`: the direction of the cylinder's axis
    /// * `radius`: the radius of the cylinder
    /// * `length`: the full length of the cylinder along its axis; the cylinder extends from
    ///   `center` to `center + direction * length`
    ///
    /// returns: Cylinder3
    pub fn new(center: Point3, direction: UnitVec3, radius: f64, length: f64) -> Self {
        Self {
            center,
            direction,
            radius,
            length,
        }
    }

    /// Create a cylinder whose axis runs between two points, with the given radius. `center` is
    /// set to `p0`, the direction points from `p0` towards `p1`, and the length is the distance
    /// between them.
    ///
    /// Returns an error if `p0` and `p1` are coincident, since the axis direction would then be
    /// undefined.
    ///
    /// # Arguments
    ///
    /// * `p0`: the point at the center of the cylinder's starting cap
    /// * `p1`: the point at the center of the cylinder's ending cap
    /// * `radius`: the radius of the cylinder
    ///
    /// returns: Result<Cylinder3>
    pub fn from_points(p0: &impl PCoords<3>, p1: &impl PCoords<3>, radius: f64) -> Result<Self> {
        let a = Point3::from(p0.coords());
        let b = Point3::from(p1.coords());
        let diff = b - a;
        let direction = UnitVec3::try_new(diff, 1e-10).ok_or("Points are coincident")?;
        Ok(Self::new(a, direction, radius, diff.norm()))
    }

    /// Returns the radius of the cylinder.
    pub fn r(&self) -> f64 {
        self.radius
    }

    /// Returns the point at the center of the cylinder's starting cap. Identical to `center`;
    /// provided for consistency with the `Segment3` endpoint API.
    pub fn a(&self) -> Point3 {
        self.center
    }

    /// Returns the point at the center of the cylinder's ending cap, at `center + direction * length`.
    pub fn b(&self) -> Point3 {
        self.center + self.direction.into_inner() * self.length
    }

    /// Returns the infinite line running through the cylinder's axis, in the direction of
    /// `direction`.
    pub fn axis(&self) -> Line3 {
        Line3::new(self.center, self.direction.into_inner())
    }

    /// Returns the circle bounding the starting cap of the cylinder, with its normal pointing
    /// outward (opposite `self.direction`).
    pub fn start_cap(&self) -> Circle3 {
        Circle3::new(self.a(), -self.direction, self.radius)
    }

    /// Returns the circle bounding the ending cap of the cylinder, with its normal pointing
    /// outward (the same direction as `self.direction`).
    pub fn end_cap(&self) -> Circle3 {
        Circle3::new(self.b(), self.direction, self.radius)
    }

    /// Returns the volume of the (solid) cylinder.
    pub fn volume(&self) -> f64 {
        PI * self.radius * self.radius * self.length
    }

    /// Returns the area of the cylinder's lateral (side) surface
    pub fn lateral_area(&self) -> f64 {
        2.0 * PI * self.radius * self.length
    }

    /// Returns `true` if `test_point` lies within (or on the boundary of) the solid cylinder.
    ///
    /// # Arguments
    ///
    /// * `test_point`: a point in world space to test
    ///
    /// returns: bool
    pub fn contains_point(&self, test_point: &impl PCoords<3>) -> bool {
        let axis = self.axis();
        let radial = axis.distance_to(test_point);
        let axial = axis.scalar_project(test_point);
        axial > -EPSILON && axial < self.length + EPSILON && radial < self.radius + EPSILON
    }

    /// Returns a new cylinder transformed by the given isometry, without modifying the original.
    pub fn transformed_by(&self, iso: &Iso3) -> Self {
        Self {
            center: iso * self.center,
            direction: UnitVec3::new_normalize(iso.rotation * self.direction.into_inner()),
            radius: self.radius,
            length: self.length,
        }
    }

    /// Returns a new cylinder occupying the same physical volume but with the axis direction
    /// reversed, without modifying the original. The starting point (`center`/`a()`) becomes the
    /// old ending point (`b()`) and vice versa.
    pub fn reversed(&self) -> Self {
        Self {
            center: self.b(),
            direction: -self.direction,
            radius: self.radius,
            length: self.length,
        }
    }

    /// Project the test position onto the cylinder and return a `SurfacePoint3` at the site of the
    /// closest point. Will return a `None` if the test point lies on the cylinder axis, where all
    /// points are equally close.
    ///
    /// Set `infinite` to `true` if you want the closest point to be found to any point at `radius`
    /// distance from the central axis, even if it's beyond the ends of the cylinder. Set `infinite`
    /// to `false` if you want the result to be clamped to the cylinder region.
    ///
    /// This query models a cylinder with only lateral walls, no end caps.
    ///
    /// # Arguments
    ///
    /// * `test_point`: the position in 3D space to test
    /// * `infinite`: set to `false` if the result should be clamped to within the cylinder ends,
    ///   otherwise set to `true`
    ///
    /// returns: Option<SurfacePoint<3>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn closest_point(
        &self,
        test_point: &impl PCoords<3>,
        infinite: bool,
    ) -> Option<SurfacePoint3> {
        let axis = self.axis();

        // The axial parameter, in this case the distance of the test point along the cylinder axis
        let axial_t = axis.scalar_project(test_point);

        // The radial vector to the test point
        let radial_t = Point3::from(test_point.coords()) - axis.at(axial_t);

        // The radial vector to the closest point
        let radial_cn = UnitVec3::try_new(radial_t, EPSILON)?;
        let radial_c = radial_cn.into_inner() * self.radius;

        // The axial parameter of the closest point
        let axial_c = if infinite {
            axial_t
        } else {
            axial_t.clamp(0.0, self.length)
        };

        Some(SurfacePoint3::new(axis.at(axial_c) + radial_c, radial_cn))
    }
}

impl ops::Mul<Cylinder3> for Iso3 {
    type Output = Cylinder3;
    fn mul(self, rhs: Cylinder3) -> Cylinder3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<&Cylinder3> for Iso3 {
    type Output = Cylinder3;
    fn mul(self, rhs: &Cylinder3) -> Cylinder3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<Cylinder3> for &Iso3 {
    type Output = Cylinder3;
    fn mul(self, rhs: Cylinder3) -> Cylinder3 {
        rhs.transformed_by(self)
    }
}

impl ops::Mul<&Cylinder3> for &Iso3 {
    type Output = Cylinder3;
    fn mul(self, rhs: &Cylinder3) -> Cylinder3 {
        rhs.transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use crate::geom3::tests::RandomGeometry3;
    use approx::assert_relative_eq;

    fn axis_cylinder() -> Cylinder3 {
        Cylinder3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            2.0,
            10.0,
        )
    }

    fn tilted_cylinder() -> Cylinder3 {
        let center = Point3::new(1.0, 2.0, 3.0);
        let direction = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 1.0));
        Cylinder3::new(center, direction, 1.5, 6.0)
    }

    fn random_cylinder() -> Cylinder3 {
        let mut rg = RandomGeometry3::new();
        let r = rg.f64(0.5, 4.0);
        let l = rg.f64(1.0, 10.0);
        let center = rg.point3(10.0);
        let direction = rg.unit_vec3();
        Cylinder3::new(center, direction, r, l)
    }

    #[test]
    fn new_stores_fields() {
        let center = Point3::new(1.0, 2.0, 3.0);
        let direction = UnitVec3::new_normalize(Vector3::z());
        let cyl = Cylinder3::new(center, direction, 2.0, 5.0);
        assert_relative_eq!(cyl.center, center, epsilon = 1e-12);
        assert_relative_eq!(
            cyl.direction.into_inner(),
            direction.into_inner(),
            epsilon = 1e-12
        );
        assert_relative_eq!(cyl.r(), 2.0, epsilon = 1e-12);
        assert_relative_eq!(cyl.length, 5.0, epsilon = 1e-12);
    }

    #[test]
    fn from_points_matches_endpoints() {
        let p0 = Point3::new(0.0, 0.0, 0.0);
        let p1 = Point3::new(0.0, 0.0, 4.0);
        let cyl = Cylinder3::from_points(&p0, &p1, 1.5).unwrap();
        assert_relative_eq!(cyl.center, p0, epsilon = 1e-12);
        assert_relative_eq!(cyl.length, 4.0, epsilon = 1e-12);
        assert_relative_eq!(cyl.direction.into_inner(), Vector3::z(), epsilon = 1e-12);
        assert_relative_eq!(cyl.a(), p0, epsilon = 1e-12);
        assert_relative_eq!(cyl.b(), p1, epsilon = 1e-12);
    }

    #[test]
    fn from_points_coincident_is_error() {
        let p = Point3::new(1.0, 1.0, 1.0);
        assert!(Cylinder3::from_points(&p, &p, 1.0).is_err());
    }

    #[test]
    fn a_is_center_b_is_offset_by_length() {
        let cyl = tilted_cylinder();
        assert_relative_eq!(cyl.a(), cyl.center, epsilon = 1e-12);
        assert_relative_eq!(
            cyl.b(),
            cyl.center + cyl.direction.into_inner() * cyl.length,
            epsilon = 1e-12
        );
        assert_relative_eq!((cyl.b() - cyl.a()).norm(), cyl.length, epsilon = 1e-10);
    }

    #[test]
    fn axis_passes_through_center() {
        let cyl = tilted_cylinder();
        let axis = cyl.axis();
        assert_relative_eq!(axis.origin, cyl.center, epsilon = 1e-12);
        assert_relative_eq!(axis.direction, cyl.direction.into_inner(), epsilon = 1e-12);
    }

    #[test]
    fn caps_are_at_endpoints_with_outward_normals() {
        let cyl = tilted_cylinder();
        let start = cyl.start_cap();
        let end = cyl.end_cap();

        assert_relative_eq!(start.center, cyl.a(), epsilon = 1e-12);
        assert_relative_eq!(end.center, cyl.b(), epsilon = 1e-12);
        assert_relative_eq!(start.r(), cyl.r(), epsilon = 1e-12);
        assert_relative_eq!(end.r(), cyl.r(), epsilon = 1e-12);
        assert_relative_eq!(
            start.normal.into_inner(),
            -cyl.direction.into_inner(),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            end.normal.into_inner(),
            cyl.direction.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn volume_matches_known_value() {
        // radius 2, length 10 -> pi * 4 * 10
        let cyl = axis_cylinder();
        assert_relative_eq!(cyl.volume(), PI * 4.0 * 10.0, epsilon = 1e-10);
    }

    #[test]
    fn lateral_area_matches_known_value() {
        let cyl = axis_cylinder();
        assert_relative_eq!(cyl.lateral_area(), 2.0 * PI * 2.0 * 10.0, epsilon = 1e-10);
    }

    #[test]
    fn zero_length_degenerates_to_disc() {
        let cyl = Cylinder3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            2.0,
            0.0,
        );
        assert_relative_eq!(cyl.volume(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(cyl.a(), cyl.b(), epsilon = 1e-12);
        assert_relative_eq!(cyl.a(), cyl.center, epsilon = 1e-12);
    }

    #[test]
    fn zero_radius_degenerates_to_segment() {
        let cyl = Cylinder3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            0.0,
            4.0,
        );
        assert_relative_eq!(cyl.volume(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(cyl.lateral_area(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn zero_length_and_radius_is_a_point() {
        let cyl = Cylinder3::new(
            Point3::new(1.0, 2.0, 3.0),
            UnitVec3::new_normalize(Vector3::z()),
            0.0,
            0.0,
        );
        assert_relative_eq!(cyl.a(), cyl.center, epsilon = 1e-12);
        assert_relative_eq!(cyl.b(), cyl.center, epsilon = 1e-12);
        assert_relative_eq!(cyl.volume(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn contains_point_center_is_inside() {
        let cyl = axis_cylinder();
        assert!(cyl.contains_point(&cyl.center));
    }

    #[test]
    fn contains_point_outside_radius_is_false() {
        let cyl = axis_cylinder();
        let outside = cyl.center + Vector3::new(cyl.r() + 1.0, 0.0, 0.0);
        assert!(!cyl.contains_point(&outside));
    }

    #[test]
    fn contains_point_beyond_length_is_false() {
        let cyl = axis_cylinder();
        let beyond = cyl.b() + cyl.direction.into_inner() * 1.0;
        assert!(!cyl.contains_point(&beyond));
    }

    #[test]
    fn contains_point_before_start_is_false() {
        let cyl = axis_cylinder();
        let before = cyl.a() - cyl.direction.into_inner() * 1.0;
        assert!(!cyl.contains_point(&before));
    }

    #[test]
    fn contains_point_on_boundary_is_true() {
        let cyl = axis_cylinder();
        assert!(cyl.contains_point(&cyl.a()));
        assert!(cyl.contains_point(&cyl.b()));
        let mid = cyl.center + cyl.direction.into_inner() * (cyl.length * 0.5);
        let rim = mid + Vector3::new(cyl.r(), 0.0, 0.0);
        assert!(cyl.contains_point(&rim));
    }

    // -------------------------------------------------------------------------
    // closest_point tests
    // -------------------------------------------------------------------------

    #[test]
    fn closest_point_on_axis_is_none() {
        let cyl = axis_cylinder();
        let mid = cyl.center + cyl.direction.into_inner() * (cyl.length * 0.5);
        assert!(cyl.closest_point(&mid, false).is_none());
        assert!(cyl.closest_point(&mid, true).is_none());
    }

    #[test]
    fn closest_point_on_axis_beyond_length_is_still_none() {
        // Ambiguity is about the radial direction, not the axial position, so points on the axis
        // beyond the cylinder's ends are just as ambiguous as ones in the middle.
        let cyl = axis_cylinder();
        let beyond = cyl.b() + cyl.direction.into_inner() * 5.0;
        assert!(cyl.closest_point(&beyond, false).is_none());
        assert!(cyl.closest_point(&beyond, true).is_none());
    }

    #[test]
    fn closest_point_within_length_matches_for_both_modes() {
        let cyl = axis_cylinder(); // center origin, axis +z, r=2, length=10
        let mid = cyl.center + cyl.direction.into_inner() * (cyl.length * 0.5);
        let test = mid + Vector3::new(5.0, 0.0, 0.0);

        let clamped = cyl.closest_point(&test, false).unwrap();
        let infinite = cyl.closest_point(&test, true).unwrap();
        let expected_point = mid + Vector3::new(2.0, 0.0, 0.0);

        assert_relative_eq!(clamped.point, expected_point, epsilon = 1e-10);
        assert_relative_eq!(infinite.point, expected_point, epsilon = 1e-10);
        assert_relative_eq!(clamped.normal.into_inner(), Vector3::x(), epsilon = 1e-10);
        assert_relative_eq!(infinite.normal.into_inner(), Vector3::x(), epsilon = 1e-10);
    }

    #[test]
    fn closest_point_beyond_end_infinite_extends_past_length() {
        let cyl = axis_cylinder(); // length 10, from z=0 to z=10
        let test = Point3::new(3.0, 0.0, 15.0); // beyond the +z end
        let sp = cyl.closest_point(&test, true).unwrap();
        assert_relative_eq!(sp.point, Point3::new(2.0, 0.0, 15.0), epsilon = 1e-10);
        assert_relative_eq!(sp.normal.into_inner(), Vector3::x(), epsilon = 1e-10);
    }

    #[test]
    fn closest_point_beyond_end_clamped_stops_at_length() {
        let cyl = axis_cylinder();
        let test = Point3::new(3.0, 0.0, 15.0);
        let sp = cyl.closest_point(&test, false).unwrap();
        assert_relative_eq!(sp.point, Point3::new(2.0, 0.0, 10.0), epsilon = 1e-10);
        assert_relative_eq!(sp.normal.into_inner(), Vector3::x(), epsilon = 1e-10);
    }

    #[test]
    fn closest_point_before_start_infinite_extends_before_zero() {
        let cyl = axis_cylinder();
        let test = Point3::new(3.0, 0.0, -15.0); // before the start
        let sp = cyl.closest_point(&test, true).unwrap();
        assert_relative_eq!(sp.point, Point3::new(2.0, 0.0, -15.0), epsilon = 1e-10);
        assert_relative_eq!(sp.normal.into_inner(), Vector3::x(), epsilon = 1e-10);
    }

    #[test]
    fn closest_point_before_start_clamped_stops_at_zero() {
        let cyl = axis_cylinder();
        let test = Point3::new(3.0, 0.0, -15.0);
        let sp = cyl.closest_point(&test, false).unwrap();
        assert_relative_eq!(sp.point, Point3::new(2.0, 0.0, 0.0), epsilon = 1e-10);
        assert_relative_eq!(sp.normal.into_inner(), Vector3::x(), epsilon = 1e-10);
    }

    #[test]
    fn closest_point_clamped_and_infinite_share_radial_direction_beyond_end() {
        // Clamping only moves the result along the axis; the radial direction (and thus the
        // normal) should be identical to the unclamped result.
        let cyl = axis_cylinder();
        let test = Point3::new(3.0, 4.0, 15.0);
        let clamped = cyl.closest_point(&test, false).unwrap();
        let infinite = cyl.closest_point(&test, true).unwrap();
        assert_relative_eq!(
            clamped.normal.into_inner(),
            infinite.normal.into_inner(),
            epsilon = 1e-10
        );
    }

    #[test]
    fn closest_point_zero_radius_is_axis_projection() {
        let cyl = Cylinder3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            0.0,
            10.0,
        );
        let test = Point3::new(3.0, 0.0, 5.0);
        let sp = cyl.closest_point(&test, false).unwrap();
        assert_relative_eq!(sp.point, Point3::new(0.0, 0.0, 5.0), epsilon = 1e-10);
    }

    #[test]
    fn stress_closest_point_infinite_is_always_at_radius_and_matches_projection() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..500 {
            let cyl = random_cylinder();
            let test = rg.point3(20.0);
            let axis = cyl.axis();
            let Some(sp) = cyl.closest_point(&test, true) else {
                continue;
            };

            assert_relative_eq!(axis.distance_to(&sp.point), cyl.r(), epsilon = 1e-8);
            // The axial position of the result must match the test point's own (unclamped)
            // projection onto the axis, since the infinite mode never moves it.
            assert_relative_eq!(
                axis.scalar_project(&sp.point),
                axis.scalar_project(&test),
                epsilon = 1e-8
            );
            // The normal must be perpendicular to the cylinder's axis.
            assert_relative_eq!(sp.normal.dot(&cyl.direction), 0.0, epsilon = 1e-8);
        }
    }

    #[test]
    fn stress_closest_point_clamped_stays_within_length_and_at_radius() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..500 {
            let cyl = random_cylinder();
            let test = rg.point3(20.0);
            let axis = cyl.axis();
            let Some(sp) = cyl.closest_point(&test, false) else {
                continue;
            };

            assert_relative_eq!(axis.distance_to(&sp.point), cyl.r(), epsilon = 1e-8);
            let axial_result = axis.scalar_project(&sp.point);
            assert!(axial_result >= -1e-8 && axial_result <= cyl.length + 1e-8);
        }
    }

    // -------------------------------------------------------------------------
    // Transformation tests
    // -------------------------------------------------------------------------

    #[test]
    fn transformed_by_identity_preserves_all() {
        let cyl = tilted_cylinder();
        let result = cyl.transformed_by(&Iso3::identity());
        assert_relative_eq!(result.r(), cyl.r(), epsilon = 1e-12);
        assert_relative_eq!(result.length, cyl.length, epsilon = 1e-12);
        assert_relative_eq!(result.center, cyl.center, epsilon = 1e-12);
        assert_relative_eq!(
            result.direction.into_inner(),
            cyl.direction.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn stress_transformed_by() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..1000 {
            let original = random_cylinder();
            let iso = rg.iso3(10.0);

            let moved = original.transformed_by(&iso);
            assert_relative_eq!(moved.r(), original.r(), epsilon = 1e-9);
            assert_relative_eq!(moved.length, original.length, epsilon = 1e-9);
            assert_relative_eq!(moved.center, iso * original.center, epsilon = 1e-9);
            assert_relative_eq!(moved.a(), iso * original.a(), epsilon = 1e-9);
            assert_relative_eq!(moved.b(), iso * original.b(), epsilon = 1e-9);
        }
    }

    #[test]
    fn iso_mul_operator_matches_transformed_by() {
        let cyl = tilted_cylinder();
        let iso = Iso3::translation(1.0, 2.0, 3.0);
        let a = cyl.transformed_by(&iso);
        let b = iso * cyl;
        assert_relative_eq!(a.center, b.center, epsilon = 1e-12);
        assert_relative_eq!(
            a.direction.into_inner(),
            b.direction.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn reversed_reverses_direction_preserves_radius_length_and_shape() {
        let cyl = tilted_cylinder();
        let reversed = cyl.reversed();
        assert_relative_eq!(reversed.r(), cyl.r(), epsilon = 1e-12);
        assert_relative_eq!(reversed.length, cyl.length, epsilon = 1e-12);
        assert_relative_eq!(
            reversed.direction.into_inner(),
            -cyl.direction.into_inner(),
            epsilon = 1e-12
        );
        // Reversing the direction swaps which endpoint is `a()` vs `b()`, but the physical
        // cylinder occupies the same space.
        assert_relative_eq!(reversed.a(), cyl.b(), epsilon = 1e-10);
        assert_relative_eq!(reversed.b(), cyl.a(), epsilon = 1e-10);
    }

    #[test]
    fn diagonal_direction_sanity() {
        // Basic sanity check that a non-axis-aligned direction behaves as expected: the length is
        // projected along the (already unit) direction vector.
        let direction = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 0.0));
        let cyl = Cylinder3::new(Point3::origin(), direction, 1.0, std::f64::consts::SQRT_2);
        assert_relative_eq!(cyl.b().coords.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(cyl.b().coords.y, 1.0, epsilon = 1e-10);
    }
}
