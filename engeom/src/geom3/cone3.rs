use crate::common::PCoords;
use crate::geom3::circle3::Circle3;
use crate::geom3::line3::Line3;
use crate::{Iso3, Point3, Result, SurfacePoint3, UnitVec3};
use std::f64::consts::PI;
use std::ops;

mod fitting;

/// Epsilon for bounds/degeneracy checks.
const EPSILON: f64 = 1e-10;

/// A cone in 3D space, defined by a tip (apex) point, a unit direction, a height, and a
/// base radius. The cone's axis runs from `tip` to `tip + direction * height`, and the lateral
/// surface is every point whose distance from the axis grows linearly from `0` at `tip` to
/// `radius` at the base.
///
/// Because the axis direction is stored as a unit vector rather than being implied by a
/// non-zero height, `height` and `radius` may independently be zero: a zero-height cone is a
/// flat disc, a zero-radius cone is a line segment, and a cone with both zero is a point.
///
/// This is one of `engeom`'s 3D geometric primitives.
///
/// > [!NOTE]
/// > I have some misgivings about this representation. The separate height/radius makes it quick
/// > to retrieve a unit speed axis line and allows for independent zero-height and zero-length
/// > cones, but those are degenerate cases, and it makes it so we need a separate representation
/// > when we do fitting in the `./fitting.rs` module.  I'm not fully sold on the cost/benefit here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone3 {
    pub tip: Point3,
    pub direction: UnitVec3,
    pub height: f64,
    pub radius: f64,
}

impl Cone3 {
    /// Create a cone from a tip (apex) point, a unit axis direction, a height, and a base
    /// radius.
    ///
    /// # Arguments
    ///
    /// * `tip`: the apex of the cone, in world space
    /// * `direction`: the direction of the cone's axis, from the tip towards the base
    /// * `height`: the distance from the tip to the base along the axis; the base lies at
    ///   `tip + direction * height`
    /// * `radius`: the radius of the circular base
    ///
    /// returns: Cone3
    pub fn new(tip: Point3, direction: UnitVec3, height: f64, radius: f64) -> Self {
        Self {
            tip,
            direction,
            height,
            radius,
        }
    }

    /// Create a cone from a tip point and a base center point, with the given base radius. The
    /// axis points from `tip` towards `base_center`, and the height is the distance between
    /// them.
    ///
    /// Returns an error if `tip` and `base_center` are coincident, since the axis direction
    /// would then be undefined.
    ///
    /// # Arguments
    ///
    /// * `tip`: the apex of the cone
    /// * `base_center`: the center point of the cone's circular base
    /// * `radius`: the radius of the circular base
    ///
    /// returns: Result<Cone3>
    pub fn from_points(
        tip: &impl PCoords<3>,
        base_center: &impl PCoords<3>,
        radius: f64,
    ) -> Result<Self> {
        let t = Point3::from(tip.coords());
        let b = Point3::from(base_center.coords());
        let diff = b - t;
        let direction = UnitVec3::try_new(diff, EPSILON).ok_or("Points are coincident")?;
        Ok(Self::new(t, direction, diff.norm(), radius))
    }

    /// Returns the radius of the cone's base.
    pub fn r(&self) -> f64 {
        self.radius
    }

    /// Returns the center point of the cone's base, at `tip + direction * height`.
    pub fn base_center(&self) -> Point3 {
        self.tip + self.direction.into_inner() * self.height
    }

    /// Returns the infinite line running through the cone's axis, from the tip in the direction
    /// of `direction`.
    pub fn axis(&self) -> Line3 {
        Line3::new(self.tip, self.direction.into_inner())
    }

    /// Returns the circle bounding the base of the cone, with its normal pointing outward (the
    /// same direction as `self.direction`).
    pub fn base(&self) -> Circle3 {
        Circle3::new(self.base_center(), self.direction, self.radius)
    }

    /// Returns the half-angle of the cone: the angle between the axis and the lateral surface,
    /// in radians.
    pub fn half_angle(&self) -> f64 {
        self.radius.atan2(self.height)
    }

    /// Returns the slant height of the cone: the distance from the tip to a point on the rim of
    /// the base, along the lateral surface.
    pub fn slant_height(&self) -> f64 {
        (self.height * self.height + self.radius * self.radius).sqrt()
    }

    /// Returns the volume of the (solid) cone.
    pub fn volume(&self) -> f64 {
        PI * self.radius * self.radius * self.height / 3.0
    }

    /// Returns the area of the cone's lateral surface, excluding the base.
    pub fn lateral_area(&self) -> f64 {
        PI * self.radius * self.slant_height()
    }

    /// Returns `true` if `test_point` lies within (or on the boundary of) the solid cone.
    ///
    /// # Arguments
    ///
    /// * `test_point`: a point in world space to test
    ///
    /// returns: bool
    pub fn contains_point(&self, test_point: &impl PCoords<3>) -> bool {
        let axis = self.axis();
        let axial = axis.scalar_project(test_point);
        if axial < -EPSILON || axial > self.height + EPSILON {
            return false;
        }

        let radial = axis.distance_to(test_point);
        let allowed = if self.height > EPSILON {
            self.radius * (axial / self.height)
        } else {
            self.radius
        };
        radial <= allowed + EPSILON
    }

    /// Returns a new cone transformed by the given isometry, without modifying the original.
    pub fn transformed_by(&self, iso: &Iso3) -> Self {
        Self {
            tip: iso * self.tip,
            direction: UnitVec3::new_normalize(iso.rotation * self.direction.into_inner()),
            height: self.height,
            radius: self.radius,
        }
    }

    /// Project the test position onto the cone's lateral surface and return a `SurfacePoint3` at
    /// the site of the closest point. If the `infinite` argument is set to `false`, the result
    /// will be clamped between the tip and the base rim, otherwise the result will be anywhere
    /// on the hourglass shape formed by an infinite cone.
    ///
    /// Returns `None` if the test point lies on the cone's axis, where all directions around the
    /// surface are equally close, or if the cone is degenerate (zero height and radius), where the
    /// lateral surface collapses to a single point with no well-defined normal.
    ///
    /// This query models a cone with only a lateral surface, no base cap.
    ///
    /// # Arguments
    ///
    /// * `test_point`:
    /// * `infinite`:
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
        if self.height <= EPSILON {
            return None;
        }

        // This starts similar to the cylinder, but we find the slant line: the line that lies
        // on the lateral wall of the cone in the plane formed by the cone axis and the test
        // point. The closest point on the infinite cone will be the projection of the test point
        // onto the slant line.
        let axis = self.axis();
        let axial_t = axis.scalar_project(test_point);

        let radial_dir = UnitVec3::try_new(
            Point3::from(test_point.coords()) - axis.at(axial_t),
            EPSILON,
        )?;
        let slant_line = Line3::from_points(
            &self.tip,
            &(self.base_center() + radial_dir.into_inner() * self.radius),
        );

        // Now we'll get the projection of the test point onto the slant line and clamp it if
        // necessary.
        let t_inf = slant_line.scalar_project(test_point);
        let t = if infinite {
            t_inf
        } else {
            t_inf.clamp(0.0, 1.0)
        };

        // Now we have to get the cone's surface normal. To do that we have to rotate the slant
        // line's direction 90° clockwise in the plane of the cone axis and the test point. To find
        // that we'll use the cross product of the radial direction and the axis.
        let normal_rot = radial_dir.into_inner().cross(&axis.direction).normalize();
        let normal = UnitVec3::try_new(
            Iso3::rotation(normal_rot * -PI / 2.0) * slant_line.direction,
            EPSILON,
        )?;

        Some(SurfacePoint3::new(slant_line.at(t), normal))
    }
}

impl ops::Mul<Cone3> for Iso3 {
    type Output = Cone3;
    fn mul(self, rhs: Cone3) -> Cone3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<&Cone3> for Iso3 {
    type Output = Cone3;
    fn mul(self, rhs: &Cone3) -> Cone3 {
        rhs.transformed_by(&self)
    }
}

impl ops::Mul<Cone3> for &Iso3 {
    type Output = Cone3;
    fn mul(self, rhs: Cone3) -> Cone3 {
        rhs.transformed_by(self)
    }
}

impl ops::Mul<&Cone3> for &Iso3 {
    type Output = Cone3;
    fn mul(self, rhs: &Cone3) -> Cone3 {
        rhs.transformed_by(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use crate::common::random_geometry::RandomGeometry3;
    use approx::assert_relative_eq;

    fn axis_cone() -> Cone3 {
        // Tip at origin, opening up +z, height 10, base radius 2 (at z=10).
        Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            10.0,
            2.0,
        )
    }

    fn tilted_cone() -> Cone3 {
        let tip = Point3::new(1.0, 2.0, 3.0);
        let direction = UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 1.0));
        Cone3::new(tip, direction, 6.0, 1.5)
    }

    fn random_cone() -> Cone3 {
        let mut rg = RandomGeometry3::new();
        let r = rg.f64(0.5, 4.0);
        let h = rg.f64(1.0, 10.0);
        let tip = rg.point(10.0);
        let direction = rg.unit_vec();
        Cone3::new(tip, direction, h, r)
    }

    #[test]
    fn new_stores_fields() {
        let tip = Point3::new(1.0, 2.0, 3.0);
        let direction = UnitVec3::new_normalize(Vector3::z());
        let cone = Cone3::new(tip, direction, 5.0, 2.0);
        assert_relative_eq!(cone.tip, tip, epsilon = 1e-12);
        assert_relative_eq!(
            cone.direction.into_inner(),
            direction.into_inner(),
            epsilon = 1e-12
        );
        assert_relative_eq!(cone.height, 5.0, epsilon = 1e-12);
        assert_relative_eq!(cone.r(), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn from_points_matches_endpoints() {
        let tip = Point3::new(0.0, 0.0, 0.0);
        let base = Point3::new(0.0, 0.0, 4.0);
        let cone = Cone3::from_points(&tip, &base, 1.5).unwrap();
        assert_relative_eq!(cone.tip, tip, epsilon = 1e-12);
        assert_relative_eq!(cone.height, 4.0, epsilon = 1e-12);
        assert_relative_eq!(cone.direction.into_inner(), Vector3::z(), epsilon = 1e-12);
        assert_relative_eq!(cone.base_center(), base, epsilon = 1e-12);
    }

    #[test]
    fn from_points_coincident_is_error() {
        let p = Point3::new(1.0, 1.0, 1.0);
        assert!(Cone3::from_points(&p, &p, 1.0).is_err());
    }

    #[test]
    fn base_center_is_offset_by_height() {
        let cone = tilted_cone();
        assert_relative_eq!(
            cone.base_center(),
            cone.tip + cone.direction.into_inner() * cone.height,
            epsilon = 1e-12
        );
    }

    #[test]
    fn axis_line_passes_through_tip() {
        let cone = tilted_cone();
        let axis = cone.axis();
        assert_relative_eq!(axis.origin, cone.tip, epsilon = 1e-12);
        assert_relative_eq!(axis.direction, cone.direction.into_inner(), epsilon = 1e-12);
    }

    #[test]
    fn base_is_at_base_center_with_outward_normal() {
        let cone = tilted_cone();
        let base = cone.base();
        assert_relative_eq!(base.center, cone.base_center(), epsilon = 1e-12);
        assert_relative_eq!(base.r(), cone.r(), epsilon = 1e-12);
        assert_relative_eq!(
            base.normal.into_inner(),
            cone.direction.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn half_angle_matches_known_value() {
        // height 1, radius 1 -> 45 degree half angle
        let cone = Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            1.0,
            1.0,
        );
        assert_relative_eq!(cone.half_angle(), PI / 4.0, epsilon = 1e-12);
    }

    #[test]
    fn slant_height_matches_known_value() {
        // height 3, radius 4 -> slant 5
        let cone = Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            3.0,
            4.0,
        );
        assert_relative_eq!(cone.slant_height(), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn volume_matches_known_value() {
        // radius 2, height 10 -> pi * 4 * 10 / 3
        let cone = axis_cone();
        assert_relative_eq!(cone.volume(), PI * 4.0 * 10.0 / 3.0, epsilon = 1e-10);
    }

    #[test]
    fn lateral_area_matches_known_value() {
        // radius 3, height 4 -> slant 5, area = pi * 3 * 5
        let cone = Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            4.0,
            3.0,
        );
        assert_relative_eq!(cone.lateral_area(), PI * 3.0 * 5.0, epsilon = 1e-10);
    }

    #[test]
    fn zero_height_degenerates_to_disc() {
        let cone = Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            0.0,
            2.0,
        );
        assert_relative_eq!(cone.volume(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(cone.base_center(), cone.tip, epsilon = 1e-12);
    }

    #[test]
    fn zero_radius_degenerates_to_segment() {
        let cone = Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            4.0,
            0.0,
        );
        assert_relative_eq!(cone.volume(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(cone.lateral_area(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn contains_point_tip_is_inside() {
        let cone = axis_cone();
        assert!(cone.contains_point(&cone.tip));
    }

    #[test]
    fn contains_point_base_center_is_inside() {
        let cone = axis_cone();
        assert!(cone.contains_point(&cone.base_center()));
    }

    #[test]
    fn contains_point_outside_radius_at_base_is_false() {
        let cone = axis_cone();
        let outside = cone.base_center() + Vector3::new(cone.r() + 1.0, 0.0, 0.0);
        assert!(!cone.contains_point(&outside));
    }

    #[test]
    fn contains_point_within_radius_at_base_is_true() {
        let cone = axis_cone();
        let inside = cone.base_center() + Vector3::new(cone.r() - 0.1, 0.0, 0.0);
        assert!(cone.contains_point(&inside));
    }

    #[test]
    fn contains_point_at_full_radius_partway_up_is_false() {
        // At half the height, the allowed radius is half the base radius.
        let cone = axis_cone();
        let halfway = cone.tip + cone.direction.into_inner() * (cone.height * 0.5);
        let at_base_radius = halfway + Vector3::new(cone.r(), 0.0, 0.0);
        assert!(!cone.contains_point(&at_base_radius));
        let at_half_radius = halfway + Vector3::new(cone.r() * 0.5, 0.0, 0.0);
        assert!(cone.contains_point(&at_half_radius));
    }

    #[test]
    fn contains_point_beyond_base_is_false() {
        let cone = axis_cone();
        let beyond = cone.base_center() + cone.direction.into_inner() * 1.0;
        assert!(!cone.contains_point(&beyond));
    }

    #[test]
    fn contains_point_before_tip_is_false() {
        let cone = axis_cone();
        let before = cone.tip - cone.direction.into_inner() * 1.0;
        assert!(!cone.contains_point(&before));
    }

    // -------------------------------------------------------------------------
    // closest_point tests
    // -------------------------------------------------------------------------

    #[test]
    fn closest_point_on_axis_is_none() {
        let cone = axis_cone();
        let on_axis = cone.tip + cone.direction.into_inner() * 5.0;
        assert!(cone.closest_point(&on_axis, false).is_none());
        assert!(cone.closest_point(&on_axis, true).is_none());
    }

    #[test]
    fn closest_point_zero_height_is_always_none() {
        // Height <= epsilon means the cone can't form a slant line at all, regardless of radius.
        let cone = Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            0.0,
            2.0,
        );
        let test = Point3::new(5.0, 0.0, 0.0);
        assert!(cone.closest_point(&test, false).is_none());
        assert!(cone.closest_point(&test, true).is_none());
    }

    #[test]
    fn closest_point_zero_radius_projects_onto_axis() {
        // A zero-radius cone degenerates to its axis segment; the "surface" point should just be
        // the perpendicular projection of the test point onto the axis.
        let cone = Cone3::new(
            Point3::origin(),
            UnitVec3::new_normalize(Vector3::z()),
            10.0,
            0.0,
        );
        let test = Point3::new(3.0, 0.0, 4.0);
        let sp = cone.closest_point(&test, false).unwrap();
        assert_relative_eq!(sp.point, Point3::new(0.0, 0.0, 4.0), epsilon = 1e-10);
    }

    #[test]
    fn closest_point_on_surface_within_bounds_matches_for_both_modes() {
        let cone = axis_cone(); // tip at origin, axis +z, height 10, radius 2
        // At axial position z=5 (t=0.5), the cone's radius is exactly 1, so this point already
        // lies on the lateral surface.
        let on_surface = Point3::new(1.0, 0.0, 5.0);
        let clamped = cone.closest_point(&on_surface, false).unwrap();
        let infinite = cone.closest_point(&on_surface, true).unwrap();
        assert_relative_eq!(clamped.point, on_surface, epsilon = 1e-10);
        assert_relative_eq!(infinite.point, on_surface, epsilon = 1e-10);
    }

    #[test]
    fn closest_point_beyond_base_clamps_to_rim_but_infinite_extends_past_it() {
        let cone = axis_cone(); // tip at origin, axis +z, height 10, radius 2
        // The slant line from the tip through the rim point (2, 0, 10) has direction (2, 0, 10);
        // this test point sits at twice that vector from the tip, i.e. exactly on the infinite
        // extension of the slant line beyond the rim.
        let test = cone.tip + Vector3::new(4.0, 0.0, 20.0);

        let clamped = cone.closest_point(&test, false).unwrap();
        let rim = cone.base_center() + Vector3::new(cone.r(), 0.0, 0.0);
        assert_relative_eq!(clamped.point, rim, epsilon = 1e-10);

        // Since the test point already lies exactly on the infinite cone's surface, the
        // unclamped query should return the test point itself.
        let infinite = cone.closest_point(&test, true).unwrap();
        assert_relative_eq!(infinite.point, test, epsilon = 1e-10);
    }

    #[test]
    fn closest_point_before_tip_clamps_to_tip_but_infinite_extends_through_apex() {
        let cone = axis_cone(); // tip at origin, axis +z, height 10, radius 2
        // A point behind (below) the tip, off-axis.
        let test = Point3::new(1.0, 0.0, -5.0);

        let clamped = cone.closest_point(&test, false).unwrap();
        assert_relative_eq!(clamped.point, cone.tip, epsilon = 1e-10);

        // In infinite mode, the query is free to project onto the slant line's negative-t
        // extension through the apex, landing on a different point than the clamped result.
        let infinite = cone.closest_point(&test, true).unwrap();
        assert!(
            (infinite.point - clamped.point).norm() > 1e-6,
            "infinite mode should extend past the tip rather than clamping to it"
        );

        // Wherever it lands, the result must still be a valid closest-point projection: the
        // test point must lie on the line through the surface point along its normal.
        let normal_line = Line3::new(infinite.point, infinite.normal.into_inner());
        assert_relative_eq!(normal_line.distance_to(&test), 0.0, epsilon = 1e-8);
    }

    #[test]
    fn closest_point_normal_is_outward_and_perpendicular_to_slant() {
        let cone = axis_cone();
        let test = Point3::new(5.0, 0.0, 3.0);
        let sp = cone.closest_point(&test, true).unwrap();

        let radial_dir = Vector3::new(1.0, 0.0, 0.0); // matches the test point's azimuth
        let tangent = cone.direction.into_inner() * cone.height + radial_dir * cone.radius;
        assert_relative_eq!(sp.normal.dot(&tangent), 0.0, epsilon = 1e-8);
        assert!(
            sp.normal.dot(&radial_dir) > 0.0,
            "normal should point away from the axis, not towards it"
        );
    }

    #[test]
    fn stress_closest_point_clamped_stays_within_bounds() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..500 {
            let cone = random_cone();
            let test = rg.point(20.0);
            let Some(sp) = cone.closest_point(&test, false) else {
                continue;
            };

            let axis = cone.axis();
            let axial = axis.scalar_project(&sp.point);
            let radial = axis.distance_to(&sp.point);

            assert!(axial >= -1e-8 && axial <= cone.height + 1e-8);
            let expected_radial = cone.radius * (axial / cone.height);
            assert_relative_eq!(radial, expected_radial, epsilon = 1e-6);
        }
    }

    #[test]
    fn stress_closest_point_infinite_test_point_is_colinear_with_surface_normal() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..500 {
            let cone = random_cone();
            let test = rg.point(20.0);
            let Some(sp) = cone.closest_point(&test, true) else {
                continue;
            };

            // Check colinearity
            assert_relative_eq!(sp.planar_distance(&test), 0.0, epsilon = 1e-8);

            // The returned point should also actually lie on the infinite (double-napped) cone
            // surface: its radial distance from the axis should match the cone's slope applied
            // to the absolute axial distance from the tip, since points on the far nappe (beyond
            // the tip) mirror the near one.
            let axis = cone.axis();
            let axial = axis.scalar_project(&sp.point);
            let radial = axis.distance_to(&sp.point);
            let expected_radial = cone.radius * axial.abs() / cone.height;
            assert_relative_eq!(radial, expected_radial, epsilon = 1e-6);
        }
    }

    // -------------------------------------------------------------------------
    // Transformation tests
    // -------------------------------------------------------------------------

    #[test]
    fn transformed_by_identity_preserves_all() {
        let cone = tilted_cone();
        let result = cone.transformed_by(&Iso3::identity());
        assert_relative_eq!(result.r(), cone.r(), epsilon = 1e-12);
        assert_relative_eq!(result.height, cone.height, epsilon = 1e-12);
        assert_relative_eq!(result.tip, cone.tip, epsilon = 1e-12);
        assert_relative_eq!(
            result.direction.into_inner(),
            cone.direction.into_inner(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn stress_transformed_by() {
        let mut rg = RandomGeometry3::new();
        for _ in 0..1000 {
            let original = random_cone();
            let iso = rg.iso3(10.0);

            let moved = original.transformed_by(&iso);
            assert_relative_eq!(moved.r(), original.r(), epsilon = 1e-9);
            assert_relative_eq!(moved.height, original.height, epsilon = 1e-9);
            assert_relative_eq!(moved.tip, iso * original.tip, epsilon = 1e-9);
            assert_relative_eq!(
                moved.base_center(),
                iso * original.base_center(),
                epsilon = 1e-9
            );
        }
    }

    #[test]
    fn iso_mul_operator_matches_transformed_by() {
        let cone = tilted_cone();
        let iso = Iso3::translation(1.0, 2.0, 3.0);
        let a = cone.transformed_by(&iso);
        let b = iso * cone;
        assert_relative_eq!(a.tip, b.tip, epsilon = 1e-12);
        assert_relative_eq!(
            a.direction.into_inner(),
            b.direction.into_inner(),
            epsilon = 1e-12
        );
    }
}
