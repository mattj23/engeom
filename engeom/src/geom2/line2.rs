use crate::AngleDir::Cw;
use crate::common::vec_f64::sort_and_dedup;
use crate::common::{Intersection, PCoords};
use crate::geom2::{Aabb2, Ray2, SurfacePoint2, UnitVec2, rot90, signed_angle};
use crate::{Iso2, Point2, Vector2};
use parry2d_f64::query::Ray;
use std::ops;

/// A parameterized line in 2D space: `P(t) = origin + t * direction`.
///
/// The direction is not required to be normalized; use `new_normalize` for unit-speed
/// parameterization where `t` equals unit length.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Line2 {
    pub origin: Point2,
    pub direction: Vector2,
}

impl Line2 {
    pub fn x_axis() -> Self {
        Self::new(Point2::origin(), Vector2::x())
    }

    pub fn y_axis() -> Self {
        Self::new(Point2::origin(), Vector2::y())
    }

    /// Create a new parallel line with a given offset along the normal direction. A positive
    /// `delta_n` moves the origin to the right of the line, while a negative `delta_n` moves the
    /// origin to the left.
    ///
    /// # Arguments
    ///
    /// * `delta_n`: The offset along the normal direction.
    ///
    /// returns: Line2
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Line2;
    /// use approx::assert_relative_eq;
    ///
    /// let l0 = Line2::new(Point2::new(0.0, 0.0), Vector2::new(0.0, 1.0));
    /// let l1 = l0.new_parallel(1.0);
    ///
    /// assert_relative_eq!(l1.origin, Point2::new(0.0, 1.0), epsilon = 1e-6);
    /// assert_relative_eq!(l1.direction, Vector2::new(0.0, 1.0), epsilon = 1e-6);
    /// ```
    pub fn new_parallel(&self, delta_n: f64) -> Self {
        let n = self.normal().into_inner();
        Self::new(self.origin + n * delta_n, self.direction)
    }

    /// Create a line from an origin point and a direction vector (stored as-is, not normalized).
    pub fn new(origin: Point2, direction: Vector2) -> Self {
        Self { origin, direction }
    }

    /// Create a line from an origin point and a direction vector, normalizing the direction so
    /// that the parameter `t` equals arc length from the origin.
    pub fn new_normalize(origin: Point2, direction: Vector2) -> Self {
        Self::new(origin, direction.normalize())
    }

    /// Create a line through two points. The direction is `p2 - p1` (not normalized).
    pub fn from_points(p1: &impl PCoords<2>, p2: &impl PCoords<2>) -> Self {
        let p1 = Point2::from(p1.coords());
        let p2 = Point2::from(p2.coords());
        Self::new(p1, p2 - p1)
    }

    /// Normalizes the direction vector in place so that `t` equals arc length from the origin.
    pub fn normalize(&mut self) {
        self.direction = self.direction.normalize();
    }

    /// Returns a new line with the same origin but a normalized direction, so that `t` equals arc
    /// length from the origin.
    pub fn normalized(&self) -> Self {
        Self::new(self.origin, self.direction.normalize())
    }

    /// Returns the point at parameter `t`: `P(t) = origin + t * direction`.
    pub fn at(&self, t: f64) -> Point2 {
        self.origin + self.direction * t
    }

    /// Moves the origin of the line by a given amount along the direction of the line. A positive
    /// `delta_t` moves the origin forward along the direction of the line, while a negative
    /// `delta_t` moves the origin backward. The line is modified in place.
    pub fn shift_along(&mut self, delta_t: f64) {
        self.origin += self.direction * delta_t;
    }

    /// Returns a new line with the origin shifted by a given amount along the direction of the
    /// line. The direction of the new line is the same as the original line. The original is left
    /// unchanged.
    pub fn new_shifted_along(&self, delta_t: f64) -> Self {
        Self {
            origin: self.origin + self.direction * delta_t,
            direction: self.direction,
        }
    }

    /// Returns the parameter `t` such that `P(t)` is the closest point on the line to `point`.
    pub fn scalar_project(&self, point: &impl PCoords<2>) -> f64 {
        (Point2::from(point.coords()) - self.origin).dot(&self.direction)
            / self.direction.norm_squared()
    }

    /// Returns the closest point on the line to `point`.
    pub fn closest_point(&self, point: &impl PCoords<2>) -> Point2 {
        self.at(self.scalar_project(point))
    }

    /// Returns the perpendicular distance from `point` to the line.
    pub fn distance_to(&self, point: &impl PCoords<2>) -> f64 {
        let pt = Point2::from(point.coords());
        (pt - self.closest_point(&pt)).norm()
    }

    /// Returns a new line with both origin and direction transformed by the given isometry.
    pub fn new_transformed_by(&self, iso: &Iso2) -> Self {
        let mut clone = self.clone();
        clone.transform_by(iso);
        clone
    }

    /// Transforms this line in place by the given isometry.
    pub fn transform_by(&mut self, iso: &Iso2) {
        self.origin = iso * self.origin;
        self.direction = iso.rotation * self.direction;
    }

    /// Returns the unit normal to this line: the direction rotated 90 degrees clockwise.
    /// By convention, this points to the "right" when traveling in the line's direction,
    /// consistent with outward normals on counter-clockwise-wound 2D geometry.
    pub fn normal(&self) -> UnitVec2 {
        UnitVec2::new_normalize(rot90(Cw) * self.direction)
    }

    /// Returns the signed perpendicular distance from `point` to this line. Positive values
    /// indicate the point is to the right of the direction of travel (on the normal side),
    /// negative values indicate the point is to the left.
    pub fn signed_distance_to(&self, point: &impl PCoords<2>) -> f64 {
        self.normal().dot(&(point.coords() - self.origin.coords))
    }

    /// Returns the intersection point with another `Line2`, or `None` if the lines are parallel.
    pub fn intersect(&self, other: &impl LineOps2) -> Option<Point2> {
        let (t, _) =
            intersection_param(&self.origin, &self.direction, &other.origin(), &other.dir())?;
        Some(self.at(t))
    }
}

impl From<Line2> for SurfacePoint2 {
    fn from(line: Line2) -> Self {
        SurfacePoint2::new_normalize(line.origin, line.direction)
    }
}

impl<T: LineOps2> From<&T> for Line2 {
    fn from(line: &T) -> Self {
        Self::new(line.origin(), line.dir())
    }
}

impl LineOps2 for Line2 {
    fn origin(&self) -> Point2 {
        self.origin
    }

    fn dir(&self) -> Vector2 {
        self.direction
    }

    fn at(&self, t: f64) -> Point2 {
        self.origin + self.direction * t
    }
}

impl ops::Mul<Line2> for Iso2 {
    type Output = Line2;
    fn mul(self, rhs: Line2) -> Line2 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<&Line2> for Iso2 {
    type Output = Line2;
    fn mul(self, rhs: &Line2) -> Line2 {
        rhs.new_transformed_by(&self)
    }
}

impl ops::Mul<Line2> for &Iso2 {
    type Output = Line2;
    fn mul(self, rhs: Line2) -> Line2 {
        rhs.new_transformed_by(self)
    }
}

impl ops::Mul<&Line2> for &Iso2 {
    type Output = Line2;
    fn mul(self, rhs: &Line2) -> Line2 {
        rhs.new_transformed_by(self)
    }
}

/// Compute the intersection parameters between two parameterized lines. Will return None if
/// the two directions are parallel to each other
pub fn intersection_param(
    a0: &Point2,
    ad: &Vector2,
    b0: &Point2,
    bd: &Vector2,
) -> Option<(f64, f64)> {
    let det: f64 = bd.x * ad.y - bd.y * ad.x;
    if det.abs() < 1e-12 {
        return None;
    }

    let dx = b0.x - a0.x;
    let dy = b0.y - a0.y;

    Some(((dy * bd.x - dx * bd.y) / det, (dy * ad.x - dx * ad.y) / det))
}

pub fn intersect_lines(a: &impl LineOps2, b: &impl LineOps2) -> Option<(f64, f64)> {
    intersection_param(&a.origin(), &a.dir(), &b.origin(), &b.dir())
}

pub fn intersect_rays(r0: &Ray, r1: &Ray) -> Option<(f64, f64)> {
    intersection_param(&r0.origin, &r0.dir, &r1.origin, &r1.dir)
}

pub trait LineOps2 {
    fn origin(&self) -> Point2;
    fn dir(&self) -> Vector2;
    fn at(&self, t: f64) -> Point2;

    fn projected_parameter(&self, p: &impl PCoords<2>) -> f64 {
        let n = p.coords() - self.origin().coords;
        self.dir().dot(&n) / self.dir().dot(&self.dir())
    }

    fn projected_point(&self, p: &impl PCoords<2>) -> Point2 {
        self.at(self.projected_parameter(p))
    }

    /// Returns the direction of the vector turned in its orthogonal direction by rotating it
    /// -90 degrees, this is typically used as a normal
    fn orthogonal(&self) -> Vector2 {
        Iso2::rotation(-std::f64::consts::PI / 2.0) * self.dir()
    }

    /// Returns the signed perpendicular distance from `point` to this line. Positive values
    /// indicate the point is to the right of the direction of travel (on the normal side),
    /// negative values indicate the point is to the left.
    fn signed_projection_dist(&self, point: &impl PCoords<2>) -> f64 {
        self.orthogonal()
            .normalize()
            .dot(&(point.coords() - self.origin().coords))
    }

    fn intersection_params(&self, other: &impl LineOps2) -> Option<(f64, f64)> {
        intersection_param(&self.origin(), &self.dir(), &other.origin(), &other.dir())
    }

    /// Creates an isometry where the origin is the same as the line's origin and the X direction
    /// is the same as the line's direction. Transforming an entity at the origin by this isometry
    /// will move it to the same position and orientation as the line. Transforming an entity by
    /// the inverse of this isometry will transfer its relationship with the line to the origin.
    fn to_iso_from_x(&self) -> Iso2 {
        let r = Iso2::rotation(signed_angle(&Vector2::x(), &self.dir()));
        let t = Iso2::translation(self.origin().x, self.origin().y);

        t * r
    }

    /// Creates an isometry where the origin is the same as the line's origin and the Y direction
    /// is the same as the line's direction. Transforming an entity at the origin by this isometry
    /// will move it to the same position and orientation as the line. Transforming an entity by
    /// the inverse of this isometry will transfer its relationship with the line to the origin.
    fn to_iso_from_y(&self) -> Iso2 {
        let r = Iso2::rotation(signed_angle(&Vector2::y(), &self.dir()));
        let t = Iso2::translation(self.origin().x, self.origin().y);

        t * r
    }
}

impl<T: LineOps2> LineOps2 for &T {
    // Blanket impl for all types that implement LineOps2
    fn origin(&self) -> Point2 {
        (**self).origin()
    }

    fn dir(&self) -> Vector2 {
        (**self).dir()
    }

    fn at(&self, t: f64) -> Point2 {
        (**self).at(t)
    }
}

impl LineOps2 for Ray2 {
    fn origin(&self) -> Point2 {
        self.origin
    }

    fn dir(&self) -> Vector2 {
        self.dir
    }

    fn at(&self, t: f64) -> Point2 {
        self.point_at(t)
    }
}

impl<T: LineOps2> Intersection<Aabb2, Vec<f64>> for T {
    fn intersection(&self, other: Aabb2) -> Vec<f64> {
        let x_inv = 1.0 / self.dir().x;
        let y_inv = 1.0 / self.dir().y;

        let t1x = (other.mins.x - self.origin().x) * x_inv;
        let t2x = (other.maxs.x - self.origin().x) * x_inv;
        let t1y = (other.mins.y - self.origin().y) * y_inv;
        let t2y = (other.maxs.y - self.origin().y) * y_inv;

        let mut result = Vec::with_capacity(2);

        if t1x.is_finite() && other.contains_local_point(&self.at(t1x)) {
            result.push(t1x);
        }
        if t2x.is_finite() && other.contains_local_point(&self.at(t2x)) {
            result.push(t2x);
        }
        if t1y.is_finite() && other.contains_local_point(&self.at(t1y)) {
            result.push(t1y);
        }
        if t2y.is_finite() && other.contains_local_point(&self.at(t2y)) {
            result.push(t2y);
        }
        sort_and_dedup(&mut result, Some(1e-12));

        result
    }
}

pub fn slab_method2(bv: &Aabb2, origin: &Point2, n_inv: &Vector2) -> bool {
    let mut t1 = (bv.mins.x - origin.x) * n_inv.x;
    let mut t2 = (bv.maxs.x - origin.x) * n_inv.x;

    let tmin = t1.min(t2);
    let tmax = t1.max(t2);

    t1 = (bv.mins.y - origin.y) * n_inv.y;
    t2 = (bv.maxs.y - origin.y) * n_inv.y;

    let tmin = tmin.max(t1.min(t2).min(tmax));
    let tmax = tmax.min(t1.max(t2).max(tmin));

    tmax >= tmin
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use test_case::test_case;

    #[test]
    fn line2_basis_x() {
        let line = Line2::new([1.0, 1.0].into(), [0.0, 1.0].into());
        let p = line.at(1.0);
        let t = line.to_iso_from_x().inverse() * p;
        assert_relative_eq!(Point2::new(1.0, 0.0), t, epsilon = 1e-12);
    }

    #[test]
    fn line2_basis_y() {
        let line = Line2::new([1.0, 1.0].into(), [1.0, 0.0].into());
        let p = line.at(1.0);
        let t = line.to_iso_from_y().inverse() * p;
        assert_relative_eq!(Point2::new(0.0, 1.0), t, epsilon = 1e-12);
    }

    /// These tests check that the intersection parameter calculation between two parameterized
    /// lines works as expected. The test cases were generated by starting with a random
    /// intersection point and selecting random orthogonal vector and parameter values that rounded
    /// to a single decimal place.
    #[test_case((11.0, 0.7, -4.2, -2.7), (-0.1, -4.7, 1.8, 0.0), (2.0, 1.5))]
    #[test_case((6.2, 3.0, 1.3, -1.9), (-1.3, 8.3, 3.1, -1.7), (-1.0, 2.0))]
    #[test_case((3.6, -3.3, -4.1, 3.2), (9.9, 0.9, 3.0, 2.0), (0.0, -2.1))]
    #[test_case((7.4, -10.1, -2.4, 3.2), (-7.0, -7.7, 4.2, 2.8), (2.5, 2.0))]
    #[test_case((-1.6, 4.3, -0.0, -2.3), (-9.8, 7.1, -4.1, 3.7), (2.0, -2.0))]
    #[test_case((1.7, -0.5, -0.5, -3.3), (3.9, -2.3, 0.5, -1.5), (-1.0, -3.4))]
    #[test_case((4.7, 0.4, -0.6, 0.9), (4.7, 1.3, -3.0, -0.0), (1.0, 0.2))]
    #[test_case((-4.8, -2.0, -1.1, -1.3), (-0.8, -21.0, -0.4, 4.8), (-2.0, 4.5))]
    #[test_case((9.1, 5.7, 4.9, 1.6), (-9.1, -15.5, 2.1, 4.5), (-2.0, 4.0))]
    #[test_case((2.8, 15.7, 0.7, 3.0), (-3.9, -7.1, 1.3, 3.6), (-4.0, 3.0))]
    #[test_case((-5.0, 6.4, -2.6, 0.6), (5.3, 2.3, -2.4, 4.0), (-3.5, 0.5))]
    #[test_case((0.4, -1.9, 2.5, 1.5), (7.6, -19.0, -3.2, 4.2), (-1.6, 3.5))]
    #[test_case((10.6, 5.9, 2.0, 0.7), (-2.0, 7.1, 2.6, -4.7), (-5.0, 1.0))]
    #[test_case((7.3, -0.1, -3.5, 0.9), (14.4, 14.6, -4.4, -3.0), (3.0, 4.0))]
    #[test_case((-5.3, 11.5, -0.8, 4.6), (-15.7, 10.0, -3.1, 2.5), (-2.5, -4.0))]
    #[test_case((4.8, 1.9, 1.0, 1.0), (4.8, 1.9, -1.1, 1.0), (0.0, 0.0))]
    fn inter_param_success(av: (f64, f64, f64, f64), bv: (f64, f64, f64, f64), p: (f64, f64)) {
        let a = Point2::new(av.0, av.1);
        let an = Vector2::new(av.2, av.3);
        let b = Point2::new(bv.0, bv.1);
        let bn = Vector2::new(bv.2, bv.3);

        let (ap, bp) = intersection_param(&a, &an, &b, &bn).unwrap();

        assert_relative_eq!(p.0, ap, epsilon = 1.0e-6);
        assert_relative_eq!(p.1, bp, epsilon = 1.0e-6);
    }

    /// These tests check that the intersection parameter calculation between two parameterized
    /// lines which are parallel returns a None
    #[test_case((-5.0, 2.8, 2.2, 1.8), (-4.2, -0.2, 6.6, 5.4))]
    #[test_case((3.3, 2.5, 4.0, 1.0), (3.2, 0.7, -20.0, -5.0))]
    #[test_case((4.2, -2.3, -0.6, 1.4), (-1.0, 0.5, -2.4, 5.6))]
    #[test_case((-1.1, 2.0, 5.0, 4.0), (4.9, -2.8, 19.5, 15.6))]
    #[test_case((2.4, -3.0, -1.8, -2.6), (0.1, 0.7, 7.2, 10.4))]
    #[test_case((1.2, 2.1, 4.3, -1.0), (-1.4, 3.9, 8.6, -2.0))]
    #[test_case((-4.8, -2.0, -0.1, -0.5), (3.0, -3.6, 0.4, 2.0))]
    #[test_case((-4.4, -0.4, 3.1, 1.1), (3.1, 4.9, -3.1, -1.1))]
    #[test_case((-1.0, -0.1, 0.1, 1.0), (1.2, -2.8, 0.3, 3.0))]
    #[test_case((4.7, -3.7, 4.0, -2.5), (2.5, -0.4, -11.2, 7.0))]
    #[test_case((-1.2, 0.4, 3.9, 0.9), (2.7, -4.6, -11.7, -2.7))]
    #[test_case((-4.8, 4.1, 4.4, -3.0), (4.3, 3.2, -6.6, 4.5))]
    #[test_case((-3.7, 4.7, -1.4, 2.6), (1.0, -4.5, -2.1, 3.9))]
    #[test_case((-0.1, 2.6, 1.4, -0.3), (-0.6, -1.5, -4.2, 0.9))]
    #[test_case((1.8, -2.0, 4.5, -2.0), (0.5, -4.7, 18.0, -8.0))]
    fn inter_parallel_fail(av: (f64, f64, f64, f64), bv: (f64, f64, f64, f64)) {
        let a = Point2::new(av.0, av.1);
        let an = Vector2::new(av.2, av.3);
        let b = Point2::new(bv.0, bv.1);
        let bn = Vector2::new(bv.2, bv.3);

        let result = intersection_param(&a, &an, &b, &bn);

        assert_eq!(None, result);
    }

    // ── Line2 tests ──────────────────────────────────────────────────────────

    fn x_axis_line2() -> Line2 {
        Line2::new(Point2::origin(), Vector2::x())
    }

    #[test]
    fn line2_at_zero_is_origin() {
        let line = x_axis_line2();
        assert_relative_eq!(line.at(0.0), Point2::origin(), epsilon = 1e-12);
    }

    #[test]
    fn line2_at_one_is_origin_plus_direction() {
        let line = Line2::new(Point2::new(1.0, 2.0), Vector2::new(0.0, 1.0));
        assert_relative_eq!(line.at(1.0), Point2::new(1.0, 3.0), epsilon = 1e-12);
        assert_relative_eq!(line.at(-1.0), Point2::new(1.0, 1.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_new_normalize_gives_unit_direction() {
        let line = Line2::new_normalize(Point2::origin(), Vector2::new(3.0, 0.0));
        assert_relative_eq!(line.direction.norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn line2_from_points_direction_is_difference() {
        let line = Line2::from_points(&Point2::new(1.0, 0.0), &Point2::new(4.0, 0.0));
        assert_relative_eq!(line.direction, Vector2::new(3.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_scalar_project_on_line() {
        let line = x_axis_line2();
        assert_relative_eq!(
            line.scalar_project(&Point2::new(5.0, 0.0)),
            5.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn line2_scalar_project_perpendicular_offset() {
        let line = x_axis_line2();
        assert_relative_eq!(
            line.scalar_project(&Point2::new(0.0, 3.0)),
            0.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn line2_closest_point_perpendicular_drop() {
        let line = x_axis_line2();
        let cp = line.closest_point(&Point2::new(4.0, 3.0));
        assert_relative_eq!(cp, Point2::new(4.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_distance_to_known_value() {
        let line = x_axis_line2();
        assert_relative_eq!(
            line.distance_to(&Point2::new(0.0, 3.0)),
            3.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn line2_normal_is_cw_unit_perpendicular() {
        // X-axis: direction = (1, 0), CW normal = (0, -1)
        let line = x_axis_line2();
        assert_relative_eq!(
            line.normal().into_inner(),
            Vector2::new(0.0, -1.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn line2_signed_distance_right_is_positive() {
        // X-axis, point below (right when traveling +x) → positive
        let line = x_axis_line2();
        assert_relative_eq!(
            line.signed_distance_to(&Point2::new(0.0, -3.0)),
            3.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn line2_signed_distance_left_is_negative() {
        let line = x_axis_line2();
        assert_relative_eq!(
            line.signed_distance_to(&Point2::new(0.0, 3.0)),
            -3.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn line2_intersect_perpendicular_lines() {
        let a = Line2::new(Point2::new(0.0, 0.0), Vector2::x());
        let b = Line2::new(Point2::new(3.0, 0.0), Vector2::y());
        let pt = a.intersect(&b).unwrap();
        assert_relative_eq!(pt, Point2::new(3.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_intersect_parallel_returns_none() {
        let a = Line2::new(Point2::new(0.0, 0.0), Vector2::x());
        let b = Line2::new(Point2::new(0.0, 1.0), Vector2::x());
        assert!(a.intersect(&b).is_none());
    }

    #[test]
    fn line2_from_surface_point_roundtrip() {
        use crate::geom2::SurfacePoint2;
        let sp = SurfacePoint2::new_normalize(Point2::new(1.0, 2.0), Vector2::new(0.0, 1.0));
        let line = Line2::from(&sp);
        assert_relative_eq!(line.origin(), sp.point, epsilon = 1e-12);
        assert_relative_eq!(line.direction, sp.normal.into_inner(), epsilon = 1e-12);
    }

    #[test]
    fn line2_to_surface_point_normal_aligns_with_direction() {
        let line = Line2::new_normalize(Point2::new(1.0, 2.0), Vector2::new(1.0, 0.0));
        let sp = SurfacePoint2::from(line);
        assert_relative_eq!(sp.point, line.origin(), epsilon = 1e-12);
        assert_relative_eq!(sp.normal.into_inner(), line.direction, epsilon = 1e-12);
    }

    #[test]
    fn line2_transform_by_preserves_points_on_line() {
        use crate::geom2::tests::Random2;
        let mut rg = Random2::new();
        for _ in 0..200 {
            let line = Line2::new(rg.point(10.0), {
                let a = rg.angle_sym_pi();
                Vector2::new(a.cos(), a.sin())
            });
            let iso = Iso2::new(rg.point(10.0).coords, rg.angle_sym_pi());
            let transformed = line.new_transformed_by(&iso);
            for t in [-2.0, 0.0, 1.0, 3.0] {
                assert_relative_eq!(iso * line.at(t), transformed.at(t), epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn line2_iso_mul_operator() {
        let line = Line2::new(Point2::new(1.0, 0.0), Vector2::x());
        let iso = Iso2::new(Vector2::new(0.0, 5.0), 0.0);
        let t1 = iso * line.clone();
        let t2 = &iso * &line;
        assert_relative_eq!(t1.origin(), Point2::new(1.0, 5.0), epsilon = 1e-12);
        assert_relative_eq!(t2.origin(), Point2::new(1.0, 5.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_aabb_intersection_horizontal_line_through_box() {
        let line = Line2::new(Point2::new(-2.0, 0.0), Vector2::x());
        let aabb = Aabb2::new(Point2::new(-1.0, -1.0), Point2::new(1.0, 1.0));

        let mut ts = line.intersection(aabb);
        ts.sort_by(|a, b| a.partial_cmp(b).unwrap());

        assert_eq!(ts.len(), 2);
        assert_relative_eq!(ts[0], 1.0, epsilon = 1e-12);
        assert_relative_eq!(ts[1], 3.0, epsilon = 1e-12);
        assert_relative_eq!(line.at(ts[0]), Point2::new(-1.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(line.at(ts[1]), Point2::new(1.0, 0.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_aabb_intersection_vertical_line_through_box() {
        let line = Line2::new(Point2::new(0.0, -2.0), Vector2::y());
        let aabb = Aabb2::new(Point2::new(-1.0, -1.0), Point2::new(1.0, 1.0));

        let ts = line.intersection(aabb);

        assert_eq!(ts.len(), 2);
        assert_relative_eq!(ts[0], 1.0, epsilon = 1e-12);
        assert_relative_eq!(ts[1], 3.0, epsilon = 1e-12);
        assert_relative_eq!(line.at(ts[0]), Point2::new(0.0, -1.0), epsilon = 1e-12);
        assert_relative_eq!(line.at(ts[1]), Point2::new(0.0, 1.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_aabb_intersection_corner_returns_single_parameter() {
        let line = Line2::new(Point2::new(-1.0, 1.0), Vector2::new(1.0, 1.0));
        let aabb = Aabb2::new(Point2::new(-1.0, -1.0), Point2::new(1.0, 1.0));

        let ts = line.intersection(aabb);

        assert_eq!(ts.len(), 1);
        assert_relative_eq!(ts[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(line.at(ts[0]), Point2::new(-1.0, 1.0), epsilon = 1e-12);
    }

    #[test]
    fn line2_aabb_intersection_miss_returns_empty() {
        let line = Line2::new(Point2::new(-2.0, 2.0), Vector2::x());
        let aabb = Aabb2::new(Point2::new(-1.0, -1.0), Point2::new(1.0, 1.0));

        let ts = line.intersection(aabb);

        assert!(ts.is_empty());
    }
}
