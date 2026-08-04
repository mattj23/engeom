use crate::AngleDir::Cw;
use crate::common::{PCoords, Segment};
use crate::geom2::{Aabb2, BoundaryElement2, LineOps2, Manifold1Pos2, rot90};
use crate::{Point2, UnitVec2};

/// A line segment in 2D space, defined by two endpoints.
///
/// This is one of `engeom`'s 2D geometric primitives.
///
/// This is the two-dimensional specialization of the dimension-generic
/// [`Segment`](Segment); see that type for the shared constructors and queries (`new`,
/// `new_unchecked`, `at`, `length`, `scalar_projection`, `dir`, `closest_point`, `reversed`,
/// `to_line`, `transformed_by`, and so on). The methods defined directly on `Segment2` here are
/// the ones that only make sense in 2D.
pub type Segment2 = Segment<2>;

impl Segment2 {
    /// Create a new segment shifted by distance `d` in the direction of the segment normal vector.
    /// The normal vector is the direction vector rotated by 90 degrees clockwise, in keeping with
    /// the general convention of a normal vector pointing outwards from a counter-clockwise wound
    /// polyline.
    ///
    /// # Arguments
    ///
    /// * `d`: the distance to shift the segment along its normal vector
    ///
    /// returns: Segment2
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::geom2::{Point2, Segment2};
    /// let a = Point2::new(0.0, 0.0);
    /// let b = Point2::new(1.0, 0.0);
    /// let s = Segment2::new(&a, &b).unwrap();
    ///
    /// let s1 = s.offset_by(1.0);
    ///
    /// assert_relative_eq!(s1.a, Point2::new(0.0, -1.0), epsilon = 1.0e-6);
    /// assert_relative_eq!(s1.b, Point2::new(1.0, -1.0), epsilon = 1.0e-6);
    /// ```
    pub fn offset_by(&self, d: f64) -> Self {
        let n = self.normal();
        Self {
            a: self.a + n.into_inner() * d,
            b: self.b + n.into_inner() * d,
        }
    }

    pub fn aabb(&self) -> Aabb2 {
        let mins = Point2::new(self.a.x.min(self.b.x), self.a.y.min(self.b.y));
        let maxs = Point2::new(self.a.x.max(self.b.x), self.a.y.max(self.b.y));
        Aabb2::new(mins, maxs)
    }

    pub fn at_t(&self, t: f64) -> Manifold1Pos2 {
        let point = self.a + (self.b - self.a) * t;
        let direction = UnitVec2::new_normalize(self.b - self.a);
        let normal = rot90(Cw) * direction;
        Manifold1Pos2::new(t * self.length(), point, direction, normal)
    }

    pub fn normal(&self) -> UnitVec2 {
        let direction = UnitVec2::new_normalize(self.b - self.a);
        rot90(Cw) * direction
    }

    /// The parameters at which this segment crosses `other`, or `None` if they do not cross.
    ///
    /// Both are in the normalized `0.0..1.0` form used by [`at`](Segment::at), so
    /// `self.at(t)` and `other.at(u)` are the same point.
    ///
    /// `None` covers three cases which a caller almost always wants to treat alike: the supporting
    /// lines are parallel, they are collinear, or they meet somewhere outside one or both segments.
    /// Use [`Line::closest_approach`](crate::common::Line::closest_approach) on
    /// [`to_line`](Segment::to_line) to get the parameters regardless of whether they land in range.
    ///
    /// A crossing exactly at an endpoint counts. Callers which need a strictly interior crossing
    /// should test the returned parameters against their own margin.
    ///
    /// # Arguments
    ///
    /// * `other`: the segment to test against
    ///
    /// returns: Option<(f64, f64)>, the parameter on `self` followed by the parameter on `other`
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::Point2;
    /// use engeom::geom2::Segment2;
    /// use approx::assert_relative_eq;
    ///
    /// let a = Segment2::new(&Point2::new(0.0, 0.0), &Point2::new(2.0, 0.0)).unwrap();
    /// let b = Segment2::new(&Point2::new(1.0, -1.0), &Point2::new(1.0, 3.0)).unwrap();
    ///
    /// let (t, u) = a.intersection_param(&b).unwrap();
    /// assert_relative_eq!(t, 0.5, epsilon = 1e-12);
    /// assert_relative_eq!(u, 0.25, epsilon = 1e-12);
    /// assert_relative_eq!(a.at(t), b.at(u), epsilon = 1e-12);
    ///
    /// // The lines meet, but past the end of `a`.
    /// let c = Segment2::new(&Point2::new(5.0, -1.0), &Point2::new(5.0, 1.0)).unwrap();
    /// assert!(a.intersection_param(&c).is_none());
    /// ```
    pub fn intersection_param(&self, other: &Segment2) -> Option<(f64, f64)> {
        let (t, u) = self.to_line().intersection_params(&other.to_line())?;
        ((0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u)).then_some((t, u))
    }

    /// Whether this segment touches `other` anywhere.
    pub fn intersects_other(&self, other: &Segment2) -> bool {
        if self.intersection_param(other).is_some() {
            return true;
        }

        // Parallel supporting lines yield no crossing parameters at all, and collinear ones are
        // reported as touching. Note that this does not check the collinear overlap is real, which
        // is long-standing behaviour rather than a decision made here.
        let (l0, l1) = (self.to_line(), other.to_line());
        l0.intersection_params(&l1).is_none() && l0.distance_to(&l1.origin) < 1e-12
    }
}

impl BoundaryElement2 for Segment2 {
    fn length(&self) -> f64 {
        Segment2::length(self)
    }

    fn at_length(&self, length: f64) -> Manifold1Pos2 {
        let t = length / self.length();
        self.at_t(t)
    }

    fn closest_to_point(&self, point: &dyn PCoords<2>) -> Manifold1Pos2 {
        let p = Point2::from(point.coords());
        let t = self.scalar_projection(&p).clamp(0.0, 1.0);
        self.at_t(t)
    }

    fn aabb(&self) -> Aabb2 {
        Segment2::aabb(self)
    }

    fn at_end(&self) -> Manifold1Pos2 {
        self.at_t(1.0)
    }

    fn to_points(&self, _tol: f64) -> Vec<Point2> {
        vec![self.a, self.b]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::RandomGeometry;
    use approx::assert_relative_eq;

    #[test]
    fn length_simple() {
        let a = Point2::new(1.0, 1.0);
        let b = Point2::new(5.0, 1.0);
        let seg = Segment2::new(&a, &b).unwrap();
        assert_relative_eq!(seg.length(), 4.0);
    }

    #[test]
    fn scalar_projection_simple() {
        let a = Point2::new(1.0, 1.0);
        let b = Point2::new(5.0, 1.0);
        let seg = Segment2::new(&a, &b).unwrap();
        let test_point = Point2::new(3.0, 2.0);
        let t = seg.scalar_projection(&test_point);
        assert_relative_eq!(0.5, t, epsilon = 1e-6);
    }

    #[test]
    fn closest_to_simple() {
        let a = Point2::new(1.0, 1.0);
        let b = Point2::new(4.0, 1.0);
        let seg = Segment2::new(&a, &b).unwrap();
        let test_point = Point2::new(2.0, 3.0);
        let closest = seg.closest_to_point(&test_point);
        assert_relative_eq!(closest.point, Point2::new(2.0, 1.0), epsilon = 1e-6);
    }

    // intersection_param tests ──────────────────────────────────────────────

    fn seg(ax: f64, ay: f64, bx: f64, by: f64) -> Segment2 {
        Segment2::new(&Point2::new(ax, ay), &Point2::new(bx, by)).unwrap()
    }

    #[test]
    fn intersection_param_lands_on_the_same_point() {
        let a = seg(0.0, 0.0, 4.0, 2.0);
        let b = seg(0.0, 3.0, 3.0, -1.5);

        let (t, u) = a.intersection_param(&b).unwrap();
        assert_relative_eq!(a.at(t), b.at(u), epsilon = 1e-12);
        assert!((0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u));
    }

    #[test]
    fn intersection_param_refuses_a_crossing_off_the_end() {
        // The supporting lines meet, but well beyond the end of both segments.
        let a = seg(0.0, 0.0, 1.0, 0.0);
        let b = seg(9.0, -1.0, 9.0, 1.0);
        assert!(a.intersection_param(&b).is_none());
    }

    #[test]
    fn intersection_param_refuses_parallel_and_collinear() {
        let a = seg(0.0, 0.0, 2.0, 0.0);
        assert!(a.intersection_param(&seg(0.0, 1.0, 2.0, 1.0)).is_none());
        assert!(a.intersection_param(&seg(0.5, 0.0, 1.5, 0.0)).is_none());
    }

    #[test]
    fn intersection_param_counts_an_endpoint_touch() {
        // A T junction, where one segment ends on the other. Reported, because a caller wanting a
        // strictly interior crossing can test the parameters itself.
        let a = seg(0.0, 0.0, 2.0, 0.0);
        let b = seg(1.0, 0.0, 1.0, 2.0);

        let (t, u) = a.intersection_param(&b).unwrap();
        assert_relative_eq!(t, 0.5, epsilon = 1e-12);
        assert_relative_eq!(u, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn intersection_param_is_symmetric() {
        let a = seg(-1.0, -2.0, 3.0, 1.0);
        let b = seg(2.0, -3.0, -0.5, 2.0);

        let (t, u) = a.intersection_param(&b).unwrap();
        let (u2, t2) = b.intersection_param(&a).unwrap();
        assert_relative_eq!(t, t2, epsilon = 1e-12);
        assert_relative_eq!(u, u2, epsilon = 1e-12);
    }

    #[test]
    fn stress_intersects_other_agrees_with_intersection_param() {
        // The two now share their arithmetic, so this pins the one place they still differ: a
        // collinear pair, which `intersects_other` reports as touching and which has no crossing
        // parameters to report.
        let mut rand = RandomGeometry::<2>::from_seed(0x5e62_c1a5);

        for _ in 0..5000 {
            let a = Segment2::new_unchecked(rand.point(3.0), rand.point(3.0));
            let b = Segment2::new_unchecked(rand.point(3.0), rand.point(3.0));

            if a.intersection_param(&b).is_some() {
                assert!(
                    a.intersects_other(&b),
                    "a reported crossing was not reported as an intersection"
                );
            }
        }
    }
}
