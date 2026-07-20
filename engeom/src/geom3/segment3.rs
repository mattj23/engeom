use crate::common::{PCoords, Segment};
use crate::geom3::Aabb3;
use crate::{Iso3, Manifold1Pos3, Point3, UnitVec3};

/// A line segment in 3D space, defined by two endpoints.
///
/// This is one of `engeom`'s 3D geometric primitives.
///
/// This is the three-dimensional specialization of the dimension-generic
/// [`Segment`](Segment); see that type for the shared constructors and queries (`new`,
/// `new_unchecked`, `at`, `length`, `scalar_projection`, `dir`, `closest_point`, `reversed`,
/// `to_line`, `transformed_by`, and so on). The methods defined directly on `Segment3` here are
/// the ones that only make sense in 3D.
pub type Segment3 = Segment<3>;

impl Segment3 {
    pub fn aabb(&self) -> Aabb3 {
        let mins = Point3::new(
            self.a.x.min(self.b.x),
            self.a.y.min(self.b.y),
            self.a.z.min(self.b.z),
        );
        let maxs = Point3::new(
            self.a.x.max(self.b.x),
            self.a.y.max(self.b.y),
            self.a.z.max(self.b.z),
        );
        Aabb3::new(mins, maxs)
    }

    pub fn at_t(&self, t: f64) -> Manifold1Pos3 {
        let point = self.a + (self.b - self.a) * t;
        let direction = UnitVec3::new_normalize(self.b - self.a);
        Manifold1Pos3::new(t * self.length(), point, direction)
    }

    /// Returns the manifold position of the point on the segment closest to `point`, clamped to
    /// the segment's endpoints.
    pub fn closest_to_point(&self, point: &impl PCoords<3>) -> Manifold1Pos3 {
        let t = self.scalar_projection(point).clamp(0.0, 1.0);
        self.at_t(t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use approx::assert_relative_eq;

    #[test]
    fn new_rejects_coincident_points() {
        let a = Point3::new(1.0, 1.0, 1.0);
        assert!(Segment3::new(&a, &a).is_err());
    }

    #[test]
    fn length_simple() {
        let a = Point3::new(1.0, 1.0, 1.0);
        let b = Point3::new(5.0, 1.0, 1.0);
        let seg = Segment3::new(&a, &b).unwrap();
        assert_relative_eq!(seg.length(), 4.0);
    }

    #[test]
    fn scalar_projection_simple() {
        let a = Point3::new(1.0, 1.0, 0.0);
        let b = Point3::new(5.0, 1.0, 0.0);
        let seg = Segment3::new(&a, &b).unwrap();
        let test_point = Point3::new(3.0, 2.0, 0.0);
        let t = seg.scalar_projection(&test_point);
        assert_relative_eq!(0.5, t, epsilon = 1e-6);
    }

    #[test]
    fn closest_to_simple() {
        let a = Point3::new(1.0, 1.0, 0.0);
        let b = Point3::new(4.0, 1.0, 0.0);
        let seg = Segment3::new(&a, &b).unwrap();
        let test_point = Point3::new(2.0, 3.0, 5.0);
        let closest = seg.closest_to_point(&test_point);
        assert_relative_eq!(closest.point, Point3::new(2.0, 1.0, 0.0), epsilon = 1e-6);
    }

    #[test]
    fn closest_to_point_clamps_to_endpoints() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let seg = Segment3::new(&a, &b).unwrap();

        let before = seg.closest_to_point(&Point3::new(-5.0, 1.0, 0.0));
        assert_relative_eq!(before.point, a, epsilon = 1e-12);

        let after = seg.closest_to_point(&Point3::new(5.0, 1.0, 0.0));
        assert_relative_eq!(after.point, b, epsilon = 1e-12);
    }

    #[test]
    fn reversed_swaps_endpoints() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 2.0, 3.0);
        let seg = Segment3::new(&a, &b).unwrap();
        let rev = seg.reversed();
        assert_relative_eq!(rev.a, b, epsilon = 1e-12);
        assert_relative_eq!(rev.b, a, epsilon = 1e-12);
        assert_relative_eq!(rev.length(), seg.length(), epsilon = 1e-12);
    }

    #[test]
    fn aabb_matches_bounding_extents() {
        let a = Point3::new(1.0, -1.0, 3.0);
        let b = Point3::new(-2.0, 4.0, 0.0);
        let seg = Segment3::new(&a, &b).unwrap();
        let aabb = seg.aabb();
        assert_relative_eq!(aabb.mins, Point3::new(-2.0, -1.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs, Point3::new(1.0, 4.0, 3.0), epsilon = 1e-12);
    }

    #[test]
    fn at_and_at_t_agree_on_position() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(4.0, 0.0, 0.0);
        let seg = Segment3::new(&a, &b).unwrap();
        for t in [0.0, 0.25, 0.5, 1.0] {
            assert_relative_eq!(seg.at(t), seg.at_t(t).point, epsilon = 1e-12);
        }
    }

    #[test]
    fn at_t_l_matches_arc_length() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(4.0, 0.0, 0.0);
        let seg = Segment3::new(&a, &b).unwrap();
        let pos = seg.at_t(0.5);
        assert_relative_eq!(pos.l, 2.0, epsilon = 1e-12);
        assert_relative_eq!(pos.direction.into_inner(), Vector3::x(), epsilon = 1e-12);
    }

    #[test]
    fn transformed_by_moves_endpoints() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let seg = Segment3::new(&a, &b).unwrap();
        let iso = Iso3::translation(1.0, 2.0, 3.0);
        let moved = seg.transformed_by(&iso);
        assert_relative_eq!(moved.a, Point3::new(1.0, 2.0, 3.0), epsilon = 1e-12);
        assert_relative_eq!(moved.b, Point3::new(2.0, 2.0, 3.0), epsilon = 1e-12);
    }

    #[test]
    fn to_line_passes_through_endpoints() {
        let a = Point3::new(1.0, 2.0, 3.0);
        let b = Point3::new(4.0, 5.0, 6.0);
        let seg = Segment3::new(&a, &b).unwrap();
        let line = seg.to_line();
        assert_relative_eq!(line.distance_to(&a), 0.0, epsilon = 1e-12);
        assert_relative_eq!(line.distance_to(&b), 0.0, epsilon = 1e-12);
    }
}
