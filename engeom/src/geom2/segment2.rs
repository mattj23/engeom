use crate::AngleDir::Cw;
use crate::common::{PCoords, Segment};
use crate::geom2::{Aabb2, BoundaryElement2, LineOps2, Manifold1Pos2, rot90};
use crate::{Iso2, Point2, UnitVec2};

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

    pub fn intersects_other(&self, other: &Segment2) -> bool {
        let l0 = self.to_line();
        let l1 = other.to_line();
        let Some((t0, t1)) = l0.intersection_params(&l1) else {
            return l0.distance_to(&l1.origin) < 1e-12;
        };

        t0 >= 0.0 && t1 >= 0.0 && t0 <= 1.0 && t1 <= 1.0
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
}
