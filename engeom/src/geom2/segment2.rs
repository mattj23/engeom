use crate::AngleDir::Cw;
use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom2::{Aabb2, BoundaryElement, Line2, ManifoldPosition2, rot90};
use crate::{Iso2, Point2, Result, TransformBy, UnitVec2, Vector2};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Segment2 {
    pub a: Point2,
    pub b: Point2,
    pub length: f64,
}

impl Segment2 {
    pub fn try_new(a: &impl PCoords<2>, b: &impl PCoords<2>) -> Result<Self> {
        let length = dist(a, b);
        if length < 1e-12 {
            Err("The two points are too close to each other".into())
        } else {
            let a = Point2::from(a.coords());
            let b = Point2::from(b.coords());
            Ok(Self { a, b, length })
        }
    }

    /// Calculate the scalar projection of a set of coordinates onto the line segment, in which
    /// 0.0 represents a point at the segment's starting point `a` and 1.0 represents a point at
    /// the segment's end point `b`.  The result can be any finite value, including negative ones.
    ///
    /// # Arguments
    ///
    /// * `other`: an element with a position in 2d space
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn scalar_projection(&self, other: &impl PCoords<2>) -> f64 {
        let dir = self.b - self.a;
        let test = other.coords() - self.a.coords();
        dir.dot(&test) / self.length.powi(2)
    }

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
    /// let s = Segment2::try_new(&a, &b).unwrap();
    ///
    /// let s1 = s.with_offset(1.0);
    ///
    /// assert_relative_eq!(s1.a, Point2::new(0.0, -1.0), epsilon = 1.0e-6);
    /// assert_relative_eq!(s1.b, Point2::new(1.0, -1.0), epsilon = 1.0e-6);
    /// ```
    pub fn with_offset(&self, d: f64) -> Self {
        let n = UnitVec2::new_normalize(self.orthogonal());
        Self {
            a: self.a + n.into_inner() * d,
            b: self.b + n.into_inner() * d,
            length: self.length,
        }
    }

    /// Create a new segment with the points reversed
    pub fn reversed(&self) -> Self {
        Self {
            a: self.b,
            b: self.a,
            length: self.length,
        }
    }

    pub fn aabb(&self) -> Aabb2 {
        let mins = Point2::new(self.a.x.min(self.b.x), self.a.y.min(self.b.y));
        let maxs = Point2::new(self.a.x.max(self.b.x), self.a.y.max(self.b.y));
        Aabb2::new(mins, maxs)
    }

    pub fn at_t(&self, t: f64) -> ManifoldPosition2 {
        let point = self.a + (self.b - self.a) * t;
        let direction = UnitVec2::new_normalize(self.b - self.a);
        let normal = rot90(Cw) * direction;
        ManifoldPosition2::new(t * self.length, point, direction, normal)
    }
}

impl TransformBy<Iso2, Segment2> for Segment2 {
    fn transform_by(&self, t: &Iso2) -> Self {
        Self {
            a: t.transform_point(&self.a),
            b: t.transform_point(&self.b),
            length: self.length,
        }
    }
}

impl Line2 for Segment2 {
    fn origin(&self) -> Point2 {
        self.a
    }

    fn dir(&self) -> Vector2 {
        self.b - self.a
    }

    fn at(&self, t: f64) -> Point2 {
        self.a + self.dir() * t
    }
}

impl BoundaryElement for Segment2 {
    fn length(&self) -> f64 {
        self.length
    }

    fn at_length(&self, length: f64) -> ManifoldPosition2 {
        let t = length / self.length;
        self.at_t(t)
    }

    fn closest_to_point(&self, point: &impl PCoords<2>) -> ManifoldPosition2 {
        let t = self.scalar_projection(point).clamp(0.0, 1.0);
        self.at_t(t)
    }

    fn aabb(&self) -> Aabb2 {
        Segment2::aabb(self)
    }

    fn at_end(&self) -> ManifoldPosition2 {
        self.at_t(1.0)
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
        let seg = Segment2::try_new(&a, &b).unwrap();
        assert_relative_eq!(seg.length(), 4.0);
    }

    #[test]
    fn scalar_projection_simple() {
        let a = Point2::new(1.0, 1.0);
        let b = Point2::new(5.0, 1.0);
        let seg = Segment2::try_new(&a, &b).unwrap();
        let test_point = Point2::new(3.0, 2.0);
        let t = seg.scalar_projection(&test_point);
        assert_relative_eq!(0.5, t, epsilon = 1e-6);
    }

    #[test]
    fn closest_to_simple() {
        let a = Point2::new(1.0, 1.0);
        let b = Point2::new(4.0, 1.0);
        let seg = Segment2::try_new(&a, &b).unwrap();
        let test_point = Point2::new(2.0, 3.0);
        let closest = seg.closest_to_point(&test_point);
        assert_relative_eq!(closest.point, Point2::new(2.0, 1.0), epsilon = 1e-6);
    }
}
