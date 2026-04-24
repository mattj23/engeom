use crate::{Iso2, Line2, Point2, Vector2};
use crate::common::{Intersection, PCoords};
use crate::common::vec_f64::sort_and_dedup;
use crate::geom2::{intersection_param, signed_angle, Aabb2, Ray2};

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
        let a = other.mins.clone();
        let b = Point2::new(other.maxs.x, other.mins.y);
        let c = other.maxs.clone();
        let d = Point2::new(other.mins.x, other.maxs.y);

        let x_min = Line2::from_points(&a, &d);
        let x_max = Line2::from_points(&b, &c);
        let y_min = Line2::from_points(&a, &b);
        let y_max = Line2::from_points(&c, &d);

        let mut result = Vec::with_capacity(2);

        if let Some((t0, t1)) = self.intersection_params(&x_min) {
            if 0.0 <= t1 && t1 <= 1.0 {
                result.push(t0);
            }
        }
        if let Some((t0, t1)) = self.intersection_params(&x_max) {
            if 0.0 <= t1 && t1 <= 1.0 {
                result.push(t0);
            }
        }
        if let Some((t0, t1)) = self.intersection_params(&y_min) {
            if 0.0 <= t1 && t1 <= 1.0 {
                result.push(t0);
            }
        }
        if let Some((t0, t1)) = self.intersection_params(&y_max) {
            if 0.0 <= t1 && t1 <= 1.0 {
                result.push(t0);
            }
        }

        sort_and_dedup(&mut result, Some(1e-12));

        result
    }
}
