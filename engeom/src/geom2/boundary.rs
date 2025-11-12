use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom2::{Aabb2, Line2, Segment2};
use crate::{Arc2, SurfacePoint2};

pub trait BoundaryElement {
    fn closest_to_point(&self, point: &impl PCoords<2>) -> f64;

    fn aabb(&self) -> Aabb2;

    fn at_start(&self) -> SurfacePoint2;

    fn at_end(&self) -> SurfacePoint2;

    fn at_length(&self, length: f64) -> SurfacePoint2;
}

impl BoundaryElement for Segment2 {
    fn closest_to_point(&self, point: &impl PCoords<2>) -> f64 {
        let proj = self.projected_parameter(point);
        if proj < 0.0 {
            0.0
        } else if proj > 1.0 {
            self.length()
        } else {
            proj * self.length()
        }
    }

    fn aabb(&self) -> Aabb2 {
        Segment2::aabb(self)
    }

    fn at_start(&self) -> SurfacePoint2 {
        SurfacePoint2::new_normalize(self.a, self.orthogonal())
    }

    fn at_end(&self) -> SurfacePoint2 {
        SurfacePoint2::new_normalize(self.b, self.orthogonal())
    }

    fn at_length(&self, length: f64) -> SurfacePoint2 {
        let fraction = length / self.length();
        SurfacePoint2::new_normalize(self.at(fraction), self.orthogonal())
    }
}
