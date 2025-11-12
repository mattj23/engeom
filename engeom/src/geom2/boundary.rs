use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom2::{Aabb2, Line2, ManifoldPosition2};
use crate::{Arc2, SurfacePoint2};

pub trait BoundaryElement {
    /// The total length of the element's manifold domain. For example, for a line segment this
    /// would be the distance from the start to the end point. For an arc it would be the total
    /// arc length.
    fn length(&self) -> f64;

    ///
    ///
    /// # Arguments
    ///
    /// * `length`:
    ///
    /// returns: ManifoldPosition2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn at_length(&self, length: f64) -> ManifoldPosition2;

    fn closest_to_point(&self, point: &impl PCoords<2>) -> ManifoldPosition2;

    fn aabb(&self) -> Aabb2;

    fn at_start(&self) -> ManifoldPosition2 {
        self.at_length(0.0)
    }

    fn at_end(&self) -> ManifoldPosition2 {
        self.at_length(self.length())
    }
}
