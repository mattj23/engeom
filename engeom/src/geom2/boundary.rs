use crate::common::PCoords;
use crate::geom2::{Aabb2, ManifoldPosition2};
use crate::Point2;
use crate::geom2::{Segment2, Arc2};

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

#[derive(Debug, Clone)]
enum BKind {
    /// A line segment, containing the end point
    Seg((f64, f64)),

    /// An arc, containing the center point, end point, and whether the arc is clockwise
    Arc((f64, f64, f64, f64, bool))
}

#[derive(Debug, Clone)]
pub struct BoundaryData2 {
    pub start: Point2,
    elements: Vec<BKind>,
}

impl BoundaryData2 {
    pub fn new(start: Point2) -> Self {
        Self { start, elements: Vec::new() }
    }

    pub fn add_seg_xy(&mut self, x: f64, y: f64) {
        self.elements.push(BKind::Seg((x, y)));
    }

    pub fn add_seg(&mut self, p: &impl PCoords<2>) {
        self.add_seg_xy(p.coords().x, p.coords().y);
    }

    pub fn add_arc_xy(&mut self, center_x: f64, center_y: f64, end_x: f64, end_y: f64, clockwise: bool) {
        self.elements.push(BKind::Arc((center_x, center_y, end_x, end_y, clockwise)));
    }

    pub fn add_arc(&mut self, center: &impl PCoords<2>, end: &impl PCoords<2>, clockwise: bool) {
        self.add_arc_xy(center.coords().x, center.coords().y, end.coords().x, end.coords().y, clockwise);
    }
}

/// Contains the geometry of a boundary, which is a collection of elements that can be queried for
///
pub struct Boundary2 {

}