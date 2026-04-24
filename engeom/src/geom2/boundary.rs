use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom2::{Aabb2, ManifoldPosition2};
use crate::geom2::{Arc2, Segment2};
use crate::{Point2, Result};
use parry2d_f64::bounding_volume::BoundingVolume;

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

    fn closest_to_point(&self, point: &dyn PCoords<2>) -> ManifoldPosition2;

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
    Arc((f64, f64, f64, f64, bool)),
}

#[derive(Debug, Clone)]
pub struct BoundaryData2 {
    pub start: Point2,
    elements: Vec<BKind>,
}

impl BoundaryData2 {
    pub fn new(start: Point2) -> Self {
        Self {
            start,
            elements: Vec::new(),
        }
    }

    pub fn add_seg_xy(&mut self, x: f64, y: f64) {
        self.elements.push(BKind::Seg((x, y)));
    }

    pub fn add_seg(&mut self, p: &impl PCoords<2>) {
        self.add_seg_xy(p.coords().x, p.coords().y);
    }

    pub fn add_arc_xy(
        &mut self,
        center_x: f64,
        center_y: f64,
        end_x: f64,
        end_y: f64,
        clockwise: bool,
    ) {
        self.elements
            .push(BKind::Arc((center_x, center_y, end_x, end_y, clockwise)));
    }

    pub fn add_arc(&mut self, center: &impl PCoords<2>, end: &impl PCoords<2>, clockwise: bool) {
        self.add_arc_xy(
            center.coords().x,
            center.coords().y,
            end.coords().x,
            end.coords().y,
            clockwise,
        );
    }

    pub fn try_to_boundary(&self) -> Result<Boundary2> {
        let mut elements: Vec<Box<dyn BoundaryElement>> = Vec::new();
        let mut last_point = self.start;

        for e in self.elements.iter() {
            match e {
                BKind::Seg((x, y)) => {
                    let end = Point2::new(*x, *y);
                    let seg = Segment2::try_new(&last_point, &end)?;
                    last_point = end;
                    elements.push(Box::new(seg));
                }
                BKind::Arc((cx, cy, ex, ey, cw)) => {
                    let end = Point2::new(*ex, *ey);
                    let center = Point2::new(*cx, *cy);
                    let arc = Arc2::try_new_ends(&last_point, &end, &center, *cw)?;
                    last_point = end;
                    elements.push(Box::new(arc));
                }
            }
        }

        Ok(Boundary2::new(elements))
    }
}

/// Contains the geometry of a boundary, which is a collection of elements that can be queried for
///
pub struct Boundary2 {
    elements: Vec<Box<dyn BoundaryElement>>,
    lengths: Vec<f64>,
}

impl Boundary2 {
    pub fn new(elements: Vec<Box<dyn BoundaryElement>>) -> Self {
        let mut lengths = vec![0.0];
        let mut total_length = 0.0;
        for element in elements.iter() {
            total_length += element.length();
            lengths.push(total_length);
        }

        Self { elements, lengths }
    }

    pub fn at_closest_to_point(&self, point: &dyn PCoords<2>) -> ManifoldPosition2 {
        let point = Point2::from(point.coords());
        let (k, _, m) = (0..self.elements.len())
            .map(|i| {
                let m = self.elements[i].closest_to_point(&point);
                let d = dist(&point, &m);

                (i, d, m)
            })
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();

        let length = self.lengths[k] + m.l;
        ManifoldPosition2::new(length, m.point, m.direction, m.normal)
    }

    pub fn at_length(&self, length: f64) -> Option<ManifoldPosition2> {
        let pre_mod = if length < 0.0 || length > self.length() {
            None
        } else {
            let search = self
                .lengths
                .binary_search_by(|a| a.partial_cmp(&length).unwrap());
            match search {
                Ok(i) => {
                    if i == 0 {
                        Some(self.elements[0].at_start())
                    } else {
                        Some(self.elements[i - 1].at_end())
                    }
                }
                Err(next_i) => {
                    let i = next_i - 1;
                    let rem_l = length - self.lengths[i];
                    Some(self.elements[i].at_length(rem_l))
                }
            }
        };

        pre_mod.map(|pr| ManifoldPosition2 { l: length, ..pr })
    }

    pub fn length(&self) -> f64 {
        *self.lengths.last().unwrap()
    }

    pub fn aabb(&self) -> Aabb2 {
        let mut aabb = self.elements[0].aabb().clone();
        for element in self.elements.iter().skip(1) {
            aabb.merge(&element.aabb())
        }
        aabb
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector2;
    use approx::assert_relative_eq;
    use faer::rand::Rng;
    use rand::RngExt;
    use std::f64::consts::PI;

    fn simple_data() -> BoundaryData2 {
        let mut data = BoundaryData2::new(Point2::new(0.0, 0.0));
        data.add_seg_xy(1.0, 0.0);
        data.add_arc_xy(1.0, 0.5, 1.0, 1.0, false);
        data.add_seg_xy(0.0, 1.0);
        data
    }

    #[test]
    fn simple_boundary_builds() {
        let data = simple_data();
        let boundary = data.try_to_boundary().unwrap();
        let aabb = boundary.aabb();
        assert_eq!(boundary.elements.len(), 3);

        assert_relative_eq!(aabb.mins.x, 0.0);
        assert_relative_eq!(aabb.mins.y, 0.0);
        assert_relative_eq!(aabb.maxs.x, 1.5);
        assert_relative_eq!(aabb.maxs.y, 1.0);
    }

    #[test]
    fn simple_boundary_length() {
        let data = simple_data();
        let boundary = data.try_to_boundary().unwrap();
        let expected = 2.0 + PI / 2.0;

        assert_relative_eq!(boundary.length(), expected, epsilon = 1e-12);
    }

    #[test]
    fn simple_boundary_at_length_start() {
        let data = simple_data();
        let boundary = data.try_to_boundary().unwrap();
        let m = boundary.at_length(0.0).unwrap();
        assert_relative_eq!(Point2::origin(), m.point, epsilon = 1e-12);
        assert_relative_eq!(Vector2::x(), m.direction, epsilon = 1e-12)
    }

    #[test]
    fn simple_boundary_at_length_middle() {
        let data = simple_data();
        let boundary = data.try_to_boundary().unwrap();
        let m = boundary.at_length(boundary.length() / 2.0).unwrap();
        assert_relative_eq!(Point2::new(1.5, 0.5), m.point, epsilon = 1e-12);
        assert_relative_eq!(Vector2::y(), m.direction, epsilon = 1e-12)
    }

    #[test]
    fn simple_boundary_at_length_end() {
        let data = simple_data();
        let boundary = data.try_to_boundary().unwrap();
        let m = boundary.at_length(boundary.length()).unwrap();
        assert_relative_eq!(Point2::new(0.0, 1.0), m.point, epsilon = 1e-12);
        assert_relative_eq!(-Vector2::x(), m.direction, epsilon = 1e-12)
    }

    #[test]
    fn simple_boundary_at_segment_end() {
        let data = simple_data();
        let boundary = data.try_to_boundary().unwrap();
        // let m = boundary.at_length(boundary.elements[0].length()).unwrap();
        let m = boundary.at_length(1.0).unwrap();
        assert_relative_eq!(Point2::new(1.0, 0.0), m.point, epsilon = 1e-12);
        assert_relative_eq!(Vector2::x(), m.direction, epsilon = 1e-12)
    }

    #[test]
    fn stress_boundary_length_closest_round_trip() {
        let data = simple_data();
        let boundary = data.try_to_boundary().unwrap();
        let mut rng = rand::rng();
        for _ in 0..1000 {
            let l = rng.random_range(0.0..boundary.length());
            let m = boundary.at_length(l).unwrap();
            let l2 = boundary.at_closest_to_point(&m.point).l;
            assert_relative_eq!(l, l2, epsilon = 1e-6);
        }
    }
}
