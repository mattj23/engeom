//! In `engeom`, 2D boundaries are a concept that represents a continuous manifold in 2D space,
//! with a positive and a negative side. Boundaries may be open or closed (currently, closed
//! boundaries are not implemented, but many of their features can be achieved with an open boundary
//! that returns to its starting point).
//!
//! A boundary is very similar to a [`Curve2`] element, except that while a `Curve2` is composed
//! entirely of line segments, a boundary may consist of different types of elements, all of which
//! must implement the [`BoundaryElement`] trait. Currently, only [`Segment2`] and [`Arc2`] are
//! implemented.
//!
//! In comparison to a `Curve2`, a boundary can represent curved geometry with near-theoretical
//! precision, and without the use of approximation through large amounts of points.  When
//! representing geometry with curved regions, boundaries will generally be more efficient and
//! capable.
//!
//! Boundaries are defined through the [`BoundaryData2`] struct, which provides a convenient and
//! efficient way to define their geometry while ensuring that the continuity constraint isn't
//! violated.  There are helper methods on `BoundaryData2` which allow for constructing complicated
//! geometry, alongside simple methods for full control over segments and arcs.
//!
//! Actual queryable geometry is represented by the [`Boundary2`] struct, which contains a vector
//! of boxed dynamic `BoundaryElement` instances. The `Boundary2` is most easily built from the
//! [`BoundaryData2::try_to_boundary`] method.  Once created, spatial queries and measurements can
//! be performed.
//!
//! Another important and useful feature of the boundary concept is the ability to perform
//! fitting of boundary geometry to measured/observed geometry, like points. A common use for this
//! ability is to try to extract semantic or design meaning from observed data. For example,
//! fitting a slot or square hole to a cross-section in order to determine which points should be
//! used to measure width. A generalized implementation of fitting to points is provided in
//! [`fit_boundary_to_points`].

mod construction;
mod fitting;
mod data;

use crate::common::PCoords;
use crate::common::points::dist;
use crate::geom2::{Aabb2, ManifoldPosition2};
use crate::geom2::{Arc2, Segment2};
use crate::{Point2, Result};
use parry2d_f64::bounding_volume::BoundingVolume;

pub use data::BoundaryData2;
pub use construction::BCursor;
pub use fitting::{BndBuildFn, fit_boundary_to_points};

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

    fn to_points(&self, tol: f64) -> Vec<Point2>;
}

/// Contains the geometry of a boundary, which is a collection of elements that can be queried for
///
pub struct Boundary2 {
    elements: Vec<Box<dyn BoundaryElement>>,
    lengths: Vec<f64>,
}

impl Boundary2 {
    pub fn try_new(elements: Vec<Box<dyn BoundaryElement>>) -> Result<Self> {
        if elements.is_empty() {
            return Err("Boundary must have at least one element".into());
        }

        let mut lengths = vec![0.0];
        let mut total_length = 0.0;
        for element in elements.iter() {
            total_length += element.length();
            lengths.push(total_length);
        }

        Ok(Self { elements, lengths })
    }

    pub fn to_points(&self, tol: f64) -> Result<Vec<Point2>> {
        let mut points = vec![self.elements[0].at_start().point];
        for e in self.elements.iter() {
            let e_points = e.to_points(tol);
            points.extend(e_points.into_iter().skip(1));
        }

        Ok(points)
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

    pub fn at_start(&self) -> ManifoldPosition2 {
        self.elements[0].at_start()
    }

    pub fn at_end(&self) -> ManifoldPosition2 {
        let mut result = self.elements[self.elements.len() - 1].at_end();
        result.l = self.length();
        result
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
        let mut data = BoundaryData2::new_open(Point2::new(0.0, 0.0));
        let mut cursor = data.get_cursor(None);
        cursor.add_seg_xy(1.0, 0.0);
        cursor.add_arc_xy(1.0, 0.5, 1.0, 1.0, false);
        cursor.add_seg_xy(0.0, 1.0);
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
