//! This module will be the place for tools that compute triangulations

use crate::na::Point;

pub mod parallel_row2;

/// This struct is a simple helper to keep track of vertices in a triangulation. Adding a point
/// returns the index that the point can be found in the list.  Once all points have been added,
/// the `take_points` method can be called to move the points vector out of the builder.
pub struct VertexBuilder<const D: usize> {
    points: Vec<Point<f64, D>>,
}

impl<const D: usize> Default for VertexBuilder<D> {
    fn default() -> Self {
        Self::new()
    }
}

impl<const D: usize> VertexBuilder<D> {
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }

    pub fn points(&self) -> &[Point<f64, D>] {
        &self.points
    }

    pub fn push(&mut self, point: Point<f64, D>) -> usize {
        let index = self.points.len();
        self.points.push(point);
        index
    }

    pub fn take_points(self) -> Vec<Point<f64, D>> {
        self.points
    }
}
