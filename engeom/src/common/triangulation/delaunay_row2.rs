//! This module has an implementation for computing the 2D delaunay triangulation between a set
//! of two parallel strips.  It's a specialized reduced algorithm that's useful for certain
//! meshing cases, typically when a triangulation problem can be transformed into another space.
//!
use crate::common::kd_tree::KdTreeSearch;
use crate::{KdTree2, Point2, Result};

pub fn build_delaunay_strip(
    a: &[DelaunayRowPoint],
    y_a: f64,
    b: &[DelaunayRowPoint],
    y_b: f64,
) -> Result<Vec<[usize; 3]>> {
    // Initialize the indices as the first points in each row. The x positions of the points need
    // to be sorted.
    let mut state = DelaunayRowState::try_new(a, y_a, b, y_b)?;
    let mut faces = Vec::new();

    // The algorithm works by maintaining two indices, one for the first row and one for the
    // second row.  At any point we will look forward for the next (up to) two possible valid
    // delaunay faces that can be constructed by advancing one of the indices.
    //
    // A valid delaunay face is one that does not circumscribe any other point in the two rows, and
    // that doesn't protrude beyond 2x the spacing between rows in either direction.
    //
    // If there are two valid faces, we'll take the one that has the smaller maximum edge length.
    // If there are no valid faces, we'll advance the index of the row smaller current X value.
    while !state.is_done() {
        let next_0 = state.next_0();
        let next_1 = state.next_1();

        match (next_0, next_1) {
            (Some((e0, f0)), Some((e1, f1))) => {
                // Pick the better of the two
                if e0 < e1 {
                    faces.push(f0);
                    state.advance_0();
                } else {
                    faces.push(f1);
                    state.advance_1();
                }
            }
            (Some((_, f0)), None) => {
                // Only the first row has a valid face
                faces.push(f0);
                state.advance_0();
            }
            (None, Some((_, f1))) => {
                // Only the second row has a valid face
                faces.push(f1);
                state.advance_1();
            }

            (None, None) => {
                // If both sides are none and either of the sides are at the end of their points,
                // we can exit now
                if !state.has_next_1() || !state.has_next_0() {
                    break;
                }

                // No valid faces, advance the index of the row with the smaller X value
                if state.point0().x < state.point1().x {
                    state.advance_0();
                } else {
                    state.advance_1();
                }
            }
        }
    }

    Ok(faces)
}

struct DelaunayRowState<'a> {
    /// The current index in the first row.
    i0: usize,

    /// The current index in the second row.
    i1: usize,

    /// The first row of points.
    row0: &'a [DelaunayRowPoint],

    /// The second row of points.
    row1: &'a [DelaunayRowPoint],

    /// The Y coordinate of the first row.
    y0: f64,

    /// The Y coordinate of the second row.
    y1: f64,

    points0: Vec<Point2>,
    points1: Vec<Point2>,

    /// A KD tree for fast point lookups.
    tree: KdTree2,
}

impl<'a> DelaunayRowState<'a> {
    pub fn point0(&self) -> Point2 {
        self.points0[self.i0]
    }

    pub fn point1(&self) -> Point2 {
        self.points1[self.i1]
    }

    fn validate_points(&self, pa: &Point2, pb: &Point2, pc: &Point2) -> Option<f64> {
        // Get the circle
        // let Ok(circle) = Circle2::from_3_points(pa, pb, pc) else {
        //     // If the circle cannot be formed, we cannot form a face.
        //     return None;
        // };
        //
        // // Check for any neighbors
        // let neighbors = self.tree.within(&circle.center, circle.r() + 1e-6);
        // if neighbors.len() > 3 {
        //     // If there are more than 3 neighbors, we cannot form a face.
        //     return None;
        // }
        //
        // let span = self.y1 - self.y0;
        // if circle.y() + circle.r() > self.y1 + 4.0 * span
        //     || circle.y() - circle.r() < self.y0 - 4.0 * span
        // {
        //     // If the circle is too high or too low, we cannot form a face.
        //     return None;
        // }
        // // Return the maximum edge length
        // Some((pb - pa).norm().max((pc - pb).norm()).max((pa - pc).norm()))

        let width = (pb.y - pa.y).abs();
        let max_edge = (pb - pa).norm().max((pc - pb).norm()).max((pa - pc).norm());
        if max_edge < 3.0 * width {
            Some(max_edge)
        } else {
            None
        }
    }

    /// Get the next possible face by advancing the index in row 0
    pub fn next_0(&self) -> Option<(f64, [usize; 3])> {
        if !self.has_next_0() {
            return None;
        }

        // Get the 2D points
        let pa = &self.points0[self.i0];
        let pb = &self.points1[self.i1];
        let pc = &self.points0[self.i0 + 1];

        let Some(max_edge) = self.validate_points(pa, pb, pc) else {
            // If the points do not form a valid face, return None.
            return None;
        };

        // Get the delaunay points
        let a = &self.row0[self.i0];
        let b = &self.row1[self.i1];
        let c = &self.row0[self.i0 + 1];

        Some((max_edge, [a.i, b.i, c.i]))
    }

    pub fn next_1(&self) -> Option<(f64, [usize; 3])> {
        if !self.has_next_1() {
            return None;
        }

        // Get the 2D points
        let pa = &self.points0[self.i0];
        let pb = &self.points1[self.i1];
        let pc = &self.points1[self.i1 + 1];

        let Some(max_edge) = self.validate_points(pa, pb, pc) else {
            // If the points do not form a valid face, return None.
            return None;
        };

        // Get the delaunay points
        let a = &self.row0[self.i0];
        let b = &self.row1[self.i1];
        let c = &self.row1[self.i1 + 1];

        Some((max_edge, [a.i, b.i, c.i]))
    }

    pub fn is_done(&self) -> bool {
        self.i0 == self.row0.len() - 1 && self.i1 == self.row1.len() - 1
    }

    /// Determine if there are more points in the first row beyond the current index.
    pub fn has_next_0(&self) -> bool {
        self.i0 < self.row0.len() - 1
    }

    /// Determine if there are more points in the second row beyond the current index.
    pub fn has_next_1(&self) -> bool {
        self.i1 < self.row1.len() - 1
    }

    pub fn advance_0(&mut self) {
        if self.has_next_0() {
            self.i0 += 1;
        }
    }

    pub fn advance_1(&mut self) {
        if self.has_next_1() {
            self.i1 += 1;
        }
    }

    pub fn try_new(
        row0: &'a [DelaunayRowPoint],
        y0: f64,
        row1: &'a [DelaunayRowPoint],
        y1: f64,
    ) -> Result<Self> {
        // Ensure the rows are sorted by X coordinate.
        if row0.is_empty() || row1.is_empty() {
            return Err("Rows must not be empty".into());
        }
        if !row0.windows(2).all(|w| w[0].x <= w[1].x) {
            return Err("Row 0 must be sorted by X coordinate".into());
        }
        if !row1.windows(2).all(|w| w[0].x <= w[1].x) {
            return Err("Row 1 must be sorted by X coordinate".into());
        }

        let points0 = row0
            .iter()
            .map(|p| Point2::new(p.x, y0))
            .collect::<Vec<_>>();
        let points1 = row1
            .iter()
            .map(|p| Point2::new(p.x, y1))
            .collect::<Vec<_>>();

        let mut all_points = Vec::with_capacity(row0.len() + row1.len());
        for point in points0.iter() {
            all_points.push(*point);
        }
        for point in points1.iter() {
            all_points.push(*point);
        }

        let tree = KdTree2::new(&all_points);

        Ok(Self {
            i0: 0,
            i1: 0,
            row0,
            row1,
            y0,
            y1,
            points0,
            points1,
            tree,
        })
    }
}

pub struct DelaunayRowPoint {
    /// The X coordinate of the point in the row.
    pub x: f64,

    /// The index of the point in the original (potentially untransformed) vertex set.
    pub i: usize,
}

impl DelaunayRowPoint {
    /// Creates a new `DelaunayRowPoint`.
    ///
    /// # Arguments
    ///
    /// * `x`: The X coordinate of the point in the row.
    /// * `i`: The index of the point in the original vertex set.
    pub fn new(x: f64, i: usize) -> Self {
        Self { x, i }
    }
}
