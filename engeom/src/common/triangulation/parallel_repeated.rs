use crate::Result;
use crate::common::PCoords;
use crate::na::Point;

/// This is a builder struct which manages the triangulation of repeated rows of vertices where the
/// triangulation pattern should be a regular join of each vertex to its neighbor with the same
/// index.
pub struct ParallelBuilder<const D: usize> {
    /// The number of points in each row
    n: usize,

    rows: usize,

    /// The builder which handles registering vertices
    vertices: Vec<Point<f64, D>>,

    /// The face index
    faces: Vec<[u32; 3]>,

    /// Set to true if you want the last row to join to the first row
    close: bool,
}

impl<const D: usize> ParallelBuilder<D> {
    pub fn new(n: usize, close: bool) -> Self {
        Self {
            n,
            rows: 0,
            vertices: Vec::new(),
            faces: Vec::new(),
            close,
        }
    }

    pub fn push(&mut self, points: &[impl PCoords<D>]) -> Result<()> {
        if points.len() != self.n {
            return Err(format!(
                "Expected {} points in the row, but got {}",
                self.n,
                points.len()
            )
            .into());
        }

        for p in points {
            self.vertices.push(Point::from(p.coords()))
        }

        self.rows += 1;

        if self.rows > 1 {
            for i in 0..(self.n - 1) {
                let b = i + (self.rows - 1) * self.n;
                let a = b - self.n;
                let c = b + 1;
                let d = c - self.n;
                self.faces.push([a as u32, b as u32, c as u32]);
                self.faces.push([a as u32, c as u32, d as u32]);
            }
        }

        Ok(())
    }

    pub fn take(mut self) -> (Vec<Point<f64, D>>, Vec<[u32; 3]>) {
        if self.close {
            for i in 0..(self.n - 1) {
                let b = i;
                let a = i + (self.rows - 1) * self.n;
                let c = b + 1;
                let d = a + 1;
                self.faces.push([a as u32, b as u32, c as u32]);
                self.faces.push([a as u32, c as u32, d as u32]);
            }
        }

        (self.vertices, self.faces)
    }
}
