use crate::raster2::{Point2I, Vector2I};

/// A struct representing a rectangular region of interest (ROI) in a raster image. The semantics
/// of this structure are similar to a bounding box, except that the maximum corner is exclusive
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RasterRoi {
    pub min: Point2I,
    pub max: Point2I,
}

impl RasterRoi {
    pub fn is_empty(&self) -> bool {
        self.min.x == self.max.x || self.min.y == self.max.y
    }

    pub fn empty() -> Self {
        Self {
            min: Point2I::new(0, 0),
            max: Point2I::new(0, 0),
        }
    }

    pub fn new(min: Point2I, max: Point2I) -> Self {
        assert!(min.x <= max.x && min.y <= max.y, "Invalid ROI bounds");
        Self { min, max }
    }

    pub fn from_bounds(min_x: i32, min_y: i32, max_x: i32, max_y: i32) -> Self {
        Self {
            min: Point2I::new(min_x, min_y),
            max: Point2I::new(max_x, max_y),
        }
    }

    pub fn extent(&self) -> Vector2I {
        self.max - self.min
    }

    pub fn contains_indices(&self, x: i32, y: i32) -> bool {
        x >= self.min.x
            && x < self.max.x
            && y >= self.min.y
            && y < self.max.y
    }

    /// Takes an x, y coordinate from the world and returns the x, y coordinate referenced from the
    /// minimum corner of the bounds
    pub fn out_to_in(&self, outside: Point2I) -> Point2I {
        Point2I::from(outside - self.min)
    }

    /// Takes an x, y coordinates referencing the minimum corner of the bounds and returns the
    /// corresponding x, y coordinates from the world
    pub fn in_to_out(&self, inside: Point2I) -> Point2I {
        Point2I::from(inside + self.min.coords)
    }

    pub fn expand_to_contain(&mut self, outside: Point2I) {
        if self.is_empty() {
            self.min = outside;
            self.max = self.min + Vector2I::new(1, 1);
        } else {
            self.min.x = self.min.x.min(outside.x);
            self.min.y = self.min.y.min(outside.y);
            self.max.x = self.max.x.max(outside.x + 1);
            self.max.y = self.max.y.max(outside.y + 1);
        }
    }
}

impl Default for RasterRoi {
    fn default() -> Self {
        Self::empty()
    }
}