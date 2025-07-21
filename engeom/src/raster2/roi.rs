use crate::raster2::{Point2I, Vector2I};

/// This is a small lightweight struct (128 bits in size) that represents the correspondence between
/// a point in a region of interest (ROI) and its parent point in the original image. The two
/// different coordinates refer to the same pixel, but one is in the indices of the ROI and the
/// other in indices of the full sized image.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RoiPoint {
    /// The local point has coordinates in the ROI itself, typically ranging  between (0, 0) and
    /// (width, height) of the ROI
    pub local: Point2I,

    /// The parent point has coordinates in the original image, typically ranging between
    /// (min_x, min_y) and (max_x, max_y) of the original ROI.  This is the point in the parent
    /// image that corresponds the `local` point in the ROI.
    pub parent: Point2I,
}

#[derive(Debug, Clone)]
pub struct RoiOverlap {
    a: RasterRoi,
    b: RasterRoi,
    i: RasterRoi,
}

impl RoiOverlap {
    pub fn a(&self) -> &RasterRoi {
        &self.a
    }

    pub fn b(&self) -> &RasterRoi {
        &self.b
    }

    pub fn new(a: RasterRoi, b: RasterRoi) -> Self {
        let i = a.intersection(&b);
        Self { a, b, i }
    }

    pub fn iter_intersection_a(&self) -> RoiIterator {
        let offset = self.i.min - self.a.min;
        RoiIterator::new(&self.i, offset)
    }

    pub fn iter_intersection_b(&self) -> RoiIterator {
        let offset = self.i.min - self.b.min;
        RoiIterator::new(&self.i, offset)
    }
}

/// A struct representing a rectangular region of interest (ROI) in a raster image. The semantics
/// of this structure are similar to a bounding box, except that the maximum corner is exclusive
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
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
        x >= self.min.x && x < self.max.x && y >= self.min.y && y < self.max.y
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

    pub fn union(&self, other: &RasterRoi) -> Self {
        if self.is_empty() {
            return *other;
        }
        if other.is_empty() {
            return *self;
        }
        Self {
            min: Point2I::new(self.min.x.min(other.min.x), self.min.y.min(other.min.y)),
            max: Point2I::new(self.max.x.max(other.max.x), self.max.y.max(other.max.y)),
        }
    }

    pub fn intersects(&self, other: &RasterRoi) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }
        self.min.x < other.max.x
            && self.max.x > other.min.x
            && self.min.y < other.max.y
            && self.max.y > other.min.y
    }

    pub fn intersection(&self, other: &RasterRoi) -> Self {
        if !self.intersects(other) {
            return Self::empty();
        }
        Self {
            min: Point2I::new(self.min.x.max(other.min.x), self.min.y.max(other.min.y)),
            max: Point2I::new(self.max.x.min(other.max.x), self.max.y.min(other.max.y)),
        }
    }

    pub fn iter_points(&self) -> RoiIterator {
        RoiIterator::new(self, Vector2I::default())
    }

    pub fn expanded(&self, padding: u32) -> RasterRoi {
        let padding = padding as i32;
        let offset = Vector2I::new(padding, padding);
        let new_min = self.min - offset;
        let new_max = self.max + offset;
        RasterRoi::new(new_min, new_max)
    }
}

impl Default for RasterRoi {
    fn default() -> Self {
        Self::empty()
    }
}

pub struct RoiIterator<'a> {
    pub roi: &'a RasterRoi,
    pub x: i32,
    pub y: i32,
    pub offset: Vector2I,
}

impl<'a> RoiIterator<'a> {
    pub fn new(roi: &'a RasterRoi, offset: Vector2I) -> Self {
        Self {
            roi,
            x: 0,
            y: 0,
            offset,
        }
    }
}

impl<'a> Iterator for RoiIterator<'a> {
    type Item = RoiPoint;

    fn next(&mut self) -> Option<Self::Item> {
        let extent = self.roi.extent();

        if self.y >= extent.y {
            return None;
        }

        let roi = Point2I::new(self.x, self.y);
        let parent = self.roi.in_to_out(roi);

        let point = RoiPoint {
            local: roi + self.offset,
            parent,
        };

        self.x += 1;

        if self.x >= extent.x {
            self.x = 0;
            self.y += 1;
        }

        Some(point)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::image::{ImageBuffer, Luma};
    use crate::raster2::RasterMask;
    use crate::raster2::roi_mask::RoiMask;
    use imageproc::definitions::Image;
    use std::collections::HashMap;

    fn make_test_setup() -> (HashMap<Point2I, usize>, Image<Luma<usize>>, RasterRoi) {
        let mut map = HashMap::new();
        let mut image = ImageBuffer::new(60, 20);
        let mut c = 1;
        let mut roi = RasterRoi::empty();
        let corner = Vector2I::new(12, 8);

        for x in 0..24 {
            for y in 0..7 {
                let inside = Point2I::new(x, y);
                map.insert(inside, c);

                let outside = inside + corner;
                roi.expand_to_contain(outside);

                image.put_pixel(outside.x as u32, outside.y as u32, Luma([c]));

                c += 1;
            }
        }
        (map, image, roi)
    }

    #[test]
    fn roi_iterator_forward() {
        let (map, image, roi) = make_test_setup();

        for p in roi.iter_points() {
            // ONLY the points in the ROI should be iterated through, and they are the only
            // ones that should be in the map
            let expected = map.get(&p.local);
            assert!(expected.is_some(), "Point {:?} not found in map", p.local);
            let expected = expected.unwrap();

            // The pixel in the image should match the value in the map
            let actual = image.get_pixel(p.parent.x as u32, p.parent.y as u32)[0];
            assert_eq!(*expected, actual, "Mismatch at point {:?}", p);
        }
    }

    #[test]
    fn roi_mask_offset_matching() {
        let (map, image, roi) = make_test_setup();

        // We'll create a new expanded roi mask that is 2 pixels larger in each direction than
        let mask = RasterMask::empty(roi.extent().x as u32 + 4, roi.extent().y as u32 + 4);
        let offset = Vector2I::new(-2, -2);
        let roi_mask = RoiMask::new_resized(mask, roi, roi.expanded(2));

        // We'll create a manual mapping of where the values should be in the expanded mask
        let mut expected_map = HashMap::new();
        for (kp, v) in map.iter() {
            // The coordinates in the expanded mask are shifted by the inverse of the offset vector
            // so for instance, (2, 2) in the original roi would be (4, 4) in the expanded mask
            expected_map.insert(kp - offset, *v);
        }

        // Check that I did that right
        assert_eq!(*expected_map.get(&Point2I::new(2, 2)).unwrap(), 1);
        assert!(!expected_map.contains_key(&Point2I::new(0, 0)));

        for p in roi_mask.iter_shared_points() {
            let expected = expected_map.get(&p.local);
            assert!(expected.is_some(), "Point {:?} not found in map", p.local);
            let expected = expected.unwrap();

            // The pixel in the image should match the value in the map
            let actual = image.get_pixel(p.parent.x as u32, p.parent.y as u32)[0];
            assert_eq!(*expected, actual, "Mismatch at point {:?}", p);
        }
    }

    #[test]
    fn roi_expand() {
        let mut roi = RasterRoi::empty();
        assert!(roi.is_empty());
        roi.expand_to_contain(Point2I::new(1, 1));
        assert!(!roi.is_empty());
        assert_eq!(roi.min, Point2I::new(1, 1));
        assert_eq!(roi.max, Point2I::new(2, 2));

        roi.expand_to_contain(Point2I::new(0, 0));
        assert_eq!(roi.min, Point2I::new(0, 0));
        assert_eq!(roi.max, Point2I::new(2, 2));

        roi.expand_to_contain(Point2I::new(3, 3));
        assert_eq!(roi.min, Point2I::new(0, 0));
        assert_eq!(roi.max, Point2I::new(4, 4));
    }

    #[test]
    fn test_raster_roi_empty() {
        let roi = RasterRoi::empty();
        assert!(roi.is_empty());
    }
}
