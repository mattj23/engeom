use crate::Result;
use crate::image::{GenericImage, GrayImage, ImageFormat, ImageReader, Luma};
use crate::raster2::index_iter::IndexIter;
use crate::raster2::{LabeledRegions, Point2I, zhang_suen_thinning};
use imageproc::distance_transform::Norm;
use imageproc::drawing::{
    draw_filled_circle_mut, draw_filled_rect_mut, draw_hollow_circle_mut, draw_hollow_rect_mut,
    draw_polygon_mut,
};
use imageproc::morphology::{dilate_mut, erode_mut};
use imageproc::rect::Rect;
use imageproc::region_labelling::Connectivity;
use parry2d_f64::utils::hashmap::HashMap;
use std::io::BufWriter;
use std::path::Path;

type IpPoint = imageproc::point::Point<i32>;

#[derive(Clone, Debug)]
pub struct RasterMask {
    pub buffer: GrayImage,
}

impl RasterMask {
    /// Create a new `RasterMask` by taking ownership of a `GrayImage`.
    pub fn new(buffer: GrayImage) -> RasterMask {
        RasterMask { buffer }
    }

    pub fn empty(width: u32, height: u32) -> RasterMask {
        let buffer = GrayImage::new(width, height);
        RasterMask { buffer }
    }

    pub fn empty_like(example: &impl GenericImage) -> RasterMask {
        let buffer = GrayImage::new(example.width(), example.height());
        RasterMask { buffer }
    }

    pub fn save(&self, path: &Path, fmt: ImageFormat) -> Result<()> {
        let file = std::fs::File::create(path)?;
        let mut writer = BufWriter::new(file);

        self.buffer
            .write_to(&mut writer, fmt)
            .map_err(|e| format!("Failed to save mask to {}: {}", path.display(), e).into())
    }

    pub fn load(path: &Path) -> Result<RasterMask> {
        let loaded = ImageReader::open(path)?.with_guessed_format()?.decode()?;
        let buffer = loaded.into_luma8();
        Ok(RasterMask { buffer })
    }

    /// Sets a point/pixel in the mask to a specified value. If the point is out of bounds, it will
    /// return an error.
    ///
    /// # Arguments
    ///
    /// * `p`: the point to set, represented as a `Point2I`.
    /// * `value`: the value to set the point to, represented as a `bool`. If true, the pixel will
    ///   be set to white (255), and if false, it will be set to black (0).
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    pub fn set_point(&mut self, p: Point2I, value: bool) -> Result<()> {
        if self.point_in_bounds(p) {
            self.set_point_unchecked(p, value);
            Ok(())
        } else {
            Err(format!(
                "Point ({}, {}) is out of bounds for mask of size {}x{}",
                p.x,
                p.y,
                self.width(),
                self.height()
            )
            .into())
        }
    }

    /// Sets a point in the mask to a specified value without checking if the point is within
    /// bounds. If the point has an index below 0 it will silently write to the pixel at the 0
    /// index, and if the point has an index above the width or height of the mask, it will
    /// panic with an out-of-bounds error.  Only use this method if you are sure the point is
    /// within bounds of the mask.
    ///
    /// # Arguments
    ///
    /// * `p`: the point to set, represented as a `Point2I`.
    /// * `value`: the value to set the point to, represented as a `bool`. If true, the pixel will
    ///   be set to white (255), and if false, it will be set to black (0).
    ///
    /// returns: ()
    pub fn set_point_unchecked(&mut self, p: Point2I, value: bool) {
        self.buffer
            .put_pixel(p.x as u32, p.y as u32, Luma([value as u8 * 255]))
    }

    /// Sets a point/pixel in the mask to a specified value if the point is within the bounds of the
    /// mask, otherwise does nothing. This is a convenience method to simplify the logic of
    /// common algorithms in which setting a pixel off the mask should just be ignored.
    ///
    /// If the point is outside the mask bounds, the function will also return false, in case that
    /// information is needed by the calling code.
    ///
    /// # Arguments
    ///
    /// * `p`: the point to set, represented as a `Point2I`.
    /// * `value`: the value to set the point to, represented as a `bool`. If true, the pixel will
    ///   be set to white (255), and if false, it will be set to black (0).
    ///
    /// returns: bool
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::raster2::{RasterMask, Point2I};
    /// let mut mask = RasterMask::empty(10, 10);
    /// let p0 = Point2I::new(5, 5);
    ///
    /// assert!(mask.set_point_if_in_bounds(p0, true));
    /// assert!(mask.get_point(p0));
    ///
    /// let p1 = Point2I::new(-5, -5);
    /// assert!(!mask.set_point_if_in_bounds(p1, true));
    /// assert!(!mask.get_point(p1));
    /// ```
    pub fn set_point_if_in_bounds(&mut self, p: Point2I, value: bool) -> bool {
        if self.point_in_bounds(p) {
            self.set_point_unchecked(p, value);
            true
        } else {
            false
        }
    }

    // pub fn get(&self, x: u32, y: u32) -> bool {
    //     self.buffer.get_pixel(x, y)[0] > 0
    // }

    /// Get the value of the mask at the specified pixel/point. If the point is out of the bounds
    /// of the mask, it is by default considered to be false.
    ///
    /// # Arguments
    ///
    /// * `p`: the point to check, represented as a `Point2I`.
    ///
    /// returns: bool
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn get_point(&self, p: Point2I) -> bool {
        self.point_in_bounds(p) && self.buffer.get_pixel(p.x as u32, p.y as u32)[0] == 255
    }

    pub fn width(&self) -> u32 {
        self.buffer.width()
    }

    pub fn height(&self) -> u32 {
        self.buffer.height()
    }

    // ==========================================================================================
    // Index Operations
    // ==========================================================================================

    /// Returns true if a point is within the bounds of the mask. This requires the x and y
    /// coordinates to be non-negative and less than the width and height of the mask, respectively.
    ///
    /// # Arguments
    ///
    /// * `p`: the point to check, represented as a `Point2I`.
    ///
    /// returns: bool
    pub fn point_in_bounds(&self, p: Point2I) -> bool {
        p.x >= 0 && p.x < self.width() as i32 && p.y >= 0 && p.y < self.height() as i32
    }

    /// Computes the neighbors of a point in the mask based on the specified connectivity. The
    /// neighbors are returned as a vector of `Point2I` objects. The connectivity can be either
    /// four-connected (4) or eight-connected (8). The function filters out any neighbors that
    /// are out of bounds of the mask, so the returned vector will have a length between 0 and 8.
    ///
    /// # Arguments
    ///
    /// * `p`: the point for which to find neighbors
    /// * `connectivity`: the type of connectivity to use, either `Connectivity::Four` or
    ///   `Connectivity::Eight`
    ///
    /// returns: Vec<OPoint<i32, Const<2>>, Global>
    pub fn get_neighbors(&self, p: Point2I, connectivity: Connectivity) -> Vec<Point2I> {
        let mut candidates = vec![
            Point2I::new(p.x - 1, p.y),
            Point2I::new(p.x + 1, p.y),
            Point2I::new(p.x, p.y - 1),
            Point2I::new(p.x, p.y + 1),
        ];

        if connectivity == Connectivity::Eight {
            candidates.extend(vec![
                Point2I::new(p.x - 1, p.y - 1),
                Point2I::new(p.x + 1, p.y - 1),
                Point2I::new(p.x - 1, p.y + 1),
                Point2I::new(p.x + 1, p.y + 1),
            ]);
        }

        // Filter out candidates that are out of bounds
        candidates
            .into_iter()
            .filter(|&np| self.point_in_bounds(np))
            .collect()
    }

    pub fn iter_all(&self) -> impl Iterator<Item = Point2I> {
        IndexIter::new(self.width() as usize, self.height() as usize)
    }

    // ==========================================================================================
    // Truth operations
    // ==========================================================================================

    pub fn iter_true(&self) -> RasterMaskTrueIterator {
        RasterMaskTrueIterator {
            mask: self,
            x: 0,
            y: 0,
        }
    }

    pub fn count_true(&self) -> usize {
        self.buffer.as_raw().iter().filter(|&&v| v > 0).count()
    }

    // ==========================================================================================
    // NOT Operations
    // ==========================================================================================
    /// Invert the mask in place, i.e., set all true values to false and vice versa.
    pub fn not_mut(&mut self) {
        for v in self.buffer.iter_mut() {
            *v = !*v;
        }
    }

    /// Create a new mask that is the inverse of the current mask.
    pub fn not(&self) -> RasterMask {
        let mut new_mask = self.clone();
        new_mask.not_mut();
        new_mask
    }

    // ==========================================================================================
    // OR (Union) Operations
    // ==========================================================================================

    /// Set the value at (x, y) to true if either of the masks has a true value at that position.
    pub fn or_mut(&mut self, other: &RasterMask) -> Result<()> {
        if self.width() != other.width() || self.height() != other.height() {
            return Err("Masks must have the same dimensions".into());
        }
        for (va, vb) in self.buffer.iter_mut().zip(other.buffer.iter()) {
            *va |= *vb;
        }

        Ok(())
    }

    /// Create a new mask that is the union of the current mask and another mask.
    pub fn or(&self, other: &RasterMask) -> Result<RasterMask> {
        if self.width() != other.width() || self.height() != other.height() {
            return Err("Masks must have the same dimensions".into());
        }
        let mut new_mask = self.clone();
        new_mask.or_mut(other)?;
        Ok(new_mask)
    }

    // ==========================================================================================
    // AND (Intersection) Operations
    // ==========================================================================================
    /// Set the value at (x, y) to true if both masks have a true value at that position.
    /// This is the same as the logical AND of the two masks, or a set intersection.
    pub fn and_mut(&mut self, other: &RasterMask) -> Result<()> {
        if self.width() != other.width() || self.height() != other.height() {
            return Err("Masks must have the same dimensions".into());
        }
        for (va, vb) in self.buffer.iter_mut().zip(other.buffer.iter()) {
            *va &= *vb;
        }

        Ok(())
    }

    /// Create a new mask that is the intersection of the current mask and another mask.
    /// This is the same as the logical AND of the two masks, or a set intersection.
    pub fn and(&self, other: &RasterMask) -> Result<RasterMask> {
        if self.width() != other.width() || self.height() != other.height() {
            return Err("Masks must have the same dimensions".into());
        }
        let mut new_mask = self.clone();
        new_mask.and_mut(other)?;
        Ok(new_mask)
    }

    // ==========================================================================================
    // AND NOT (Difference) Operations
    // =========================================================================================
    /// Set the value at (x, y) to false if both masks have a true value at that position. This is
    /// the same as a logical NOT operation on the second mask followed by a logical AND with the
    /// first mask, which is the equivalent to a set difference operation, or subtracting items
    /// from the second mask from the first mask.
    pub fn and_not_mut(&mut self, other: &RasterMask) -> Result<()> {
        if self.width() != other.width() || self.height() != other.height() {
            return Err("Masks must have the same dimensions".into());
        }
        for (va, vb) in self.buffer.iter_mut().zip(other.buffer.iter()) {
            *va &= !*vb;
        }

        Ok(())
    }

    /// Create a new mask that is the logical AND NOT of the current mask and another mask. This
    /// is the same as a logical NOT operation on the second mask followed by a logical AND with
    /// first mask, also the equivalent to a set difference operation, or subtracting items in the
    /// second mask from the first mask.
    pub fn and_not(&self, other: &RasterMask) -> Result<RasterMask> {
        if self.width() != other.width() || self.height() != other.height() {
            return Err("Masks must have the same dimensions".into());
        }
        let mut new_mask = self.clone();
        new_mask.and_not_mut(other)?;
        Ok(new_mask)
    }

    // ==========================================================================================
    // Morphological Operations
    // =========================================================================================
    pub fn erode_mut(&mut self, norm: Norm, k: u8) {
        erode_mut(&mut self.buffer, norm, k);
    }

    pub fn eroded(&self, norm: Norm, k: u8) -> RasterMask {
        let mut new_mask = self.clone();
        new_mask.erode_mut(norm, k);
        new_mask
    }

    pub fn dilate_mut(&mut self, norm: Norm, k: u8) {
        dilate_mut(&mut self.buffer, norm, k);
    }

    pub fn dilated(&self, norm: Norm, k: u8) -> RasterMask {
        let mut new_mask = self.clone();
        new_mask.dilate_mut(norm, k);
        new_mask
    }

    pub fn zhang_suen_thin(&mut self) {
        zhang_suen_thinning(self);
    }

    pub fn erode_alternating_norms_mut(&mut self, count: usize) {
        for i in 0..count {
            if i % 2 == 0 {
                erode_mut(&mut self.buffer, Norm::L1, 1);
            } else {
                erode_mut(&mut self.buffer, Norm::LInf, 1);
            }
        }
    }

    pub fn eroded_alternating_norms(&self, count: usize) -> RasterMask {
        let mut new_mask = self.clone();
        new_mask.erode_alternating_norms_mut(count);
        new_mask
    }

    pub fn dilate_alternating_norms_mut(&mut self, count: usize) {
        for i in 0..count {
            if i % 2 == 0 {
                dilate_mut(&mut self.buffer, Norm::L1, 1);
            } else {
                dilate_mut(&mut self.buffer, Norm::LInf, 1);
            }
        }
    }

    pub fn dilated_alternating_norms(&self, count: usize) -> RasterMask {
        let mut new_mask = self.clone();
        new_mask.dilate_alternating_norms_mut(count);
        new_mask
    }

    /// Performs a flood fill operation starting from a set of points. The flood fill will expand
    /// to all false pixels that are connected to the starting points, based on the specified
    /// connectivity (4 or 8).
    ///
    /// The function will return a new `RasterMask` where only the pixels that were filled in the
    /// operation are set to true, and all other pixels are false. The starting points will be
    /// included in the filled area.
    ///
    /// If a starting point is out of bounds but has valid neighbors, the flood fill will begin
    /// from those neighbors. If a starting point is `true` in the original mask, it will be
    /// discarded, and it will be used as a starting point, even if it has valid neighbors which
    /// are `false`.
    ///
    /// # Arguments
    ///
    /// * `points`: a slice of `Point2I` representing the starting points for the flood fill.
    /// * `connectivity`: the type of connectivity to use for the flood fill, either
    ///   `Connectivity::Four` or `Connectivity::Eight`.
    ///
    /// returns: RasterMask
    pub fn get_flood_fill_from_points(
        &self,
        points: &[Point2I],
        connectivity: Connectivity,
    ) -> RasterMask {
        // We'll check the original set of points given to us to verify that they are all set to
        // false in the original mask, and if they are not we'll remove them now. We will not check
        // if the points are in bounds, because we will allow them to be out of bounds in case they
        // have valid neighbors which are in bounds.
        let mut stack = points
            .iter()
            .filter_map(|p| if !self.get_point(*p) { Some(*p) } else { None })
            .collect::<Vec<_>>();
        let mut output = RasterMask::empty_like(&self.buffer);

        // Because we have prefiltered the stack to only contain points which are false in the
        // original mask, and because we don't care if the points are in bounds or not, we can
        // safely pop each point from the stack and set it to true in the output mask. When
        // filtering the neighbors, we must check that they are false in the original mask before
        // putting them into the stack in order to prevent this assumption from being violated
        // as the process continues.
        while let Some(p) = stack.pop() {
            // Set the working point in the output mask to true
            output.set_point_unchecked(p, true);

            // Now iterate through all neighbors of the working point and push any which are both
            // false in the original mask and the output mask onto the stack. We are guaranteed
            // from the `get_neighbors` function that all neighbors will be in bounds, and that
            // any point already set to true in the output mask has already been checked for
            // neighbors (including the point adjacent to this working point).
            for n in self.get_neighbors(p, connectivity) {
                if !self.get_point(n) && !output.get_point(n) {
                    stack.push(n);
                }
            }
        }

        output
    }

    /// Gets the result of a flood fill operation starting from the borders of the mask. The flood
    /// fill will expand to all false pixels that are connected to the borders, based on the
    /// specified connectivity (4 or 8).
    ///
    /// This is a convenient way to extract the exterior of a mask, where the exterior is defined
    /// as false pixels that touch the mask border. This is useful for processes that need to
    /// distinguish between interior empty spaces and exterior ones, for things like hole counting
    /// or filling.
    ///
    /// # Arguments
    ///
    /// * `connectivity`: the type of connectivity to use for the flood fill, either
    ///   `Connectivity::Four` or `Connectivity::Eight`.
    ///
    /// returns: RasterMask
    pub fn get_flood_fill_from_borders(&self, connectivity: Connectivity) -> RasterMask {
        let mut border_points = Vec::new();

        for i in 0..self.width() {
            border_points.push(Point2I::new(i as i32, 0));
            border_points.push(Point2I::new(i as i32, (self.height() - 1) as i32));
        }

        for i in 0..self.height() {
            border_points.push(Point2I::new(0, i as i32));
            border_points.push(Point2I::new((self.width() - 1) as i32, i as i32));
        }

        self.get_flood_fill_from_points(&border_points, connectivity)
    }

    // ==========================================================================================
    // Drawing Operations
    // =========================================================================================
    pub fn draw_rect_mut(&mut self, min: Point2I, max: Point2I, value: bool, filled: bool) {
        let size = max - min;
        let r = Rect::at(min.x, min.y).of_size(size.x as u32, size.y as u32);
        let color = if value { Luma([255]) } else { Luma([0]) };
        if filled {
            draw_filled_rect_mut(&mut self.buffer, r, color);
        } else {
            draw_hollow_rect_mut(&mut self.buffer, r, color);
        }
    }

    pub fn draw_circle_mut(&mut self, center: Point2I, radius: i32, value: bool, filled: bool) {
        let color = if value { Luma([255]) } else { Luma([0]) };
        if filled {
            draw_filled_circle_mut(&mut self.buffer, (center.x, center.y), radius, color);
        } else {
            draw_hollow_circle_mut(&mut self.buffer, (center.x, center.y), radius, color);
        }
    }

    pub fn draw_polygon_mut(&mut self, points: &[Point2I], value: bool, filled: bool) -> Result<()> {
        let color = if value { Luma([255]) } else { Luma([0]) };
        if filled {
            let ipoints = polygon_ipoints(points);
            if ipoints.len() < 3 {
                return Err("Cannot draw a polygon with less than 3 points".into());
            }
            draw_polygon_mut(&mut self.buffer, &ipoints, color);
        } else {
            todo!("drawing a hollow polygon is not implemented yet");
        }

        Ok(())
    }

    // ==========================================================================================
    // Misc Convenience Operations
    // =========================================================================================
    pub fn connected_regions(&self, connectivity: Connectivity) -> LabeledRegions {
        LabeledRegions::from_connected(&self.buffer, connectivity, Luma([0]))
    }

    pub fn convex_hull(&self) -> Vec<Point2I> {
        let result = imageproc::geometry::convex_hull(self);
        result.into_iter().map(|p| Point2I::new(p.x, p.y)).collect()
    }

    /// Generates a list of vertices and triangles from the mask, where each true pixel's
    /// coordinates are treated as a vertex, and each triangle is formed by connecting three
    /// vertices.
    ///
    /// This function is a helper in generating meshes from rasters, in that it will generate
    /// the triangular mesh structure from the connectivity of the mask pixels, but the actual
    /// vertices will be calculated through some operation on the pixel index on some other
    /// data structure, such as a ScalarRaster or a UV map.
    pub fn triangle_structure(&self) -> (Vec<Point2I>, Vec<[u32; 3]>) {
        let mut vertices = Vec::new();
        let mut by_index = HashMap::new();
        for p in self.iter_true() {
            vertices.push(p);
            by_index.insert(p, (vertices.len() - 1) as u32);
        }

        let mut faces = Vec::<[u32; 3]>::new();
        for (p, i) in by_index.iter() {
            let p_right = Point2I::new(p.x + 1, p.y);
            let p_up = Point2I::new(p.x, p.y + 1);
            let p_up_right = Point2I::new(p.x + 1, p.y + 1);
            if let Some(i_up_right) = by_index.get(&p_up_right) {
                // If we do have the upper right corner, we will form the two faces as possible
                // going diagonally to the corner
                if let Some(i_right) = by_index.get(&p_right) {
                    faces.push([*i, *i_right, *i_up_right]);
                }

                if let Some(i_up) = by_index.get(&p_up) {
                    faces.push([*i, *i_up_right, *i_up]);
                }
            } else {
                if let (Some(i_right), Some(i_up)) = (by_index.get(&p_right), by_index.get(&p_up)) {
                    // If we do not have the upper right corner, but do have the upper and right
                    // corners independently, we will form the face canted in the other
                    // direction.
                    faces.push([*i, *i_right, *i_up]);
                }
            }

            // Lastly, we'll check if the point below is missing, if so we'll fill in the
            // diagonal face if the right and right lower points exist
            let p_down = Point2I::new(p.x, p.y - 1);
            if by_index.contains_key(&p_down) {
                continue;
            }

            let p_down_right = Point2I::new(p.x + 1, p.y - 1);
            if let (Some(i_down_right), Some(i_right)) =
                (by_index.get(&p_down_right), by_index.get(&p_right))
            {
                faces.push([*i, *i_down_right, *i_right]);
            }
        }
        (vertices, faces)
    }
}

fn polygon_ipoints(points: &[Point2I]) -> Vec<IpPoint> {
    let mut working = points.to_vec();
    while working.len() > 2 && working[0] == working[working.len() - 1] {
        // Remove the last point if it is the same as the first point
        working.pop();
    }

    working.iter().map(|p| IpPoint::new(p.x, p.y)).collect()
}

pub struct RasterMaskTrueIterator<'a> {
    mask: &'a RasterMask,
    x: u32,
    y: u32,
}

impl<'a> RasterMaskTrueIterator<'a> {
    pub fn new(mask: &'a RasterMask) -> Self {
        RasterMaskTrueIterator { mask, x: 0, y: 0 }
    }
}

impl<'a> Iterator for RasterMaskTrueIterator<'a> {
    type Item = Point2I;

    fn next(&mut self) -> Option<Self::Item> {
        while self.y < self.mask.buffer.height() {
            if self.x >= self.mask.buffer.width() {
                self.x = 0;
                self.y += 1;
            }
            if self.y >= self.mask.buffer.height() {
                return None;
            }

            let current = Point2I::new(self.x as i32, self.y as i32);

            if self.mask.get_point(current) {
                self.x += 1;
                return Some(current);
            }
            self.x += 1;
        }
        None
    }
}

impl From<&RasterMask> for Vec<IpPoint> {
    fn from(mask: &RasterMask) -> Vec<IpPoint> {
        let mut points = Vec::new();
        for p in mask.iter_true() {
            points.push(IpPoint::new(p.x, p.y));
        }
        points
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_true() -> Result<()> {
        let mut mask = RasterMask::empty(4, 2);
        mask.set_point(Point2I::new(0, 0), true)?;
        mask.set_point(Point2I::new(1, 1), true)?;
        mask.set_point(Point2I::new(2, 0), true)?;

        assert_eq!(mask.count_true(), 3);
        Ok(())
    }

    #[test]
    fn not_mut() -> Result<()> {
        let mut mask = RasterMask::empty(4, 2);
        mask.set_point(Point2I::new(0, 0), true)?;
        mask.set_point(Point2I::new(1, 1), true)?;
        mask.not_mut();

        assert!(!mask.get_point(Point2I::new(0, 0)));
        assert!(!mask.get_point(Point2I::new(1, 1)));
        assert_eq!(mask.count_true(), 6);

        Ok(())
    }

    #[test]
    fn create_empty() {
        let mask = RasterMask::empty(10, 10);
        assert_eq!(mask.width(), 10);
        assert_eq!(mask.height(), 10);
    }

    #[test]
    fn value_set_get() -> Result<()> {
        let mut mask = RasterMask::empty(5, 5);
        mask.set_point(Point2I::new(2, 2), true)?;
        assert!(mask.get_point(Point2I::new(2, 2)));
        assert!(!mask.get_point(Point2I::new(1, 1)));
        Ok(())
    }

    #[test]
    fn buffer_is_row_major_order() -> Result<()> {
        // Verify that the buffer is in row-major order by checking the first few pixels
        let mut mask = RasterMask::empty(5, 3);
        mask.set_point(Point2I::new(0, 0), true)?;
        mask.set_point(Point2I::new(2, 0), true)?;

        let true_indices = mask
            .buffer
            .as_raw()
            .into_iter()
            .enumerate()
            .filter_map(|(i, &v)| if v > 0 { Some(i) } else { None })
            .collect::<Vec<_>>();

        assert_eq!(true_indices, vec![0, 2]);
        Ok(())
    }

    #[test]
    fn iter_true() -> Result<()> {
        let mut mask = RasterMask::empty(4, 3);

        mask.set_point(Point2I::new(1, 0), true)?;
        mask.set_point(Point2I::new(3, 1), true)?;
        mask.set_point(Point2I::new(2, 2), true)?;

        let mut iter = mask.iter_true();
        assert_eq!(iter.next(), Some(Point2I::new(1, 0)));
        assert_eq!(iter.next(), Some(Point2I::new(3, 1)));
        assert_eq!(iter.next(), Some(Point2I::new(2, 2)));
        assert_eq!(iter.next(), None);

        Ok(())
    }

    #[test]
    fn exterior_flood_fill() {
        let mut mask = RasterMask::empty(20, 10);
        mask.draw_rect_mut(Point2I::new(5, 0), Point2I::new(15, 10), true, true);

        let expected = mask.not();

        mask.draw_rect_mut(Point2I::new(9, 3), Point2I::new(13, 7), false, true);
        let exterior = mask.get_flood_fill_from_borders(Connectivity::Eight);

        for p in exterior.iter_all() {
            assert_eq!(
                exterior.get_point(p),
                expected.get_point(p),
                "Point {:?} does not match expected value",
                p
            );
        }
    }
}
