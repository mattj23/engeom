//! This module contains an abstraction for mapping between a 2D cartesian space of real numbers
//! and a raster space of pixels.

use crate::image::GrayImage;
use crate::na::{DMatrix, Scalar};
use crate::raster2::{Point2I, RasterMask};
use crate::{Iso2, Point2, Vector2};
use num_traits::Zero;

#[derive(Clone, Copy)]
pub struct RasterMapping {
    origin: Point2,

    /// The shape of the corresponding image as (rows, cols)
    shape: (usize, usize),

    /// The size of each pixel in the physical length units of the cartesian space (e.g. mm,
    /// inches, meters, etc.)
    px_size: f64,

    /// A transform applied to a point in the cartesian space before it is mapped to the image
    forward: Iso2,

    /// The inverse of the forward transform
    inverse: Iso2,
}

impl RasterMapping {
    /// Create a new `nalgebra` dmatrix of the same shape as the raster mapping, filled with zeros.
    pub fn make_zero_matrix<T: Scalar + Zero>(&self) -> DMatrix<T> {
        DMatrix::zeros(self.shape.0, self.shape.1)
    }

    /// Create a new `RasterMask` of the same shape as the raster mapping, filled with zeros.
    /// This is an 8-bit grayscale image that can be used for fast raster masking operations.
    pub fn make_mask(&self) -> RasterMask {
        RasterMask::empty(self.shape.1 as u32, self.shape.0 as u32)
    }

    /// Create a new `ImageMapping` given the origin, shape, pixel size, and an optional transform
    /// which is applied to a point in the cartesian space before it is mapped to the image.
    ///
    /// # Arguments
    ///
    /// * `origin`: The point in the cartesian space which is mapped to the top-left corner of the
    ///   image.
    /// * `shape`: The `(rows, cols)` of the image, starting at the origin. Points in the first
    ///   quadrant (positive y, positive x) will be mapped into the image rows and columns,
    ///   respectively.
    /// * `px_size`: The size of each pixel in the cartesian space units (e.g. mm, inches, meters,
    ///   etc.)
    /// * `transform`: An optional transform which is applied to a points in the cartesian space
    ///   before they are mapped to the image.  If `None` is provided, the identity transform is
    ///   used.
    ///
    /// returns: RasterMapping
    pub fn new(
        origin: Point2,
        shape: (usize, usize),
        px_size: f64,
        transform: Option<Iso2>,
    ) -> Self {
        let forward = transform.unwrap_or(Iso2::identity());
        let inverse = forward.inverse();

        Self {
            origin,
            shape,
            px_size,
            forward,
            inverse,
        }
    }

    pub fn origin(&self) -> Point2 {
        self.origin
    }

    pub fn shape(&self) -> (usize, usize) {
        self.shape
    }

    pub fn px_size(&self) -> f64 {
        self.px_size
    }

    pub fn forward(&self) -> Iso2 {
        self.forward
    }

    pub fn inverse(&self) -> Iso2 {
        self.inverse
    }

    pub fn rows(&self) -> usize {
        self.shape.0
    }

    pub fn cols(&self) -> usize {
        self.shape.1
    }

    /// Return the row and column of the image which is mapped to the given point in the cartesian
    /// space.
    ///
    /// # Arguments
    ///
    /// * `p`:
    ///
    /// returns: (usize, usize)
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn row_col_of(&self, p: &Point2) -> (usize, usize) {
        let p = self.forward * p;
        let row = ((p.y - self.origin.y) / self.px_size) as usize;
        let col = ((p.x - self.origin.x) / self.px_size) as usize;
        (row, col)
    }

    /// Return the image point which is mapped to the given point in the cartesian space.  This
    /// performs the same operation as `row_col_of` but returns a `f64` point with subpixel
    /// accuracy rather than returning the integer row and column.  Use this when you are not
    /// indexing into an image but are rather working with coordinates in the image space.
    ///
    /// # Arguments
    ///
    /// * `p`: The point in the outer, cartesian space to convert to an image point.
    ///
    /// returns: OPoint<f64, Const<2>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn image_point_of(&self, p: &Point2) -> Point2 {
        let p = self.forward * p;
        let y = (p.y - self.origin.y) / self.px_size;
        let x = (p.x - self.origin.x) / self.px_size;
        Point2::new(x, y)
    }

    pub fn image_index_of(&self, p: &Point2) -> Point2I {
        let img_point = self.image_point_of(p);
        Point2I::new(img_point.x as i32, img_point.y as i32)
    }

    /// Return the point in the cartesian space which is mapped to the given row and column of the
    /// image.
    ///
    /// # Arguments
    ///
    /// * `row`:
    /// * `col`:
    ///
    /// returns: OPoint<f64, Const<2>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn point_of_row_col(&self, row: usize, col: usize) -> Point2 {
        let x = self.origin.x + col as f64 * self.px_size;
        let y = self.origin.y + row as f64 * self.px_size;
        self.inverse * Point2::new(x, y)
    }

    /// Return the point in the cartesian space which is mapped to the given image point.  This is
    /// the inverse of `image_point_of`.
    ///
    /// # Arguments
    ///
    /// * `p`:
    ///
    /// returns: OPoint<f64, Const<2>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn point_of_image_point(&self, img_point: &Point2) -> Point2 {
        let p = self.origin + (img_point.coords * self.px_size);
        self.inverse * p
    }

    /// Takes a pixel index point (represented as `Point2I`) and returns the corresponding point
    /// in the cartesian UV space.
    ///
    /// # Arguments
    ///
    /// * `img_point`: The pixel index point in the image space, represented as `Point2I`.
    ///
    /// returns: OPoint<f64, Const<2>>
    pub fn point_of_image_point_i(&self, img_point: Point2I) -> Point2 {
        let c = Vector2::new(img_point.x as f64, img_point.y as f64);
        let p = self.origin + (c * self.px_size);
        self.inverse * p
    }

    pub fn iter_row_col(&self) -> RowColIter {
        RowColIter::new(self)
    }
}

pub struct RowColIter<'a> {
    mapping: &'a RasterMapping,
    row: usize,
    col: usize,
}

impl<'a> RowColIter<'a> {
    pub fn new(mapping: &'a RasterMapping) -> Self {
        Self {
            mapping,
            row: 0,
            col: 0,
        }
    }
}

impl<'a> Iterator for RowColIter<'a> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.row >= self.mapping.rows() {
            return None;
        }
        let (row, col) = (self.row, self.col);
        self.col += 1;
        if self.col >= self.mapping.cols() {
            self.col = 0;
            self.row += 1;
        }
        Some((row, col))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::Vector2;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_forward_mapping() {
        let origin = Point2::new(-1.0, -1.0);
        let transform = Iso2::new(Vector2::new(2.0, 3.0), PI / 2.0);
        let mapping = RasterMapping::new(origin, (100, 100), 0.10, Some(transform));
        let (row, col) = mapping.row_col_of(&Point2::new(1.0, 0.0));
        assert_eq!(col, 30);
        assert_eq!(row, 50);

        let p2 = mapping.image_point_of(&Point2::new(1.0, 0.0));
        assert_relative_eq!(p2.x, 30.0, epsilon = 0.0001);
        assert_relative_eq!(p2.y, 50.0, epsilon = 0.0001);
    }

    #[test]
    fn test_inverse_mapping() {
        let origin = Point2::new(-1.0, -1.0);
        let transform = Iso2::new(Vector2::new(2.0, 3.0), PI / 2.0);
        let mapping = RasterMapping::new(origin, (100, 100), 0.10, Some(transform));
        let p = mapping.point_of_row_col(50, 30);
        assert_relative_eq!(p.x, 1.0, epsilon = 0.0001);
        assert_relative_eq!(p.y, 0.0, epsilon = 0.0001);

        let p2 = mapping.point_of_image_point(&Point2::new(30.0, 50.0));
        assert_relative_eq!(p2.x, 1.0, epsilon = 0.0001);
        assert_relative_eq!(p2.y, 0.0, epsilon = 0.0001);
    }

    #[test]
    fn test_iterator() {
        let mapping = RasterMapping::new(Point2::new(0.0, 0.0), (5, 10), 1.0, None);

        let mut expected = Vec::new();
        for row in 0..mapping.rows() {
            for col in 0..mapping.cols() {
                expected.push((row, col));
            }
        }

        let mut actual = Vec::new();
        for (row, col) in mapping.iter_row_col() {
            actual.push((row, col));
        }

        assert_eq!(expected, actual);
    }
}
