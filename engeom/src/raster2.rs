//! This module contains tools for working with 2D raster data, such as images, depth maps, and
//! other information which can be represented as a grid of values. It has some image processing
//! tools, but it is not specifically an image processing library.

mod area_average;
mod ball_rolling;
mod index_iter;
mod inpaint;
mod kernel;
mod mapping;
mod mask_ops;
mod raster_mask;
mod region_labeling;
mod roi;
mod roi_mask;
mod scalar_raster;
mod visualize;
mod zhang_suen;

use crate::Result;
use crate::common::{PointNI, VectorNI};
use crate::image::{ImageBuffer, Luma};
use crate::na::{DMatrix, Scalar};
pub use ball_rolling::*;
pub use index_iter::{IndexIter, SizeForIndex};
pub use inpaint::inpaint;
pub use kernel::*;
pub use mapping::RasterMapping;
pub use raster_mask::{RasterMask, RasterMaskTrueIterator};
pub use region_labeling::*;
pub use roi::{RasterRoi, RoiOverlap};
pub use scalar_raster::*;
pub use visualize::*;
pub use zhang_suen::*;

pub type Point2I = PointNI<2>;
pub type Vector2I = VectorNI<2>;

pub trait Point2IIndexAccess<T> {
    /// Get the value at the given point in the raster.
    ///
    /// # Arguments
    ///
    /// * `point`: The point to access.
    ///
    /// returns: T
    fn get_at(&self, point: Point2I) -> Option<T>;

    fn set_at(&mut self, point: Point2I, value: T) -> Result<()>;
}

impl<T: Scalar + Copy> Point2IIndexAccess<T> for DMatrix<T> {
    fn get_at(&self, point: Point2I) -> Option<T> {
        if point.x < 0
            || point.y < 0
            || point.x >= self.ncols() as i32
            || point.y >= self.nrows() as i32
        {
            None
        } else {
            let x = point.x as usize;
            let y = point.y as usize;
            Some(self[(y, x)])
        }
    }

    fn set_at(&mut self, point: Point2I, value: T) -> Result<()> {
        if point.x < 0
            || point.y < 0
            || point.x >= self.ncols() as i32
            || point.y >= self.nrows() as i32
        {
            Err("Point out of bounds".into())
        } else {
            let x = point.x as usize;
            let y = point.y as usize;
            self[(y, x)] = value;
            Ok(())
        }
    }
}
impl Point2IIndexAccess<u16> for ImageBuffer<Luma<u16>, Vec<u16>> {
    fn get_at(&self, point: Point2I) -> Option<u16> {
        if point.x < 0
            || point.y < 0
            || point.x >= self.width() as i32
            || point.y >= self.height() as i32
        {
            None
        } else {
            let x = point.x as u32;
            let y = point.y as u32;
            Some(self.get_pixel(x, y)[0])
        }
    }

    fn set_at(&mut self, point: Point2I, value: u16) -> Result<()> {
        if point.x < 0
            || point.y < 0
            || point.x >= self.width() as i32
            || point.y >= self.height() as i32
        {
            return Err("Point out of bounds".into());
        }
        let x = point.x as u32;
        let y = point.y as u32;
        self.put_pixel(x, y, Luma([value]));
        Ok(())
    }
}

pub trait ToMatrixIndices {
    fn mat_idx(&self) -> (usize, usize);
}

impl ToMatrixIndices for Point2I {
    fn mat_idx(&self) -> (usize, usize) {
        (self.y as usize, self.x as usize)
    }
}

/// Find the minimum and maximum _finite_ values in a DMatrix.
///
/// # Arguments
///
/// * `matrix`: the DMatrix<f64> to search for min and max values.
///
/// returns: (f64, f64)
pub fn d_matrix_min_max(matrix: &DMatrix<f64>) -> (f64, f64) {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;

    for value in matrix.iter() {
        if value.is_finite() {
            if *value < min {
                min = *value;
            }
            if *value > max {
                max = *value;
            }
        }
    }

    (min, max)
}

pub fn d_matrix_mean_stdev(matrix: &DMatrix<f64>) -> (f64, f64) {
    let mut sum = 0.0;
    let mut count = 0.0;

    for value in matrix.iter() {
        if value.is_finite() {
            sum += *value;
            count += 1.0;
        }
    }

    if count == 0.0 {
        return (f64::NAN, f64::NAN);
    }

    let mean = sum / count;

    let mut variance_sum = 0.0;
    for value in matrix.iter() {
        if value.is_finite() {
            variance_sum += (*value - mean).powi(2);
        }
    }

    let stdev = (variance_sum / count).sqrt();

    (mean, stdev)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::na::DMatrix;

    #[test]
    fn d_matrix_min_max_with_nans() {
        let matrix = DMatrix::from_row_slice(
            3,
            3,
            &[
                1.0,
                2.0,
                f64::INFINITY,
                4.0,
                5.0,
                6.0,
                f64::NEG_INFINITY,
                8.0,
                9.0,
            ],
        );
        let (min, max) = d_matrix_min_max(&matrix);
        assert_eq!(min, 1.0);
        assert_eq!(max, 9.0);
    }

    #[test]
    fn matrix_point_index_get() -> Result<()> {
        let mut matrix = DMatrix::from_row_slice(3, 3, &[1i32, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(matrix.get_at(Point2I::new(0, 0)), Some(1));
        assert_eq!(matrix.get_at(Point2I::new(1, 0)), Some(2));
        assert_eq!(matrix.get_at(Point2I::new(2, 0)), Some(3));
        assert_eq!(matrix.get_at(Point2I::new(0, 1)), Some(4));
        assert_eq!(matrix.get_at(Point2I::new(0, 2)), Some(7));

        // Test out of bounds
        assert_eq!(matrix.get_at(Point2I::new(-1, -1)), None);
        assert_eq!(matrix.get_at(Point2I::new(3, 3)), None);

        Ok(())
    }

    #[test]
    fn matrix_point_index_set() -> Result<()> {
        let mut matrix = DMatrix::zeros(3, 3);

        // Test set the first column to 1, 2, 3
        matrix.set_at(Point2I::new(0, 0), 1)?;
        matrix.set_at(Point2I::new(0, 1), 2)?;
        matrix.set_at(Point2I::new(0, 2), 3)?;

        // Test set the rest of the first row to 4, 5
        matrix.set_at(Point2I::new(1, 0), 4)?;
        matrix.set_at(Point2I::new(2, 0), 5)?;

        // Now check the values of the first row
        assert_eq!(matrix[(0, 0)], 1);
        assert_eq!(matrix[(0, 1)], 4);
        assert_eq!(matrix[(0, 2)], 5);

        // Check the rest of the first column
        assert_eq!(matrix[(1, 0)], 2);
        assert_eq!(matrix[(2, 0)], 3);

        Ok(())
    }
}
