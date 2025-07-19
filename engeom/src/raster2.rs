//! This module contains tools for working with 2D raster data, such as images, depth maps, and
//! other information which can be represented as a grid of values. It has some image processing
//! tools, but it is not specifically an image processing library.

mod area_average;
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
mod zhang_suen;

use crate::Result;
use crate::common::{PointNI, VectorNI};
use crate::image::{GenericImage, ImageBuffer, Luma};
use crate::na::DMatrix;
pub use inpaint::inpaint;
pub use kernel::*;
pub use mapping::RasterMapping;
pub use raster_mask::{RasterMask, RasterMaskTrueIterator};
pub use region_labeling::*;
pub use roi::RasterRoi;
pub use scalar_raster::*;
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

    fn set_at(&mut self, point: Point2I, value: T) -> Result<()> ;
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
