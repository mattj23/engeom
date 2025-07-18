//! This module contains tools for working with 2D raster data, such as images, depth maps, and
//! other information which can be represented as a grid of values. It has some image processing
//! tools, but it is not specifically an image processing library.

mod area_average;
mod inpaint;
mod kernel;
mod mapping;
mod mask_ops;
mod raster_mask;
mod scalar_raster;
mod zhang_suen;

use crate::na::DMatrix;
pub use inpaint::inpaint;
pub use kernel::*;
pub use mapping::RasterMapping;
pub use raster_mask::{RasterMask, RasterMaskTrueIterator};
pub use scalar_raster::*;
pub use zhang_suen::*;

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
