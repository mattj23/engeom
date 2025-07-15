//! This module contains tools for working with 2D raster data, such as images, depth maps, and
//! other information which can be represented as a grid of values. It has some image processing
//! tools, but it is not specifically an image processing library.

mod mapping;
mod scalar_raster;
mod mask_ops;
mod inpaint;
mod area_average;
mod kernel;

pub use mapping::RasterMapping;
pub use scalar_raster::*;
pub use inpaint::inpaint;
pub use mask_ops::*;
pub use kernel::*;