//! This module contains tools for working with 2D raster data, such as images, depth maps, and
//! other information which can be represented as a grid of values. It has some image processing
//! tools, but it is not specifically an image processing library.

mod area_average;
mod inpaint;
mod kernel;
mod mapping;
mod mask_ops;
mod scalar_raster;

pub use inpaint::inpaint;
pub use kernel::*;
pub use mapping::RasterMapping;
pub use mask_ops::*;
pub use scalar_raster::*;
