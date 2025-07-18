//! Uses the `imageproc` crate to perform region labeling on a raster mask.

use faer::prelude::default;
use imageproc::definitions::Image;
use imageproc::region_labelling::connected_components;
pub use imageproc::region_labelling::Connectivity;
use crate::image::{GenericImage, Luma};
use crate::raster2::index_iter::SizeForIndex;
use crate::raster2::roi::RasterRoi;

#[derive(Clone)]
pub struct LabeledRegions {
    buffer: Image<Luma<u32>>,
    roi: Vec<RasterRoi>,
}

impl LabeledRegions {
    pub fn buffer(&self) -> &Image<Luma<u32>> {
        &self.buffer
    }

    pub fn roi(&self) -> &[RasterRoi] {
        &self.roi
    }

    pub fn from_connected_components<I>(image: &I, conn: Connectivity, background: I::Pixel) -> Self
    where
        I: GenericImage,
        I::Pixel: Eq,
    {
        let result = connected_components(image, conn, background);
        let mut regions: Vec<RasterRoi> = Vec::new();
        for i in result.iter_indices() {
            let vi = result.get_pixel(i.x as u32, i.y as u32)[0] as usize;
            if vi == 0 {
                continue;
            }

            // If we don't have an roi up to this value, pad the regions vector with empty ROIs
            // until we reach the value
            while regions.len() <= vi {
                regions.push(default());
            }

            regions[vi].expand_to_contain(i);
        }

        Self { buffer: result, roi: regions }
    }
}