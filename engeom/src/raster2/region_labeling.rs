//! Uses the `imageproc` crate to perform region labeling on a raster mask.

use crate::image::{GenericImage, Luma};
use crate::raster2::index_iter::SizeForIndex;
use crate::raster2::roi::{RasterRoi, RoiOverlap};
use crate::raster2::roi_mask::RoiMask;
use crate::raster2::{Point2I, RasterMask};
use faer::prelude::default;
use imageproc::definitions::Image;
pub use imageproc::region_labelling::Connectivity;
use imageproc::region_labelling::connected_components;

/// This is a very lightweight temporary structure which contains information about a single
/// labeled region in a `LabeledRegions` object. It contains the label, the region of interest (ROI)
/// for that label, and the count of pixels in that region. It also holds a reference to the
/// `LabeledRegions` object so that it can access the backing image buffer if needed.
///
/// This is meant to be an intermediate object that you can use to get some basic information about
/// an identified region, and can then be used to generate heavier objects while also maintaining
/// tracking information back to the original pixel positions.
pub struct Region<'a> {
    labeled_regions: &'a LabeledRegions,
    label: u32,
    roi: RasterRoi,
    count: usize,
}

impl<'a> Region<'a> {
    /// Create a new `RoiMask` for this region. This will contain a `RasterMask` with only the
    /// pixels that have this label. The mask will be the size of the ROI plus a border of
    /// `padding` extra blank pixels around the ROI.  The `RoiMask` also contains the mappings that
    /// associate the pixels in the mask with positions in the original `LabeledRegions` buffer,
    /// and thus with pixels in the original image from which the `LabeledRegions` was created.
    ///
    /// # Arguments
    ///
    /// * `padding`: an empty border this many pixels wide will be added around the ROI when
    ///   creating the mask. Use this if you need extra space for morphological or other operations
    ///   but still need to keep track of which positions in the mask correspond to the original
    ///   image.
    ///
    /// returns: RoiMask
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn create_roi_mask(&self, padding: u32) -> RoiMask {
        let height = self.roi.extent().y as u32 + 2 * padding;
        let width = self.roi.extent().x as u32 + 2 * padding;
        let mut mask = RasterMask::empty(width, height);
        let mask_roi = self.roi.expanded(padding);
        let overlay = RoiOverlap::new(mask_roi, self.roi);

        for p in overlay.iter_intersection_a() {
            if self.label == self.labeled_regions.label_at(p.parent) {
                // This should be safe because we are iterating over a known set of
                // indices calculated from the bounds of this mask.
                mask.set_point(p.local, true).unwrap();
            }
        }

        RoiMask::new_resized(mask, self.roi, mask_roi)
    }
}

impl<'a> Region<'a> {
    /// Get the `u32` label associated with this region in the `LabeledRegions` buffer.
    pub fn label(&self) -> u32 {
        self.label
    }

    /// Get the minimum sized bounding box on the original buffer which contains all the pixels
    /// in this region.
    pub fn roi(&self) -> &RasterRoi {
        &self.roi
    }

    /// Get the number of pixels which were counted in this region.
    pub fn count(&self) -> usize {
        self.count
    }
}

#[derive(Clone)]
pub struct LabeledRegions {
    buffer: Image<Luma<u32>>,
    raw_roi: Vec<RasterRoi>,
    raw_counts: Vec<usize>,
}

impl LabeledRegions {
    /// Get a reference to the backing image buffer holding all labeled pixels.
    pub fn buffer(&self) -> &Image<Luma<u32>> {
        &self.buffer
    }

    pub fn label_at(&self, p: Point2I) -> u32 {
        if p.x < 0
            || p.y < 0
            || p.x >= self.buffer.width() as i32
            || p.y >= self.buffer.height() as i32
        {
            return 0; // Out of bounds, return background label
        }
        self.buffer.get_pixel(p.x as u32, p.y as u32)[0]
    }

    /// Try to get a `Region` object for a specific label. If the label is 0 or out of bounds,
    /// it returns `None`. Otherwise, it returns a `Region` object containing the label, the
    /// region of interest (ROI) for that label, and the count of pixels in that region.
    ///
    /// The returned `Region` object also contains an immutable reference to the parent
    /// `LabeledRegions` object (this object) so that it can extract a working sub-region mask.
    ///
    /// # Arguments
    ///
    /// * `label`: the numeric label of the region to try to retrieve. Labels start at 1, as 0 is
    ///   the background label.
    ///
    /// returns: Option<Region>
    pub fn get_region(&self, label: u32) -> Option<Region> {
        if label == 0 || label as usize >= self.raw_roi.len() {
            return None; // Label 0 is reserved for the background
        }

        let roi = &self.raw_roi[label as usize];
        let count = self.raw_counts[label as usize];

        Some(Region {
            labeled_regions: self,
            label,
            roi: *roi,
            count,
        })
    }

    pub fn iter(&self) -> RegionIterator {
        RegionIterator {
            labeled: self,
            current: 1, // Start at 1 to skip the background label
        }
    }

    pub fn region_count(&self) -> usize {
        self.raw_roi.len() - 1 // Subtract 1 to exclude the background region labeled with 0
    }

    /// This function constructs a `LabeledRegions` object from the `connected_components` function
    /// in the `imageproc` crate. It creates a backing single channel image buffer where each pixel
    /// is a `u32` value representing the label of the region it belongs to, beginning at 1 and
    /// incrementing for each new region found. The background pixels are labeled with 0.
    ///
    /// # Arguments
    ///
    /// * `image`:
    /// * `conn`:
    /// * `background`:
    ///
    /// returns: LabeledRegions
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_connected<I>(image: &I, conn: Connectivity, background: I::Pixel) -> Self
    where
        I: GenericImage,
        I::Pixel: Eq,
    {
        let result = connected_components(image, conn, background);

        let mut raw_regions: Vec<RasterRoi> = Vec::new();
        let mut raw_counts: Vec<usize> = Vec::new();
        for i in result.iter_indices() {
            let vi = result.get_pixel(i.x as u32, i.y as u32)[0] as usize;
            if vi == 0 {
                continue;
            }

            // If we don't have an roi up to this value, pad the regions and counts vector with
            // default values until we reach the current index.
            while raw_regions.len() <= vi {
                raw_regions.push(default());
                raw_counts.push(0);
            }

            raw_regions[vi].expand_to_contain(i);
            raw_counts[vi] += 1;
        }

        Self {
            buffer: result,
            raw_roi: raw_regions,
            raw_counts,
        }
    }
}

pub struct RegionIterator<'a> {
    labeled: &'a LabeledRegions,
    current: usize,
}

impl<'a> Iterator for RegionIterator<'a> {
    type Item = Region<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let label = self.current as u32;
        self.current += 1;
        self.labeled.get_region(label)
    }
}
