use crate::raster2::roi::{RoiIterator, RoiOverlap, RoiPoint};
use crate::raster2::{RasterMask, RasterRoi};

/// A struct that combines an owned working `RasterMask` with two regions of interest on another
/// raster image.  One is the full size of the mask, and the other represents an original overlap of
/// the mask with the other image.
///
/// This is a convenience tool to keep track of operations where you want to create and work with a
/// mask that represents part of a larger image, but you may need to change its size or shape to
/// perform morphological operations, such as those that require an empty border.  In such a case
/// you can no longer guarantee that all pixels of the mask overlap with the original image.
pub struct RoiMask {
    pub mask: RasterMask,
    overlay: RoiOverlap,
}

impl RoiMask {
    pub fn new(mask: RasterMask, roi: RasterRoi) -> Self {
        Self::new_resized(mask, roi, roi)
    }

    pub fn new_resized(mask: RasterMask, original: RasterRoi, mask_roi: RasterRoi) -> Self {
        let overlay = RoiOverlap::new(mask_roi, original);
        Self { mask, overlay }
    }

    /// This will iterate through the points in the mask that are shared with the ROI in the
    /// original image. The point in the `roi` field will correspond to the points in the mask as
    /// expected, and the points in the `parent` field will correspond to the points in the
    /// original image.  However, if the mask has been expanded or contracted, the points which are
    /// visited by the iterator will not necessarily be the full set of points in the mask and/or
    /// the full set of points in the original ROI.
    pub fn iter_shared_points(&self) -> RoiIterator {
        self.overlay.iter_intersection_a()
    }

    pub fn iter_true(&self) -> impl Iterator<Item = RoiPoint> {
        self.iter_shared_points()
            .filter(|p| self.mask.get_point(p.local))
    }
}
