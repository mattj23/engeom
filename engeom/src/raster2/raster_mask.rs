use crate::Result;
use crate::image::{GenericImage, GrayImage, Luma};
use crate::raster2::zhang_suen_thinning;
use imageproc::distance_transform::Norm;
use imageproc::morphology::{dilate_mut, erode_mut};

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

    pub fn set(&mut self, x: u32, y: u32, value: bool) {
        self.buffer.put_pixel(x, y, Luma([value as u8 * 255]));
    }

    pub fn get(&self, x: u32, y: u32) -> bool {
        self.buffer.get_pixel(x, y)[0] > 0
    }

    pub fn width(&self) -> u32 {
        self.buffer.width()
    }

    pub fn height(&self) -> u32 {
        self.buffer.height()
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

    pub fn dilate_alternating_norms_mut(&mut self, count: usize) {
        for i in 0..count {
            if i % 2 == 0 {
                dilate_mut(&mut self.buffer, Norm::L1, 1);
            } else {
                dilate_mut(&mut self.buffer, Norm::LInf, 1);
            }
        }
    }
}

pub struct RasterMaskTrueIterator<'a> {
    mask: &'a RasterMask,
    x: u32,
    y: u32,
}

impl<'a> Iterator for RasterMaskTrueIterator<'a> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        while self.y < self.mask.buffer.height() {
            if self.x >= self.mask.buffer.width() {
                self.x = 0;
                self.y += 1;
            }
            if self.y >= self.mask.buffer.height() {
                return None;
            }
            if self.mask.get(self.x, self.y) {
                let current = (self.x, self.y);
                self.x += 1;
                return Some(current);
            }
            self.x += 1;
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_true() {
        let mut mask = RasterMask::empty(4, 2);
        mask.set(0, 0, true);
        mask.set(1, 1, true);
        mask.set(2, 0, true);

        assert_eq!(mask.count_true(), 3);
    }

    #[test]
    fn not_mut() {
        let mut mask = RasterMask::empty(4, 2);
        mask.set(0, 0, true);
        mask.set(1, 1, true);

        mask.not_mut();

        assert!(!mask.get(0, 0));
        assert!(!mask.get(1, 1));
        assert_eq!(mask.count_true(), 6);
    }

    #[test]
    fn create_empty() {
        let mask = RasterMask::empty(10, 10);
        assert_eq!(mask.width(), 10);
        assert_eq!(mask.height(), 10);
    }

    #[test]
    fn value_set_get() {
        let mut mask = RasterMask::empty(5, 5);
        mask.set(2, 2, true);
        assert!(mask.get(2, 2));
        assert!(!mask.get(1, 1));
    }

    #[test]
    fn buffer_is_row_major_order() {
        // Verify that the buffer is in row-major order by checking the first few pixels
        let mut mask = RasterMask::empty(5, 3);
        mask.set(0, 0, true);
        mask.set(2, 0, true);

        let true_indices = mask
            .buffer
            .as_raw()
            .into_iter()
            .enumerate()
            .filter_map(|(i, &v)| if v > 0 { Some(i) } else { None })
            .collect::<Vec<_>>();

        assert_eq!(true_indices, vec![0, 2]);
    }

    #[test]
    fn iter_true() {
        let mut mask = RasterMask::empty(4, 3);
        mask.set(1, 0, true);
        mask.set(3, 1, true);
        mask.set(2, 2, true);

        let mut iter = mask.iter_true();
        assert_eq!(iter.next(), Some((1, 0)));
        assert_eq!(iter.next(), Some((3, 1)));
        assert_eq!(iter.next(), Some((2, 2)));
        assert_eq!(iter.next(), None);
    }
}
