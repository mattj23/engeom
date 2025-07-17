//! This module has an implementation of the Zhang-Suen thinning algorithm for binary images, which
//! is very similar to a medial axis transform. It will progressively erode a positive region of
//! pixels until it reaches roughly a single pixel wide representation of the region.

use crate::image::GrayImage;
use crate::raster2::{MaskOperations, MaskValue};

/// Perform Zhang-Suen thinning on a binary image mask. This algorithm is used to reduce binary
/// regions in the image to a single point wide skeleton, similar to a medial axis transform.
///
/// # Arguments
///
/// * `mask`: the binary image mask to be thinned. The mask should be a `GrayImage` where white
///   pixels will be the ones getting eroded.
///
/// returns: ()
pub fn zhang_suen_thinning(mask: &mut GrayImage) {
    let mut zhang_suen = ZhangSuen::new(mask);
    while zhang_suen.zhang_suen_iter() {}
}

struct ZhangSuen<'a> {
    mask: &'a mut GrayImage,
}

impl<'a> ZhangSuen<'a> {
    fn new(mask: &'a mut GrayImage) -> Self {
        Self { mask }
    }

    /// Returns the number of transitions from not exists to exists moving around the neighbors
    /// of the pixel at x, y
    fn zhang_suen_a(&self, x: i32, y: i32) -> i32 {
        let neighbor_values = zhang_suen_neighbors(x, y)
            .iter()
            .map(|(x, y)| self.has_pixel(*x, *y))
            .collect::<Vec<_>>();

        let mut count = 0;
        for i in 0..neighbor_values.len() {
            let j = (i + 1) % neighbor_values.len();
            if !neighbor_values[i] && neighbor_values[j] {
                count += 1;
            }
        }

        count
    }

    /// Returns the total number of existing neighbors of the pixel at x, y
    fn zhang_suen_b(&self, x: i32, y: i32) -> i32 {
        zhang_suen_neighbors(x, y)
            .iter()
            .map(|(x, y)| self.has_pixel(*x, *y))
            .filter(|x| *x)
            .count() as i32
    }

    /// Returns true if at least one of the three given neighbors of the pixel at x, y does not
    /// exist
    fn zhang_suen_c(&self, x: i32, y: i32, a: usize, b: usize, c: usize) -> bool {
        let neighbor_values = zhang_suen_neighbors(x, y)
            .iter()
            .map(|(x, y)| self.has_pixel(*x, *y))
            .collect::<Vec<_>>();

        !neighbor_values[a - 2] || !neighbor_values[b - 2] || !neighbor_values[c - 2]
    }

    fn delete_pixel(&mut self, x: i32, y: i32) {
        self.mask.set_unmasked(x as u32, y as u32);
    }

    fn has_pixel(&self, x: i32, y: i32) -> bool {
        if x < 0 || y < 0 || x >= self.mask.width() as i32 || y >= self.mask.height() as i32 {
            return false;
        }
        self.mask.get_pixel(x as u32, y as u32).is_masked()
    }

    pub fn zhang_suen_iter(&mut self) -> bool {
        // https://rosettacode.org/wiki/Zhang-Suen_thinning_algorithm

        // Step 1
        let mut to_delete_1 = Vec::new();
        for (x, y, v) in self.mask.enumerate_pixels() {
            if !v.is_masked() {
                continue;
            }

            let x = x as i32;
            let y = y as i32;

            let a = self.zhang_suen_a(x, y);
            let b = self.zhang_suen_b(x, y);

            if !(2..=6).contains(&b) || a != 1 {
                continue;
            }

            if !self.zhang_suen_c(x, y, 2, 4, 6) {
                continue;
            }

            if !self.zhang_suen_c(x, y, 4, 6, 8) {
                continue;
            }

            to_delete_1.push((x, y));
        }

        for (x, y) in to_delete_1.iter() {
            self.delete_pixel(*x, *y);
        }

        // Step 2
        let mut to_delete_2 = Vec::new();
        for (x, y, v) in self.mask.enumerate_pixels() {
            if !v.is_masked() {
                continue;
            }
            let x = x as i32;
            let y = y as i32;

            let a = self.zhang_suen_a(x, y);
            let b = self.zhang_suen_b(x, y);

            if !(2..=6).contains(&b) || a != 1 {
                continue;
            }

            if !self.zhang_suen_c(x, y, 2, 4, 8) {
                continue;
            }

            if !self.zhang_suen_c(x, y, 2, 6, 8) {
                continue;
            }

            to_delete_2.push((x, y));
        }

        for (x, y) in to_delete_2.iter() {
            self.delete_pixel(*x, *y);
        }

        !to_delete_1.is_empty() || !to_delete_2.is_empty()
    }
}

fn zhang_suen_neighbors(x: i32, y: i32) -> [(i32, i32); 8] {
    [
        (x, y - 1),
        (x + 1, y - 1),
        (x + 1, y),
        (x + 1, y + 1),
        (x, y + 1),
        (x - 1, y + 1),
        (x - 1, y),
        (x - 1, y - 1),
    ]
}
