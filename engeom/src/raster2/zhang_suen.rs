//! This module has an implementation of the Zhang-Suen thinning algorithm for binary images, which
//! is very similar to a medial axis transform. It will progressively erode a positive region of
//! pixels until it reaches roughly a single pixel wide representation of the region.

use crate::raster2::Point2I;
use crate::raster2::raster_mask::RasterMask;
use itertools::Itertools;

/// Perform Zhang-Suen thinning on a binary image mask. This algorithm is used to reduce binary
/// regions in the image to a single point wide skeleton, similar to a medial axis transform.
///
/// # Arguments
///
/// * `mask`: the binary image mask to be thinned. The mask should be a `GrayImage` where white
///   pixels will be the ones getting eroded.
///
/// returns: ()
pub fn zhang_suen_thinning(mask: &mut RasterMask) {
    let mut zhang_suen = ZhangSuen::new(mask);
    while zhang_suen.zhang_suen_iter() {}
}

struct ZhangSuen<'a> {
    mask: &'a mut RasterMask,
}

impl<'a> ZhangSuen<'a> {
    fn new(mask: &'a mut RasterMask) -> Self {
        Self { mask }
    }

    /// Returns the number of transitions from not exists to exists moving around the neighbors
    /// of the pixel at x, y
    fn zhang_suen_a(&self, p: Point2I) -> i32 {
        let neighbor_values = zhang_suen_neighbors(p)
            .iter()
            .map(|p| self.mask.get_point(*p))
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
    fn zhang_suen_b(&self, p: Point2I) -> i32 {
        zhang_suen_neighbors(p)
            .iter()
            .map(|p| self.mask.get_point(*p))
            .filter(|x| *x)
            .count() as i32
    }

    /// Returns true if at least one of the three given neighbors of the pixel at x, y does not
    /// exist
    fn zhang_suen_c(&self, p: Point2I, a: usize, b: usize, c: usize) -> bool {
        let neighbor_values = zhang_suen_neighbors(p)
            .iter()
            .map(|p| self.mask.get_point(*p))
            .collect::<Vec<_>>();

        !neighbor_values[a - 2] || !neighbor_values[b - 2] || !neighbor_values[c - 2]
    }

    // fn delete_pixel(&mut self, x: i32, y: i32) {
    //     self.mask.set(x as u32, y as u32, false);
    // }

    // fn has_pixel(&self, p: Point2I) -> bool {
    //     self.mask.get_point(p)
    // }

    pub fn zhang_suen_iter(&mut self) -> bool {
        // https://rosettacode.org/wiki/Zhang-Suen_thinning_algorithm

        // Step 1
        let mut to_delete_1 = Vec::new();
        for p in self.mask.iter_true() {
            let a = self.zhang_suen_a(p);
            let b = self.zhang_suen_b(p);

            if !(2..=6).contains(&b) || a != 1 {
                continue;
            }

            if !self.zhang_suen_c(p, 2, 4, 6) {
                continue;
            }

            if !self.zhang_suen_c(p, 4, 6, 8) {
                continue;
            }

            to_delete_1.push(p);
        }

        for p in to_delete_1.iter() {
            self.mask.set_point_if_in_bounds(*p, false);
        }

        // Step 2
        let mut to_delete_2 = Vec::new();
        for p in self.mask.iter_true() {
            let a = self.zhang_suen_a(p);
            let b = self.zhang_suen_b(p);

            if !(2..=6).contains(&b) || a != 1 {
                continue;
            }

            if !self.zhang_suen_c(p, 2, 4, 8) {
                continue;
            }

            if !self.zhang_suen_c(p, 2, 6, 8) {
                continue;
            }

            to_delete_2.push(p);
        }

        for p in to_delete_2.iter() {
            self.mask.set_point_if_in_bounds(*p, false);
        }

        !to_delete_1.is_empty() || !to_delete_2.is_empty()
    }
}

fn zhang_suen_neighbors(p: Point2I) -> [Point2I; 8] {
    [
        Point2I::new(p.x, p.y - 1),
        Point2I::new(p.x + 1, p.y - 1),
        Point2I::new(p.x + 1, p.y),
        Point2I::new(p.x + 1, p.y + 1),
        Point2I::new(p.x, p.y + 1),
        Point2I::new(p.x - 1, p.y + 1),
        Point2I::new(p.x - 1, p.y),
        Point2I::new(p.x - 1, p.y - 1),
        // (x, y - 1),
        // (x + 1, y - 1),
        // (x + 1, y),
        // (x + 1, y + 1),
        // (x, y + 1),
        // (x - 1, y + 1),
        // (x - 1, y),
        // (x - 1, y - 1),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Result;
    use crate::raster2::raster_mask::RasterMask;

    fn sample_to_mask(input: Vec<&str>) -> Result<RasterMask> {
        let height = input.len();
        let width = input[0].len();
        let mut mask = RasterMask::empty(width as u32, height as u32);

        for (y, line) in input.iter().enumerate() {
            for (x, c) in line.chars().enumerate() {
                if c == '1' {
                    mask.set_point(Point2I::new(x as i32, y as i32), true)?;
                }
            }
        }

        Ok(mask)
    }

    #[test]
    fn sample_io() -> Result<()> {
        // sample from https://rosettacode.org/wiki/Zhang-Suen_thinning_algorithm

        let input = vec![
            "00000000000000000000000000000000",
            "01111111110000000111111110000000",
            "01110001111000001111001111000000",
            "01110000111000001110000111000000",
            "01110001111000001110000000000000",
            "01111111110000001110000000000000",
            "01110111100000001110000111000000",
            "01110011110011101111001111011100",
            "01110001111011100111111110011100",
            "00000000000000000000000000000000",
        ];

        let expected = vec![
            "00000000000000000000000000000000",
            "00111111100000000011111100000000",
            "00100000100000000110000000000000",
            "01000000100000000100000000000000",
            "01000000100000001000000000000000",
            "01111110100000001000000000000000",
            "00000001000000000100000000000000",
            "00000000100001000110000110001000",
            "00000000010000000001111000000000",
            "00000000000000000000000000000000",
        ];

        let mut input = sample_to_mask(input)?;
        let expected = sample_to_mask(expected)?;

        input.buffer.save("D:/temp/k/zs-input.png")?;
        expected.buffer.save("D:/temp/k/zs-expected.png")?;

        zhang_suen_thinning(&mut input);
        input.buffer.save("D:/temp/k/zs-output.png")?;

        assert_eq!(expected.count_true(), input.count_true());

        for p in expected.iter_all() {
            assert_eq!(
                input.get_point(p),
                expected.get_point(p),
                "Pixel mismatch at {:?}",
                p
            );
        }

        Ok(())
    }

    #[test]
    fn test_zhang_suen() {
        let mut mask = RasterMask::empty(50, 50);
        // Create a simple cross shape
        mask.draw_rect_mut(Point2I::new(20, 5), Point2I::new(30, 45), true, true);
        mask.draw_rect_mut(Point2I::new(5, 20), Point2I::new(45, 30), true, true);

        mask.buffer.save("D:/temp/k/zs-input.png").unwrap();

        zhang_suen_thinning(&mut mask);
        mask.buffer.save("D:/temp/k/zs-output.png").unwrap();
    }
}
