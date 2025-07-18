use super::ScalarImage;
use crate::image::{GenericImageView, GrayImage, Luma, SubImage};
use imageproc::distance_transform::Norm::{L1, LInf};
use imageproc::morphology::dilate_mut;

pub struct AreaAverage<'a> {
    pub depth: SubImage<&'a ScalarImage<u16>>,
    pub mask: SubImage<&'a GrayImage>,
    pub x0: u32,
    pub y0: u32,
    pub width: u32,
    pub height: u32,
    pub radius_px: u32,
}

impl<'a> AreaAverage<'a> {
    pub fn from(
        depth: &'a ScalarImage<u16>,
        mask: &'a GrayImage,
        x: u32,
        y: u32,
        radius_px: u32,
    ) -> Self {
        // Find the working limits of the box
        let half_size = radius_px as i32 + 2;
        let x_min = (x as i32 - half_size).max(0) as u32;
        let x_max = (x as i32 + half_size).min(depth.width() as i32 - 1) as u32;
        let y_min = (y as i32 - half_size).max(0) as u32;
        let y_max = (y as i32 + half_size).min(depth.height() as i32 - 1) as u32;

        // Find the point coordinate in the box
        let x0 = (x as i32 - x_min as i32) as u32;
        let y0 = (y as i32 - y_min as i32) as u32;

        let width = x_max - x_min + 1;
        let height = y_max - y_min + 1;
        let depth_view = depth.view(x_min, y_min, width, height);
        let mask_view = mask.view(x_min, y_min, width, height);
        Self {
            depth: depth_view,
            mask: mask_view,
            x0,
            y0,
            width,
            height,
            radius_px,
        }
    }

    pub fn get_average(&self) -> u16 {
        let area_mask = self.get_area_mask();

        let mut count = 0.0;
        let mut sum = 0.0;
        for (x, y, v) in area_mask.enumerate_pixels() {
            if v[0] == 255 {
                count += 1.0;
                sum += self.depth.get_pixel(x, y)[0] as f32;
            }
        }

        (sum / count) as u16
    }

    fn get_area_mask(&self) -> GrayImage {
        let mut fill = GrayImage::new(self.width, self.height);
        // Set everything to zero
        for (_, _, v) in fill.enumerate_pixels_mut() {
            *v = Luma([0]);
        }

        // Set the center to one
        fill.put_pixel(self.x0, self.y0, Luma([255]));

        for i in 0..self.radius_px {
            // Perform the dilation, alternating between L1 and LInf norms to get an octagon
            if i % 2 == 0 {
                dilate_mut(&mut fill, LInf, 1);
            } else {
                dilate_mut(&mut fill, L1, 1);
            }

            // Now wipe out the pixels which are not in the mask
            // TODO: Should these be saved to a list so they can be accessed faster? or is the
            // allocation overhead too much?
            for (x, y, v) in fill.enumerate_pixels_mut() {
                if self.mask.get_pixel(x, y)[0] < 255 {
                    *v = Luma([0]);
                }
            }
        }

        fill
    }
}
