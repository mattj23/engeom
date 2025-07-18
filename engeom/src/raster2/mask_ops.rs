// //! These are operations that can be conceptually performed on GrayImages which serve as masks
// //!
//
// use crate::image::{GrayImage, Luma};
// use imageproc::distance_transform::Norm;
// use imageproc::morphology::{dilate_mut, erode_mut};
//
// pub trait MaskValue {
//     fn set_masked(&mut self);
//     fn set_unmasked(&mut self);
//     fn is_masked(&self) -> bool;
//     fn is_unmasked(&self) -> bool {
//         !self.is_masked()
//     }
// }
//
// impl MaskValue for Luma<u8> {
//     fn set_masked(&mut self) {
//         self[0] = 255;
//     }
//
//     fn set_unmasked(&mut self) {
//         self[0] = 0;
//     }
//
//     fn is_masked(&self) -> bool {
//         self[0] == 255
//     }
// }
//
// pub trait MaskOperations {
//     /// Returns an inverted version of the image mask
//     fn inverted(&self) -> Self;
//
//     /// Returns a mask which is all the pixels in the mask which are part of an unmasked region
//     /// which touches the edge of the image.  This is useful for finding the exterior of a mask.
//     /// For example, if a mask consists of a circle in the center of the image with black holes
//     /// in it, this method will return a white image with a solid black circle in the center.
//     fn get_unmasked_exterior(&self) -> Self;
//
//     /// Return a mask which is the difference of the current mask with positive pixels in the
//     /// other mask removed.
//     fn subtract(&self, other: &Self) -> Self;
//
//     /// Return a mask which is the union of the current mask with the other mask.
//     fn union(&self, other: &Self) -> Self;
//
//     /// Return a mask which is the intersection of the current mask with the other mask.
//     /// This is the same as the logical AND of the two masks.
//     fn intersection(&self, other: &Self) -> Self;
//
//     /// Returns true if the pixel at the given coordinates exists AND is masked, otherwise false.
//     /// This is NOT the same as !is_pixel_unmasked(x, y) because both will return false if the pixel
//     /// does not exist.
//     fn is_pixel_masked(&self, x: i32, y: i32) -> bool;
//
//     /// Returns true if the pixel at the given coordinates exists AND is unmasked, otherwise false.
//     /// This is NOT the same as !is_pixel_masked(x, y) because both will return false if the pixel
//     /// does not exist.
//     fn is_pixel_unmasked(&self, x: i32, y: i32) -> bool;
//
//     /// Set a pixel to be masked (i.e. set to 255)
//     fn set_masked(&mut self, x: u32, y: u32);
//
//     /// Set a pixel to be unmasked (i.e. set to 0)
//     fn set_unmasked(&mut self, x: u32, y: u32);
//
//     /// Erode the mask with alternating L1 and LInf norms
//     fn erode_alternating(&mut self, count: usize);
//
//     /// Dilate the mask with alternating L1 and LInf norms
//     fn dilate_alternating(&mut self, count: usize);
//
//     fn eroded_alternating(&self, count: usize) -> Self;
//     fn dilated_alternating(&self, count: usize) -> Self;
//
//     fn invert(&mut self);
// }
//
// impl MaskOperations for GrayImage {
//     fn inverted(&self) -> GrayImage {
//         let mut inverted = self.clone();
//         for (_, _, v) in inverted.enumerate_pixels_mut() {
//             if v.is_masked() {
//                 v.set_unmasked();
//             } else {
//                 v.set_masked();
//             }
//         }
//         inverted
//     }
//
//     fn get_unmasked_exterior(&self) -> GrayImage {
//         let mut output = GrayImage::new(self.width(), self.height());
//         let mut stack = Vec::new();
//
//         // Push any top and bottom row pixels which are unmasked onto the stack
//         for i in 0..self.width() {
//             if self.get_pixel(i, 0).is_unmasked() {
//                 stack.push((i, 0));
//                 output.set_masked(i, 0);
//             }
//             if self.get_pixel(i, self.height() - 1).is_unmasked() {
//                 stack.push((i, self.height() - 1));
//                 output.set_masked(i, self.height() - 1);
//             }
//         }
//
//         // Push any left or right column pixels which are unmasked onto the stack
//         for i in 0..self.height() {
//             if self.get_pixel(0, i).is_unmasked() {
//                 stack.push((0, i));
//                 output.set_masked(0, i);
//             }
//             if self.get_pixel(self.width() - 1, i).is_unmasked() {
//                 stack.push((self.width() - 1, i));
//                 output.set_masked(self.width() - 1, i);
//             }
//         }
//
//         // Now, for each pixel in the stack, set the output pixel to unmasked and push any
//         // unmasked neighbors onto the stack
//         while let Some((x, y)) = stack.pop() {
//             for (xn, yn) in &[
//                 (x as i32 - 1, y as i32),
//                 (x as i32 + 1, y as i32),
//                 (x as i32, y as i32 - 1),
//                 (x as i32, y as i32 + 1),
//             ] {
//                 if output.is_pixel_unmasked(*xn, *yn)
//                     && self.get_pixel(*xn as u32, *yn as u32).is_unmasked()
//                 {
//                     stack.push((*xn as u32, *yn as u32));
//                     output.set_masked(*xn as u32, *yn as u32);
//                 }
//             }
//         }
//
//         output
//     }
//
//     fn subtract(&self, other: &Self) -> Self {
//         let mut output = self.clone();
//         for (x, y, v) in output.enumerate_pixels_mut() {
//             if other.get_pixel(x, y)[0] == 255 {
//                 *v = Luma([0]);
//             }
//         }
//         output
//     }
//
//     fn union(&self, other: &Self) -> Self {
//         let mut output = self.clone();
//         for (x, y, v) in output.enumerate_pixels_mut() {
//             if other.get_pixel(x, y)[0] == 255 {
//                 *v = Luma([255]);
//             }
//         }
//         output
//     }
//
//     fn intersection(&self, other: &Self) -> Self {
//         let mut output = GrayImage::new(self.width(), self.height());
//         for (x, y, v) in output.enumerate_pixels_mut() {
//             if self.get_pixel(x, y)[0] == 255 && other.get_pixel(x, y)[0] == 255 {
//                 *v = Luma([255]);
//             }
//         }
//         output
//     }
//
//     /// Returns true if the pixel at the given coordinates exists AND is masked, otherwise false.
//     /// Masked means that the pixel value is 255, and conceptually it represents a real pixel with
//     /// a known value.
//     fn is_pixel_masked(&self, x: i32, y: i32) -> bool {
//         if x < 0 || x >= self.width() as i32 || y < 0 || y >= self.height() as i32 {
//             false
//         } else {
//             self.get_pixel(x as u32, y as u32).is_masked()
//         }
//     }
//
//     /// Returns true if the pixel at the given coordinates exists AND is unmasked, otherwise false.
//     /// Unmasked means that the pixel value is 0, and conceptually it represents a pixel that is
//     /// not part of the dataset (similar to null or NaN).
//     fn is_pixel_unmasked(&self, x: i32, y: i32) -> bool {
//         if x < 0 || x >= self.width() as i32 || y < 0 || y >= self.height() as i32 {
//             false
//         } else {
//             self.get_pixel(x as u32, y as u32).is_unmasked()
//         }
//     }
//
//     fn set_masked(&mut self, x: u32, y: u32) {
//         self.get_pixel_mut(x, y).set_masked();
//     }
//
//     fn set_unmasked(&mut self, x: u32, y: u32) {
//         self.get_pixel_mut(x, y).set_unmasked();
//     }
//
//     fn erode_alternating(&mut self, count: usize) {
//         for i in 0..count {
//             if i % 2 == 0 {
//                 erode_mut(self, Norm::L1, 1);
//             } else {
//                 erode_mut(self, Norm::LInf, 1);
//             }
//         }
//     }
//
//     fn dilate_alternating(&mut self, count: usize) {
//         for i in 0..count {
//             if i % 2 == 0 {
//                 dilate_mut(self, Norm::L1, 1);
//             } else {
//                 dilate_mut(self, Norm::LInf, 1);
//             }
//         }
//     }
//
//     fn eroded_alternating(&self, count: usize) -> Self {
//         let mut working = self.clone();
//         working.erode_alternating(count);
//         working
//     }
//
//     fn dilated_alternating(&self, count: usize) -> Self {
//         let mut working = self.clone();
//         working.dilate_alternating(count);
//         working
//     }
//
//     fn invert(&mut self) {
//         for (_, _, v) in self.enumerate_pixels_mut() {
//             if v.is_masked() {
//                 v.set_unmasked();
//             } else {
//                 v.set_masked();
//             }
//         }
//     }
// }
