//! This module contains a scalar field raster type for representing 2d raster fields backed by
//! an image buffer. It can be used for depth maps, intensity maps, and other discretized scalar
//! fields in a way that allows for image processing algorithms and operations to be applied
//! without losing a connection to a spatial coordinate system.

use crate::Result;
use crate::image::imageops::{FilterType, resize};
use crate::image::{GrayImage, ImageBuffer, ImageFormat, ImageReader, Luma, Rgba, RgbaImage};
use crate::na::DMatrix;
use crate::raster2::area_average::AreaAverage;
use crate::raster2::raster_mask::RasterMask;
use crate::raster2::{FastApproxKernel, RasterKernel, inpaint};
use colorgrad::Gradient;
use imageproc::distance_transform::Norm::L1;
use imageproc::morphology::{dilate_mut, erode_mut};
use rayon::prelude::*;
use std::fs::File;
use std::io::{BufReader, BufWriter, Cursor, Read, Write};
use std::path::Path;

pub type ScalarImage<T> = ImageBuffer<Luma<T>, Vec<T>>;

/// A `ScalarRaster` is a special 2D raster type that encodes scalar floating point values into a
/// `u16` image buffer, with a companion `u8` mask image that indicates which pixels are considered
/// "valid" or "real".  This is a lossy representation of the scalar field, but it has two
/// distinct advantages:
///
/// 1. It allows for the entire spectrum of image processing algorithms to be applied directly to
///    the raster and its values, such as filtering, dilation, erosion, inpainting, convolutions,
///    thresholding, blob detection, etc.
/// 2. It allows for the raster to be serialized and deserialized to a compact binary format using
///    any sort of image serialization format, including ones with efficient compression like PNG.
///
/// Conceptually, a `ScalarRaster` differs from the concept of a plain image buffer or a plain
/// floating point matrix in that it serves as a dual representation of a discretized scalar field
/// that still has a connection to a spatial coordinate system.  It does this in several ways:
///
/// - Positions in the 2D field represented by the `ScalarRaster` have dual meaning. They are
///   addressable by pixel coordinates, but pixel coordinates are anchored to physical world space
///   through the `px_size` field, which is the physical world length of each pixel in the raster.
///   For any coordinate in the raster a relative position can be calculated by multiplying the
///   pixel coordinates by `px_size`. Through the use of a `RasterMapping`, this can be converted
///   to absolute coordinates.
///
/// - The `px_size` field also allows for operations performed on the raster to be aware of the
///   physical scale of their operations. For example, when convolving with a Gaussian kernel, the
///   kernel size can be specified in physical units (for example, mm) based entirely on
///   information encoded directly in the `ScalarRaster`.
///
/// - When the `ScalarRaster` is scaled, the `px_size` field is updated to reflect the scaling
///   operation so that the distance between features in the raster stays the same. If the raster
///   is shrunk to half of its original size, the `px_size` is doubled, as each pixel now represents
///   twice as much of the physical world as it did before.
///
/// - The `min_z` and `max_z` fields allow for a direct mapping between floating point and `u16`
///   values, so that operations can use either the integer representation or the floating point
///   value.
///
/// The `ScalarRaster` exists to be a convenient way to work with rasterized scalar fields by
/// tracking spatial information for the client code while providing a consistent interface to
/// both image processing and floating point operations.
///
/// Some notes on its use:
///
/// - Take note that the `z_min` and `z_max` values are encoding limits for the `u16` values in the
///   `values` field; they are not necessarily the minimum and maximum values being stored in the
///   raster. They _are_ the minimum and maximum _possible_ values that can be encoded in the
///   raster, so they should be set high enough to accommodate the range of values you expect to
///   produce in any mutating operation on the raster.  The encoding is a direct linear mapping, so
///   a floating point value of `z_min` will be encoded as `0`, and a floating point value of
///   `z_max` will be encoded as `65535`.  Any floating point value outside of this range will be
///   clamped to the nearest limit.
///
/// - Take care to calculate the lost resolution of the `u16` datatype. If you know the resolution
///   you need, you can see what the representable range is by multiplying it by 2^16. For instance,
///   if the pixel values represent depth and your application requires a resolution of 1µm, then
///   the full range of what you can represent as a depth map is 65.536mm. Or, if you know what the
///   maximum and minimum values you need to represent are, you can calculate the resolution by
///   dividing the range by 65536. For example, if you need to represent depths between 0 and 10m,
///   then the resolution is 10m / 65536 = 0.1525mm. Values that are closer to each other than this
///   number will be compressed to the same `u16` value, and thus will not be distinguishable.
///   For most engineering applications this has more resolution than the certainty of any
///   measurement equipment, but for applications that require high precision over a very large
///   range, you may be better served by using a floating point matrix directly.
///
/// - The `ScalarRaster` has built in serialization and deserialization methods that allow it to
///   be turned into bytes and/or saved directly to a file with all of its internal information
///   preserved. This involves serializing the `u16` and `u8` values as grayscale images with an
///   encoding of your choice, and then stuffing them into a single binary file along with the
///   pixel size and encoding z limits.
///
/// - The fact that the `ScalarRaster` tracks pixel size through rescaling operations allows for
///   fast approximate convolutions to be performed on the raster. This works by having an entity
///   that provides kernels based on a physical size, and then allowing the `ScalarRaster` to
///   shrink itself until a kernel of an optimal pixel count is reached. The convolution is
///   performed on both the lower resolution raster and the mask, and when the image is upscaled
///   the mask is used to determine which pixels are valid in the output. Invalid pixels are
///   reconvolved at a larger resolution to fill in the gaps. This compensates for the tendency
///   to create impractically large kernels when needing to base them on physical dimensions, in
///   exchange for an approximate result.
#[derive(Clone, Debug)]
pub struct ScalarRaster {
    pub values: ScalarImage<u16>,
    pub mask: RasterMask,
    pub px_size: f64,
    pub min_z: f64,
    pub max_z: f64,
}

impl ScalarRaster {
    pub fn serialized_bytes(&self, fmt: ImageFormat) -> Vec<u8> {
        let mut value_bytes = Vec::new();
        let mut mask_bytes = Vec::new();
        self.values
            .write_to(&mut Cursor::new(&mut value_bytes), fmt)
            .unwrap();
        self.mask
            .buffer
            .write_to(&mut Cursor::new(&mut mask_bytes), fmt)
            .unwrap();

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.px_size.to_le_bytes());
        bytes.extend_from_slice(&self.min_z.to_le_bytes());
        bytes.extend_from_slice(&self.max_z.to_le_bytes());
        bytes.extend_from_slice(&(value_bytes.len() as u64).to_le_bytes());
        bytes.extend(value_bytes);
        bytes.extend_from_slice(&(mask_bytes.len() as u64).to_le_bytes());
        bytes.extend(mask_bytes);

        bytes
    }

    pub fn from_serialized_bytes(bytes: &[u8]) -> Result<Self> {
        let mut reader = Cursor::new(bytes);
        Self::from_reader(&mut reader)
    }

    fn from_reader(reader: &mut impl Read) -> Result<Self> {
        let mut bytes_8 = [0u8; 8];
        reader.read_exact(&mut bytes_8)?;
        let px_size = f64::from_le_bytes(bytes_8);

        reader.read_exact(&mut bytes_8)?;
        let min_z = f64::from_le_bytes(bytes_8);

        reader.read_exact(&mut bytes_8)?;
        let max_z = f64::from_le_bytes(bytes_8);

        reader.read_exact(&mut bytes_8)?;
        let value_len = u64::from_le_bytes(bytes_8) as usize;
        let mut value_bytes = vec![0u8; value_len];
        reader.read_exact(&mut value_bytes)?;

        reader.read_exact(&mut bytes_8)?;
        let mask_len = u64::from_le_bytes(bytes_8) as usize;
        let mut mask_bytes = vec![0u8; mask_len];
        reader.read_exact(&mut mask_bytes)?;

        let values = ImageReader::new(Cursor::new(value_bytes))
            .with_guessed_format()?
            .decode()?;

        let mask_buffer = ImageReader::new(Cursor::new(mask_bytes))
            .with_guessed_format()?
            .decode()?;

        if values.width() != mask_buffer.width() || values.height() != mask_buffer.height() {
            return Err("Values and mask images must have the same dimensions".into());
        }

        let mask = RasterMask::new(mask_buffer.into_luma8());

        Ok(Self {
            values: values.into_luma16(),
            mask,
            px_size,
            min_z,
            max_z,
        })
    }

    pub fn save_combined(&self, path: &Path, fmt: ImageFormat) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writer.write_all(&self.serialized_bytes(fmt))?;
        Ok(())
    }

    pub fn load_combined(path: &Path) -> Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        Self::from_reader(&mut reader)
    }

    pub fn f_at(&self, x: i32, y: i32) -> f64 {
        self.u_at(x, y).map(|u| self.u_to_f(u)).unwrap_or(f64::NAN)
    }

    pub fn u_at(&self, x: i32, y: i32) -> Option<u16> {
        if self.mask.get(x as u32, y as u32) {
            None
        } else {
            Some(self.values.get_pixel(x as u32, y as u32)[0])
        }
    }

    pub fn width(&self) -> u32 {
        self.values.width()
    }

    pub fn height(&self) -> u32 {
        self.values.height()
    }

    fn f_to_u(&self, val: f64) -> u16 {
        let f = (val.clamp(self.min_z, self.max_z) - self.min_z) / (self.max_z - self.min_z);
        (f * (u16::MAX - u16::MIN) as f64) as u16 + u16::MIN
    }

    fn u_to_f(&self, val: u16) -> f64 {
        let f = (val - u16::MIN) as f64 / (u16::MAX - u16::MIN) as f64;
        f * (self.max_z - self.min_z) + self.min_z
    }

    pub fn set_f_at(&mut self, x: i32, y: i32, val: f64) {
        if val.is_nan() {
            self.mask.set(x as u32, y as u32, false);
        } else {
            self.mask.set(x as u32, y as u32, true);
            let u = self.f_to_u(val);
            self.values.put_pixel(x as u32, y as u32, Luma([u]));
        }
    }

    pub fn set_u_at(&mut self, x: i32, y: i32, val: Option<u16>) {
        if let Some(v) = val {
            self.mask.set(x as u32, y as u32, true);
            self.values.put_pixel(x as u32, y as u32, Luma([v]));
        } else {
            self.mask.set(x as u32, y as u32, false);
            self.values.put_pixel(x as u32, y as u32, Luma([0]));
        }
    }

    pub fn empty_like(other: &Self) -> Self {
        Self::empty(
            other.width(),
            other.height(),
            other.px_size,
            other.min_z,
            other.max_z,
        )
    }

    pub fn empty(width: u32, height: u32, px_size: f64, min_z: f64, max_z: f64) -> Self {
        let values = ScalarImage::new(width, height);
        let mask = RasterMask::new(GrayImage::new(width, height));

        ScalarRaster {
            values,
            mask,
            px_size,
            min_z,
            max_z,
        }
    }

    /// This function will render the depth map to an image file using a color gradient map from
    /// the `colorgrad` crate.  Optionally, you can provide a tuple of `(min_z, max_z)` to clip
    /// the values to a specific range, otherwise the z_min and z_max of the raster will be used.
    ///
    /// This function will output RGBA pixels, where the alpha channel is set to 0 for invalid
    /// pixels. When specifying the output path, be sure to use the extension of an image format
    /// that supports transparency, such as PNG.
    ///
    /// This function is useful for visualization and debugging of raw value maps, but it is not a
    /// direct replacement for tools like `matplotlib`'s `imshow` function, which can do things
    /// like output a scale bar and perform intelligent scaling.
    ///
    /// # Arguments
    ///
    /// * `path`: the path to the output image file
    /// * `gradient`: a `Gradient` trait object that provides the color mapping
    /// * `min_max`: an optional tuple of `(min_z, max_z)` to clip the values to a specific range
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    pub fn render_with_cmap(
        &self,
        path: &Path,
        gradient: &dyn Gradient,
        min_max: Option<(f64, f64)>,
    ) -> Result<()> {
        let (min_z, max_z) = min_max.unwrap_or((self.min_z, self.max_z));
        let mut img = RgbaImage::new(self.width(), self.height());

        for (x, y, p) in img.enumerate_pixels_mut() {
            if !self.mask.get(x, y) {
                *p = Rgba([0, 0, 0, 0]); // Transparent pixel
            } else {
                let value = self.u_to_f(self.values.get_pixel(x, y)[0]);
                let f = (value - min_z) / (max_z - min_z);

                let color = gradient.at(f as f32).to_rgba8();
                *p = Rgba(color);
            }
        }

        img.save(path).map_err(|e| e.into())
    }

    /// Create a new `ScalarRaster` from a matrix of floating point values. The physical world
    /// length of each pixel should be specified in `px_size`, and the minimum and maximum z values
    /// for the `u16` encoding should be specified in `min_z` and `max_z`. The values in the matrix
    /// will be clamped to the range `[min_z, max_z]` before being converted to `u16` values.
    ///
    /// # Arguments
    ///
    /// * `matrix`: a 2D matrix of floating point values to be converted into the raster
    /// * `px_size`: a physical world length of each pixel in the raster, in whatever units you're
    ///   working with.
    /// * `min_z`: The lower bound of the z range which will be mapped to the `u16` range. This
    ///   value will be a 0 in the resulting raster.
    /// * `max_z`: The upper bound of the z range which will be mapped to the `u16` range. This
    ///   value will be a `u16::MAX` in the resulting raster.
    ///
    /// returns: ScalarRaster
    pub fn from_matrix(matrix: &DMatrix<f64>, px_size: f64, min_z: f64, max_z: f64) -> Self {
        let mut value = ScalarImage::new(matrix.ncols() as u32, matrix.nrows() as u32);
        let mut mask = RasterMask::empty_like(&value);

        // TODO: make generic later
        let type_min = u16::MIN;
        let type_max = u16::MAX;
        let type_range = (type_max - type_min) as f64;

        for i in 0..matrix.nrows() {
            for j in 0..matrix.ncols() {
                let v = matrix[(i, j)];
                if !v.is_finite() {
                    continue;
                }

                // The original value in the cell
                let v = v.clamp(min_z, max_z);

                // Find the fraction of the full scale range that this value represents
                let f = (v - min_z) / (max_z - min_z);

                // Scale it to the type range and convert to the target type
                let r = (f * type_range) as u16 + type_min;

                // Set the pixel value and mask
                value.put_pixel(j as u32, i as u32, Luma([r]));
                mask.set(j as u32, i as u32, true);
            }
        }

        ScalarRaster {
            values: value,
            mask,
            px_size,
            min_z,
            max_z,
        }
    }

    /// Calculates the mean and standard deviation of the valid pixels in the depth map, using
    /// their floating point values. The results are returned as a tuple of `(mean, stdev)`.
    pub fn mean_stdev(&self) -> (f64, f64) {
        let mut sum = 0.0;
        let mut count = 0.0;

        for (x, y, pixel) in self.values.enumerate_pixels() {
            if self.mask.get(x, y) {
                sum += self.u_to_f(pixel.0[0]);
                count += 1.0;
            }
        }

        if count == 0.0 {
            return (f64::NAN, f64::NAN);
        }

        let mean = sum / count;

        let mut variance_sum = 0.0;
        for (x, y, pixel) in self.values.enumerate_pixels() {
            if self.mask.get(x, y) {
                let value = self.u_to_f(pixel.0[0]);
                variance_sum += (value - mean).powi(2);
            }
        }

        let stdev = (variance_sum / count).sqrt();
        (mean, stdev)
    }

    /// Converts the `ScalarRaster` to a `DMatrix<f64>` where each pixel is represented by its
    /// floating point value. Pixels that are masked (i.e., not valid) will be represented as `NaN`
    /// in the matrix.
    pub fn to_matrix(&self) -> DMatrix<f64> {
        let mut matrix =
            DMatrix::from_element(self.height() as usize, self.width() as usize, f64::NAN);

        for (x, y) in self.mask.iter_true() {
            matrix[(y as usize, x as usize)] = self.u_to_f(self.values.get_pixel(x, y).0[0]);
        }
        matrix
    }

    /// This function converts the `ScalarRaster` to a tuple of two `DMatrix<f64>` objects, the
    /// first being the value matrix and the second being the mask matrix. These matrices are
    /// prepared for the normalized convolution operation, according to the following rules:
    ///
    /// - The value matrix will have all NaN values replaced with 0.0.
    /// - For every position in the raster's mask that has a 0 value, the corresponding position
    ///   in the mask matrix will be set to 0.0.  Every position in the mask matrix that has a
    ///   255 value will be set to 1.0.
    ///
    /// This output is specifically made to be compatible with the requirements of the
    /// `RasterKernel::convolve` method, but is a general purpose conversion that re-conceptualizes
    /// the raster as a set of values and certainties.
    pub fn to_value_and_mask_matrices(&self) -> (DMatrix<f64>, DMatrix<f64>) {
        let mut matrix = DMatrix::zeros(self.height() as usize, self.width() as usize);
        let mut mask = DMatrix::zeros(matrix.nrows(), matrix.ncols());

        for (x, y) in self.mask.iter_true() {
            matrix[(y as usize, x as usize)] = self.u_to_f(self.values.get_pixel(x, y).0[0]);
            mask[(y as usize, x as usize)] = 1.0;
        }

        (matrix, mask)
    }

    /// Performs an inpainting operation on the depth map using Alexander Telea's method. Provide
    /// a mask image which indicates the pixels to be inpainted, and a radius in pixels which is
    /// the maximum distance to search for valid pixels to use in the inpainting operation.
    ///
    /// The operation will be performed in place, modifying the `values` and `mask` fields of the
    /// `ScalarRaster`. The `mask` will be updated to include the pixels that were inpainted.
    ///
    /// # Arguments
    ///
    /// * `mask`: a `RasterMask` indicating which pixels to inpaint. The mask should have the same
    ///   dimensions as the `ScalarRaster`
    /// * `inpaint_radius`: the radius in pixels to search for valid pixels to use in the
    ///   inpainting operation
    ///
    /// returns: ()
    pub fn inpaint(&mut self, mask: &RasterMask, inpaint_radius: usize) {
        let filled = inpaint(&self.values, &mask.buffer, &self.mask, inpaint_radius);
        self.values = filled;
        self.mask = self
            .mask
            .clone()
            .or(mask)
            .expect("Mask and raster must have same dimensions");
    }

    pub fn erode(&mut self, pixels: u8) {
        erode_mut(&mut self.mask.buffer, L1, pixels);
        for (x, y, v) in self.values.enumerate_pixels_mut() {
            if self.mask.get(x, y) {
                *v = Luma([0]);
            }
        }
    }

    pub fn delete_by_mask(&mut self, mask: &RasterMask) {
        for (x, y, v) in self.values.enumerate_pixels_mut() {
            if mask.get(x, y) {
                *v = Luma([0]);
                self.mask.set(x, y, false);
            }
        }
    }

    /// Clear a row of pixels by both setting the values to zero and the mask to zero
    fn clear_row(&mut self, y: u32) {
        for x in 0..self.values.width() {
            self.values.put_pixel(x, y, Luma([0]));
            self.mask.set(x, y, false);
        }
    }

    /// Clear a column of pixels by both setting the valeus to zero and the mask to zero
    fn clear_col(&mut self, x: u32) {
        for y in 0..self.values.height() {
            self.values.put_pixel(x, y, Luma([0]));
            self.mask.set(x, y, false);
        }
    }

    /// Clear the boundary of the depth map by setting the first and last rows and columns to zero
    /// in both the depth and mask images.
    pub fn clear_boundary(&mut self) {
        self.clear_row(0);
        self.clear_row(self.values.height() - 1);
        self.clear_col(0);
        self.clear_col(self.values.width() - 1);
    }

    /// Get a mask of just the border pixels by using the L1 erode operation on the current mask
    /// and comparing the difference.
    pub fn get_border_mask(&self, pixels: u8) -> GrayImage {
        let mut border_mask = self.mask.clone();
        erode_mut(&mut border_mask.buffer, L1, pixels);
        diff_img(&self.mask.buffer, &border_mask.buffer)
    }

    pub fn remove_border_outliers(&mut self) -> usize {
        let radius = (1.0 / self.px_size).ceil() as i32;
        let border_mask = self.get_border_mask(5);
        let outliers = border_mask
            .enumerate_pixels()
            .filter(|(_, _, v)| v[0] != 0)
            .map(|(x, y, _)| (x, y, self.is_px_outlier(x, y, radius)))
            .filter(|(_, _, outlier)| *outlier)
            .collect::<Vec<_>>();

        let mut count = 0;
        for (x, y, yes) in outliers.iter() {
            if *yes {
                self.mask.set(*x, *y, false);
                self.values.put_pixel(*x, *y, Luma([0]));
                count += 1;
            }
        }

        count
    }

    fn is_px_outlier(&self, x: u32, y: u32, radius: i32) -> bool {
        if self.mask.get(x, y) {
            false
        } else {
            let mut values = Vec::new();
            for i in -radius..=radius {
                for j in -radius..=radius {
                    let x_loc = x as i32 + i;
                    let y_loc = y as i32 + j;

                    if x_loc < 0
                        || y_loc < 0
                        || x_loc >= self.values.width() as i32
                        || y_loc >= self.values.height() as i32
                    {
                        continue;
                    }

                    if self.mask.get(x_loc as u32, y_loc as u32) {
                        continue;
                    }

                    values.push(self.values.get_pixel(x_loc as u32, y_loc as u32)[0] as f32);
                }
            }

            // Find the mean
            let mean = values.iter().sum::<f32>() / values.len() as f32;

            // Find the standard deviation
            let variance =
                values.iter().map(|v| (v - mean).powi(2)).sum::<f32>() / values.len() as f32;
            let std_dev = variance.sqrt();

            // Find if we're more than 2 standard deviations away from the mean
            let this_px = self.values.get_pixel(x, y)[0] as f32;
            (this_px - mean).abs() > 3.0 * std_dev
        }
    }

    pub fn create_scaled(&self, scale: f64) -> Self {
        let width = (self.values.width() as f64 * scale) as u32;
        let height = (self.values.height() as f64 * scale) as u32;

        let values = resize(&self.values, width, height, FilterType::Nearest);
        let buffer = resize(&self.mask.buffer, width, height, FilterType::Nearest);
        let mask = RasterMask::new(buffer);

        Self {
            values,
            mask,
            px_size: self.px_size / scale,
            min_z: self.min_z,
            max_z: self.max_z,
        }
    }

    pub fn create_shrunk(&self, shrink_factor: u32) -> Self {
        let width = self.values.width() / shrink_factor;
        let height = self.values.height() / shrink_factor;

        let values = resize(&self.values, width, height, FilterType::Nearest);
        let buffer = resize(&self.mask.buffer, width, height, FilterType::Nearest);
        let mask = RasterMask::new(buffer);

        Self {
            values,
            mask,
            px_size: self.px_size * shrink_factor as f64,
            min_z: self.min_z,
            max_z: self.max_z,
        }
    }

    /// Given a reference value image, this function will return a new ScalarRaster which is has the
    /// values of the reference subtracted from the values of this image.
    pub fn subtract(&self, reference: &Self) -> Result<Self> {
        let mut corrected = self.clone();
        for (x, y, v) in corrected.values.enumerate_pixels_mut() {
            if !self.mask.get(x, y) {
                *v = Luma([0]);
            } else {
                let ref_val = reference.u_to_f(reference.values.get_pixel(x, y)[0]);
                let self_val = self.u_to_f(v[0]);
                let f = self_val - ref_val;
                if f.is_nan() {
                    *v = Luma([0]);
                    continue;
                }
                if f < self.min_z || f > self.max_z {
                    return Err(format!(
                        "Value {} out of bounds for raster with min_z {} and max_z {}",
                        f, self.min_z, self.max_z
                    )
                    .into());
                }
                *v = Luma([self.f_to_u(self_val - ref_val)]);
            }
        }

        Ok(corrected)
    }

    /// Find the average value of the pixels in the area neighborhood around the given pixel. The
    /// neighborhood is defined by a radius in pixels, but the neighborhood is found by alternating
    /// dilation of the L1 and LInf norms from the seed pixel constrained by the image mask. This
    /// produces an octagonal neighborhood which cannot jump across mask regions, and should
    /// better approximate a local neighborhood reachable along the object surface.
    ///
    /// It is, however, very expensive to do in bulk.
    ///
    /// # Arguments
    ///
    /// * `x`: the x coordinate of the pixel to find the area average around
    /// * `y`: the y coordinate of the pixel to find the area average around
    /// * `radius_px`: the radius in pixels to search for valid pixels to use in the area average
    ///
    /// returns: u16
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn area_average_px(&self, x: u32, y: u32, radius_px: u32) -> u16 {
        let area = AreaAverage::from(&self.values, &self.mask.buffer, x, y, radius_px);
        area.get_average()
    }

    /// Performs a straightforward area average of every valid pixel in entire image, using the
    /// rayon library for parallelization.
    fn area_average(&self, radius_mm: f32) -> ScalarImage<u16> {
        let radius = (radius_mm / self.px_size as f32).ceil() as u32;
        let mut blurred = ScalarImage::new(self.width(), self.height());

        // Find all of unmasked pixels in the image and collect them into a list for the parallel
        // iterator to use. I don't know if this is more efficient than having the individual
        // threads look at the mask to reduce the vector collection, or worse because of the
        // extra allocation.
        let xys = self
            .mask
            .buffer
            .enumerate_pixels()
            .filter(|(_, _, p)| (*p)[0] == 255)
            .map(|(x, y, _)| (x, y))
            .collect::<Vec<_>>();

        // Now we use rayon to iterate through the unmasked pixels in parallel, collecting the
        // result into a vector.
        let result = xys
            .par_iter()
            .map(|(x, y)| (x, y, self.area_average_px(*x, *y, radius)))
            .collect::<Vec<_>>();

        // Finally we write the results into the buffer
        for (x, y, v) in result {
            blurred.put_pixel(*x, *y, Luma([v]));
        }

        blurred
    }

    /// Performs an area average of every valid pixel in the image, but preferentially uses a
    /// resized image which has (presumably) had the same operation performed as the initial source
    /// of values. This will then fill in any missing values by directly computing them on the
    /// full sized image, which can be quite costly.
    fn area_average_with_resized(&self, radius_mm: f32, resized: &Self) -> ScalarImage<u16> {
        let mut blurred = ScalarImage::new(self.width(), self.height());
        let mut transferred = GrayImage::new(self.mask.width(), self.mask.height());
        let radius = (radius_mm / self.px_size as f32).ceil() as u32;
        let scale_x = resized.mask.width() as f32 / self.mask.width() as f32;
        let scale_y = resized.mask.height() as f32 / self.mask.height() as f32;

        // All the pixels in the resized image which have at a full mask value can be copied
        // directly and put into the corresponding pixels in the blurred image
        for (x, y, v) in blurred.enumerate_pixels_mut() {
            if !self.mask.get(x, y) {
                continue;
            }
            let sx = (x as f32 * scale_x).floor() as u32;
            let sy = (y as f32 * scale_y).floor() as u32;
            if !resized.mask.get(sx, sy) {
                continue;
            }

            *v = Luma([resized.values.get_pixel(sx, sy)[0]]);
            transferred.put_pixel(x, y, Luma([255]));
        }

        // Now we'll manually transfer the pixels which are not in the resized image
        let xys = self
            .mask
            .buffer
            .enumerate_pixels()
            .filter(|(x, y, p)| (*p)[0] != 0 && transferred.get_pixel(*x, *y)[0] == 0)
            .map(|(x, y, _)| (x, y))
            .collect::<Vec<_>>();

        let result = xys
            .par_iter()
            .map(|(x, y)| (x, y, self.area_average_px(*x, *y, radius)))
            .collect::<Vec<_>>();

        for (x, y, v) in result {
            blurred.put_pixel(*x, *y, Luma([v]));
        }

        blurred
    }

    /// Convolve this scalar raster with the given kernel, returning a new scalar raster. This is
    /// identical to using the `RasterKernel`'s `convolve` method, but is provided for convenience.
    ///
    /// # Arguments
    ///
    /// * `kernel`: the kernel to convolve with
    /// * `skip_unmasked`: if true, skip pixels in the raster which are unmasked (i.e. have a mask
    ///   value of 0). Unmasked pixels don't contribute to the convolution of other pixels (the
    ///   kernel is deweighted by the mask), but they can still have a valid convolution result if
    ///   the kernel overlaps with a masked pixel. If this is true, then the convolution will skip
    ///   those pixels entirely. Use this when the convolution doesn't have a meaningful result on
    ///   pixels that aren't considered part of the data set.
    ///
    /// returns: ScalarRaster
    pub fn convolve(&self, kernel: &RasterKernel, skip_unmasked: bool, keep_zlim: bool) -> Self {
        kernel.convolve(self, skip_unmasked, keep_zlim)
    }

    /// Performs a generic convolution operation on the depth image using the given kernel, but
    /// with a fast approximate method that will shrink the image and kernel down to a smaller
    /// target size, perform the convolution on the smaller raster, and then re-expand the
    /// result back to the full size of the original image. This is useful for operations that
    /// (1) need to occur in the scale of the world space represented by the raster, (2) tend to
    /// generate very large kernels which would be too expensive to compute at the full size, and
    /// (3) can tolerate some approximation in the result.
    ///
    /// # Arguments
    ///
    /// * `kernel`:
    ///
    /// returns: ScalarRaster
    pub fn convolve_fast(
        &self,
        kernel: &dyn FastApproxKernel,
        skip_unmasked: bool,
        keep_zlim: bool,
        target_size: usize,
    ) -> Result<Self> {
        let full_kernel = kernel.make(self.px_size)?;
        let shrink_factor = ((full_kernel.size as f64) / target_size as f64).floor();

        if shrink_factor + f64::EPSILON <= 1.0 {
            // If the shrink factor is less than or equal to 1, then we can just use the full
            // kernel and convolve the full image.
            return Ok(full_kernel.convolve(self, skip_unmasked, keep_zlim));
        }

        let shrunk_values = self.create_shrunk(shrink_factor as u32);
        let shrunk_kernel = kernel.make(shrunk_values.px_size)?;

        let shrunk_result = shrunk_kernel.convolve(&shrunk_values, skip_unmasked, keep_zlim);

        // The individual scale_x and scale_y factors are important because the scale may
        // not be equal after conversion to discrete pixel dimensions.
        let scaled_result = shrunk_result.create_scaled(shrink_factor);
        let scale_x = scaled_result.width() as f32 / self.width() as f32;
        let scale_y = scaled_result.height() as f32 / self.height() as f32;

        let mut full_result = ScalarRaster::empty_like(self);

        // All the pixels in the resized image which are at a full mask value can be copied
        // directly and put into the corresponding pixels in the convolved image
        let mut need_calc = Vec::new();
        for (x, y, v) in full_result.values.enumerate_pixels_mut() {
            // If we're skipping masked pixels, then we can skip this pixel if the self mask
            // value is not full.
            if skip_unmasked && !self.mask.get(x, y) {
                continue;
            }

            let sx = (x as f32 * scale_x).floor() as u32;
            let sy = (y as f32 * scale_y).floor() as u32;
            let scaled_mask = scaled_result.mask.buffer.get_pixel(sx, sy)[0];

            // If the scaled mask value is zero, then there was no valid convolution result
            // for this pixel, so we can skip it.
            if scaled_mask == 0 {
                continue;
            }

            // If the scaled result mask is not full, then we need to calculate this pixel
            // at full size later.
            if scaled_mask < 255 {
                need_calc.push((x, y));
                continue;
            }

            // Otherwise, we can just copy the value from the scaled result to the full result
            *v = Luma([scaled_result.values.get_pixel(
                (x as f32 * scale_x).floor() as u32,
                (y as f32 * scale_y).floor() as u32,
            )[0]]);

            // And set the mask
            full_result.mask.set(x, y, true);
        }

        let (target_matrix, target_mask) = self.to_value_and_mask_matrices();

        let final_calcs = need_calc
            .par_iter()
            .map(|(x, y)| {
                (
                    x,
                    y,
                    full_kernel.convolved_pixel_mat(
                        *x as usize,
                        *y as usize,
                        &target_matrix,
                        &target_mask,
                    ),
                )
            })
            .collect::<Vec<_>>();

        for (x, y, v) in final_calcs {
            full_result.set_f_at(*x as i32, *y as i32, v);
        }

        Ok(full_result)
    }

    /// Performs a smoothing operation on the depth image by averaging the pixels in a roughly
    /// circular neighborhood reachable along the object surface within approximately a given
    /// radius.
    pub fn area_blurred(&self, radius_mm: f32) -> Self {
        let optimal_scale = get_shrink_factor(radius_mm, self.px_size as f32, 16);

        // If the optimal scale is greater than 1 we will do a resizing operation and then
        // use the hybrid resized area average method. Otherwise, we will just do a normal
        // area average.
        let blurred = if optimal_scale > 1 {
            let small = self.create_shrunk(optimal_scale);
            let small_blurred = small.area_blurred(radius_mm);
            self.area_average_with_resized(radius_mm, &small_blurred)
        } else {
            self.area_average(radius_mm)
        };

        Self {
            values: blurred,
            mask: self.mask.clone(),
            px_size: self.px_size,
            min_z: self.min_z,
            max_z: self.max_z,
        }
    }

    /// Performs a high-pass filtering by taking an area blur at the given radius and then
    /// subtracting the result from the original image.
    pub fn area_filtered(&self, radius_mm: f32) -> Result<Self> {
        let blurred = self.area_blurred(radius_mm);
        self.subtract(&blurred)
    }

    pub fn count_valid_pixels(&self) -> usize {
        self.mask.count_true()
    }
}

/// Naive way of finding the optimal scale factor for a given target radius. Replace this with
/// direct division when you figure out the ceil/floor implications
fn get_shrink_factor(radius_mm: f32, mm_per_px: f32, target: u32) -> u32 {
    let mut shrink_factor = 1;
    let radius = (radius_mm / mm_per_px).ceil() as u32;
    while radius / shrink_factor > target {
        shrink_factor += 1;
    }

    shrink_factor
}

// pub fn mm_to_px_depth(value: f64) -> u16 {
//     // The zero value (0mm) is 32768
//     // The max value is +2mm, the min value is -2mm
//     let f = (value - (-DEPTH_MAX)) / (DEPTH_MAX - (-DEPTH_MAX));
//     (f * 65535.0) as u16
// }
//
// pub fn px_depth_to_mm(value: u16) -> f64 {
//     let f = value as f64 / 65535.0;
//     f * (DEPTH_MAX - (-DEPTH_MAX)) + (-DEPTH_MAX)
// }

pub fn diff_img(img1: &GrayImage, img2: &GrayImage) -> GrayImage {
    let mut diff = GrayImage::new(img1.width(), img1.height());

    for (p1, p2, p3) in diff.enumerate_pixels_mut() {
        let v1 = img1.get_pixel(p1, p2)[0];
        let v2 = img2.get_pixel(p1, p2)[0];

        if v1 >= v2 {
            *p3 = Luma([v1 - v2]);
        } else {
            *p3 = Luma([0]);
        }
    }

    diff
}

pub fn fill_cycle(mask: &GrayImage, count: usize) -> GrayImage {
    let mut dilated = mask.clone();

    for _ in 0..count {
        dilate_mut(&mut dilated, L1, 1);
    }

    for _ in 0..count {
        erode_mut(&mut dilated, L1, 1);
    }

    dilated
}

pub fn erode_cycle(mask: &GrayImage, count: usize) -> GrayImage {
    let mut dilated = mask.clone();

    for _ in 0..count {
        erode_mut(&mut dilated, L1, 1);
    }

    for _ in 0..count {
        dilate_mut(&mut dilated, L1, 1);
    }

    dilated
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::na::DMatrix;
    use approx::assert_relative_eq;

    #[test]
    fn scalar_raster_round_trip() {
        let expected = DMatrix::from_row_slice(
            3,
            4,
            &[
                -1.0, -2.0, -3.0, -4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
            ],
        );
        let raster = ScalarRaster::from_matrix(&expected, 1.0, -10.0, 20.0);
        let matrix = raster.to_matrix();

        for i in 0..expected.nrows() {
            for j in 0..expected.ncols() {
                assert_relative_eq!(expected[(i, j)], matrix[(i, j)], epsilon = 5e-4)
            }
        }
    }

    #[test]
    fn round_trip_serialization() {
        let m = 300;
        let n = 200;
        let mut values = DMatrix::from_fn(m, n, |i, j| {
            if (i + j) % 2 == 0 {
                (i * j) as f64
            } else {
                0.0
            }
        });
        values /= values.sum();

        let expected = ScalarRaster::from_matrix(&values, 1.0, -2.0, 2.0);

        let bytes = expected.serialized_bytes(ImageFormat::Png);
        let deserialized = ScalarRaster::from_serialized_bytes(&bytes).unwrap();

        for i in 0..m {
            for j in 0..n {
                let expected_value = expected.f_at(j as i32, i as i32);
                let actual_value = deserialized.f_at(j as i32, i as i32);
                if expected_value.is_nan() {
                    assert!(actual_value.is_nan());
                    continue;
                }
                assert_relative_eq!(expected_value, actual_value, epsilon = 1e-4);
            }
        }
    }
}
