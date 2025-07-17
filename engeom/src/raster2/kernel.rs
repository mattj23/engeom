//! This module has tools for working with kernels and convolutions on 2D raster data.

use crate::Result;
use crate::na::DMatrix;
use crate::raster2::{MaskOperations, ScalarRaster};
use rayon::iter::IntoParallelRefIterator;
use std::ops::Range;

pub trait FastApproxKernel {
    fn make(&self, pixel_size: f64) -> Result<RasterKernel>;
}

pub struct FastApproxGaussian {
    sigma: f64,
}

impl FastApproxGaussian {
    /// Creates a new `FastApproxGaussian` with the specified standard deviation (sigma) in world
    /// distance units.
    ///
    /// # Arguments
    ///
    /// * `sigma`: The standard deviation of the Gaussian distribution in pixels.
    /// returns: FastApproxGaussian
    pub fn new(sigma: f64) -> Self {
        Self { sigma }
    }
}

impl FastApproxKernel for FastApproxGaussian {
    /// Generates a Gaussian kernel based on the specified standard deviation (sigma) in world
    /// distance units. The kernel size is determined by the formula `size = ceil(sigma * 6)`, which
    /// ensures it is odd and large enough to capture the Gaussian distribution effectively.
    ///
    /// The kernel values are generated using the Gaussian function and then normalized so that
    /// all values sum to 1.0.
    ///
    /// # Arguments
    ///
    /// * `pixel_size`: The size of a pixel in world distance units, used to scale the sigma value.
    ///
    /// returns: DMatrix<f64>
    fn make(&self, pixel_size: f64) -> Result<RasterKernel> {
        let sigma = self.sigma / pixel_size;
        Ok(RasterKernel::gaussian(sigma))
    }
}

/// Represents a convolution kernel for raster data, which is a square matrix of an odd numbered
/// size full of f64 values which ideally sum to 1.0.  The `RasterKernel` struct can manage to
/// perform convolution operations on a `ScalarRaster` using the kernel values. It can be created
/// from any square, odd-sized matrix of f64 values, but some convenience methods are provided
/// to create typical kernels.
pub struct RasterKernel {
    pub values: DMatrix<f64>,
    pub size: usize,
}

impl RasterKernel {
    /// Generate a symmetrical Gaussian kernel with a given standard deviation (sigma) measured in
    /// pixels. The kernel size is determined by the formula `size = ceil(sigma * 6)`, ensuring it
    /// is odd and large enough to capture the Gaussian distribution effectively.
    ///
    /// The kernel values are generated using the Gaussian function and then normalized so that
    /// all values sum to 1.0.
    ///
    /// # Arguments
    ///
    /// * `sigma`: The standard deviation of the Gaussian distribution in pixels; this controls the
    ///   size and spread of the kernel.
    ///
    /// returns: Result<RasterKernel, Box<dyn Error, Global>>
    pub fn gaussian(sigma: f64) -> Self {
        let values = gaussian_kernel_matrix(sigma);
        let size = values.nrows();
        Self { values, size }
    }

    /// Creates a new `RasterKernel` from a 2D matrix of f64 values.
    pub fn new(values: DMatrix<f64>) -> Result<Self> {
        if values.nrows() == 0 || values.ncols() == 0 {
            return Err("Kernel must have non-zero dimensions".into());
        }
        if values.nrows() != values.ncols() {
            return Err("Kernel must be square (nrows == ncols)".into());
        }
        if values.nrows() % 2 == 0 {
            return Err("Kernel size must be odd (nrows and ncols must be odd)".into());
        }
        let size = values.nrows();

        Ok(Self { values, size })
    }

    // pub fn convolve(&self, raster: &ScalarRaster, skip_unmasked: bool) -> ScalarRaster {
    //     let mut result = ScalarRaster::empty_like(&raster);
    //
    //     for y in 0..raster.height() as i32 {
    //         for x in 0..raster.width() as i32 {
    //             if skip_unmasked && raster.mask.is_pixel_unmasked(x, y) {
    //                 result.mask.set_unmasked(x as u32, y as u32);
    //             } else {
    //                 let v = self.convolved_pixel(raster, x, y);
    //                 result.set_f_at(x, y, v);
    //             }
    //         }
    //     }
    //
    //     result
    // }
    pub fn convolve(&self, raster: &ScalarRaster, skip_unmasked: bool) -> ScalarRaster {
        let matrix = raster.to_matrix();
        let mut result = ScalarRaster::empty_like(&raster);

        for y in 0..raster.height() as i32 {
            for x in 0..raster.width() as i32 {
                if skip_unmasked && raster.mask.is_pixel_unmasked(x, y) {
                    result.mask.set_unmasked(x as u32, y as u32);
                } else {
                    let v = self.convolved_pixel_mat(raster, x, y, &matrix);
                    result.set_f_at(x, y, v);
                }
            }
        }

        result
    }
    pub fn convolved_pixel_mat(
        &self,
        raster: &ScalarRaster,
        x: i32,
        y: i32,
        matrix: &DMatrix<f64>,
    ) -> f64 {
        let pose = KernelPose::new(&self.values, raster, x, y);

        // Copy the kernel ignoring NAN values
        let mut kernel_copy = DMatrix::zeros(self.size, self.size);
        for c in pose.all() {
            if !c.mvf().is_nan() {
                kernel_copy[(c.ky, c.kx)] = c.kv();
            }
        }
        let total = kernel_copy.sum();
        if total == 0.0 {
            return f64::NAN; // Avoid division by zero
        }

        kernel_copy = kernel_copy / total;

        let mut sum = 0.0;
        for c in pose.all() {
            let mvf = matrix[(c.my as usize, c.mx as usize)];
            if !mvf.is_nan() {
                sum += mvf * kernel_copy[(c.ky, c.kx)];
            }
        }

        if sum.is_nan() {
            panic!("sum is nan");
        }

        sum
    }

    pub fn convolved_pixel(&self, raster: &ScalarRaster, x: i32, y: i32) -> f64 {
        let pose = KernelPose::new(&self.values, raster, x, y);

        // Copy the kernel ignoring NAN values
        let mut kernel_copy = DMatrix::zeros(self.size, self.size);
        for c in pose.all() {
            if !c.mvf().is_nan() {
                kernel_copy[(c.ky, c.kx)] = c.kv();
            }
        }
        let total = kernel_copy.sum();
        if total == 0.0 {
            return f64::NAN; // Avoid division by zero
        }

        kernel_copy = kernel_copy / total;

        let mut sum = 0.0;
        for c in pose.all() {
            if !c.mvf().is_nan() {
                sum += c.mvf() * kernel_copy[(c.ky, c.kx)];
            }
        }

        if sum.is_nan() {
            panic!("sum is nan");
        }

        sum
    }
}

/// A `KernelPose` represents the position of a convolution kernel on a `ScalarRaster`.  The pose
/// provides a method to iterate over all the kernel coordinates that overlap with valid indices
/// in the reference raster, taking into account the kernel's size and the raster's bounds.
struct KernelPose<'a> {
    /// A reference to the kernel matrix
    kernel: &'a DMatrix<f64>,

    /// A reference to the matrix being convolved
    reference: &'a ScalarRaster,

    /// The x (column) coordinate of the kernel's center in the reference matrix
    x: i32,

    /// The y (row) coordinate of the kernel's center in the reference matrix
    y: i32,

    x_bounds: KernelBounds,

    y_bounds: KernelBounds,

    radius: i32,
}

impl<'a> KernelPose<'a> {
    pub fn new(kernel: &'a DMatrix<f64>, reference: &'a ScalarRaster, x: i32, y: i32) -> Self {
        let radius = kernel.ncols() as i32 / 2;

        Self {
            kernel,
            reference,
            x,
            y,
            x_bounds: KernelBounds::from(x, radius, reference.width() as i32),
            y_bounds: KernelBounds::from(y, radius, reference.height() as i32),
            radius,
        }
    }

    pub fn all(&self) -> Vec<KCoords> {
        let mut result = Vec::new();
        for kx in self.x_bounds.kvals() {
            for ky in self.y_bounds.kvals() {
                let mx = self.x + kx as i32 - self.radius;
                let my = self.y + ky as i32 - self.radius;
                result.push(KCoords::new(kx, ky, mx, my, self.kernel, self.reference));
            }
        }
        result
    }
}

/// Represents a dual set of coordinates where kx, ky are the coordinates in the kernel and mx, my
/// are the corresponding coordinates in the ScalarRaster.  It is used to perform the mapping
/// between the kernel and the raster data for a single pixel during convolution operations.
struct KCoords<'a> {
    kx: usize,
    ky: usize,
    mx: i32,
    my: i32,

    kernel: &'a DMatrix<f64>,
    reference: &'a ScalarRaster,
}

impl<'a> KCoords<'a> {
    fn new(
        kx: usize,
        ky: usize,
        mx: i32,
        my: i32,
        kernel: &'a DMatrix<f64>,
        reference: &'a ScalarRaster,
    ) -> Self {
        Self {
            kx,
            ky,
            mx,
            my,
            kernel,
            reference,
        }
    }

    /// Returns the value of the kernel at the coordinates (kx, ky)
    pub fn kv(&self) -> f64 {
        self.kernel[(self.ky, self.kx)]
    }

    /// Returns the f64 value of the reference raster at the coordinates (mx, my). If the mask is
    /// not set, this will return NaN.
    pub fn mvf(&self) -> f64 {
        self.reference.f_at(self.mx, self.my)
    }

    /// Returns the u16 value of the reference raster at the coordinates (mx, my). If the mask is
    /// not set, this will return None.
    pub fn mv(&self) -> Option<u16> {
        self.reference.u_at(self.mx, self.my)
    }
}

/// This struct handles the boundaries of a single direction of the kernel in relation to the
/// reference raster being convolved. It makes iteration over the kernel coordinates transparent
/// to calling code when the kernel bounds might be outside the bounds of the raster.
struct KernelBounds {
    /// The number of valid kernel coordinates in this dimension starting from `k0`.
    count: usize,

    /// The starting index of the kernel in the kernel coordinates, such that both `k0` and `m0` are
    /// the first valid indices in the kernel and matrix respectively for the given dimension at a
    /// specific center coordinate of the kernel.
    ///
    /// For instance, if the kernel is 5x5 and kernel pose was placed at 0,0 in the raster, then
    /// the first valid matrix coordinates would be 0 and the first valid kernel coordinates
    /// would be 2.
    k0: usize,

    /// The starting index of the kernel in the matrix coordinates, such that both `k0` and `m0` are
    /// the first valid indices in the kernel and matrix respectively for the given dimension at a
    /// specific center coordinate of the kernel.
    ///
    /// For instance, if the kernel is 5x5 and kernel pose was placed at 0,0 in the raster, then
    /// the first valid matrix coordinates would be 0 and the first valid kernel coordinates
    /// would be 2.
    m0: usize,
}

impl KernelBounds {
    fn new(count: usize, k0: usize, m0: usize) -> Self {
        Self { count, k0, m0 }
    }

    fn from(x: i32, radius: i32, m_size: i32) -> Self {
        let m0 = x - radius; // The start of the kernel in the matrix coordinates ignoring bounds
        let m_lower = m0.max(0); // The start of the kernel in the matrix clipped to 0
        let k0 = m_lower - m0; // The start of the kernel in the kernel coordinates

        let m1 = x + radius; // The end of the kernel in the matrix coordinates ignoring bounds
        let m_upper = m1.min(m_size - 1); // The end of the kernel in the matrix clipped to the matrix size

        let count = (m_upper - m_lower + 1) as usize;
        Self::new(count, k0 as usize, m_lower as usize)
    }

    /// Returns a range of kernel coordinates for this dimension
    pub fn kvals(&self) -> Range<usize> {
        self.k0..self.k0 + self.count
    }
}

fn gaussian_kernel_matrix(sigma: f64) -> DMatrix<f64> {
    let size = (sigma * 6.0).ceil() as usize | 1; // Ensure size is odd

    let radius = size as f64 / 2.0;
    let mut values = DMatrix::zeros(size, size);

    for y in 0..size {
        for x in 0..size {
            let dx = x as f64 - radius;
            let dy = y as f64 - radius;
            let value = (-((dx * dx + dy * dy) / (2.0 * sigma * sigma))).exp();
            values[(y, x)] = value;
        }
    }

    // Normalize the kernel
    let total: f64 = values.sum();
    if total == 0.0 {
        // This shouldn't happen?
        panic!("RasterKernel::gaussian failed underflow: total sum is zero");
    }
    values / total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raster2::ScalarRaster;
    use approx::assert_relative_eq;

    fn obviously_correct_gaussian(target: &DMatrix<f64>) -> DMatrix<f64> {
        // Do an obviously correct convolution with a Gaussian kernel on a target where NANs are
        // the same as an unmasked value.
        let kernel = gaussian_kernel_matrix(1.0);
        let tail_n = (kernel.nrows() - 1) as i32 / 2i32;
        let n_rows = target.nrows() as i32;
        let n_cols = target.ncols() as i32;
        let mut result = DMatrix::zeros(target.nrows(), target.ncols());

        // i is rows and j is columns
        for i in 0..target.nrows() {
            for j in 0..target.ncols() {
                let mut kernel_sum = 0.0;
                let mut pixel_sum = 0.0;

                for ki in 0..kernel.nrows() {
                    for kj in 0..kernel.ncols() {
                        let mi = i as i32 + (ki as i32 - tail_n);
                        let mj = j as i32 + (kj as i32 - tail_n);

                        if mi >= 0 && mi < n_rows && mj >= 0 && mj < n_cols {
                            if target[(mi as usize, mj as usize)].is_nan() {
                                continue; // Skip masked pixels
                            } else {
                                pixel_sum += target[(mi as usize, mj as usize)] * kernel[(ki, kj)];
                                kernel_sum += kernel[(ki, kj)];
                            }
                        }
                    }
                }
                result[(i, j)] = pixel_sum / kernel_sum;
            }
        }

        result
    }

    fn striped_target() -> DMatrix<f64> {
        // Create a striped target matrix with some values
        let mut target = DMatrix::zeros(400, 600);
        for i in 0..target.nrows() {
            for j in 0..target.ncols() {
                target[(i, j)] = ((i + j) % 50) as f64 / 100.0;
            }
        }
        for i in 0..target.nrows() {
            for j in 0..target.ncols() {
                if (i + j) % 10 == 0 {
                    target[(i, j)] = f64::NAN; // Mask some pixels
                }
            }
        }
        target
    }

    #[test]
    fn gaussian_kernel() {
        let target = striped_target();
        let expected = obviously_correct_gaussian(&target);

        let max = target.max() * 2.0;
        let raster = ScalarRaster::from_matrix(&target, 1.0, -max, max);

        let kernel = RasterKernel::gaussian(1.0);
        let result = kernel.convolve(&raster, false);
        let result_matrix = result.to_matrix();

        assert_eq!(expected.nrows(), result_matrix.nrows());
        assert_eq!(expected.ncols(), result_matrix.ncols());

        for i in 0..expected.nrows() {
            for j in 0..expected.ncols() {
                assert_relative_eq!(expected[(i, j)], result_matrix[(i, j)], epsilon = 2e-4);
            }
        }
    }
}
