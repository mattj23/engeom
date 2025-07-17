//! This module has tools for working with kernels and convolutions on 2D raster data.

use crate::Result;
use crate::na::DMatrix;
use crate::raster2::{MaskOperations, MaskValue, ScalarRaster};
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
#[derive(Debug, Clone)]
pub struct RasterKernel {
    pub values: DMatrix<f64>,
    pub size: usize,
    pub tail_n: i32,
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
        Self::new(values).expect(
            "Failed to create Gaussian kernel, something is wrong with the \
        kernel matrix generation",
        )
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
        let tail_n = (size as i32 - 1) / 2; // The tail length is half the size minus one

        Ok(Self {
            values,
            size,
            tail_n,
        })
    }

    /// This function performs a normalized convolution (a convolution in which certain values are
    /// suppressed or skipped) of a kernel against a `ScalarRaster` entity. The values of the
    /// `ScalarRaster` are converted to a matrix of floating point pixel values alongside a mask
    /// matrix that replicates the meaning of the `ScalarRaster`'s mask. The matrices are then
    /// convolved with the kernel, and the result is re-encoded into a new `ScalarRaster`.
    ///
    /// The pixel size value of the resulting `ScalarRaster` is the same as the input raster,
    /// but the minimum and maximum z encoding values are recalculated based on the limits of the
    /// convolved product.
    ///
    /// The `skip_unmasked` argument determines whether to skip pixels that aren't in the input
    /// raster's mask. If it is set to true, all pixels that aren't in the input raster's mask will
    /// also not be in the output raster.  If set to false, the function will attempt to calculate
    /// all pixels in the raster, but those that are beyond the kernel's radius from the nearest
    /// will not have a valid convolution result and so will be excluded in the output raster mask.
    ///
    /// # Arguments
    ///
    /// * `raster`: The `ScalarRaster` containing the pixel values to be convolved.
    /// * `skip_unmasked`: If true, pixels that are not in the input raster's mask will be skipped
    ///   and no convolved value will be calculated for them. The output raster's mask will not
    ///   include these pixels.
    ///
    /// returns: ScalarRaster
    pub fn convolve(&self, raster: &ScalarRaster, skip_unmasked: bool) -> ScalarRaster {
        let (matrix, mask) = raster.to_value_and_mask_matrices();
        let convolved = self.convolve_matrix_and_mask(&matrix, &mask, skip_unmasked);
        let z_min = convolved.min();
        let z_max = convolved.max();

        ScalarRaster::from_matrix(&convolved, raster.px_size, z_min, z_max)
    }

    /// This function performs a normalized convolution (a convolution in which certain values are
    /// suppressed or skipped) of a kernel against a matrix of floating point pixel values.
    ///
    /// The function takes two `DMatrix<f64>` arguments: `matrix` and `mask`. The `matrix` contains
    /// the scalar pixel values to be convolved, while the `mask` contains the mask values that
    /// determines which pixels are valid for the convolution. The mask should be a matrix of the
    /// same size as the `matrix`, otherwise the function will panic.
    ///
    /// The values in `matrix` and in `mask` must be prepared according to the following rules:
    ///
    /// - Replace any NAN elements in the value matrix with 0.0
    /// - Create a mask matrix where all valid elements of the value matrix have a value of 1.0
    ///   in the corresponding position, and all invalid elements have a value of 0.0.
    ///
    /// During the convolution, the kernel's product with the mask will be summed alongside the
    /// typical kernel product with the pixel values.  The kernel/pixel product will then be
    /// divided by the kernel/mask product to normalize the result.
    ///
    /// In theory the `mask` matrix can be thought of as a matrix of weights, and the values in the
    /// mask can be any floating point value rather than just 0.0 and 1.0, allowing convolutions to
    /// include measures of certainty instead of just binary masking.  However, the `skip_unmasked`
    /// argument will then have a different meaning, and shouldn't be used.
    ///
    /// The `skip_unmasked` argument determines whether to skip pixels that have a value of 0.0
    /// (within the floating point epsilon) in the mask matrix. If `skip_unmasked` is true, pixels
    /// without a valid mask value won't be calculated, and the corresponding pixel in the
    /// result matrix will be a NAN. This is used for skipping pixels in cases where there isn't a
    /// conceptual meaning to a missing pixel...such as in a depth map where a missing pixel might
    /// mean there is no physical entity at that location, and so we don't care about the
    /// convolved value at that position.
    ///
    /// If `skip_unmasked` is false, all pixels will be calculated, though not all may have valid
    /// results and so some pixels may still be NAN in the result matrix. Any pixels that is beyond
    /// the kernel's radius from the nearest pixel with a nonzero mask will end up as a NAN in the
    /// final result.
    ///
    /// # Arguments
    ///
    /// * `matrix`: A `DMatrix<f64>` containing the pixel values of the matrix to be convolved.
    ///   Follow the guidelines above to prepare this matrix.
    /// * `mask`: A `DMatrix<f64>` containing the mask values, where valid pixels have a value of
    ///   1.0 and invalid pixels have a value of 0.0. This matrix must be the same size as `matrix`
    ///   and should also follow the guidelines above.
    /// * `skip_unmasked`: If true, pixels with a value of 0.0 in the mask matrix will be skipped
    ///   and no convolved value will be calculated for them. The value in the result matrix will
    ///   be a NAN where the mask value is 0.0.  If false, all pixels will be calculated, but not
    ///   all will have valid results.
    ///
    /// returns: Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>
    pub fn convolve_matrix_and_mask(
        &self,
        matrix: &DMatrix<f64>,
        mask: &DMatrix<f64>,
        skip_unmasked: bool,
    ) -> DMatrix<f64> {
        let mut result = DMatrix::zeros(matrix.nrows(), matrix.ncols());

        for y in 0..matrix.nrows() {
            for x in 0..matrix.ncols() {
                if skip_unmasked && mask[(y, x)].abs() < f64::EPSILON {
                    result[(y, x)] = f64::NAN;
                } else {
                    result[(y, x)] = self.convolved_pixel_mat(x, y, matrix, mask);
                }
            }
        }

        result
    }

    /// This function computes the value of a pixel at position `(x, y)` in a matrix by convolving
    /// it against the kernel values and a mask. This is the fundamental building block of the
    /// normalized convolutions performed by the `RasterKernel` struct.
    ///
    /// The function takes both a `DMatrix<f64>` representing the pixel values of the matrix and a
    /// `DMatrix<f64>` representing the mask values. These matrices must be prepared according to
    /// the following rules:
    ///
    /// - Replace any NAN elements in the value matrix with 0.0
    /// - Create a mask matrix where all valid elements of the value matrix have a value of 1.0
    ///   in the corresponding position, and all invalid elements have a value of 0.0.
    ///
    /// There is no similar function to operate directly on a `ScalarRaster` for performance
    /// reasons.  The following optimizations are what make this function fast enough to not be
    /// unusable:
    ///
    /// - The mask matrix is prepared in advance, and the valid elements of the value matrix have
    ///   a value of 1.0 in the corresponding position in the mask matrix. This allows the
    ///   normalization information to be collected in a single branchless multiplication operation
    /// - The invalid elements of the value matrix are set to zero, so there are no NANs to catch
    ///   with branches; the mask matrix replaces any if/then logic that would otherwise be
    ///   required.
    /// - The mask and value matrices are both pre-computed and stored externally, since
    ///   convolving the entire matrix will entail a LOT of repeated retrievals of the same values.
    ///   By computing these matrices in advance we only have to calculate the floating point
    ///   values a single time per convolution.
    ///
    /// # Arguments
    ///
    /// * `x`: The x-coordinate (column index, j) of the pixel in the matrix to be convolved.
    /// * `y`: The y-coordinate (row index, i) of the pixel in the matrix to be convolved.
    /// * `matrix`: The `DMatrix<f64>` containing the pixel values of the matrix to be convolved.
    /// * `mask`: A `DMatrix<f64>` containing the mask values, where valid pixels have a value of
    ///   1.0 and invalid pixels have a value of 0.0.
    ///
    /// returns: f64
    pub fn convolved_pixel_mat(
        &self,
        x: usize,
        y: usize,
        matrix: &DMatrix<f64>,
        mask: &DMatrix<f64>,
    ) -> f64 {
        let mut kernel_sum = 0.0;
        let mut pixel_sum = 0.0;

        // This initial arithmatic is to determine the mutually valid range of the kernel as it
        // lies over the pixel at (x, y) in the matrix.
        //
        // The `*_mi` and `*_mj` variables represent the valid range of indices in the matrix which
        // overlap with the kernel, while the `*_ki` and `*_kj` variables represent the
        // corresponding window indices in the kernel itself. Ultimately, all values will be
        // positive integers.
        // ----------------------------------------------------------------------------------------

        // For the row dimension (`i`/`y`), the first valid index in the matrix is either zero or
        // the `y` coordinate minus the tail length of the kernel, whichever is larger. The first
        // valid index in the kernel has the same offset from the center of the kernel as the first
        // valid index in the matrix has from the `y` coordinate.
        //
        // The last valid index in the matrix is either the `y` coordinate plus the tail length of
        // the kernel or the last index in the matrix, whichever is smaller. The total number of
        // valid indices in the kernel is the difference between the last and first valid indices
        // in the matrix, plus one to include the last index itself.
        let min_mi = (y as i32 - self.tail_n).max(0) as usize;
        let min_ki = min_mi - (y as i32 - self.tail_n) as usize;
        let max_mi = (y + self.tail_n as usize).min(matrix.nrows() - 1);
        let count_ki = max_mi - min_mi + 1;

        // The same logic applies to the column dimension (`j`/`x`), where the first and last valid
        // indices in the matrix are determined by the `x` coordinate, the tail length of the
        // kernel, and the limits of the matrix. The start of the kernel and the total count of
        // the kernel are calculated from the matrix limits.
        let min_mj = (x as i32 - self.tail_n).max(0) as usize;
        let min_kj = min_mj - (x as i32 - self.tail_n) as usize;
        let max_mj = (x + self.tail_n as usize).min(matrix.ncols() - 1);
        let count_kj = max_mj - min_mj + 1;

        // Now we will iterate over the valid indices in each direction, calculating both the
        // matrix and the kernel indices. From there we can accumulate the pixel and kernel sums.
        for i in 0..count_ki {
            for j in 0..count_kj {
                let ki = min_ki + i;
                let kj = min_kj + j;
                let mi = min_mi + i;
                let mj = min_mj + j;
                pixel_sum += matrix[(mi, mj)] * self.values[(ki, kj)];
                kernel_sum += mask[(mi, mj)] * self.values[(ki, kj)];
            }
        }

        if kernel_sum == 0.0 {
            f64::NAN
        } else {
            pixel_sum / kernel_sum
        }
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
        panic!("kernel total is zero, cannot normalize");
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
