use crate::common::vec_f64::mean_and_stdev;
use parry3d_f64::na::Isometry;

/// A container for the results of an alignment operation, including the full transformation, the
/// various component transformations, and the residuals of the alignment.
pub struct Alignment<R, const D: usize> {
    full: Isometry<f64, R, D>,
    alignment: Isometry<f64, R, D>,
    local_origin: Isometry<f64, R, D>,
    offset: Isometry<f64, R, D>,
    residuals: Vec<f64>,
}

impl<R, const D: usize> Alignment<R, D> {
    pub(crate) fn new(
        full: Isometry<f64, R, D>,
        alignment: Isometry<f64, R, D>,
        local_origin: Isometry<f64, R, D>,
        offset: Isometry<f64, R, D>,
        residuals: Vec<f64>,
    ) -> Self {
        Self {
            full,
            alignment,
            local_origin,
            offset,
            residuals,
        }
    }

    /// Gets the full transformation, which is a single transformation that brings the test
    /// entity geometry directly to the target geometry. This is a composite of the local origin's
    /// inverse, the alignment, and the work offset ($W * A * L^{-1}$).
    pub fn full(&self) -> &Isometry<f64, R, D> {
        &self.full
    }

    /// Gets the alignment transformation, which is the transformation about the origin that is
    /// produced by the tx, ty, tz, rx, ry, rz parameters adjusted by the alignment algorithm.
    pub fn alignment(&self) -> &Isometry<f64, R, D> {
        &self.alignment
    }

    /// Gets the local origin, which is a transformation from world origin to the local origin
    /// used by the alignment algorithm.
    pub fn local_origin(&self) -> &Isometry<f64, R, D> {
        &self.local_origin
    }

    /// Gets the offset transform, which is a transformation applied to the test entity geometry
    /// after the alignment transformation. It is typically used to counteract the local origin, or
    /// to provide an initial guess.
    pub fn offset(&self) -> &Isometry<f64, R, D> {
        &self.offset
    }

    /// Gets the residuals of the alignment, which are the difference between the target geometry
    /// and the test entity geometry after the alignment transformation. The order is the same as
    /// the order of the target geometry.
    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }

    /// Calculates the average residual of the alignment.
    pub fn residual_mean(&self) -> f64 {
        self.residuals.iter().sum::<f64>() / self.residuals.len() as f64
    }

    /// Calculates the mean and standard deviation of the residuals of the alignment.
    pub fn residual_mean_std_dev(&self) -> (f64, f64) {
        mean_and_stdev(&self.residuals).unwrap()
    }
}

/// The result of an alignment operation, including the transform and the residuals
pub struct Align<R, const D: usize> {
    transform: Isometry<f64, R, D>,
    residuals: Vec<f64>,
}

impl<R, const D: usize> Align<R, D> {
    pub fn new(transform: Isometry<f64, R, D>, residuals: Vec<f64>) -> Self {
        Self {
            transform,
            residuals,
        }
    }

    pub fn transform(&self) -> &Isometry<f64, R, D> {
        &self.transform
    }

    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }

    pub fn avg_residual(&self) -> f64 {
        self.residuals.iter().sum::<f64>() / self.residuals.len() as f64
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DistMode {
    ToPoint,
    ToPlane,
}
