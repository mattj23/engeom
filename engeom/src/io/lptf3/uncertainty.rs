//! Uncertainty models for laser triangulation scanners
//!

use crate::io::lptf3::Lptf3UncertaintyModel;

pub struct DiffTanModel {
    detector_z_dist: f64,
    detector_y_dist: f64,
    scale: f64,
}

impl DiffTanModel {
    /// Create a model of triangulation uncertainty based on scaling the derivative of `tan(theta)`
    /// where `theta` is the angle between the direction of the laser beam and the direction of
    /// the detector at every measured point.
    ///
    /// This is the inherent triangulation depth uncertainty for a small angular uncertainty at
    /// the detector, so uncertainty scales nonlinearly with distance from the detector.
    ///
    /// Pay close attention to the meaning of the arguments, as they are not necessarily the same
    /// as the physical parameters of the scanner, but specify where the emitter is located from
    /// inside the scanner's data coordinate system, which assumes that Z+ is pointing in the
    /// direction of the emitter and Y+ is pointing orthogonal to the emitter plane.
    ///
    /// # Arguments
    ///
    /// * `detector_z_dist`: The distance from Z=0 in the scanner's data coordinate system to the
    ///   _optical center_ of the detector when modeled as a pinhole camera.
    /// * `detector_y_dist`: The distance from Y=0 in the scanner's data coordinate system to the
    ///   _optical center_ of the detector when modeled as a pinhole camera.
    /// * `scale`: A scale factor that is applied to the uncertainty value, which should then
    ///   produce the standard deviation of the expected distribution of z values if the same point
    ///   were to be repeatedly sampled.
    ///
    /// returns: DiffTanModel
    pub fn new(detector_z_dist: f64, detector_y_dist: f64, scale: f64) -> Self {
        Self {
            detector_z_dist,
            detector_y_dist,
            scale,
        }
    }
}

impl Lptf3UncertaintyModel for DiffTanModel {
    fn value(&self, _x: f64, z: f64) -> f64 {
        let theta = (self.detector_z_dist - z).atan2(self.detector_y_dist);
        self.scale / theta.cos().powi(2)
    }
}
