use parry3d_f64::na::Isometry;

/// The result of an alignment operation, including the transform and the residuals
pub struct Alignment<R, const D: usize> {
    transform: Isometry<f64, R, D>,
    residuals: Vec<f64>,
}

impl<R, const D: usize> Alignment<R, D> {
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

pub enum DistMode {
    ToPoint,
    ToPlane,
}
