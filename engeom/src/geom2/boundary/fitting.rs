//! This module generalizes fitting of boundary geometry to sample points

use crate::common::points::dist;
use crate::geom2::Boundary2;
use crate::na::{DVector, Dyn, Matrix, Owned, U1, Vector};
use crate::{Point2, Result};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

pub type BuildFn = Box<dyn Fn(&DVector<f64>) -> Result<Boundary2>>;
const DELTA: f64 = 1e-6;

pub fn fit_boundary_to_points(
    points: &[Point2],
    builder: &BuildFn,
    initial: DVector<f64>,
) -> Result<DVector<f64>> {
    let problem = BoundaryFit::try_new(points, builder, initial)?;
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        Ok(result.params)
    } else {
        Err(format!("Fitting failed: {:?}", report.termination).into())
    }
}

struct BoundaryFit<'a> {
    points: &'a [Point2],
    params: DVector<f64>,
    builder: &'a BuildFn,
    current: Option<Boundary2>,
    residuals: Option<DVector<f64>>,
}

impl<'a> BoundaryFit<'a> {
    fn try_new(points: &'a [Point2], builder: &'a BuildFn, initial: DVector<f64>) -> Result<Self> {
        // Check to make sure that the initial value doesn't fail
        let _ = builder(&initial)?;

        let mut problem = BoundaryFit {
            points,
            params: DVector::zeros(initial.len()),
            builder,
            current: None,
            residuals: None,
        };

        problem.set_params(&initial);

        Ok(problem)
    }
}

impl LeastSquaresProblem<f64, Dyn, Dyn> for BoundaryFit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, Dyn>;
    type ParameterStorage = Owned<f64, Dyn>;

    fn set_params(&mut self, x: &Vector<f64, Dyn, Self::ParameterStorage>) {
        self.params = x.clone();

        if let Ok(boundary) = (self.builder)(&self.params) {
            let mut res = DVector::zeros(self.points.len());
            for i in 0..self.points.len() {
                let p = boundary.at_closest_to_point(&self.points[i]).point;
                res[i] = dist(&p, &self.points[i]);
            }
            self.residuals = Some(res);
            self.current = Some(boundary);
        } else {
            self.residuals = None;
            self.current = None;
        }
    }

    fn params(&self) -> Vector<f64, Dyn, Self::ParameterStorage> {
        self.params.clone()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        self.residuals.clone()
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, Dyn, Self::JacobianStorage>> {
        let Some(residuals) = &self.residuals else {
            return None;
        };

        let mut jac = Matrix::<f64, Dyn, Dyn, Self::JacobianStorage>::zeros(
            residuals.len(), self.params.len(),
        );

        for k in 0..self.params.len() {
            let mut params = self.params.clone();
            params[k] += DELTA;
            let Ok(disturbed) = (self.builder)(&params) else {
                continue;
            };

            for i in 0..self.points.len() {
                let p = &self.points[i];
                let d = dist(p, &disturbed.at_closest_to_point(p).point);
                jac[(i, k)] = (d - residuals[i]) / DELTA;
            }
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use crate::common::points::{fill_gaps, to_points};
    use crate::geom2::BoundaryData2;
    use approx::assert_relative_eq;
    use super::*;

    #[test]
    fn simple_triangle() {
        let corners = to_points(&[[1.0, 1.0], [3.0, 2.0], [2.0, 4.0], [1.0, 1.0]]);
        let points = fill_gaps(&corners, 0.1);

        let builder: BuildFn = Box::new(|params: &DVector<f64>| {
            let mut bdata = BoundaryData2::new(Point2::new(params[0], params[1]));
            bdata.add_seg_xy(params[2], params[3]);
            bdata.add_seg_xy(params[4], params[5]);
            bdata.add_seg_xy(params[0], params[1]);
            bdata.try_to_boundary()
        });

        let initial = DVector::from(vec![0.0, 0.0, 4.0, 0.0, 1.0, 7.0]);
        let result = fit_boundary_to_points(&points, &builder, initial).unwrap();
        let expected = DVector::from(vec![1.0, 1.0, 3.0, 2.0, 2.0, 4.0]);
        assert_relative_eq!(result, expected, epsilon = 1.0e-6);
    }

}