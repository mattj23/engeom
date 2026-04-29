//! This module generalizes fitting of boundary geometry to sample points

use crate::common::points::dist;
use crate::geom2::Boundary2;
use crate::na::{DVector, Dyn, Matrix, Owned, Vector, U1};
use crate::{Point2, Result, SurfacePoint2};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use crate::common::VecDot;

pub type BndBuildFn = Box<dyn Fn(&DVector<f64>) -> Result<Boundary2>>;
const DELTA: f64 = 1e-6;

// =============================================================================================
// Fitting a boundary to discrete points
// =============================================================================================

/// Given a set of points, a function which takes a parameter vector and produces a [`Boundary2`],
/// and an initial guess, this function will attempt to perform a Levenberg-Marquardt best fit of
/// the parameters to produce a `Boundary2` that minimizes the residuals of the points projected
/// onto the boundary.
///
/// You will need to provide a `DVector` with a reasonable initial guess, and a [`BndBuildFn`] that
/// accepts a `DVector` of that size and produces a boundary. The function is allowed to return a
/// `Err(Box<dyn Error>)` if the parameters produce invalid geometry, but the initial guess values
/// provided may not fail or this function will exit before attempting the minimization.
///
/// Because the residuals are the values of the points projected onto the nearest point on the
/// boundary, be aware that there is no inherent pressure preventing the boundary geometry from
/// growing beyond the sample points _unless_ it exists within the nature of the geometry. For
/// example, if the points lie on a line and the boundary geometry is a segment, there is nothing
/// that will constrain the segment from growing beyond the ends of the points.
///
/// For this reason, it is best to write the builder function so that any constraints are inherently
/// contained within it. For example, write the builder function so that the segment ends are
/// only parameterized in one direction.
///
/// Refer to the examples for more information.
///
/// # Arguments
///
/// * `points`: A slice of `Point2` that the boundary will be fit to. The residuals are calculated
///   as the distance from each point to the closest point on the boundary, so the points will
///   "pull" the closest area of the boundary towards themselves, but won't have any effect on
///   the boundary farther away.
/// * `builder`: A function that takes a `DVector` of parameters and produces a `Boundary2`.
/// * `initial`: An initial guess for the parameters. The `builder` function must return an
///   `Ok(Boundary2)` with this input for the fitting algorithm to initialize. The algorithm will
///   use the size of the initial guess to determine the number of parameters, so make sure it
///   matches what the `builder` function is expecting.
/// * `ignore_ends`: if `true`, points that project onto the ends of an open boundary will have
///   residuals of zero
///
/// returns: Result<Matrix<f64, Dyn, Const<1>, VecStorage<f64, Dyn, Const<1>>>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
/// // This example will fit a boundary consisting of three line segments to a bunch of points
/// // arranged in a triangle.
/// // ------------------------------------------------------------------------------------------
/// use engeom::Point2;
/// use engeom::na::DVector;
/// use engeom::common::points::{to_points, fill_gaps};
/// use engeom::geom2::{BndBuildFn, fit_boundary_to_points, BoundaryData2};
/// use approx::assert_relative_eq;
///
/// // We'll create the initial corners and then use the `fill_gaps` helper to generate points
/// // between them. Normally this would be your observed data.
/// let corners = to_points(&[[1.0, 1.0], [3.0, 2.0], [2.0, 4.0], [1.0, 1.0]]);
/// let points = fill_gaps(&corners, 0.1);
///
/// // Here we define the function which creates the boundary from six parameters. This is an
/// // extremely simple parameterization that just encodes the three corners as x,y pairs
/// let builder: BndBuildFn = Box::new(|params: &DVector<f64>| {
///    let mut bdata = BoundaryData2::new_open(Point2::new(params[0], params[1]));
///    let mut cursor = bdata.get_cursor(None);
///    cursor.add_seg_xy(params[2], params[3]);
///    cursor.add_seg_xy(params[4], params[5]);
///    cursor.add_seg_xy(params[0], params[1]);
///    bdata.try_to_boundary()
/// });
///
/// // We'll provide an initial guess which has the corners in roughly the correct areas, just way
/// // too large for the actual data and not aligned with the edges. Then we'll run the fitting
/// // algorithm.
/// let initial = DVector::from(vec![0.0, 0.0, 4.0, 0.0, 1.0, 7.0]);
/// let result = fit_boundary_to_points(&points, &builder, initial, false).unwrap();
///
/// // Finally we'll verify that the corners match the ones we originally provided.
/// let expected = DVector::from(vec![1.0, 1.0, 3.0, 2.0, 2.0, 4.0]);
/// assert_relative_eq!(result, expected, epsilon = 1.0e-6);
///
/// ```
pub fn fit_boundary_to_points(
    points: &[Point2],
    builder: &BndBuildFn,
    initial: DVector<f64>,
    ignore_ends: bool,
) -> Result<DVector<f64>> {
    let fitting = BoundaryToPoints::new(points, ignore_ends);
    let problem = BoundaryFit::try_new(&fitting, builder, initial)?;
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        Ok(result.params)
    } else {
        Err(format!("Fitting failed: {:?}", report.termination).into())
    }
}

struct BoundaryToPoints<'a> {
    points: &'a [Point2],
    ignore_ends: bool,
}

impl<'a> BoundaryToPoints<'a> {
    fn new(points: &'a [Point2], ignore_ends: bool) -> Self {
        BoundaryToPoints { points, ignore_ends }
    }
}

impl BoundaryFittable for BoundaryToPoints<'_> {
    fn residuals_and_weights(&self, boundary: &Boundary2) -> (DVector<f64>, DVector<f64>) {
        let bounds = f64::EPSILON..(boundary.length() - f64::EPSILON);
        let mut res = DVector::zeros(self.points.len());
        let mut weights = DVector::zeros(self.points.len());
        weights.fill(1.0);

        for i in 0..self.points.len() {
            let (_, m) = boundary.at_closest_to_point(&self.points[i]);
            if self.ignore_ends {
                weights[i] = if bounds.contains(&m.l) { 1.0 } else { 0.0 };
            }
            res[i] = dist(&m, &self.points[i]);
        }

        (res, weights)
    }

    fn residual_only(&self, sample_i: usize, boundary: &Boundary2) -> f64 {
        let (_, m) = boundary.at_closest_to_point(&self.points[sample_i]);
        dist(&m, &self.points[sample_i])
    }
}

// =============================================================================================
// Fitting a boundary to discrete surface points
// =============================================================================================

/// Given a set of surface points, a function which takes a parameter vector and produces a
/// [`Boundary2`], and an initial guess, this function will attempt to perform a Levenberg-Marquardt
/// best fit of the parameters to produce a `Boundary2` that minimizes the residuals of the surface
/// points projected onto the boundary and weighted by their dot product with the boundary normal
/// at the projection site.
///
/// You will need to provide a `DVector` with a reasonable initial guess, and a [`BndBuildFn`] that
/// accepts a `DVector` of that size and produces a boundary. The function is allowed to return a
/// `Err(Box<dyn Error>)` if the parameters produce invalid geometry, but the initial guess values
/// provided may not fail or this function will exit before attempting the minimization.
///
/// Refer to the examples for more information.
pub fn fit_boundary_to_surface_points(
    points: &[SurfacePoint2],
    builder: &BndBuildFn,
    initial: DVector<f64>,
    weight_mode: VecDot,
) -> Result<DVector<f64>> {
    let fitting = BoundaryToSurfacePoints::new(points, weight_mode);
    let problem = BoundaryFit::try_new(&fitting, builder, initial)?;
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        Ok(result.params)
    } else {
        Err(format!("Fitting failed: {:?}", report.termination).into())
    }
}

struct BoundaryToSurfacePoints<'a> {
    points: &'a [SurfacePoint2],
    weight_mode: VecDot,
}

impl<'a> BoundaryToSurfacePoints<'a> {
    fn new(points: &'a [SurfacePoint2], weight_mode: VecDot) -> Self {
        BoundaryToSurfacePoints { points, weight_mode }
    }
}

impl BoundaryFittable for BoundaryToSurfacePoints<'_> {
    fn residuals_and_weights(&self, boundary: &Boundary2) -> (DVector<f64>, DVector<f64>) {
        let mut res = DVector::zeros(self.points.len());
        let mut weights = DVector::zeros(self.points.len());
        weights.fill(1.0);

        for i in 0..self.points.len() {
            let (_, m) = boundary.at_closest_to_point(&self.points[i]);
            let dot = self.points[i].normal.dot(&m.normal);
            weights[i] = match self.weight_mode {
                VecDot::AsIs => dot,
                VecDot::Abs => dot.abs(),
                VecDot::ClampPos => dot.max(0.0),
            };

            res[i] = dist(&m, &self.points[i]);
        }

        (res, weights)
    }

    fn residual_only(&self, sample_i: usize, boundary: &Boundary2) -> f64 {
        let (_, m) = boundary.at_closest_to_point(&self.points[sample_i]);
        dist(&m, &self.points[sample_i])
    }
}

// =============================================================================================
// `BoundaryFit` and `BoundaryFittable` together offer a generic mechanism for performing the
// fitting of boundaries to a set of samples of a finite count.
// =============================================================================================

pub trait BoundaryFittable {
    /// Given a boundary, this should return two equally sized `DVector`s of the residuals and
    /// the weights. The weights are not calculated independently during the jabcobian, but simply
    /// are reused from this step. Return the _un-weighted_ residuals, and the weights that will
    /// be used to scale them. For any element `i` in the residual vector, the value of
    /// `residual_only(i, ...)` should be identical.
    ///
    /// Return in the order: `(residuals, weights)`
    fn residuals_and_weights(&self, boundary: &Boundary2) -> (DVector<f64>, DVector<f64>);

    /// Given a boundary, this should return the residual for sample `sample_i`. This step is used
    /// for the numerical jacobians. It should produce values identical to the ones from
    /// `residuals_and_weights`
    fn residual_only(&self, sample_i: usize, boundary: &Boundary2) -> f64;
}

struct BoundaryFit<'a> {
    fitting: &'a dyn BoundaryFittable,
    params: DVector<f64>,
    builder: &'a BndBuildFn,
    current: Option<Boundary2>,
    residuals: Option<DVector<f64>>,
    weights: Option<DVector<f64>>,
}

impl<'a> BoundaryFit<'a> {
    fn try_new(
        fitting: &'a dyn BoundaryFittable,
        builder: &'a BndBuildFn,
        initial: DVector<f64>,
    ) -> Result<Self> {
        // Check to make sure that the initial value doesn't fail
        let _ = builder(&initial)?;

        let mut problem = BoundaryFit {
            fitting,
            params: DVector::zeros(initial.len()),
            builder,
            current: None,
            residuals: None,
            weights: None,
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
            let (residuals, weights) = self.fitting.residuals_and_weights(&boundary);
            self.residuals = Some(residuals);
            self.weights = Some(weights);
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
        let Some(weights) = &self.weights else {
            return None;
        };

        let mut jac = Matrix::<f64, Dyn, Dyn, Self::JacobianStorage>::zeros(
            residuals.len(),
            self.params.len(),
        );

        for k in 0..self.params.len() {
            let mut params = self.params.clone();
            params[k] += DELTA;
            let Ok(disturbed) = (self.builder)(&params) else {
                continue;
            };

            for i in 0..residuals.len() {
                let d = self.fitting.residual_only(i, &disturbed);
                jac[(i, k)] = weights[i] * (d - residuals[i]) / DELTA;
            }
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::{fill_gaps, to_points};
    use crate::geom2::BoundaryData2;
    use approx::assert_relative_eq;

    #[test]
    fn simple_triangle() {
        let corners = to_points(&[[1.0, 1.0], [3.0, 2.0], [2.0, 4.0], [1.0, 1.0]]);
        let points = fill_gaps(&corners, 0.1);

        let builder: BndBuildFn = Box::new(|params: &DVector<f64>| {
            let mut bdata = BoundaryData2::new_open(Point2::new(params[0], params[1]));
            let mut cursor = bdata.get_cursor(None);
            cursor.add_seg_xy(params[2], params[3]);
            cursor.add_seg_xy(params[4], params[5]);
            cursor.add_seg_xy(params[0], params[1]);
            bdata.try_to_boundary()
        });

        let initial = DVector::from(vec![0.0, 0.0, 4.0, 0.0, 1.0, 7.0]);
        let result = fit_boundary_to_points(&points, &builder, initial, false).unwrap();
        let expected = DVector::from(vec![1.0, 1.0, 3.0, 2.0, 2.0, 4.0]);
        assert_relative_eq!(result, expected, epsilon = 1.0e-6);
    }
}
