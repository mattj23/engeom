use crate::common::points::mean_point;
use crate::common::TolZone;
use crate::geom2::align2::jacobian::{copy_jacobian, point_surface_jacobian};
use crate::geom2::align2::RcParams2;
use crate::geom2::curve2::Curve2;
use crate::geom2::domain_tolerance_map::TolZoneMap;
use crate::geom2::{Align2, Iso2, Point2, SurfacePoint2};
use crate::Result;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry2d_f64::na::{Dyn, Matrix, Owned, Vector, U1, U3};

struct PointsToCurveZones<'a> {
    /// The original, unmodified points to be aligned.
    points: &'a [Point2],

    /// The curve which the points are being aligned to
    curve: &'a Curve2,

    /// The parameters which are being optimized, expressed as an isometry around a rotation center
    params: RcParams2,

    /// Storage for the internal points after they've been moved. This is constantly being updated
    /// whenever the parameters are changed.
    moved: Vec<Point2>,

    /// Storage for the closest point stations on the curve. This is constantly being updated
    /// whenever the parameters are changed.
    closest: Vec<SurfacePoint2>,

    /// The tolerance zone map
    tol_zone_map: &'a dyn TolZoneMap,

    /// Stored residuals for each point
    residuals: Vec<f64>,

    /// Out of zone
    scale: Vec<f64>,
}

impl<'a> PointsToCurveZones<'a> {
    fn new(
        points: &'a [Point2],
        curve: &'a Curve2,
        initial: &Iso2,
        tol_zone_map: &'a dyn TolZoneMap,
    ) -> Self {
        let mp = mean_point(points);
        let params = RcParams2::from_initial(initial, &mp);

        let mut item = Self {
            points,
            curve,
            params,
            moved: Vec::with_capacity(points.len()),
            closest: Vec::with_capacity(points.len()),
            tol_zone_map,
            residuals: Vec::with_capacity(points.len()),
            scale: Vec::with_capacity(points.len()),
        };
        item.move_points();
        item
    }

    fn move_points(&mut self) {
        self.moved.clear();
        self.closest.clear();
        self.residuals.clear();
        self.scale.clear();
        for p in self.points {
            let m = self.params.transform() * *p;
            let station = self.curve.at_closest_to_point(&m);
            let tol_zone = self.tol_zone_map.get(station.length_along());

            // We have to check if the point is within the tolerance zone so that the
            // residual can be scaled later. This check involves calculating the scalar projection
            // which is also used in the residual calculation, so it makes sense to do it once.
            let sp = station.surface_point();
            let residual = sp.scalar_projection(&m);

            self.scale.push(if tol_zone.contains(residual) {
                0.1
            } else {
                1.0
            });
            self.residuals.push(residual);
            self.closest.push(station.surface_point());
            self.moved.push(m);
        }
    }
}

/* LeastSquaresProblem
  ---------------------------------------------------------------------------------------------

  The way that the `levenberg_marquardt` crate works is by defining a LeastSquaresProblem trait
  through which the solver will interact with your implementation of the problem.  The solver will
  set the parameters, get the residuals, get the jacobian, and get the parameters. These
  operations will happen non-concurrently, allowing the implementation to be stateful and for
  concurrency (if desired) to be included in the problem implementation itself.

  The solver will work by calling `set_params` with parameters and then calling `residuals`
  and `jacobian` to extract these values.  Internally it will use these values to calculate the
  next step in the minimization process.  This will continue until the termination condition is
  met.

  To speed up the calculation of the residuals and the jacobian, we will pre-calculate any
  necessary common values when the parameters are set. Any stateful information which will be
  accessed more than once, and especially by both the residuals and the jacobian, should be
  calculated and stored in the implementation struct.  Then, when the `residuals` and `jacobian`
  functions are called, they can access this pre-calculated information.
*/

impl<'a> LeastSquaresProblem<f64, Dyn, U3> for PointsToCurveZones<'a> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U3>;
    type ParameterStorage = Owned<f64, U3>;

    fn set_params(&mut self, x: &Vector<f64, U3, Self::ParameterStorage>) {
        self.params.set(x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, U3, Self::ParameterStorage> {
        *self.params.x()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for (i, (r, s)) in self.residuals.iter().zip(self.scale.iter()).enumerate() {
            res[i] = r * s;
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U3, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, U3, Self::JacobianStorage>::zeros(self.points.len());

        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = point_surface_jacobian(p, c, &self.params) * self.scale[i];
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

/// Simple Levenberg-Marquardt minimization of the transform between a set of 2d points and a
/// `Curve2` polyline. This takes into account all points at all distances, projecting every point
/// in the given set onto its closest station on the curve, and minimizing the squared residual
/// distances.  Like any LM minimization, it will terminate on the first local minima it finds from
/// its starting position.
///
/// # Arguments
///
/// * `points`: A slice containing `Point2` entities to be aligned to the curve
/// * `curve`: A reference to the `Curve2` entity which will be the stationary reference
/// * `tol_zone_map`:
/// * `initial`: An initial starting position for the LM algorithm. Use Iso2::identity() to start
///     from the current point position.
///
/// returns: Result<Alignment<Unit<Complex<f64>>, 2>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn points_to_curve_zones(
    points: &[Point2],
    curve: &Curve2,
    tol_zone_map: &dyn TolZoneMap,
    initial: &Iso2,
) -> Result<Align2> {
    let problem = PointsToCurveZones::new(points, curve, initial, tol_zone_map);
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        let residuals = result.residuals().unwrap().as_slice().to_vec();
        Ok(Align2::new(*result.params.transform(), residuals))
    } else {
        let text = format!("Failed to align points to curve: {:?}", report.termination);
        Err(text.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::align2::iso2_from_param;
    use approx::assert_relative_eq;
    use parry3d_f64::na::Vector3;
    use std::f64::consts::PI;

    fn to_pts(v: &[(f64, f64)]) -> Vec<Point2> {
        v.iter().map(|(x, y)| Point2::new(*x, *y)).collect()
    }

    fn closed_curve(p: &[(f64, f64)]) -> Curve2 {
        let points = to_pts(p);
        Curve2::from_points(&points, 1e-8, true).unwrap()
    }

    #[test]
    fn test_align_points_to_curve() {
        let curve = closed_curve(&[(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (0.0, 1.0)]);

        let points = to_pts(&[
            (1.0, 0.0),
            (2.0, 0.0),
            (3.0, 0.0),
            (5.0, 0.25),
            (5.0, 0.75),
            (1.0, 1.0),
            (2.0, 1.0),
            (3.0, 1.0),
            (0.0, 0.25),
            (0.0, 0.75),
        ]);

        let shift = iso2_from_param(&Vector3::new(0.05, -0.05, 10.0 * PI / 180.0));
        let moved = points.iter().map(|p| shift * *p).collect::<Vec<_>>();

        let alignment = points_to_curve(&moved, &curve, &Iso2::identity()).unwrap();
        let result = alignment.transform().to_matrix();

        let shift_i = shift.inverse();

        assert_relative_eq!(result, shift_i.to_matrix(), epsilon = 1e-10);
    }
}
