use crate::Result;
use crate::common::dist;
use crate::geom2::Alignment2;
use crate::geom2::Point2;
use crate::geom2::align2::jacobian::{copy_jacobian, point_surf_jacobian2};
use crate::geom2::align2::{
    AlignOptions2, AlignParams2, AlignSurfMatch2, AlignValues2, SurfaceTarget2,
};
use crate::na::{Dyn, Matrix, Owned, U1, U3, Vector};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

/// Performs a Levenberg-Marquardt minimization to align a set of 2D points to a surface target.
///
/// The points are repeatedly projected onto their closest position on the target as the solver
/// moves them, so the correspondences are re-established at every step rather than being fixed up
/// front.
///
/// # Arguments
///
/// * `points`: the 2D points to be aligned, in their own local coordinate system
/// * `target`: an entity that implements the [`SurfaceTarget2`] trait, such as a `Curve2` or a
///   `Boundary2`, which will be the stationary target of the points.
/// * `params`: the alignment parameters, see [`AlignParams2`] for details
/// * `opts`: options controlling weighting and filtering, see [`AlignOptions2`]
///
/// returns: Result<Alignment<Unit<Complex<f64>>, 2>, Box<dyn Error, Global>>
pub fn points_to_surface2(
    points: &[Point2],
    target: &impl SurfaceTarget2,
    params: AlignParams2,
    opts: &AlignOptions2,
) -> Result<Alignment2> {
    let problem = PointsToSurface2::new(points, target, params, opts);

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        let residuals = result.residuals().unwrap().as_slice().to_vec();
        let c = result.params.current_values();
        Ok(Alignment2::new(
            c.transform,
            c.align,
            result.params.local,
            result.params.offset,
            residuals,
        ))
    } else {
        Err(format!(
            "Failed to align points to surface: {:?}",
            report.termination
        )
        .into())
    }
}

struct PointsToSurface2<'a, T: SurfaceTarget2> {
    points: &'a [Point2],
    target: &'a T,
    params: AlignParams2,

    /// The alignment values for the current parameters. Cached here so that `jacobian()` doesn't
    /// have to recompute what `move_points` already worked out.
    current: AlignValues2,

    /// The test points after being moved by the current transform, in the target's space.
    moved: Vec<Point2>,

    /// The closest match on the target for each moved point.
    closest: Vec<AlignSurfMatch2>,

    /// The signed distance from each moved point to its match.
    residuals: Vec<f64>,

    /// The weight applied to each point's residual and jacobian row.
    weights: Vec<f64>,

    ignore_off: bool,
}

impl<'a, T: SurfaceTarget2> PointsToSurface2<'a, T> {
    fn new(
        points: &'a [Point2],
        target: &'a T,
        params: AlignParams2,
        opts: &AlignOptions2,
    ) -> Self {
        let current = params.current_values();
        let mut x = Self {
            points,
            target,
            params,
            current,
            moved: vec![Point2::origin(); points.len()],
            closest: vec![AlignSurfMatch2::default(); points.len()],
            residuals: vec![0.0; points.len()],
            weights: vec![1.0; points.len()],
            ignore_off: opts.ignore_off,
        };

        x.move_points();
        x
    }

    /// Moves the points by the current transform, finds the closest position on the target to
    /// each, and recomputes the residuals and weights.
    fn move_points(&mut self) {
        self.current = self.params.current_values();
        let transform = self.current.transform;

        for i in 0..self.points.len() {
            let m = transform * self.points[i];
            let c = self.target.align_surf_closest_to(&m);

            // The residual is the distance between the test point and the closest point on the
            // target, adjusted for the direction of the scalar projection.
            self.residuals[i] = dist(&m, &c.point) * c.dn(&m).signum();

            self.weights[i] = if self.ignore_off {
                c.weight * f64::from(c.is_on)
            } else {
                c.weight
            };

            self.moved[i] = m;
            self.closest[i] = c;
        }
    }
}

impl<T: SurfaceTarget2> LeastSquaresProblem<f64, Dyn, U3> for PointsToSurface2<'_, T> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U3>;
    type ParameterStorage = Owned<f64, U3>;

    fn set_params(&mut self, x: &Vector<f64, U3, Self::ParameterStorage>) {
        self.params.set_storage(*x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, U3, Self::ParameterStorage> {
        self.params.get_storage()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for i in 0..self.points.len() {
            // The weights are currently binary, so folding them in directly is the same as
            // folding in their square root. Step 5 introduces continuous weights, at which point
            // this must become `sqrt(w)` to keep the effective weight on `r^2` equal to `w`.
            res[i] = self.residuals[i] * self.weights[i];
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U3, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, U3, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = point_surf_jacobian2(p, c, &self.current) * self.weights[i];
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::{mean_point, transform_points};
    use crate::geom2::align2::Dof3;
    use crate::geom2::{Boundary2, BoundaryData2, BoundaryEditor};
    use crate::{Curve2, Iso2};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    fn to_pts(v: &[(f64, f64)]) -> Vec<Point2> {
        v.iter().map(|(x, y)| Point2::new(*x, *y)).collect()
    }

    fn closed_curve(p: &[(f64, f64)]) -> Curve2 {
        Curve2::from_points(&to_pts(p), 1e-8, true).unwrap()
    }

    /// A closed CCW "stadium" boundary: two horizontal segments joined by semicircular arcs at
    /// each end. Exercises a `Boundary2` target containing genuine `Arc2` elements rather than
    /// only segments.
    fn stadium_boundary() -> Boundary2 {
        let mut data = BoundaryData2::new_closed();
        let mut cursor = data.get_cursor(None);
        // Bottom edge, left to right
        cursor.add_seg_xy(2.0, 0.0);
        // Right cap, counter-clockwise around (2.0, 0.5)
        cursor.add_arc_xy(2.0, 0.5, 2.0, 1.0, false);
        // Top edge, right to left
        cursor.add_seg_xy(0.0, 1.0);
        // Left cap, counter-clockwise around (0.0, 0.5)
        cursor.add_arc_xy(0.0, 0.5, 0.0, 0.0, false);
        data.try_to_boundary().unwrap()
    }

    #[test]
    fn curve_recovers_disturbance() -> Result<()> {
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

        let disturb = Iso2::new(crate::Vector2::new(0.05, -0.05), 10.0 * PI / 180.0);
        let moved = transform_points(&points, &disturb);

        let params = AlignParams2::new_at_origin(None);
        let result = points_to_surface2(&moved, &curve, params, &AlignOptions2::default())?;

        assert_relative_eq!(
            result.full().to_matrix(),
            disturb.inverse().to_matrix(),
            epsilon = 1e-10
        );
        Ok(())
    }

    #[test]
    fn curve_recovers_disturbance_at_center() -> Result<()> {
        // Same as above, but rotating about the centroid of the moved points rather than the
        // world origin. The recovered full transform must be identical either way.
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

        let disturb = Iso2::new(crate::Vector2::new(0.05, -0.05), 10.0 * PI / 180.0);
        let moved = transform_points(&points, &disturb);

        let params = AlignParams2::new_at_center(mean_point(&moved), None);
        let result = points_to_surface2(&moved, &curve, params, &AlignOptions2::default())?;

        assert_relative_eq!(
            result.full().to_matrix(),
            disturb.inverse().to_matrix(),
            epsilon = 1e-10
        );
        Ok(())
    }

    #[test]
    fn boundary_with_arcs_recovers_disturbance() -> Result<()> {
        let boundary = stadium_boundary();
        let points = boundary.to_points(1e-4)?;

        let disturb = Iso2::new(crate::Vector2::new(0.03, -0.02), 5.0 * PI / 180.0);
        let moved = transform_points(&points, &disturb);

        let params = AlignParams2::new_at_center(mean_point(&moved), None);
        let result = points_to_surface2(&moved, &boundary, params, &AlignOptions2::default())?;

        // The sampled points sit on the theoretical boundary, so the alignment should recover the
        // disturbance almost exactly. The tolerance is looser than the curve case because the
        // points are a chordal approximation of the arcs, not exact positions on them.
        assert_relative_eq!(
            result.full().to_matrix(),
            disturb.inverse().to_matrix(),
            epsilon = 1e-6
        );
        Ok(())
    }

    #[test]
    fn locked_tx_is_not_recovered() -> Result<()> {
        let curve = closed_curve(&[(0.0, 0.0), (5.0, 0.0), (5.0, 1.0), (0.0, 1.0)]);
        let points = to_pts(&[
            (1.0, 0.0),
            (2.0, 0.0),
            (3.0, 0.0),
            (1.0, 1.0),
            (2.0, 1.0),
            (3.0, 1.0),
            (0.0, 0.25),
            (0.0, 0.75),
            (5.0, 0.25),
            (5.0, 0.75),
        ]);

        // A pure x-translation, which the alignment is forbidden from undoing.
        let disturb = Iso2::translation(0.3, 0.0);
        let moved = transform_points(&points, &disturb);

        let dof = Dof3::new(false, true, true);
        let params = AlignParams2::new_at_origin(Some(dof));
        let result = points_to_surface2(&moved, &curve, params, &AlignOptions2::default())?;

        // With `local` and `offset` both identity, the full transform is exactly the alignment
        // transform, so a locked tx must leave the x translation at precisely zero.
        assert_eq!(result.full().translation.vector.x, 0.0);

        // ...and the disturbance genuinely was not recovered.
        let recovered = transform_points(&moved, result.full());
        let max_dev = recovered
            .iter()
            .zip(points.iter())
            .map(|(a, e)| (a - e).norm())
            .fold(0.0_f64, f64::max);
        assert!(
            max_dev > 0.1,
            "locked tx should have prevented recovery, max deviation was {max_dev}"
        );

        Ok(())
    }
}
