//! This module contains common implementations for computing the values of the Jacobian matrix
//! for different 2D alignment Levenberg-Marquardt problems.

use crate::common::{PCoords, SPCoords};
use crate::geom2::align2::{AlignStorage2, AlignValues2};
use parry2d_f64::na::{Dim, Matrix, RawStorageMut, Storage, U3};

/// This is a helper function to calculate the partial derivatives of the parameters for a residual
/// distance between a test point and a surface point on a target entity.
///
/// This is the 2D counterpart of `geom3::align3::jacobian::point_surf_jacobian`, and is
/// deliberately structured identically to it.
///
/// # Arguments
///
/// * `p`: the test point, transformed into the target entity's coordinate system
/// * `c`: the closest point on the target surface
/// * `align`: the current alignment values, allowing for fast access of the direction vectors
///   associated with the different partial differentials.
///
/// returns: Matrix<f64, Const<3>, Const<1>, ArrayStorage<f64, 3, 1>>
pub fn point_surf_jacobian2(
    p: &impl PCoords<2>,
    c: &impl SPCoords<2>,
    align: &AlignValues2,
) -> AlignStorage2 {
    let mut result = AlignStorage2::zeros();

    // We'll grab the sign of the scalar projection, allowing us to know if we're outside or inside
    // the target surface.
    let sign = c.scalar_projection(p).signum();

    // First, we want to calculate the deviation direction. If the magnitude of the deviation is
    // close to zero, we'll use the surface normal instead, as the two will converge as the point
    // approaches the surface.
    let dev = p.coords() - c.coords();
    let dir = if dev.norm_squared() < 1e-16 {
        c.normal().into_inner()
    } else {
        dev.normalize() * sign
    };

    // The translations will be the dot product of the deviation direction and the partial
    // differential translation directions. For instance, if the translation is going across the
    // surface, there's no real penalty because we're staying equidistant.
    result[0] = val_or_zero(align.dtx.dot(&dir), align.dof.tx);
    result[1] = val_or_zero(align.dty.dot(&dir), align.dof.ty);

    // The rotation will be the dot product of the deviation direction and the partial differential
    // rotation direction.
    //
    // This is evaluated at the test point `p` because that is what the residual's derivative
    // literally calls for: it is `p` that the parameters move, while `c` is a fixed position on
    // the stationary target. Evaluating at `c` instead, as `align3` does, gives an identical
    // result, and it is worth recording why rather than leaving it to look like a discrepancy.
    //
    // `pre_rot`'s rotation part is the inverse of `offset`'s rotation, and `m_drz` is `offset`'s
    // rotation times `SK2`, so the translation cancels in the difference and
    //
    //     drz(p) - drz(c) = post_rot * SK2 * post_rot^-1 * (p - c)
    //
    // Conjugating a skew-symmetric matrix by a rotation leaves it skew-symmetric (in 2D, since
    // rotations commute with `SK2`, it is just `SK2` again), and `v' * S * v == 0` for every
    // skew-symmetric `S` and every vector `v`. Since `dir` is parallel to `p - c`, the difference
    // contributes exactly zero. The same argument covers 3D, where the Euler partials are built
    // as `Q * SK * Q'` and so are skew-symmetric too.
    result[2] = val_or_zero(align.drz(p).dot(&dir), align.dof.rz);

    result
}

/// The counterpart to [`point_surf_jacobian2`] for the case where it is the *target* entity whose
/// transform is being optimized, rather than the test point's.
///
/// This is what a multi-body adjustment needs. When two measured curves are aligned to each other,
/// a correspondence between them constrains both bodies: moving the test curve slides `p`, and
/// moving the reference curve slides `c`. Each contributes a block to the same jacobian row, and
/// this function supplies the second one.
///
/// The residual is the same signed point-to-point distance that [`point_surf_jacobian2`]
/// differentiates, so the two differ only in which point the parameters move and therefore in
/// sign: displacing the target by `v` changes the distance by `-dir . v` where `dir` is the unit
/// deviation direction.
///
/// The rotation partial is evaluated at `c`, which is the point this function's parameters
/// actually move. As in the forward case, the choice turns out not to matter: the skew-symmetry
/// argument in [`point_surf_jacobian2`] applies unchanged here, since the difference between
/// evaluating at `p` and at `c` is still orthogonal to a `dir` that is still parallel to `p - c`.
/// `stress_surf_rev_against_numeric` passes either way. `c` is used because it is the point the
/// derivative is actually taken with respect to, not because the other form is wrong.
///
/// This is the 2D counterpart of `geom3::align3::jacobian::point_surf_jacobian_rev`, and is
/// deliberately structured identically to it.
///
/// # Arguments
///
/// * `p`: the test point, in the common coordinate system the residual is measured in
/// * `c`: the closest point on the target surface, in that same coordinate system, meaning it has
///   already been moved by the target's current transform
/// * `align`: the current alignment values of the **target** entity
///
/// returns: Matrix<f64, Const<3>, Const<1>, ArrayStorage<f64, 3, 1>>
pub fn point_surf_jacobian2_rev(
    p: &impl PCoords<2>,
    c: &impl SPCoords<2>,
    align: &AlignValues2,
) -> AlignStorage2 {
    let mut result = AlignStorage2::zeros();

    let sign = c.scalar_projection(p).signum();

    let dev = p.coords() - c.coords();
    let dir = if dev.norm_squared() < 1e-16 {
        c.normal().into_inner()
    } else {
        dev.normalize() * sign
    };

    // The target moving away from the test point closes the same distance that the test point
    // moving toward the target would, hence the negation.
    let dir = -dir;

    result[0] = val_or_zero(align.dtx.dot(&dir), align.dof.tx);
    result[1] = val_or_zero(align.dty.dot(&dir), align.dof.ty);

    result[2] = val_or_zero(align.drz(c).dot(&dir), align.dof.rz);

    result
}

fn val_or_zero(value: f64, condition: bool) -> f64 {
    if condition { value } else { 0.0 }
}

/// Generic helper to copy the contents of a single row into a larger jacobian matrix of either
/// fixed or dynamic row count
///
/// # Arguments
///
/// * `j`: The source Jacobian row to copy
/// * `matrix`: The destination matrix to copy into
/// * `row`: The row index in the destination matrix to copy into
///
/// returns: ()
pub fn copy_jacobian<R, S>(j: &AlignStorage2, matrix: &mut Matrix<f64, R, U3, S>, row: usize)
where
    R: Dim,
    S: RawStorageMut<f64, R, U3> + Storage<f64, R, U3>,
{
    matrix.row_mut(row).copy_from_slice(j.as_slice());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::dist;
    use crate::common::random_geometry::RandomGeometry2;
    use crate::geom2::align2::{AlignOrigin2, AlignParams2, AlignSurfMatch2};
    use crate::geom2::{Point2, SurfacePoint2, UnitVec2, Vector2};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    const NUMERIC_EPS: f64 = 1e-7;

    /// The residual that `point_surf_jacobian2` differentiates: the signed distance from the test
    /// point to its match, with the sign taken from which side of the target's normal it lies on.
    ///
    /// This must stay in step with the residual computed in `points_to_surface.rs`; the jacobian
    /// is only correct with respect to this exact expression.
    fn residual(p: &Point2, c: &AlignSurfMatch2) -> f64 {
        dist(p, &c.point) * c.scalar_projection(p).signum()
    }

    /// Builds a match sitting at `p + offset`, with its normal aimed back along the offset either
    /// toward the test point or away from it.
    fn make_match(p: &Point2, offset: Vector2, flip: bool) -> AlignSurfMatch2 {
        let c_point = p + offset;
        let sign = if flip { -1.0 } else { 1.0 };
        let normal = UnitVec2::new_normalize((p - c_point) * sign);
        AlignSurfMatch2::new(c_point, normal, true, 1.0)
    }

    /// A central-difference estimate of the partial derivative of the residual with respect to one
    /// parameter, with the correspondence `c` held fixed.
    ///
    /// Holding `c` fixed is deliberate and matches what the analytic jacobian models: the solver
    /// re-establishes correspondences between steps, but within a single linearization the match
    /// is a fixed position on the stationary target.
    fn numeric(params: &AlignParams2, p_local: &Point2, c: &AlignSurfMatch2, index: usize) -> f64 {
        let stored = params.storage();
        let mut w0 = params.clone();
        let mut w1 = params.clone();
        w0.set_index(index, stored[index] - NUMERIC_EPS);
        w1.set_index(index, stored[index] + NUMERIC_EPS);

        let p0 = w0.compute_values().transform * p_local;
        let p1 = w1.compute_values().transform * p_local;

        (residual(&p1, c) - residual(&p0, c)) / (2.0 * NUMERIC_EPS)
    }

    #[test]
    fn partials_at_identity() {
        // The test point sits one unit in +x from its match, so the deviation direction is +x.
        // Translating in +x pulls it further away, translating in +y slides it across, and
        // rotating about the origin moves it perpendicular to the deviation.
        let params = AlignParams2::from_origin(None);
        let p = Point2::new(1.0, 0.0);
        let c = make_match(&p, Vector2::new(-1.0, 0.0), false);

        let j = point_surf_jacobian2(&p, &c, &params.compute_values());

        assert_relative_eq!(j.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(j.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(j.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn flipped_normal_flips_the_sign() {
        // Same geometry, but with the target normal pointing away from the test point, which puts
        // the point on the negative side and reverses every partial.
        let params = AlignParams2::from_origin(None);
        let p = Point2::new(1.0, 0.0);
        let c = make_match(&p, Vector2::new(-1.0, 0.0), true);

        let j = point_surf_jacobian2(&p, &c, &params.compute_values());

        assert_relative_eq!(j.x, -1.0, epsilon = 1e-10);
    }

    #[test]
    fn locked_dof_zeroes_its_column() {
        let dof = crate::geom2::align2::Dof3::new(false, true, false);
        let params = AlignParams2::from_origin(Some(dof));
        let p = Point2::new(1.0, 2.0);
        let c = make_match(&p, Vector2::new(-0.5, -0.3), false);

        let j = point_surf_jacobian2(&p, &c, &params.compute_values());

        assert_eq!(j.x, 0.0);
        assert_eq!(j.z, 0.0);
        assert_ne!(j.y, 0.0);
    }

    #[test]
    fn stress_against_numeric() {
        // The replacement for the old `stress_point_surf_against_numeric`, rewritten against the
        // signed point-to-point residual that `points_to_surface2` actually minimizes rather than
        // the point-to-plane residual the previous generation used.
        //
        // The deviation magnitudes here are deliberately large (up to 2 units). The distinction
        // between evaluating the rotation partial at the test point and at the match is of order
        // `|p - c|`, so a test that only sampled near-converged configurations would pass either
        // way and prove nothing.
        let mut rg = RandomGeometry2::new();
        for _ in 0..10000 {
            let local = rg.iso2(10.0);
            let offset = rg.iso2(10.0);
            let params = AlignParams2::new(AlignOrigin2::Local(local), Some(offset), None)
                .with_storage(rg.vector(PI));

            let p_local = rg.point(10.0);
            let p = params.compute_values().transform * p_local;

            // A deviation with a guaranteed non-trivial magnitude, so the match never degenerates
            // onto the test point.
            let dev = rg.unit_vec().into_inner() * rg.f64(0.1, 2.0);
            let c = make_match(&p, dev, rg.bool());

            let expected = [
                numeric(&params, &p_local, &c, 0),
                numeric(&params, &p_local, &c, 1),
                numeric(&params, &p_local, &c, 2),
            ];

            let j = point_surf_jacobian2(&p, &c, &params.compute_values());

            assert_relative_eq!(j.x, expected[0], epsilon = 1e-5);
            assert_relative_eq!(j.y, expected[1], epsilon = 1e-5);
            assert_relative_eq!(j.z, expected[2], epsilon = 1e-5);
        }
    }

    // ============================================================================================
    // The reverse point-to-surface jacobian, checked against finite differences
    // ============================================================================================

    /// The residual `point_surf_jacobian2_rev` differentiates: the signed point-to-point distance
    /// from a fixed test point to a target point which the target's own transform moves.
    ///
    /// `c_local` and `n_local` describe the match in the target's own coordinates, so that
    /// perturbing the target's parameters moves it as it would during a real solve.
    fn surf_rev_residual(
        params: &AlignParams2,
        p: &Point2,
        c_local: &Point2,
        n_local: &UnitVec2,
    ) -> f64 {
        let sp = SurfacePoint2::new(*c_local, *n_local).transformed_by(&params.compute_transform());
        dist(p, &sp.point) * sp.scalar_projection(p).signum()
    }

    /// A finite-difference estimate of the partial derivative of that residual with respect to one
    /// of the target's three parameters.
    fn surf_rev_numeric(
        params: &AlignParams2,
        p: &Point2,
        c_local: &Point2,
        n_local: &UnitVec2,
        index: usize,
    ) -> f64 {
        let mut lo = params.clone();
        lo.set_index(index, params.storage()[index] - NUMERIC_EPS);
        let mut hi = params.clone();
        hi.set_index(index, params.storage()[index] + NUMERIC_EPS);

        let r_lo = surf_rev_residual(&lo, p, c_local, n_local);
        let r_hi = surf_rev_residual(&hi, p, c_local, n_local);
        (r_hi - r_lo) / (2.0 * NUMERIC_EPS)
    }

    #[test]
    fn stress_surf_rev_against_numeric() {
        // The analytic reverse jacobian against a central finite difference, over random local
        // origins, working offsets, parameter vectors, and geometry. This is the gate on the
        // reverse form: the sign convention and the requirement that the rotation partial be
        // evaluated at `c` are both easy to get wrong by inspection and obvious here.
        //
        // Seeded so a failure is reproducible and so this can never join the flaky-by-RNG set.
        let mut rg = RandomGeometry2::from_seed(0x5eed_a11c);

        for _ in 0..2000 {
            let params =
                AlignParams2::new(AlignOrigin2::Local(rg.iso2(5.0)), Some(rg.iso2(5.0)), None)
                    .with_storage(AlignStorage2::new(
                        rg.f64_sym(3.0),
                        rg.f64_sym(3.0),
                        rg.f64_sym(0.8),
                    ));

            let c_local = rg.point(5.0);
            let n_local = rg.unit_vec();

            // Keep the test point well away from the match, so the deviation direction is well
            // conditioned and the residual is differentiable.
            let sp =
                SurfacePoint2::new(c_local, n_local).transformed_by(&params.compute_transform());
            let p = sp.point + sp.normal.into_inner() * rg.f64_sym(4.0).abs().max(0.5);

            let analytic = point_surf_jacobian2_rev(&p, &sp, &params.compute_values());

            for i in 0..3 {
                let numeric = surf_rev_numeric(&params, &p, &c_local, &n_local, i);
                assert_relative_eq!(analytic[i], numeric, epsilon = 1e-6);
            }
        }
    }

    #[test]
    fn surf_rev_is_the_negative_of_the_forward_form_for_translations() {
        // Sliding the target by `v` and sliding the test point by `-v` change the distance
        // identically, so the translation partials of the two forms must be exact negatives. The
        // rotation partials are not related this way, because they pivot about different points.
        let params = AlignParams2::from_origin(None);
        let values = params.compute_values();

        let c = Point2::new(1.0, 2.0);
        let n = UnitVec2::new_normalize(Vector2::new(0.0, 1.0));
        let p = Point2::new(1.0, 3.5);
        let sp = SurfacePoint2::new(c, n);

        let fwd = point_surf_jacobian2(&p, &sp, &values);
        let rev = point_surf_jacobian2_rev(&p, &sp, &values);

        for i in 0..2 {
            assert_relative_eq!(fwd[i], -rev[i], epsilon = 1e-12);
        }
    }

    #[test]
    fn surf_rev_respects_locked_dof() {
        let dof = crate::geom2::align2::Dof3::new(true, false, false);
        let params = AlignParams2::from_origin(Some(dof));
        let values = params.compute_values();

        let c = Point2::new(1.0, 2.0);
        let n = UnitVec2::new_normalize(Vector2::new(0.0, 1.0));
        let p = Point2::new(2.0, 3.0);

        let rev = point_surf_jacobian2_rev(&p, &SurfacePoint2::new(c, n), &values);

        for i in [1usize, 2] {
            assert_eq!(
                rev[i], 0.0,
                "locked parameter {i} should have a zero column"
            );
        }
        assert_ne!(rev[0], 0.0);
    }
}
