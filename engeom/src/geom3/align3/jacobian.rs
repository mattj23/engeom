//! Utilities for building Jacobian rows for 3d geometric optimization problems.
//!
//! A **Jacobian** is a matrix of partial derivatives. In this module, each row describes how a
//! distance-like residual changes when the underlying transform parameters are perturbed. These
//! derivatives are especially useful in iterative solvers such as least-squares optimization,
//! where the Jacobian helps the solver decide how to update position and orientation parameters.
//!
//! The functionality in this module focuses on common geometric distance models used throughout the
//! crate:
//!
//! - Point-to-plane approximations, where a point is compared against a nearby surface point and
//!   normal
//! - A reverse point-to-plane form for cases where the reference entity is transformed instead of
//!   the test point, this is used for multi-mesh alignments
//! - Point-to-point distance derivatives
//! - Helpers for inserting a single Jacobian row into a larger matrix
//!
//! These helpers are intended to make it easy to assemble Jacobians for alignment, registration,
//! and fitting routines in a consistent 6-parameter rigid-body setting.

use super::*;
use crate::geom3::align3::params::AlignValues3;
use parry2d_f64::na::Dim;
use parry3d_f64::na::{Matrix, RawStorageMut, Storage, U6};

/// This is a helper function to calculate the partial derivatives of the parameters for a residual
/// distance between a test point and a surface point on a target entity.
///
/// # Arguments
///
/// * `p`: the test point, transformed into the target entity's coordinate system
/// * `c`: the closest point on the target surface
/// * `align`: the current alignment parameters, allowing for fast access of the direction vectors
///   associated with the different partial differentials.
///
/// returns: Matrix<f64, Const<6>, Const<1>, ArrayStorage<f64, 6, 1>>
pub fn point_surf_jacobian(
    p: &impl PCoords<3>,
    c: &impl SPCoords<3>,
    align: &AlignValues3,
) -> AlignStorage3 {
    let mut result = AlignStorage3::zeros();

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

    // The transformations will be the dot product of the deviation direction and the partial
    // differential translation directions. For instance, if the translation is going across the
    // surface, there's no real penalty because we're staying equidistant.
    result[0] = val_or_zero(align.dtx.dot(&dir), align.dof.tx);
    result[1] = val_or_zero(align.dty.dot(&dir), align.dof.ty);
    result[2] = val_or_zero(align.dtz.dot(&dir), align.dof.tz);

    // The rotations will be the dot product of the deviation direction and the partial differential
    // rotation directions.
    //
    // These are evaluated at the closest point `c` rather than the test point `p`, even though
    // the residual's derivative nominally calls for `p` (it is `p` that the parameters move, `c`
    // being a fixed position on the stationary target). The two are exactly equal, so the choice
    // is immaterial:
    //
    // `pre_rot`'s rotation part is the inverse of `offset`'s rotation and each `m_dr*` is
    // `offset`'s rotation times an Euler partial, so the translation cancels in the difference
    // and `dr(p) - dr(c) = post_rot * E * post_rot^-1 * (p - c)`. `euler_partials` builds every
    // `E` as `Q * SK * Q'` for some rotation `Q`, which is skew-symmetric, and conjugating by
    // `post_rot` keeps it so. Since `v' * S * v == 0` for any skew-symmetric `S`, and `dir` is
    // parallel to `p - c`, the difference contributes exactly zero.
    //
    // The 2D counterpart in `geom2::align2::jacobian` uses `p` instead, purely because it reads
    // more directly as the derivative; its `stress_against_numeric` test confirms both forms
    // agree against a finite-difference estimate far from convergence.
    result[3] = val_or_zero(align.drx(c).dot(&dir), align.dof.rx);
    result[4] = val_or_zero(align.dry(c).dot(&dir), align.dof.ry);
    result[5] = val_or_zero(align.drz(c).dot(&dir), align.dof.rz);

    result
}

/// The counterpart to [`point_surf_jacobian`] for the case where it is the *target* entity whose
/// transform is being optimized, rather than the test point's.
///
/// This is what a multi-body adjustment needs. When two measured meshes are aligned to each other,
/// a correspondence between them constrains both bodies: moving the test mesh slides `p`, and
/// moving the reference mesh slides `c`. Each contributes a block to the same jacobian row, and
/// this function supplies the second one.
///
/// The residual is the same signed point-to-point distance that [`point_surf_jacobian`]
/// differentiates, so the two differ only in which point the parameters move and therefore in
/// sign: displacing the target by `v` changes the distance by `-dir . v` where `dir` is the unit
/// deviation direction.
///
/// The rotation partials are evaluated at `c`, which is the point this function's parameters
/// actually move. As in the forward case, though, the choice turns out not to matter: the
/// skew-symmetry argument in [`point_surf_jacobian`] applies unchanged here, since the difference
/// between evaluating at `p` and at `c` is still orthogonal to a `dir` that is still parallel to
/// `p - c`. `stress_surf_rev_against_numeric` passes either way. `c` is used because it reads as
/// the more honest derivative, not because the other form is wrong.
///
/// # Arguments
///
/// * `p`: the test point, in the common coordinate system the residual is measured in
/// * `c`: the closest point on the target surface, in that same coordinate system, meaning it has
///   already been moved by the target's current transform
/// * `align`: the current alignment values of the **target** entity
///
/// returns: Matrix<f64, Const<6>, Const<1>, ArrayStorage<f64, 6, 1>>
pub fn point_surf_jacobian_rev(
    p: &impl PCoords<3>,
    c: &impl SPCoords<3>,
    align: &AlignValues3,
) -> AlignStorage3 {
    let mut result = AlignStorage3::zeros();

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
    result[2] = val_or_zero(align.dtz.dot(&dir), align.dof.tz);

    result[3] = val_or_zero(align.drx(c).dot(&dir), align.dof.rx);
    result[4] = val_or_zero(align.dry(c).dot(&dir), align.dof.ry);
    result[5] = val_or_zero(align.drz(c).dot(&dir), align.dof.rz);

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
pub fn copy_jacobian<R, S>(j: &AlignStorage3, matrix: &mut Matrix<f64, R, U6, S>, row: usize)
where
    R: Dim,
    S: RawStorageMut<f64, R, U6> + Storage<f64, R, U6>,
{
    matrix.row_mut(row).copy_from_slice(j.as_slice());
}

#[cfg(test)]
mod tests {
    //! For lack of a more clever way of testing the jacobians, I've written a bunch of numerical
    //! versions of the jacobian calculations. The speed of the LM solver in the real implementation
    //! depends pretty heavily on how the real jacobian calculations work from making the assumption
    //! that an infinitesimally small motion of a test point against the closest surface of a
    //! reference won't change the location of the closest point, allowing the partial differentials
    //! to be approximated just from the deviation and surface directions.
    //!
    //! Since I had to figure these out on my own, I wanted to check them against something that
    //! I had relatively high confidence in, and the numerical method seemed like the most
    //! straightforward, obviously correct option.  I use the numerical alternative to check against
    //! a few different categories of cases.
    use super::*;
    use crate::UnitVec3;
    use crate::common::random_geometry::RandomGeometry3;
    use crate::geom3::{Point3, SurfacePoint3, Vector3};
    use approx::assert_relative_eq;
    use parry3d_f64::na::{Dyn, Owned};

    const NUMERIC_EPSILON: f64 = 1e-8;

    // ============================================================================================
    // The reverse point-to-surface jacobian, checked against finite differences
    // ============================================================================================

    /// The residual `point_surf_jacobian_rev` differentiates: the signed point-to-point distance
    /// from a fixed test point to a target point which the target's own transform moves.
    ///
    /// `c_local` and `n_local` describe the match in the target's own coordinates, so that
    /// perturbing the target's parameters moves it exactly as it would during a real solve.
    fn surf_rev_residual(
        params: &AlignParams3,
        p: &Point3,
        c_local: &Point3,
        n_local: &UnitVec3,
    ) -> f64 {
        let t = params.compute_transform();
        let c = t * c_local;
        let n = t.rotation * *n_local;
        let sp = SurfacePoint3::new(c, n);
        crate::common::dist(p, &c) * sp.scalar_projection(p).signum()
    }

    /// A finite-difference estimate of the partial derivative of that residual with respect to
    /// one of the target's six parameters.
    fn surf_rev_numeric(
        params: &AlignParams3,
        p: &Point3,
        c_local: &Point3,
        n_local: &UnitVec3,
        index: usize,
    ) -> f64 {
        let mut lo = params.clone();
        lo.set_index(index, params.storage()[index] - NUMERIC_EPSILON);
        let mut hi = params.clone();
        hi.set_index(index, params.storage()[index] + NUMERIC_EPSILON);

        let r_lo = surf_rev_residual(&lo, p, c_local, n_local);
        let r_hi = surf_rev_residual(&hi, p, c_local, n_local);
        (r_hi - r_lo) / (2.0 * NUMERIC_EPSILON)
    }

    #[test]
    fn stress_surf_rev_against_numeric() {
        // The analytic reverse jacobian against a central finite difference, over random local
        // origins, working offsets, parameter vectors, and geometry. This is the gate on the
        // reverse form: the sign convention and the requirement that the rotation partials be
        // evaluated at `c` are both easy to get wrong by inspection and obvious here.
        let mut rg = RandomGeometry3::from_seed(0x5eed_a11c);

        for _ in 0..2000 {
            let params =
                AlignParams3::new(AlignOrigin3::Local(rg.iso3(5.0)), Some(rg.iso3(5.0)), None)
                    .with_storage(AlignStorage3::from_iterator((0..6).map(|i| {
                        if i < 3 {
                            rg.f64_sym(3.0)
                        } else {
                            rg.f64_sym(0.8)
                        }
                    })));

            let c_local = rg.point(5.0);
            let n_local = rg.unit_vec();

            // Keep the test point well away from the match, so the deviation direction is
            // well conditioned and the residual is differentiable.
            let t = params.compute_transform();
            let c = t * c_local;
            let n = t.rotation * n_local;
            let p = c + n.into_inner() * rg.f64_sym(4.0).abs().max(0.5);

            let analytic =
                point_surf_jacobian_rev(&p, &SurfacePoint3::new(c, n), &params.compute_values());

            for i in 0..6 {
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
        let params = AlignParams3::from_origin(None);
        let values = params.compute_values();

        let c = Point3::new(1.0, 2.0, 3.0);
        let n = UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0));
        let p = Point3::new(1.0, 2.0, 4.5);
        let sp = SurfacePoint3::new(c, n);

        let fwd = point_surf_jacobian(&p, &sp, &values);
        let rev = point_surf_jacobian_rev(&p, &sp, &values);

        for i in 0..3 {
            assert_relative_eq!(fwd[i], -rev[i], epsilon = 1e-12);
        }
    }

    #[test]
    fn surf_rev_respects_locked_dof() {
        let dof = Dof6::new(true, false, true, false, true, false);
        let params = AlignParams3::from_origin(Some(dof));
        let values = params.compute_values();

        let c = Point3::new(1.0, 2.0, 3.0);
        let n = UnitVec3::new_normalize(Vector3::new(0.0, 1.0, 1.0));
        let p = Point3::new(2.0, 3.0, 5.0);

        let rev = point_surf_jacobian_rev(&p, &SurfacePoint3::new(c, n), &values);

        for i in [1usize, 3, 5] {
            assert_eq!(
                rev[i], 0.0,
                "locked parameter {i} should have a zero column"
            );
        }
    }

    #[test]
    fn test_jacobian_copy() {
        let x = AlignStorage3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        let mut target = Matrix::<f64, Dyn, U6, Owned<f64, Dyn, U6>>::zeros(10);
        // Copy the jacobian into the target matrix
        copy_jacobian(&x, &mut target, 4);

        assert_eq!(target[(4, 0)], 1.0);
        assert_eq!(target[(4, 1)], 2.0);
        assert_eq!(target[(4, 2)], 3.0);
        assert_eq!(target[(4, 3)], 4.0);
        assert_eq!(target[(4, 4)], 5.0);
        assert_eq!(target[(4, 5)], 6.0);
    }
}
