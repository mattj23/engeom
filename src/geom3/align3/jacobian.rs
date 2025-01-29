//! This module contains common helpers for computing the values of the jacobian matrix for
//! different Levenberg-Marquardt alignments.

use super::*;
use crate::geom3::{Point3, SurfacePoint3};
use parry2d_f64::na::Dim;
use parry3d_f64::na::{Matrix, RawStorageMut, Storage, U6};

/// This is a helper function for computing the partial derivatives of the parameters (a single
/// row of the Jacobian matrix) for a distance function approximated by a point and its closest
/// point on a plane.  This is a reasonable approximation for distances measured between points and
/// the surface of a mesh, or points and a point/normal cloud where the points act locally like
/// planes.
///
/// This approximation assumes that the vector from the point to the closest point on the plane is
/// very close to (if not exactly) the normal of the plane.
///
/// # Arguments
///
/// * `p`: the test point (a sample point in the data being optimized)
/// * `c`: the reference surface point (a point and normal in the model) closest to `p`
/// * `rc`: a point which is the center of rotation for the test points
///
/// returns: Matrix<f64, Const<6>, Const<1>, ArrayStorage<f64, 6, 1>>
pub fn point_plane_jacobian(p: &Point3, c: &SurfacePoint3, params: &RcParams3) -> T3Storage {
    let s = c.scalar_projection(p).signum();

    // The point with relation to the current center of rotation
    let from_rc = Point3::from(p - params.current_rc());

    point_plane_core(s, c, from_rc, params)
}

/// This is a helper function for computing the partial derivatives of the parameters for a
/// distance function approximated by a point and its closest point on a plane (represented as
/// a `SurfacePoint3`), but where the parameters being evaluated are the transform of the reference
/// plane, not the test point.
///
/// This is very similar to `point_plane_jacobian`, but is used in multi-entity simultaneous
/// alignments where not only are the test points being transformed, but the reference points are
/// potentially being transformed by a transform on that entity.
///
/// # Arguments
///
/// * `p`: the test point (a sample point in the data being optimized)
/// * `c`: the reference surface point (a point and normal in the model) closest to `p`
/// * `rc`: a point which is the center of rotation **for the reference points**
///
/// returns: Matrix<f64, Const<6>, Const<1>, ArrayStorage<f64, 6, 1>>
pub fn point_plane_jacobian_rev(p: &Point3, c: &SurfacePoint3, params: &RcParams3) -> T3Storage {
    let s = c.scalar_projection(p).signum();

    // The point with relation to the current center of rotation
    let from_rc = Point3::from(c.point - params.current_rc());

    point_plane_core(-s, c, from_rc, params)
}

fn point_plane_core(s: f64, c: &SurfacePoint3, from_rc: Point3, params: &RcParams3) -> T3Storage {
    let mut result = T3Storage::zeros();
    let n = c.normal.into_inner() * s;

    result[0] = n.x;
    result[1] = n.y;
    result[2] = n.z;

    result[3] = n.dot(&(params.rotations().rd.x * from_rc).coords);
    result[4] = n.dot(&(params.rotations().rd.y * from_rc).coords);
    result[5] = n.dot(&(params.rotations().rd.z * from_rc).coords);

    result
}

pub fn point_point_jacobian(p: &Point3, c: &Point3, params: &RcParams3) -> T3Storage {
    let mut result = T3Storage::zeros();
    let m = p - c;
    if m.norm_squared() < 1e-16 {
        result
    } else {
        let n = m.normalize();
        result.x = n.x;
        result.y = n.y;
        result.z = n.z;

        let from_rc = Point3::from(p - params.current_rc());
        result.w = n.dot(&(params.rotations().rd.x * from_rc).coords);
        result.a = n.dot(&(params.rotations().rd.y * from_rc).coords);
        result.b = n.dot(&(params.rotations().rd.z * from_rc).coords);

        result
    }
}

/// Generic helper to copy the contents of a single row into a larger jacobian matrix of either
/// fixed or dynamic row count
///
/// # Arguments
///
/// * `j`:
/// * `matrix`:
/// * `row`:
///
/// returns: ()
pub fn copy_jacobian<R, S>(j: &T3Storage, matrix: &mut Matrix<f64, R, U6, S>, row: usize)
where
    R: Dim,
    S: RawStorageMut<f64, R, U6> + Storage<f64, R, U6>,
{
    matrix.row_mut(row).copy_from_slice(j.as_slice());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Point3;
    use approx::assert_relative_eq;
    use parry3d_f64::na::{Dyn, Owned};
    use std::f64::consts::PI;

    const NUMERIC_EPSILON: f64 = 1e-8;

    fn point_plane_numeric(
        params: &RcParams3,
        p: &Point3,
        closest: &SurfacePoint3,
        index: usize,
    ) -> f64 {
        let mut params = params.clone();
        let t_i = params.transform().inverse();
        let mut x = *params.x();
        x[index] += NUMERIC_EPSILON;
        params.set(&x);
        let t = params.transform() * t_i;

        let moved = t * *p;
        let d0 = closest.scalar_projection(&p).abs();
        let d1 = closest.scalar_projection(&moved).abs();
        (d1 - d0) / NUMERIC_EPSILON
    }

    fn point_plane_numeric_rev(
        params: &RcParams3,
        p: &Point3,
        closest: &SurfacePoint3,
        index: usize,
    ) -> f64 {
        let mut params = params.clone();
        let t_i = params.transform().inverse();
        let mut x = *params.x();
        x[index] += NUMERIC_EPSILON;
        params.set(&x);
        let t = params.transform() * t_i;

        // Get the original value
        let d0 = closest.scalar_projection(p);

        // Move the reference point
        let moved = closest.transformed(&t);
        let d1 = moved.scalar_projection(p);

        (d1 - d0) / NUMERIC_EPSILON
    }

    #[test]
    fn test_point_plane_translation() {
        let p = Point3::new(1.0, 2.0, 3.0);
        let c = Point3::new(0.0, 0.0, 0.0);
        let sp = SurfacePoint3::new_normalize(c, p - c);

        let initial = Iso3::from_parts(
            Translation3::new(8.0, -5.0, -6.0),
            UnitQuaternion::from_euler_angles(-0.2, 0.3, 0.5),
        );

        let rc = Point3::new(-1.0, -2.0, 3.0);
        let params = RcParams3::from_initial(&initial, &(initial.inverse() * rc));

        let test = point_plane_jacobian(&p, &sp, &params);

        assert_relative_eq!(
            point_plane_numeric(&params, &p, &sp, 0),
            test.x,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            point_plane_numeric(&params, &p, &sp, 1),
            test.y,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            point_plane_numeric(&params, &p, &sp, 2),
            test.z,
            epsilon = 1e-6
        );
    }

    /// This is the simplest possible test of the jacobians for the rotations, and involves a
    /// parameter set that is currently an identity transform and a center of rotation which is
    /// currently at the origin.
    #[test]
    fn test_point_plane_simple_rotations() {
        let p = Point3::new(1.0, 2.0, 3.0);
        let c = Point3::new(1.0, 3.0, 4.0);
        let sp = SurfacePoint3::new_normalize(c, p - c);
        let rc = Point3::origin();

        let params = RcParams3::from_initial(&Iso3::identity(), &rc);
        let test = point_plane_jacobian(&p, &sp, &params);

        let expected_w = point_plane_numeric(&params, &p, &sp, 3);
        let expected_a = point_plane_numeric(&params, &p, &sp, 4);
        let expected_b = point_plane_numeric(&params, &p, &sp, 5);
        assert_relative_eq!(expected_w, test.w, epsilon = 1e-6);
        assert_relative_eq!(expected_a, test.a, epsilon = 1e-6);
        assert_relative_eq!(expected_b, test.b, epsilon = 1e-6);
    }

    /// This test is a little more complicated, and involves a reference center of rotation which
    /// is not at the origin. However the parameter set is still an identity transform.
    #[test]
    fn test_point_plane_centered_rotations() {
        let p = Point3::new(1.0, 2.0, 3.0);
        let c = Point3::new(1.0, 3.0, 4.0);
        let sp = SurfacePoint3::new_normalize(c, p - c);
        let rc = Point3::new(-1.2, -3.5, -0.75);

        let params = RcParams3::from_initial(&Iso3::identity(), &rc);
        let test = point_plane_jacobian(&p, &sp, &params);

        let expected_w = point_plane_numeric(&params, &p, &sp, 3);
        let expected_a = point_plane_numeric(&params, &p, &sp, 4);
        let expected_b = point_plane_numeric(&params, &p, &sp, 5);
        assert_relative_eq!(expected_w, test.w, epsilon = 1e-6);
        assert_relative_eq!(expected_a, test.a, epsilon = 1e-6);
        assert_relative_eq!(expected_b, test.b, epsilon = 1e-6);
    }

    #[test]
    fn test_point_plane_initial_rz_rotation() {
        let p1 = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(1.0, 0.0, 0.0);
        let sp = SurfacePoint3::new_normalize(c, Vector3::new(0.0, 1.0, 0.0));
        let rc1 = Point3::origin();

        let initial = Iso3::from_parts(
            Translation3::identity(),
            UnitQuaternion::from_euler_angles(0.0, 0.0, -PI / 2.0),
        );

        let rc0 = initial.inverse() * rc1;
        let p0 = initial.inverse() * p1;

        let params = RcParams3::from_initial(&initial, &rc0);
        assert_relative_eq!(*params.current_rc(), rc1, epsilon = 1e-6);
        assert_relative_eq!(params.transform() * p0, p1, epsilon = 1e-6);

        let test = point_plane_jacobian(&p1, &sp, &params);
        //
        let expected_w = point_plane_numeric(&params, &p1, &sp, 3);
        let expected_a = point_plane_numeric(&params, &p1, &sp, 4);
        let expected_b = point_plane_numeric(&params, &p1, &sp, 5);
        assert_relative_eq!(expected_w, test.w, epsilon = 1e-6);
        assert_relative_eq!(expected_a, test.a, epsilon = 1e-6);
        assert_relative_eq!(expected_b, test.b, epsilon = 1e-6);
    }

    #[test]
    fn test_point_plane_rotation() {
        // These are the constructs after being transformed by the initial transform
        let p1 = Point3::new(10.0, 20.0, 30.0);
        let r1 = Point3::new(5.0, 8.0, 18.0);
        let c = Point3::new(12.0, 18.0, 28.0);
        let sp = SurfacePoint3::new_normalize(c, p1 - c);

        let initial = Iso3::from_parts(
            Translation3::new(8.0, -5.0, -6.0),
            UnitQuaternion::from_euler_angles(-0.2, 0.3, 0.5),
        );

        let p0 = initial.inverse() * p1;
        let r0 = initial.inverse() * r1;

        let params = RcParams3::from_initial(&initial, &r0);
        assert_relative_eq!(p1, params.transform() * p0, epsilon = 1e-8);
        assert_relative_eq!(r1, params.transform() * r0, epsilon = 1e-8);

        let test = point_plane_jacobian(&p1, &sp, &params);
        let expected_w = point_plane_numeric(&params, &p1, &sp, 3);
        let expected_a = point_plane_numeric(&params, &p1, &sp, 4);
        let expected_b = point_plane_numeric(&params, &p1, &sp, 5);
        assert_relative_eq!(expected_w, test.w, epsilon = 1e-6);
        assert_relative_eq!(expected_a, test.a, epsilon = 1e-6);
        assert_relative_eq!(expected_b, test.b, epsilon = 1e-6);
    }

    fn rev_test() -> (RcParams3, Point3, SurfacePoint3, T3Storage) {
        let p = Point3::new(2.0, 3.0, 4.0);
        let cp = Point3::new(1.0, 2.0, 3.0);
        let cn = Vector3::new(1.0, 1.0, 1.0);
        let c = SurfacePoint3::new_normalize(cp, cn);
        let rc = Point3::new(-0.5, -0.75, 3.25);

        let initial = Iso3::from_parts(
            Translation3::new(4.0, 3.0, 2.0),
            UnitQuaternion::from_euler_angles(0.2, -0.3, 0.5),
        );

        let rc0 = initial.inverse() * rc;
        let params = RcParams3::from_initial(&initial, &rc0);

        let test = point_plane_jacobian_rev(&p, &c, &params);
        (params, p, c, test)
    }

    #[test]
    fn test_point_plane_rev_translation() {
        let (params, p, c, test) = rev_test();
        assert_relative_eq!(
            point_plane_numeric_rev(&params, &p, &c, 0),
            test.x,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            point_plane_numeric_rev(&params, &p, &c, 1),
            test.y,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            point_plane_numeric_rev(&params, &p, &c, 2),
            test.z,
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_point_plane_rev_rotation() {
        let (params, p, c, test) = rev_test();
        assert_relative_eq!(
            point_plane_numeric_rev(&params, &p, &c, 3),
            test.w,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            point_plane_numeric_rev(&params, &p, &c, 4),
            test.a,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            point_plane_numeric_rev(&params, &p, &c, 5),
            test.b,
            epsilon = 1e-6
        );
    }

    #[test]
    fn test_jacobian_copy() {
        let x = T3Storage::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
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
