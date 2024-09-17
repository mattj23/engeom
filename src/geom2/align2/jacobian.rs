//! This module contains common implementations for computing the values of the Jacobian matrix
//! for different 2D alignment Levenberg-Marquardt problems.

use crate::geom2::align2::{RcParams2, T2Storage};
use crate::geom2::{Point2, SurfacePoint2, Vector2};
use parry2d_f64::na::{Dim, Matrix, RawStorageMut, Storage, U3};

/// This is a helper function for computing the partial derivatives of the parameters (a single row
/// of the Jacobian matrix) for a distance function approximated by a 2d test point and its closest
/// point on a 2d surface.  This is a reasonable approximation for distances measured between points
/// and the surface of a continuous 2d domain, such as a curve or the boundary of a shape. It will
/// be an exact solution for the distance between a point and a straight line.
///
/// This approximation assumes that the vector between the test point and the closest point on the
/// surface is very close to (or exactly) the normal of the surface at that point. That is, the
/// test point lies on or very close to the line defined by the surface point and its normals.
///
/// # Arguments
///
/// * `p`:
/// * `s`:
/// * `params`:
///
/// returns: Matrix<f64, Const<3>, Const<1>, ArrayStorage<f64, 3, 1>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn point_surface_jacobian(p: &Point2, s: &SurfacePoint2, params: &RcParams2) -> T2Storage {
    // The vector of motion of the test point under rotation is the vector from the rotation center
    // turned 90 degrees.
    let from_rc = p - params.current_rc();
    let v_rot = Vector2::new(-from_rc.y, from_rc.x);

    T2Storage::new(s.normal.x, s.normal.y, s.normal.dot(&v_rot))
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
pub fn copy_jacobian<R, S>(j: &T2Storage, matrix: &mut Matrix<f64, R, U3, S>, row: usize)
where
    R: Dim,
    S: RawStorageMut<f64, R, U3> + Storage<f64, R, U3>,
{
    matrix.row_mut(row).copy_from_slice(j.as_slice());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{Iso2, Vector2};
    use approx::assert_relative_eq;
    use parry2d_f64::na::{Translation2, UnitComplex};
    use rand::distributions::Uniform;
    use rand::prelude::Distribution;
    use rand::thread_rng;
    use std::f64::consts::{FRAC_PI_2, PI};

    const NUMERIC_EPS: f64 = 1e-8;

    /// Perform the numeric approximation of the partial derivative of a single parameter between
    /// a test point and a surface point.
    ///
    /// The test point is a point which is assumed to be already transformed by the current
    /// parameters.
    ///
    /// # Arguments
    ///
    /// * `params`:
    /// * `p`:
    /// * `s`:
    /// * `index`:
    ///
    /// returns: f64
    fn point_surf_numeric(params: &RcParams2, p: &Point2, s: &SurfacePoint2, index: usize) -> f64 {
        let mut params = params.clone();
        let t_i = params.transform().inverse();
        let mut x = *params.x();

        x[index] += NUMERIC_EPS;
        params.set(&x);
        let t = params.transform() * t_i;

        let moved = t * *p;
        let d0 = s.scalar_projection(&p);
        let d1 = s.scalar_projection(&moved);
        (d1 - d0) / NUMERIC_EPS
    }

    /// Generate a point/surface point pair for testing. The surface point will be an offset from
    /// the point by the given translation vector, and the surface point's normal will be aimed
    /// directly at the test point, and flipped if d is negative.
    ///
    /// # Arguments
    ///
    /// * `px`: The x coordinate of the test point
    /// * `py`: The y coordinate of the test point
    /// * `tx`: The x distance to offset the surface point from the test point
    /// * `ty`: The y distance to offset the surface point from the test point
    /// * `d`: Set to -1 if the surface normal should be flipped away from the test point
    ///
    /// returns: (OPoint<f64, Const<2>>, SurfacePoint<2>)
    fn make_surf_test_pair(px: f64, py: f64, tx: f64, ty: f64, d: f64) -> (Point2, SurfacePoint2) {
        let p = Point2::new(px, py);
        let sp = p + Vector2::new(tx, ty);
        let s = SurfacePoint2::new_normalize(sp, (p - sp) * d);
        (p, s)
    }

    #[test]
    fn test_point_surf_numeric_translation0() {
        let (p, s) = make_surf_test_pair(0.0, 0.0, 1.0, 0.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));

        let tx = point_surf_numeric(&params, &p, &s, 0);
        let ty = point_surf_numeric(&params, &p, &s, 1);

        assert_relative_eq!(tx, -1.0, epsilon = 1e-6);
        assert_relative_eq!(ty, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_numeric_translation1() {
        let (p, s) = make_surf_test_pair(0.0, 0.0, 1.0, 0.0, -1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));

        let tx = point_surf_numeric(&params, &p, &s, 0);
        let ty = point_surf_numeric(&params, &p, &s, 1);

        assert_relative_eq!(tx, 1.0, epsilon = 1e-6);
        assert_relative_eq!(ty, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_numeric_rot0() {
        // Point is at (1, 0), surface is at (1, 1), rotation center is at (0, 0), so the
        // rotation brings the test point closer to the surface by the radius from the rotation
        // center.
        let (p, s) = make_surf_test_pair(1.0, 0.0, 0.0, 1.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));
        let tr = point_surf_numeric(&params, &p, &s, 2);

        assert_relative_eq!(tr, -1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_numeric_rot1() {
        // Point is at (2, 0), surface is at (2, 1), rotation center is at (0, 0), so the
        // rotation brings the test point closer to the surface by the radius from the rotation
        // center.
        let (p, s) = make_surf_test_pair(2.0, 0.0, 0.0, 1.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));
        let tr = point_surf_numeric(&params, &p, &s, 2);

        assert_relative_eq!(tr, -2.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_numeric_rot2() {
        // Point is at (1, 0), surface is at (2, 0), rotation center is at (0, 0), so the
        // rotation moves the point parallel to the surface
        let (p, s) = make_surf_test_pair(1.0, 0.0, 1.0, 0.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));
        let tr = point_surf_numeric(&params, &p, &s, 2);

        assert_relative_eq!(tr, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_numeric_rot_c0() {
        // Test when the rotation center is not at the origin
        // Point is at (1, 0), surface is at (1, 1), rotation center is at (2, 0), so the
        // rotation brings the test point farther from the surface by the radius from the rotation
        // center.
        let (p, s) = make_surf_test_pair(1.0, 0.0, 0.0, 1.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(2.0, 0.0));
        let tr = point_surf_numeric(&params, &p, &s, 2);

        assert_relative_eq!(tr, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_numeric_complex() {
        // Assume the original point is at (1, 0) and the original rotation center was at (0, 0).
        // The current parameters shift it by (1, 1) and rotate 90, so the current point is at
        // (1, 2) and the current rotation center is at (1, 1).  The surface is now at (2, 2).
        let (p, s) = make_surf_test_pair(1.0, 2.0, 1.0, 0.0, 1.0);
        let initial = Iso2::from_parts(Translation2::new(1.0, 1.0), UnitComplex::new(FRAC_PI_2));

        // Just check we did the math right
        let original = initial.inverse() * p;
        assert_relative_eq!(original, Point2::new(1.0, 0.0));

        // Create the parameters and check that the current rotation center is where we expect it
        let params = RcParams2::from_initial(&initial, &Point2::new(0.0, 0.0));
        assert_relative_eq!(*params.current_rc(), Point2::new(1.0, 1.0));

        // Now check the numeric approximation
        let tx = point_surf_numeric(&params, &p, &s, 0);
        let ty = point_surf_numeric(&params, &p, &s, 1);
        let tr = point_surf_numeric(&params, &p, &s, 2);

        assert_relative_eq!(tx, -1.0, epsilon = 1e-6);
        assert_relative_eq!(ty, 0.0, epsilon = 1e-6);
        assert_relative_eq!(tr, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_0() {
        let (p, s) = make_surf_test_pair(0.0, 0.0, 1.0, 0.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));
        let j = point_surface_jacobian(&p, &s, &params);

        assert_relative_eq!(j.x, -1.0, epsilon = 1e-6);
        assert_relative_eq!(j.y, 0.0, epsilon = 1e-6);
        assert_relative_eq!(j.z, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_1() {
        let (p, s) = make_surf_test_pair(1.0, 0.0, 0.0, 1.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));
        let j = point_surface_jacobian(&p, &s, &params);

        assert_relative_eq!(j.x, 0.0, epsilon = 1e-6);
        assert_relative_eq!(j.y, -1.0, epsilon = 1e-6);
        assert_relative_eq!(j.z, -1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_2() {
        let (p, s) = make_surf_test_pair(1.0, 0.0, 0.0, 1.0, -1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(0.0, 0.0));
        let j = point_surface_jacobian(&p, &s, &params);

        assert_relative_eq!(j.x, 0.0, epsilon = 1e-6);
        assert_relative_eq!(j.y, 1.0, epsilon = 1e-6);
        assert_relative_eq!(j.z, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_rot_c0() {
        // Test when the rotation center is not at the origin
        // Point is at (1, 0), surface is at (1, 1), rotation center is at (2, 0), so the
        // rotation brings the test point farther from the surface by the radius from the rotation
        // center.
        let (p, s) = make_surf_test_pair(1.0, 0.0, 0.0, 1.0, 1.0);
        let params = RcParams2::from_initial(&Iso2::identity(), &Point2::new(2.0, 0.0));
        let j = point_surface_jacobian(&p, &s, &params);

        assert_relative_eq!(j.x, 0.0, epsilon = 1e-6);
        assert_relative_eq!(j.y, -1.0, epsilon = 1e-6);
        assert_relative_eq!(j.z, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_point_surf_complex() {
        // Assume the original point is at (1, 0) and the original rotation center was at (0, 0).
        // The current parameters shift it by (1, 1) and rotate 90, so the current point is at
        // (1, 2) and the current rotation center is at (1, 1).  The surface is now at (2, 2).
        let (p, s) = make_surf_test_pair(1.0, 2.0, 1.0, 0.0, 1.0);
        let initial = Iso2::from_parts(Translation2::new(1.0, 1.0), UnitComplex::new(FRAC_PI_2));

        // Just check we did the math right
        let original = initial.inverse() * p;
        assert_relative_eq!(original, Point2::new(1.0, 0.0));

        // Create the parameters and check that the current rotation center is where we expect it
        let params = RcParams2::from_initial(&initial, &Point2::new(0.0, 0.0));
        assert_relative_eq!(*params.current_rc(), Point2::new(1.0, 1.0));

        // Now check the numeric approximation
        let j = point_surface_jacobian(&p, &s, &params);

        assert_relative_eq!(j.x, -1.0, epsilon = 1e-6);
        assert_relative_eq!(j.y, 0.0, epsilon = 1e-6);
        assert_relative_eq!(j.z, 1.0, epsilon = 1e-6);
    }

    fn random_iso2() -> Iso2 {
        let mut rng = thread_rng();
        let v = Vector2::new(
            Uniform::from(-10.0..10.0).sample(&mut rng),
            Uniform::from(-10.0..10.0).sample(&mut rng),
        );
        let r = Uniform::from(-PI..PI).sample(&mut rng);
        Iso2::from_parts(Translation2::from(v), UnitComplex::new(r))
    }

    fn random_point() -> Point2 {
        let mut rng = thread_rng();
        Point2::new(
            Uniform::from(-10.0..10.0).sample(&mut rng),
            Uniform::from(-10.0..10.0).sample(&mut rng),
        )
    }

    fn random_dir() -> f64 {
        let mut rng = thread_rng();
        Uniform::from(-1.0..1.0_f64).sample(&mut rng).signum()
    }

    #[test]
    fn stress_point_surf_against_numeric() {
        for _ in 0..10000 {
            let p = random_point();
            let (p, s) = make_surf_test_pair(p.x, p.y, 1.0, 0.0, random_dir());
            let rc = random_point();
            let params = RcParams2::from_initial(&random_iso2(), &rc);

            let tx = point_surf_numeric(&params, &p, &s, 0);
            let ty = point_surf_numeric(&params, &p, &s, 1);
            let tr = point_surf_numeric(&params, &p, &s, 2);

            let j = point_surface_jacobian(&p, &s, &params);

            assert_relative_eq!(j.x, tx, epsilon = 1e-6);
            assert_relative_eq!(j.y, ty, epsilon = 1e-6);
            assert_relative_eq!(j.z, tr, epsilon = 1e-6);
        }
    }
}
