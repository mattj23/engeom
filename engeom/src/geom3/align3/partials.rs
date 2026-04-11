//! This module is for working with the partial derivatives of motion in the alignment problem.

#[cfg(test)]
mod tests {
    //! These tests are here as I worked out the partial derivatives of motion using numerical
    //! approximation to guide myself.  I started with very simple cases and gradually added
    //! complexity to make sure that my understanding was solid.
    //!
    //! The [`AlignParams3`] struct holds the parameters being optimized in a 3D alignment problem,
    //! expressed as an euler angle rotation + translation problem around an arbitrary rotation
    //! centerpoint. Additionally, a working transformation exists which can alter the neutral
    //! relationship between the test and target entities, and a general 6-dof constraint can be
    //! applied to suppress parameters.
    //!
    //! The

    use std::f64::consts::PI;
    use crate::geom3::align3::Dof6;
    use crate::geom3::align3::params::{AlignOrigin, AlignParams3};
    use crate::{Iso3, Point3, Vector3};
    use approx::assert_relative_eq;
    use crate::na::{Translation, Translation3, UnitQuaternion};

    const ANGLE_EPSILON: f64 = 1e-8;
    const TRANS_EPSILON: f64 = 1e-8;

    /// This function computes the vector of motion of a point in the original coordinate system
    /// of the test entity for a finitely approximated infinitesimal change in one of the
    /// parameters. The parameters by index are: (0: tx, 1: ty, 2: tz, 3: rx, 4: ry, 5: rz)
    fn finite_diff(params: &AlignParams3, point: &Point3, index: usize) -> Vector3 {
        let mut w0 = params.clone();
        let mut w1 = params.clone();
        let stored = params.get_storage();

        let eps = if index < 3 {
            TRANS_EPSILON
        } else {
            ANGLE_EPSILON
        };

        w0.set_index(index, stored[index] - eps);
        w1.set_index(index, stored[index] + eps);

        let p0 = w0.current_values().transform * point;
        let p1 = w1.current_values().transform * point;

        (p1 - p0) / (2.0 * eps)
    }

    #[test]
    fn partials_of_translations_at_zero() {
        let params = AlignParams3::new_local(AlignOrigin::Origin, None);
        let current = params.current_values();

        assert_relative_eq!(current.dtx, Vector3::x_axis(), epsilon = 1e-12);
        assert_relative_eq!(current.dty, Vector3::y_axis(), epsilon = 1e-12);
        assert_relative_eq!(current.dtz, Vector3::z_axis(), epsilon = 1e-12);
    }

    #[test]
    fn partials_of_translations_with_rotations() {
        let params = AlignParams3::new_local(AlignOrigin::Origin,None)
            .with_rx(0.1)
            .with_ry(0.2)
            .with_rz(0.3);
        let test_point = Point3::new(1.0, 2.0, 3.0);

        let exp_x = finite_diff(&params, &test_point, 0);
        let exp_y = finite_diff(&params, &test_point, 1);
        let exp_z = finite_diff(&params, &test_point, 2);
        let c = params.current_values();

        // Because translations are applied after rotations, the direction vectors don't change
        // when rotation parameters are applied
        assert_relative_eq!(c.dtx, Vector3::x_axis(), epsilon = 1e-12);
        assert_relative_eq!(c.dty, Vector3::y_axis(), epsilon = 1e-12);
        assert_relative_eq!(c.dtz, Vector3::z_axis(), epsilon = 1e-12);

        // Sanity check that I've thought through this correctly
        assert_relative_eq!(exp_x, c.dtx, epsilon = 1e-8);
        assert_relative_eq!(exp_y, c.dty, epsilon = 1e-8);
        assert_relative_eq!(exp_z, c.dtz, epsilon = 1e-8);
    }

    #[test]
    fn partials_of_translations_with_working_transform() {
        // The working transformation is applied before the parameter transform, and its role is to
        // effectively do the same thing as transforming the test entity by the working transform
        // before doing the alignment. This means that the directions of X, Y, and Z may be
        // different.

        // let working = Iso3::from_parts(
        //     Translation3::new(2.0, 3.0, 4.0),
        //     UnitQuaternion::from_euler_angles(PI / 3.0, PI / 4.0, PI / 5.0),
        // );
        //
        // let params = AlignParams3::new();
        // let test_point = Point3::new(1.0, 2.0, 3.0);
        //
        // let exp_x = finite_diff(&params, &test_point, 0);
        // let exp_y = finite_diff(&params, &test_point, 1);
        // let exp_z = finite_diff(&params, &test_point, 2);
        // let c = params.current_values();
        //
        // assert_relative_eq!(exp_x, c.dtx, epsilon = 1e-8);
        // assert_relative_eq!(exp_y, c.dty, epsilon = 1e-8);
        // assert_relative_eq!(exp_z, c.dtz, epsilon = 1e-8);
    }
}
