//! This module contains an abstraction for handling the parameters of 2d alignments re-expressed
//! as a transformation around a rotation center point.

use crate::geom2::align2::{iso2_from_param, param_from_iso2, T2Storage};
use crate::geom2::{Iso2, Point2};

/// Manages the parameters for a 2D alignment problem with a rotation center point.
///
/// # Overview
/// This struct is used to re-express the transformation problem as a translation in the world
/// coordinate system but a rotation around a different origin.
///
/// # The Forward Transformation
/// The rotation center point is
/// initially specified in the world coordinate system.  The
///
#[derive(Clone)]
pub struct RcParams2 {
    /// The rotation center point in the same coordinate system as the test entity
    rc: Point2,

    /// The current parameters
    x: T2Storage,

    /// The current active transformation
    transform: Iso2,

    /// The current inverse transformation
    inverse: Iso2,

    /// The current rotation-only matrix of the transformation
    rotation: Iso2,

    /// The currently active center of rotation, computed by transforming the rotation center point
    /// `rc` by the current transformation `transform`
    current_rc: Point2,
}

impl RcParams2 {
    ///
    ///
    /// # Arguments
    ///
    /// * `initial`:
    /// * `rc`:
    ///
    /// returns: RcParams2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn from_initial(initial: &Iso2, rc: &Point2) -> Self {
        // Get the initial transformation as a transformation about the rotation center point
        let about_rc = as_iso_about_center(rc, initial);
        let x = param_from_iso2(&about_rc);
        let current_rc = initial * rc;

        let mut item = Self {
            rc: *rc,
            x,
            transform: Iso2::identity(),
            inverse: Iso2::identity(),
            rotation: Iso2::rotation(initial.rotation.angle()),
            current_rc,
        };

        item.compute();
        item
    }

    /// Returns the unmoved rotation center point in the original global coordinate system.
    pub fn rc(&self) -> &Point2 {
        &self.rc
    }

    /// Returns the current rotation center point as it's been transformed by the currently active
    /// transformation.
    pub fn current_rc(&self) -> &Point2 {
        &self.current_rc
    }

    ///
    ///
    /// # Arguments
    ///
    /// * `x`:
    ///
    /// returns: ()
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn set(&mut self, x: &T2Storage) {
        self.x = *x;
        self.compute();
    }

    pub fn x(&self) -> &T2Storage {
        &self.x
    }

    pub fn transform(&self) -> &Iso2 {
        &self.transform
    }

    pub fn inverse(&self) -> &Iso2 {
        &self.inverse
    }

    pub fn rotation(&self) -> &Iso2 {
        &self.rotation
    }

    fn compute(&mut self) {
        let t = iso2_from_param(&self.x);
        self.rotation = Iso2::rotation(t.rotation.angle());
        self.transform = as_iso_about_origin(&self.rc, &t);
        self.inverse = self.transform.inverse();
        self.current_rc = self.transform * self.rc;
    }
}

/// Given an isometry in the global coordinate system, returns the isometry which is equivalent but
/// with its rotations about the given rotation center point.
///
/// # Arguments
///
/// * `rc`: The rotation center point
/// * `t`: The isometry in the global coordinate system
///
/// returns: Isometry<f64, Unit<Complex<f64>>, 2>
///
/// # Examples
///
/// ```
///
/// ```
fn as_iso_about_center(rc: &Point2, t: &Iso2) -> Iso2 {
    let fwd = Iso2::translation(rc.x, rc.y);
    let back = Iso2::translation(-rc.x, -rc.y);

    back * t * fwd
}

/// Given an isometry about the given rotation center point, returns the isometry in the global
/// coordinate system which is equivalent.
///
/// # Arguments
///
/// * `rc`: The rotation center point
/// * `t`: The isometry about the rotation center point
///
/// returns: Isometry<f64, Unit<Complex<f64>>, 2>
///
/// # Examples
///
/// ```
///
/// ```
fn as_iso_about_origin(rc: &Point2, t: &Iso2) -> Iso2 {
    let fwd = Iso2::translation(rc.x, rc.y);
    let back = Iso2::translation(-rc.x, -rc.y);

    fwd * t * back
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::Iso2;
    use crate::geom2::Point2;
    use approx::assert_relative_eq;
    use parry2d_f64::na::{Translation2, UnitComplex};
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn test_iso_about_origin_rot() {
        let rc = Point2::new(10.0, 20.0);
        let t = Iso2::rotation(FRAC_PI_2);

        let t2 = as_iso_about_origin(&rc, &t);
        let p0 = Point2::new(11.0, 20.0);
        let p1 = Point2::new(10.0, 21.0);

        let r0 = t2 * p0;
        let r1 = t2 * p1;
        assert_eq!(r0, Point2::new(10.0, 21.0));
        assert_eq!(r1, Point2::new(9.0, 20.0));
    }

    #[test]
    fn test_iso_about_origin_shift() {
        let rc = Point2::new(10.0, 20.0);
        let t = Iso2::translation(1.0, 1.0);

        let t2 = as_iso_about_origin(&rc, &t);
        let p0 = Point2::new(11.0, 20.0);
        let p1 = Point2::new(10.0, 21.0);

        let r0 = t2 * p0;
        let r1 = t2 * p1;
        assert_eq!(r0, Point2::new(12.0, 21.0));
        assert_eq!(r1, Point2::new(11.0, 22.0));
    }

    #[test]
    fn test_iso_about_point_rot() {
        // We'll calculate the isometry forward, assuming that test_iso_about_origin_XX is correct,
        // and then we'll test by using the as_iso_about_center function to reverse it
        let rc = Point2::new(10.0, 20.0);
        let t_rc = Iso2::rotation(FRAC_PI_2);

        let t_global = as_iso_about_origin(&rc, &t_rc);

        let t_rc2 = as_iso_about_center(&rc, &t_global);
        assert_relative_eq!(t_rc2, t_rc);
    }

    #[test]
    fn test_iso_about_point_shift() {
        // We'll calculate the isometry forward, assuming that test_iso_about_origin_XX is correct,
        // and then we'll test by using the as_iso_about_center function to reverse it
        let rc = Point2::new(10.0, 20.0);
        let t_rc = Iso2::translation(1.0, 1.0);

        let t_global = as_iso_about_origin(&rc, &t_rc);

        let t_rc2 = as_iso_about_center(&rc, &t_global);
        assert_relative_eq!(t_rc2, t_rc);
    }
}
