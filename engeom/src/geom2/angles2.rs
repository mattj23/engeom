//! Tools for working with angles specifically in 2D

use crate::common::AngleDir;
use crate::geom2::{Iso2, Vector2};
use std::f64::consts::{FRAC_PI_2, PI};

/// Convenience function for creating a rotation matrix for a rotation of 90 degrees in the
/// specified direction.  A counter-clockwise rotation corresponds with a positive angle, and a
/// clockwise rotation corresponds with a negative angle.
///
/// # Arguments
///
/// * `dir`: The direction of the rotation
///
/// returns: Isometry<f64, Unit<Complex<f64>>, 2>
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::AngleDir::{Ccw, Cw};
/// use engeom::geom2::{Point2, rot90};
///
/// let p_ccw = rot90(Ccw) * Point2::new(1.0, 0.0);
/// let p_cw = rot90(Cw) * Point2::new(1.0, 0.0);
///
/// assert_relative_eq!(p_ccw, Point2::new(0.0, 1.0));
/// assert_relative_eq!(p_cw, Point2::new(0.0, -1.0));
/// ```
pub fn rot90(dir: AngleDir) -> Iso2 {
    match dir {
        AngleDir::Ccw => Iso2::rotation(FRAC_PI_2),
        AngleDir::Cw => Iso2::rotation(-FRAC_PI_2),
    }
}

/// Convenience function for creating a rotation matrix for a rotation of 270 degrees (-90) in the
/// specified direction.  A counter-clockwise rotation corresponds with a positive angle, and a
/// clockwise rotation corresponds with a negative angle.
///
/// # Arguments
///
/// * `dir`:
///
/// returns: Isometry<f64, Unit<Complex<f64>>, 2>
pub fn rot270(dir: AngleDir) -> Iso2 {
    match dir {
        AngleDir::Ccw => Iso2::rotation(-FRAC_PI_2),
        AngleDir::Cw => Iso2::rotation(FRAC_PI_2),
    }
}

/// Compute the signed angle from vector `v1` to vector `v2`. The sign of the angle is positive if
/// the shortest rotation from `v1` to `v2` is counter-clockwise, and negative if it is clockwise.
/// The angle is measured in radians, and will be in the range [-pi, pi].
///
/// Alternately, you can think of the angle as the amount of rotation needed to rotate `v1` to
/// match the direction of `v2`.
///
/// # Arguments
///
/// * `v1`: The reference vector
/// * `v2`: The vector to which the angle is measured
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::Vector2;
/// use engeom::geom2::signed_angle;
///
/// let v = Vector2::new(1.0, 0.0);
///
/// let a0 = signed_angle(&v, &Vector2::new(1.0, 1.0));
/// assert_relative_eq!(a0.to_degrees(), 45.0);
///
/// let a1 = signed_angle(&v, &Vector2::new(1.0, -1.0));
/// assert_relative_eq!(a1.to_degrees(), -45.0);
///
/// let a2 = signed_angle(&v, &Vector2::new(-1.0, 1.0));
/// assert_relative_eq!(a2.to_degrees(), 135.0);
///
/// let a3 = signed_angle(&v, &Vector2::new(-1.0, -1.0));
/// assert_relative_eq!(a3.to_degrees(), -135.0);
/// ```
pub fn signed_angle(v1: &Vector2, v2: &Vector2) -> f64 {
    (v1.x * v2.y - v1.y * v2.x).atan2(v1.x * v2.x + v1.y * v2.y)
}

/// Returns the angle between two vectors, in radians, in the direction specified.  This can be
/// thought of as the amount of rotation needed to rotate `v1` to match the direction of `v2` when
/// rotating only in `direction`.  The angle will be in the range [0, pi].
///
/// # Arguments
///
/// * `v1`: The reference vector
/// * `v2`: The vector to which the angle is measured
/// * `direction`: The direction in which to measure the angle
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use approx::assert_relative_eq;
/// use engeom::geom2::{Vector2, directed_angle};
/// use engeom::AngleDir::{Ccw, Cw};
///
/// let v = Vector2::new(1.0, 0.0);
///
/// let v0 = Vector2::new(1.0, 1.0);
/// assert_relative_eq!(directed_angle(&v, &v0, Ccw).to_degrees(), 45.0);
/// assert_relative_eq!(directed_angle(&v, &v0, Cw).to_degrees(), 315.0);
/// ```
pub fn directed_angle(v1: &Vector2, v2: &Vector2, direction: AngleDir) -> f64 {
    use AngleDir::{Ccw, Cw};
    let a = signed_angle(v1, v2)
        * match direction {
            Ccw => 1.0,
            Cw => -1.0,
        };
    if a < 0.0 {
        a + 2.0 * PI
    } else {
        a
    }
}
