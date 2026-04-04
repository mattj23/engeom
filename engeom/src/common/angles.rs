//! This module contains common constructs for working with angles

use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

pub const ANGLE_TOL: f64 = 1.0e-12;

/// Lists the two possible directions of rotation, clockwise and counter-clockwise.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum AngleDir {
    /// Clockwise rotation (negative angular direction).
    Cw,
    /// Counter-clockwise rotation (positive angular direction).
    Ccw,
}

impl AngleDir {
    /// Returns `-1.0` for `Cw` and `1.0` for `Ccw`, suitable for use as a sign multiplier.
    pub fn to_sign(self) -> f64 {
        match self {
            AngleDir::Cw => -1.0,
            AngleDir::Ccw => 1.0,
        }
    }

    /// Returns `Cw` if `sign < 0.0`, otherwise `Ccw`.
    pub fn from_sign(sign: f64) -> Self {
        if sign < 0.0 {
            AngleDir::Cw
        } else {
            AngleDir::Ccw
        }
    }

    /// Returns the opposite direction: `Cw` <-> `Ccw`.
    pub fn opposite(self) -> Self {
        match self {
            AngleDir::Cw => AngleDir::Ccw,
            AngleDir::Ccw => AngleDir::Cw,
        }
    }
}

/// Calculates the angle between two angles in a given rotational direction. The angle returned
/// is the angle in the given rotational direction (clockwise or counter-clockwise) which `radians0`
/// would need to be rotated to align with `radians1`. The result will always be positive, in the
/// range [0, 2π].
///
/// # Arguments
///
/// * `radians0`: the starting angle, in radians
/// * `radians1`: the ending angle, in radians
/// * `angle_dir`: the rotational direction to consider
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::{angle_in_direction, AngleDir};
/// use std::f64::consts::{PI, FRAC_PI_2};
/// use approx::assert_relative_eq;
///
/// // Going CCW from 0 to π/2 is a quarter turn
/// let ccw = angle_in_direction(0.0, FRAC_PI_2, AngleDir::Ccw);
/// assert_relative_eq!(ccw, FRAC_PI_2, epsilon = 1.0e-10);
///
/// // Going CW from 0 to π/2 is three-quarter turns (the long way around)
/// let cw = angle_in_direction(0.0, FRAC_PI_2, AngleDir::Cw);
/// assert_relative_eq!(cw, 3.0 * FRAC_PI_2, epsilon = 1.0e-10);
/// ```
pub fn angle_in_direction(radians0: f64, radians1: f64, angle_dir: AngleDir) -> f64 {
    let t0 = angle_signed_pi(radians0);
    let t1 = angle_signed_pi(radians1);
    match angle_dir {
        AngleDir::Cw => {
            let t1 = if t1 > t0 { t1 - 2.0 * PI } else { t1 };
            t0 - t1
        }
        AngleDir::Ccw => {
            let t1 = if t1 < t0 { t1 + 2.0 * PI } else { t1 };
            t1 - t0
        }
    }
}

/// Returns the signed shortest angular distance from `radians0` to `radians1`.
///
/// A positive result means the shortest path is counter-clockwise; a negative result means the
/// shortest path is clockwise. The magnitude is always in the range [0, π].
///
/// # Arguments
///
/// * `radians0`: the starting angle, in radians
/// * `radians1`: the ending angle, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::shortest_angle_between;
/// use std::f64::consts::{PI, FRAC_PI_2, FRAC_PI_4};
/// use approx::assert_relative_eq;
///
/// // CCW from 0 to π/2 is the short way, so result is positive
/// let a = shortest_angle_between(0.0, FRAC_PI_2);
/// assert_relative_eq!(a, FRAC_PI_2, epsilon = 1.0e-10);
///
/// // CW from 0 to -π/4 is the short way, so result is negative
/// let b = shortest_angle_between(0.0, -FRAC_PI_4);
/// assert_relative_eq!(b, -FRAC_PI_4, epsilon = 1.0e-10);
/// ```
pub fn shortest_angle_between(radians0: f64, radians1: f64) -> f64 {
    let cw = angle_in_direction(radians0, radians1, AngleDir::Cw);
    let ccw = angle_in_direction(radians0, radians1, AngleDir::Ccw);
    if cw < ccw { -cw } else { ccw }
}

/// Re-expresses an angle, specified in radians, in the range [-π, π).  If the angle was already
/// in the range [-π, π), it is returned unchanged.
///
/// > **Note:** π (and any angle equivalent to π) maps to -π, not π, since the range is
/// > half-open. If your code depends on π being returned for inputs at exactly ±π, use a
/// > different normalization.
///
/// # Arguments
///
/// * `radians`: the angle to re-express, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::angle_signed_pi;
/// use std::f64::consts::PI;
/// use approx::assert_relative_eq;
/// let new_angle = angle_signed_pi(2.5 * PI);
/// assert_relative_eq!(new_angle, PI / 2.0, epsilon = 1.0e-10);
/// ```
pub fn angle_signed_pi(radians: f64) -> f64 {
    (radians + PI).rem_euclid(2.0 * PI) - PI
}

/// Re-expresses an angle, specified in radians, in the range [0, 2pi].  If the angle was already
/// in the range [0, 2pi], it is returned unchanged.
///
/// # Arguments
///
/// * `radians`: The angle to re-express, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::angle_to_2pi;
/// use std::f64::consts::PI;
/// use approx::assert_relative_eq;
/// let new_angle = angle_to_2pi(-PI);
/// assert_relative_eq!(new_angle, PI, epsilon = 1.0e-10);
/// ```
pub fn angle_to_2pi(radians: f64) -> f64 {
    let mut angle = radians % (2.0 * PI);
    if angle < 0.0 {
        angle += 2.0 * PI;
    }
    angle
}

/// Returns the signed complement of an angle with respect to a full rotation (2π).
///
/// Given an angle, this returns the remaining arc needed to complete a full circle, with the same
/// sign convention: a positive input returns a negative complement (since the remaining arc goes
/// the other way), and vice versa. The result is in the range (-2π, 2π].
///
/// This is useful when you have an arc angle and need the arc that covers the rest of the circle.
///
/// # Arguments
///
/// * `radians`: the angle to complement, in radians
///
/// returns: f64
///
/// # Examples
///
/// ```
/// use engeom::common::signed_compliment_2pi;
/// use std::f64::consts::{PI, FRAC_PI_2};
/// use approx::assert_relative_eq;
///
/// // A quarter turn CCW leaves three-quarter turns in the CW direction
/// let a = signed_compliment_2pi(FRAC_PI_2);
/// assert_relative_eq!(a, -3.0 * FRAC_PI_2, epsilon = 1.0e-10);
///
/// // A quarter turn CW leaves three-quarter turns in the CCW direction
/// let b = signed_compliment_2pi(-FRAC_PI_2);
/// assert_relative_eq!(b, 3.0 * FRAC_PI_2, epsilon = 1.0e-10);
/// ```
pub fn signed_compliment_2pi(radians: f64) -> f64 {
    if radians >= 0.0 {
        (-2.0 * PI) + radians
    } else {
        2.0 * PI + radians
    }
}

/// An `AngleInterval` represents a continuous range of angles, specified by a starting angle and
/// a positive (counter-clockwise) included length.  This is similar to an interval on a number
/// line, but with the added complexity that angles wrap.
///
/// When defining an `AngleInterval`, remember that all directions are positive (counter-clockwise).
/// To represent an interval with a negative length (for instance, starting at 0 and going to -π),
/// the interval must be defined as starting at -π and having a length of π. Some of the original
/// information is lost in this representation.
#[derive(Copy, Clone, Debug)]
pub struct AngleInterval {
    /// The starting angle of the interval, in radians. Will always take a value in the range
    /// [0, 2pi].
    start: f64,

    /// The inclusive angle of the interval, in radians. Will always take a value in the range
    /// [0, 2pi].
    angle: f64,
}

impl AngleInterval {
    /// Create a new `AngleInterval` with the given starting angle and included angle.  The
    /// included angle *may* be positive or negative, but if it is negative, the start and end
    /// will be reversed and the included angle will be inverted, and the directional information
    /// will be lost.
    ///
    /// In all cases, both the start angle and included angle will be converted and clamped into
    /// the range [0, 2pi].
    ///
    /// # Arguments
    ///
    /// * `start`: The starting angle of the interval, in radians
    /// * `angle`: The included angle of the interval, in radians. The end of the interval is
    ///   essentially `start + angle`. If `angle` is negative, the start and end will be swapped and
    ///   the included angle will be made positive.
    ///
    /// returns: AngleInterval
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new(start: f64, angle: f64) -> Self {
        if angle < 0.0 {
            let start = angle_to_2pi(start + angle);
            Self {
                start,
                angle: angle.abs().min(2.0 * PI),
            }
        } else {
            let start = angle_to_2pi(start);
            Self {
                start,
                angle: angle.min(2.0 * PI),
            }
        }
    }

    pub fn start(&self) -> f64 {
        self.start
    }

    pub fn angle(&self) -> f64 {
        self.angle
    }

    /// Returns true if the interval contains the given angle
    ///
    /// # Arguments
    ///
    /// * `angle`:
    ///
    /// returns: bool
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn contains(&self, angle: f64) -> bool {
        let angle = angle_to_2pi(angle);
        if angle >= self.start - ANGLE_TOL {
            angle <= self.start + self.angle + ANGLE_TOL
        } else {
            angle + 2.0 * PI <= self.start + self.angle + ANGLE_TOL
        }
    }

    /// Returns true if the interval intersects with the other interval.  An intersection occurs
    /// if either interval contains the start of the other interval.
    ///
    /// # Arguments
    ///
    /// * `other`:
    ///
    /// returns: bool
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn intersects(&self, other: &Self) -> bool {
        // In order for there to be an intersection, one of the intervals must contain the start
        // of the other interval.
        self.contains(other.start) || other.contains(self.start)
    }

    pub fn at_fraction(&self, f: f64) -> f64 {
        self.start + self.angle * f
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::linear_space;
    use crate::{Circle2, Iso2};
    use approx::assert_relative_eq;
    use rand::{Rng, RngExt};
    use test_case::test_case;

    #[test]
    fn angle_dir_to_sign() {
        assert_eq!(AngleDir::Cw.to_sign(), -1.0);
        assert_eq!(AngleDir::Ccw.to_sign(), 1.0);
    }

    #[test]
    fn angle_dir_from_sign() {
        assert_eq!(AngleDir::from_sign(-1.0), AngleDir::Cw);
        assert_eq!(AngleDir::from_sign(1.0), AngleDir::Ccw);
        assert_eq!(AngleDir::from_sign(0.0), AngleDir::Ccw);
    }

    #[test]
    fn angle_dir_opposite() {
        assert_eq!(AngleDir::Cw.opposite(), AngleDir::Ccw);
        assert_eq!(AngleDir::Ccw.opposite(), AngleDir::Cw);
    }

    #[test]
    fn angle_dir_roundtrip() {
        assert_eq!(AngleDir::from_sign(AngleDir::Cw.to_sign()), AngleDir::Cw);
        assert_eq!(AngleDir::from_sign(AngleDir::Ccw.to_sign()), AngleDir::Ccw);
    }

    #[test_case(90.0, -270.0)]
    #[test_case(180.0, -180.0)]
    #[test_case(270.0, -90.0)]
    #[test_case(-91.0, 269.0)]
    #[test_case(-181.0, 179.0)]
    #[test_case(-271.0, 89.0)]
    fn test_signed_compliment_0(angle: f64, compliment: f64) {
        let test = signed_compliment_2pi(angle.to_radians());
        assert_relative_eq!(test, compliment.to_radians(), epsilon = 1.0e-10);
    }

    #[test]
    fn stress_test_angle_to_2pi() {
        let mut rnd = rand::rng();
        for _ in 0..1000 {
            let angle = rnd.random_range(-8.0 * PI..8.0 * PI);
            let test = angle_to_2pi(angle);
            assert!(
                (0.0..2.0 * PI).contains(&test),
                "Failed Angle to 2pi: angle={}, test={}",
                angle,
                test
            );

            assert_relative_eq!(f64::sin(angle), f64::sin(test), epsilon = 1.0e-10);
            assert_relative_eq!(f64::cos(angle), f64::cos(test), epsilon = 1.0e-10);
        }
    }

    #[test]
    fn stress_test_angle_signed_pi() {
        let mut rnd = rand::rng();
        for _ in 0..1000 {
            let angle = rnd.random_range(-8.0 * PI..8.0 * PI);
            let test = angle_signed_pi(angle);
            assert!(
                (-PI..PI).contains(&test),
                "Failed Angle Signed Pi: angle={}, test={}",
                angle,
                test
            );

            assert_relative_eq!(f64::sin(angle), f64::sin(test), epsilon = 1.0e-10);
            assert_relative_eq!(f64::cos(angle), f64::cos(test), epsilon = 1.0e-10);
        }
    }

    #[test]
    fn stress_test_angle_in_direction_counterclockwise() {
        let mut rnd = rand::rng();
        let c = Circle2::new(0.0, 0.0, 1.0);

        for _ in 0..1000 {
            let start = rnd.random_range(-2.0 * PI..2.0 * PI);
            let end = rnd.random_range(-2.0 * PI..2.0 * PI);

            let v0 = c.point_at_angle(start);
            let v1 = c.point_at_angle(end);

            let test = angle_in_direction(start, end, AngleDir::Ccw);

            assert!(
                (0.0..2.0 * PI).contains(&test),
                "Failed Angle in Direction CW: start={}, end={}, test={} not in [0, 2pi]",
                start,
                end,
                test
            );

            // By rotating the start vector by the test angle, we should get the end vector
            let rot = Iso2::rotation(test);
            let test0 = rot * v0;

            assert_relative_eq!(test0.x, v1.x, epsilon = 1.0e-10);
            assert_relative_eq!(test0.y, v1.y, epsilon = 1.0e-10);
        }
    }
    #[test]
    fn stress_test_angle_in_direction_clockwise() {
        let mut rnd = rand::rng();
        let c = Circle2::new(0.0, 0.0, 1.0);

        for _ in 0..1000 {
            let start = rnd.random_range(-2.0 * PI..2.0 * PI);
            let end = rnd.random_range(-2.0 * PI..2.0 * PI);

            let v0 = c.point_at_angle(start);
            let v1 = c.point_at_angle(end);

            let test = angle_in_direction(start, end, AngleDir::Cw);

            assert!(
                (0.0..2.0 * PI).contains(&test),
                "Failed Angle in Direction CW: start={}, end={}, test={} not in [0, 2pi]",
                start,
                end,
                test
            );

            // By rotating the start vector by the test angle, we should get the end vector
            let rot = Iso2::rotation(-test);
            let test0 = rot * v0;

            assert_relative_eq!(test0.x, v1.x, epsilon = 1.0e-10);
            assert_relative_eq!(test0.y, v1.y, epsilon = 1.0e-10);
        }
    }

    #[test]
    fn test_angle_includes() {
        let mut rnd = rand::rng();
        for _ in 0..1000 {
            let start = rnd.random_range(-2.0 * PI..2.0 * PI);
            let angle = rnd.random_range(-2.0 * PI..2.0 * PI);
            let interval = AngleInterval::new(start, angle);

            for da in linear_space(0.0, angle, 100).values() {
                let test = start + da;
                assert!(
                    interval.contains(test),
                    "Failed Include {:?}, start={}, da={}, angle={}, test={}",
                    interval,
                    start,
                    da,
                    angle,
                    test
                );
            }

            let compliment = signed_compliment_2pi(angle);
            if compliment.abs() > 0.1 {
                let to_check = linear_space(0.0, compliment, 100);
                for da in to_check.values()[1..to_check.len() - 2].iter() {
                    let test = start + da;
                    assert!(
                        !interval.contains(test),
                        "Failed Exclude {:?}, start={}, da={}, angle={}, test={}",
                        interval,
                        start,
                        da,
                        angle,
                        test
                    );
                }
            }
        }
    }
}
