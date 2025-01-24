//! This module contains common constructs for working with angles

use std::f64::consts::PI;

pub const ANGLE_TOL: f64 = 1.0e-12;

/// Enumerates the two possible directions of rotation, clockwise and counter-clockwise
#[derive(Copy, Clone, Debug)]
pub enum AngleDir {
    Cw,
    Ccw,
}

impl AngleDir {
    pub fn to_sign(self) -> f64 {
        match self {
            AngleDir::Cw => -1.0,
            AngleDir::Ccw => 1.0,
        }
    }

    pub fn from_sign(sign: f64) -> Self {
        if sign < 0.0 {
            AngleDir::Cw
        } else {
            AngleDir::Ccw
        }
    }

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
/// range [0, 2pi].
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
///
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

/// Re-expresses an angle, specified in radians, in the range [-pi, pi].  If the angle was already
/// in the range [-pi, pi], it is returned unchanged.
///
/// # Arguments
///
/// * `radians`:
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
    let mut angle = radians % (2.0 * PI);
    if angle > PI {
        angle -= 2.0 * PI;
    } else if angle < -PI {
        angle += 2.0 * PI;
    }
    angle
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

/// Returns the signed compliment of an angle, specified in radians, in the range [-2pi, 2pi].
///
/// # Arguments
///
/// * `radians`:
///
/// returns: f64
///
/// # Examples
///
/// ```
///
/// ```
pub fn signed_compliment_2pi(radians: f64) -> f64 {
    if radians >= 0.0 {
        (-2.0 * PI) + radians
    } else {
        2.0 * PI + radians
    }
}

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
    use rand::{thread_rng, Rng};
    use test_case::test_case;

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
        let mut rnd = thread_rng();
        for _ in 0..1000 {
            let angle = rnd.gen_range(-8.0 * PI..8.0 * PI);
            let test = angle_to_2pi(angle);
            assert!(
                test >= 0.0 && test < 2.0 * PI,
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
        let mut rnd = thread_rng();
        for _ in 0..1000 {
            let angle = rnd.gen_range(-8.0 * PI..8.0 * PI);
            let test = angle_signed_pi(angle);
            assert!(
                test >= -PI && test < PI,
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
        let mut rnd = thread_rng();
        let c = Circle2::new(0.0, 0.0, 1.0);

        for _ in 0..1000 {
            let start = rnd.gen_range(-2.0 * PI..2.0 * PI);
            let end = rnd.gen_range(-2.0 * PI..2.0 * PI);

            let v0 = c.point_at_angle(start);
            let v1 = c.point_at_angle(end);

            let test = angle_in_direction(start, end, AngleDir::Ccw);

            assert!(
                test >= 0.0 && test < 2.0 * PI,
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
        let mut rnd = thread_rng();
        let c = Circle2::new(0.0, 0.0, 1.0);

        for _ in 0..1000 {
            let start = rnd.gen_range(-2.0 * PI..2.0 * PI);
            let end = rnd.gen_range(-2.0 * PI..2.0 * PI);

            let v0 = c.point_at_angle(start);
            let v1 = c.point_at_angle(end);

            let test = angle_in_direction(start, end, AngleDir::Cw);

            assert!(
                test >= 0.0 && test < 2.0 * PI,
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
        let mut rnd = thread_rng();
        for _ in 0..1000 {
            let start = rnd.gen_range(-2.0 * PI..2.0 * PI);
            let angle = rnd.gen_range(-2.0 * PI..2.0 * PI);
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
