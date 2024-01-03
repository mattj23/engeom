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

/// Re-expresses an angle, specified in radians, in the range [0, 2pi].  If the angle was already
/// in the range [0, 2pi], it is returned unchanged.
///
/// # Arguments
///
/// * `angle`: The angle to re-express, in radians
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
pub fn angle_to_2pi(angle: f64) -> f64 {
    let mut angle = angle % (2.0 * PI);
    if angle < 0.0 {
        angle += 2.0 * PI;
    }
    angle
}

/// Returns the signed compliment of an angle, specified in radians, in the range [-2pi, 2pi].
///
/// # Arguments
///
/// * `angle`:
///
/// returns: f64
///
/// # Examples
///
/// ```
///
/// ```
pub fn signed_compliment_2pi(angle: f64) -> f64 {
    if angle >= 0.0 {
        (-2.0 * PI) + angle
    } else {
        2.0 * PI + angle
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
    fn test_angle_includes() {
        let mut rnd = thread_rng();
        for _ in 0..1000 {
            let start = rnd.gen_range(-2.0 * PI..2.0 * PI);
            let angle = rnd.gen_range(-2.0 * PI..2.0 * PI);
            let interval = AngleInterval::new(start, angle);

            for da in linear_space(0.0, angle, 100) {
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
                for da in to_check[1..to_check.len() - 2].iter() {
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
