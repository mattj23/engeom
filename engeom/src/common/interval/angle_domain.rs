use crate::common::{ANGLE_TOL, angle_to_2pi};
use std::f64::consts::PI;

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
    use crate::common::interval::angle_domain::AngleInterval;
    use crate::common::{linear_space, signed_compliment_2pi};
    use rand::RngExt;
    use std::f64::consts::PI;

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
