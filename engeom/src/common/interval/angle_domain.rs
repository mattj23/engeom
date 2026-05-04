use crate::common::interval::IntervalOps;
use crate::common::{ANGLE_TOL, Interval, angle_to_2pi};
use std::f64::consts::PI;
use itertools::Itertools;

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
    /// [0, 2π].
    min: f64,

    /// The ending angle of the interval, in radians. If `max` is less than `min`, it's because the
    /// interval wraps beyond 2π.
    /// [0, 2π].
    max: f64,

    is_full: bool,
}

impl IntervalOps for AngleInterval {
    fn min(&self) -> f64 {
        self.min
    }

    fn max(&self) -> f64 {
        self.max
    }

    fn contains_value(&self, x: f64) -> bool {
        if self.is_full {
            return true;
        }

        let test = angle_to_2pi(x);

        if self.is_wrapping() {
            test >= self.min || test <= self.max
        } else {
            test <= self.max && test >= self.min
        }
    }

    fn contains_other(&self, other: Self) -> bool {
        if self.is_full {
            return true;
        }

        if other.is_full {
            return false;
        }

        // In order to fully contain the other, it must contain both the min and the max, however
        // that is not sufficient, because if one of the two wraps we might be spanning the ends
        // but disjoint in the middle
        if !self.contains_value(other.min) || !self.contains_value(other.max) {
            return false;
        }

        // To get the wrapping scenario, one or both of the intervals must wrap.
        match (self.is_wrapping(), other.is_wrapping()) {
            (false, false) => true, // neither wrap, so the contains above was sufficient
            (true, false) => {
                // We're wrapping and the other isn't, we exist between [min, 2π] and [0, max]. The
                // other doesn't wrap, so it must be entirely in one of those two ranges, or it
                // isn't contained
                let (a, b) = self.wrapping_sub_intervals().unwrap();
                a.contains_other(other) || b.contains_other(other)
            }
            (false, true) => {
                // We don't wrap, but the other does. This is an indication that we don't contain
                // the other, since by definition if it was contained inside this interval it also
                // wouldn't wrap
                false
            }
            (true, true) => {
                // We both wrap. We both have sub-intervals that exist between [min, 2π] and
                // [0, max]. In order for this interval to contain the other, our two intervals must
                // contain the other's two
                let (a0, b0) = self.wrapping_sub_intervals().unwrap();
                let (a1, b1) = other.wrapping_sub_intervals().unwrap();

                a0.contains_other(a1) && b0.contains_other(b1)
            }
        }
    }

    fn extent(&self) -> f64 {
        if self.is_full {
            return 2.0 * PI;
        }

        if self.is_wrapping() {
            2.0 * PI - (self.min - self.max)
        } else {
            self.max - self.min
        }
    }

    fn overlaps(&self, other: Self) -> bool {
        self.contains_value(other.min) || other.contains_value(self.max)
    }

    fn intersection(&self, other: Self) -> (Option<Interval>, Option<Interval>) {
        todo!()
    }

    fn clamp_value(&self, x: f64) -> f64 {
        todo!()
    }

    fn center(&self) -> f64 {
        todo!()
    }

    fn is_empty(&self) -> bool {
        todo!()
    }

    fn new_containing(&self, other: &Self) -> Self {
        todo!()
    }

    fn new_full() -> Self {
        todo!()
    }

    fn wraps() -> bool {
        todo!()
    }
}

impl AngleInterval {
    fn is_wrapping(&self) -> bool {
        self.max < self.min
    }

    fn wrapping_sub_intervals(&self) -> Option<(AngleInterval, AngleInterval)> {
        if self.is_wrapping() {
            let a0 = AngleInterval { min: self.min, max: 2.0 * PI, is_full: false};
            let a1 = AngleInterval { min: 0.0, max: self.max, is_full: false };
            Some((a0, a1))
        } else {
            None
        }
    }

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
        match angle {
            a if a.abs() >= 2.0 * PI => {
                // If the angle encircles more than the entire domain, we just return a full
                // interval that covers the entire domain
                Self {
                    min: 0.0,
                    max: 2.0 * PI,
                    is_full: true,
                }
            }
            a if a < 0.0 => {
                // The angle is not the full domain but is negative.
                let min = angle_to_2pi(start + a);
                let max = angle_to_2pi(start);
                Self {
                    min,
                    max,
                    is_full: false,
                }
            }
            _ => {
                let min = angle_to_2pi(start);
                let max = angle_to_2pi(start + angle);

                Self {
                    min,
                    max,
                    is_full: false,
                }
            }
        }
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
