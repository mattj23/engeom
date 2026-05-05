use crate::common::interval::IntervalOps;
use crate::common::{Interval, angle_to_2pi};
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
        self.contains_value(other.min) || other.contains_value(self.min)
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

    fn offset(&self, x: f64) -> Self {
        Self::new_start_angle(self.min + x, self.extent())
    }

    fn new_full() -> Self {
        Self {
            min: 0.0,
            max: 2.0 * PI,
            is_full: true,
        }
    }

    fn wraps() -> bool {
        todo!()
    }
}

impl AngleInterval {
    /// Get whether the interval spans the wrapping line. If the interval does not wrap, then the
    /// min is lower than the max. If the interval does wrap, it means that the interval goes from
    /// min -> wrap-val, and from 0.0 -> max, and that max is less than min.
    fn is_wrapping(&self) -> bool {
        self.max < self.min
    }

    /// If the interval is wrapping (crosses the wrapping line and re-emerges at zero), this
    /// method will return the two sub intervals that together represent the full interval: the one
    /// that goes from `self.min` to 2π, and the one that goes from 0.0 to `self.max`.
    fn wrapping_sub_intervals(&self) -> Option<(AngleInterval, AngleInterval)> {
        if self.is_wrapping() {
            let a0 = AngleInterval {
                min: self.min,
                max: 2.0 * PI,
                is_full: false,
            };
            let a1 = AngleInterval {
                min: 0.0,
                max: self.max,
                is_full: false,
            };
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
    pub fn new_start_angle(start: f64, angle: f64) -> Self {
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
    use crate::common::{angle_to_2pi, linear_space, signed_compliment_2pi};
    use crate::{IntervalOps, Iso2, Vector2};
    use approx::assert_relative_eq;
    use rand::{RngExt, rng};
    use std::f64::consts::PI;
    use num_traits::Signed;

    fn v_at(theta: f64) -> Vector2 {
        Iso2::rotation(theta) * Vector2::x()
    }

    #[test]
    fn stress_creation() {
        let mut rng = rng();
        for _ in 0..10000 {
            let start = rng.random_range(-PI..PI);
            let span = rng.random_range(-(2.0 * PI)..(2.0 * PI));
            let interval = AngleInterval::new_start_angle(start, span);

            let v0 = v_at(start);
            let v1 = v_at(start + span);
            let (c_min, c_max) = if span.is_sign_positive() {
                (v0, v1)
            } else {
                (v1, v0)
            };

            assert_relative_eq!(v_at(interval.min), c_min, epsilon = 1e-6);
            assert_relative_eq!(v_at(interval.max), c_max, epsilon = 1e-6);
            assert_relative_eq!(interval.extent(), span.abs(), epsilon = 1.0e-8);
        }
    }

    #[test]
    fn offset_full_sweep() {
        let t0 = 0.0;
        let t1 = PI / 4.0;
        let original = AngleInterval::new_start_angle(t0, t1);
        assert_relative_eq!(original.extent(), t1, epsilon = 1.0e-8);
        assert_relative_eq!(original.min(), t0, epsilon = 1.0e-8);
        assert_relative_eq!(original.max(), t1, epsilon = 1.0e-8);

        let v0 = Vector2::x();
        let v1 = Iso2::rotation(t1) * v0;

        for t in linear_space(0.0, 2.0 * PI, 10000).iter() {
            let v0a = Iso2::rotation(*t) * v0;
            let v1a = Iso2::rotation(*t) * v1;

            let offset = original.offset(*t);

            assert_relative_eq!(original.extent(), offset.extent(), epsilon = 1.0e-8);
            assert_relative_eq!(v_at(offset.min()), v0a, epsilon = 1.0e-8);
            assert_relative_eq!(v_at(offset.max()), v1a, epsilon = 1.0e-8);
        }
    }

    fn reduce(val: f64) -> f64 {
        val - (1e-6 * val.signum())
    }

    #[test]
    fn stress_contains_value() {
        let mut rnd = rand::rng();
        for _ in 0..1000 {
            let start = rnd.random_range(-2.0 * PI..2.0 * PI);
            let angle = rnd.random_range(-2.0 * PI..2.0 * PI);
            let interval = AngleInterval::new_start_angle(start, angle);

            for da in linear_space(1e-6 * angle.signum(), reduce(angle), 100).values() {
                let test = start + da;
                assert!(
                    interval.contains_value(test),
                    "Failed Include {:?}, start={}, da={}, angle={}, test={} ({})",
                    interval,
                    start,
                    da,
                    angle,
                    test,
                    angle_to_2pi(test)
                );
            }

            let compliment = signed_compliment_2pi(angle);
            if compliment.abs() > 0.1 {
                let to_check = linear_space(1e-6 * compliment.signum(), reduce(compliment), 100);
                for da in to_check.values()[1..to_check.len() - 2].iter() {
                    let test = start + da;
                    assert!(
                        !interval.contains_value(test),
                        "Failed Exclude {:?}, start={}, da={}, angle={}, test={} ({})",
                        interval,
                        start,
                        da,
                        angle,
                        test,
                        angle_to_2pi(test)
                    );
                }
            }
        }
    }
}
