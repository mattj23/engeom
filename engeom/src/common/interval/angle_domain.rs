use crate::common::interval::IntervalOps;
use crate::common::{angle_ccw_to, angle_in_direction, angle_to_2pi, shortest_angle_between};
use std::f64::consts::PI;

const ANGLE_TOL: f64 = f64::EPSILON * 3.0 * PI;

use crate::AngleDir::Ccw;
#[cfg(test)]
use approx::{AbsDiffEq, RelativeEq};

/// An `AngleInterval` represents a continuous range of angles, specified by a starting angle and
/// a positive (counter-clockwise) included length.  This is similar to an interval on a number
/// line, but with the added complexity that angles wrap.
///
/// When defining an `AngleInterval`, remember that all directions are positive (counter-clockwise).
/// To represent an interval with a negative length (for instance, starting at 0 and going to -π),
/// the interval must be defined as starting at -π and having a length of π. Some of the original
/// information is lost in this representation.
#[derive(Copy, Clone, Debug, PartialEq)]
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
        angle_in_direction(self.min, x, Ccw) < self.extent()
    }

    fn contains_other(&self, other: Self) -> bool {
        if self.is_full {
            return true;
        }
        if self.is_empty() {
            return false;
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
        // This is used to determine `is_empty()`, so it cannot reference it

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

    fn intersection(&self, other: Self) -> (Option<Self>, Option<Self>) {
        if self.is_full {
            return (Some(other), None);
        }

        if other.is_full {
            return (Some(*self), None);
        }

        // If they're the same, return self. This is to distinguish from the case where they're
        // not identical, but compliments
        if (self.min - other.min).abs() < ANGLE_TOL && (self.max - other.max).abs() < ANGLE_TOL {
            return (Some(*self), None);
        }

        // First we check if one contains the start of the other
        let s_ovr = self.contains_value(other.min);
        let o_ovr = other.contains_value(self.min);

        match (s_ovr, o_ovr) {
            // They can only contain each other's starts if they were identical (which we guarded
            // for above), or if they overlap each other's ends via wrapping. In that case they
            // have two intersection regions
            (true, true) => {
                let a = ccw_between(self.min, other.max);
                let b = ccw_between(other.min, self.max);
                (Some(a), Some(b))
            }

            // We go from other.min to self.max
            (true, false) => (Some(ccw_between(other.min, self.max)), None),

            // We go from self.min to other.max
            (false, true) => (Some(ccw_between(self.min, other.max)), None),

            // No overlap
            (false, false) => (None, None),
        }
    }

    fn clamp_value(&self, x: f64) -> f64 {
        self.center() + shortest_angle_between(self.center(), x).signum() * self.extent() / 2.0
    }

    fn center(&self) -> f64 {
        angle_to_2pi(self.min + self.extent() / 2.0)
    }

    fn is_empty(&self) -> bool {
        if self.is_full {
            return false;
        }

        self.extent() < f64::EPSILON * 4.0 * PI
    }

    fn new_containing(&self, other: &Self) -> Self {
        if self.is_full || other.is_full {
            return Self::new_full();
        }

        // If they're the same, return self. This is to distinguish from the case where they're
        // not identical, but compliments
        if (self.min - other.min).abs() < ANGLE_TOL && (self.max - other.max).abs() < ANGLE_TOL {
            return *self;
        }

        // First we check if one contains the start of the other
        let s_ovr = self.contains_value(other.min);
        let o_ovr = other.contains_value(self.min);

        match (s_ovr, o_ovr) {
            // They can only contain each other's starts if they were identical (which we guarded
            // for above), or if they overlap each other's ends via wrapping. In that case they
            // span the whole domain.
            (true, true) => Self::new_full(),

            // We go from self.min to other.max
            (true, false) => Self::new_start_angle(self.min, angle_ccw_to(self.min, other.max)),

            // We go from other.min to self.max
            (false, true) => Self::new_start_angle(other.min, angle_ccw_to(other.min, self.max)),

            // No overlap, we need to determine which side is closest
            (false, false) => {
                let s_to_o = angle_ccw_to(self.max, other.min);
                let o_to_s = angle_ccw_to(other.max, self.min);
                match s_to_o < o_to_s {
                    // We go from self.min to other.max
                    true => Self::new_start_angle(self.min, angle_ccw_to(self.min, other.max)),
                    // We go from other.min to self.max
                    false => Self::new_start_angle(other.min, angle_ccw_to(other.min, self.max)),
                }
            }
        }
    }

    fn offset(&self, x: f64) -> Self {
        Self::new_start_angle(self.min + x, self.extent())
    }

    fn expanded(&self, half_width: f64) -> Self {
        let extent = self.extent() + 2.0 * half_width.abs();
        Self::new_start_angle(self.center() - extent / 2.0, extent)
    }

    fn dilated(&self, half_width: f64) -> Self {
        // An infinite interval can't be dilated
        if self.is_full {
            return Self::new_full();
        }

        let extent = (self.extent() - 2.0 * half_width.abs()).max(0.0);
        Self::new_start_angle(self.center() - extent / 2.0, extent)
    }

    fn new_full() -> Self {
        Self {
            min: 0.0,
            max: 2.0 * PI,
            is_full: true,
        }
    }

    fn wraps() -> bool {
        true
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

fn ccw_between(min: f64, max: f64) -> AngleInterval {
    AngleInterval::new_start_angle(min, angle_ccw_to(min, max))
}

#[cfg(test)]
impl AbsDiffEq for AngleInterval {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.min.abs_diff_eq(&other.min, epsilon)
            && self.max.abs_diff_eq(&other.max, epsilon)
            && self.is_full == other.is_full
    }
}

#[cfg(test)]
impl RelativeEq for AngleInterval {
    fn default_max_relative() -> Self::Epsilon {
        1e-8
    }

    fn relative_eq(
        &self,
        other: &Self,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        self.min.relative_eq(&other.min, epsilon, max_relative)
            && self.max.relative_eq(&other.max, epsilon, max_relative)
            && self.is_full == other.is_full
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::common::interval::angle_domain::AngleInterval;
    use crate::common::{angle_to_2pi, linear_space, signed_compliment_2pi};
    use crate::{IntervalOps, Iso2, Vector2};
    use approx::assert_relative_eq;
    use rand::{RngExt, rng};
    use std::f64::consts::PI;

    fn v_at(theta: f64) -> Vector2 {
        Iso2::rotation(theta) * Vector2::x()
    }

    fn by_deg(deg0: i32, deg1: i32) -> AngleInterval {
        AngleInterval {
            min: (deg0 as f64).to_radians(),
            max: (deg1 as f64).to_radians(),
            is_full: false,
        }
    }

    fn sweep_test(
        a: AngleInterval,
        b: AngleInterval,
        expected: AngleInterval,
        combine: &impl Fn(AngleInterval, AngleInterval) -> AngleInterval,
    ) {
        for t in linear_space(0.0, 2.0 * PI, 1000).iter() {
            let a_t = a.offset(*t);
            let b_t = b.offset(*t);
            let expected_t = expected.offset(*t);

            let test0 = combine(a_t, b_t);
            let test1 = combine(b_t, a_t);
            assert_relative_eq!(expected_t, test0, epsilon = 1e-10);
            assert_relative_eq!(expected_t, test1, epsilon = 1e-10);
        }
    }

    fn sweep_single(
        a: AngleInterval,
        e: AngleInterval,
        modify: &impl Fn(AngleInterval) -> AngleInterval,
    ) {
        for t in linear_space(0.0, 2.0 * PI, 1000).iter() {
            let a_t = a.offset(*t);
            let e_t = e.offset(*t);
            let test0 = modify(a_t);
            assert_relative_eq!(e_t, test0, epsilon = 1e-10);
        }
    }

    #[test]
    fn sweep_expand() {
        let a = by_deg(10, 20);
        let ex = by_deg(5, 25);
        sweep_single(a, ex, &|x| x.expanded(5.0f64.to_radians()));
    }

    #[test]
    fn sweep_dilate() {
        let ex = by_deg(10, 20);
        let a = by_deg(0, 30);
        sweep_single(a, ex, &|x| x.dilated(10.0f64.to_radians()));
    }

    #[test]
    fn sweep_intersection_wrap() {
        let a = by_deg(40, 320);
        let b = by_deg(310, 50);
        let e0 = by_deg(40, 50);
        let e1 = by_deg(310, 320);

        for t in linear_space(0.0, 2.0 * PI, 1000).iter() {
            let a_t = a.offset(*t);
            let b_t = b.offset(*t);
            let e0_t = e0.offset(*t);
            let e1_t = e1.offset(*t);

            let (t0_0, t0_1) = a_t.intersection(b_t);
            assert_relative_eq!(e0_t, t0_0.unwrap(), epsilon = 1e-10);
            assert_relative_eq!(e1_t, t0_1.unwrap(), epsilon = 1e-10);

            let (t1_0, t1_1) = b_t.intersection(a_t);
            assert_relative_eq!(e1_t, t1_0.unwrap(), epsilon = 1e-10);
            assert_relative_eq!(e0_t, t1_1.unwrap(), epsilon = 1e-10);
        }
    }

    #[test]
    fn sweep_intersection_overlap() {
        let a = by_deg(5, 50);
        let b = by_deg(30, 80);
        let expected = by_deg(30, 50);
        sweep_test(a, b, expected, &|x, y| {
            let (r0, r1) = x.intersection(y);
            assert!(
                r1.is_none(),
                "Expected only one intersection, got two: {:?} and {:?}",
                r0,
                r1
            );
            r0.unwrap()
        });
    }

    #[test]
    fn sweep_intersection_no_overlap() {
        let a = by_deg(5, 20);
        let b = by_deg(50, 80);
        for t in linear_space(0.0, 2.0 * PI, 1000).iter() {
            let a_t = a.offset(*t);
            let b_t = b.offset(*t);
            assert_eq!(a_t.intersection(b_t), (None, None));
            assert_eq!(b_t.intersection(a_t), (None, None));
        }
    }

    #[test]
    fn sweep_new_contains_overlap() {
        let a = by_deg(5, 50);
        let b = by_deg(30, 80);
        let expected = by_deg(5, 80);
        sweep_test(a, b, expected, &|x, y| x.new_containing(&y))
    }

    #[test]
    fn sweep_new_contains_no_overlap() {
        let a = by_deg(5, 20);
        let b = by_deg(30, 80);
        let expected = by_deg(5, 80);
        sweep_test(a, b, expected, &|x, y| x.new_containing(&y))
    }

    #[test]
    fn sweep_new_contains_equal() {
        let a = by_deg(5, 20);
        let b = by_deg(5, 20);
        let expected = by_deg(5, 20);
        sweep_test(a, b, expected, &|x, y| x.new_containing(&y))
    }

    #[test]
    fn sweep_new_contains_wrap() {
        let a = by_deg(20, 340);
        let b = by_deg(320, 40);
        let expected = AngleInterval::new_full();
        sweep_test(a, b, expected, &|x, y| x.new_containing(&y))
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
