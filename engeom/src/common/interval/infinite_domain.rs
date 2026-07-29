use crate::common::interval::IntervalOps;
use serde::{Deserialize, Serialize};

/// An interval on an infinite scalar domain, such as the interval [0, 10] on the real number line.
/// Intervals can be thought of as 1d shapes, and are subject to boolean operations and tests
/// such as intersections, unions, containment, and composition.
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct Interval {
    /// The minimum value of the interval, inclusive.
    pub min: f64,

    /// The maximum value of the interval, inclusive.
    pub max: f64,
}

impl IntervalOps for Interval {
    fn min(&self) -> f64 {
        self.min
    }

    fn max(&self) -> f64 {
        self.max
    }

    fn contains_value(&self, x: f64) -> bool {
        x >= self.min && x < self.max
    }

    fn contains_other(&self, other: Self) -> bool {
        self.contains_value(other.min) && self.contains_value(other.max)
    }

    fn extent(&self) -> f64 {
        self.max - self.min
    }

    fn overlaps(&self, other: Self) -> bool {
        self.contains_value(other.min) || other.contains_value(self.min)
    }

    fn intersection(&self, other: Self) -> (Option<Self>, Option<Self>) {
        if self.overlaps(other) {
            let a = Interval::new(self.min.max(other.min), self.max.min(other.max));
            (Some(a), None)
        } else {
            (None, None)
        }
    }

    fn clamp_value(&self, x: f64) -> f64 {
        x.clamp(self.min, self.max)
    }

    fn center(&self) -> f64 {
        (self.min + self.max) / 2.0
    }

    fn is_empty(&self) -> bool {
        // Intentionally `!(min < max)` rather than `min >= max`: a NaN bound makes `min < max`
        // false, so this treats a NaN interval as empty. `min >= max` would report it non-empty.
        #[allow(clippy::neg_cmp_op_on_partial_ord)]
        {
            !(self.min < self.max)
        }
    }

    fn new_containing(&self, other: &Self) -> Self {
        let min = self.min.min(other.min);
        let max = self.max.max(other.max);
        Interval::new_unchecked(min, max)
    }

    fn offset(&self, x: f64) -> Self {
        Self::new_unchecked(self.min + x, self.max + x)
    }

    fn expanded(&self, half_width: f64) -> Self {
        Self::new_unchecked(self.min - half_width.abs(), self.max + half_width.abs())
    }

    fn dilated(&self, half_width: f64) -> Self {
        let extent = (self.extent() - 2.0 * half_width).max(0.0);
        Self::new_unchecked(self.center() - extent / 2.0, self.center() + extent / 2.0)
    }

    fn new_full() -> Self {
        Interval::new_unchecked(f64::NEG_INFINITY, f64::INFINITY)
    }

    fn wraps() -> bool {
        false
    }
}

impl Interval {
    /// Create a new interval with the given minimum and maximum values.  The minimum and maximum
    /// values will be swapped if the minimum is greater than the maximum.
    ///
    /// An interval may contain infinite values, but not NaN values.  If either the minimum or
    /// maximum value is NaN, the creation of the interval will trigger a panic at runtime.  To
    /// avoid this, ensure that the minimum and maximum values are finite or use the `try_new`
    /// method instead.
    ///
    /// If you absolutely know that the minimum and maximum values are finite, you can use the
    /// `new_unchecked` method instead, which will not panic at runtime.
    ///
    /// # Arguments
    ///
    /// * `min`: the minimum value of the interval, inclusive
    /// * `max`: the maximum value of the interval, inclusive
    ///
    /// returns: Interval
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new(0.0, 1.0);
    /// assert_eq!(interval.min, 0.0);
    /// assert_eq!(interval.max, 1.0);
    /// ```
    pub fn new(min: f64, max: f64) -> Self {
        assert!(!min.is_nan());
        assert!(!max.is_nan());
        Self {
            min: min.min(max),
            max: min.max(max),
        }
    }

    /// Create a new interval with the given minimum and maximum values.  The minimum and maximum
    /// values will be swapped if the minimum is greater than the maximum.
    ///
    /// An interval may contain infinite values, but not NaN values.  If either the minimum or
    /// maximum value is NaN, the creation of the interval will return an error result.
    ///
    /// # Arguments
    ///
    /// * `min`: the minimum value of the interval, inclusive
    /// * `max`: the maximum value of the interval, inclusive
    ///
    /// returns: Result<Interval>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// if let Ok(interval) = Interval::try_new(0.0, 1.0) {
    ///     assert_eq!(interval.min, 0.0);
    ///     assert_eq!(interval.max, 1.0);
    /// } else {
    ///     panic!("Interval::try_new returned an error");
    /// }
    /// ```
    pub fn try_new(min: f64, max: f64) -> crate::Result<Self> {
        if min.is_nan() || max.is_nan() {
            Err("Interval::try_new received a NaN value".into())
        } else {
            Ok(Self {
                min: min.min(max),
                max: min.max(max),
            })
        }
    }

    /// Create a new interval with the given minimum and maximum values and perform no validation.
    /// The minimum and maximum values will *not* be swapped like `new` or `try_new` if the minimum
    /// is greater than the maximum, nor will either value be checked for NaN.
    ///
    /// If either the minimum or maximum value is NaN, or if the minimum is greater than the
    /// maximum, the behavior of the resulting `Interval` is undefined.
    ///
    /// # Arguments
    ///
    /// * `min`: the minimum value of the interval, inclusive
    /// * `max`: the maximum value of the interval, inclusive
    ///
    /// returns: Interval
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new_unchecked(0.0, 1.0);
    /// assert_eq!(interval.min, 0.0);
    /// assert_eq!(interval.max, 1.0);
    /// ```
    pub fn new_unchecked(min: f64, max: f64) -> Self {
        Self { min, max }
    }
    //
    // /// Returns the length of the interval.
    // ///
    // /// returns: f64
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use engeom::common::Interval;
    // /// let interval = Interval::new(0.0, 1.0);
    // /// assert_eq!(interval.length(), 1.0);
    // /// ```
    // pub fn length(&self) -> f64 {
    //     self.max - self.min
    // }
    //
    // /// Returns true if the interval contains the given value.
    // ///
    // /// # Arguments
    // ///
    // /// * `x`: the value to test for containment
    // ///
    // /// returns: bool
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use engeom::common::Interval;
    // /// let interval = Interval::new(0.0, 1.0);
    // /// assert!(interval.contains(0.5));
    // /// assert!(!interval.contains(1.5));
    // /// ```
    // pub fn contains(&self, x: f64) -> bool {
    //     x >= self.min && x <= self.max
    // }
    //
    // /// Returns true if the interval contains the other interval.
    // ///
    // /// # Arguments
    // ///
    // /// * `other`: the interval to test for containment
    // ///
    // /// returns: bool
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use engeom::common::Interval;
    // /// let interval = Interval::new(0.0, 1.0);
    // /// let other = Interval::new(0.25, 0.75);
    // /// assert!(interval.contains_interval(&other));
    // /// ```
    // pub fn contains_interval(&self, other: &Interval) -> bool {
    //     self.contains(other.min) && self.contains(other.max)
    // }
    //
    // /// Returns true if the interval overlaps with the other interval.  An overlap occurs if either
    // /// interval contains the start of the other interval.
    // ///
    // /// # Arguments
    // ///
    // /// * `other`: the interval to test for overlap
    // ///
    // /// returns: bool
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use engeom::common::Interval;
    // /// let interval = Interval::new(0.0, 1.0);
    // /// let other = Interval::new(0.5, 1.5);
    // /// assert!(interval.overlaps(&other));
    // /// ```
    // pub fn overlaps(&self, other: &Interval) -> bool {
    //     self.contains(other.min) || other.contains(self.min)
    // }
    //
    // /// Returns the intersection of the interval with the other interval.  If the intervals do not
    // /// overlap, the intersection will be None.
    // ///
    // /// # Arguments
    // ///
    // /// * `other`: the interval to intersect with
    // ///
    // /// returns: Option<Interval>
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use engeom::common::Interval;
    // /// let interval = Interval::new(0.0, 1.0);
    // /// let other = Interval::new(0.5, 1.5);
    // /// if let Some(intersection) = interval.intersection(&other) {
    // ///    assert_eq!(intersection, Interval::new(0.5, 1.0));
    // /// } else {
    // ///    panic!("interval.intersection returned None");
    // /// }
    // /// ```
    // pub fn intersection(&self, other: &Interval) -> Option<Interval> {
    //     if self.overlaps(other) {
    //         Some(Interval::new(
    //             self.min.max(other.min),
    //             self.max.min(other.max),
    //         ))
    //     } else {
    //         None
    //     }
    // }
    //
    // /// Clamps a value to the interval.
    // ///
    // /// # Arguments
    // ///
    // /// * `x`: the value to clamp
    // ///
    // /// returns: f64
    // ///
    // /// # Examples
    // ///
    // /// ```
    // /// use engeom::common::Interval;
    // /// let interval = Interval::new(0.0, 1.0);
    // ///
    // /// // Clamping a value within the interval returns the value
    // /// assert_eq!(interval.clamp(0.5), 0.5);
    // ///
    // /// // Clamping a value below the minimum returns the minimum
    // /// assert_eq!(interval.clamp(-1.0), 0.0);
    // ///
    // /// // Clamping a value above the maximum returns the maximum
    // /// assert_eq!(interval.clamp(2.0), 1.0);
    // /// ```
    // pub fn clamp(&self, x: f64) -> f64 {
    //     x.min(self.max).max(self.min)
    // }
    //
    // /// Returns the center of the interval, which is the average of the minimum and maximum values.
    // pub fn center(&self) -> f64 {
    //     (self.min + self.max) / 2.0
    // }
    //
    // pub fn width(&self) -> f64 {
    //     self.max - self.min
    // }
    //
    // pub fn is_empty(&self) -> bool {
    //     !(self.min < self.max)
    // }
    //
    // /// Creates an interval that fully contains both parent intervals.
    // ///
    // /// # Arguments
    // ///
    // /// * `a`: the first interval to join
    // /// * `b`: the second interval to join
    // ///
    // /// returns: Interval
    // pub fn new_contains(a: &Interval, b: &Interval) -> Interval {
    //     let min = a.min.min(b.min);
    //     let max = a.max.max(b.max);
    //     Interval::new_unchecked(min, max)
    // }
    //
    // /// Return a new interval which has been expanded by the half-width in each direction.
    // ///
    // /// # Arguments
    // ///
    // /// * `half_width`: The amount of distance to expand each side of the interval.
    // ///
    // /// returns: Interval
    // pub fn new_expanded(&self, half_width: f64) -> Interval {
    //     Interval::new_unchecked(self.min - half_width.abs(), self.max + half_width.abs())
    // }
    //
    // /// The opposite of `new_expanded`, this shrinks an interval by the half-width in each
    // /// direction. The total size, however, will not fall below zero. If the interval is shrunk by
    // /// more than its size, it will be left as a zero-width interval at the center.
    // ///
    // /// # Arguments
    // ///
    // /// * `half_width`:
    // ///
    // /// returns: Interval
    // ///
    // /// # Examples
    // ///
    // /// ```
    // ///
    // /// ```
    // pub fn new_dilated(&self, half_width: f64) -> Interval {
    //     let min = self.min + half_width.abs();
    //     let max = self.max - half_width.abs();
    //     if min > max {
    //         Interval::new_unchecked(self.center(), self.center())
    //     } else {
    //         Interval::new_unchecked(min, max)
    //     }
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_interval() {
        let interval = Interval::new(0.0, 1.0);
        assert_eq!(interval.min, 0.0);
        assert_eq!(interval.max, 1.0);
    }

    #[test]
    fn new_interval_values_flipped() {
        let interval = Interval::new(1.0, 0.0);
        assert_eq!(interval.min, 0.0);
        assert_eq!(interval.max, 1.0);
    }

    #[test]
    fn try_new_interval() {
        if let Ok(interval) = Interval::try_new(0.0, 1.0) {
            assert_eq!(interval.min, 0.0);
            assert_eq!(interval.max, 1.0);
        } else {
            panic!("Interval::try_new returned an error");
        }
    }

    #[test]
    #[should_panic]
    fn new_nan() {
        Interval::new(f64::NAN, 1.0);
    }

    #[test]
    fn new_unchecked_nan() {
        Interval::new_unchecked(f64::NAN, 1.0);
    }

    #[test]
    fn new_unchecked_swapped() {
        // This case tests that the minimum and maximum values are not fixed in the unchecked version, which should
        // naively accept whatever the user provides.
        let interval = Interval::new_unchecked(1.0, 0.0);
        assert_eq!(interval.min, 1.0);
        assert_eq!(interval.max, 0.0);
    }

    #[test]
    fn new_unchecked() {
        let interval = Interval::new_unchecked(0.0, 1.0);
        assert_eq!(interval.min, 0.0);
        assert_eq!(interval.max, 1.0);
    }

    #[test]
    fn interval_extent() {
        let interval = Interval::new(0.0, 1.0);
        assert_eq!(interval.extent(), 1.0);
    }

    #[test]
    fn interval_contains() {
        let interval = Interval::new(0.0, 1.0);
        assert!(interval.contains_value(0.5));
        assert!(!interval.contains_value(1.5));
    }

    #[test]
    fn interval_contains_interval() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.25, 0.75);
        assert!(interval.contains_other(other));
    }

    #[test]
    fn interval_doesnt_contain_interval() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.25, 1.25);
        assert!(!interval.contains_other(other));
    }

    #[test]
    fn interval_overlaps() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.5, 1.5);
        assert!(interval.overlaps(other));
    }

    #[test]
    fn interval_doesnt_overlap() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(1.5, 2.5);
        assert!(!interval.overlaps(other));
    }

    #[test]
    fn interval_intersection() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.5, 1.5);
        if let Some(intersection) = interval.intersection(other).0 {
            assert_eq!(intersection, Interval::new(0.5, 1.0));
        } else {
            panic!("interval.intersection returned None");
        }
    }
}
