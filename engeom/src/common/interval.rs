//! This module contains an abstraction which represents an interval on a continuous scalar domain,
//! such as the interval [0, 1] on the real number line.  Intervals allow for the testing of
//! intersections between ranges of values.

use crate::Result;
use serde::{Deserialize, Serialize};

/// An interval on a continuous scalar domain, such as the interval [0, 10] on the real number line.
/// Intervals can be thought of as 1d shapes, and are subject to boolean operations and tests
/// such as intersections, unions, containment, and composition.
#[derive(Debug, Copy, Clone, PartialEq, Serialize, Deserialize)]
pub struct Interval {
    /// The minimum value of the interval, inclusive.
    pub min: f64,

    /// The maximum value of the interval, inclusive.
    pub max: f64,
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
    pub fn try_new(min: f64, max: f64) -> Result<Self> {
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

    /// Returns the length of the interval.
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new(0.0, 1.0);
    /// assert_eq!(interval.length(), 1.0);
    /// ```
    pub fn length(&self) -> f64 {
        self.max - self.min
    }

    /// Returns true if the interval contains the given value.
    ///
    /// # Arguments
    ///
    /// * `x`: the value to test for containment
    ///
    /// returns: bool
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new(0.0, 1.0);
    /// assert!(interval.contains(0.5));
    /// assert!(!interval.contains(1.5));
    /// ```
    pub fn contains(&self, x: f64) -> bool {
        x >= self.min && x <= self.max
    }

    /// Returns true if the interval contains the other interval.
    ///
    /// # Arguments
    ///
    /// * `other`: the interval to test for containment
    ///
    /// returns: bool
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new(0.0, 1.0);
    /// let other = Interval::new(0.25, 0.75);
    /// assert!(interval.contains_interval(&other));
    /// ```
    pub fn contains_interval(&self, other: &Interval) -> bool {
        self.contains(other.min) && self.contains(other.max)
    }

    /// Returns true if the interval overlaps with the other interval.  An overlap occurs if either
    /// interval contains the start of the other interval.
    ///
    /// # Arguments
    ///
    /// * `other`: the interval to test for overlap
    ///
    /// returns: bool
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new(0.0, 1.0);
    /// let other = Interval::new(0.5, 1.5);
    /// assert!(interval.overlaps(&other));
    /// ```
    pub fn overlaps(&self, other: &Interval) -> bool {
        self.contains(other.min) || other.contains(self.min)
    }

    /// Returns the intersection of the interval with the other interval.  If the intervals do not
    /// overlap, the intersection will be None.
    ///
    /// # Arguments
    ///
    /// * `other`: the interval to intersect with
    ///
    /// returns: Option<Interval>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new(0.0, 1.0);
    /// let other = Interval::new(0.5, 1.5);
    /// if let Some(intersection) = interval.intersection(&other) {
    ///    assert_eq!(intersection, Interval::new(0.5, 1.0));
    /// } else {
    ///    panic!("interval.intersection returned None");
    /// }
    /// ```
    pub fn intersection(&self, other: &Interval) -> Option<Interval> {
        if self.overlaps(other) {
            Some(Interval::new(
                self.min.max(other.min),
                self.max.min(other.max),
            ))
        } else {
            None
        }
    }

    /// Clamps a value to the interval.
    ///
    /// # Arguments
    ///
    /// * `x`: the value to clamp
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::Interval;
    /// let interval = Interval::new(0.0, 1.0);
    ///
    /// // Clamping a value within the interval returns the value
    /// assert_eq!(interval.clamp(0.5), 0.5);
    ///
    /// // Clamping a value below the minimum returns the minimum
    /// assert_eq!(interval.clamp(-1.0), 0.0);
    ///
    /// // Clamping a value above the maximum returns the maximum
    /// assert_eq!(interval.clamp(2.0), 1.0);
    /// ```
    pub fn clamp(&self, x: f64) -> f64 {
        x.min(self.max).max(self.min)
    }
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
    fn interval_length() {
        let interval = Interval::new(0.0, 1.0);
        assert_eq!(interval.length(), 1.0);
    }

    #[test]
    fn interval_contains() {
        let interval = Interval::new(0.0, 1.0);
        assert!(interval.contains(0.5));
        assert!(!interval.contains(1.5));
    }

    #[test]
    fn interval_contains_interval() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.25, 0.75);
        assert!(interval.contains_interval(&other));
    }

    #[test]
    fn interval_doesnt_contain_interval() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.25, 1.25);
        assert!(!interval.contains_interval(&other));
    }

    #[test]
    fn interval_overlaps() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.5, 1.5);
        assert!(interval.overlaps(&other));
    }

    #[test]
    fn interval_doesnt_overlap() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(1.5, 2.5);
        assert!(!interval.overlaps(&other));
    }

    #[test]
    fn interval_intersection() {
        let interval = Interval::new(0.0, 1.0);
        let other = Interval::new(0.5, 1.5);
        if let Some(intersection) = interval.intersection(&other) {
            assert_eq!(intersection, Interval::new(0.5, 1.0));
        } else {
            panic!("interval.intersection returned None");
        }
    }
}
