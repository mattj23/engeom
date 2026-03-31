//! This module provides [`Averager`], a lightweight accumulator for computing the mean or weighted
//! mean of a stream of numeric values.

use num_traits::{One, Zero};
use std::ops::{Add, Div};

/// A lightweight accumulator for computing the arithmetic mean or weighted mean of a stream of
/// values.
///
/// `Averager` stores only two values internally — a running sum and a running count — making it
/// suitable for use when you need to track several independent averages simultaneously without
/// allocating a collection for each one.
///
/// # Type parameters
///
/// `T` must support addition, division, and the `Zero`/`One` identity constants from
/// [`num_traits`]. Common choices are `f32`, `f64`, and integer types.
///
/// # Examples
///
/// Unweighted average:
///
/// ```
/// use engeom::common::Averager;
///
/// let mut a = Averager::<f64>::new();
/// a.add(1.0);
/// a.add(2.0);
/// a.add(3.0);
/// assert_eq!(a.average(), Some(2.0));
/// ```
///
/// Weighted average:
///
/// ```
/// use engeom::common::Averager;
///
/// let mut a = Averager::<f64>::new();
/// a.add_weight(0.0, 1.0);  // contributes weight 1
/// a.add_weight(10.0, 3.0); // contributes weight 3 → result = (0 + 30) / 4 = 7.5
/// assert_eq!(a.average(), Some(7.5));
/// ```
pub struct Averager<T>
where
    T: Add + Div + Zero + One + Copy,
{
    sum: T,
    count: T,
}

impl<T> Default for Averager<T>
where
    T: Add<Output = T> + Div<Output = T> + Zero + One + Copy,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Averager<T>
where
    T: Add<Output = T> + Div<Output = T> + Zero + One + Copy,
{
    /// Creates a new, empty `Averager` with a sum and count of zero.
    pub fn new() -> Self {
        Self {
            sum: T::zero(),
            count: T::zero(),
        }
    }

    /// Adds a value with an implicit weight of one.
    pub fn add(&mut self, value: T) {
        self.sum = self.sum + value;
        self.count = self.count + T::one();
    }

    /// Adds a value with an explicit `weight`, contributing `value * weight` to the sum and
    /// `weight` to the total count.
    ///
    /// Mixing calls to [`add`](Self::add) and `add_weight` is valid: `add(v)` is equivalent to
    /// `add_weight(v, T::one())`.
    pub fn add_weight(&mut self, value: T, weight: T) {
        self.sum = self.sum + (value * weight);
        self.count = self.count + weight;
    }

    /// Returns the mean of all values added so far, or `None` if no values have been added.
    pub fn average(&self) -> Option<T> {
        if self.count.is_zero() {
            return None;
        }
        Some(self.sum / self.count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn empty_returns_none() {
        let a = Averager::<f64>::new();
        assert!(a.average().is_none());
    }

    #[test]
    fn default_returns_none() {
        let a = Averager::<f64>::default();
        assert!(a.average().is_none());
    }

    #[test]
    fn single_value_returns_that_value() {
        let mut a = Averager::<f64>::new();
        a.add(7.0);
        assert_relative_eq!(a.average().unwrap(), 7.0);
    }

    #[test]
    fn two_values_returns_mean() {
        let mut a = Averager::<f64>::new();
        a.add(3.0);
        a.add(7.0);
        assert_relative_eq!(a.average().unwrap(), 5.0);
    }

    #[test]
    fn multiple_values_returns_mean() {
        let mut a = Averager::<f64>::new();
        for v in [1.0, 2.0, 3.0, 4.0, 5.0] {
            a.add(v);
        }
        assert_relative_eq!(a.average().unwrap(), 3.0);
    }

    #[test]
    fn weighted_single_value_returns_that_value() {
        let mut a = Averager::<f64>::new();
        a.add_weight(4.0, 10.0);
        assert_relative_eq!(a.average().unwrap(), 4.0);
    }

    #[test]
    fn weighted_average_two_values() {
        let mut a = Averager::<f64>::new();
        a.add_weight(0.0, 1.0);
        a.add_weight(10.0, 3.0); // weight 3x heavier → expected (0*1 + 10*3) / 4 = 7.5
        assert_relative_eq!(a.average().unwrap(), 7.5);
    }

    #[test]
    fn equal_weights_match_unweighted() {
        let mut weighted = Averager::<f64>::new();
        let mut unweighted = Averager::<f64>::new();
        for v in [2.0, 5.0, 8.0] {
            weighted.add_weight(v, 1.0);
            unweighted.add(v);
        }
        assert_relative_eq!(
            weighted.average().unwrap(),
            unweighted.average().unwrap(),
            epsilon = 1e-10
        );
    }

    #[test]
    fn integer_average() {
        let mut a = Averager::<i32>::new();
        a.add(4);
        a.add(6);
        assert_eq!(a.average().unwrap(), 5);
    }
}
