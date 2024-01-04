//! This module contains an abstraction for working with a discrete domain of scalar f64 values,
//! where the values are always ordered and only finite values are allowed.

use crate::common::vecf64::{are_all_finite, are_in_ascending_order};
use crate::{min_max, Result};
use std::error::Error;
use std::ops::Deref;

/// Generate a discrete domain of values which are linearly spaced between `start` and `end` and
/// which have a total count of `n`. The first value will be `start` and the last value will be
/// `end`.
///
/// # Arguments
///
/// * `start`: the starting value of the domain, inclusive
/// * `end`: the ending value of the domain, inclusive
/// * `n`: the total number of discrete, evenly spaced values in the domain
///
/// returns: DiscreteDomain
///
/// # Examples
///
/// ```
/// use engeom::common::linear_space;
/// let domain = linear_space(0.0, 1.0, 3);
/// assert_eq!(domain.values(), vec![0.0, 0.5, 1.0]);
/// ```
pub fn linear_space(start: f64, end: f64, n: usize) -> DiscreteDomain {
    let mut values = Vec::with_capacity(n);
    let step = (end - start) / (n - 1) as f64;
    for i in 0..n {
        values.push(start + i as f64 * step);
    }
    DiscreteDomain { values }
}

/// A discrete domain of scalar f64 values, in which all values are guaranteed to be finite and
/// in ascending order.
#[derive(Debug, Default)]
pub struct DiscreteDomain {
    values: Vec<f64>,
}

impl DiscreteDomain {
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    pub fn linear(start: f64, end: f64, n: usize) -> Self {
        let mut values = Vec::with_capacity(n);
        let (start, end) = min_max(start, end);
        let step = (end - start) / (n - 1) as f64;
        for i in 0..n {
            values.push(start + i as f64 * step);
        }
        DiscreteDomain { values }
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Try to push a value onto the end of the domain. The value must be finite and greater than
    /// the last value in the domain (unless the domain is empty).  If the value is not finite or
    /// is less than the last value in the domain, an error is returned.
    ///
    /// # Arguments
    ///
    /// * `value`: a finite value to add to the domain, must be greater than the last value in the
    /// domain (unless the domain is empty)
    ///
    /// returns: Result<(), Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::DiscreteDomain;
    /// let mut domain = DiscreteDomain::default();
    /// domain.push(1.0).unwrap();
    /// domain.push(2.0).unwrap();
    ///
    /// assert_eq!(domain.values(), vec![1.0, 2.0]);
    /// ```
    pub fn push(&mut self, value: f64) -> Result<()> {
        if !value.is_finite() {
            return Err(Box::from(
                "Cannot add a non-finite value to a discrete domain",
            ));
        }
        if !self.is_empty() && value < self.values[self.values.len() - 1] {
            return Err(Box::from(
                "Cannot add a value to a discrete domain that is less than the last value",
            ));
        }
        self.values.push(value);
        Ok(())
    }

    pub fn iter(&self) -> impl Iterator<Item = &f64> {
        self.values.iter()
    }
}

impl Deref for DiscreteDomain {
    type Target = [f64];

    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

impl TryFrom<Vec<f64>> for DiscreteDomain {
    type Error = Box<dyn Error>;

    fn try_from(values: Vec<f64>) -> Result<Self> {
        if !are_all_finite(&values) {
            return Err(Box::from(
                "Cannot create a discrete domain from a vector containing NaN or infinite values",
            ));
        }

        if !are_in_ascending_order(&values) {
            return Err(Box::from(
                "Cannot create a discrete domain from a vector that is not in ascending order",
            ));
        }

        Ok(Self { values })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iterate_values() {
        let mut working = Vec::new();
        for v in linear_space(0.0, 1.0, 3).iter() {
            working.push(*v);
        }

        assert_eq!(working, vec![0.0, 0.5, 1.0]);
    }

    #[test]
    fn try_linear_space() {
        let domain = linear_space(0.0, 1.0, 3);
        assert_eq!(domain.values(), vec![0.0, 0.5, 1.0]);
    }

    #[test]
    fn push_value() {
        let mut domain = DiscreteDomain::default();
        domain.push(1.0).unwrap();
        domain.push(2.0).unwrap();

        assert_eq!(domain.values(), vec![1.0, 2.0]);
    }

    #[test]
    fn try_from() {
        let domain = DiscreteDomain::try_from(vec![1.0, 2.0]).unwrap();
        assert_eq!(domain.values(), vec![1.0, 2.0]);
    }

    #[test]
    fn try_from_with_nan() {
        let result = DiscreteDomain::try_from(vec![1.0, f64::NAN]);
        assert!(result.is_err());
    }

    #[test]
    fn try_from_with_infinity() {
        let result = DiscreteDomain::try_from(vec![1.0, f64::INFINITY]);
        assert!(result.is_err());
    }

    #[test]
    fn try_from_with_descending_order() {
        let result = DiscreteDomain::try_from(vec![2.0, 1.0]);
        assert!(result.is_err());
    }
}
