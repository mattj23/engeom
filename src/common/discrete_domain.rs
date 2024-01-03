//! This module contains an abstraction for working with a discrete domain of scalar f64 values,
//! where the values are always ordered and only finite values are allowed.

use std::error::Error;
use crate::Result;
use crate::common::vecf64::{are_all_finite, are_in_ascending_order};

/// A discrete domain of scalar f64 values, in which all values are guaranteed to be finite and
/// in ascending order.
pub struct DiscreteDomain {
    values: Vec<f64>,
}

impl DiscreteDomain {
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn push(&mut self, value: f64) -> Result<()> {
        if !value.is_finite() {
            return Err(Box::from("Cannot add a non-finite value to a discrete domain"));
        }
        if !self.values.is_empty() && value < self.values[self.values.len() - 1] {
            return Err(Box::from("Cannot add a value to a discrete domain that is less than the last value"));
        }
        self.values.push(value);
        Ok(())
    }

}

impl Default for DiscreteDomain {
    fn default() -> Self {
        Self { values: Vec::new() }
    }
}

impl TryFrom<Vec<f64>> for DiscreteDomain {
    type Error = Box<dyn Error>;

    fn try_from(values: Vec<f64>) -> Result<Self> {
        if !are_all_finite(&values) {
            return Err(Box::from("Cannot create a discrete domain from a vector containing NaN or infinite values"));
        }

        if !are_in_ascending_order(&values) {
            return Err(Box::from("Cannot create a discrete domain from a vector that is not in ascending order"));
        }

        Ok(Self { values })
    }
}
