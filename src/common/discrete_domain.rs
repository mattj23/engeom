//! This module contains an abstraction for working with a discrete domain of scalar f64 values,
//! where the values are always ordered and only finite values are allowed.

pub struct DiscreteDomain {
    values: Vec<f64>,
}

impl Default for DiscreteDomain {
    fn default() -> Self {
        Self { values: Vec::new() }
    }
}

impl From<Vec<f64>> for DiscreteDomain {
    fn from(values: Vec<f64>) -> Self {
        Self { values }
    }
}
