// This module has an abstraction for working with a tolerance zone in reference to a scalar
// value.
use crate::Result;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TolZone {
    /// The lower bound of the tolerance zone
    pub lower: f64,

    /// The upper bound of the tolerance zone
    pub upper: f64,
}

impl TolZone {
    /// Create a new tolerance zone with the given nominal value and bounds, without checking
    /// that the bounds are valid. You must ensure that `lower` <= `upper`.
    pub fn new_unchecked(lower: f64, upper: f64) -> Self {
        Self { lower, upper }
    }

    /// Create a new tolerance zone with the given nominal value and bounds, checking that the
    /// bounds are valid. Returns an error if the bounds are not valid. The bounds are valid if
    /// `lower` <= `upper`.
    pub fn new(lower: f64, upper: f64) -> Result<Self> {
        if lower <= upper {
            Ok(Self { lower, upper })
        } else {
            Err("Invalid tolerance zone bounds".into())
        }
    }

    /// Returns true if the given value is within the tolerance zone
    pub fn contains(&self, x: f64) -> bool {
        x >= self.lower && x <= self.upper
    }

    /// Returns the size of the tolerance zone
    pub fn size(&self) -> f64 {
        self.upper - self.lower
    }
}
