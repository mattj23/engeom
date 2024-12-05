// This module has an abstraction for working with a tolerance zone in reference to a scalar
// value.
use crate::Result;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Tolerance {
    /// The lower bound of the tolerance zone
    pub lower: f64,

    /// The upper bound of the tolerance zone
    pub upper: f64,
}

impl Tolerance {
    /// Create a new tolerance zone with the given nominal value and bounds, without checking
    /// that the bounds are valid. You must ensure that `lower` <= `upper`.
    pub fn new_unchecked(lower: f64, upper: f64) -> Self {
        Self { lower, upper }
    }

    /// Create a new tolerance zone with the given nominal value and bounds, checking that the
    /// bounds are valid. Returns an error if the bounds are not valid. The bounds are valid if
    /// `lower` <= `upper`.
    pub fn try_new(lower: f64, upper: f64) -> Result<Self> {
        if lower <= upper {
            Ok(Self { lower, upper })
        } else {
            Err("Invalid tolerance zone bounds".into())
        }
    }

    /// Create a new tolerance zone with the given nominal value and half width. The tolerance
    /// zone will extend from `center - half_width.abs()` to `center + half_width.abs()`.
    ///
    /// # Arguments
    ///
    /// * `center`: the nominal value of the tolerance zone
    /// * `half_width`: the half width of the tolerance zone. If a negative value is given, it
    /// will be converted to a positive value.
    ///
    /// returns: TolZone
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::metrology::Tolerance;
    /// let zone = Tolerance::symmetrical(0.0, 1.0);
    /// assert_eq!(zone.lower, -1.0);
    /// assert_eq!(zone.upper, 1.0);
    ///
    /// let zone = Tolerance::symmetrical(2.0, 1.0);
    /// assert_eq!(zone.lower, 1.0);
    /// assert_eq!(zone.upper, 3.0);
    /// ```
    pub fn symmetrical(center: f64, half_width: f64) -> Self {
        Self {
            lower: center - half_width.abs(),
            upper: center + half_width.abs(),
        }
    }

    /// Returns true if the given value is within the tolerance zone
    pub fn conforms(&self, x: f64) -> bool {
        x >= self.lower && x <= self.upper
    }

    /// Returns the size of the tolerance zone (upper - lower)
    pub fn size(&self) -> f64 {
        self.upper - self.lower
    }

    /// Returns the center of the tolerance zone, which is the value equidistant from the lower
    /// and upper bounds
    pub fn center(&self) -> f64 {
        (self.upper + self.lower) / 2.0
    }
}
