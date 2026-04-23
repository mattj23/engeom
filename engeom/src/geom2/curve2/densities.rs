//! This module has some special density-like functions for 2D curves.

use crate::{Result, Curve2, Iso2, Series1};

impl Curve2 {
    /// Creates a series of the length density of the curve with respect to a given direction. This
    /// can be used for finding lines/edges in a curve.
    ///
    /// Effectively, the curve will be transformed by the inverse of the `facing` isometry. The
    /// algorithm will sweep from the minimum `x` value to the maximum `x` value in increments of
    /// `bin_width/2`, summing up the total `y` distance traveled by segments within `bin_width/2`
    /// of the x value being interrogated.
    ///
    /// The result will be a `Series1` where the domain is the x direction of the `facing` isometry
    /// and the range (y values) is the linear distance of the curve traveled in the y direction
    /// for the bin at each x value.
    ///
    /// # Arguments
    ///
    /// * `facing`:
    /// * `bin_width`:
    ///
    /// returns: Result<Series1, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn directional_density(&self, facing: &Iso2, bin_width: f64) -> Result<Series1> {

        todo!()
    }

}
