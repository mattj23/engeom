//! This module contains a set of common mathematical functions which can represent continuous
//! functions over the x domain.  Excluded from this set are polynomials, which are provided in
//! the `polynomial` module.

use crate::func1::Func1;
use crate::Result;

/// A `Line1` is a linear function in one dimension, defined by slope `m` and intercept `b`.  The
/// `Line1` type is an alias for a `Polynomial<2>` in order to take advantage of the generic
/// least squares algorithm.
pub type Line1 = super::polynomial::Polynomial<2>;

impl Line1 {
    /// Create a new `Line1` with the given slope and intercept.
    ///
    /// # Arguments
    ///
    /// * `m`: the slope of the line
    /// * `b`: the y intercept of the line
    ///
    /// returns: Polynomial<2>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::func1::{Func1, Line1};
    /// let line = Line1::new_mxb(2.0, 1.0);
    /// let y = line.f(2.0);
    /// assert_eq!(y, 5.0);
    /// ```
    pub fn new_mxb(m: f64, b: f64) -> Self {
        Self::new([b, m])
    }

    /// Try to create a new `Line1` from two points.  If the two points are too close together, the
    /// slope would be infinite, and an `Err` is returned instead.
    ///
    /// # Arguments
    ///
    /// * `x0`: the x value of the first point
    /// * `y0`: the y value of the first point
    /// * `x1`: the x value of the second point
    /// * `y1`: the y value of the second point
    ///
    /// returns: Result<Polynomial<2>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::func1::{Func1, Line1};
    /// let line1 = Line1::try_from_points(0.0, 1.0, 1.0, 3.0);
    /// assert!(line1.is_ok());
    ///
    /// let line2 = Line1::try_from_points(0.0, 1.0, 0.0, 3.0);
    /// assert!(line2.is_err());
    /// ```
    pub fn try_from_points(x0: f64, y0: f64, x1: f64, y1: f64) -> Result<Self> {
        if (x1 - x0).abs() < 1e-12 {
            return Err(
                "x1 and x0 are too close together and will create an infinite slope".into(),
            );
        }
        let m = (y1 - y0) / (x1 - x0);
        let b = y0 - m * x0;
        Ok(Self::new_mxb(m, b))
    }

    /// Returns the slope (m) of the line.
    pub fn m(&self) -> f64 {
        self.c[1]
    }

    /// Returns the y intercept (b) of the line.
    pub fn b(&self) -> f64 {
        self.c[0]
    }
}

/// A `Gaussian1D` is a Gaussian function in one dimension, defined by mean `mean` and standard
/// deviation `sigma`.  The gaussian is a pure implementation of the function, and has a maximum
/// value of 1.0 at `mean`.
pub struct Gaussian1 {
    pub mean: f64,
    pub sigma: f64,
}

impl Gaussian1 {
    /// Create a new `Gaussian1` with the given mean and standard deviation.  The value of the
    /// function at the mean will be 1.0.
    ///
    /// # Arguments
    ///
    /// * `mean`: the mean of the gaussian function
    /// * `sigma`: the standard deviation of the gaussian function
    ///
    /// returns: Gaussian1
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::func1::{Func1, Gaussian1};
    /// let gaussian = Gaussian1::new(0.0, 1.0);
    /// assert_relative_eq!(gaussian.f(0.0), 1.0);
    /// assert_relative_eq!(gaussian.f(1.0), 0.6065306597126334);
    /// ```
    pub fn new(mean: f64, sigma: f64) -> Self {
        Self { mean, sigma }
    }
}

impl Func1 for Gaussian1 {
    fn f(&self, x: f64) -> f64 {
        (-0.5 * ((x - self.mean).powi(2) / self.sigma.powi(2))).exp()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::func1::Func1;
    use approx::assert_relative_eq;

    #[test]
    fn test_line1() {
        let line = Line1::new_mxb(2.0, 1.0);
        assert_relative_eq!(line.f(0.0), 1.0);
        assert_relative_eq!(line.f(1.0), 3.0);
        assert_relative_eq!(line.f(2.0), 5.0);
    }
}
