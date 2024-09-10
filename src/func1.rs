//! This module and its submodules provide an abstraction for numerical functions over a continuous
//! one dimensional domain of scalar f64 values.

use crate::common::DiscreteDomain;

mod common_functions;
mod polynomial;
mod series1;

pub use common_functions::{Gaussian1, Line1};
pub use polynomial::{Cubic, Polynomial, Quadratic, Quartic, Quintic};
pub use series1::Series1;

/// A function over a continuous one dimensional domain of scalar f64 values.
pub trait Func1 {
    /// Evaluate the function at the given value.
    ///
    /// # Arguments
    ///
    /// * `x`: the value at which to evaluate the function
    ///
    /// returns: f64
    fn f(&self, x: f64) -> f64;

    /// Evaluate the function at each value in the given discrete domain.  A naive default
    /// implementation is provided, but implementations which can make use of the ordering
    /// guarantees of the discrete domain should override this method.
    ///
    /// # Arguments
    ///
    /// * `xs`: the discrete domain over which to evaluate the function
    ///
    /// returns: Vec<f64, Global>
    fn fs(&self, xs: &DiscreteDomain) -> Vec<f64> {
        xs.iter().map(|x| self.f(*x)).collect()
    }
}
