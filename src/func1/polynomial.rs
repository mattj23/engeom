//! This module contains a generic implementation of a polynomial function of a single variable (x)
//! and the trait implementations for a `Func1`.  Common polynomials (quadratic, cubic, quartic,
//! quintic) are provided as type aliases.  A generic least squares fitting algorithm for any
//! polynomial is also provided.

use super::Func1;
use parry2d_f64::na::DMatrix;

pub type Quadratic = Polynomial<3>;
pub type Cubic = Polynomial<4>;
pub type Quartic = Polynomial<5>;
pub type Quintic = Polynomial<6>;

/// A polynomial function of a single variable (x), where the total number of exponents is
/// specified by the generic parameter `K`.  Because x^0 counts as an exponent, the highest
/// exponent is K-1.  For example, a quadratic polynomial will be defined as a `Polynomial<3>`,
/// where the exponents are x^0, x^1, and x^2.
///
/// A polynomial can be constructed with either a known set of coefficients, or by fitting a set
/// of data points using a least squares algorithm capable of using weights.
#[derive(Debug, Clone)]
pub struct Polynomial<const K: usize> {
    /// The coefficients of the polynomial, ordered from lowest exponent to highest exponent.
    pub c: [f64; K],
}

impl<const K: usize> Polynomial<K> {
    /// Construct a new polynomial with the given coefficients.  The coefficients are ordered from
    /// lowest exponent to highest exponent.  For example, a quadratic polynomial will be defined
    /// as a `Polynomial<3>`, where the coefficients are ordered as c[0], c[1], and c[2], in which
    /// c[0] is the coefficient for x^0, c[1] is the coefficient for x^1, and c[2] is the
    /// coefficient for x^2.
    ///
    /// # Arguments
    ///
    /// * `c`: a slice of f64 values representing the coefficients of the polynomial
    ///
    /// returns: Polynomial<{ K }>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::func1::{Func1, Polynomial};
    /// let p = Polynomial::new([1.0, 2.0, 3.0]);
    /// let y = p.f(2.0);
    /// assert_eq!(y, 17.0);
    /// ```
    pub fn new(c: [f64; K]) -> Self {
        Self { c }
    }

    /// Perform a least squares fit of a polynomial to the given data points, optionally using
    /// weights.  This is a generic implementation of the standard least squares calculation for
    /// a polynomial of any degree.  The residuals are calculated as the difference in the y values
    /// of the sample points and the model values at the same x value.
    ///
    /// The weights are optional, but if provided, must be the same length as the sample points.
    /// The sample points are specified as two separate slices of the same length where x and y
    /// values of each sample point occupy the same index in their respective slices.
    ///
    /// X and y slices of different lengths will cause a panic.
    ///
    /// # Arguments
    ///
    /// * `xs`: a slice of f64 values representing the x values of the sample points
    /// * `ys`: a slice of f64 values representing the y values of the sample points
    /// * `weights`: an optional slice of f64 values representing the weights of the sample points
    ///
    /// returns: Polynomial<{ K }>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::func1::{Func1, Polynomial};
    /// let xs = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
    /// let ys = vec![1.0, 0.25, 0.0, 0.25, 1.0];
    /// let p = Polynomial::<3>::least_squares(&xs, &ys, None);
    ///
    /// assert_eq!(p.c[0], 0.0);
    /// assert_eq!(p.c[1], 0.0);
    /// assert_eq!(p.c[2], 1.0);
    /// ```
    pub fn least_squares(xs: &[f64], ys: &[f64], weights: Option<&[f64]>) -> Self {
        assert_eq!(xs.len(), ys.len());
        let w = Weights::new(weights);

        let mut sums = vec![0.0; 2 * K + 1];
        let mut rhs: DMatrix<f64> = DMatrix::zeros(K, 1);
        for i in 0..xs.len() {
            let w = w.get(i);

            for k in 0..K {
                let wxk = w * xs[i].powi(k as i32);
                rhs[(k, 0)] += wxk * ys[i];
                sums[k] += wxk;
            }
            for k in K + 1..=2 * K {
                let xk = xs[i].powi(k as i32);
                sums[k] += w * xk;
            }
        }

        let mut matrix = DMatrix::zeros(K, K);
        for r in 0..K {
            for c in 0..K {
                matrix[(r, c)] = sums[r + c];
            }
        }

        let p = matrix.try_inverse().unwrap() * rhs;
        let mut c = [0.0; K];
        for i in 0..K {
            c[i] = p[(i, 0)];
        }
        Self::new(c)
    }
}

impl<const K: usize> Func1 for Polynomial<K> {
    fn f(&self, x: f64) -> f64 {
        let mut y = 0.0;
        for i in 0..K {
            y += self.c[i] * x.powi(i as i32);
        }
        y
    }
}

/// Private helper struct for handling weights in a least squares fit.
struct Weights<'a> {
    values: Option<&'a [f64]>,
}

impl<'a> Weights<'a> {
    fn new(values: Option<&'a [f64]>) -> Self {
        Self { values }
    }

    fn get(&self, i: usize) -> f64 {
        if let Some(values) = self.values {
            values[i]
        } else {
            1.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::DiscreteDomain;
    use approx::assert_relative_eq;

    #[test]
    fn least_squares_quadratic() {
        let origin = Quadratic::new([0.0, 0.0, 1.0]);
        let x = DiscreteDomain::linear(-1.0, 1.0, 5);
        let y = origin.fs(&x);

        let fit = Quadratic::least_squares(&x, &y, None);

        assert_relative_eq!(origin.c[0], fit.c[0], epsilon = 1e-10);
        assert_relative_eq!(origin.c[1], fit.c[1], epsilon = 1e-10);
        assert_relative_eq!(origin.c[2], fit.c[2], epsilon = 1e-10);
    }

    #[test]
    fn test_quartic() {
        let origin = Quartic::new([5.0, -2.0, 1.0, 4.0, -3.0]);
        let x = DiscreteDomain::linear(-1.0, 1.0, 1000);
        let y = origin.fs(&x);

        let fit = Quartic::least_squares(&x, &y, None);

        assert_relative_eq!(origin.c[0], fit.c[0], epsilon = 1e-10);
        assert_relative_eq!(origin.c[1], fit.c[1], epsilon = 1e-10);
        assert_relative_eq!(origin.c[2], fit.c[2], epsilon = 1e-10);
        assert_relative_eq!(origin.c[3], fit.c[3], epsilon = 1e-10);
        assert_relative_eq!(origin.c[4], fit.c[4], epsilon = 1e-10);
    }
}
