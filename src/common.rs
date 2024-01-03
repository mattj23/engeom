mod angles;

pub use angles::{signed_compliment_2pi, AngleInterval, AngleDir, angle_to_2pi};

/// General purpose options for resampling data over a discrete domain.
pub enum Resample {
    /// Resample by a given number of points, evenly spaced over the domain
    ByCount(usize),

    /// Resample with a specific spacing between points, understanding that if the spacing does not
    /// divide evenly into the domain the end points may not be centered in the original domain
    BySpacing(f64),

    /// Resample with a maximum spacing between points. The number of points will be chosen
    /// automatically such that the entire domain is covered (as if `BySpacing` was used) but the
    /// spacing between points will not exceed the given value.
    ByMaxSpacing(f64),
}


/// General purpose options for smoothing data over a discrete domain.
pub enum Smoothing {
    /// A Gaussian filter with the given standard deviation, where the filter size is truncated to
    /// 3 standard deviations
    Gaussian(f64),

    /// A quadratic fit filter with the given window size. A quadratic polynomial is fit to items
    /// within the window, and the item is replaced with the value of the polynomial at the same
    /// position
    Quadratic(f64),

    /// A cubic fit filter with the given window size. A cubic polynomial is fit to items within
    /// the window, and the item is replaced with the value of the polynomial at the same position
    Cubic(f64),
}

/// Generate a vec of domain values which are linearly spaced between `start` and `end` and which
/// have a count of `count`. The first value will be `start` and the last value will be `end`.
///
/// # Arguments
///
/// * `start`: the starting value of the domain, inclusive
/// * `end`: the ending value of the domain, inclusive
/// * `count`: the total number of discrete, evenly spaced values in the domain
///
/// returns: Vec<f64, Global>
///
/// # Examples
///
/// ```
/// use engeom::common::linear_space;
/// let domain = linear_space(0.0, 1.0, 3);
/// assert_eq!(domain, vec![0.0, 0.5, 1.0]);
/// ```
pub fn linear_space(start: f64, end: f64, count: usize) -> Vec<f64> {
    let mut result = Vec::with_capacity(count);
    let step = (end - start) / (count - 1) as f64;
    for i in 0..count {
        result.push(start + i as f64 * step);
    }
    result
}

