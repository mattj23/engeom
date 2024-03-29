mod angles;
mod discrete_domain;
mod interval;
pub mod vecf64;

pub use angles::{angle_to_2pi, signed_compliment_2pi, AngleDir, AngleInterval};
pub use discrete_domain::{linear_space, DiscreteDomain};
pub use interval::Interval;

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
