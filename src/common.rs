pub mod align;
mod angles;
mod convert_2d_3d;
mod discrete_domain;
pub mod indices;
mod interval;
pub mod points;
pub mod surface_point;
pub mod svd_basis;
pub mod vec_f64;

pub use align::DistMode;
pub use angles::{
    angle_in_direction, angle_signed_pi, angle_to_2pi, signed_compliment_2pi, AngleDir,
    AngleInterval,
};
pub use convert_2d_3d::{To2D, To3D};
pub use discrete_domain::{linear_space, DiscreteDomain};
pub use interval::Interval;
pub use parry3d_f64::query::SplitResult;
pub use surface_point::{SurfacePoint, SurfacePointCollection};

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

/// General purpose options for fitting data to a model
#[derive(Debug, Clone, Copy)]
pub enum BestFit {
    /// Use all samples and perform a least squares minimization
    All,

    /// De-weight samples based on their standard deviation from the mean
    Gaussian(f64),
}

/// A trait for projecting an entity to another entity
pub trait Project<TEntity, TResult> {
    fn project(&self, entity: TEntity) -> TResult;
}

/// A trait for intersecting an entity with another entity
pub trait Intersection<TOther, TResult> {
    fn intersection(&self, other: TOther) -> TResult;
}

/// A trait for transforming an entity by another entity
pub trait TransformBy<T, TOut> {
    fn transform_by(&self, transform: &T) -> TOut;
}
