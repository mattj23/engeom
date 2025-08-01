pub mod align;
mod angles;
pub mod average;
mod convert_2d_3d;
mod discrete_domain;
mod domain_map;
mod index_mask;
pub mod indices;
mod interval;
pub mod kd_tree;
pub mod points;
pub mod poisson_disk;
pub mod surface_point;
pub mod svd_basis;
pub mod triangulation;
pub mod vec_f64;
mod voxel_downsample;
pub mod domain_window;

use crate::na::{Point, SVector};
pub use align::DistMode;
pub use angles::{
    AngleDir, AngleInterval, angle_in_direction, angle_signed_pi, angle_to_2pi,
    signed_compliment_2pi,
};
pub use convert_2d_3d::{To2D, To3D};
pub use discrete_domain::{DiscreteDomain, linear_space};
pub use domain_map::DomainMap;
pub use index_mask::IndexMask;
pub use interval::Interval;
pub use parry3d_f64::query::SplitResult;
pub use surface_point::{SurfacePoint, SurfacePointCollection};
pub use voxel_downsample::voxel_downsample;

/// A type alias for signed integer indices for rasters
pub type PointNI<const D: usize> = Point<i32, D>;

/// A type alias for a vector of signed integers, typically used for manipulating indices in
/// D-dimensional raster spaces.
pub type VectorNI<const D: usize> = SVector<i32, D>;

/// General purpose option for starting the selection of a set of items, either from everything,
/// nothing, a specific set of indices, or a bitmask.
#[derive(Debug, Clone)]
pub enum Selection {
    /// Start with no items selected. This is used to indicate that the selection should start with
    /// nothing selected, and then items can be selected or modified.
    None,

    /// Select all items in the set. This is used to indicate that the selection should start with
    /// everything selected, and then items can be deselected or modified.
    All,

    /// A specific set of indices to select. This is passed as a vector of indices and not as
    /// a reference to a slice because the selection will need to be able to own and modify
    /// the indices.
    Indices(Vec<usize>),

    /// A bitmask which indicates which items are selected. This is passed not as a reference
    /// because the selection will need to be able to own and modify the mask.
    Mask(IndexMask),
}

/// General purpose option for selecting or deselecting items from a set
#[derive(Debug, Clone, Copy)]
pub enum SelectOp {
    /// The items identified by the operation should be added to the existing selection
    Add,

    /// The items identified by the operation should be removed from the existing selection
    Remove,

    /// The items identified by the operation should be retained in the selection, while
    /// the rest of the selection is cleared
    KeepOnly,
}

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
    /// Use all samples and perform a least-squares minimization
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

/// A generic trait for points or point-like structures in D-dimensional space which provides a
/// generic way to access the coordinates of the point as a vector.
pub trait PCoords<const D: usize> {
    /// Returns the coordinates of the point as a vector.
    fn coords(&self) -> SVector<f64, D>;
}

impl<const D: usize> PCoords<D> for [f64; D] {
    fn coords(&self) -> SVector<f64, D> {
        SVector::from_column_slice(self)
    }
}

impl<const D: usize> PCoords<D> for SurfacePoint<D> {
    fn coords(&self) -> SVector<f64, D> {
        self.point.coords
    }
}

impl<const D: usize> PCoords<D> for Point<f64, D> {
    fn coords(&self) -> SVector<f64, D> {
        self.coords
    }
}

impl<const D: usize> PCoords<D> for SVector<f64, D> {
    fn coords(&self) -> SVector<f64, D> {
        *self
    }
}
