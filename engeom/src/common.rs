pub mod align;
mod angles;
pub mod average;
mod convert_2d_3d;
pub mod cubic_spline;
mod discrete_domain;
mod domain_map;
pub mod domain_window;
mod index_mask;
pub mod indices;
pub mod interval;
pub mod kd_tree;
pub mod points;
pub mod poisson_disk;
pub mod ransac_tools;
pub mod surface_point;
pub mod svd_basis;
pub mod triangulation;
pub mod vec_f64;
mod voxel_downsample;

use crate::na::{Point, SVector, Unit};
pub use align::DistMode;
pub use angles::*;
pub use average::Averager;
pub use convert_2d_3d::{To2D, To3D};
pub use discrete_domain::{DiscreteDomain, linear_space};
pub use domain_map::DomainMap;
pub use index_mask::IndexMask;
pub use interval::{AngleInterval, Interval, IntervalMergeDomain, IntervalOps};
pub use parry3d_f64::query::SplitResult;
pub use points::*;
pub use surface_point::{SurfacePoint, SurfacePointCollection};
pub use voxel_downsample::voxel_downsample;

/// A type alias for signed integer indices for rasters
pub type PointNI<const D: usize> = Point<i32, D>;

/// A type alias for a vector of signed integers, typically used for manipulating indices in
/// D-dimensional raster spaces.
pub type VectorNI<const D: usize> = SVector<i32, D>;

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

pub trait SPCoords<const D: usize>: PCoords<D> {
    fn normal(&self) -> Unit<SVector<f64, D>>;

    fn scalar_projection(&self, p: &impl PCoords<D>) -> f64 {
        self.normal().dot(&(p.coords() - self.coords()))
    }
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

impl<const D: usize> SPCoords<D> for SurfacePoint<D> {
    fn normal(&self) -> Unit<SVector<f64, D>> {
        self.normal
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
