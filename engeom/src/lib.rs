extern crate core;

use std::error::Error;

pub mod airfoil2;
pub mod common;
pub mod errors;
pub mod func1;
pub mod geom2;
pub mod geom3;
pub mod io;
pub mod metrology;
pub mod raster2;
pub mod raster3;
pub mod sensors;
pub mod stats;

#[cfg(feature = "three_d")]
pub mod td;
pub mod utility;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;
pub type ResultCode<T> = std::result::Result<T, usize>;

// Re-export some commonly used crates for convenience
pub use alum;
pub use colorgrad;
pub use imageproc;
pub use imageproc::image;
pub use levenberg_marquardt;
pub use parry2d_f64 as parry2d;
pub use parry3d_f64 as parry3d;
pub use parry3d_f64::na;
pub use rayon;
pub use serde;
pub use serde_json;
// Re-exported so a caller can reach `tol_compress::reorder` without taking the dependency and
// matching versions themselves. Writing a tcmesh renumbers vertices, and that is the module which
// says how, so anyone holding per-vertex data outside the file needs it.
pub use tol_compress;

// Re-export the `three_d` crate if the feature is enabled
#[cfg(feature = "three_d")]
pub use three_d;

// Common one-dimensional functions
pub use func1::{Func1, Gaussian1, Line1, Polynomial, Series1};

// Extremely common angle tools
pub use common::{AngleDir, AngleInterval, IndexMask, IntervalOps};

// Nalgebra exports
pub type DVector = na::DVector<f64>;
pub type DMatrix = na::DMatrix<f64>;

// Extremely common 2D types
pub use geom2::{
    Arc2, Circle2, Curve2, CurveGroup2, CurveStation2, Iso2, KdTree2, Line2, Point2, SurfacePoint2,
    SvdBasis2, UnitVec2, Vector2,
};

// Extremely common 3D types
pub use geom3::{
    Aabb3, CloudIndex3, Curve3, CurveStation3, HalfEdgeMesh3, Iso3, KdTree3, Line3, Manifold1Pos3,
    Mesh3, MeshData3, Plane3, Point3, PointCloud3, PointCloudOverlap, Sphere3, SurfacePoint3,
    SvdBasis3, UnitVec3, VOXEL_COHERENCE_ATTR, VOXEL_COUNT_ATTR, Vector3,
};

// Extremely common conversion tools
pub use common::{To2D, To3D};

/// General purpose option for how to handle the result of a dot product between directional
/// vectors.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum VecDot {
    /// Use the dot product as is, allowing it to range from any negative to any positive number
    AsIs,

    /// Use the absolute value of the dot product
    Abs,

    /// Clamp the dot product to a positive value. Values below 0 are raised to zero.
    ClampPos,
}

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

#[cfg(test)]
pub mod tests {
    use crate::io::{read_tc_curve2_from, read_tc_mesh_from};

    use crate::{Curve2, Mesh3};

    /// Load a mesh with the stanford bunny reconstruction at resolution 4. The vertices are within
    /// 0.00000189 of the original ply file.
    pub fn stanford_bun_4() -> Mesh3 {
        embedded_mesh(include_bytes!("../tests/data/stanford_bun_4.tcmesh"))
    }

    /// Load a mesh with the stanford bunny reconstruction at resolution 2. The vertices are within
    /// 0.00000189 of the original ply file.
    pub fn stanford_bun_2() -> Mesh3 {
        embedded_mesh(include_bytes!("../tests/data/stanford_bun_2.tcmesh"))
    }

    /// Load a mesh with the stanford bunny reconstruction at resolution 3. The vertices are within
    /// 0.00000189 of the original ply file.
    pub fn stanford_bun_3() -> Mesh3 {
        embedded_mesh(include_bytes!("../tests/data/stanford_bun_3.tcmesh"))
    }

    /// Load a mesh of a small engine blade. The mesh has 21795 vertices and 43586 faces. Dimensions
    /// are in millimeters. The mesh has been processed externally to remove mesh errors and should
    /// be watertight. Stored at a 1 micron tolerance.
    pub fn engine_blade() -> Mesh3 {
        embedded_mesh(include_bytes!("../tests/data/engine-blade.tcmesh"))
    }

    fn embedded_mesh(bytes: &[u8]) -> Mesh3 {
        let data = read_tc_mesh_from(&mut { bytes }).unwrap();
        Mesh3::from_data(data, false).unwrap()
    }

    pub fn airfoil_curve() -> Curve2 {
        let bytes = include_bytes!("../tests/data/airfoil-0.tccurve2");
        let mut working = bytes.as_slice();
        read_tc_curve2_from(&mut working).unwrap()
    }

    /// Get the path to the test data directory.
    pub fn get_test_data_path() -> std::path::PathBuf {
        let mut path = std::env::current_dir().unwrap();
        path.push("tests");
        path.push("data");
        path
    }

    /// Get the path to a test file in the test data directory. This does not check if the file
    /// exists.
    pub fn get_test_file_path(str: &str) -> std::path::PathBuf {
        let mut path = get_test_data_path();
        path.push(str);
        path
    }
}
