use std::error::Error;

pub mod airfoil;
pub mod common;
pub mod errors;
pub mod func1;
pub mod geom2;
pub mod geom3;
pub mod io;
pub mod metrology;
pub mod stats;
pub mod utility;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

// Common one dimensional functions
pub use func1::{Func1, Gaussian1, Line1, Polynomial, Series1};

// Extremely common angle tools
pub use common::{AngleDir, AngleInterval};

// Extremely common 2D types
pub use geom2::{
    Arc2, Circle2, Curve2, CurveStation2, Iso2, Point2, SurfacePoint2, SvdBasis2, UnitVec2, Vector2,
};

// Extremely common 3D types
pub use geom3::{
    Curve3, CurveStation3, Iso3, Mesh, Plane3, Point3, SurfacePoint3, UnitVec3, Vector3, SvdBasis3,
};

// Extremely common conversion tools
pub use common::{To2D, To3D, TransformBy};

#[cfg(test)]
mod tests {}
