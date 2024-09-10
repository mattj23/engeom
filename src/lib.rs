use std::error::Error;

pub mod airfoil;
pub mod common;
pub mod errors;
pub mod func1;
pub mod geom2;
pub mod geom3;
pub mod stats;
pub mod utility;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

// Extremely common angle tools
pub use common::{AngleDir, AngleInterval};

// Extremely common 2D types
pub use geom2::{
    Arc2, Circle2, Curve2, CurveStation2, Iso2, Point2, SurfacePoint2, SvdBasis2, Vector2,
};

// Extremely common 3D types
pub use geom3::{Iso3, Point3, SurfacePoint3, Vector3};

// Extremely common conversion tools
pub use common::{To2D, To3D};

#[cfg(test)]
mod tests {}
