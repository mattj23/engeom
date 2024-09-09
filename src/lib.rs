use std::error::Error;

pub mod common;
pub mod func1;
pub mod geom2;
pub mod geom3;
pub mod errors;
pub mod stats;
pub mod utility;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

// Extremely common angle tools
pub use common::{AngleDir, AngleInterval};

// Extremely common 2D types
pub use geom2::{Iso2, Point2, SvdBasis2, Vector2, SurfacePoint2, Curve2, CurveStation2, Circle2,
                Arc2};

// Extremely common 3D types
pub use geom3::{Iso3, Point3, Vector3, SurfacePoint3};

// Extremely common conversion tools
pub use common::{To2D, To3D};

#[cfg(test)]
mod tests {}
