use std::error::Error;

pub mod common;
pub mod func1;
pub mod geom2;
pub mod geom3;
mod errors;

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

// Extremely common angle tools
pub use common::{AngleDir, AngleInterval};

// Extremely common 2D types
pub use geom2::{Iso2, Point2, SvdBasis2, Vector2, SurfacePoint2, Curve2, CurveStation2};

#[cfg(test)]
mod tests {}
