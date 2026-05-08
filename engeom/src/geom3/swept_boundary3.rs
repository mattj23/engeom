use crate::Iso3;
use crate::geom2::Boundary2;

pub struct ExtrudedBoundary3 {
    shape: Boundary2,
    start: Iso3,
    length: f64,
}

pub struct PolarBoundary3 {}
