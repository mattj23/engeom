use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Debug)]
pub enum InvalidGeometry {
    NotEnoughPoints,
    GeometricOpFailed,
}

impl Display for InvalidGeometry {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Error for InvalidGeometry {}
