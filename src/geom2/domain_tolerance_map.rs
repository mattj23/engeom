/// Contains an abstraction for working with a tolerance zone mapped onto a
/// discrete domain.
use crate::common::TolZone;

pub trait TolZoneMap {
    fn get(&self, x: f64) -> TolZone;
}

/// A tolerance zone map that returns a constant tolerance zone for all values of x.
pub struct ConstantTolZone {
    tol_zone: TolZone,
}

impl ConstantTolZone {
    pub fn new(tol_zone: TolZone) -> Self {
        Self { tol_zone }
    }
}

impl TolZoneMap for ConstantTolZone {
    fn get(&self, _x: f64) -> TolZone {
        self.tol_zone
    }
}
