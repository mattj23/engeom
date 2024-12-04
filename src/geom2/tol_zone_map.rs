use crate::common::{DiscreteDomain, TolZone};
/// Contains an abstraction for working with a tolerance zone mapped onto a
/// discrete domain.
use crate::Result;

pub trait TolZoneMap {
    fn get(&self, x: f64) -> Option<TolZone>;
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
    fn get(&self, _x: f64) -> Option<TolZone> {
        Some(self.tol_zone)
    }
}

/// A tolerance zone map that returns a tolerance zone based on a 1D discrete domain. Each tolerance
/// zone is associated with a value in the domain and extends to the next value in the domain, with
/// the last value in the domain extending to infinity.
pub struct DiscreteDomainTolZone {
    domain: DiscreteDomain,
    tol_zones: Vec<TolZone>,
}

impl DiscreteDomainTolZone {
    /// Create a new tolerance zone map from a list of values and tolerance zones. The values must
    /// be in increasing order and the tolerance zones must be valid. Each `TolZone` is considered
    /// to start at the corresponding value in the domain and extend to the next value in the domain.
    /// The last value in the domain is considered to extend to infinity.
    pub fn try_new(values: &[(f64, TolZone)]) -> Result<Self> {
        let mut domain = DiscreteDomain::default();
        let mut tol_zones = Vec::with_capacity(values.len());
        for (value, tol_zone) in values {
            domain.push(*value)?;
            tol_zones.push(*tol_zone);
        }
        Ok(Self { domain, tol_zones })
    }
}

impl TolZoneMap for DiscreteDomainTolZone {
    fn get(&self, x: f64) -> Option<TolZone> {
        if self.domain.is_empty() {
            None
        } else {
            if let Some(i) = self.domain.index_of(x) {
                Some(self.tol_zones[i])
            } else {
                Some(self.tol_zones[self.tol_zones.len() - 1])
            }
        }
    }
}
