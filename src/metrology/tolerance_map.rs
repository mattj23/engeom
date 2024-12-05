use crate::common::DiscreteDomain;
use crate::metrology::Tolerance;
use crate::Result;

pub trait ToleranceMap {
    fn get(&self, x: f64) -> Option<Tolerance>;
}

/// A tolerance zone map that returns a constant tolerance zone for all values of x.
pub struct ConstantTolMap {
    tol_zone: Tolerance,
}

impl ConstantTolMap {
    pub fn new(tol_zone: Tolerance) -> Self {
        Self { tol_zone }
    }
}

impl ToleranceMap for ConstantTolMap {
    fn get(&self, _x: f64) -> Option<Tolerance> {
        Some(self.tol_zone)
    }
}

/// A tolerance zone map that returns a tolerance zone based on a 1D discrete domain. Each tolerance
/// zone is associated with a value in the domain and extends to the next value in the domain, with
/// the last value in the domain extending to infinity.
pub struct DiscreteDomainTolMap {
    pub domain: DiscreteDomain,
    pub tol_zones: Vec<Tolerance>,
}

impl DiscreteDomainTolMap {
    /// Create a new tolerance zone map from a discrete domain and a list of tolerance zones. The
    /// number of values in the domain must be the same as the number of tolerance zones. If the
    /// number of values and tolerance zones are not the same, an error is returned.
    ///
    /// # Arguments
    ///
    /// * `domain`: a discrete domain where each value is associated with a tolerance zone by index
    /// * `zones`: a vector of tolerance zones, one for each value in the domain
    ///
    /// returns: Result<DiscreteDomainTolZone, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::common::DiscreteDomain;
    /// use engeom::metrology::{Tolerance, ToleranceMap, DiscreteDomainTolMap};
    /// use approx::assert_relative_eq;
    ///
    /// let domain = DiscreteDomain::try_from(vec![1.0, 2.0, 3.0]).unwrap();
    /// let zones = vec![
    ///     Tolerance::symmetrical(0.0, 1.0),
    ///     Tolerance::symmetrical(0.0, 2.0),
    ///     Tolerance::symmetrical(0.0, 3.0)
    /// ];
    ///
    /// let map = DiscreteDomainTolMap::try_new(domain, zones).unwrap();
    ///
    /// let zone = map.get(2.5).unwrap();
    /// assert_relative_eq!(zone.size(), 4.0);
    /// ```
    pub fn try_new(domain: DiscreteDomain, tol_zones: Vec<Tolerance>) -> Result<Self> {
        if domain.len() != tol_zones.len() {
            Err("The number of values and tolerance zones must be the same".into())
        } else {
            Ok(Self { domain, tol_zones })
        }
    }
}

impl ToleranceMap for DiscreteDomainTolMap {
    fn get(&self, x: f64) -> Option<Tolerance> {
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
