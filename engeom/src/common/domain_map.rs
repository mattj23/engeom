use serde::{Deserialize, Serialize};

/// This is an extremely simple linear mapping from one domain to another. It's constructed so that
/// 0.0 in domain A maps to `x0` in domain B, and 1.0 in domain A maps to `x0 + m` in domain B.
///
/// The purpose of this struct is to unambiguously handle a mapping between two domains so that
/// the calling code doesn't need to think about it.
#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct DomainMap {
    pub x0: f64,
    pub m: f64,
}

impl DomainMap {
    pub fn new(x0: f64, m: f64) -> Self {
        Self { x0, m }
    }

    /// Convert a value from A to B
    pub fn to(&self, a: f64) -> f64 {
        self.x0 + self.m * a
    }

    /// Convert a value from B to A
    pub fn from(&self, b: f64) -> f64 {
        (b - self.x0) / self.m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_domain_map_to_fwd() {
        let dm = DomainMap::new(0.5, 1.0);
        assert_relative_eq!(dm.to(0.1), 0.6);
    }

    #[test]
    fn test_domain_map_from_fwd() {
        let dm = DomainMap::new(0.5, 1.0);
        assert_relative_eq!(dm.from(0.6), 0.1);
    }

    #[test]
    fn test_domain_map_to_rev() {
        let dm = DomainMap::new(0.5, -1.0);
        assert_relative_eq!(dm.to(0.1), 0.4);
    }

    #[test]
    fn test_domain_map_from_rev() {
        let dm = DomainMap::new(0.5, -1.0);
        assert_relative_eq!(dm.from(0.4), 0.1);
    }
}
