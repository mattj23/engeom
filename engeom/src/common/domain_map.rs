use serde::{Deserialize, Serialize};

/// This is an extremely simple linear mapping from one domain to another. It's constructed so that
/// 0.0 in the outer domain maps to x0 in the inner domain, and 1.0 in the outer domain maps to
/// x0 + m in the inner domain.
///
/// The purpose of this struct is to unambiguously handle a mapping between two domains so that
/// the calling code doesn't need to think about it.
/// TODO: Clarify the naming and provide some examples
#[derive(Debug, Copy, Clone, Serialize, Deserialize)]
pub struct DomainMap {
    pub x0: f64,
    pub m: f64,
}

impl DomainMap {
    pub fn new(x0: f64, m: f64) -> Self {
        Self { x0, m }
    }

    pub fn to(&self, x: f64) -> f64 {
        self.x0 + self.m * x
    }

    pub fn from(&self, x: f64) -> f64 {
        (x - self.x0) / self.m
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
