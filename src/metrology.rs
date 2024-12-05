pub mod line_profiles;
mod surface_deviation;
mod tolerance;
mod tolerance_map;

pub use tolerance::Tolerance;
pub use tolerance_map::{ConstantTolMap, DiscreteDomainTolMap, ToleranceMap};

pub type SurfaceDeviation2 = surface_deviation::SurfaceDeviation<2>;
pub type SurfaceDeviationSet2 = surface_deviation::SurfaceDeviationSet<2>;
pub type SurfaceDeviation3 = surface_deviation::SurfaceDeviation<3>;
pub type SurfaceDeviationSet3 = surface_deviation::SurfaceDeviationSet<3>;
