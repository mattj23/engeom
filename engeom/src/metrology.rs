mod dimension;
pub mod line_profiles;
mod surface_deviation;
mod tolerance;
mod tolerance_map;

pub use tolerance::Tolerance;
pub use tolerance_map::{ConstantTolMap, DiscreteDomainTolMap, ToleranceMap};

use crate::{Iso3, To2D, To3D};
pub use dimension::Measurement;

pub type SurfaceDeviation2 = surface_deviation::SurfaceDeviation<2>;
pub type SurfaceDeviationSet2 = surface_deviation::SurfaceDeviationSet<2>;
pub type SurfaceDeviation3 = surface_deviation::SurfaceDeviation<3>;
pub type SurfaceDeviationSet3 = surface_deviation::SurfaceDeviationSet<3>;

pub type Distance2 = dimension::Distance<2>;
pub type Distance3 = dimension::Distance<3>;

// Conversions between 2D and 3D distances
impl Distance2 {
    /// Convert a 2D distance to a 3D distance using an isometry transformation. The 2D distance
    /// is converted to 3D by adding a zero z-component, and then it is transformed by the provided
    /// isometry to move it to some other location in 3D space.
    ///
    /// # Arguments
    ///
    /// * `iso`: The isometry transformation to apply to the 2D distance
    ///
    /// returns: Distance<3>
    pub fn to_3d(&self, iso: &Iso3) -> Distance3 {
        let a = iso * self.a.to_3d();
        let b = iso * self.b.to_3d();
        let direction = iso * self.direction.to_3d();
        Distance3::new(a, b, Some(direction))
    }
}

impl Distance3 {
    /// Convert a 3D distance to a 2D distance using an isometry transformation. The 3D distance
    /// is first transformed by the provided isometry, and then it is converted to 2D by dropping
    /// the z-component.
    ///
    /// # Arguments
    ///
    /// * `iso`: The isometry transformation to apply to the 3D distance before dropping the
    ///   z-component
    ///
    /// returns: Distance<2>
    pub fn to_2d(&self, iso: &Iso3) -> Distance2 {
        let a = (iso * self.a).to_2d();
        let b = (iso * self.b).to_2d();
        let direction = (iso * self.direction).to_2d();
        Distance2::new(a, b, Some(direction))
    }
}
