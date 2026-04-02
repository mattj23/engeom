//! Representations for manifold positions and orientations in 3D space.

use crate::common::PCoords;
use crate::na::SVector;
use crate::{Point3, SurfacePoint3, UnitVec3};

#[derive(Clone, Debug)]
pub struct Manifold1Pos3 {
    /// The position of the point along the manifold length
    pub l: f64,

    /// The position of the point in full 3D cartesian space
    pub point: Point3,

    /// The direction of the positive manifold length in full 3D cartesian space at the current
    /// position
    pub direction: UnitVec3,
}

impl Manifold1Pos3 {
    pub fn new(l: f64, point: Point3, direction: UnitVec3) -> Self {
        Self {
            l,
            point,
            direction,
        }
    }

    pub fn as_surface_point(&self) -> SurfacePoint3 {
        SurfacePoint3::new(self.point, self.direction)
    }
}

impl PCoords<3> for Manifold1Pos3 {
    fn coords(&self) -> SVector<f64, 3> {
        self.point.coords
    }
}
