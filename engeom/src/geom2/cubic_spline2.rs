use crate::AngleDir::Cw;
use crate::UnitVec2;
use crate::common::cubic_spline::{CubicSpline, CubicSplineLookup};
use crate::geom2::rot90;

pub type CubicSpline2 = CubicSpline<2>;
pub type CubicSplineLookup2 = CubicSplineLookup<2>;

impl CubicSpline2 {
    /// Return the unit normal vector of the curve at parameter `t`. The normal vector is the
    /// tangent vector rotated clockwise by 90 degrees.
    pub fn normal(&self, t: f64) -> UnitVec2 {
        rot90(Cw) * self.tangent(t)
    }
}
