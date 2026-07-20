use crate::common::cubic_spline::CubicSpline;
use crate::geom3::Aabb3;

/// A cubic Bézier curve in 3D space, defined by four control points.
///
/// This is one of `engeom`'s 3D geometric primitives.
///
/// This is the three-dimensional specialization of the dimension-generic
/// [`CubicSpline`](CubicSpline); see that type for the shared constructors and queries (`new`,
/// `position`, `derivative`, `tangent`, `curvature`, `split`, `polyline`, and so on). `aabb` is a
/// 3D-specific convenience wrapper over the generic
/// [`compute_bounds`](CubicSpline::compute_bounds).
pub type CubicSpline3 = CubicSpline<3>;

impl CubicSpline3 {
    /// Returns the axis-aligned bounding box of the curve, computed on demand.
    pub fn aabb(&self) -> Aabb3 {
        let (lo, hi) = self.compute_bounds();
        Aabb3::new(lo, hi)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point3;
    use approx::assert_relative_eq;

    #[test]
    fn aabb_matches_compute_bounds() {
        let c = CubicSpline::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, -1.0, 1.0),
            Point3::new(3.0, 0.0, 0.0),
        );
        let (lo, hi) = c.compute_bounds();
        let aabb = c.aabb();
        assert_relative_eq!(aabb.mins, lo, epsilon = 1e-12);
        assert_relative_eq!(aabb.maxs, hi, epsilon = 1e-12);
    }
}
