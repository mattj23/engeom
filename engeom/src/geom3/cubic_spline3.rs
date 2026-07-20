use crate::common::cubic_spline::CubicSpline;

/// A cubic Bézier curve in 3D space, defined by four control points.
///
/// This is one of `engeom`'s 3D geometric primitives.
///
/// This is the three-dimensional specialization of the dimension-generic
/// [`CubicSpline`](CubicSpline); see that type for the shared constructors and queries (`new`,
/// `position`, `derivative`, `tangent`, `curvature`, `split`, `polyline`, and so on).
pub type CubicSpline3 = CubicSpline<3>;
