//! This module contains the trait and match type used to describe a 2D alignment target: a
//! stationary entity (such as a [`Curve2`] or [`Boundary2`]) that a moving set of points can be
//! aligned to by repeatedly projecting each point onto its closest position on the target.

use crate::UnitVec2;
use crate::common::{PCoords, SPCoords};
use crate::geom2::{Boundary2, Curve2, CurveStation2, Point2, Vector2};
use crate::na::{SVector, Unit};

/// The result of projecting a single point onto a [`SurfaceTarget2`], used as the correspondence
/// for that point during an alignment.
#[derive(Debug, Clone)]
pub struct AlignSurfMatch2 {
    /// The closest point on the target to the query point
    pub point: Point2,

    /// The outward-facing normal of the target at `point`
    pub normal: UnitVec2,

    /// Whether the closest point actually lies on the target's interior, as opposed to having
    /// clamped to an open target's endpoint. Always `true` for a closed target. See the target's
    /// own `find_align_match` implementation for the exact clamping rule.
    pub is_on: bool,

    /// A scalar weight for the correspondence, in `[0, 1]`, independent of `is_on`. A target may
    /// use this to de-weight regions of itself that it considers less reliable.
    ///
    /// This is a statement of intent ("care about this correspondence less"), distinct from
    /// [`AlignSurfMatch2::sigma`], which is a statement about measurement noise.
    pub weight: f64,

    /// The measurement uncertainty of the target at `point`, as a standard deviation in the units
    /// of the geometry. Zero means the target is treated as exact, which is the default and is
    /// correct for nominal/theoretical geometry.
    ///
    /// A target built from measured data should report the uncertainty interpolated to `point`,
    /// since the match rarely lands exactly on a vertex. The alignment combines this with the
    /// test point's own uncertainty in quadrature, `sqrt(test^2 + target^2)`, which is the
    /// variance of the difference of two independent measurements.
    ///
    /// This is treated as **isotropic**. Real scanner uncertainty is usually one-dimensional
    /// (depth along the sensor axis), and the statistically correct contribution to a residual
    /// measured along direction `d` would be `sigma * |u . d|` for an uncertainty axis `u`.
    /// Nothing currently records `u`, so the isotropic treatment stands in for it. Note that the
    /// approximation only ever under-trusts a point: on a surface at grazing incidence, depth
    /// noise displaces the point along the surface rather than through it, so its true
    /// normal-direction uncertainty is smaller than the scalar suggests.
    pub sigma: f64,
}

impl AlignSurfMatch2 {
    pub fn new(point: Point2, normal: UnitVec2, is_on: bool, weight: f64) -> Self {
        Self {
            point,
            normal,
            is_on,
            weight,
            sigma: 0.0,
        }
    }

    /// Returns a copy of this match carrying the given target-side measurement uncertainty. See
    /// [`AlignSurfMatch2::sigma`] for the semantics and the isotropy caveat.
    pub fn with_sigma(&self, sigma: f64) -> Self {
        Self {
            sigma,
            ..self.clone()
        }
    }
}

impl PCoords<2> for AlignSurfMatch2 {
    fn coords(&self) -> SVector<f64, 2> {
        self.point.coords()
    }
}

impl SPCoords<2> for AlignSurfMatch2 {
    fn normal(&self) -> Unit<SVector<f64, 2>> {
        self.normal
    }
}

impl Default for AlignSurfMatch2 {
    fn default() -> Self {
        Self {
            point: Point2::origin(),
            normal: Vector2::x_axis(),
            is_on: false,
            weight: 0.0,
            sigma: 0.0,
        }
    }
}

/// A stationary 2D entity that a set of points can be aligned to, by projecting each point onto
/// its closest position on the target.
///
/// This takes `&Point2` rather than `&impl PCoords<2>` so that the trait remains object-safe,
/// allowing `Box<dyn SurfaceTarget2>` / `Vec<Box<dyn SurfaceTarget2>>` for cases (such as a future
/// multi-entity alignment) where the set of targets isn't known at compile time. Every call site
/// already holds a concrete, already-transformed point, so nothing is lost in generality.
///
/// A target derived from measured rather than nominal geometry should populate
/// [`AlignSurfMatch2::sigma`] via [`AlignSurfMatch2::with_sigma`], interpolating its own
/// uncertainty to the match position. Neither `Curve2` nor `Boundary2` does so today: a `Curve2`
/// is a bare polyline with no attributes to draw from, and a `Boundary2` is analytic geometry for
/// which measurement uncertainty isn't meaningful.
pub trait SurfaceTarget2: Sync + Send {
    fn find_align_match(&self, p: &Point2) -> AlignSurfMatch2;
}

impl SurfaceTarget2 for Curve2 {
    fn find_align_match(&self, p: &Point2) -> AlignSurfMatch2 {
        let station = self.at_closest_to_point(p);
        let sp = station.surface_point();

        let is_on = self.is_closed() || !is_curve_endpoint(&station, self);
        AlignSurfMatch2::new(sp.point, sp.normal, is_on, 1.0)
    }
}

/// Returns true if `station` sits exactly at one of the two open ends of `curve`.
///
/// `Curve2::at_closest_to_point` clamps to the nearest edge, so a projection past either end of
/// an open curve always lands with `fraction` exactly `0.0` or `1.0` at the first or last edge
/// respectively; there is no tolerance band to worry about.
fn is_curve_endpoint(station: &CurveStation2, curve: &Curve2) -> bool {
    let last_edge = curve.count() - 2;
    (station.index() == 0 && station.fraction() == 0.0)
        || (station.index() == last_edge && station.fraction() == 1.0)
}

impl SurfaceTarget2 for Boundary2 {
    fn find_align_match(&self, p: &Point2) -> AlignSurfMatch2 {
        let (_, m) = self.at_closest_to_point(p);

        // As with `Curve2`, every `BoundaryElement2::closest_to_point` implementation clamps
        // exactly to its own length bounds, so a projection past either end of an open boundary
        // lands with the cumulative length exactly `0.0` or exactly `self.length()`.
        let is_on = self.is_closed() || (m.l > 0.0 && m.l < self.length());

        AlignSurfMatch2::new(m.point, m.normal, is_on, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{BoundaryData2, BoundaryEditor};
    use approx::assert_relative_eq;

    fn ccw_square_curve(closed: bool) -> Curve2 {
        let points =
            [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)].map(|(x, y)| Point2::new(x, y));
        Curve2::from_points(&points, 1e-8, closed).unwrap()
    }

    #[test]
    fn curve_closed_normal_points_outward() {
        let curve = ccw_square_curve(true);
        let m = curve.find_align_match(&Point2::new(1.0, -5.0));
        assert_relative_eq!(
            m.normal.into_inner(),
            Vector2::new(0.0, -1.0),
            epsilon = 1e-8
        );
        assert!(m.is_on);
    }

    #[test]
    fn curve_closed_seam_is_on() {
        // The closest point to (0.0, -5.0) is the seam vertex (0,0)/(0,2)... rather pick a point
        // whose closest station is exactly the shared start/end vertex of a closed curve.
        let curve = ccw_square_curve(true);
        let m = curve.find_align_match(&Point2::new(-5.0, 0.0));
        assert!(m.is_on);
    }

    #[test]
    fn curve_open_end_is_not_on() {
        let curve = ccw_square_curve(false);
        // Past the first vertex (0,0) of the open curve
        let m = curve.find_align_match(&Point2::new(-5.0, -5.0));
        assert!(!m.is_on);

        // Past the last vertex (0,2) of the open curve
        let m = curve.find_align_match(&Point2::new(-5.0, 7.0));
        assert!(!m.is_on);
    }

    #[test]
    fn curve_open_interior_is_on() {
        let curve = ccw_square_curve(false);
        let m = curve.find_align_match(&Point2::new(1.0, -5.0));
        assert!(m.is_on);
    }

    fn ccw_square_boundary(closed: bool) -> Boundary2 {
        let mut data = if closed {
            BoundaryData2::new_closed()
        } else {
            BoundaryData2::new_open(Point2::new(0.0, 0.0))
        };

        let mut cursor = data.get_cursor(None);
        cursor.add_seg_xy(2.0, 0.0);
        cursor.add_seg_xy(2.0, 2.0);
        cursor.add_seg_xy(0.0, 2.0);
        if closed {
            cursor.add_seg_xy(0.0, 0.0);
        }

        data.try_to_boundary().unwrap()
    }

    #[test]
    fn boundary_closed_normal_points_outward() {
        let boundary = ccw_square_boundary(true);
        let m = boundary.find_align_match(&Point2::new(1.0, -5.0));
        assert_relative_eq!(
            m.normal.into_inner(),
            Vector2::new(0.0, -1.0),
            epsilon = 1e-8
        );
        assert!(m.is_on);
    }

    #[test]
    fn boundary_open_end_is_not_on() {
        let boundary = ccw_square_boundary(false);
        let m = boundary.find_align_match(&Point2::new(-5.0, -5.0));
        assert!(!m.is_on);

        let m = boundary.find_align_match(&Point2::new(-5.0, 7.0));
        assert!(!m.is_on);
    }

    #[test]
    fn boundary_open_interior_is_on() {
        let boundary = ccw_square_boundary(false);
        let m = boundary.find_align_match(&Point2::new(1.0, -5.0));
        assert!(m.is_on);
    }
}
