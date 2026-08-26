//! This module contains support functionality for curve-to-curve alignment, including the core
//! building blocks consumed by the other curve alignment routines in the `align2` module tree. All
//! public items are re-exported directly from `align2`, so callers interact with them through that
//! parent module. If you are implementing your own alignment routines, you may find that you want
//! to use some functionality contained here.
//!
//! This is the 2D counterpart of `geom3::align3::mesh`, and is deliberately structured similarly
//! to it. The differences that are not simply a matter of dimension are called out below.
//!
//! # What this module provides
//!
//! ## Sample point container
//!
//! [`CurveSurfPoint`] is an owned sample taken from a `Curve2`, the counterpart of `MeshSurfPoint`
//! in 3D. It exists because `CurveStation2` borrows its parent curve and so cannot be stored
//! across the transforms a multi-body solve applies.
//!
//! ## Alignment curve container
//!
//! [`AlignmentCurve`] wraps a `&Curve2` together with optional per-vertex uncertainty values, an
//! initial transform, and a list of [`CurveWeight`] providers. It is used by the multi-curve
//! bundle adjustment solver where each entity needs its own configuration.
//!
//! ## Residual weighting
//!
//! The [`CurveWeight`] trait and its two built-in implementations let callers scale alignment
//! residuals by region:
//!
//! - [`LengthRangeWeight`] - applies a weight to points whose distance along the curve falls
//!   within a given range, leaving all other points at weight `1.0`. This is the 2D stand-in for
//!   `FaceIndexWeight`: a curve is parameterized by arc length, so a contiguous span is the
//!   natural way to name a region of it.
//! - [`NearCurveWeight`] - applies a weight to points that are close to (and similarly oriented
//!   as) a second reference curve, based on distance and normal angle thresholds.
//!
//! # Uncertainty
//!
//! Unlike `Mesh3`, a `Curve2` carries no per-vertex standard deviation of its own, so there is no
//! equivalent of the `Mesh3::point_stdev` fallback that 3D enjoys. An explicit slice supplied to
//! [`AlignmentCurve`] is the only source of target-side uncertainty in 2D, and a curve without one
//! contributes a sigma of zero.

use crate::common::points::dist;
use crate::geom2::curve2::CurveStation2;
use crate::geom2::{Curve2, Iso2, SurfacePoint2};

// ===============================================================================================
// Sample point container
// ===============================================================================================

/// An owned sample point taken from a `Curve2`, holding both where the sample sits on the curve
/// and the surface point itself.
///
/// This is the 2D counterpart of `MeshSurfPoint`. Where that struct locates a sample by face index
/// and barycentric coordinates, this one locates it by edge index and the fraction along that
/// edge, which is the same idea in one dimension fewer.
///
/// `length_along` is redundant with `index` and `fraction` but is kept because arc length is the
/// natural coordinate for a curve, and recovering it later would require the parent curve, which a
/// [`CurveWeight`] does not have. Every field is filled in by [`CurveSurfPoint::from_station`], so
/// the three stay consistent by construction.
#[derive(Debug, Clone, Copy)]
pub struct CurveSurfPoint {
    /// The index of the vertex in the underlying polyline which directly precedes this point.
    pub index: usize,

    /// The fraction of the distance between the vertex at `index` and the next one at which this
    /// point is located, in the range `[0, 1]`.
    pub fraction: f64,

    /// The distance along the curve from its start to this point.
    pub length_along: f64,

    /// The surface point (position and normal) corresponding to this location on the curve.
    pub sp: SurfacePoint2,
}

impl CurveSurfPoint {
    /// Captures a station on a curve as an owned sample point.
    pub fn from_station(station: &CurveStation2) -> Self {
        Self {
            index: station.index(),
            fraction: station.fraction(),
            length_along: station.length_along(),
            sp: station.surface_point(),
        }
    }
}

/// Interpolates a curve's per-vertex standard deviation to the position of a sample, using the
/// edge index and fraction the sample already carries.
///
/// The interpolation is linear in the standard deviation rather than in the variance. Combining
/// the two vertex variances as `sum(w_i^2 * var_i)` would be right if the interpolated point were
/// an *average of independent measurements of the same quantity*, but it isn't: it is a position
/// on a curve which is uncertain by roughly the local noise level everywhere. That form would
/// claim the midpoint of an edge is `1/sqrt(2)` times as uncertain as its endpoints, which is not
/// true of a measured profile. Linear interpolation keeps sigma continuous across edges and equal
/// to the vertex value at each vertex, which is the behavior wanted here.
///
/// # Panics
///
/// Panics if `stdev` is shorter than the curve's vertex count. Callers reaching this through a
/// public alignment entry point are protected by that entry point's argument validation.
pub(crate) fn interpolated_stdev(cp: &CurveSurfPoint, stdev: &[f64]) -> f64 {
    let a = stdev[cp.index];
    let b = stdev[cp.index + 1];
    a + (b - a) * cp.fraction
}

// ===============================================================================================
// Alignment curve container
// ===============================================================================================

/// A container structure that holds all the information necessary to align this curve against a
/// reference. This provides a unified interface for all additional options used to refine the
/// alignment process, such as the uncertainty of the curve vertex points, an initial alignment,
/// and methods of applying weights to the sample points.
pub struct AlignmentCurve<'a> {
    pub curve: &'a Curve2,
    pub uncertainty: Option<&'a [f64]>,
    pub initial: Option<&'a Iso2>,
    pub weights: Option<&'a [Box<dyn CurveWeight + Sync>]>,
}

impl<'a> AlignmentCurve<'a> {
    /// Creates a new `AlignmentCurve` instance.
    ///
    /// # Arguments
    ///
    /// * `curve`: The curve to align.
    /// * `uncertainty`: Optional uncertainty values for the curve vertices, should be in the form
    ///   of a slice of f64 values the same length as the number of vertices in the curve. The
    ///   values should represent standard deviations of distance the vertex would be from the
    ///   current position upon repeated measurements. A normal distribution is assumed for the
    ///   sake of calculating relative probabilities.
    /// * `initial`: Optional initial transformation for the curve. If not specified, the identity
    ///   transformation will be used.
    /// * `weights`: An optional list of weight providing entities that will be used to calculate
    ///   weights of the alignment points _once_ upon initialization. These weights will be
    ///   combined and will then scale the residual calculated at the associated alignment point.
    pub fn new(
        curve: &'a Curve2,
        uncertainty: Option<&'a [f64]>,
        initial: Option<&'a Iso2>,
        weights: Option<&'a [Box<dyn CurveWeight + Sync>]>,
    ) -> Self {
        Self {
            curve,
            uncertainty,
            initial,
            weights,
        }
    }

    /// The starting pose of this curve, or the identity if none was given.
    pub fn transform(&self) -> Iso2 {
        *self.initial.unwrap_or(&Iso2::identity())
    }

    /// The standard deviation of the curve at a sample point, or zero if this curve has no
    /// uncertainty values.
    pub fn sigma_at(&self, cp: &CurveSurfPoint) -> f64 {
        self.uncertainty
            .map_or(0.0, |stdev| interpolated_stdev(cp, stdev))
    }

    /// The combined weight of every [`CurveWeight`] provider at a sample point, or `1.0` if this
    /// curve has none.
    pub fn weight_at(&self, cp: &CurveSurfPoint) -> f64 {
        self.weights
            .map_or(1.0, |ws| ws.iter().map(|w| w.weight(cp)).product())
    }
}

// ===============================================================================================
// Residual weighting
// ===============================================================================================

/// This is a trait for a generic curve weight providing entity. When given a [`CurveSurfPoint`],
/// it should return a weight that will be applied to the residual at that point during the
/// alignment process. This allows for flexible weighting strategies, such as proximity to another
/// curve, or lying within a specific span of the curve, etc.
pub trait CurveWeight {
    /// Returns the weight of the curve point.
    ///
    /// # Arguments
    ///
    /// * `point`: The curve point for which to compute the weight.
    ///
    /// returns: f64
    fn weight(&self, point: &CurveSurfPoint) -> f64;
}

/// Applies a weight to sample points whose distance along the curve falls within a given range.
/// Points outside the range are left at a weight of `1.0`.
#[derive(Clone)]
pub struct LengthRangeWeight {
    start: f64,
    end: f64,
    weight: f64,
}

impl LengthRangeWeight {
    /// Creates a new `LengthRangeWeight` instance.
    ///
    /// # Arguments
    ///
    /// * `start`: The distance along the curve at which the range begins, inclusive.
    /// * `end`: The distance along the curve at which the range ends, inclusive.
    /// * `weight`: The weight to apply to the points within the range.
    pub fn new(start: f64, end: f64, weight: f64) -> Self {
        Self { start, end, weight }
    }

    pub fn to_boxed_trait(self) -> Box<dyn CurveWeight + Sync> {
        Box::new(self)
    }
}

impl CurveWeight for LengthRangeWeight {
    fn weight(&self, point: &CurveSurfPoint) -> f64 {
        if point.length_along >= self.start && point.length_along <= self.end {
            self.weight
        } else {
            1.0
        }
    }
}

/// Applies a weight to sample points which are both close to and similarly oriented as a second
/// reference curve. Points which fail either test are left at a weight of `1.0`.
#[derive(Clone)]
pub struct NearCurveWeight {
    curve: Curve2,
    weight: f64,
    max_dist: f64,
    max_angle: f64,
}

impl NearCurveWeight {
    /// Creates a new `NearCurveWeight` instance.
    ///
    /// # Arguments
    ///
    /// * `curve`: The curve to compute the weight against.
    /// * `weight`: The weight to apply to the curve points.
    /// * `max_dist`: The maximum distance to consider for the weight.
    /// * `max_angle`: The maximum angle between normals to consider for the weight.
    pub fn new(curve: Curve2, weight: f64, max_dist: f64, max_angle: f64) -> Self {
        Self {
            curve,
            weight,
            max_dist,
            max_angle,
        }
    }

    pub fn to_boxed_trait(self) -> Box<dyn CurveWeight + Sync> {
        Box::new(self)
    }
}

impl CurveWeight for NearCurveWeight {
    fn weight(&self, point: &CurveSurfPoint) -> f64 {
        let nearest = self.curve.at_closest_to_point(&point.sp.point);
        let nearest = nearest.surface_point();

        // If the distance is greater than the maximum distance, return 1.0
        if dist(&nearest.point, &point.sp.point) > self.max_dist {
            return 1.0;
        }

        // Check the angle between the normals
        if nearest.normal.angle(&point.sp.normal) > self.max_angle {
            return 1.0;
        }

        self.weight
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{Point2, Vector2};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// A counter-clockwise square with 2-unit sides, so the cumulative vertex lengths are
    /// 0, 2, 4, 6 and (when closed) 8.
    fn square(closed: bool) -> Curve2 {
        let points =
            [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)].map(|(x, y)| Point2::new(x, y));
        Curve2::from_points(&points, 1e-8, closed).unwrap()
    }

    fn sample_at(curve: &Curve2, length: f64) -> CurveSurfPoint {
        CurveSurfPoint::from_station(&curve.at_length(length).unwrap())
    }

    #[test]
    fn a_sample_records_where_it_came_from() {
        let curve = square(false);
        // One quarter of the way along the second edge, which runs from (2,0) to (2,2).
        let cp = sample_at(&curve, 2.5);

        assert_eq!(cp.index, 1);
        assert_relative_eq!(cp.fraction, 0.25, epsilon = 1e-12);
        assert_relative_eq!(cp.length_along, 2.5, epsilon = 1e-12);
        assert_relative_eq!(cp.sp.point, Point2::new(2.0, 0.5), epsilon = 1e-12);
    }

    #[test]
    fn a_sample_carries_the_outward_normal() {
        // The curve is wound counter-clockwise, so the normal on the bottom edge points down.
        let curve = square(true);
        let cp = sample_at(&curve, 1.0);

        assert_relative_eq!(
            cp.sp.normal.into_inner(),
            Vector2::new(0.0, -1.0),
            epsilon = 1e-12
        );
    }

    #[test]
    fn uncertainty_is_interpolated_along_the_edge() {
        let curve = square(false);
        let stdev = [0.1, 0.5, 0.2, 0.4];
        let ac = AlignmentCurve::new(&curve, Some(&stdev), None, None);

        // At a vertex the value is the vertex's own.
        assert_relative_eq!(ac.sigma_at(&sample_at(&curve, 0.0)), 0.1, epsilon = 1e-12);
        assert_relative_eq!(ac.sigma_at(&sample_at(&curve, 2.0)), 0.5, epsilon = 1e-12);

        // Halfway along the first edge it is the average of the two ends.
        assert_relative_eq!(ac.sigma_at(&sample_at(&curve, 1.0)), 0.3, epsilon = 1e-12);

        // A quarter of the way along the second edge, from 0.5 toward 0.2.
        assert_relative_eq!(ac.sigma_at(&sample_at(&curve, 2.5)), 0.425, epsilon = 1e-12);
    }

    #[test]
    fn a_curve_without_uncertainty_contributes_no_sigma() {
        let curve = square(false);
        let ac = AlignmentCurve::new(&curve, None, None, None);

        assert_eq!(ac.sigma_at(&sample_at(&curve, 1.0)), 0.0);
    }

    #[test]
    fn the_transform_defaults_to_the_identity() {
        let curve = square(false);
        let ac = AlignmentCurve::new(&curve, None, None, None);
        assert_relative_eq!(
            ac.transform().to_matrix(),
            Iso2::identity().to_matrix(),
            epsilon = 1e-12
        );

        let initial = Iso2::new(Vector2::new(1.0, 2.0), 0.3);
        let ac = AlignmentCurve::new(&curve, None, Some(&initial), None);
        assert_relative_eq!(
            ac.transform().to_matrix(),
            initial.to_matrix(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn a_length_range_weight_applies_only_inside_its_span() {
        let curve = square(false);
        let w = LengthRangeWeight::new(1.0, 3.0, 0.25);

        assert_relative_eq!(w.weight(&sample_at(&curve, 0.5)), 1.0);
        assert_relative_eq!(w.weight(&sample_at(&curve, 2.0)), 0.25);
        assert_relative_eq!(w.weight(&sample_at(&curve, 5.0)), 1.0);

        // The range is inclusive at both ends.
        assert_relative_eq!(w.weight(&sample_at(&curve, 1.0)), 0.25);
        assert_relative_eq!(w.weight(&sample_at(&curve, 3.0)), 0.25);
    }

    #[test]
    fn weights_from_several_providers_multiply() {
        let curve = square(false);
        let weights: Vec<Box<dyn CurveWeight + Sync>> = vec![
            LengthRangeWeight::new(0.0, 3.0, 0.5).to_boxed_trait(),
            LengthRangeWeight::new(2.0, 6.0, 0.5).to_boxed_trait(),
        ];
        let ac = AlignmentCurve::new(&curve, None, None, Some(&weights));

        // Only the first span.
        assert_relative_eq!(ac.weight_at(&sample_at(&curve, 1.0)), 0.5, epsilon = 1e-12);
        // Both spans overlap here.
        assert_relative_eq!(ac.weight_at(&sample_at(&curve, 2.5)), 0.25, epsilon = 1e-12);
        // Only the second span.
        assert_relative_eq!(ac.weight_at(&sample_at(&curve, 5.0)), 0.5, epsilon = 1e-12);
    }

    #[test]
    fn a_curve_without_weight_providers_is_unweighted() {
        let curve = square(false);
        let ac = AlignmentCurve::new(&curve, None, None, None);
        assert_relative_eq!(ac.weight_at(&sample_at(&curve, 1.0)), 1.0);
    }

    #[test]
    fn a_near_curve_weight_needs_both_proximity_and_agreement() {
        let curve = square(false);

        // A short segment sitting 0.1 below the bottom edge, wound in the same direction so its
        // normal also points down.
        let near = Curve2::from_points(
            &[Point2::new(0.0, -0.1), Point2::new(2.0, -0.1)],
            1e-8,
            false,
        )
        .unwrap();
        let w = NearCurveWeight::new(near, 0.4, 0.5, PI / 4.0);

        // On the bottom edge: close enough and pointing the same way.
        assert_relative_eq!(w.weight(&sample_at(&curve, 1.0)), 0.4, epsilon = 1e-12);

        // On the top edge: far away, so no weight is applied.
        assert_relative_eq!(w.weight(&sample_at(&curve, 5.0)), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn a_near_curve_weight_rejects_a_disagreeing_normal() {
        let curve = square(false);

        // The same nearby segment, but wound the other way, so its normal points up while the
        // bottom edge of the square points down.
        let flipped = Curve2::from_points(
            &[Point2::new(2.0, -0.1), Point2::new(0.0, -0.1)],
            1e-8,
            false,
        )
        .unwrap();
        let w = NearCurveWeight::new(flipped, 0.4, 0.5, PI / 4.0);

        assert_relative_eq!(w.weight(&sample_at(&curve, 1.0)), 1.0, epsilon = 1e-12);
    }
}
