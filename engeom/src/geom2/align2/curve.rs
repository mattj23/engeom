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
//! ## Alignment body container
//!
//! [`AlignmentGroup`] wraps a `&CurveGroup2` together with optional per-vertex uncertainty
//! values, an initial transform, and a list of [`CurveWeight`] providers. It is used by the
//! multi-curve bundle adjustment solver where each rigid body needs its own configuration. A
//! body is always a group; a lone curve rides along as a group of one via `CurveGroup2::from`.
//!
//! ## Residual weighting
//!
//! The [`CurveWeight`] trait and its two built-in implementations let callers scale alignment
//! residuals by region:
//!
//! - [`LengthRangeWeight`] - applies a weight to points on one named member whose distance
//!   along that member falls within a given range, leaving all other points at weight `1.0`.
//!   This is the 2D stand-in for `FaceIndexWeight`: a member curve is parameterized by arc
//!   length, so a member plus a span of it is the natural way to name a region of a body.
//! - [`NearCurveWeight`] - applies a weight to points that are close to (and similarly oriented
//!   as) a second reference curve, based on distance and normal angle thresholds.
//!
//! # Uncertainty
//!
//! Unlike `Mesh3`, a curve carries no per-vertex standard deviation of its own, so there is no
//! equivalent of the `Mesh3::point_stdev` fallback that 3D enjoys. An explicit slice supplied to
//! [`AlignmentGroup`] is the only source of target-side uncertainty in 2D, and a body without
//! one contributes a sigma of zero.
//!
//! The slice covers the whole group, with the members' vertices concatenated in member order;
//! `CurveGroup2::vertex_offset` locates each member's span. Note that a closed member stores its
//! seam vertex twice, so it counts one more value than it has corners.

use crate::common::points::dist;
use crate::common::vec_f64::mean_and_stdev;
use crate::geom2::align2::target::is_curve_endpoint;
use crate::geom2::curve2::CurveStation2;
use crate::geom2::{Curve2, CurveGroup2, Iso2, SurfacePoint2};
use std::f64::consts::PI;

// ===============================================================================================
// Sample point container
// ===============================================================================================

/// An owned sample point taken from a curve within a [`crate::geom2::CurveGroup2`], holding
/// which member it sits on, where on that member it sits, and the surface point itself.
///
/// This is the 2D counterpart of `MeshSurfPoint`. Where that struct locates a sample by face index
/// and barycentric coordinates, this one locates it by edge index and the fraction along that
/// edge, which is the same idea in one dimension fewer. What has no 3D counterpart is `member`: a
/// multi-component mesh is still one `Mesh3` with one face index space, while a group of curves
/// keeps its members separate, so the member index is part of the sample's address.
///
/// `index`, `fraction`, and `length_along` are all local to the member curve, not to the group.
///
/// `length_along` is redundant with `index` and `fraction` but is kept because arc length is the
/// natural coordinate for a curve, and recovering it later would require the parent curve, which a
/// [`CurveWeight`] does not have. Every field is filled in by [`CurveSurfPoint::from_station`], so
/// they stay consistent by construction.
#[derive(Debug, Clone, Copy)]
pub struct CurveSurfPoint {
    /// The index of the member curve this point belongs to, within its owning group. A sample
    /// taken from a lone curve uses `0`, matching the group of one that `CurveGroup2::from`
    /// builds.
    pub member: usize,

    /// The index of the vertex in the member's underlying polyline which directly precedes this
    /// point.
    pub index: usize,

    /// The fraction of the distance between the vertex at `index` and the next one at which this
    /// point is located, in the range `[0, 1]`.
    pub fraction: f64,

    /// The distance along the member curve from its start to this point.
    pub length_along: f64,

    /// The surface point (position and normal) corresponding to this location on the curve.
    pub sp: SurfacePoint2,
}

impl CurveSurfPoint {
    /// Captures a station on a curve as an owned sample point, recording which member of its
    /// group the curve is.
    pub fn from_station(station: &CurveStation2, member: usize) -> Self {
        Self {
            member,
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
// Alignment body container
// ===============================================================================================

/// A container structure that holds all the information necessary to align this body against a
/// reference. This provides a unified interface for all additional options used to refine the
/// alignment process, such as the uncertainty of the curve vertex points, an initial alignment,
/// and methods of applying weights to the sample points.
///
/// A body is a [`CurveGroup2`]: one or more disjoint curves that move together rigidly, such as
/// the loops and open segments a planar section of a mesh produces. To align a lone curve, wrap
/// it in a group of one with `CurveGroup2::from`.
pub struct AlignmentGroup<'a> {
    pub group: &'a CurveGroup2,
    pub uncertainty: Option<&'a [f64]>,
    pub initial: Option<&'a Iso2>,
    pub weights: Option<&'a [Box<dyn CurveWeight + Sync>]>,
}

impl<'a> AlignmentGroup<'a> {
    /// Creates a new `AlignmentGroup` instance.
    ///
    /// # Arguments
    ///
    /// * `group`: The group of curves to align as one rigid body.
    /// * `uncertainty`: Optional uncertainty values for the curve vertices, should be in the form
    ///   of a slice of f64 values whose length is the group's total vertex count, with the
    ///   members' vertices concatenated in member order. The values should represent standard
    ///   deviations of distance the vertex would be from the current position upon repeated
    ///   measurements. A normal distribution is assumed for the sake of calculating relative
    ///   probabilities.
    /// * `initial`: Optional initial transformation for the body. If not specified, the identity
    ///   transformation will be used.
    /// * `weights`: An optional list of weight providing entities that will be used to calculate
    ///   weights of the alignment points _once_ upon initialization. These weights will be
    ///   combined and will then scale the residual calculated at the associated alignment point.
    pub fn new(
        group: &'a CurveGroup2,
        uncertainty: Option<&'a [f64]>,
        initial: Option<&'a Iso2>,
        weights: Option<&'a [Box<dyn CurveWeight + Sync>]>,
    ) -> Self {
        Self {
            group,
            uncertainty,
            initial,
            weights,
        }
    }

    /// The starting pose of this body, or the identity if none was given.
    pub fn transform(&self) -> Iso2 {
        *self.initial.unwrap_or(&Iso2::identity())
    }

    /// The standard deviation of the body at a sample point, interpolated within the owning
    /// member's span of the concatenated uncertainty values, or zero if this body has none.
    pub fn sigma_at(&self, cp: &CurveSurfPoint) -> f64 {
        self.uncertainty.map_or(0.0, |stdev| {
            let o = self.group.vertex_offset(cp.member);
            let n = self.group.curves()[cp.member].count();
            interpolated_stdev(cp, &stdev[o..o + n])
        })
    }

    /// The combined weight of every [`CurveWeight`] provider at a sample point, or `1.0` if this
    /// body has none.
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

/// Applies a weight to sample points on one named member whose distance along that member falls
/// within a given range. Points on other members, or outside the range, are left at a weight of
/// `1.0`.
///
/// The member index is part of the address because `length_along` is a member-local coordinate:
/// a span of arc length means nothing without saying which member it is measured along.
#[derive(Clone)]
pub struct LengthRangeWeight {
    member: usize,
    start: f64,
    end: f64,
    weight: f64,
}

impl LengthRangeWeight {
    /// Creates a new `LengthRangeWeight` instance.
    ///
    /// # Arguments
    ///
    /// * `member`: The index of the member curve the range lies on.
    /// * `start`: The distance along the member at which the range begins, inclusive.
    /// * `end`: The distance along the member at which the range ends, inclusive.
    /// * `weight`: The weight to apply to the points within the range.
    pub fn new(member: usize, start: f64, end: f64, weight: f64) -> Self {
        Self {
            member,
            start,
            end,
            weight,
        }
    }

    pub fn to_boxed_trait(self) -> Box<dyn CurveWeight + Sync> {
        Box::new(self)
    }
}

impl CurveWeight for LengthRangeWeight {
    fn weight(&self, point: &CurveSurfPoint) -> f64 {
        if point.member == self.member
            && point.length_along >= self.start
            && point.length_along <= self.end
        {
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

// ===============================================================================================
// Alignment point sampling
// ===============================================================================================

/// Parameters controlling [`generate_alignment_points`].
///
/// This is the 2D counterpart of `GAPParams` in `align3`, named differently only so that the two
/// can be glob-imported together. The two flatness criteria that 3D applies to a sampled
/// neighborhood (`out_of_plane_ratio` and `centroid_ratio`) have no 2D counterpart: a curve is
/// already an ordered one-dimensional path, so there is no neighborhood to fit a basis to and no
/// mesh edge to fall off. What replaces them is [`CAPParams::max_corner_angle`], which rejects
/// samples sitting on a sharp turn, and [`CAPParams::end_margin`], which keeps samples away from
/// the ends of an open curve.
#[derive(Debug, Clone, Copy)]
pub struct CAPParams {
    /// The arc length spacing between samples taken along the test curve. It also sets the
    /// distance at which the turn of the curve is examined on either side of a sample.
    pub sample_spacing: f64,

    /// The maximum permissible angle, in radians, between the normal at a sample and the normals
    /// one sample spacing to either side of it, on both the test and the reference curve.
    ///
    /// This is what rejects samples on or beside a corner, where the closest point on the other
    /// curve jumps rather than slides as the bodies move. A value of `PI / 3.0` (60 degrees) is a
    /// reasonable starting point; leave it large for a noisy curve.
    pub max_corner_angle: f64,

    /// The distance from either end of an open curve within which samples are discarded, in the
    /// units of the geometry.
    ///
    /// A projection past the end of an open curve clamps to the endpoint, which produces a
    /// confident-looking correspondence to a place the curve does not actually go. Both the test
    /// and the reference side are kept clear of this. A closed curve has no ends and ignores it.
    pub end_margin: f64,

    /// An optional value that, if provided, will result in a filtering operation of the final
    /// selected candidates based on the distance between them and the reference curve. A value of
    /// `Some(3.0)` will filter out candidates that are more than 3 standard deviations above the
    /// mean candidate distance to the reference curve. This can be used to remove outlying areas
    /// of the test curve as the two curves begin to converge towards alignment.
    pub filter_distances: Option<f64>,
}

impl CAPParams {
    /// Creates a new set of parameters for the curve sampling algorithm.
    ///
    /// # Arguments
    ///
    /// * `sample_spacing`: The arc length spacing between samples on the test curve. Smaller
    ///   values mean a denser sampling, but also more points to check and more influence from
    ///   smaller features in the curve.
    /// * `max_corner_angle`: The maximum permissible angle between the normal at a sample and the
    ///   normals one sample spacing to either side of it.
    /// * `end_margin`: The distance from either end of an open curve within which samples are
    ///   discarded.
    /// * `filter_distances`: An optional number of standard deviations above the mean
    ///   candidate-to-reference distance beyond which candidates are discarded.
    ///
    /// returns: CAPParams
    pub fn new(
        sample_spacing: f64,
        max_corner_angle: f64,
        end_margin: f64,
        filter_distances: Option<f64>,
    ) -> Self {
        Self {
            sample_spacing,
            max_corner_angle,
            end_margin,
            filter_distances,
        }
    }

    /// Creates a new set of default parameters for the curve sampling algorithm, requiring only
    /// the sample spacing to be specified.
    ///
    /// # Arguments
    ///
    /// * `sample_spacing`: The arc length spacing between samples on the test curve.
    ///
    /// returns: CAPParams
    pub fn defaults(sample_spacing: f64) -> Self {
        Self {
            sample_spacing,
            max_corner_angle: PI / 3.0,
            end_margin: 2.0 * sample_spacing,
            filter_distances: Some(3.0),
        }
    }
}

/// A sampling algorithm that finds a set of good alignment points on a test body which can be
/// used to align it with a reference body. Pay close attention to the parameters.
///
/// The method walks every member of the test group at a uniform arc length spacing and keeps the
/// samples which look like they will behave during a solve. Each check runs against the member
/// curve that owns the position being checked, so a group mixing closed loops and open segments
/// is treated at every position like the single curve that position belongs to. A sample is
/// retained when:
///
/// 1. It is at least `end_margin` from either end of its own member, if that member is open.
/// 2. The normals one sample spacing to either side of it, along its own member, are within
///    `max_corner_angle` of its own.
/// 3. After being moved by `iso`, its closest point on any member of the reference group is not
///    a clamped endpoint of an open member.
/// 4. The normal at that closest point faces the same way as the moved sample's normal.
/// 5. The owning reference member satisfies the same corner and end-margin checks at the closest
///    point.
///
/// At the very end, if a sigma value is provided in `filter_distances`, the mean and standard
/// deviation of the distance from each candidate to its corresponding projected point are
/// computed, and any candidates further than `sigma` standard deviations above the mean are
/// discarded. The filter is pooled across the whole body rather than run per member, because it
/// measures body-to-body convergence.
///
/// # Arguments
///
/// * `test`: The body which will be sampled for alignment points: the resulting points will be
///   on this group's members, in the group's own coordinates.
/// * `reference`: The body which is used as a reference for the alignment. The resulting points
///   will be good points to use when aligning the test body to this reference body.
/// * `iso`: An initial transform that will be applied to the test body's points before
///   projecting them onto the reference. This would represent an initial guess of the alignment
///   that is to follow, and takes the test body's coordinates into the reference's.
/// * `params`: The parameters for the sampling algorithm.
///
/// returns: Vec<CurveSurfPoint, Global>
pub fn generate_alignment_points(
    test: &CurveGroup2,
    reference: &CurveGroup2,
    iso: &Iso2,
    params: &CAPParams,
) -> Vec<CurveSurfPoint> {
    let mut candidates = Vec::new();

    for (member, test_curve) in test.curves().iter().enumerate() {
        for cp in sample_curve(test_curve, params.sample_spacing, member) {
            if !clear_of_the_ends(test_curve, cp.length_along, params.end_margin) {
                continue;
            }
            if !turns_gently(test_curve, &cp, params) {
                continue;
            }

            // Move the sample into the reference body's frame and find where it lands, on
            // whichever member is nearest.
            let moved = cp.sp.transformed_by(iso);
            let (ref_m, station) = reference.at_closest_to_point(&moved.point);
            let ref_curve = &reference.curves()[ref_m];

            // A projection past the end of an open member clamps to the endpoint, which is not a
            // real correspondence.
            if !ref_curve.is_closed() && is_curve_endpoint(&station, ref_curve) {
                continue;
            }

            let ref_cp = CurveSurfPoint::from_station(&station, ref_m);
            if ref_cp.sp.normal.dot(&moved.normal) < 0.0 {
                continue;
            }
            if !clear_of_the_ends(ref_curve, ref_cp.length_along, params.end_margin) {
                continue;
            }
            if !turns_gently(ref_curve, &ref_cp, params) {
                continue;
            }

            candidates.push((dist(&moved.point, &ref_cp.sp.point), cp));
        }
    }

    // Lastly, filter out candidates with distances more than a specified number of standard
    // deviations above the mean distance.
    //
    // A zero spread is left alone rather than filtered. With no variation there is no outlier to
    // find, and the threshold would land exactly on the mean and so discard every candidate. Two
    // coincident curves are the case that reaches this, which is rare in measured data but is
    // precisely what a synthetic test looks like.
    if let Some(sigma) = params.filter_distances
        && candidates.len() > 10
    {
        let distances = candidates.iter().map(|(d, _)| *d).collect::<Vec<_>>();
        if let Some((mean, st_dev)) = mean_and_stdev(&distances)
            && st_dev > 0.0
        {
            let threshold = mean + sigma * st_dev;
            candidates.retain(|(d, _)| *d < threshold);
        }
    }

    candidates.into_iter().map(|(_, cp)| cp).collect()
}

/// Walks a member curve at a uniform arc length spacing, offset by half a spacing so that
/// samples do not land on the vertices of a regularly shaped polyline, where the normal is
/// ambiguous. The samples are stamped with the member index they belong to.
fn sample_curve(curve: &Curve2, spacing: f64, member: usize) -> Vec<CurveSurfPoint> {
    if !spacing.is_finite() || spacing <= 0.0 || curve.length() <= 0.0 {
        return Vec::new();
    }

    let count = (curve.length() / spacing).floor() as usize;
    (0..count)
        .filter_map(|k| {
            let l = (k as f64 + 0.5) * spacing;
            curve
                .at_length(l)
                .map(|s| CurveSurfPoint::from_station(&s, member))
        })
        .collect()
}

/// Whether a position along a curve is far enough from either end. A closed curve has no ends.
fn clear_of_the_ends(curve: &Curve2, length_along: f64, margin: f64) -> bool {
    curve.is_closed() || (length_along >= margin && length_along <= curve.length() - margin)
}

/// Whether the curve turns by no more than `max_corner_angle` within one sample spacing on either
/// side of a point. This is what keeps samples off corners, where the closest point on the other
/// curve jumps rather than slides as the bodies move.
fn turns_gently(curve: &Curve2, cp: &CurveSurfPoint, params: &CAPParams) -> bool {
    for delta in [-params.sample_spacing, params.sample_spacing] {
        let raw = cp.length_along + delta;
        let l = if curve.is_closed() {
            raw.rem_euclid(curve.length())
        } else {
            raw
        };

        // Off the end of an open curve, which `clear_of_the_ends` is responsible for.
        let Some(station) = curve.at_length(l) else {
            continue;
        };

        if station.normal().angle(&cp.sp.normal) > params.max_corner_angle {
            return false;
        }
    }
    true
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
        CurveSurfPoint::from_station(&curve.at_length(length).unwrap(), 0)
    }

    /// Wraps a fixture curve as a body of one for the container tests.
    fn group_of(curve: &Curve2) -> CurveGroup2 {
        CurveGroup2::from(curve.clone())
    }

    #[test]
    fn a_sample_records_where_it_came_from() {
        let curve = square(false);
        // One quarter of the way along the second edge, which runs from (2,0) to (2,2).
        let cp = sample_at(&curve, 2.5);

        assert_eq!(cp.member, 0);
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
        let group = group_of(&curve);
        let stdev = [0.1, 0.5, 0.2, 0.4];
        let ac = AlignmentGroup::new(&group, Some(&stdev), None, None);

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
        let group = group_of(&curve);
        let ac = AlignmentGroup::new(&group, None, None, None);

        assert_eq!(ac.sigma_at(&sample_at(&curve, 1.0)), 0.0);
    }

    #[test]
    fn the_transform_defaults_to_the_identity() {
        let curve = square(false);
        let group = group_of(&curve);
        let ac = AlignmentGroup::new(&group, None, None, None);
        assert_relative_eq!(
            ac.transform().to_matrix(),
            Iso2::identity().to_matrix(),
            epsilon = 1e-12
        );

        let initial = Iso2::new(Vector2::new(1.0, 2.0), 0.3);
        let ac = AlignmentGroup::new(&group, None, Some(&initial), None);
        assert_relative_eq!(
            ac.transform().to_matrix(),
            initial.to_matrix(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn a_length_range_weight_applies_only_inside_its_span() {
        let curve = square(false);
        let w = LengthRangeWeight::new(0, 1.0, 3.0, 0.25);

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
            LengthRangeWeight::new(0, 0.0, 3.0, 0.5).to_boxed_trait(),
            LengthRangeWeight::new(0, 2.0, 6.0, 0.5).to_boxed_trait(),
        ];
        let group = group_of(&curve);
        let ac = AlignmentGroup::new(&group, None, None, Some(&weights));

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
        let group = group_of(&curve);
        let ac = AlignmentGroup::new(&group, None, None, None);
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

    // ============================================================================================
    // Alignment point sampling
    // ============================================================================================

    /// A closed counter-clockwise rectangle 10 by 6, centered on the origin. Its corners sit at
    /// arc lengths 0, 10, 16 and 26 out of a perimeter of 32.
    fn rect_curve() -> Curve2 {
        let points =
            [(-5.0, -3.0), (5.0, -3.0), (5.0, 3.0), (-5.0, 3.0)].map(|(x, y)| Point2::new(x, y));
        Curve2::from_points(&points, 1e-8, true).unwrap()
    }

    /// A straight open line along x, so it has ends but no corners at all.
    fn open_line() -> Curve2 {
        Curve2::from_points(
            &[Point2::new(0.0, 0.0), Point2::new(10.0, 0.0)],
            1e-8,
            false,
        )
        .unwrap()
    }

    #[test]
    fn samples_are_spaced_along_the_curve_and_miss_the_vertices() {
        let curve = rect_curve();
        let samples = sample_curve(&curve, 0.5, 0);

        // A perimeter of 32 at a spacing of 0.5 gives 64 samples, offset by half a spacing so
        // that none of them land on a corner.
        assert_eq!(samples.len(), 64);
        assert_relative_eq!(samples[0].length_along, 0.25, epsilon = 1e-12);
        assert_relative_eq!(samples[1].length_along, 0.75, epsilon = 1e-12);
        for s in &samples {
            for corner in [0.0, 10.0, 16.0, 26.0, 32.0] {
                assert!((s.length_along - corner).abs() > 1e-9);
            }
        }
    }

    #[test]
    fn a_degenerate_spacing_samples_nothing() {
        let curve = rect_curve();
        assert!(sample_curve(&curve, 0.0, 0).is_empty());
        assert!(sample_curve(&curve, -1.0, 0).is_empty());
        assert!(sample_curve(&curve, f64::NAN, 0).is_empty());
    }

    #[test]
    fn the_corner_check_rejects_samples_beside_a_corner() {
        let curve = rect_curve();
        let params = CAPParams::defaults(0.5);

        // Well inside the bottom edge, with both neighbors on the same edge.
        assert!(turns_gently(&curve, &sample_at(&curve, 5.0), &params));

        // A quarter spacing short of the corner at 10.0, so the neighbor one spacing ahead has
        // already turned the 90 degree corner.
        assert!(!turns_gently(&curve, &sample_at(&curve, 9.75), &params));
        assert!(!turns_gently(&curve, &sample_at(&curve, 10.25), &params));

        // A closed curve wraps, so the seam at 0.0 is treated like any other corner.
        assert!(!turns_gently(&curve, &sample_at(&curve, 0.25), &params));
        assert!(!turns_gently(&curve, &sample_at(&curve, 31.75), &params));
    }

    #[test]
    fn the_end_margin_only_applies_to_open_curves() {
        let closed = rect_curve();
        let open = open_line();

        // A closed curve has no ends, so nothing is ever near one.
        assert!(clear_of_the_ends(&closed, 0.0, 2.0));
        assert!(clear_of_the_ends(&closed, closed.length(), 2.0));

        assert!(!clear_of_the_ends(&open, 0.5, 2.0));
        assert!(!clear_of_the_ends(&open, 9.5, 2.0));
        assert!(clear_of_the_ends(&open, 5.0, 2.0));
        // The margin itself is inclusive at both ends.
        assert!(clear_of_the_ends(&open, 2.0, 2.0));
        assert!(clear_of_the_ends(&open, 8.0, 2.0));
    }

    #[test]
    fn generated_points_avoid_the_corners_of_a_closed_curve() {
        let test = rect_curve();
        let reference = rect_curve();
        let params = CAPParams::defaults(0.5);

        let points = generate_alignment_points(
            &group_of(&test),
            &group_of(&reference),
            &Iso2::identity(),
            &params,
        );

        // Two of the 64 samples beside each of the four corners are rejected on the test curve,
        // and the reference curve rejects the same ones, leaving 56.
        assert_eq!(points.len(), 56);
        for p in &points {
            for corner in [0.0, 10.0, 16.0, 26.0, 32.0] {
                assert!(
                    (p.length_along - corner).abs() > 0.5,
                    "a sample at {} is too close to the corner at {corner}",
                    p.length_along
                );
            }
        }
    }

    #[test]
    fn coincident_curves_still_produce_correspondences() {
        // Every candidate distance is zero here, so the mean is zero and the spread is zero. If
        // the distance filter were applied anyway its threshold would land on the mean and throw
        // away every candidate, leaving the solve with nothing to work with.
        let test = rect_curve();
        let reference = rect_curve();
        let params = CAPParams::defaults(0.5);
        assert!(params.filter_distances.is_some());

        let points = generate_alignment_points(
            &group_of(&test),
            &group_of(&reference),
            &Iso2::identity(),
            &params,
        );
        assert!(!points.is_empty());
    }

    #[test]
    fn the_distance_filter_drops_the_furthest_candidates() {
        // Two straight lines fanning apart, so the distance from the test curve to the reference
        // grows steadily along its length and the far end is what the filter should remove.
        // Straight lines keep the corner check out of the way entirely.
        let test = open_line();
        let reference = Curve2::from_points(
            &[Point2::new(0.0, 0.0), Point2::new(10.0, 2.0)],
            1e-8,
            false,
        )
        .unwrap();

        let unfiltered = CAPParams::new(0.25, PI / 3.0, 0.5, None);
        let filtered = CAPParams::new(0.25, PI / 3.0, 0.5, Some(1.0));

        let test = group_of(&test);
        let reference = group_of(&reference);
        let a = generate_alignment_points(&test, &reference, &Iso2::identity(), &unfiltered);
        let b = generate_alignment_points(&test, &reference, &Iso2::identity(), &filtered);

        assert!(
            b.len() < a.len(),
            "the filter should have dropped the far samples: {} unfiltered, {} filtered",
            a.len(),
            b.len()
        );

        // What survives is the near end of the curve, not an arbitrary subset of it.
        let worst_kept = b.iter().map(|p| p.length_along).fold(f64::MIN, f64::max);
        let worst_dropped = a.iter().map(|p| p.length_along).fold(f64::MIN, f64::max);
        assert!(worst_kept < worst_dropped);
    }

    #[test]
    fn samples_whose_normals_disagree_are_rejected() {
        // The same line twice, but the reference is wound the other way so its normal points the
        // opposite direction everywhere. Nothing should survive.
        let test = open_line();
        let reference = Curve2::from_points(
            &[Point2::new(10.0, 0.0), Point2::new(0.0, 0.0)],
            1e-8,
            false,
        )
        .unwrap();

        let params = CAPParams::new(0.5, PI / 3.0, 1.0, None);
        let points = generate_alignment_points(
            &group_of(&test),
            &group_of(&reference),
            &Iso2::identity(),
            &params,
        );

        assert!(points.is_empty());
    }

    #[test]
    fn samples_projecting_past_an_open_end_are_rejected() {
        // The reference covers only the first half of the test curve, so samples on the second
        // half clamp to its endpoint and must be discarded rather than matched there.
        let test = open_line();
        let reference =
            Curve2::from_points(&[Point2::new(0.0, 0.0), Point2::new(5.0, 0.0)], 1e-8, false)
                .unwrap();

        let params = CAPParams::new(0.5, PI / 3.0, 0.5, None);
        let points = generate_alignment_points(
            &group_of(&test),
            &group_of(&reference),
            &Iso2::identity(),
            &params,
        );

        assert!(!points.is_empty());
        for p in &points {
            assert!(
                p.length_along < 5.0,
                "a sample at {} matched past the end of the reference",
                p.length_along
            );
        }
    }

    #[test]
    fn samples_carry_their_member_and_respect_each_members_ends() {
        // A test body of two parallel open lines just off a common reference line. Samples must
        // come from both members, stamped with the right index, and each member's end margin is
        // measured along that member.
        let above = Curve2::from_points(
            &[Point2::new(0.0, 0.1), Point2::new(10.0, 0.1)],
            1e-8,
            false,
        )
        .unwrap();
        let below = Curve2::from_points(
            &[Point2::new(0.0, -0.1), Point2::new(10.0, -0.1)],
            1e-8,
            false,
        )
        .unwrap();
        let test = CurveGroup2::new(vec![above, below]).unwrap();
        let reference = group_of(&open_line());

        let params = CAPParams::new(0.5, PI / 3.0, 0.5, None);
        let points = generate_alignment_points(&test, &reference, &Iso2::identity(), &params);

        // Each 10-unit member yields 20 half-offset samples, of which the end margin removes the
        // first and last, leaving 18 per member.
        let from_0 = points.iter().filter(|p| p.member == 0).count();
        let from_1 = points.iter().filter(|p| p.member == 1).count();
        assert_eq!(from_0, 18);
        assert_eq!(from_1, 18);
        for p in &points {
            assert!(p.length_along >= 0.5 && p.length_along <= 9.5);
        }
    }

    #[test]
    fn reference_checks_run_against_the_owning_member() {
        // The reference body mixes a short open segment far away (member 0) with a closed
        // rectangle (member 1), and the test body coincides with the rectangle. Every match
        // lands on member 1, and the checks must be judged against it: if the endpoint clamp or
        // end margin were wrongly evaluated against the short open member 0, nearly every sample
        // would be rejected, because their arc lengths run far past its end.
        let far_segment = Curve2::from_points(
            &[Point2::new(0.0, 50.0), Point2::new(1.0, 50.0)],
            1e-8,
            false,
        )
        .unwrap();
        let reference = CurveGroup2::new(vec![far_segment, rect_curve()]).unwrap();
        let test = group_of(&rect_curve());

        let params = CAPParams::defaults(0.5);
        let points = generate_alignment_points(&test, &reference, &Iso2::identity(), &params);

        // The same count the single-curve corner test establishes: 64 samples minus two beside
        // each of the four corners.
        assert_eq!(points.len(), 56);
    }
}
