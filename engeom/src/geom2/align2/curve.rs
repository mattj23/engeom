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
//! # What was deliberately left out
//!
//! There is no per-vertex uncertainty and no residual weighting machinery here, even though the
//! 3D mirror (`geom3::align3::mesh`) has both. In 3D they have a producer: `Mesh3::point_stdev`
//! carries real measurement uncertainty from scan data. In 2D nothing produces such data today -
//! a `Curve2` has no stdev attribute, and a mesh section does not interpolate its parent's - so
//! the plumbing was removed rather than carried speculatively. If a producer appears, mirror the
//! machinery back from the 3D side.

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
/// natural coordinate for a curve, and recovering it later would require the parent curve, which
/// consumers of an owned sample do not have. Every field is filled in by
/// [`CurveSurfPoint::from_station`], so they stay consistent by construction.
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

// ===============================================================================================
// Alignment point sampling
// ===============================================================================================

/// Curve-alignment-point parameters, controlling [`generate_alignment_points`].
///
/// This is the 2D counterpart of `GAPParams` ("generate-alignment-points parameters") in
/// `align3`, named differently only so that the two can be glob-imported together. The two flatness criteria that 3D applies to a sampled
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
