//! An alignment target built over an unordered set of measured points rather than a curve.
//!
//! This is the 2D counterpart of `geom3::align3::cloud`, and is deliberately structured similarly
//! to it. The one structural difference is the input: 3D has a `PointCloud3` which carries its own
//! normals, standard deviations, and voxel coherence, while 2D has no such container, so the parts
//! are passed in as parallel slices.
//!
//! # Why the match is not simply the nearest point
//!
//! Taking the nearest sample outright would leave a residual floor that no amount of solving can
//! remove: a test point landing between two samples can get no closer than the gap between them.
//! On a set reduced to a coarse spacing that error is the same order as the misalignment being
//! solved for, so it would not be a refinement but a floor.
//!
//! Instead the match point is the query projected onto the tangent line at the nearest neighbor,
//! which makes the residual a point-to-line distance and removes the along-surface component of
//! the quantization error. This is the cheap form of what Generalized ICP does with a full
//! covariance per point.
//!
//! # A known mismatch with the robust weighting, which this target inherits rather than causes
//!
//! `points_to_surface2` sets its residual degrees of freedom to 2, on the grounds that the
//! residual is a Euclidean point-to-target distance in the plane. That describes the arithmetic
//! rather than the distribution. What sets the degrees of freedom is how many independent
//! dimensions of noise survive into the residual, and against a surface the answer is usually one:
//! noise which slides a test point along the local tangent does not change its distance to that
//! surface at all, because the closest point slides along to follow it. Only the normal component
//! contributes.
//!
//! **This is not specific to this target.** A curve match landing inside a segment is the
//! perpendicular distance to that segment's line, which is just as one-dimensional. The residual
//! only regains a dimension where the projection clamps:
//!
//! | match lands | `Curve2` | `CloudTarget2` |
//! |---|---|---|
//! | segment interior | 1-D, point to line | 1-D |
//! | vertex, or past an open end | 2-D, point to point | 1-D, line extrapolated |
//!
//! So a curve is mostly one-dimensional with a minority of two-dimensional matches at its vertices
//! and ends, while this target is uniformly one-dimensional because an infinite line never clamps.
//! That last row is the one worth watching: where a curve would clamp to an end and report a true
//! point-to-point distance, this reports a confident distance to a line extrapolated into
//! territory nothing was measured in. `max_extrapolation` is what marks those.
//!
//! The practical effect of the mismatch is that MAGSAC++ treats residuals as drawn from a wider
//! distribution than they are, so robust refinement down-weights outliers less aggressively than
//! it should. The geometry is unaffected. Fixing it means changing the weighting in
//! `points_to_surface2`, and `MagsacWeight` will not accept a dof of 1 as written, so it is
//! recorded here rather than tuned around. It is a pre-existing property of the solver, and any
//! fix should be judged against the curve path as much as this one.

use crate::common::kd_tree::KdTreeSearch;
use crate::geom2::align2::{AlignSurfMatch2, SurfaceTarget2};
use crate::geom2::{KdTree2, Point2};
use crate::{Result, UnitVec2};

/// An alignment target over an unordered set of measured points, each carrying a normal.
///
/// See the module documentation for why the match is a projection onto the tangent line rather
/// than the nearest point itself.
pub struct CloudTarget2<'a> {
    index: KdTree2,

    points: &'a [Point2],

    /// Held rather than fetched per query, since `find_align_match` runs once per test point per
    /// solver iteration.
    normals: &'a [UnitVec2],
    stdev: Option<&'a [f64]>,

    max_extrapolation: f64,
}

impl<'a> CloudTarget2<'a> {
    /// Build an alignment target over a set of points which carry per-point normals.
    ///
    /// # Arguments
    ///
    /// * `points`: the measured positions to align to.
    /// * `normals`: one normal per point, in the same order. A normal supplies both the tangent
    ///   line the match lands on and the sign of the residual, and neither can be recovered from
    ///   positions alone.
    /// * `stdev`: an optional measurement standard deviation per point, in the same order, which
    ///   becomes the target-side uncertainty of a match.
    /// * `max_extrapolation`: how far *along the surface* a query may sit from the nearest sample
    ///   and still get a match reported as on-surface. See below, because this is not a distance
    ///   threshold in the way it first appears.
    ///
    /// # What `max_extrapolation` measures, and what it does not
    ///
    /// It bounds the **along-surface** distance from the nearest sample, not the total distance to
    /// it. The normal component is deliberately excluded, because that component *is* the residual
    /// the alignment exists to remove: a test point sitting a long way off the surface but directly
    /// above a sample is exactly the case a coarse alignment must be able to see, and gating on
    /// total distance would discard the whole point set on the first iteration of any solve that
    /// started far away.
    ///
    /// What the along-surface distance does say is whether the tangent line is fiction. Beyond the
    /// end of the point set, or across a gap in it, the nearest sample's line is an extrapolation
    /// into territory nothing was measured in, and a confident match there is worse than no match.
    /// Set this at a small multiple of the set's sample spacing.
    ///
    /// Matches beyond the bound are still returned, but with `is_on` false, so they take effect
    /// only when the solve sets `ignore_off_target`.
    ///
    /// # Errors
    ///
    /// Returns an error if `max_extrapolation` is not finite and positive, if the point set is
    /// empty, or if the normals or standard deviations disagree in length with the points.
    pub fn try_new(
        points: &'a [Point2],
        normals: &'a [UnitVec2],
        stdev: Option<&'a [f64]>,
        max_extrapolation: f64,
    ) -> Result<Self> {
        if !max_extrapolation.is_finite() || max_extrapolation <= 0.0 {
            return Err(format!(
                "max_extrapolation must be finite and positive, got {max_extrapolation}"
            )
            .into());
        }
        if points.is_empty() {
            return Err("a point set alignment target must have at least one point".into());
        }
        if normals.len() != points.len() {
            return Err(format!(
                "there are {} points but {} normals",
                points.len(),
                normals.len()
            )
            .into());
        }
        if let Some(s) = stdev
            && s.len() != points.len()
        {
            return Err(format!(
                "there are {} points but {} standard deviations",
                points.len(),
                s.len()
            )
            .into());
        }

        Ok(Self {
            index: KdTree2::try_new(points)?,
            points,
            normals,
            stdev,
            max_extrapolation,
        })
    }

    /// The points this target was built over.
    pub fn points(&self) -> &'a [Point2] {
        self.points
    }
}

impl SurfaceTarget2 for CloudTarget2<'_> {
    fn find_align_match(&self, p: &Point2) -> AlignSurfMatch2 {
        let (i, distance) = self.index.nearest_one(p);

        let normal = self.normals[i];
        let nearest = self.points[i];

        // Split the offset from the nearest sample into its normal and along-surface parts. The
        // normal part is the residual the solver is here to remove; the along-surface part is how
        // far the tangent line is being extrapolated, and is the one that decides whether to
        // trust it.
        let offset = normal.dot(&(p - nearest));
        let lateral = (distance * distance - offset * offset).max(0.0).sqrt();

        let projected = p - normal.into_inner() * offset;

        let matched =
            AlignSurfMatch2::new(projected, normal, lateral <= self.max_extrapolation, 1.0);

        match self.stdev {
            Some(stdev) => matched.with_sigma(stdev[i]),
            None => matched,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::dist;
    use crate::geom2::align2::{AlignOptions2, AlignParams2, points_to_surface2};
    use crate::geom2::{Curve2, Iso2, Vector2};
    use approx::assert_relative_eq;

    /// A horizontal run of samples along `y = 0` at unit spacing, all facing +y.
    fn flat_run() -> (Vec<Point2>, Vec<UnitVec2>) {
        let points: Vec<Point2> = (-5..=5).map(|i| Point2::new(i as f64, 0.0)).collect();
        let normals = vec![UnitVec2::new_unchecked(Vector2::y()); points.len()];
        (points, normals)
    }

    #[test]
    fn a_nonsense_extrapolation_bound_is_rejected() {
        let (points, normals) = flat_run();
        assert!(CloudTarget2::try_new(&points, &normals, None, 0.0).is_err());
        assert!(CloudTarget2::try_new(&points, &normals, None, -1.0).is_err());
        assert!(CloudTarget2::try_new(&points, &normals, None, f64::NAN).is_err());
    }

    #[test]
    fn mismatched_or_empty_inputs_are_rejected() {
        let (points, normals) = flat_run();

        assert!(CloudTarget2::try_new(&[], &[], None, 1.0).is_err());
        assert!(CloudTarget2::try_new(&points, &normals[..3], None, 1.0).is_err());

        let short = vec![0.1; 3];
        assert!(CloudTarget2::try_new(&points, &normals, Some(&short), 1.0).is_err());
    }

    #[test]
    fn the_match_is_a_projection_onto_the_tangent_line() {
        // A query above the midpoint between two samples. Taking the nearest sample outright
        // would report a match half a unit to the side; the tangent line projection puts it
        // directly below the query, which is what removes the sampling floor.
        let (points, normals) = flat_run();
        let target = CloudTarget2::try_new(&points, &normals, None, 1.0).unwrap();

        let query = Point2::new(0.5, 2.0);
        let m = target.find_align_match(&query);

        assert_relative_eq!(m.point, Point2::new(0.5, 0.0), epsilon = 1e-12);
        assert_relative_eq!(dist(&query, &m.point), 2.0, epsilon = 1e-12);
        assert!(m.is_on);
    }

    #[test]
    fn the_gate_is_on_along_surface_distance_not_total_distance() {
        let (points, normals) = flat_run();
        let target = CloudTarget2::try_new(&points, &normals, None, 1.0).unwrap();

        // Far away along the normal, but directly above a sample. This is the case a coarse
        // alignment has to be able to see, so it must not be gated out.
        let above = target.find_align_match(&Point2::new(0.0, 50.0));
        assert!(above.is_on, "a large normal offset must not be gated out");

        // Just off the end of the run, where the tangent line is pure extrapolation.
        let beside = target.find_align_match(&Point2::new(9.0, 0.0));
        assert!(!beside.is_on, "extrapolating past the set must be gated");
    }

    #[test]
    fn the_target_reports_its_own_uncertainty() {
        let (points, normals) = flat_run();
        let stdev: Vec<f64> = (0..points.len()).map(|i| 0.01 * (i + 1) as f64).collect();
        let target = CloudTarget2::try_new(&points, &normals, Some(&stdev), 1.0).unwrap();

        // Nearest to the sample at x = -5, which is index 0.
        let m = target.find_align_match(&Point2::new(-5.0, 1.0));
        assert_relative_eq!(m.sigma, 0.01, epsilon = 1e-12);

        // ...and to the sample at x = 5, the last one.
        let m = target.find_align_match(&Point2::new(5.0, 1.0));
        assert_relative_eq!(m.sigma, 0.01 * points.len() as f64, epsilon = 1e-12);
    }

    #[test]
    fn the_tangent_line_projection_is_what_removes_the_sampling_floor() {
        // The claim the whole design rests on. A coarsely sampled point set is aligned against,
        // and the projection has to do better than a target which reports the nearest sample
        // outright, because that one cannot resolve below its own spacing.
        //
        // The sampling is deliberately irregular. With a uniformly sampled symmetric shape the
        // quantization error of a nearest-sample match cancels around the shape and the pose comes
        // out fine anyway, which would make this test prove nothing. Uneven spacing is both the
        // harder case and the one measured data actually presents.
        struct NearestSample<'a> {
            points: &'a [Point2],
            normals: &'a [UnitVec2],
            index: KdTree2,
        }

        impl SurfaceTarget2 for NearestSample<'_> {
            fn find_align_match(&self, p: &Point2) -> AlignSurfMatch2 {
                let (i, _) = self.index.nearest_one(p);
                AlignSurfMatch2::new(self.points[i], self.normals[i], true, 1.0)
            }
        }

        // A coarsely and irregularly sampled ellipse, which is what measured data looks like:
        // the gaps between samples are uneven, so the quantization error of a nearest-sample
        // match does not cancel out around the shape the way a uniform sampling would let it.
        let mut rg = crate::common::random_geometry::RandomGeometry2::from_seed(0x5a3_f100);
        let on_ellipse = |a: f64| Point2::new(8.0 * a.cos(), 3.0 * a.sin());
        let normal_at = |a: f64| {
            // Outward normal of the ellipse at the given parameter.
            UnitVec2::new_normalize(Vector2::new(3.0 * a.cos(), 8.0 * a.sin()))
        };

        let n = 60;
        let angles: Vec<f64> = (0..n)
            .map(|i| {
                let base = std::f64::consts::TAU * (i as f64) / (n as f64);
                base + rg.f64_sym(0.4) * std::f64::consts::TAU / (n as f64)
            })
            .collect();
        let sampled: Vec<Point2> = angles.iter().map(|&a| on_ellipse(a)).collect();
        let normals: Vec<UnitVec2> = angles.iter().map(|&a| normal_at(a)).collect();

        // Test points taken densely and at an unrelated phase, so the true answer is the identity
        // but no test point sits on top of a sample.
        let dense: Vec<Point2> = (0..500)
            .map(|i| on_ellipse(std::f64::consts::TAU * (i as f64 + 0.37) / 500.0))
            .collect();

        let disturb = Iso2::new(Vector2::new(0.08, -0.06), 0.01);
        let moved: Vec<Point2> = dense.iter().map(|p| disturb * p).collect();

        let opts = AlignOptions2::default();
        let params = || AlignParams2::from_origin(None);

        let projected = CloudTarget2::try_new(&sampled, &normals, None, 2.0).unwrap();
        let nearest = NearestSample {
            points: &sampled,
            normals: &normals,
            index: KdTree2::try_new(&sampled).unwrap(),
        };

        // How far the test points are left from where they belong, which is the quantity that
        // matters rather than any single component of the transform.
        let rms = |t: &Iso2| -> f64 {
            let residual = t * disturb;
            (dense
                .iter()
                .map(|p| dist(p, &(residual * p)).powi(2))
                .sum::<f64>()
                / dense.len() as f64)
                .sqrt()
        };

        let with_projection = points_to_surface2(&moved, &projected, params(), &opts).unwrap();
        let with_nearest = points_to_surface2(&moved, &nearest, params(), &opts).unwrap();

        let e_projection = rms(with_projection.alignment().full_transform());
        let e_nearest = rms(with_nearest.alignment().full_transform());

        assert!(
            e_projection < 0.1 * e_nearest,
            "the tangent line projection should resolve well below the sample spacing: \
             projected {e_projection}, nearest sample {e_nearest}"
        );
    }

    #[test]
    fn a_point_set_target_recovers_a_disturbed_pose() {
        // End to end through the single-body solver, against a point set sampled off a curve.
        let curve = Curve2::from_points(
            &[
                Point2::new(-5.0, -3.0),
                Point2::new(5.0, -3.0),
                Point2::new(5.0, 3.0),
                Point2::new(-5.0, 3.0),
            ],
            1e-8,
            true,
        )
        .unwrap();

        let spacing = 0.25;
        let count = (curve.length() / spacing).floor() as usize;
        let stations: Vec<_> = (0..count)
            .filter_map(|k| curve.at_length((k as f64 + 0.5) * spacing))
            .collect();
        let sampled: Vec<Point2> = stations.iter().map(|s| s.point()).collect();
        let normals: Vec<UnitVec2> = stations.iter().map(|s| s.normal()).collect();

        let target = CloudTarget2::try_new(&sampled, &normals, None, 1.0).unwrap();

        let disturb = Iso2::new(Vector2::new(0.15, -0.1), 0.02);
        let moved: Vec<Point2> = sampled.iter().map(|p| disturb * p).collect();

        let outcome = points_to_surface2(
            &moved,
            &target,
            AlignParams2::from_origin(None),
            &AlignOptions2::default(),
        )
        .unwrap();

        // The recovered transform should undo the disturbance.
        let composed = outcome.alignment().full_transform() * disturb;
        assert_relative_eq!(
            composed.to_matrix(),
            Iso2::identity().to_matrix(),
            epsilon = 1e-4
        );
    }
}
