//! Aligning points to a point cloud, by way of the tangent plane at each nearest neighbor.

use crate::common::kd_tree::KdTreeSearch;
use crate::geom3::align3::{AlignSurfMatch3, SurfaceTarget3};
use crate::{CloudIndex3, Point3, PointCloud3, Result, UnitVec3, VOXEL_COHERENCE_ATTR};

/// A [`PointCloud3`] presented as something a set of points can be aligned to.
///
/// # Why this is a configured wrapper and `Mesh3` is not
///
/// `Mesh3` implements [`SurfaceTarget3`] directly, because a mesh *is* a surface and the closest
/// point on it is well defined everywhere. A cloud is not a surface, only samples of one, so
/// something has to stand in for the surface between the samples and that stand-in needs a
/// parameter saying how far it may be trusted. Hence a wrapper carrying `max_extrapolation`.
///
/// # The match is a tangent plane, not the nearest point
///
/// Returning the nearest cloud point itself would make the solver's residual a point-to-point
/// distance, which carries a quantization error of roughly half the sample spacing even at a
/// perfect alignment: a test point landing between two cloud points can get no closer than the gap
/// between them. On a cloud reduced to a coarse grid that error is the same order as the
/// misalignment being solved for, so it would not be a refinement but a floor.
///
/// Instead the match point is the query projected onto the tangent plane at the nearest neighbor,
/// which makes the residual a point-to-plane distance and removes the in-plane component of the
/// quantization error. This is the cheap form of what Generalized ICP does with a full covariance
/// per point.
///
/// # A known mismatch with the robust weighting, which this target inherits rather than causes
///
/// `points_to_surface3` sets `RESIDUAL_DOF` to 3, on the grounds that the residual is a Euclidean
/// point-to-target distance in space. That describes the arithmetic rather than the distribution.
/// What sets the degrees of freedom is how many independent dimensions of noise survive into the
/// residual, and against a surface the answer is usually one: noise which moves a test point within
/// the local tangent plane does not change its distance to that surface at all, because the closest
/// point slides along to follow it. Only the normal component contributes.
///
/// **This is not specific to a cloud target.** A mesh match landing inside a triangle is the
/// perpendicular distance to that triangle's plane, which is just as one-dimensional. The residual
/// only regains dimensions where the projection clamps:
///
/// | match lands | `Mesh3` | `CloudTarget3` |
/// |---|---|---|
/// | face interior | 1-D, point to plane | 1-D |
/// | edge | 2-D, point to line | 1-D, plane extended |
/// | vertex, or past the boundary | 3-D, point to point | 1-D, plane extrapolated |
///
/// So a mesh is mostly one-dimensional with a minority of higher-dimensional matches at edges and
/// boundaries, while this target is uniformly one-dimensional because an infinite plane never
/// clamps. That last row is the one worth watching: where a mesh would clamp to a boundary and
/// report a true point-to-point distance, this reports a confident distance to a plane
/// extrapolated into territory nothing was measured in. `max_extrapolation` is what marks those.
///
/// The practical effect of the mismatch is that MAGSAC++ treats residuals as drawn from a wider
/// distribution than they are, so robust refinement down-weights outliers less aggressively than it
/// should. The geometry is unaffected. Fixing it means changing the weighting in
/// `points_to_surface3`, and `MagsacWeight` will not accept a dof of 1 as written, so it is recorded
/// here rather than tuned around. It is a pre-existing property of the solver, and any fix should be
/// judged against the mesh path as much as this one.
pub struct CloudTarget3<'a> {
    index: CloudIndex3<'a>,

    /// Held rather than fetched per query, since `find_align_match` runs once per test point per
    /// solver iteration and these would otherwise be an `Option` unwrap each time.
    normals: &'a [UnitVec3],
    stdev: Option<&'a [f64]>,
    coherence: Option<&'a [f64]>,

    max_extrapolation: f64,
}

impl<'a> CloudTarget3<'a> {
    /// Build an alignment target over a cloud, which must carry per-point normals.
    ///
    /// # Arguments
    ///
    /// * `cloud`: the cloud to align to. It must have normals, because a normal supplies both the
    ///   tangent plane the match lands on and the sign of the residual, and neither can be recovered
    ///   from positions alone. Estimate them with `CloudIndex3::estimate_normals` if the cloud does
    ///   not already carry them.
    /// * `max_extrapolation`: how far *laterally* a query may sit from the nearest cloud point and
    ///   still get a match reported as on-surface. See below, because this is not a distance
    ///   threshold in the way it first appears.
    ///
    /// # What `max_extrapolation` measures, and what it does not
    ///
    /// It bounds the **in-plane** distance from the nearest cloud point, not the total distance to
    /// it. The normal component is deliberately excluded, because that component *is* the residual
    /// the alignment exists to remove: a test point sitting a long way off the surface but directly
    /// above a sample is exactly the case a coarse alignment must be able to see, and gating on
    /// total distance would discard the whole point set on the first iteration of any solve that
    /// started far away.
    ///
    /// What the in-plane distance does say is whether the tangent plane is fiction. Beyond the edge
    /// of the cloud, or across a gap in it, the nearest sample's plane is an extrapolation into
    /// territory nothing was measured in, and a confident match there is worse than no match. Set
    /// this at a small multiple of the cloud's sample spacing.
    ///
    /// Matches beyond the bound are still returned, but with `is_on` false, so they take effect only
    /// when the solve sets `ignore_off_target`.
    pub fn try_new(cloud: &'a PointCloud3, max_extrapolation: f64) -> Result<Self> {
        if !max_extrapolation.is_finite() || max_extrapolation <= 0.0 {
            return Err(format!(
                "max_extrapolation must be finite and positive, got {max_extrapolation}"
            )
            .into());
        }

        let normals = cloud.point_normals().ok_or(
            "A cloud used as an alignment target must carry per-point normals, which supply both \
             the tangent plane and the sign of the residual. See CloudIndex3::estimate_normals.",
        )?;

        Ok(Self {
            index: cloud.compute_index()?,
            normals,
            stdev: cloud.point_stdev(),
            coherence: cloud
                .point_attr(VOXEL_COHERENCE_ATTR)
                .and_then(|a| a.as_scalar()),
            max_extrapolation,
        })
    }

    /// The cloud this target was built over.
    pub fn cloud(&self) -> &'a PointCloud3 {
        self.index.cloud()
    }
}

impl SurfaceTarget3 for CloudTarget3<'_> {
    fn find_align_match(&self, p: &Point3) -> AlignSurfMatch3 {
        let (i, distance) = self.index.nearest_one(p);

        let normal = self.normals[i];
        let nearest = self.index.points()[i];

        // Split the offset from the nearest sample into its normal and in-plane parts. The normal
        // part is the residual the solver is here to remove; the in-plane part is how far the
        // tangent plane is being extrapolated, and is the one that decides whether to trust it.
        let offset = normal.dot(&(p - nearest));
        let lateral = (distance * distance - offset * offset).max(0.0).sqrt();

        let projected = p - normal.into_inner() * offset;

        // A voxel-reduced cloud reports how well each cell's normals agreed. Where they did not,
        // the averaged point is a blend of surfaces facing different ways and its tangent plane
        // means little, so the cloud's own assessment becomes the match weight.
        let weight = self.coherence.map_or(1.0, |c| c[i]);

        let matched =
            AlignSurfMatch3::new(projected, normal, lateral <= self.max_extrapolation, weight);

        match self.stdev {
            Some(stdev) => matched.with_sigma(stdev[i]),
            None => matched,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::align3::{AlignOptions3, AlignParams3, points_to_surface3};
    use crate::{Iso3, Mesh3, Vector3};
    use approx::assert_relative_eq;

    fn engine_blade() -> Mesh3 {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/engine-blade.tcmesh");
        let data = crate::io::read_tc_mesh_file(&path).expect("failed to load engine blade");
        Mesh3::from_data(data, false).expect("failed to accelerate engine blade")
    }

    /// RMS distance the sample is left from where it started, which is the geometric measure of how
    /// well a pose was recovered. A norm over the six parameters would mix radians with millimetres.
    fn rms_displacement(points: &[Point3], t: &Iso3) -> f64 {
        let sum: f64 = points.iter().map(|p| (t * p - p).norm_squared()).sum();
        (sum / points.len() as f64).sqrt()
    }

    #[test]
    fn cloud_target_requires_normals() {
        let bare = PointCloud3::new(vec![Point3::origin(), Point3::new(1.0, 0.0, 0.0)]);
        assert!(CloudTarget3::try_new(&bare, 1.0).is_err());
    }

    #[test]
    fn cloud_target_rejects_a_nonsense_extrapolation() {
        let cloud = engine_blade().sample_dense(1.0, None).unwrap();
        assert!(CloudTarget3::try_new(&cloud, 0.0).is_err());
        assert!(CloudTarget3::try_new(&cloud, f64::NAN).is_err());
    }

    /// The match must be the query projected onto the tangent plane, not the nearest sample. On a
    /// flat sheet of samples a query above a gap should come back at its own footprint, not at the
    /// nearest sample's position.
    #[test]
    fn cloud_target_projects_onto_the_tangent_plane() {
        let mut cloud = PointCloud3::new(vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ]);
        cloud
            .set_point_normals(Some(vec![UnitVec3::new_normalize(Vector3::z()); 4]))
            .unwrap();

        let target = CloudTarget3::try_new(&cloud, 10.0).unwrap();

        // Directly between all four samples, 2 units up.
        let query = Point3::new(0.5, 0.5, 2.0);
        let m = target.find_align_match(&query);

        // Point-to-point would land on a corner and report a distance of sqrt(0.5 + 2^2).
        // Point-to-plane lands directly below the query and reports exactly the height.
        assert_relative_eq!(m.point, Point3::new(0.5, 0.5, 0.0), epsilon = 1e-12);
        assert_relative_eq!((query - m.point).norm(), 2.0, epsilon = 1e-12);
        assert!(m.is_on);
    }

    /// `max_extrapolation` bounds the in-plane distance only. A query far along the normal is the
    /// residual the solver exists to remove and must stay on-surface; a query far to the side is
    /// off the end of the cloud and must not.
    #[test]
    fn cloud_target_gates_on_lateral_distance_not_total_distance() {
        let mut cloud = PointCloud3::new(vec![Point3::origin(), Point3::new(0.5, 0.0, 0.0)]);
        cloud
            .set_point_normals(Some(vec![UnitVec3::new_normalize(Vector3::z()); 2]))
            .unwrap();

        let target = CloudTarget3::try_new(&cloud, 1.0).unwrap();

        // 50 units up, directly above a sample. Enormous total distance, no lateral distance.
        let above = target.find_align_match(&Point3::new(0.0, 0.0, 50.0));
        assert!(above.is_on, "a large normal offset must not be gated out");

        // 5 units to the side, on the surface. Small total distance, all of it lateral.
        let beside = target.find_align_match(&Point3::new(5.0, 0.0, 0.0));
        assert!(!beside.is_on, "extrapolating past the cloud must be gated");
    }

    /// A reduced cloud's coherence becomes the match weight, so voxels which straddled an edge
    /// speak for themselves.
    #[test]
    fn cloud_target_takes_its_weight_from_voxel_coherence() {
        let mesh = Mesh3::create_box(20.0, 20.0, 20.0, false);
        let cloud = mesh
            .sample_dense(0.5, None)
            .unwrap()
            .reduce_by_voxel(2.0)
            .unwrap();

        let target = CloudTarget3::try_new(&cloud, 4.0).unwrap();

        // Every weight should be the cloud's own coherence value, so at least one match on an edge
        // voxel must come back weighted below one.
        let weights: Vec<f64> = cloud
            .points()
            .iter()
            .map(|p| target.find_align_match(p).weight)
            .collect();

        let lowest = weights.iter().copied().fold(f64::INFINITY, f64::min);
        let highest = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        assert!(highest > 0.99, "a flat face should be weighted near one");
        assert!(
            lowest < 0.9,
            "an edge voxel should be weighted down, got {lowest}"
        );
    }

    /// Returns the nearest cloud point itself, which is what `CloudTarget3` deliberately does not
    /// do. Used only to measure what that decision is worth.
    struct NearestPointTarget<'a> {
        index: CloudIndex3<'a>,
        normals: &'a [UnitVec3],
    }

    impl SurfaceTarget3 for NearestPointTarget<'_> {
        fn find_align_match(&self, p: &Point3) -> AlignSurfMatch3 {
            let (i, _) = self.index.nearest_one(p);
            AlignSurfMatch3::new(self.index.points()[i], self.normals[i], true, 1.0)
        }
    }

    /// The tangent-plane projection is the whole reason this type is not three lines long, so the
    /// alternative is measured rather than argued about.
    ///
    /// A point-to-point residual cannot fall below the sample spacing even at a perfect pose,
    /// because a test point landing between samples can get no closer than the gap between them.
    /// A point-to-plane residual has no such floor on a smooth surface.
    #[test]
    fn the_tangent_plane_projection_is_what_removes_the_sampling_floor() {
        let mesh = engine_blade();
        let test_points = mesh.sample_poisson(1.0, None).unwrap().points().to_vec();

        let center = mesh.aabb().center();
        let delta = Iso3::new(Vector3::new(1.5, -0.8, 0.6), Vector3::new(0.0, 0.0, 0.04));
        let moved: Vec<Point3> = test_points.iter().map(|p| delta * p).collect();

        let spacing = 2.0;
        let cloud = mesh.sample_dense(spacing, None).unwrap();
        let normals = cloud.point_normals().unwrap();

        let opts = AlignOptions3::default();
        let recover = |out: &crate::geom3::AlignOutcome3| {
            rms_displacement(&test_points, &(out.alignment().full_transform() * delta))
        };

        let plane = CloudTarget3::try_new(&cloud, spacing * 3.0).unwrap();
        let plane_err = recover(
            &points_to_surface3(
                &moved,
                &plane,
                AlignParams3::from_center(center, None),
                &opts,
            )
            .expect("tangent plane solve failed"),
        );

        let nearest = NearestPointTarget {
            index: cloud.compute_index().unwrap(),
            normals,
        };
        let nearest_err = recover(
            &points_to_surface3(
                &moved,
                &nearest,
                AlignParams3::from_center(center, None),
                &opts,
            )
            .expect("nearest point solve failed"),
        );

        println!();
        println!("at {spacing} mm sample spacing:");
        println!("  tangent plane:  {plane_err:.3e} mm");
        println!("  nearest point:  {nearest_err:.3e} mm");
        println!(
            "  factor:         {:.0}x",
            nearest_err / plane_err.max(1e-15)
        );

        assert!(
            plane_err < nearest_err / 100.0,
            "the tangent plane projection should be far better than a nearest point match, \
             got {plane_err} against {nearest_err}"
        );
    }

    /// The headline check: recover a known pose against clouds of several sample spacings, and
    /// report each against a mesh target solving the identical problem.
    ///
    /// The gap is the cost of standing samples in for a surface. It should grow with the spacing,
    /// because the tangent plane at the nearest sample is a worse stand-in for the real surface the
    /// further apart the samples are, and that trend is the thing worth knowing before using a
    /// reduced cloud as an alignment target. The numbers are printed rather than asserted tightly,
    /// because pinning them would turn a measurement into a tripwire.
    #[test]
    fn cloud_target_recovers_a_pose_and_reports_its_cost_against_a_mesh() {
        let mesh = engine_blade();

        let test_points = mesh.sample_poisson(1.0, None).unwrap().points().to_vec();
        assert!(test_points.len() > 200);

        let center = mesh.aabb().center();
        let delta = Iso3::new(Vector3::new(1.5, -0.8, 0.6), Vector3::new(0.0, 0.0, 0.04));
        let moved: Vec<Point3> = test_points.iter().map(|p| delta * p).collect();
        let start = rms_displacement(&test_points, &delta);

        let opts = AlignOptions3::default();
        let recover = |out: &crate::geom3::AlignOutcome3| {
            rms_displacement(&test_points, &(out.alignment().full_transform() * delta))
        };

        let mesh_out = points_to_surface3(
            &moved,
            &mesh,
            AlignParams3::from_center(center, None),
            &opts,
        )
        .expect("mesh solve failed");
        let mesh_err = recover(&mesh_out);

        println!();
        println!("start displacement: {start:.4} mm");
        println!("mesh target:        {mesh_err:.3e} mm");
        println!();
        println!(
            "{:>10}  {:>10}  {:>12}  {:>10}",
            "spacing", "points", "RMS mm", "vs mesh"
        );

        let mut previous = 0.0;
        for spacing in [0.25, 1.0, 2.0] {
            let cloud = mesh.sample_dense(spacing, None).unwrap();
            let target = CloudTarget3::try_new(&cloud, spacing * 3.0).unwrap();

            let out = points_to_surface3(
                &moved,
                &target,
                AlignParams3::from_center(center, None),
                &opts,
            )
            .expect("cloud solve failed");
            let err = recover(&out);

            println!(
                "{:>10.2}  {:>10}  {:>12.3e}  {:>9.1}x",
                spacing,
                cloud.point_count(),
                err,
                err / mesh_err.max(1e-15)
            );

            assert!(
                err < start / 10.0,
                "spacing {spacing} did not recover the pose: {err} from {start}"
            );

            assert!(
                err >= previous * 0.25,
                "error fell sharply as samples got sparser, which is backwards"
            );
            previous = err;
        }
    }
}
