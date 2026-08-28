//! Tools for asking how well a set of points constrains a 3D alignment, and for choosing a
//! subset of them that preserves that conditioning.
//!
//! The analysis itself is dimension-independent and lives in
//! [`crate::common::align::information::AlignInfo`], whose module documentation explains
//! the information matrix and the two uses it serves:
//!
//! - **Diagnosis.** Is this alignment well conditioned, and if not, along what motion is it free
//!   to slide? See [`AlignInfo3::marginal_precision`] and
//!   [`AlignInfo3::weak_directions`]. In 3D the weak direction of an alignment is rarely
//!   axis-aligned: points on a cylinder slide along a helix, points on a plane slide in any
//!   in-plane direction and spin about the normal, and a shallow dish is weakest in a coupled
//!   tilt-and-shift.
//! - **Pruning.** Which points are actually carrying the alignment, and what is the smallest
//!   subset that would constrain it just as well? See [`AlignInfo3::leverage`] and
//!   [`AlignInfo3::select_d_optimal`]. Pruning is the reason this exists: a simultaneous
//!   alignment of fifteen or twenty scans is expensive in proportion to the number of alignment
//!   points, and most of those points are redundant with their neighbors.
//!
//! This module supplies the 3D ingredients: the [`AlignDofs`] implementation on [`Dof6`] and the
//! [`AlignInfo3::from_points`] constructor, which turns points and a [`SurfaceTarget3`]
//! into jacobian rows. The `geom2::align2::information` module does the same for 2D.
//!
//! # A caveat about units
//!
//! The translation columns of `H` are dimensionless (a dot product of two unit vectors) while the
//! rotation columns carry units of length, because a rotation partial scales with the distance
//! from the local origin. `H` is therefore heterogeneous, and two consequences follow.
//!
//! [`AlignInfo3::marginal_precision`] is safe: each entry concerns a single parameter, so
//! its units are internally consistent. [`AlignInfo3::weak_directions`] mixes translation
//! and rotation in one vector and is only meaningful when the two are comparably scaled, which in
//! practice means the local origin should sit near the middle of the point set (as
//! `AlignParams3::from_center` arranges). Placing the origin far away inflates the rotation block
//! and will make rotations look artificially well constrained.

use crate::Point3;
use crate::common::align::information::{AlignDofs, AlignInfo};
use crate::geom3::align3::jacobian::point_surf_jacobian;
use crate::geom3::align3::{AlignParams3, Dof6, SurfaceTarget3};

/// The information content of a set of test points with respect to a 3D alignment, over the six
/// parameters `tx`, `ty`, `tz`, `rx`, `ry`, `rz`.
///
/// Built with [`AlignInfo3::from_points`], which projects the points onto a target and
/// takes the jacobian row of each correspondence, or with [`AlignInfo3::from_rows`] when
/// the rows have already been computed. See the module documentation for what it is good for.
pub type AlignInfo3 = AlignInfo<Dof6, 6>;

impl AlignDofs<6> for Dof6 {
    fn free_indices(&self) -> Vec<usize> {
        [self.tx, self.ty, self.tz, self.rx, self.ry, self.rz]
            .iter()
            .enumerate()
            .filter_map(|(i, &free)| free.then_some(i))
            .collect()
    }
}

impl AlignInfo3 {
    /// Builds the information content of a set of points against a target, in the pose described
    /// by `params`.
    ///
    /// Each point is moved by the current transform, projected onto the target, and turned into a
    /// jacobian row exactly as [`crate::geom3::align3::points_to_surface3`] would do. A point is
    /// weighted by the target's own `weight` for the correspondence, divided by the variance the
    /// target reports there, so that a match on a noisier part of the target carries
    /// proportionally less information.
    ///
    /// This deliberately ignores robust weighting. The question being asked is a geometric one
    /// ("how well does this arrangement of points pin down the pose"), and the answer should not
    /// depend on which points happen to be outliers in some particular fit.
    ///
    /// # Arguments
    ///
    /// * `points`: the test points, in their own local coordinate system
    /// * `target`: the stationary entity the points would be aligned to
    /// * `params`: the alignment parameterization, which fixes the local origin the rotation
    ///   partials are taken about and which degrees of freedom are free
    pub fn from_points(
        points: &[Point3],
        target: &impl SurfaceTarget3,
        params: &AlignParams3,
    ) -> Self {
        let current = params.compute_values();
        let mut rows = Vec::with_capacity(points.len());
        let mut weights = Vec::with_capacity(points.len());

        for p in points {
            let m = current.transform * p;
            let c = target.find_align_match(&m);

            rows.push(point_surf_jacobian(&m, &c, &current));
            weights.push(if c.sigma > 0.0 {
                c.weight / (c.sigma * c.sigma)
            } else {
                c.weight
            });
        }

        Self::assemble(rows, weights, params.dof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::align3::{AlignStorage3, AlignSurfMatch3};
    use crate::na::Matrix6;
    use crate::{Iso3, Mesh3, UnitVec3, Vector3};
    use approx::assert_relative_eq;

    /// A flat target in the z = 0 plane, so every match reports the same upward normal. Its
    /// information content is known in closed form, which is what makes it useful here.
    struct GroundPlane;

    impl SurfaceTarget3 for GroundPlane {
        fn find_align_match(&self, p: &Point3) -> AlignSurfMatch3 {
            AlignSurfMatch3::new(
                Point3::new(p.x, p.y, 0.0),
                UnitVec3::new_unchecked(Vector3::z()),
                true,
                1.0,
            )
        }
    }

    /// Points scattered over a patch of the z = 1 plane, one unit above `GroundPlane`.
    fn plane_points() -> Vec<Point3> {
        let mut points = Vec::new();
        for i in -3..=3 {
            for j in -3..=3 {
                points.push(Point3::new(i as f64, j as f64, 1.0));
            }
        }
        points
    }

    fn box_mesh() -> Mesh3 {
        Mesh3::create_box(10.0, 5.0, 2.0, false)
    }

    fn box_points(mesh: &Mesh3) -> Vec<Point3> {
        mesh.sample_poisson(0.5, None)
            .iter()
            .map(|p| p.point())
            .collect()
    }

    // ============================================================================================
    // Conditioning against geometry whose answer is known in advance
    // ============================================================================================

    #[test]
    fn a_plane_leaves_three_motions_unconstrained() {
        // Points on a plane can only ever resist motion along its normal. Sliding in x or y, or
        // spinning about z, moves every point along the surface and changes no residual, so
        // exactly three of the six degrees of freedom are free.
        let info = AlignInfo3::from_points(
            &plane_points(),
            &GroundPlane,
            &AlignParams3::from_origin(None),
        );

        let weak = info.weak_directions();
        assert_eq!(weak.len(), 6);

        // Three eigenvalues at zero, three well away from it.
        let unconstrained: Vec<f64> = weak.iter().take(3).map(|(v, _)| *v).collect();
        for v in &unconstrained {
            assert!(v.abs() < 1e-9, "expected an unconstrained motion, got {v}");
        }
        assert!(
            weak[3].0 > 1.0,
            "expected the remaining motions to be constrained, got {}",
            weak[3].0
        );

        // ...and those three motions live entirely in the tx / ty / rz subspace. Any individual
        // eigenvector may mix them, so the test is that they have no tz / rx / ry component.
        for (_, direction) in weak.iter().take(3) {
            for i in [2usize, 3, 4] {
                assert!(
                    direction[i].abs() < 1e-9,
                    "unconstrained direction leaked into parameter {i}: {direction:?}"
                );
            }
        }
    }

    #[test]
    fn a_singular_problem_reports_rather_than_lying() {
        let info = AlignInfo3::from_points(
            &plane_points(),
            &GroundPlane,
            &AlignParams3::from_origin(None),
        );

        let err = info.marginal_precision().unwrap_err().to_string();
        assert!(err.contains("singular"), "unexpected message: {err}");
        assert!(info.leverage().is_err());
    }

    #[test]
    fn locking_the_free_motions_makes_a_plane_well_conditioned() {
        // The same points, but with the three motions they cannot see locked out. What remains is
        // fully determined, so the analysis should now succeed.
        let dof = Dof6::new(false, false, true, true, true, false);
        let params = AlignParams3::from_origin(Some(dof));
        let info = AlignInfo3::from_points(&plane_points(), &GroundPlane, &params);

        assert_eq!(info.free_dof_count(), 3);

        let precision = info.marginal_precision().unwrap();
        assert_eq!(precision[0], None, "tx was locked");
        assert_eq!(precision[1], None, "ty was locked");
        assert_eq!(precision[5], None, "rz was locked");
        for i in [2usize, 3, 4] {
            assert!(
                precision[i].unwrap() > 0.0,
                "parameter {i} should be constrained"
            );
        }
    }

    #[test]
    fn a_box_constrains_everything() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let params = AlignParams3::from_origin(None);
        let info = AlignInfo3::from_points(&points, &mesh, &params);

        assert_eq!(info.free_dof_count(), 6);
        let precision = info.marginal_precision().unwrap();
        for (i, p) in precision.iter().enumerate() {
            assert!(
                p.unwrap() > 0.0,
                "parameter {i} should be constrained on a closed box"
            );
        }
    }

    // ============================================================================================
    // Leverage
    // ============================================================================================

    #[test]
    fn leverage_sums_to_the_free_dof_count() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let info = AlignInfo3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 6.0, epsilon = 1e-9);
    }

    #[test]
    fn leverage_sums_to_the_free_dof_count_when_constrained() {
        // The identity `sum(h_i) == k` has to hold on the free subspace too, not just at k = 6.
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let dof = Dof6::new(false, true, true, true, false, true);
        let info = AlignInfo3::from_points(&points, &mesh, &AlignParams3::from_origin(Some(dof)));

        assert_eq!(info.free_dof_count(), 4);
        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 4.0, epsilon = 1e-9);
    }

    // ============================================================================================
    // Greedy selection
    // ============================================================================================

    #[test]
    fn selection_prefers_isolated_points_over_a_redundant_cluster() {
        // This is the property that separates greedy selection from independent scoring. A tight
        // cluster of near-identical points all contribute the same row, so after the first one is
        // taken the rest add nothing; the isolated points elsewhere on the box each bring
        // something new and should be chosen first.
        let mesh = box_mesh();

        let mut points = Vec::new();

        // Twenty near-duplicates in one spot on the top face.
        for i in 0..20 {
            points.push(Point3::new(1.0 + i as f64 * 1e-4, 1.0, 1.0));
        }
        let cluster = 0..20;

        // Six isolated points spread over other faces.
        let isolated = [
            Point3::new(-4.0, -2.0, 1.0),
            Point3::new(4.0, 2.0, -1.0),
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(-5.0, 1.0, 0.5),
            Point3::new(0.0, 2.5, -0.5),
            Point3::new(2.0, -2.5, 0.5),
        ];
        points.extend_from_slice(&isolated);

        let info = AlignInfo3::from_points(&points, &mesh, &AlignParams3::from_origin(None));
        let picked = info.select_d_optimal(7, None);

        assert_eq!(picked.len(), 7);

        // Every isolated point should be taken before a second member of the cluster is.
        let from_cluster = picked.iter().filter(|&&i| cluster.contains(&i)).count();
        assert_eq!(
            from_cluster, 1,
            "expected exactly one cluster member among {picked:?}"
        );
    }

    #[test]
    fn selection_is_a_prefix_so_it_can_be_truncated() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let info = AlignInfo3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        let long = info.select_d_optimal(20, None);
        let short = info.select_d_optimal(8, None);

        assert_eq!(&long[..8], &short[..]);
    }

    #[test]
    fn a_selected_subset_stays_well_conditioned() {
        // The point of pruning: a small chosen subset should hold the alignment up nearly as well
        // as the whole set, and far better than the same number of arbitrary points.
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let params = AlignParams3::from_origin(None);
        let info = AlignInfo3::from_points(&points, &mesh, &params);

        let n = 40;
        let picked = info.select_d_optimal(n, None);
        let chosen: Vec<Point3> = picked.iter().map(|&i| points[i]).collect();
        let arbitrary: Vec<Point3> = points.iter().take(n).copied().collect();

        // The comparison is on the smallest eigenvalue, which measures how close the subset comes
        // to leaving some motion unconstrained. It is used rather than `marginal_precision`
        // because an arbitrary subset can be outright singular (the first `n` Poisson samples
        // tend to land on a single face of the box), and there is no finite precision to compare
        // in that case.
        let weakest = |pts: &[Point3]| -> f64 {
            AlignInfo3::from_points(pts, &mesh, &params).weak_directions()[0].0
        };

        let chosen_weakest = weakest(&chosen);
        let arbitrary_weakest = weakest(&arbitrary);

        assert!(
            chosen_weakest > arbitrary_weakest * 10.0,
            "the selected subset should be far better conditioned than an arbitrary one: \
             selected {chosen_weakest}, arbitrary {arbitrary_weakest}"
        );

        // ...and unlike the arbitrary subset, it should constrain the pose outright.
        assert!(
            AlignInfo3::from_points(&chosen, &mesh, &params)
                .marginal_precision()
                .is_ok(),
            "the selected subset should leave no motion unconstrained"
        );
    }

    #[test]
    fn selection_respects_the_available_count() {
        let mesh = box_mesh();
        let points = box_points(&mesh);
        let info = AlignInfo3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        assert!(info.select_d_optimal(0, None).is_empty());
        assert_eq!(
            info.select_d_optimal(points.len() + 50, None).len(),
            points.len()
        );
    }

    #[test]
    fn zero_weight_points_are_never_selected() {
        let rows = vec![
            AlignStorage3::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            AlignStorage3::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            AlignStorage3::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
        ];
        let weights = vec![1.0, 0.0, 1.0];
        let dof = Dof6::new(true, true, true, false, false, false);

        let info = AlignInfo3::from_rows(rows, weights, dof).unwrap();
        let picked = info.select_d_optimal(3, None);

        assert!(!picked.contains(&1), "a zero-weight point was selected");
        assert_eq!(picked.len(), 2);
    }

    // ============================================================================================
    // Construction and validation
    // ============================================================================================

    #[test]
    fn information_matrix_is_expanded_with_zeros_for_locked_dof() {
        let rows = vec![AlignStorage3::new(1.0, 2.0, 3.0, 0.0, 0.0, 0.0)];
        let dof = Dof6::new(true, true, true, false, false, false);
        let info = AlignInfo3::from_rows(rows, vec![1.0], dof).unwrap();

        let h = info.information();
        // The free block is the outer product of (1, 2, 3) with itself.
        assert_relative_eq!(h[(0, 0)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(h[(0, 1)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(h[(2, 2)], 9.0, epsilon = 1e-12);
        // ...and every locked row and column is zero.
        for i in 3..6 {
            for j in 0..6 {
                assert_eq!(h[(i, j)], 0.0);
                assert_eq!(h[(j, i)], 0.0);
            }
        }
    }

    #[test]
    fn mismatched_rows_and_weights_are_rejected() {
        let rows = vec![AlignStorage3::zeros(); 3];
        let err = AlignInfo3::from_rows(rows, vec![1.0; 2], Dof6::all())
            .unwrap_err()
            .to_string();
        assert!(err.contains("weights"), "unexpected message: {err}");
    }

    #[test]
    fn invalid_weights_are_rejected() {
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let rows = vec![AlignStorage3::zeros(); 2];
            let result = AlignInfo3::from_rows(rows, vec![1.0, bad], Dof6::all());
            assert!(result.is_err(), "weight {bad} should have been rejected");
        }
    }

    #[test]
    fn a_fully_locked_problem_degrades_gracefully() {
        let dof = Dof6::new(false, false, false, false, false, false);
        let info =
            AlignInfo3::from_rows(vec![AlignStorage3::zeros(); 4], vec![1.0; 4], dof).unwrap();

        assert_eq!(info.free_dof_count(), 0);
        assert!(info.weak_directions().is_empty());
        assert!(info.select_d_optimal(4, None).is_empty());
        assert_eq!(info.information(), Matrix6::zeros());
    }

    #[test]
    fn the_local_origin_is_respected() {
        // Rotation partials are taken about the local origin, so moving it changes the rotation
        // block of the information matrix. This guards against the origin being quietly ignored.
        let mesh = box_mesh();
        let points = box_points(&mesh);

        let at_origin = AlignInfo3::from_points(&points, &mesh, &AlignParams3::from_origin(None));
        let offset = AlignInfo3::from_points(
            &points,
            &mesh,
            &AlignParams3::from_center(Point3::new(50.0, 0.0, 0.0), None),
        );

        let delta = (at_origin.information() - offset.information()).norm();
        assert!(
            delta > 1.0,
            "the local origin had no effect (delta {delta})"
        );
    }

    #[test]
    fn a_transformed_pose_is_analyzed_where_it_currently_sits() {
        // The analysis describes the pose described by `params`, not the pose the points were
        // authored in, so a working offset has to be honored.
        let mesh = box_mesh();
        let points = box_points(&mesh);

        let here = AlignInfo3::from_points(&points, &mesh, &AlignParams3::from_origin(None));

        let displaced = AlignParams3::new(
            crate::geom3::align3::AlignOrigin3::Origin,
            Some(Iso3::translation(0.0, 0.0, 20.0)),
            None,
        );
        let there = AlignInfo3::from_points(&points, &mesh, &displaced);

        let delta = (here.information() - there.information()).norm();
        assert!(
            delta > 1e-6,
            "the working offset had no effect (delta {delta})"
        );
    }
}
