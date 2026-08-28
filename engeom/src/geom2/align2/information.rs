//! Tools for asking how well a set of points constrains a 2D alignment, and for choosing a
//! subset of them that preserves that conditioning.
//!
//! The analysis itself is dimension-independent and lives in
//! [`crate::common::align::information::AlignInfo`], whose module documentation explains
//! the information matrix and the two uses it serves:
//!
//! - **Diagnosis.** Is this alignment well conditioned, and if not, along what motion is it free
//!   to slide? See [`AlignInfo2::marginal_precision`] and
//!   [`AlignInfo2::weak_directions`]. In 2D the weak direction of an alignment is rarely
//!   axis-aligned: points on a circle can spin about its center, points on a straight line slide
//!   along it, and a shallow arc is weakest in a coupled turn-and-shift.
//! - **Pruning.** Which points are actually carrying the alignment, and what is the smallest
//!   subset that would constrain it just as well? See [`AlignInfo2::leverage`] and
//!   [`AlignInfo2::select_d_optimal`].
//!
//! This module supplies the 2D ingredients: the [`AlignDofs`] implementation on [`Dof3`] and the
//! [`AlignInfo2::from_points`] constructor, which turns points and a [`SurfaceTarget2`]
//! into jacobian rows. The `geom3::align3::information` module does the same for 3D.
//!
//! # A caveat about units
//!
//! The translation columns of `H` are dimensionless (a dot product of two unit vectors) while the
//! rotation column carries units of length, because a rotation partial scales with the distance
//! from the local origin. `H` is therefore heterogeneous, and two consequences follow.
//!
//! [`AlignInfo2::marginal_precision`] is safe: each entry concerns a single parameter, so
//! its units are internally consistent. [`AlignInfo2::weak_directions`] mixes translation
//! and rotation in one vector and is only meaningful when the two are comparably scaled, which in
//! practice means the local origin should sit near the middle of the point set (as
//! `AlignParams2::from_center` arranges). Placing the origin far away inflates the rotation term
//! and will make rotation look artificially well constrained.

use crate::Point2;
use crate::common::align::information::{AlignDofs, AlignInfo};
use crate::geom2::align2::jacobian::point_surf_jacobian;
use crate::geom2::align2::{AlignParams2, Dof3, SurfaceTarget2};

/// The information content of a set of test points with respect to a 2D alignment, over the three
/// parameters `tx`, `ty`, `rz`.
///
/// Built with [`AlignInfo2::from_points`], which projects the points onto a target and
/// takes the jacobian row of each correspondence, or with [`AlignInfo2::from_rows`] when
/// the rows have already been computed. See the module documentation for what it is good for.
pub type AlignInfo2 = AlignInfo<Dof3, 3>;

impl AlignDofs<3> for Dof3 {
    fn free_indices(&self) -> Vec<usize> {
        [self.tx, self.ty, self.rz]
            .iter()
            .enumerate()
            .filter_map(|(i, &free)| free.then_some(i))
            .collect()
    }
}

impl AlignInfo2 {
    /// Builds the information content of a set of points against a target, in the pose described
    /// by `params`.
    ///
    /// Each point is moved by the current transform, projected onto the target, and turned into a
    /// jacobian row exactly as [`crate::geom2::align2::points_to_surface2`] would do. A point is
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
    ///   partial is taken about and which degrees of freedom are free
    pub fn from_points(
        points: &[Point2],
        target: &impl SurfaceTarget2,
        params: &AlignParams2,
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
    use crate::geom2::align2::{AlignOrigin2, AlignStorage2, AlignSurfMatch2};
    use crate::na::Matrix3;
    use crate::{Curve2, Iso2, UnitVec2, Vector2};
    use approx::assert_relative_eq;

    /// A flat target along `y = 0`, so every match reports the same upward normal. Its
    /// information content is known in closed form, which is what makes it useful here.
    struct GroundLine;

    impl SurfaceTarget2 for GroundLine {
        fn find_align_match(&self, p: &Point2) -> AlignSurfMatch2 {
            AlignSurfMatch2::new(
                Point2::new(p.x, 0.0),
                UnitVec2::new_unchecked(Vector2::y()),
                true,
                1.0,
            )
        }
    }

    /// Points spread along `y = 1`, one unit above `GroundLine`.
    fn line_points() -> Vec<Point2> {
        (-6..=6).map(|i| Point2::new(i as f64, 1.0)).collect()
    }

    /// A closed counter-clockwise rectangle 10 by 5, centered on the origin. Being closed, it
    /// constrains every degree of freedom.
    fn rect_curve() -> Curve2 {
        let points =
            [(-5.0, -2.5), (5.0, -2.5), (5.0, 2.5), (-5.0, 2.5)].map(|(x, y)| Point2::new(x, y));
        Curve2::from_points(&points, 1e-8, true).unwrap()
    }

    /// Points walked around the rectangle at a uniform arc length spacing. The perimeter is 30,
    /// so a spacing of 0.3 gives 100 of them, and the first twenty all sit on the bottom edge.
    fn rect_points(curve: &Curve2) -> Vec<Point2> {
        let spacing = 0.3;
        let n = (curve.length() / spacing).floor() as usize;
        (0..n)
            .filter_map(|k| {
                curve
                    .at_length((k as f64 + 0.5) * spacing)
                    .map(|s| s.point())
            })
            .collect()
    }

    // ============================================================================================
    // Conditioning against geometry whose answer is known in advance
    // ============================================================================================

    #[test]
    fn a_straight_line_leaves_one_motion_unconstrained() {
        // Points against a straight line can only ever resist motion across it. Sliding along the
        // line moves every point parallel to the surface and changes no residual, so exactly one
        // of the three degrees of freedom is free. Rotation is not free, because turning the set
        // lifts points away from the line in proportion to how far out they sit.
        let info = AlignInfo2::from_points(
            &line_points(),
            &GroundLine,
            &AlignParams2::from_origin(None),
        );

        let weak = info.weak_directions();
        assert_eq!(weak.len(), 3);

        assert!(
            weak[0].0.abs() < 1e-9,
            "expected an unconstrained motion, got {}",
            weak[0].0
        );
        assert!(
            weak[1].0 > 1.0,
            "expected the remaining motions to be constrained, got {}",
            weak[1].0
        );

        // ...and that motion is sliding along the line, with no ty or rz component at all.
        let direction = weak[0].1;
        assert_relative_eq!(direction[0].abs(), 1.0, epsilon = 1e-9);
        assert!(direction[1].abs() < 1e-9);
        assert!(direction[2].abs() < 1e-9);
    }

    #[test]
    fn a_singular_problem_reports_rather_than_lying() {
        let info = AlignInfo2::from_points(
            &line_points(),
            &GroundLine,
            &AlignParams2::from_origin(None),
        );

        let err = info.marginal_precision().unwrap_err().to_string();
        assert!(err.contains("singular"), "unexpected message: {err}");
        assert!(info.leverage().is_err());
    }

    #[test]
    fn locking_the_free_motion_makes_a_line_well_conditioned() {
        // The same points, but with the motion they cannot see locked out. What remains is fully
        // determined, so the analysis should now succeed.
        let dof = Dof3::new(false, true, true);
        let params = AlignParams2::from_origin(Some(dof));
        let info = AlignInfo2::from_points(&line_points(), &GroundLine, &params);

        assert_eq!(info.free_dof_count(), 2);

        let precision = info.marginal_precision().unwrap();
        assert_eq!(precision[0], None, "tx was locked");
        for i in [1usize, 2] {
            assert!(
                precision[i].unwrap() > 0.0,
                "parameter {i} should be constrained"
            );
        }
    }

    #[test]
    fn a_closed_curve_constrains_everything() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let params = AlignParams2::from_origin(None);
        let info = AlignInfo2::from_points(&points, &curve, &params);

        assert_eq!(info.free_dof_count(), 3);
        let precision = info.marginal_precision().unwrap();
        for (i, p) in precision.iter().enumerate() {
            assert!(
                p.unwrap() > 0.0,
                "parameter {i} should be constrained on a closed curve"
            );
        }
    }

    // ============================================================================================
    // Leverage
    // ============================================================================================

    #[test]
    fn leverage_sums_to_the_free_dof_count() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let info = AlignInfo2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 3.0, epsilon = 1e-9);
    }

    #[test]
    fn leverage_sums_to_the_free_dof_count_when_constrained() {
        // The identity `sum(h_i) == k` has to hold on the free subspace too, not just at k = 3.
        let curve = rect_curve();
        let points = rect_points(&curve);
        let dof = Dof3::new(false, true, true);
        let info = AlignInfo2::from_points(&points, &curve, &AlignParams2::from_origin(Some(dof)));

        assert_eq!(info.free_dof_count(), 2);
        let total: f64 = info.leverage().unwrap().iter().sum();
        assert_relative_eq!(total, 2.0, epsilon = 1e-9);
    }

    // ============================================================================================
    // Greedy selection
    // ============================================================================================

    #[test]
    fn selection_prefers_isolated_points_over_a_redundant_cluster() {
        // This is the property that separates greedy selection from independent scoring. A tight
        // cluster of near-identical points all contribute the same row, so after the first one is
        // taken the rest add nothing; the isolated points elsewhere on the rectangle each bring
        // something new and should be chosen first.
        let curve = rect_curve();

        let mut points = Vec::new();

        // Twenty near-duplicates in one spot on the bottom edge.
        for i in 0..20 {
            points.push(Point2::new(1.0 + i as f64 * 1e-4, -2.5));
        }
        let cluster = 0..20;

        // Four isolated points spread over the other three edges.
        let isolated = [
            Point2::new(5.0, 0.0),
            Point2::new(-5.0, 1.0),
            Point2::new(0.0, 2.5),
            Point2::new(-3.0, 2.5),
        ];
        points.extend_from_slice(&isolated);

        let info = AlignInfo2::from_points(&points, &curve, &AlignParams2::from_origin(None));
        // Four picks, which is how many distinct positions are on offer once the cluster has
        // contributed the one row it has to give. Twenty of the twenty-four candidates are
        // cluster members, so picking arbitrarily would take about three of them.
        let picked = info.select_d_optimal(4, None);
        assert_eq!(picked.len(), 4);

        let from_cluster = picked.iter().filter(|&&i| cluster.contains(&i)).count();
        assert_eq!(
            from_cluster, 1,
            "expected exactly one cluster member among {picked:?}"
        );

        // Asked for more than the geometry can offer, it does go back to the cluster: with three
        // degrees of freedom and every distinct position already taken, a near-duplicate is the
        // best remaining candidate. That is the selection working, not failing.
        let more = info.select_d_optimal(5, None);
        assert_eq!(
            more.iter().filter(|&&i| cluster.contains(&i)).count(),
            2,
            "expected the fifth pick to fall back to the cluster: {more:?}"
        );
    }

    #[test]
    fn selection_is_a_prefix_so_it_can_be_truncated() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let info = AlignInfo2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        let long = info.select_d_optimal(20, None);
        let short = info.select_d_optimal(8, None);

        assert_eq!(&long[..8], &short[..]);
    }

    #[test]
    fn a_selected_subset_stays_well_conditioned() {
        // The point of pruning: a small chosen subset should hold the alignment up nearly as well
        // as the whole set, and far better than the same number of arbitrary points.
        let curve = rect_curve();
        let points = rect_points(&curve);
        let params = AlignParams2::from_origin(None);
        let info = AlignInfo2::from_points(&points, &curve, &params);

        let n = 20;
        let picked = info.select_d_optimal(n, None);
        let chosen: Vec<Point2> = picked.iter().map(|&i| points[i]).collect();
        let arbitrary: Vec<Point2> = points.iter().take(n).copied().collect();

        // The comparison is on the smallest eigenvalue, which measures how close the subset comes
        // to leaving some motion unconstrained. It is used rather than `marginal_precision`
        // because an arbitrary subset can be outright singular (the first `n` samples walk along
        // a single edge of the rectangle), and there is no finite precision to compare then.
        let weakest = |pts: &[Point2]| -> f64 {
            AlignInfo2::from_points(pts, &curve, &params).weak_directions()[0].0
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
            AlignInfo2::from_points(&chosen, &curve, &params)
                .marginal_precision()
                .is_ok(),
            "the selected subset should leave no motion unconstrained"
        );
    }

    #[test]
    fn selection_respects_the_available_count() {
        let curve = rect_curve();
        let points = rect_points(&curve);
        let info = AlignInfo2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        assert!(info.select_d_optimal(0, None).is_empty());
        assert_eq!(
            info.select_d_optimal(points.len() + 50, None).len(),
            points.len()
        );
    }

    #[test]
    fn zero_weight_points_are_never_selected() {
        let rows = vec![
            AlignStorage2::new(1.0, 0.0, 0.0),
            AlignStorage2::new(0.0, 1.0, 0.0),
            AlignStorage2::new(0.0, 0.0, 1.0),
        ];
        let weights = vec![1.0, 0.0, 1.0];

        let info = AlignInfo2::from_rows(rows, weights, Dof3::all()).unwrap();
        let picked = info.select_d_optimal(3, None);

        assert!(!picked.contains(&1), "a zero-weight point was selected");
        assert_eq!(picked.len(), 2);
    }

    // ============================================================================================
    // Construction and validation
    // ============================================================================================

    #[test]
    fn information_matrix_is_expanded_with_zeros_for_locked_dof() {
        let rows = vec![AlignStorage2::new(1.0, 2.0, 0.0)];
        let dof = Dof3::new(true, true, false);
        let info = AlignInfo2::from_rows(rows, vec![1.0], dof).unwrap();

        let h = info.information();
        // The free block is the outer product of (1, 2) with itself.
        assert_relative_eq!(h[(0, 0)], 1.0, epsilon = 1e-12);
        assert_relative_eq!(h[(0, 1)], 2.0, epsilon = 1e-12);
        assert_relative_eq!(h[(1, 1)], 4.0, epsilon = 1e-12);
        // ...and every locked row and column is zero.
        for j in 0..3 {
            assert_eq!(h[(2, j)], 0.0);
            assert_eq!(h[(j, 2)], 0.0);
        }
    }

    #[test]
    fn mismatched_rows_and_weights_are_rejected() {
        let rows = vec![AlignStorage2::zeros(); 3];
        let err = AlignInfo2::from_rows(rows, vec![1.0; 2], Dof3::all())
            .unwrap_err()
            .to_string();
        assert!(err.contains("weights"), "unexpected message: {err}");
    }

    #[test]
    fn invalid_weights_are_rejected() {
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            let rows = vec![AlignStorage2::zeros(); 2];
            let result = AlignInfo2::from_rows(rows, vec![1.0, bad], Dof3::all());
            assert!(result.is_err(), "weight {bad} should have been rejected");
        }
    }

    #[test]
    fn a_fully_locked_problem_degrades_gracefully() {
        let dof = Dof3::new(false, false, false);
        let info =
            AlignInfo2::from_rows(vec![AlignStorage2::zeros(); 4], vec![1.0; 4], dof).unwrap();

        assert_eq!(info.free_dof_count(), 0);
        assert!(info.weak_directions().is_empty());
        assert!(info.select_d_optimal(4, None).is_empty());
        assert_eq!(info.information(), Matrix3::zeros());
    }

    #[test]
    fn the_local_origin_is_respected() {
        // The rotation partial is taken about the local origin, so moving it changes the rotation
        // entry of the information matrix. This guards against the origin being quietly ignored.
        let curve = rect_curve();
        let points = rect_points(&curve);

        let at_origin = AlignInfo2::from_points(&points, &curve, &AlignParams2::from_origin(None));
        let offset = AlignInfo2::from_points(
            &points,
            &curve,
            &AlignParams2::from_center(Point2::new(50.0, 0.0), None),
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
        let curve = rect_curve();
        let points = rect_points(&curve);

        let here = AlignInfo2::from_points(&points, &curve, &AlignParams2::from_origin(None));

        let displaced = AlignParams2::new(
            AlignOrigin2::Origin,
            Some(Iso2::translation(0.0, 20.0)),
            None,
        );
        let there = AlignInfo2::from_points(&points, &curve, &displaced);

        let delta = (here.information() - there.information()).norm();
        assert!(
            delta > 1e-6,
            "the working offset had no effect (delta {delta})"
        );
    }
}
