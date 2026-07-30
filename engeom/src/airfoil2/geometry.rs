//! This module has the composed internal tools for doing the geometric analysis of airfoil
//! sections.  These tools pull together the lower-level tools in the sibling modules to make
//! more convenient ways to extract geometry for most use cases.  The user can always still use
//! the individual pieces manually.

use crate::airfoil2::camber::extract_inscribed_circles;
use crate::airfoil2::edges::{
    AfEdgeFit, fit_blended_round_edge, fit_full_round_edge, fit_open_edge, fit_rounded_square_edge,
    fit_sharp_edge, fit_spline_max_k, fit_square_edge,
};
use crate::airfoil2::inscribed::Inscribed;
use crate::airfoil2::{
    AfEdge, AfEdgeGeometry, AfEdgeSearch, AfGeometry, OrientFwdAft, OrientUpperLower, SectionInput,
};
use crate::geom2::LineOps2;
use crate::{AngleDir, Curve2, Point2, Result};

pub fn geometry_only_analysis(
    section: &Curve2,
    general_tol: f64,
    fwd_aft: OrientFwdAft,
    upper_lower: OrientUpperLower,
    le_search: AfEdgeSearch,
    te_search: AfEdgeSearch,
) -> Result<AfGeometry> {
    // Create the unified input struct and then start by capturing the unambiguous inscribed
    // circles. The orientation of the inscribed circles is completely unknown.
    let input = SectionInput::new(section, general_tol);
    let unoriented = extract_inscribed_circles(&input)?;

    // Now we will orient the inscribed circles in the forward/aft direction
    let fwd_oriented = fwd_aft.apply(unoriented)?;

    // Now we orient them in the upper/lower direction
    let oriented = upper_lower.apply(fwd_oriented)?;

    // Now we capture the leading and trailing edges
    let (oriented, trailing) = run_edge_fit(te_search, &input, oriented, false)?;
    let (oriented, leading) = run_edge_fit(le_search, &input, oriented, true)?;

    // Now we extract the camber curve
    let camber = make_camber_curve(leading.point, trailing.point, &oriented, general_tol * 1e-2)?;
    let (lower, upper) = split_sides(leading, trailing, section, &oriented)?;

    Ok(AfGeometry {
        leading,
        trailing,
        camber,
        upper,
        lower,
        circles: oriented,
    })
}

fn split_sides(
    leading: AfEdge,
    trailing: AfEdge,
    section: &Curve2,
    oriented: &[Inscribed],
) -> Result<(Curve2, Curve2)> {
    // Now we split the two sides
    let (side0, side1) = match (leading.geometry, trailing.geometry) {
        (AfEdgeGeometry::Open, AfEdgeGeometry::Open) => {
            // This shouldn't be possible
            return Err("Both leading and trailing edges cannot be open".into());
        }
        (AfEdgeGeometry::Open, _) => {
            // Leading is open, trailing is closed
            let l = section.at_closest_to_point(&trailing.point).length_along();
            section.split_open_at_length(l)?
        }
        (_, AfEdgeGeometry::Open) => {
            // Leading is closed, trailing is open
            let l = section.at_closest_to_point(&leading.point).length_along();
            section.split_open_at_length(l)?
        }
        (_, _) => {
            // Both sides are closed
            let l0 = section.at_closest_to_point(&leading.point).length_along();
            let l1 = section.at_closest_to_point(&trailing.point).length_along();
            section.split_closed_at_lengths(l0, l1)?
        }
    };

    // To identify the sides, we'll check an inscribed circle and see which one is closer to `p0`,
    // which would indicate the lower/pressure side.
    let test_circle = &oriented[oriented.len() / 2];
    let (lower, upper) =
        if side0.dist_to_point(&test_circle.p0) < side1.dist_to_point(&test_circle.p0) {
            (side0, side1)
        } else {
            (side1, side0)
        };

    // Ensure that the upper and lower skins are oriented in the correct direction so that their
    // surface normals face outward
    let line_l = lower.at_closest_to_point(&test_circle.c).direction_line();
    let lower = if line_l.winding_direction(&test_circle.c) == AngleDir::Cw {
        lower.reversed()
    } else {
        lower
    };

    let line_u = upper.at_closest_to_point(&test_circle.c).direction_line();
    let upper = if line_u.winding_direction(&test_circle.c) == AngleDir::Cw {
        upper.reversed()
    } else {
        upper
    };

    Ok((lower, upper))
}

fn make_camber_curve(
    le_point: Point2,
    te_point: Point2,
    oriented_circles: &[Inscribed],
    tol: f64,
) -> Result<Curve2> {
    let mut points = vec![le_point];
    for c in oriented_circles.iter() {
        points.push(c.center());
    }
    points.push(te_point);
    Curve2::from_points(&points, tol, false)
}

fn run_edge_fit(
    search: AfEdgeSearch,
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<(Vec<Inscribed>, AfEdge)> {
    let result = match search {
        AfEdgeSearch::Auto => fit_auto_edge(input, circles, at_front),
        AfEdgeSearch::Open => fit_open_edge(input, circles, at_front),
        AfEdgeSearch::Sharp => fit_sharp_edge(input, circles, at_front),
        AfEdgeSearch::Square => fit_square_edge(input, circles, at_front),
        AfEdgeSearch::RoundedSquare => fit_rounded_square_edge(input, circles, at_front),
        AfEdgeSearch::FullRound => fit_full_round_edge(input, circles, at_front),
        AfEdgeSearch::BlendedRound => fit_blended_round_edge(input, circles, at_front),
        AfEdgeSearch::SplineMaxK => fit_spline_max_k(input, circles, at_front).map(|r| r.0),
    }?;

    Ok((result.circles, result.edge))
}

struct FitCandidate {
    complexity: usize,
    fit: AfEdgeFit,
}

impl FitCandidate {
    fn new(complexity: usize, fit: AfEdgeFit) -> Self {
        Self { complexity, fit }
    }
}

fn fit_auto_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // Candidates have a complexity score based on the type of edge fit, with the goal being to
    // use the lowest complexity method which adequately captures the geometry. For a method to be
    // selected, it must have an average residual at least some tolerance better than the best
    // method with a lower complexity score.
    const COMPLEXITY_TOL: f64 = 0.05; // Starting with 5%

    let mut candidates = Vec::new();
    if let Ok(sharp) = fit_sharp_edge(input, circles.clone(), at_front) {
        candidates.push(FitCandidate::new(0, sharp));
    }
    if let Ok(square) = fit_square_edge(input, circles.clone(), at_front) {
        candidates.push(FitCandidate::new(0, square));
    }
    if let Ok(full_round) = fit_full_round_edge(input, circles.clone(), at_front) {
        candidates.push(FitCandidate::new(0, full_round));
    }

    if let Ok(rounded_square) = fit_rounded_square_edge(input, circles.clone(), at_front) {
        candidates.push(FitCandidate::new(1, rounded_square));
    }

    if let Ok(blended_round) = fit_blended_round_edge(input, circles.clone(), at_front) {
        candidates.push(FitCandidate::new(1, blended_round));
    }

    if let Ok(spline_max) = fit_spline_max_k(input, circles.clone(), at_front) {
        candidates.push(FitCandidate::new(2, spline_max.0));
    }

    // Drop any fit whose quality metric is not finite (e.g. an empty residual set yields a NaN
    // mean), since it can neither be ranked nor meaningfully compared against other tiers.
    candidates.retain(|c| c.fit.avg_residual.is_finite());
    if candidates.is_empty() {
        return Err("No edge fits with a valid residual for auto edge search".into());
    }

    // Walk the complexity tiers from lowest to highest. A higher-complexity tier is only adopted
    // when its best residual improves on the best residual achieved by *any* lower-complexity tier
    // by at least `COMPLEXITY_TOL` (a relative fraction), so that added complexity has to pay for
    // itself. `best_lower` accumulates the best residual seen strictly below the current tier and
    // is independent of which tier ended up selected.
    let max_complexity = candidates.iter().map(|c| c.complexity).max().unwrap();

    let mut selected: Option<&FitCandidate> = None;
    let mut best_lower = f64::INFINITY;
    for tier in 0..=max_complexity {
        let best_in_tier = candidates
            .iter()
            .filter(|c| c.complexity == tier)
            .min_by(|a, b| a.fit.avg_residual.total_cmp(&b.fit.avg_residual));

        if let Some(best) = best_in_tier {
            // The first (lowest) tier present is always taken as the baseline; later tiers must
            // beat the best lower residual by the tolerance. Note that a lower tier achieving a
            // residual of exactly zero locks the selection, since nothing can beat zero by a
            // relative margin.
            if selected.is_none() || best.fit.avg_residual < best_lower * (1.0 - COMPLEXITY_TOL) {
                selected = Some(best);
            }
            best_lower = best_lower.min(best.fit.avg_residual);
        }
    }

    Ok(selected
        .expect("non-empty candidates guarantees a selection")
        .fit
        .clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::airfoil_curve;

    fn af_curve_known_ccw() -> Curve2 {
        let donor = airfoil_curve();
        Curve2::from_points_ccw(donor.points(), donor.tol(), true).unwrap()
    }

    fn basic_geom(section: &Curve2) -> Result<(Vec<Inscribed>, AfEdge, AfEdge)> {
        let input = SectionInput::new(section, 1e-3);
        let unoriented = extract_inscribed_circles(&input)?;
        let partial = OrientFwdAft::TmaxFwd.apply(unoriented)?;
        let oriented = OrientUpperLower::Curvature.apply(partial)?;

        // Now we capture the leading and trailing edges
        let (oriented, trailing) = run_edge_fit(AfEdgeSearch::FullRound, &input, oriented, false)?;
        let (oriented, leading) = run_edge_fit(AfEdgeSearch::BlendedRound, &input, oriented, true)?;
        Ok((oriented, leading, trailing))
    }

    fn verify_upper_lower(upper: &Curve2, lower: &Curve2, expected: &Curve2) -> Result<()> {
        for p in upper.points() {
            let local = upper.at_closest_to_point(p);
            let check = expected.at_closest_to_point(p);
            assert!(
                local.normal().dot(&check.normal()).is_sign_positive(),
                "Upper point is not facing the right way"
            );
        }

        for p in lower.points() {
            let local = lower.at_closest_to_point(p);
            let check = expected.at_closest_to_point(p);
            assert!(
                local.normal().dot(&check.normal()).is_sign_positive(),
                "Lower point is not facing the right way"
            );
        }

        Ok(())
    }

    #[test]
    fn side_orientation_ccw() -> Result<()> {
        let expected = af_curve_known_ccw();
        let (oriented, leading, trailing) = basic_geom(&expected)?;
        let (lower, upper) = split_sides(leading, trailing, &expected, &oriented)?;
        verify_upper_lower(&upper, &lower, &expected)?;
        Ok(())
    }

    #[test]
    fn side_orientation_cw() -> Result<()> {
        let expected = af_curve_known_ccw();
        let reversed = expected.reversed();

        let (oriented, leading, trailing) = basic_geom(&reversed)?;
        let (lower, upper) = split_sides(leading, trailing, &reversed, &oriented)?;

        verify_upper_lower(&upper, &lower, &expected)?;
        Ok(())
    }

    #[test]
    fn auto_edges() -> Result<()> {
        let curve = airfoil_curve();
        let input = SectionInput::new(&curve, 1e-3);

        let circles = extract_inscribed_circles(&input)?;
        let circles = OrientFwdAft::TmaxFwd.apply(circles)?;
        let circles = OrientUpperLower::Curvature.apply(circles)?;

        let result = fit_auto_edge(&input, circles, true)?;
        let _result = fit_auto_edge(&input, result.circles, false)?;

        Ok(())
    }
}
