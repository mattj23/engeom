//! This module has the composed internal tools for doing the geometric analysis of airfoil
//! sections.  These tools pull together the lower-level tools in the sibling modules to make
//! more convenient ways to extract geometry for most use cases.  The user can always still use
//! the individual pieces manually.

use crate::airfoil2::camber::extract_inscribed_circles;
use crate::airfoil2::edges::{
    AfEdgeFit, fit_blended_round_edge, fit_full_round_edge, fit_rounded_square_edge,
    fit_sharp_edge, fit_square_edge,
};
use crate::airfoil2::inscribed::Inscribed;
use crate::airfoil2::{
    AfEdge, AfEdgeGeometry, AfEdgeSearch, AfGeometry, OrientFwdAft, OrientUpperLower, SectionInput,
};
use crate::{Curve2, Point2, Result};

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
    let (oriented, leading) = run_edge_fit(le_search, &input, oriented, true)?;
    let (oriented, trailing) = run_edge_fit(te_search, &input, oriented, false)?;

    // Now we extract the camber curve
    let camber = make_camber_curve(leading.point, trailing.point, &oriented, general_tol * 1e-2)?;

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
    let (lower, upper) =
        if side0.dist_to_point(&oriented[0].p0) < side1.dist_to_point(&oriented[0].p0) {
            (side0, side1)
        } else {
            (side1, side0)
        };

    Ok(AfGeometry {
        leading,
        trailing,
        camber,
        upper,
        lower,
        circles: oriented,
    })
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
        AfEdgeSearch::Open => todo!(),
        AfEdgeSearch::Sharp => fit_sharp_edge(input, circles, at_front),
        AfEdgeSearch::Square => fit_square_edge(input, circles, at_front),
        AfEdgeSearch::RoundedSquare => fit_rounded_square_edge(input, circles, at_front),
        AfEdgeSearch::FullRound => fit_full_round_edge(input, circles, at_front),
        AfEdgeSearch::BlendedRound => fit_blended_round_edge(input, circles, at_front),
    }?;

    Ok((result.circles, result.edge))
}

fn fit_auto_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    let mut candidates = Vec::new();

    if let Ok(sharp) = fit_sharp_edge(input, circles.clone(), at_front) {
        candidates.push(sharp);
    }
    if let Ok(square) = fit_square_edge(input, circles.clone(), at_front) {
        candidates.push(square);
    }
    if let Ok(rounded_square) = fit_rounded_square_edge(input, circles.clone(), at_front) {
        candidates.push(rounded_square);
    }
    if let Ok(full_round) = fit_full_round_edge(input, circles.clone(), at_front) {
        candidates.push(full_round);
    }
    if let Ok(blended_round) = fit_blended_round_edge(input, circles.clone(), at_front) {
        candidates.push(blended_round);
    }
    if candidates.len() == 0 {
        return Err("No edge fits succeeded for auto edge search".into());
    }

    candidates.sort_by(|a, b| {
        a.avg_residual
            .partial_cmp(&b.avg_residual)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(candidates[0].clone())
}
