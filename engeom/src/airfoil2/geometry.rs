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
    AfEdge, AfEdgeSearch, AfGeometry, OrientFwdAft, OrientUpperLower, SectionInput,
};
use crate::{Curve2, Result};

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

    todo!()
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
    let sharp = fit_sharp_edge(input, circles.clone(), at_front)?;
    let square = fit_square_edge(input, circles.clone(), at_front)?;
    let rounded_square = fit_rounded_square_edge(input, circles.clone(), at_front)?;
    let full_round = fit_full_round_edge(input, circles.clone(), at_front)?;
    let blended_round = fit_blended_round_edge(input, circles.clone(), at_front)?;

    let mut candidates = vec![sharp, full_round, blended_round];
    if rounded_square.avg_residual < square.avg_residual * 0.9 {
        candidates.push(rounded_square);
    } else {
        candidates.push(square);
    }

    candidates.sort_by(|a, b| {
        a.avg_residual
            .partial_cmp(&b.avg_residual)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(candidates[0].clone())
}
