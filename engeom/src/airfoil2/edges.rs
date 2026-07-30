use crate::AngleDir::Ccw;
use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::airfoil2::{AfEdge, AfEdgeGeometry, SectionInput};
use crate::common::cubic_spline::{SplineBuildFn, fit_spline_to_points};
use crate::common::{PCoords, dist, mid_point};
use crate::geom2::{
    BndBuildFn, BoundaryData2, BoundaryEditor, BoundaryElement2, CubicSpline2, LineOps2, Segment2,
    fit_boundary_to_points, signed_angle,
};
use crate::{Arc2, Circle2, Curve2, DVector, Iso2, Line2, Point2, Result, SurfacePoint2, Vector2};
use std::f64::consts::PI;

/// This is the result of an airfoil edge fitting operation. It contains the detected edge point
/// and geometry, any fitting point residuals left from the result geometry, and an updated stack
/// of inscribed circles, which may have been refined during the edge fitting.
#[derive(Clone, Debug)]
pub struct AfEdgeFit {
    /// The fitted edge: a canonical edge location and its geometric description.
    pub edge: AfEdge,

    /// Point-to-boundary residuals from the fitting optimization, one entry per section point
    /// used in the fit.
    pub residuals: DVector,

    /// The inscribed circle stack as it stood at the end of the fit. Some edge fitters
    /// (e.g. blended-round) extend the stack with additional refined circles near the edge.
    pub circles: Vec<Inscribed>,

    /// The mean of `residuals`, cached for cheap comparison between candidate fits.
    pub avg_residual: f64,
}

impl AfEdgeFit {
    /// Build a new `AfEdgeFit` and compute `avg_residual` as the mean of `residuals`.
    pub fn new(result: AfEdge, residuals: DVector, circles: Vec<Inscribed>) -> Self {
        let avg_residual = residuals.mean();
        Self {
            edge: result,
            residuals,
            circles,
            avg_residual,
        }
    }
}

/// Fit an open trailing or leading edge on the airfoil section data.  An open edge is one where
/// the airfoil shape is incomplete, which can happen for a number of different reasons. Most
/// commonly, this will be when either the measurement system does not capture all the way around
/// the airfoil (such as when only part of the airfoil is of interest) or when the nominal geometry
/// is not a full airfoil (sometimes occurring near the root or tip where the airfoil is blending
/// into other geometry).
///
/// Because there is no edge geometry to fit, the edge point is instead found by projecting the end
/// of the camber line onto a chord spanning the gap, from the first point of the section curve to
/// the last. To keep that projection from being dominated by the direction of a camber line which
/// stopped well short of the gap, the camber line is first advanced into the gap by repeatedly
/// stepping halfway to the nearer lip and pushing a new refined inscribed circle, until the
/// projected edge point stops moving by more than the general tolerance.
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, process the front (leading edge) end of the airfoil; when `false`,
///   process the rear (trailing edge) end.
///
/// # Returns
///
/// An [`AfEdgeFit`] whose edge geometry is [`AfEdgeGeometry::Open`], with an empty residual vector
/// (there is no fitted geometry to measure against) and an inscribed circle stack extended toward
/// the gap.
pub fn fit_open_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // Each step halves the distance remaining to the gap, so this is a generous upper bound on the
    // number of steps needed to converge for any sane tolerance.
    const MAX_ITERATIONS: usize = 20;

    let mut working = EdgeWork::new(input, circles, at_front)?;
    let cap = open_gap_segment(input)?;

    let mut point = end_cap_intersection(&working.stack, &cap)?;
    let mut drift = f64::INFINITY;

    for _ in 0..MAX_ITERATIONS {
        if drift <= input.general_tol {
            break;
        }

        let sp = end_camber_ray(&working.stack)?;
        let max_dist = sp
            .scalar_projection(&cap.a)
            .min(sp.scalar_projection(&cap.b));
        if max_dist < 0.0 {
            return Err("The open gap does not lie ahead of the end of the camber line".into());
        }

        // Probe halfway to the nearer lip of the gap. The test line is oriented so that the `p0`
        // and `p1` of the resulting circle come out in the same order as the rest of the stack,
        // since `refine_and_push` does not orient what it's given.
        let probe = Line2::from(&sp.shifted(max_dist * 0.5).normal_rotated_90(Ccw));
        let test_line = if probe.direction.dot(&working.last()?.contact_dir()) < 0.0 {
            probe.reversed()
        } else {
            probe
        };

        // Failing to find a crossing line means the camber has been advanced as far as the section
        // data supports, so we stop and keep the best edge point found so far.
        let Some(inscribed) = input.try_inscribed(&test_line) else {
            break;
        };
        working.stack.refine_and_push(inscribed, input);

        let next = end_cap_intersection(&working.stack, &cap)?;
        drift = dist(&point, &next);
        point = next;
    }

    let c = working.take_circles();
    let edge = AfEdge::new(point, AfEdgeGeometry::Open);

    // The residuals are deliberately empty. An open edge has no fitted geometry, so it cannot be
    // ranked against the closed-edge fits by residual; callers which compare fits must resolve the
    // open case separately rather than letting it into the comparison.
    Ok(AfEdgeFit::new(edge, DVector::zeros(0), c))
}

/// Determine whether the section is open at the end of the inscribed circle stack being processed.
///
/// Whether an edge is open is a property of the section topology rather than something which can be
/// measured against a candidate geometry, so this is the test used to resolve the open case before
/// any edge fitting is attempted.
///
/// A section is open at this end when it is not a closed curve *and* both lips of the gap lie beyond
/// the end clipping line of the stack. Requiring both lips is what picks out the correct end: the
/// leading and trailing clipping lines face in opposite directions, so at most one end can pass. A
/// section with a gap in the middle of a surface passes at neither end and is left to the normal
/// edge fitting, which matches the input contract documented on
/// [`extract_inscribed_circles`](crate::airfoil2::camber::extract_inscribed_circles).
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, test the front (leading edge) end of the airfoil; when `false`, test
///   the rear (trailing edge) end.
///
/// # Returns
///
/// `true` if the section is open at the requested end. Errors if the stack has fewer than two
/// inscribed circles, which is too few to establish an end clipping line.
pub fn is_open_at_end(input: &SectionInput, circles: &[Inscribed], at_front: bool) -> Result<bool> {
    if input.section.is_closed() {
        return Ok(false);
    }

    let mut stack = InscribedVec::new(circles.to_vec());
    if at_front {
        stack.reverse_order_only();
    }
    let clip = stack.end_clip_line()?;

    let a = input.section.at_front().point();
    let b = input.section.at_back().point();
    Ok(clip.scalar_project(&a) > 0.0 && clip.scalar_project(&b) > 0.0)
}

/// Fit a square (flat) trailing or leading edge to airfoil section data.
///
/// The edge geometry consists of two corner points connected by a straight flat face, with the
/// edge point placed at the midpoint between them. The corners are optimized by fitting an open
/// three-segment polyline (from `p0` through two free corner points to `p1`) to the section
/// points that lie beyond the clipping line of the last inscribed circle.
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, process the front (leading edge) end of the airfoil; when `false`,
///   process the rear (trailing edge) end.
///
/// # Returns
///
/// An [`AfEdgeFit`] whose edge geometry is [`AfEdgeGeometry::Square`], containing the two fitted
/// corner points and the midpoint as the edge location.
pub fn fit_square_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0;
    let p1 = working.last()?.p1;

    // To find the initial, we'll find how far the farthest point is from the end clipping line
    // and add that vector to p0 and p1
    let d = working.clip.direction * working.clip_max_scalar();
    let i0 = p0 + d;
    let i1 = p1 + d;
    let initial = DVector::from(vec![i0.x, i0.y, i1.x, i1.y]);

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        let mut bdata = BoundaryData2::new_open(p0);
        bdata.add_seg_xy(params[0], params[1]);
        bdata.add_seg_xy(params[2], params[3]);
        bdata.add_seg(&p1);
        bdata.try_to_boundary()
    });

    let result = fit_boundary_to_points(&working.fit_points, &builder, initial, false)?;
    let corner0 = Point2::new(result.params[0], result.params[1]);
    let corner1 = Point2::new(result.params[2], result.params[3]);
    let point = end_intersection(input, &working.last()?.c, &mid_point(&corner0, &corner1))?;

    let c = working.take_circles();
    let edge = AfEdge::new(point, AfEdgeGeometry::Square(corner0, corner1));
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

/// Fit a rounded-square trailing or leading edge to airfoil section data.
///
/// Like [`fit_square_edge`], but the two corners connecting the flat face to the airfoil surfaces
/// are replaced by circular arc fillets of a single optimized radius. The boundary is an open
/// path from `p0` through two filleted corners to `p1`, fit to the section points beyond the
/// clipping line of the last inscribed circle. The edge point is the midpoint between the two
/// unfilleted corner positions.
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, process the front (leading edge) end of the airfoil; when `false`,
///   process the rear (trailing edge) end.
///
/// # Returns
///
/// An [`AfEdgeFit`] whose edge geometry is [`AfEdgeGeometry::RoundedSquare`], containing the two
/// unfilleted corner positions, the fillet radius, and the midpoint as the edge location.
pub fn fit_rounded_square_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0;
    let p1 = working.last()?.p1;

    // To find the initial, we'll find how far the farthest point is from the end clipping line
    // and add that vector to p0 and p1
    let d = working.clip.direction * working.clip_max_scalar();
    let i0 = p0 + d;
    let i1 = p1 + d;
    let initial = DVector::from(vec![i0.x, i0.y, i1.x, i1.y, dist(&p0, &p1) / 8.0]);

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        let mut bdata = BoundaryData2::new_open(p0);
        bdata.add_corner_fillets(
            &[
                Point2::new(params[0], params[1]),
                Point2::new(params[2], params[3]),
                p1,
            ],
            params[4],
        )?;
        bdata.try_to_boundary()
    });

    let result = fit_boundary_to_points(&working.fit_points, &builder, initial, false)?;
    let corner0 = Point2::new(result.params[0], result.params[1]);
    let corner1 = Point2::new(result.params[2], result.params[3]);
    let radius = result.params[4] as f32;
    let point = end_intersection(input, &working.last()?.c, &mid_point(&corner0, &corner1))?;

    let c = working.take_circles();
    let edge = AfEdge::new(
        point,
        AfEdgeGeometry::RoundedSquare(corner0, corner1, radius),
    );
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

/// Fit a sharp corner trailing or leading edge to airfoil section data.
///
/// The edge geometry is a single apex point where the two airfoil surfaces meet. An open
/// two-segment polyline (from `p0` through one free corner point to `p1`) is fit to the section
/// points beyond the clipping line of the last inscribed circle. The apex is both the corner
/// and the reported edge point.
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, process the front (leading edge) end of the airfoil; when `false`,
///   process the rear (trailing edge) end.
///
/// # Returns
///
/// An [`AfEdgeFit`] whose edge geometry is [`AfEdgeGeometry::Sharp`], containing the fitted apex
/// point as both the corner and the edge location.
pub fn fit_sharp_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0;
    let p1 = working.last()?.p1;

    // To find the initial, we'll find how far the farthest point is from the end clipping line
    // and add that vector to p0 and p1
    let c = working.clip.at(working.clip_max_scalar());
    let initial = DVector::from(vec![c.x, c.y]);

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        let mut bdata = BoundaryData2::new_open(p0);
        bdata.add_seg_xy(params[0], params[1]);
        bdata.add_seg(&p1);
        bdata.try_to_boundary()
    });

    let result = fit_boundary_to_points(&working.fit_points, &builder, initial, false)?;
    let corner = Point2::new(result.params[0], result.params[1]);
    let point = end_intersection(input, &working.last()?.c, &corner)?;

    let c = working.take_circles();
    let edge = AfEdge::new(point, AfEdgeGeometry::Sharp(corner));
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

/// Fit a full-round (semicircular) trailing or leading edge to airfoil section data.
///
/// The edge is a single circular arc spanning from `p0` to `p1`, parameterized by a free circle
/// center and radius. The arc is constrained to wind in the correct direction relative to the
/// airfoil interior.
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, process the front (leading edge) end of the airfoil; when `false`,
///   process the rear (trailing edge) end.
///
/// # Returns
///
/// An [`AfEdgeFit`] whose edge geometry is [`AfEdgeGeometry::FullRound`], containing the fitted
/// circle center and radius, with the edge location at the outermost point on the camber axis.
pub fn fit_full_round_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0;
    let p1 = working.last()?.p1;
    let wind_line = Line2::new(p0, working.clip.direction);
    let wind_dir = wind_line.winding_direction(&working.clip.origin);

    let initial_r = dist(&p0, &p1) / 2.0;
    let initial_c = working.clip.at(working.clip_max_scalar() - initial_r * 0.7);
    let initial = DVector::from(vec![initial_c.x, initial_c.y, initial_r]);

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        let mut bdata = BoundaryData2::new_open(p0);
        let c = Point2::new(params[0], params[1]);
        bdata.add_full_round(&c, params[2], &p1, wind_dir)?;
        bdata.try_to_boundary()
    });

    let result = fit_boundary_to_points(&working.fit_points, &builder, initial, false)?;
    let circle = Circle2::new(result.params[0], result.params[1], result.params[2]);
    let point = end_intersection(input, &working.last()?.c, &circle)?;

    let c = working.take_circles();
    let edge = AfEdge::new(
        point,
        AfEdgeGeometry::FullRound(circle.center, circle.r() as f32),
    );
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

/// Fit a blended-round trailing or leading edge to airfoil section data.
///
/// The edge is a three-arc boundary that smoothly blends the airfoil surfaces into a central
/// leading/trailing edge circle. Starting from `p0`, the boundary consists of:
///   1. A blend arc tangent to the suction surface and internally tangent to the edge circle.
///   2. The central edge circle arc.
///   3. A blend arc tangent to the pressure surface and internally tangent to the edge circle.
///
/// The blend arcs are derived from the surface tangents at the contact points of the last
/// inscribed circle, shifted by the edge circle radius so that each blend arc is simultaneously
/// tangent to its surface tangent line and internally tangent to the edge circle.
///
/// After fitting, the inscribed circle stack is refined by inserting an additional circle
/// seeded from the fitted edge circle geometry and then dynamically refining to the original
/// camber line tolerance criteria.
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, process the front (leading edge) end of the airfoil; when `false`,
///   process the rear (trailing edge) end.
///
/// # Returns
///
/// An [`AfEdgeFit`] whose edge geometry is [`AfEdgeGeometry::BlendedRound`], containing the
/// fitted edge circle center and radius, with the edge location at the outermost point on the
/// camber axis. The returned inscribed circle stack includes one additional refined circle near
/// the edge.
pub fn fit_blended_round_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    let mut working = EdgeWork::new(input, circles, at_front)?;
    let (p0, p1) = working.last_points()?;

    // Find the winding direction, which will be the same for all arcs that get built
    let wind_line = Line2::new(p0, working.clip.direction);
    let wind_dir = wind_line.winding_direction(&working.clip.origin);

    // Prepare the initial guess
    let initial_r = dist(&p0, &p1) / 4.0;
    let initial_c = working.clip.at(working.clip_max_scalar() - initial_r * 0.7);
    let initial = DVector::from(vec![initial_c.x, initial_c.y, initial_r]);

    // Prep geometry that gets moved into the builder function
    let clip = working.clip;
    let (t0, t1) = working.last_tangents()?;

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        // First we need to actually compute the arcs
        // -----------------------------------------------------------
        let center = Point2::new(params[0], params[1]);
        let radius = params[2];
        let (arc0, arc1) = end_arcs(&t0, &t1, &clip, &center, radius);

        // Now we construct the boundary
        // -----------------------------------------------------------
        let mut bdata = BoundaryData2::new_open(p0);
        bdata.add_arc(&arc0.center, &arc0.at_end(), wind_dir.is_cw());
        bdata.add_arc(&center, &arc1.at_end(), wind_dir.is_cw());
        bdata.add_arc(&arc1.center, &p1, wind_dir.is_cw());
        bdata.try_to_boundary()
    });

    // Perform the fitting and get the best fit circle
    let result = fit_boundary_to_points(&working.fit_points, &builder, initial, false)?;
    let circle = Circle2::new(result.params[0], result.params[1], result.params[2]);

    // Now let's refine the inscribed circles to get closer to the edge circle
    refine_from_edge_circle(&mut working, input, &circle)?;
    let point = end_intersection(input, &working.last()?.c, &circle)?;

    let c = working.take_circles();
    let edge = AfEdge::new(
        point,
        AfEdgeGeometry::BlendedRound(circle.center, circle.r() as f32),
    );
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

/// Fit a cubic spline to the end of the section data and extract the point of maximum curvature.
///
/// After fitting, the inscribed circle stack is refined by inserting an additional circle
/// seeded from the fitted edge circle geometry and then dynamically refining to the original
/// camber line tolerance criteria.
///
/// # Arguments
///
/// * `input` - The section geometry and search tolerances.
/// * `circles` - The inscribed circle stack produced by the camber line fitting step.
/// * `at_front` - When `true`, process the front (leading edge) end of the airfoil; when `false`,
///   process the rear (trailing edge) end.
///
/// # Returns
///
/// An [`AfEdgeFit`] whose edge geometry is [`AfEdgeGeometry::SplineMaxK`], containing the
/// osculating edge circle center and radius, with the edge location at the outermost point on the
/// camber axis. The returned inscribed circle stack includes one additional refined circle near
/// the edge.
pub fn fit_spline_max_k(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<(AfEdgeFit, CubicSpline2)> {
    let mut working = EdgeWork::new(input, circles, at_front)?;
    let (t0, t1) = working.last_tangents()?;
    let clip = working.clip;

    let mut stack = InscribedVec::new(vec![working.last()?.clone()]);

    // Get the first attempt at a split
    let mut fittings = vec![spline_fit(&working.fit_points, &clip, &t0, &t1, &stack)?];

    for _ in 0..8 {
        // If the distance between the last circle and the current one is less than 5% of the
        // circle radius, we can stop now
        let last_circle = stack.last().unwrap();
        let last_fit = &fittings.last().unwrap();
        let inscribed_dist = dist(&last_circle.center(), &last_fit.circle.center);
        if inscribed_dist < 0.05 * last_fit.circle.r() {
            break;
        }

        // Find the halfway point and get a new inscribed circle
        let l = half_line(last_circle, &clip.direction, fittings.last().unwrap());
        stack.refine_and_push(input.try_inscribed(&l).unwrap(), input);

        let clip = stack.end_clip_line()?;
        let (t0, t1) = stack.last_tangents()?;
        let Ok(fit) = spline_fit(&working.fit_points, &clip, &t0, &t1, &stack) else {
            break;
        };

        // In order for this to be valid, the arc length of the spline must be at least 1/3 of the
        // arc length of the minimum circle
        if fit.spline.arc_length() < 2.0 * PI * fit.circle.r() / 3.0 {
            break;
        }

        let this_circle = fit.circle;
        let last_circle = last_fit.circle;
        fittings.push(fit);

        let circle_error = dist(&this_circle, &last_circle);
        let radius_error = (this_circle.r() - last_circle.r()).abs();
        if circle_error < 0.01 * this_circle.r() && radius_error < 0.01 * this_circle.r() {
            break;
        }
    }

    let fit = fittings.last().unwrap();
    refine_from_edge_circle(&mut working, input, &fit.circle)?;

    let expected_point = fit.spline.position(fit.t_max_k);
    let point = end_intersection(input, &fit.circle.center, &expected_point)?;

    let c = working.take_circles();
    let edge = AfEdge::new(
        point,
        AfEdgeGeometry::SplineMaxK(fit.circle.center, fit.circle.r() as f32),
    );
    Ok((
        AfEdgeFit::new(edge, fit.residuals.clone(), c),
        fit.spline.clone(),
    ))
}

fn half_line(inscribed: &Inscribed, clip_dir: &Vector2, last_result: &SfResult) -> Line2 {
    let l0 = Line2::new(inscribed.center(), *clip_dir);
    let l1 = last_result.end_line();
    let l = l0.slerp(&l1, 0.5).rotated(PI / 2.0);
    if l.direction.dot(&inscribed.contact_dir()) < 0.0 {
        l.reversed()
    } else {
        l
    }
}

struct SfResult {
    spline: CubicSpline2,
    t_max_k: f64,
    circle: Circle2,
    residuals: DVector,
    #[allow(dead_code)] // computed and stored but not yet consumed
    avg_residual: f64,
}

impl SfResult {
    fn end_line(&self) -> Line2 {
        Line2::new_normalize(
            self.circle.center,
            self.spline.position(self.t_max_k) - self.circle.center,
        )
    }
}

fn spline_fit(
    all_points: &[Point2],
    clip: &Line2,
    t0: &Line2,
    t1: &Line2,
    stack: &InscribedVec,
) -> Result<SfResult> {
    let clip_v = clip.normal().into_inner();
    let fit_points = all_points
        .iter()
        .filter_map(|p| {
            if clip.scalar_project(p) > 0.0 {
                Some(*p)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();

    let t0 = *t0;
    let t1 = *t1;

    // First stage simplified spline fit, constrained to existing tangent direction
    let build0: SplineBuildFn<2> = Box::new(move |params: &DVector| {
        let i0 = Iso2::from(clip_v * params[0]);
        let i1 = Iso2::from(clip_v * params[1]);

        let ts0 = i0 * t0;
        let ts1 = i1 * t1;
        let p1 = ts0.at(params[2]);
        let p2 = ts1.at(params[3]);
        Ok(CubicSpline2::new(ts0.origin, p1, p2, ts1.origin()))
    });

    let ti = fit_points
        .iter()
        .map(|p| clip.scalar_project(p))
        .fold(0.0f64, |a, b| a.max(b));

    let initial0 = DVector::from(vec![0.0, 0.0, ti, ti]);
    let result0 = fit_spline_to_points(&fit_points, &build0, initial0)?.params;

    // First stage simplified spline fit, constrained to existing tangent direction
    let build1: SplineBuildFn<2> = Box::new(move |params: &DVector| {
        let i0 = Iso2::from(clip_v * params[0]);
        let i1 = Iso2::from(clip_v * params[1]);

        let ts0 = (i0 * t0).rotated(params[4]);
        let ts1 = (i1 * t1).rotated(params[5]);
        let p1 = ts0.at(params[2]);
        let p2 = ts1.at(params[3]);
        Ok(CubicSpline2::new(ts0.origin, p1, p2, ts1.origin()))
    });
    let initial1 = DVector::from(vec![
        result0[0], result0[1], result0[2], result0[3], 0.0, 0.0,
    ]);
    let result1 = fit_spline_to_points(&fit_points, &build1, initial1)?;

    let spline = build1(&result1.params)?;

    // We need to find the curvature local maximum closest to the center
    let local_maxima = spline
        .find_curvature_maxima()
        .into_iter()
        .filter(|k| k.t > f64::EPSILON && k.t < 1.0 - f64::EPSILON)
        .collect::<Vec<_>>();
    if local_maxima.len() != 1 {
        return Err(format!(
            "Expected one curvature maximum, found {}",
            local_maxima.len()
        )
        .into());
    }
    let max_k = local_maxima[0];
    let circle = spline
        .curvature_circle(max_k.t)
        .ok_or("Failed to find curvature circle at maximum curvature")?;

    // Now we have to compute the residuals
    let query = spline.into_query();
    let back_ref = if stack.len() > 1 {
        let p0s = stack.iter().map(|c| c.p0).collect::<Vec<_>>();
        let p1s = stack.iter().map(|c| c.p1).collect::<Vec<_>>();
        Some((
            Curve2::from_points(&p0s, 1e-12, false)?,
            Curve2::from_points(&p1s, 1e-12, false)?,
        ))
    } else {
        None
    };
    let mut residuals = Vec::new();
    for p in all_points.iter() {
        if clip.scalar_project(p) > 0.0 {
            residuals.push(query.project_point(p).value)
        } else {
            if let Some((c0, c1)) = &back_ref {
                residuals.push(c0.dist_to_point(p).min(c1.dist_to_point(p)))
            } else {
                residuals.push(query.project_point(p).value)
            }
        }
    }
    let residuals = DVector::from_vec(residuals);
    let avg_residual = residuals.mean();

    Ok(SfResult {
        spline: query.into_spline(),
        t_max_k: max_k.t,
        circle,
        residuals,
        avg_residual,
    })
}

// =============================================================================================
// Common tools for the edge implementations
// =============================================================================================

/// Convenience wrapper for finding the LE/TE point given two known points on the camberline. The
/// first point `e0` should be the last inscribed circle center. The second point, `e1` should
/// come from the fitting geometry and would be where the theoretical camber end would be if the
/// fitting boundary fit exactly to the observed points.  The resulting point should be very close
/// to `e1` but will be the result of intersecting the line `e0->e1` with the section.
fn end_intersection(
    input: &SectionInput,
    e0: &impl PCoords<2>,
    e1: &impl PCoords<2>,
) -> Result<Point2> {
    let line = Line2::new_normalize(Point2::from(e0.coords()), e1.coords() - e0.coords());
    let ts = input.section.intersections_with_line(&line);
    let t = ts
        .iter()
        .filter_map(|(t, _)| if t.is_sign_positive() { Some(t) } else { None })
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .ok_or("Failed to find intersection between end line and section")?;
    Ok(line.at(*t))
}

/// Build the chord which spans the open gap in a section, running from the first point of the
/// section curve to the last. Returns an error if the section is closed, in which case there is no
/// gap to span.
fn open_gap_segment(input: &SectionInput) -> Result<Segment2> {
    if input.section.is_closed() {
        return Err("Cannot fit an open edge on a closed section".into());
    }
    Segment2::new(
        &input.section.at_front().point(),
        &input.section.at_back().point(),
    )
}

/// Build the outward facing ray at the end of the inscribed circle stack. The origin is the center
/// of the last inscribed circle, which is definitionally on the camber line, and the direction
/// comes from the end clipping line, which is already normalized and already faces away from the
/// interior of the airfoil.
fn end_camber_ray(stack: &InscribedVec) -> Result<SurfacePoint2> {
    let clip = stack.end_clip_line()?;
    let last = stack
        .last()
        .ok_or("Cannot build an end camber ray from an empty inscribed circle stack")?;
    Ok(SurfacePoint2::new_normalize(last.center(), clip.direction))
}

/// The open section counterpart to `end_intersection`. Because the section perimeter does not close
/// across the gap there is nothing there to intersect, so the end of the camber line is projected
/// onto the chord spanning the gap instead of onto the section itself.
fn end_cap_intersection(stack: &InscribedVec, cap: &Segment2) -> Result<Point2> {
    let ray = Line2::from(&end_camber_ray(stack)?);
    let (t0, t1) = ray
        .intersection_params(&cap.to_line())
        .ok_or("The end of the camber line is parallel to the open gap")?;

    // `ray` has a unit direction and `cap.to_line()` spans the gap over `t` in [0, 1], so this
    // rejects a gap which lies behind the camber end or off to one side of it.
    if t0 <= 0.0 || !(0.0..=1.0).contains(&t1) {
        return Err("The end of the camber line does not project onto the open gap".into());
    }
    Ok(ray.at(t0))
}

fn refine_from_edge_circle(
    working: &mut EdgeWork,
    input: &SectionInput,
    edge_circle: &Circle2,
) -> Result<()> {
    let (a, b) = working
        .last()?
        .c
        .outer_tangents_to(edge_circle)
        .ok_or("Failed to find outer tangents between last inscribed circle and edge circle")?;

    let fake = Inscribed::new(*edge_circle, a.b, b.b);
    let fake_camber_line = (if fake.camber_point().normal.dot(&working.clip.direction) < 0.0 {
        fake.camber_point().reversed()
    } else {
        fake.camber_point()
    })
    .shifted(-edge_circle.r() * 0.1);

    let test_line = Line2::new(fake_camber_line.point, fake.contact_dir());
    let test_line = if test_line.direction.dot(&working.last()?.contact_dir()) < 0.0 {
        test_line.reversed()
    } else {
        test_line
    };

    let inscribed = input
        .try_inscribed(&test_line)
        .ok_or("Failed to find crossing line for blended round edge refinement")?;

    working.stack.refine_and_push(inscribed, input);
    Ok(())
}

fn end_arcs(t0: &Line2, t1: &Line2, clip: &Line2, center: &Point2, radius: f64) -> (Arc2, Arc2) {
    let t0s = t0.offset_by(t0.signed_projection_dist(&clip.origin).signum() * radius);
    let t1s = t1.offset_by(t1.signed_projection_dist(&clip.origin).signum() * radius);

    // We get the blend arcs. Arc0 goes from p0 to the leading edge circle, and arc1 goes from
    // p1 to the leading edge circle. We need to keep the order of endpoints right when we
    // actuallyh build the boundary
    (
        blend_arc(&t0s, center, radius),
        blend_arc(&t1s, center, radius),
    )
}

fn blend_arc(shifted_tangent: &Line2, le_center: &Point2, le_radius: f64) -> Arc2 {
    let base_circle = Circle2::from_tangent_and_point(shifted_tangent, le_center);
    let v0 = shifted_tangent.origin - base_circle.center;
    let v1 = le_center - base_circle.center;
    let theta0 = base_circle.angle_of_point(&shifted_tangent.origin);
    let theta = signed_angle(&v0, &v1);
    Arc2::new(
        base_circle.center,
        le_radius + base_circle.r(),
        theta0,
        theta,
    )
}

struct EdgeWork {
    stack: InscribedVec,
    fit_points: Vec<Point2>,
    clip: Line2,
    at_front: bool,
}

impl EdgeWork {
    /// Get a reference to the last inscribed circle
    fn last(&self) -> Result<&Inscribed> {
        self.stack
            .last()
            .ok_or_else(|| "No inscribed circles found".into())
    }

    /// Get the contact points `p0` and `p1` of the last inscribed circle, returned in that order
    fn last_points(&self) -> Result<(Point2, Point2)> {
        let c = self.last()?;
        Ok((c.p0, c.p1))
    }

    /// Get the tangent lines at the contact points `p0` and `p1` (in that order) of the last
    /// inscribed circle
    fn last_tangents(&self) -> Result<(Line2, Line2)> {
        self.stack.last_tangents()
    }

    /// The maximum scalar projection of the fit points on the clipping line
    fn clip_max_scalar(&self) -> f64 {
        self.fit_points
            .iter()
            .map(|p| self.clip.scalar_project(p))
            .fold(0.0f64, |a, b| a.max(b))
    }

    fn new(input: &SectionInput, circles: Vec<Inscribed>, at_front: bool) -> Result<Self> {
        if circles.len() < 2 {
            return Err("Not enough inscribed circles to partition edge geometry".into());
        }
        let mut stack = InscribedVec::new(circles);
        // If we're supposed to be detecting the edge at the front of the airfoil, we'll reverse the
        // stack of circles so that we're working at the back
        if at_front {
            stack.reverse_order_only();
        }
        let clip = stack.end_clip_line()?;

        // Now we want to isolate the data beyond the end of the last circle
        let mut fit_points = Vec::new();
        for p in input.section.points() {
            if clip.scalar_project(p).is_sign_positive() {
                fit_points.push(*p);
            }
        }

        Ok(Self {
            stack,
            fit_points,
            clip,
            at_front,
        })
    }

    fn take_circles(self) -> Vec<Inscribed> {
        let mut stack = self.stack;
        if self.at_front {
            stack.reverse_order_only();
        }
        stack.take_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::fill_gaps;
    use crate::{Circle2, Curve2, Result, curve2};
    use approx::assert_relative_eq;

    macro_rules! inscribed_vec {
        ( $(($cx:expr, $cy:expr, $r:expr, $p0x:expr, $p0y:expr, $p1x:expr, $p1y:expr)),* $(,)? ) => {
            vec![
                $(
                    Inscribed::new(
                        Circle2::new($cx as f64, $cy as f64, $r as f64),
                        Point2::new($p0x as f64, $p0y as f64),
                        Point2::new($p1x as f64, $p1y as f64),
                    )
                ),*
            ]
        };
    }

    #[test]
    fn blended_round() -> Result<()> {
        #[rustfmt::skip]
        let circles = inscribed_vec!((2.2811413775, -0.0216920927, 0.9251595807, 2.5597579707, -0.9039013076, 2.6696158336, 0.8179553873), (2.2960721657, -0.0226447166, 0.9197547076, 2.5730610498, -0.8996999700, 2.6822771205, 0.8120974686));
        #[rustfmt::skip]
        let curve = curve2!((0.0400644720, -1.2493577703), (0.1994998688, -1.2339772293), (2.1250000000, -0.9921567416), (2.1595998950, -0.9871817834), (2.5077756896, -0.9203181920), (2.8791537371, -0.8030307176), (3.2023916716, -0.6572063115), (3.2337585140, -0.6419920329), (3.2591962842, -0.6275713815), (3.3117449009, -0.5909157412), (3.3591746750, -0.5478412753), (3.4007068109, -0.4990552652), (3.4356593521, -0.4453587760), (3.4634583787, -0.3876335024), (3.4836474315, -0.3268272920), (3.4958950069, -0.2639385808), (3.5000000000, -0.2000000000), (3.4958950069, -0.1360614192), (3.4836474315, -0.0731727080), (3.4634583787, -0.0123664976), (3.4567933318, 0.0033220401), (3.4356593521, 0.0453587760), (3.2843521374, 0.3034888770), (3.0543300863, 0.5621893823), (2.7763917433, 0.7685540426), (2.4622149156, 0.9139122774), (2.2845275866, 0.9586678530), (2.1595998950, 0.9871817834), (2.1250000000, 0.9921567416), (0.1994998688, 1.2339772293), (0.0400644720, 1.2493577703))?;
        let expected = Circle2::new(3.0, -0.2, 0.5);

        let input = SectionInput::new(&curve, 1e-3);
        let result = fit_blended_round_edge(&input, circles, false)?;

        assert!(matches!(
            result.edge.geometry,
            AfEdgeGeometry::BlendedRound(_, _)
        ));
        let c = match result.edge.geometry {
            AfEdgeGeometry::BlendedRound(p, r) => Circle2::new(p.x, p.y, r as f64),
            _ => unreachable!(),
        };

        assert_relative_eq!(c.center, expected.center, epsilon = 1e-2);
        assert_relative_eq!(c.r(), expected.r(), epsilon = 1e-2);
        Ok(())
    }

    #[test]
    fn open_edge() -> Result<()> {
        // A constant thickness slab with a sharp apex on the left and an open gap on the right,
        // between the first point (3, -1) and the last point (3, 1). The camber line is y = 0 and
        // every inscribed circle has a radius of 1.
        #[rustfmt::skip]
        let curve = curve2!((3.0, -1.0), (0.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (3.0, 1.0))?;
        #[rustfmt::skip]
        let circles = inscribed_vec!((1.0, 0.0, 1.0, 1.0, -1.0, 1.0, 1.0), (1.5, 0.0, 1.0, 1.5, -1.0, 1.5, 1.0));

        let input = SectionInput::new(&curve, 1e-3);
        let result = fit_open_edge(&input, circles, false)?;

        assert!(matches!(result.edge.geometry, AfEdgeGeometry::Open));

        // The section is symmetric about y = 0, so the camber line meets the gap chord at its
        // midpoint. The circle centers are only located to `resolve_tol`, which bounds how
        // precisely the end of the camber line can be aimed.
        assert_relative_eq!(
            result.edge.point,
            Point2::new(3.0, 0.0),
            epsilon = input.resolve_tol
        );

        // The camber line should have been advanced into the gap
        assert!(
            result.circles.len() > 2,
            "Expected the inscribed circle stack to grow, got {} circles",
            result.circles.len()
        );

        // Every appended circle must keep the p0/p1 ordering of the seed circles
        for c in result.circles.iter() {
            assert!(
                c.p0.y < c.p1.y,
                "Inscribed circle contact points are out of order: {:?} then {:?}",
                c.p0,
                c.p1
            );
        }

        Ok(())
    }

    #[test]
    fn open_edge_rejects_closed_section() -> Result<()> {
        #[rustfmt::skip]
        let curve = curve2!(tol: 1e-6, closed: true; (3.0, -1.0), (0.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (3.0, 1.0))?;
        #[rustfmt::skip]
        let circles = inscribed_vec!((1.0, 0.0, 1.0, 1.0, -1.0, 1.0, 1.0), (1.5, 0.0, 1.0, 1.5, -1.0, 1.5, 1.0));

        let input = SectionInput::new(&curve, 1e-3);
        assert!(fit_open_edge(&input, circles, false).is_err());
        Ok(())
    }

    #[test]
    fn sharp_corner() -> Result<()> {
        #[rustfmt::skip]
        let circles = inscribed_vec!((1.9305063262, -0.0000000382, 0.3382035549, 2.0374556693, -0.3208481102, 2.0374556923, 0.3208481026), (1.9727817705, 0.0000000367, 0.3248348872, 2.0755035923, -0.3081654692, 2.0755035702, 0.3081654766));
        #[rustfmt::skip]
        let curve = curve2!((0.0000000000, -1.0000000000), (3.0000000000, 0.0000000000), (0.0000000000, 1.0000000000))?;
        let expected = Point2::new(3.0, 0.0);

        let input = SectionInput::new(&curve, 1e-3);
        let result = fit_sharp_edge(&input, circles, false)?;

        assert!(matches!(result.edge.geometry, AfEdgeGeometry::Sharp(_)));
        let p = match result.edge.geometry {
            AfEdgeGeometry::Sharp(p) => p,
            _ => unreachable!(),
        };

        assert_relative_eq!(p, expected, epsilon = 1e-4);
        Ok(())
    }

    #[test]
    fn rounded_square() -> Result<()> {
        #[rustfmt::skip]
        let circles = inscribed_vec!((2.1599430104, 0.0000000260, 0.8848039683, 2.3095192376, -0.8720694131, 2.3095192290, 0.8720694145), (2.2705435064, -0.0000000254, 0.8661069410, 2.4169589903, -0.8536414845, 2.4169589987, 0.8536414830));
        #[rustfmt::skip]
        let curve = curve2!((0.0400644720, -1.2493577703), (0.1994998688, -1.2339772293), (3.0398999738, -0.7467954459), (3.0711318967, -0.7396669633), (3.1011958358, -0.7286031558), (3.1295981421, -0.7137856908), (3.1558724505, -0.6954578706), (3.1795873375, -0.6739206377), (3.2003534055, -0.6495276326), (3.2178296760, -0.6226793880), (3.2317291893, -0.5938167512), (3.2418237158, -0.5634136460), (3.2479475035, -0.5319692904), (3.2500000000, -0.5000000000), (3.2500000000, 0.5000000000), (3.2479475035, 0.5319692904), (3.2418237158, 0.5634136460), (3.2317291893, 0.5938167512), (3.2178296760, 0.6226793880), (3.2003534055, 0.6495276326), (3.1795873375, 0.6739206377), (3.1558724505, 0.6954578706), (3.1295981421, 0.7137856908), (3.1011958358, 0.7286031558), (3.0711318967, 0.7396669633), (3.0398999738, 0.7467954459), (0.1994998688, 1.2339772293), (0.0400644720, 1.2493577703))?;
        let exp0 = Point2::new(3.2500000000, -0.7107675827);
        let exp1 = Point2::new(3.2500000000, 0.7107675827);

        let input = SectionInput::new(&curve, 1e-3);
        let result = fit_rounded_square_edge(&input, circles, false)?;

        assert!(matches!(
            result.edge.geometry,
            AfEdgeGeometry::RoundedSquare(_, _, _)
        ));
        let (c0, c1, r) = match result.edge.geometry {
            AfEdgeGeometry::RoundedSquare(c0, c1, r) => (c0, c1, r),
            _ => unreachable!(),
        };

        assert_relative_eq!(c0, exp0, epsilon = 1e-4);
        assert_relative_eq!(c1, exp1, epsilon = 1e-4);
        assert_relative_eq!(r, 0.25, epsilon = 1e-4);
        Ok(())
    }

    #[test]
    fn square_end() -> Result<()> {
        #[rustfmt::skip]
        let circles = inscribed_vec!((1.0008511537, -0.0000000228, 0.8681375436, 1.1085303607, -0.8614337049, 1.1085303664, 0.8614337042), (1.1093683467, 0.0000000253, 0.8546776428, 1.2153780644, -0.8480777420, 1.2153780582, 0.8480777427));
        #[rustfmt::skip]
        let curve = curve2!((0.0000000000, -1.0000000000), (2.0000000000, -0.7500000000), (2.0000000000, 0.7500000000), (0.0000000000, 1.0000000000))?;
        let curve = Curve2::from_points(&fill_gaps(curve.points(), 0.1), 1e-6, false)?;

        let input = SectionInput::new(&curve, 1e-3);
        let result = fit_square_edge(&input, circles, false)?;

        assert!(matches!(result.edge.geometry, AfEdgeGeometry::Square(_, _)));
        let (corner0, corner1) = match result.edge.geometry {
            AfEdgeGeometry::Square(c0, c1) => (c0, c1),
            _ => unreachable!(),
        };
        assert_relative_eq!(corner0, Point2::new(2.0, -0.75), epsilon = 1e-6);
        assert_relative_eq!(corner1, Point2::new(2.0, 0.75), epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn full_round() -> Result<()> {
        #[rustfmt::skip]
        let circles = inscribed_vec!((1.4945495152, 0.0000000315, 1.0625288507, 1.6272841220, -1.0542054591, 1.6272841142, 1.0542054601), (1.6273656215, -0.0000000309, 1.0459370259, 1.7580275152, -1.0377436084, 1.7580275228, 1.0377436075));
        #[rustfmt::skip]
        let curve = curve2!((0.0400644720, -1.2493577703), (0.1994998688, -1.2339772293), (2.1595998950, -0.9871817834), (2.2845275866, -0.9586678530), (2.4047833431, -0.9144126230), (2.5183925683, -0.8551427630), (2.6234898019, -0.7818314825), (2.7183493501, -0.6956825506), (2.8014136219, -0.5981105305), (2.8713187041, -0.4907175520), (2.9269167573, -0.3752670049), (2.9672948630, -0.2536545839), (2.9917900138, -0.1278771617), (3.0000000000, 0.0000000000), (2.9917900138, 0.1278771617), (2.9672948630, 0.2536545839), (2.9269167573, 0.3752670049), (2.8713187041, 0.4907175520), (2.8014136219, 0.5981105305), (2.7183493501, 0.6956825506), (2.6234898019, 0.7818314825), (2.5183925683, 0.8551427630), (2.4047833431, 0.9144126230), (2.2845275866, 0.9586678530), (2.1595998950, 0.9871817834), (0.1994998688, 1.2339772293), (0.0400644720, 1.2493577703))?;

        let input = SectionInput::new(&curve, 1e-3);
        let result = fit_full_round_edge(&input, circles, false)?;

        assert!(matches!(
            result.edge.geometry,
            AfEdgeGeometry::FullRound(_, _)
        ));
        let (c, r) = match result.edge.geometry {
            AfEdgeGeometry::FullRound(x, y) => (x, y as f64),
            _ => unreachable!(),
        };
        assert_relative_eq!(r, 1.0, epsilon = 1e-6);
        assert_relative_eq!(c, Point2::new(2.0, 0.0), epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn spline_max_k() -> Result<()> {
        #[rustfmt::skip]
        let circles = inscribed_vec!((2.1997529412, -0.3448105433, 0.3280168196, 2.3049108998, -0.6555143985, 2.3925708564, -0.0794498685), (2.2103311462, -0.3464202883, 0.3231008514, 2.3139130936, -0.6524675985, 2.4002593391, -0.0850365191));
        let expected = Circle2::new(2.4966211850, -0.3899148214, 0.1416283396);

        let spline = CubicSpline2::new(
            Point2::new(0.0, -1.0),
            Point2::new(3.0, -0.7),
            Point2::new(4.0, -0.5),
            Point2::new(0.0, 1.0),
        );
        let as_points = spline.polyline(1e-4);
        let curve = Curve2::from_points(&as_points, 1e-6, false)?;
        let input = SectionInput::new(&curve, 1e-3);
        let (result, _) = fit_spline_max_k(&input, circles, false)?;

        assert!(matches!(
            result.edge.geometry,
            AfEdgeGeometry::SplineMaxK(_, _)
        ));
        let (c, r) = match result.edge.geometry {
            AfEdgeGeometry::SplineMaxK(c, r) => (c, r as f64),
            _ => unreachable!(),
        };
        assert_relative_eq!(r, expected.r(), epsilon = 1e-3);
        assert_relative_eq!(c, expected.center, epsilon = 1e-3);
        Ok(())
    }
}
