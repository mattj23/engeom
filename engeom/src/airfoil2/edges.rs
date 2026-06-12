use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::airfoil2::{AfEdge, AfEdgeGeometry, SectionInput};
use crate::common::{Intersection, dist, mid_point};
use crate::geom2::{
    BndBuildFn, BoundaryData2, BoundaryEditor, BoundaryElement2, LineOps2, fit_boundary_to_points,
    signed_angle,
};
use crate::{Arc2, Circle2, DVector, Line2, Point2, Result};

#[derive(Clone, Debug)]
pub struct AfEdgeFit {
    pub edge: AfEdge,
    pub residuals: DVector,
    pub circles: Vec<Inscribed>,
}

impl AfEdgeFit {
    pub fn new(result: AfEdge, residuals: DVector, circles: Vec<Inscribed>) -> Self {
        Self {
            edge: result,
            residuals,
            circles,
        }
    }
}

pub fn square_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0.clone();
    let p1 = working.last()?.p1.clone();

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
    let point = mid_point(&corner0, &corner1);

    let c = working.take_circles();
    let edge = AfEdge::new(point, AfEdgeGeometry::Square(corner0, corner1));
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

pub fn rounded_square_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0.clone();
    let p1 = working.last()?.p1.clone();

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
    let point = mid_point(&corner0, &corner1);

    let c = working.take_circles();
    let edge = AfEdge::new(
        point,
        AfEdgeGeometry::RoundedSquare(corner0, corner1, radius),
    );
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

pub fn sharp_corner_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0.clone();
    let p1 = working.last()?.p1.clone();

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
    let c = working.take_circles();
    let edge = AfEdge::new(corner, AfEdgeGeometry::Sharp(corner));
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

pub fn full_round_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0.clone();
    let p1 = working.last()?.p1.clone();
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

    // Find the end point as an intersection with the circle
    let end_line = Line2::new(circle.center, circle.center - working.last()?.center());
    let ts = circle
        .intersection(&end_line)
        .iter()
        .filter(|t| t.is_sign_positive())
        .cloned()
        .collect::<Vec<_>>();
    let t = ts
        .first()
        .ok_or("Failed to find intersection between round and end line")?;
    let point = end_line.at(*t);

    let c = working.take_circles();
    let edge = AfEdge::new(
        point,
        AfEdgeGeometry::FullRound(circle.center, circle.r() as f32),
    );
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

pub fn blended_round_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let (p0, p1) = working.last_points()?;

    // Find the winding direction, which will be the same for all arcs that get built
    let wind_line = Line2::new(p0, working.clip.direction);
    let wind_dir = wind_line.winding_direction(&working.clip.origin);

    // Prepare the initial guess
    let initial_r = dist(&p0, &p1) / 4.0;
    let initial_c = working.clip.at(working.clip_max_scalar() - initial_r * 0.7);
    let initial = DVector::from(vec![initial_c.x, initial_c.y, initial_r]);

    // Prep geometry that gets moved into the builder function
    let clip = working.clip.clone();
    let (t0, t1) = working.last_tangents()?;

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        // First we need to actually compute the arcs
        // -----------------------------------------------------------
        let center = Point2::new(params[0], params[1]);
        let radius = params[2];

        // We shift the tangencies in towards the center by the radius value
        let t0s = t0.new_parallel(t0.signed_projection_dist(&clip.origin).signum() * radius);
        let t1s = t1.new_parallel(t1.signed_projection_dist(&clip.origin).signum() * radius);

        // We get the blend arcs. Arc0 goes from p0 to the leading edge circle, and arc1 goes from
        // p1 to the leading edge circle. We need to keep the order of endpoints right when we
        // actuallyh build the boundary
        let arc0 = blend_arc(&t0s, &center, radius);
        let arc1 = blend_arc(&t1s, &center, radius);

        // Now we construct the boundary
        // -----------------------------------------------------------
        let mut bdata = BoundaryData2::new_open(p0);
        bdata.add_arc(&arc0.center(), &arc0.at_end(), wind_dir.is_cw());
        bdata.add_arc(&center, &arc1.at_end(), wind_dir.is_cw());
        bdata.add_arc(&arc1.center(), &p1, wind_dir.is_cw());
        bdata.try_to_boundary()
    });

    let result = fit_boundary_to_points(&working.fit_points, &builder, initial, false)?;
    let circle = Circle2::new(result.params[0], result.params[1], result.params[2]);

    // Find the end point as an intersection with the circle
    let end_line = Line2::new(circle.center, circle.center - working.last()?.center());
    let ts = circle
        .intersection(&end_line)
        .iter()
        .filter(|t| t.is_sign_positive())
        .cloned()
        .collect::<Vec<_>>();
    let t = ts
        .first()
        .ok_or("Failed to find intersection between round and end line")?;
    let point = end_line.at(*t);

    let c = working.take_circles();
    let edge = AfEdge::new(
        point,
        AfEdgeGeometry::BlendedRound(circle.center, circle.r() as f32),
    );
    Ok(AfEdgeFit::new(edge, result.residuals, c))
}

fn blend_arc(shifted_tangent: &Line2, le_center: &Point2, le_radius: f64) -> Arc2 {
    let base_circle = Circle2::new_tangent_and_point(shifted_tangent, le_center);
    let v0 = shifted_tangent.origin - base_circle.center;
    let v1 = le_center - base_circle.center;
    let theta0 = base_circle.angle_of_point(&shifted_tangent.origin);
    let theta = signed_angle(&v0, &v1);
    Arc2::circle_angles(
        base_circle.center,
        le_radius + base_circle.r(),
        theta0,
        theta,
    )
}

// =============================================================================================
// Common tools for the edge implementations
// =============================================================================================
struct EdgeWork<'a> {
    input: &'a SectionInput<'a>,
    stack: InscribedVec,
    fit_points: Vec<Point2>,
    clip: Line2,
    at_front: bool,
}

impl<'a> EdgeWork<'a> {
    /// Get a reference to the last inscribed circle
    fn last(&self) -> Result<&Inscribed> {
        self.stack
            .last()
            .ok_or_else(|| "No inscribed circles found".into())
    }

    /// Get the contact points `p0` and `p1` of the last inscribed circle, returned in that order
    fn last_points(&self) -> Result<(Point2, Point2)> {
        let c = self.last()?;
        Ok((c.p0.clone(), c.p1.clone()))
    }

    /// Get the tangent lines at the contact points `p0` and `p1` (in that order) of the last
    /// inscribed circle
    fn last_tangents(&self) -> Result<(Line2, Line2)> {
        let c = self.last()?.c.clone();
        let (p0, p1) = self.last_points()?;
        let t0 = c.at_closest_to_point(&p0).direction_line();
        let t1 = c.at_closest_to_point(&p1).direction_line();

        Ok((
            if t0.direction.dot(&self.clip.direction) < 0.0 {
                t0.new_reversed()
            } else {
                t0
            },
            if t1.direction.dot(&self.clip.direction) < 0.0 {
                t1.new_reversed()
            } else {
                t1
            },
        ))
    }

    /// The maximum scalar projection of the fit points on the clipping line
    fn clip_max_scalar(&self) -> f64 {
        let d_max = self
            .fit_points
            .iter()
            .map(|p| self.clip.scalar_project(p))
            .fold(0.0f64, |a, b| a.max(b));
        d_max
    }

    fn new(input: &'a SectionInput<'a>, circles: Vec<Inscribed>, at_front: bool) -> Result<Self> {
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
            input,
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
    use crate::{Circle2, Curve2, Line2, Result, curve2};
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
        let result = blended_round_edge(&input, circles, false)?;

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
    fn sharp_corner() -> Result<()> {
        #[rustfmt::skip]
        let circles = inscribed_vec!((1.9305063262, -0.0000000382, 0.3382035549, 2.0374556693, -0.3208481102, 2.0374556923, 0.3208481026), (1.9727817705, 0.0000000367, 0.3248348872, 2.0755035923, -0.3081654692, 2.0755035702, 0.3081654766));
        #[rustfmt::skip]
        let curve = curve2!((0.0000000000, -1.0000000000), (3.0000000000, 0.0000000000), (0.0000000000, 1.0000000000))?;
        let expected = Point2::new(3.0, 0.0);

        let input = SectionInput::new(&curve, 1e-3);
        let result = sharp_corner_edge(&input, circles, false)?;

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
        let result = rounded_square_edge(&input, circles, false)?;

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
        let result = square_edge(&input, circles, false)?;

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
        let result = full_round_edge(&input, circles, false)?;

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
}
