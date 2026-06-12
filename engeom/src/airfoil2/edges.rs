use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::airfoil2::{AfEdge, AfEdgeGeometry, SectionInput};
use crate::common::{Intersection, dist, mid_point};
use crate::geom2::{BndBuildFn, BoundaryData2, BoundaryEditor, LineOps2, fit_boundary_to_points};
use crate::{Circle2, DVector, Line2, Point2, Result};

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
    fn last(&self) -> Result<&Inscribed> {
        self.stack
            .last()
            .ok_or_else(|| "No inscribed circles found".into())
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
