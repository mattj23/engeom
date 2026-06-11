use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::airfoil2::{AfEdge, AfEdgeGeometry, SectionInput};
use crate::common::{Intersection, dist, mid_point};
use crate::geom2::{BndBuildFn, BoundaryData2, BoundaryEditor, LineOps2, fit_boundary_to_points};
use crate::{Circle2, DVector, Line2, Point2, Result};
use num_traits::real::Real;

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

pub fn full_round_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<AfEdgeFit> {
    // TODO: Can we refine the stack of inscribed circles?
    let working = EdgeWork::new(input, circles, at_front)?;
    let p0 = working.last()?.p0.clone();
    let p1 = working.last()?.p1.clone();

    let initial_r = dist(&p0, &p1) / 2.0;
    let initial_c = working.clip.at(working.clip_max_scalar() - initial_r);
    let initial = DVector::from(vec![initial_c.x, initial_c.y, initial_r]);

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        let mut bdata = BoundaryData2::new_open(p0);
        let c = Point2::new(params[0], params[1]);
        bdata.add_full_round(&c, params[2], &p1)?;
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
    use crate::{Circle2, Curve2, Line2, Result};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    fn core_circle() -> Vec<Inscribed> {
        let c0 = Circle2::new(0.0, 0.0, 1.25);
        let c1 = Circle2::new(2.0, 0.0, 1.0);
        let (s0, s1) = c0.outer_tangents_to(&c1).unwrap();
        vec![
            Inscribed::new(c0, s1.a, s0.a),
            Inscribed::new(c1, s1.b, s0.b),
        ]
    }

    #[test]
    fn square_end() -> Result<()> {
        let c = core_circle();
        let l0 = Line2::from_points(&c[0].p0, &c[1].p0);
        let l1 = Line2::from_points(&c[0].p1, &c[1].p1);
        let c0 = l0.at(1.5);
        let c1 = l1.at(1.5);
        let points = vec![c[0].p0, c0, c1, c[0].p1];
        let points = fill_gaps(&points, 0.1);
        let curve = Curve2::from_points(&points, 1e-6, false)?;
        let input = SectionInput::new(&curve, 1e-3);
        let result = square_edge(&input, c, false)?;

        assert!(matches!(result.edge.geometry, AfEdgeGeometry::Square(_, _)));
        let (corner0, corner1) = match result.edge.geometry {
            AfEdgeGeometry::Square(c0, c1) => (c0, c1),
            _ => unreachable!(),
        };
        assert_relative_eq!(corner0, c0, epsilon = 1e-6);
        assert_relative_eq!(corner1, c1, epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn full_round() -> Result<()> {
        let c = core_circle();
        let l0 = Line2::from_points(&c[0].p0, &c[1].p0);
        let l1 = Line2::from_points(&c[0].p1, &c[1].p1);

        let corner = l0.intersect(&l1).unwrap();
        let mut bdata = BoundaryData2::new_open(c[0].p0);
        bdata.add_seg(&c[1].p0);
        bdata.add_corner_fillets(&vec![corner, c[1].p1], 0.75)?;
        bdata.add_seg(&c[0].p1);
        let boundary = bdata.try_to_boundary()?;

        let expected =
            Circle2::tangent_to_corner(&corner, &(c[1].p0 - corner), &(c[1].p1 - corner), 0.75)?;

        // Quick check that the expected circle is on the boundary
        for a in vec![-PI / 4.0, 0.0, PI / 4.0] {
            let test_point = expected.point_at_angle(a);
            assert_relative_eq!(
                dist(&boundary.at_closest_to_point(&test_point).1, &test_point),
                0.0,
                epsilon = 1e-6
            );
        }
        let points = boundary.to_points(0.01)?;
        let curve = Curve2::from_points(&points, 1e-6, false)?;
        let input = SectionInput::new(&curve, 1e-3);
        let result = full_round_edge(&input, c, false)?;

        assert!(matches!(
            result.edge.geometry,
            AfEdgeGeometry::FullRound(_, _)
        ));
        let (c, r) = match result.edge.geometry {
            AfEdgeGeometry::FullRound(x, y) => (x, y as f64),
            _ => unreachable!(),
        };
        assert_relative_eq!(c, expected.center, epsilon = 1e-6);
        assert_relative_eq!(r, expected.r(), epsilon = 1e-6);
        Ok(())
    }
}
