use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::airfoil2::{AfEdge, AfEdgeGeometry, SectionInput};
use crate::common::mid_point;
use crate::geom2::{BndBuildFn, BoundaryData2, BoundaryEditor, fit_boundary_to_points};
use crate::{DVector, Point2, Result};

#[derive(Clone, Debug)]
pub struct AfEdgeFit {
    pub result: AfEdge,
    pub residuals: Vec<f64>,
    pub circles: Vec<Inscribed>,
}

pub fn edge_prep(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<(InscribedVec, Vec<Point2>)> {
    if circles.len() < 2 {
        return Err("Not enough inscribed circles to partition edge geometry".into());
    }

    let mut stack = InscribedVec::new(circles);
    // If we're supposed to be detecting the edge at the front of the airfoil, we'll reverse the
    // stack of circles so that we're working at the back
    if at_front {
        stack.reverse_order_only();
    }
    let test_line = stack.end_clip_line()?;

    // Now we want to isolate the data beyond the end of the last circle
    let mut points = Vec::new();
    for p in input.section.points() {
        if test_line.scalar_project(p).is_sign_positive() {
            points.push(*p);
        }
    }

    Ok((stack, points))
}

pub fn edge_finalize(mut circles: InscribedVec, at_front: bool) -> Vec<Inscribed> {
    if at_front {
        circles.reverse_order_only();
    }

    circles.take_vec()
}

pub fn square_edge(
    input: &SectionInput,
    circles: Vec<Inscribed>,
    at_front: bool,
) -> Result<(AfEdge, Vec<Inscribed>)> {
    // TODO: Can we refine the stack of inscribed circles?
    let (stack, fit_points) = edge_prep(input, circles, at_front)?;

    let last = stack.last().unwrap();
    let p0 = last.p0.clone();
    let p1 = last.p1.clone();

    // To find the initial, we'll find how far the farthest point is from the end clipping line
    // and add that vector to p0 and p1
    let l = stack.end_clip_line()?;
    let d_max = fit_points
        .iter()
        .map(|p| l.scalar_project(p))
        .fold(0.0f64, |a, b| a.max(b));
    let i0 = p0 + l.direction * d_max;
    let i1 = p1 + l.direction * d_max;
    let initial = DVector::from(vec![i0.x, i0.y, i1.x, i1.y]);

    let builder: BndBuildFn = Box::new(move |params: &DVector| {
        let mut bdata = BoundaryData2::new_open(p0);
        bdata.add_seg_xy(params[0], params[1]);
        bdata.add_seg_xy(params[2], params[3]);
        bdata.add_seg(&p1);
        bdata.try_to_boundary()
    });

    let result = fit_boundary_to_points(&fit_points, &builder, initial, false)?;
    let corner0 = Point2::new(result[0], result[1]);
    let corner1 = Point2::new(result[2], result[3]);
    let point = mid_point(&corner0, &corner1);

    let c = edge_finalize(stack, at_front);

    Ok((
        AfEdge::new(point, AfEdgeGeometry::Square(corner0, corner1)),
        c,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::fill_gaps;
    use crate::{Circle2, Curve2, Line2, Result};
    use approx::assert_relative_eq;

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
        let (edge, _circles) = square_edge(&input, c, false)?;

        assert!(matches!(edge.geometry, AfEdgeGeometry::Square(_, _)));
        let (corner0, corner1) = match edge.geometry {
            AfEdgeGeometry::Square(c0, c1) => (c0, c1),
            _ => unreachable!(),
        };
        assert_relative_eq!(corner0, c0, epsilon = 1e-6);
        assert_relative_eq!(corner1, c1, epsilon = 1e-6);
        Ok(())
    }
}
