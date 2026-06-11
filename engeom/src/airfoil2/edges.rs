use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::airfoil2::{AfEdge, AfEdgeGeometry, SectionInput};
use crate::geom2::{BndBuildFn, BoundaryData2, BoundaryEditor, fit_boundary_to_points};
use crate::{DVector, Point2, Result};

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
) -> Result<(AfEdge, AfEdgeGeometry, Vec<Inscribed>)> {
    let (mut stack, fit_points) = edge_prep(input, circles, at_front)?;

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

    let c = edge_finalize(stack, at_front);
    todo!()
}
