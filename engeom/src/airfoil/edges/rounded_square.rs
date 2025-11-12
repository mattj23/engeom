use crate::common::PCoords;
use crate::common::points::{dist, mid_point};
use crate::geom2::{BoundaryElement, Line2, Segment2};
use crate::{Arc2, Curve2, Point2, Result};

pub fn best_fit_rounded_square(edge_curve: &Curve2, te_intr: &Point2) -> Result<()> {
    let root_seg = Segment2::try_new(&edge_curve.at_front(), &edge_curve.at_back())?;
    let root_center = mid_point(&root_seg.a, &root_seg.b);

    let c0 = te_intr + (root_seg.a - &root_center);
    let c1 = te_intr + (root_seg.b - &root_center);
    let r0 = dist(&root_seg.a, &root_seg.b) / 4.0;

    todo!()
}

struct RoundedSquareEdge {
    seg0: Segment2,
    arc0: Arc2,
    seg1: Segment2,
    arc1: Arc2,
    seg2: Segment2,
}

impl RoundedSquareEdge {
    pub fn distance_to(&self, point: &impl PCoords<2>) -> f64 {
        let mut best = dist_to(&self.seg0, point);
        check_update(&mut best, dist_to(&self.arc0, point));
        check_update(&mut best, dist_to(&self.seg1, point));
        check_update(&mut best, dist_to(&self.arc1, point));
        check_update(&mut best, dist_to(&self.seg2, point));
        best
    }
}

fn check_update(best: &mut f64, candidate: f64) {
    if candidate.abs() < best.abs() {
        *best = candidate;
    }
}

fn dist_to(element: &impl BoundaryElement, p: &impl PCoords<2>) -> f64 {
    let closest = element.closest_to_point(p);
    let d = dist(&closest, p);
    if closest.normal_scalar_projection(p) < 0.0 {
        -d
    } else {
        d
    }
}
