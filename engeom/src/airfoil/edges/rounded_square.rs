use crate::common::PCoords;
use crate::common::points::{dist, mid_point};
use crate::geom2::{BoundaryElement, Segment2};
use crate::{Arc2, Circle2, Curve2, Point2, Result, Vector2};
use serde::{Deserialize, Serialize};
use std::path::Path;

pub fn best_fit_rounded_square(edge_curve: &Curve2, te_intr: &Point2) -> Result<()> {
    let root_seg = Segment2::try_new(&edge_curve.at_front(), &edge_curve.at_back())?;
    let root_center = mid_point(&root_seg.a, &root_seg.b);

    let c0 = te_intr + (root_seg.a - &root_center);
    let c1 = te_intr + (root_seg.b - &root_center);
    let r0 = dist(&root_seg.a, &root_seg.b) / 4.0;

    todo!()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TempOutput {
    points: Vec<Point2>,
    edge: RoundedSquareEdge,
    corner0: Point2,
    corner1: Point2,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
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

    pub fn build(
        start: &Point2,
        corner0: &Point2,
        corner1: &Point2,
        end: &Point2,
        r0: f64,
        r1: f64,
    ) -> Result<Self> {
        // First create the square end and the two circles
        let square = SquareEnd::new(start, corner0, corner1, end)?;
        let arc0 = arc_from_circle(&square.seg0, &square.seg1, r0)?;
        let arc1 = arc_from_circle(&square.seg1, &square.seg2, r1)?;
        let seg0 = Segment2::try_new(start, &arc0.a())?;
        let seg1 = Segment2::try_new(&arc0.b(), &arc1.a())?;
        let seg2 = Segment2::try_new(&arc1.b(), end)?;

        Ok(RoundedSquareEdge {
            seg0,
            arc0,
            seg1,
            arc1,
            seg2,
        })
    }
}

/// Try to create the arc connecting two corner segments with a given radius. The arc will start
/// on `seg0` and end on `seg1`, being tangent to both segments.  The endpoint of `seg0` and the
/// startpoint of `seg1` are assumed to be the corner point, and should be the same.
///
/// # Arguments
///
/// * `seg0`:
/// * `seg1`:
/// * `radius`:
///
/// returns: Result<Arc2, Box<dyn Error, Global>>
fn arc_from_circle(seg0: &Segment2, seg1: &Segment2, radius: f64) -> Result<Arc2> {
    let corner = &seg0.b;
    let v0 = seg0.a - corner;
    let v1 = seg1.b - corner;
    let circle = Circle2::tangent_to_corner(corner, &v0, &v1, radius)?;

    let m0 = seg0.closest_to_point(&circle.center);
    let m1 = seg1.closest_to_point(&circle.center);
    let mc = circle
        .project_point_to_perimeter(&corner)
        .ok_or("Failed to project corner point to circle perimeter")?;

    Ok(Arc2::three_points(m0.point, mc, m1.point))
}

struct SquareEnd {
    /// Goes from start to corner0
    seg0: Segment2,

    /// Goes from corner0 to corner1
    seg1: Segment2,

    /// Goes from corner1 to end
    seg2: Segment2,
}

impl SquareEnd {
    pub fn corner0(&self) -> &Point2 {
        &self.seg0.b
    }

    pub fn corner1(&self) -> &Point2 {
        &self.seg1.b
    }

    pub fn start(&self) -> &Point2 {
        &self.seg0.a
    }

    pub fn end(&self) -> &Point2 {
        &self.seg2.b
    }

    pub fn corner0_v0(&self) -> Vector2 {
        self.start() - self.corner0()
    }

    pub fn corner0_v1(&self) -> Vector2 {
        self.corner1() - self.corner0()
    }

    pub fn corner1_v0(&self) -> Vector2 {
        self.corner0() - self.corner1()
    }

    pub fn corner1_v1(&self) -> Vector2 {
        self.end() - self.corner1()
    }

    pub fn try_circle0(&self, r0: f64) -> Result<Circle2> {
        Circle2::tangent_to_corner(self.corner0(), &self.corner0_v0(), &self.corner0_v1(), r0)
    }

    pub fn try_circle1(&self, r1: f64) -> Result<Circle2> {
        Circle2::tangent_to_corner(self.corner1(), &self.corner1_v0(), &self.corner1_v1(), r1)
    }

    pub fn new(start: &Point2, corner0: &Point2, corner1: &Point2, end: &Point2) -> Result<Self> {
        let seg0 = Segment2::try_new(start, corner0)?;
        let seg1 = Segment2::try_new(corner0, corner1)?;
        let seg2 = Segment2::try_new(corner1, end)?;

        Ok(SquareEnd { seg0, seg1, seg2 })
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
