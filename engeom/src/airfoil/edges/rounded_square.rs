use crate::common::PCoords;
use crate::common::points::{dist, mid_point};
use crate::geom2::{BoundaryElement, Segment2};
use crate::na::{Dyn, Matrix, Owned, U1, U6, Vector, Vector6};
use crate::{Arc2, Circle2, Curve2, Point2, Result};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use serde::{Deserialize, Serialize};

pub fn best_fit_rounded_square(edge_curve: &Curve2, te_intr: &Point2) -> Result<(Arc2, Arc2, f64)> {
    let root_seg = Segment2::try_new(&edge_curve.at_front(), &edge_curve.at_back())?;
    let root_center = mid_point(&root_seg.a, &root_seg.b);

    let c0 = te_intr + (root_seg.a - root_center);
    let c1 = te_intr + (root_seg.b - root_center);
    let r0 = dist(&root_seg.a, &root_seg.b) / 4.0;

    // Make default params
    let default_x = [c0.x, c0.y, c1.x, c1.y, r0, r0];

    let points = edge_curve.clone_points();
    let fit = RoundedSquareEdgeFit::new(&points, default_x)?;
    let (result, report) = LevenbergMarquardt::new().minimize(fit);

    if !report.termination.was_successful() {
        return Err("Rounded square edge fitting did not converge successfully".into());
    }

    let edge = RoundedSquareEdge::from_params(&root_seg.a, &root_seg.b, &result.params())?;
    let max_residual = *result
        .residuals()
        .unwrap()
        .map(|x| x.abs())
        .iter()
        .max_by(|a, b| a.total_cmp(b))
        .unwrap();

    Ok((edge.arc0, edge.arc1, max_residual))
}

struct RoundedSquareEdgeFit<'a> {
    points: &'a [Point2],
    start: Point2,
    end: Point2,
    x: Vector6<f64>,
}

impl<'a> RoundedSquareEdgeFit<'a> {
    pub fn new(points: &'a [Point2], initial_x: [f64; 6]) -> Result<Self> {
        let start = points[0];
        let end = points[points.len() - 1];
        let x = Vector6::<f64>::from_column_slice(&initial_x);

        Ok(RoundedSquareEdgeFit {
            points,
            start,
            end,
            x,
        })
    }
}

impl LeastSquaresProblem<f64, Dyn, U6> for RoundedSquareEdgeFit<'_> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.x = *x;
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.x
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let current = RoundedSquareEdge::from_params(&self.start, &self.end, &self.x).ok()?;
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());

        for (i, p) in self.points.iter().enumerate() {
            res[i] = current.distance_to(p);
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let current = RoundedSquareEdge::from_params(&self.start, &self.end, &self.x).ok()?;
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());
        let delta = 1e-6;

        let mut deltas = Vec::with_capacity(6);
        for i in 0..6 {
            let mut delta_x = self.x;
            delta_x[i] += delta;
            let perturbed =
                RoundedSquareEdge::from_params(&self.start, &self.end, &delta_x).ok()?;
            deltas.push(perturbed);
        }

        for (i, p) in self.points.iter().enumerate() {
            let base_dist = current.distance_to(p);
            for j in 0..6 {
                let perturbed_dist = deltas[j].distance_to(p);
                jac[(i, j)] = (perturbed_dist - base_dist) / delta;
            }
        }

        Some(jac)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RoundedSquareEdge {
    pub seg0: Segment2,
    pub arc0: Arc2,
    pub seg1: Segment2,
    pub arc1: Arc2,
    pub seg2: Segment2,
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

    pub fn from_params(start: &Point2, end: &Point2, x: &Vector6<f64>) -> Result<Self> {
        let corner0 = Point2::new(x[0], x[1]);
        let corner1 = Point2::new(x[2], x[3]);
        let r0 = x[4];
        let r1 = x[5];

        Self::build(start, &corner0, &corner1, end, r0, r1)
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
        .project_point_to_perimeter(corner)
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
