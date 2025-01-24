use crate::common::points::dist;
use crate::common::{signed_compliment_2pi, BestFit, Intersection};
use crate::geom2::line2::Segment2;
use crate::geom2::{directed_angle, signed_angle, Iso2, Line2, Point2, Vector2};
use crate::geom3::Vector3;
use crate::stats::{compute_mean, compute_st_dev};
use crate::AngleDir::{Ccw, Cw};
use crate::AngleInterval;
use crate::Result;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use parry2d_f64::na::{Dyn, Matrix, Owned, Vector, U1, U3};
use parry2d_f64::shape::Ball;
use serde::{Deserialize, Serialize};
use std::f64::consts::FRAC_PI_2;

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Circle2 {
    pub center: Point2,
    pub ball: Ball,
}

impl Circle2 {
    /// Create a new circle from the x and y coordinates of its center and its radius
    ///
    /// # Arguments
    ///
    /// * `x`: the x coordinate of the circle's center
    /// * `y`: the y coordinate of the circle's center
    /// * `r`: the radius of the circle
    ///
    /// returns: Circle2
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Circle2, Point2};
    /// let circle = Circle2::new(1.0, 2.0, 3.0);
    ///
    /// assert_eq!(circle.x(), 1.0);
    /// assert_eq!(circle.y(), 2.0);
    /// assert_eq!(circle.r(), 3.0);
    /// ```
    pub fn new(x: f64, y: f64, r: f64) -> Circle2 {
        Circle2 {
            center: Point2::new(x, y),
            ball: Ball::new(r),
        }
    }

    /// Create a new circle from a center point and a radius
    ///
    /// # Arguments
    ///
    /// * `center`: the center point of the circle
    /// * `r`: the radius of the circle
    ///
    /// returns: Circle2
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Circle2, Point2};
    /// let circle = Circle2::from_point(Point2::new(1.0, 2.0), 3.0);
    ///
    /// assert_eq!(circle.x(), 1.0);
    /// assert_eq!(circle.y(), 2.0);
    /// assert_eq!(circle.r(), 3.0);
    /// ```
    pub fn from_point(center: Point2, r: f64) -> Circle2 {
        Circle2 {
            center,
            ball: Ball::new(r),
        }
    }

    /// Attempt to create a fitting circle from the given points and an initial guess. The fitting
    /// is an unconstrained Levenberg-Marquardt minimization of the sum of squared errors between
    /// the points and the boundary of the circle..
    ///
    /// The initial guess is used to provide an initial estimate of the circle's center and radius,
    /// for best results this should at least be in the general vicinity of the test points.
    ///
    /// The mode parameter controls the fitting algorithm. The `BestFit::All` mode will weight all
    /// points equally, while the `BestFit::Gaussian(sigma)` mode will assign zero weights to
    /// points beyond `sigma` standard deviations from the mean.
    ///
    /// # Arguments
    ///
    /// * `points`: the points to be fit to the circle
    /// * `guess`: an initial guess for the circle's center and radius
    /// * `mode`: the fitting mode to use
    ///
    /// returns: Result<Circle2, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Circle2, Point2};
    /// use engeom::common::BestFit::All;
    /// use approx::assert_relative_eq;
    ///
    /// let points = vec![
    ///     Point2::new(-1.0, 0.0),
    ///     Point2::new(0.0, 1.0),
    ///     Point2::new(1.0, 0.0),
    ///     Point2::new(0.0, -1.0),
    /// ];
    ///
    /// let guess = Circle2::new(-1.0, 1.0, 0.1);
    /// let circle = Circle2::fitting_circle(&points, &guess, All).unwrap();
    /// assert_relative_eq!(circle.x(), 0.0);
    /// assert_relative_eq!(circle.y(), 0.0);
    /// assert_relative_eq!(circle.r(), 1.0);
    /// ```
    pub fn fitting_circle(points: &[Point2], guess: &Circle2, mode: BestFit) -> Result<Circle2> {
        fit_circle(points, &guess, mode)
    }

    /// Attempt to create a fitting circle from three points. Will return an `Err` if the points
    /// are collinear.
    ///
    /// # Arguments
    ///
    /// * `p0`: the first point
    /// * `p1`: the second point
    /// * `p2`: the third point
    ///
    /// returns: Result<Circle2, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Circle2, Point2};
    /// use approx::assert_relative_eq;
    ///
    /// let p0 = Point2::new(-1.0, 0.0);
    /// let p1 = Point2::new(0.0, 1.0);
    /// let p2 = Point2::new(1.0, 0.0);
    ///
    /// let circle = Circle2::from_3_points(p0, p1, p2).unwrap();
    /// assert_relative_eq!(circle.x(), 0.0);
    /// assert_relative_eq!(circle.y(), 0.0);
    /// assert_relative_eq!(circle.r(), 1.0);
    ///
    /// ```
    pub fn from_3_points(p0: Point2, p1: Point2, p2: Point2) -> Result<Circle2> {
        let temp = p1.x.powi(2) + p1.y.powi(2);
        let bc = (p0.x.powi(2) + p0.y.powi(2) - temp) / 2.0;
        let cd = (temp - p2.x.powi(2) - p2.y.powi(2)) / 2.0;
        let det = (p0.x - p1.x) * (p1.y - p2.y) - (p1.x - p2.x) * (p0.y - p1.y);

        if det.abs() < 1.0e-6 {
            Err("Points are collinear".into())
        } else {
            let cx = (bc * (p1.y - p2.y) - cd * (p0.y - p1.y)) / det;
            let cy = ((p0.x - p1.x) * cd - (p1.x - p2.x) * bc) / det;

            let radius = ((cx - p0.x).powi(2) + (cy - p0.y).powi(2)).sqrt();
            Ok(Circle2 {
                center: Point2::new(cx, cy),
                ball: Ball::new(radius),
            })
        }
    }

    /// Returns the x coordinate of the circle's center
    pub fn x(&self) -> f64 {
        self.center.x
    }

    /// Returns the y coordinate of the circle's center
    pub fn y(&self) -> f64 {
        self.center.y
    }

    /// Returns the radius of the circle
    pub fn r(&self) -> f64 {
        self.ball.radius
    }

    /// Returns the point on the circle's perimeter at the given angle (referenced as a
    /// counter-clockwise angle where zero is the x-axis), measured in radians
    ///
    /// # Arguments
    ///
    /// * `angle`: the CCW angle from the x-axis in radians
    ///
    /// returns: OPoint<f64, Const<2>>
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::FRAC_PI_2;
    /// use engeom::{Circle2, Point2};
    /// use approx::assert_relative_eq;
    ///
    /// let c = Circle2::new(0.0, 0.0, 1.0);
    /// let p = c.point_at_angle(FRAC_PI_2);
    ///
    /// assert_relative_eq!(p.x, 0.0);
    /// assert_relative_eq!(p.y, 1.0);
    /// ```
    pub fn point_at_angle(&self, angle: f64) -> Point2 {
        let v = Vector2::new(self.ball.radius, 0.0);
        let t = Iso2::rotation(angle);
        self.center + (t * v)
    }

    /// Given a point, project it onto the perimeter of the circle. If the point is within 1.0e-10
    /// of the center, then `None` is returned.
    ///
    /// # Arguments
    ///
    /// * `point`:
    ///
    /// returns: Option<OPoint<f64, Const<2>>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn project_point_to_perimeter(&self, point: &Point2) -> Option<Point2> {
        let v = point - self.center;
        if v.norm() < 1.0e-10 {
            None
        } else {
            let n = v.normalize();
            Some(self.center + (n * self.ball.radius))
        }
    }

    /// Determines the angle of the given point relative to the circle's center, measured in
    /// radians where zero is the x-axis and positive angles are counter-clockwise.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to measure the angle of
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    /// use std::f64::consts::FRAC_PI_2;
    /// use engeom::{Circle2, Point2};
    /// use approx::assert_relative_eq;
    ///
    /// let c = Circle2::new(0.0, 0.0, 1.0);
    /// let p = Point2::new(0.0, 1.0);
    /// let a = c.angle_of_point(&p);
    ///
    /// assert_relative_eq!(a, FRAC_PI_2);
    /// ```
    pub fn angle_of_point(&self, point: &Point2) -> f64 {
        let v = point - self.center;
        v.y.atan2(v.x)
    }

    pub fn intersection_interval(&self, other: Circle2) -> Option<AngleInterval> {
        let ints = self.intersections_with(&other);
        if ints.is_empty() {
            return None;
        }
        if ints.len() == 1 {
            let start = self.angle_of_point(&ints[0]);
            return Some(AngleInterval::new(start, 0.0));
        }

        // There are two possible intervals, one going from the first point to the second in the
        // clockwise direction, and one going via its signed compliment in the counter-clockwise
        // direction.  The valid interval is the one which also contains the center of the other
        // circle.
        let v0 = ints[0] - self.center;
        let v1 = ints[1] - self.center;
        let a = signed_angle(&v0, &v1);
        let ac = signed_compliment_2pi(a);
        let s = self.angle_of_point(&ints[0]);

        let i0 = AngleInterval::new(s, a);
        let i1 = AngleInterval::new(s, ac);

        if i0.contains(self.angle_of_point(&other.center)) {
            Some(i0)
        } else {
            Some(i1)
        }
    }

    pub fn intersections_with(&self, other: &Circle2) -> Vec<Point2> {
        const TOL: f64 = 1.0e-10;

        let mut result = Vec::new();

        let d = dist(&self.center, &other.center);
        if d < TOL {
            // Circles are concentric
            return result;
        }

        let r_sum = self.ball.radius + other.ball.radius;
        if d > r_sum {
            // Circles are too far apart
            return result;
        }

        let v = (other.center - self.center).normalize();
        let a = (self.ball.radius.powi(2) - other.ball.radius.powi(2) + d.powi(2)) / (2.0 * d);
        let p2 = self.center + (v * a);

        if (d - r_sum).abs() < TOL {
            // Circles are touching
            result.push(p2);
            return result;
        }

        let h = (self.ball.radius.powi(2) - a.powi(2)).sqrt();
        let n = Iso2::rotation(FRAC_PI_2) * v;
        result.push(p2 + (n * h));
        result.push(p2 - (n * h));
        result
    }

    /// Create a full arc of the circle, starting at zero and extending for 2π radians.
    pub fn to_arc(&self) -> Arc2 {
        Arc2 {
            circle: *self,
            angle0: 0.0,
            angle: 2.0 * std::f64::consts::PI,
        }
    }

    /// Create a partial arc of the circle, starting at `angle0` and extending for `angle` radians.
    ///
    /// # Arguments
    ///
    /// * `angle0`:
    /// * `angle`:
    ///
    /// returns: Arc2
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn to_partial_arc(&self, angle0: f64, angle: f64) -> Arc2 {
        Arc2 {
            circle: *self,
            angle0,
            angle,
        }
    }

    /// Computes the distance from the test point to the outer perimeter of the circle. If the
    /// point lies within the circle boundary the distance will be negative.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to check against the circle
    ///
    /// returns: f64
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Circle2, Point2};
    /// let c = Circle2::new(0.0, 0.0, 1.0);
    ///
    /// let d0 = c.distance_to(&Point2::new(0.0, 0.0));
    /// assert_eq!(d0, -1.0);
    ///
    /// let d1 = c.distance_to(&Point2::new(1.0, 0.0));
    /// assert_eq!(d1, 0.0);
    ///
    /// let d2 = c.distance_to(&Point2::new(2.0, 0.0));
    /// assert_eq!(d2, 1.0);
    /// ```
    pub fn distance_to(&self, point: &Point2) -> f64 {
        dist(&self.center, point) - self.ball.radius
    }
}

/// Computes the intersection parameters of a line and a circle. The intersection parameters are
/// the values of the line's parameter at the intersection points. There will be zero, one, or two
/// intersection parameters.  The intersection points can be recovered from the parameters by
/// calling `line.at(t)`.
///
/// # Arguments
///
/// * `line`:
/// * `circle`:
///
/// returns: Vec<f64, Global>
///
/// # Examples
///
/// ```
///
/// ```
pub fn intersection_line_circle(line: &dyn Line2, circle: &Circle2) -> Vec<f64> {
    // Get the parameter of the circle center onto the line
    let tc = line.projected_parameter(&circle.center);

    let d = dist(&circle.center, &line.projected_point(&circle.center));

    if (d - circle.ball.radius).abs() < 1.0e-10 {
        // If the distance from the center to the line is equal to the radius, then there is one
        // single intersection point at the tangency point
        vec![tc]
    } else if d > circle.ball.radius {
        // If the distance from the center to the line is greater than the radius, then there
        // are no intersection points
        Vec::new()
    } else {
        // There are two intersection points. The distance from the tangency point to the
        // intersection points is the height of a right triangle with hypotenuse `r` and base
        // `d`.
        let h = (circle.ball.radius.powi(2) - d.powi(2)).sqrt();

        // We can't simply add and subtract `h` from `tc` because the line may not be
        // normalized, so we need to scale the value `h` by the norm of the line's direction
        // vector.
        let th = h / line.dir().norm();

        vec![tc - th, tc + th]
    }
}

impl Intersection<&Segment2, Vec<Point2>> for Circle2 {
    fn intersection(&self, other: &Segment2) -> Vec<Point2> {
        let ts = intersection_line_circle(other, self);
        ts.iter()
            .filter_map(|&t| {
                if t >= -1.0e-10 && t <= 1.0 + 1.0e-10 {
                    Some(other.at(t))
                } else {
                    None
                }
            })
            .collect()
    }
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Arc2 {
    pub circle: Circle2,
    pub angle0: f64,
    pub angle: f64,
}

impl Arc2 {
    pub fn circle_angles(center: Point2, radius: f64, angle0: f64, angle: f64) -> Self {
        Self {
            circle: Circle2::from_point(center, radius),
            angle0,
            angle,
        }
    }

    pub fn circle_point_angle(center: Point2, radius: f64, point: Point2, angle: f64) -> Self {
        let circle = Circle2::from_point(center, radius);
        let angle0 = circle.angle_of_point(&point);
        Self {
            circle,
            angle0,
            angle,
        }
    }

    pub fn three_points(p0: Point2, p1: Point2, p2: Point2) -> Self {
        let circle = Circle2::from_3_points(p0, p1, p2).unwrap();
        let angle0 = circle.angle_of_point(&p0);
        let v0 = p0 - circle.center;
        let v2 = p2 - circle.center;

        let det = (p1.x - p0.x) * (p1.y + p0.y)
            + (p2.x - p1.x) * (p2.y + p1.y)
            + (p0.x - p2.x) * (p0.y + p2.y);
        let angle = if det < 0.0 {
            directed_angle(&v0, &v2, Ccw)
        } else {
            -directed_angle(&v0, &v2, Cw)
        };

        Self {
            circle,
            angle0,
            angle,
        }
    }

    pub fn length(&self) -> f64 {
        self.circle.ball.radius * self.angle.abs()
    }

    pub fn center(&self) -> Point2 {
        self.circle.center
    }

    pub fn radius(&self) -> f64 {
        self.circle.ball.radius
    }

    pub fn point_at_angle(&self, angle: f64) -> Point2 {
        self.circle.point_at_angle(self.angle0 + angle)
    }

    pub fn point_at_fraction(&self, fraction: f64) -> Point2 {
        self.point_at_angle(self.angle * fraction)
    }

    pub fn point_at_length(&self, length: f64) -> Point2 {
        self.point_at_fraction(length / self.length())
    }
}

type Residuals = Matrix<f64, Dyn, U1, Owned<f64, Dyn, U1>>;

///
///
/// # Arguments
///
/// * `points`:
/// * `initial`:
/// * `mode`:
///
/// returns: Result<Circle2, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn fit_circle(points: &[Point2], initial: &Circle2, mode: BestFit) -> Result<Circle2> {
    let problem = CircleFit::new(points, mode, initial);
    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        Ok(result.circle)
    } else {
        let text = format!("Failed to fit circle: {:?}", report.termination);
        Err(text.into())
    }
}

struct CircleFit<'a> {
    /// The points to be fit to the circle.
    points: &'a [Point2],

    /// The best fitting mode
    mode: BestFit,

    /// The parameters being fit
    x: Vector3,

    /// The current active circle
    circle: Circle2,

    /// The active base residuals
    base_residuals: Residuals,

    /// The active weights
    weights: Residuals,
}

impl<'a> CircleFit<'a> {
    fn new(points: &'a [Point2], mode: BestFit, initial: &Circle2) -> Self {
        let x = Vector3::new(initial.center.x, initial.center.y, initial.r());
        let circle = *initial;

        // Compute the residuals
        let mut base_residuals = Residuals::zeros(points.len());
        compute_residuals_mut(points, &circle, &mut base_residuals);

        // Compute the weights
        let mut weights = Residuals::zeros(points.len());
        compute_weights_mut(&base_residuals, &mut weights, mode);

        Self {
            points,
            mode,
            x,
            circle,
            base_residuals,
            weights,
        }
    }
}

fn compute_residuals_mut(points: &[Point2], circle: &Circle2, residuals: &mut Residuals) {
    for (i, p) in points.iter().enumerate() {
        residuals[i] = circle.distance_to(p)
    }
}

fn compute_weights_mut(residuals: &Residuals, weights: &mut Residuals, mode: BestFit) {
    match mode {
        BestFit::All => {
            weights.fill(1.0);
        }
        BestFit::Gaussian(sigma) => {
            let mean = compute_mean(residuals.as_slice()).expect("Empty slice");
            let std = compute_st_dev(residuals.as_slice()).expect("Empty slice");

            for (i, r) in residuals.iter().enumerate() {
                // How many standard deviations are we from the mean?
                let d = (r - mean).abs() / std;
                if d > sigma {
                    weights[i] = 0.0;
                } else {
                    weights[i] = 1.0;
                }
            }
        }
    }
}

impl<'a> LeastSquaresProblem<f64, Dyn, U3> for CircleFit<'a> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U3>;
    type ParameterStorage = Owned<f64, U3>;

    fn set_params(&mut self, x: &Vector<f64, U3, Self::ParameterStorage>) {
        self.x = *x;
        self.circle = Circle2::new(x[0], x[1], x[2]);
        compute_residuals_mut(self.points, &self.circle, &mut self.base_residuals);
        compute_weights_mut(&self.base_residuals, &mut self.weights, self.mode);
    }

    fn params(&self) -> Vector<f64, U3, Self::ParameterStorage> {
        self.x
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Residuals::zeros(self.points.len());
        for i in 0..self.points.len() {
            res[i] = self.base_residuals[i] * self.weights[i];
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U3, Self::JacobianStorage>> {
        let mut jac = Matrix::<f64, Dyn, U3, Self::JacobianStorage>::zeros(self.points.len());
        for (i, p) in self.points.iter().enumerate() {
            // Find the vector from the center of the circle to the point
            let v = p - self.circle.center;

            // Normalize it
            let n = v.normalize();

            // Fill in the jacobian for this row
            jac[(i, 0)] = -n.x * self.weights[i];
            jac[(i, 1)] = -n.y * self.weights[i];
            jac[(i, 2)] = -1.0 * self.weights[i];
        }

        Some(jac)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom2::{Line2, Ray2};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;
    use test_case::test_case;

    #[test_case((0.0, 0.0), (1.0, 0.0), vec![-1.0, 1.0])]
    #[test_case((0.0, 0.0), (2.0, 0.0), vec![-0.5, 0.5])]
    #[test_case((1.0, 0.0), (1.0, 0.0), vec![-2.0, 0.0])]
    #[test_case((0.0, 1.0), (1.0, 0.0), vec![0.0])]
    #[test_case((0.0, 1.5), (1.0, 0.0), Vec::<f64>::new())]
    fn simple_line_intersection(s: (f64, f64), d: (f64, f64), e: Vec<f64>) {
        let c = Circle2::new(0.0, 0.0, 1.0);
        let l = Ray2::new(Point2::new(s.0, s.1), Vector2::new(d.0, d.1));

        // let result = c.intersection(&l);
        let result = intersection_line_circle(&l, &c);
        assert_eq!(result.len(), e.len());
        for (r, e) in result.iter().zip(e.iter()) {
            assert_relative_eq!(*e, *r, epsilon = 1.0e-10);
        }
    }

    #[test_case((0.0, 0.0, 1.0), 0.0, (1.0, 0.0))]
    #[test_case((0.0, 0.0, 1.0), 90.0, (0.0, 1.0))]
    #[test_case((0.0, 0.0, 1.0), 180.0, (-1.0, 0.0))]
    #[test_case((0.0, 0.0, 1.0), 360.0, (1.0, 0.0))]
    #[test_case((1.0, 1.0, 1.0), 0.0, (2.0, 1.0))]
    #[test_case((1.0, 1.0, 1.0), 90.0, (1.0, 2.0))]
    #[test_case((1.0, 1.0, 1.0), 180.0, (0.0, 1.0))]
    #[test_case((1.0, 1.0, 1.0), 360.0, (2.0, 1.0))]
    fn test_circle_point(c: (f64, f64, f64), a: f64, r: (f64, f64)) {
        let circle = Circle2::new(c.0, c.1, c.2);
        let point = circle.point_at_angle(a * std::f64::consts::PI / 180.0);
        assert_relative_eq!(r.0, point.x, epsilon = 1.0e-10);
        assert_relative_eq!(r.1, point.y, epsilon = 1.0e-10);
    }

    #[test]
    fn test_intersection_simple() {
        let c0 = Circle2::new(0.0, 0.0, 1.0);
        let c1 = Circle2::new(1.0, 0.0, 1.0);
        let result = c0.intersections_with(&c1);
        assert_eq!(result.len(), 2);
        assert_relative_eq!(result[0].x, 0.5, epsilon = 1.0e-10);
        assert_relative_eq!(result[0].y, 0.8660254037844386, epsilon = 1.0e-10);
        assert_relative_eq!(result[1].x, 0.5, epsilon = 1.0e-10);
        assert_relative_eq!(result[1].y, -0.8660254037844386, epsilon = 1.0e-10);
    }

    #[test]
    fn three_point_arc_ccw() {
        let p0 = Point2::new(1.0, 0.0);
        let p1 = Point2::new(0.0, 1.0);
        let p2 = Point2::new(0.0, -1.0);
        let arc = Arc2::three_points(p0, p1, p2);

        assert_relative_eq!(Point2::origin(), arc.center());
        assert_relative_eq!(1.0, arc.radius());
        assert_relative_eq!(0.0, arc.angle0);
        assert_relative_eq!(3.0 * PI / 2.0, arc.angle);
    }

    #[test]
    fn three_point_arc_cw() {
        let p2 = Point2::new(1.0, 0.0);
        let p1 = Point2::new(0.0, 1.0);
        let p0 = Point2::new(0.0, -1.0);
        let arc = Arc2::three_points(p0, p1, p2);

        assert_relative_eq!(Point2::origin(), arc.center());
        assert_relative_eq!(1.0, arc.radius());
        assert_relative_eq!(-PI / 2.0, arc.angle0);
        assert_relative_eq!(-3.0 * PI / 2.0, arc.angle);
    }
}
