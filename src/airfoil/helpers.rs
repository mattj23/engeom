//! Helper functions for implementing airfoil geometry algorithms.

use crate::airfoil::InscribedCircle;
use crate::common::points::dist;
use crate::geom2::polyline2::SpanningRay;
use crate::geom2::Line2;
use crate::{Circle2, Curve2, Point2, Result};

/// Reverse the order of the inscribed circles in the container. This will additionally
/// reverse the order of the points on the circles so that the upper and lower points are
/// swapped.
///
/// # Arguments
///
/// * `stations`: the inscribed circles in the airfoil section
///
/// returns: ()
pub fn reverse_inscribed_circles(stations: &mut [InscribedCircle]) {
    stations.reverse();
    stations.iter_mut().for_each(|i| i.reverse_in_place());
}

/// Takes a container of inscribed circles and generates a curve from the circle centers. The
/// curve tolerance (distance at which points are considered identical) must be specified
/// externally, but typically for camber curves a small value can be used.  The curve will be
/// oriented in the same order as the inscribed circles, which may or may not be the correct
/// orientation based on whether the circles have been previously ordered.
///
/// # Arguments
///
/// * `stations`: The list of inscribed circles to generate the curve from.
/// * `tol`: The tolerance value for the curve itself (the distance at which two points on the
/// curve are considered equal), this can be a small value for camber curves (~1e-4 for inches.
/// 1e-3 for millimeters) without causing problems.
///
/// returns: Result<Curve2, Box<dyn Error, Global>>
pub fn curve_from_inscribed_circles(stations: &[InscribedCircle], tol: f64) -> Result<Curve2> {
    let camber_points = stations.iter().map(|c| c.circle.center).collect::<Vec<_>>();
    Curve2::from_points(&camber_points, tol, false)
}

/// Find the inscribed circle with the largest diameter in the container of inscribed circles. Will
/// return a reference to the circle with the largest diameter, or None if the container is empty.
///
/// # Arguments
///
/// * `stations`: a slice of inscribed circles to search for the largest diameter
///
/// returns: Option<&InscribedCircle>
pub fn find_tmax_circle(stations: &[InscribedCircle]) -> Option<&InscribedCircle> {
    let mut max_diameter = 0.0;
    let mut max_circle = None;
    for c in stations {
        let diameter = c.circle.ball.radius * 2.0;
        if diameter > max_diameter {
            max_diameter = diameter;
            max_circle = Some(c);
        }
    }

    max_circle
}

/// Calculates the position and radius of an inscribed circle based on a spanning ray and its
/// curve. The inscribed circle center will be located on the ray, somewhere between 0 and the ray
/// length, and situated such that the circle is tangent to the curve at two points in opposite
/// directions of the ray.
///
/// This is found by evaluating the distance from points on the spanning ray to points on the curve,
/// looking for the point on the ray which is furthest from any point on the curve. The distance
/// from any point on the ray to the closest point on the curve will never be farther than the
/// distance from the point to the ray ends, but it may be *closer* when the local area of the
/// curve is not orthogonal to the ray.
///
/// To efficiently find the point of maximum distance from the section curve, this algorithm uses
/// a maximum distance binary search looking at the distance from points along the ray to their
/// nearest points on the curve and adjusting the search bounds accordingly.  It makes the
/// assumption that we are far enough from the leading and trailing edges that there are no local
/// maxima along the ray.
///
/// # Arguments
///
/// * `curve`: the airfoil section curve which the inscribed circle is being found for
/// * `ray`: a spanning ray on the curve which the inscribed circle center will be coincident with
/// * `tol`: a tolerance value which will terminate the search when the distance between the
/// limit in one direction and the limit in the other direction is less than this value. This is
/// typically the inner_tol value, or the airfoil parameters tolerance * 1e-2
///
/// returns: InscribedCircle
pub fn inscribed_from_spanning_ray(curve: &Curve2, ray: &SpanningRay, tol: f64) -> InscribedCircle {
    // Here, positive and negative refer to the directions of the limits.  The positive direction
    // starts at the ray's full length, while the negative direction starts at its origin.
    let mut positive = InscribedCircleSearchState::new(1.0, ray.at(1.0));
    let mut negative = InscribedCircleSearchState::new(0.0, ray.at(0.0));

    // `working` is a point on the ray which will be updated during the search, and used to find the
    // distance to the curve.
    let mut working;

    // While the distance between the positive and negative search bounds is greater than the
    // tolerance, continue to search for the inscribed circle center.
    while (positive.fraction - negative.fraction) * ray.dir().norm() > tol {
        // We will update the working point to be right in the middle of the positive and negative
        // direction limits.
        let fraction = (positive.fraction + negative.fraction) * 0.5;
        working = ray.at(fraction);

        // Now we find the closest position on the curve to the working point, and calculate the
        // distance and direction to that point. The direction will be used to determine which side
        // of the limits we will adjust.
        let closest = curve.at_closest_to_point(&working);
        let to_closest = closest.point() - working; // The direction vector to the closest point
        let distance = dist(&working, &closest.point());

        // If the direction vector to the closest point is in the positive direction of the ray,
        // then we will adjust the positive limit.  Otherwise, we will adjust the negative limit.
        if to_closest.dot(&ray.dir()) > 0.0 {
            positive.update(fraction, distance, closest.point());
        } else {
            negative.update(fraction, distance, closest.point());
        }
    }

    // Finally, we will put the center of the inscribed circle at the midpoint of the positive and
    // negative limits, splitting the difference one last time, and we will set the radius to be
    // the average of the positive and negative distances. By this point the difference will be
    // below the tolerance value.
    let circle = Circle2::from_point(
        ray.at((positive.fraction + negative.fraction) * 0.5),
        (positive.distance + negative.distance) * 0.5,
    );

    InscribedCircle::new(ray.clone(), positive.point, negative.point, circle)
}

/// A struct representing one side of the binary search state for the inscribed circle.
struct InscribedCircleSearchState {
    /// The fraction of the spanning ray length beyond which we know the inscribed circle center
    /// is not located.  This value will start at 0.0 for the low side and 1.0 for the high side,
    /// and will be incrementally adjusted until it converges somewhere near the middle.
    fraction: f64,

    /// The distance to the closest point on the curve from the point on the ray at the specified
    /// fraction of its length.
    distance: f64,

    /// The point on the curve which is closest to the point on the ray at the specified fraction
    point: Point2,
}

impl InscribedCircleSearchState {
    fn new(fraction: f64, point: Point2) -> InscribedCircleSearchState {
        InscribedCircleSearchState {
            fraction,
            distance: 0.0,
            point,
        }
    }

    fn update(&mut self, fraction: f64, distance: f64, point: Point2) {
        self.distance = distance;
        self.fraction = fraction;
        self.point = point;
    }
}