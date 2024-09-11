//! Helper functions for implementing airfoil geometry algorithms.

use crate::airfoil::InscribedCircle;
use crate::common::points::dist;
use crate::common::Intersection;
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

/// This function will attempt to extract a smaller curve from the larger airfoil section curve
/// which has the data just for the end of the airfoil. To do so, it will attempt to segment the
/// section into two portions, with the cuts happening at the start and end of the spanning ray of
/// the inscribed circle.  The sub-curve which has less than 25% of the total section perimeter
/// will be returned, or `None` if no such sub-curve can be found.
///
/// # Arguments
///
/// * `section`: The airfoil section to be segmented
/// * `last_station`: The last inscribed circle in the provided airfoil section
///
/// returns: Option<Curve2>
pub fn extract_edge_sub_curve(section: &Curve2, last_station: &InscribedCircle) -> Option<Curve2> {
    let split0 = section.at_closest_to_point(&last_station.spanning_ray.origin());
    let split1 = section.at_closest_to_point(
        &(last_station.spanning_ray.origin() + last_station.spanning_ray.dir()),
    );

    let portion0 = section.between_lengths(split0.length_along(), split1.length_along());
    let portion1 = section.between_lengths(split1.length_along(), split0.length_along());

    let perimeter = section.length();

    if let Some(p0) = portion0 {
        if p0.length() < perimeter * 0.25 {
            return Some(p0);
        }
    }

    if let Some(p1) = portion1 {
        if p1.length() < perimeter * 0.25 {
            return Some(p1);
        }
    }

    None
}

/// This struct provides a convenient way to store a list of inscribed circles after they've been
/// oriented and edit the end of the collection depending on the orientation.
pub struct OrientedCircles {
    circles: Vec<InscribedCircle>,
    reversed: bool,
}

impl OrientedCircles {
    pub fn create(reversed: bool) -> Self {
        OrientedCircles {
            circles: Vec::new(),
            reversed,
        }
    }

    pub fn new(circles: Vec<InscribedCircle>, reversed: bool) -> Self {
        OrientedCircles { circles, reversed }
    }

    /// Get the full camber curve from the inscribed circles in this container. The curve will be
    /// oriented such that the first point is away from the working edge (based on the `reversed`
    /// flag) and the end of the curve points towards the working edge.
    pub fn get_full_curve(&self) -> Result<Curve2> {
        let curve = curve_from_inscribed_circles(&self.circles, 1e-4)?;
        if self.reversed {
            Ok(curve.reversed())
        } else {
            Ok(curve)
        }
    }

    /// Discard any inscribed circles beyond the specified point.
    pub fn discard_sections_beyond_point(&mut self, p: &Point2) {
        let c = self.get_full_curve().unwrap();

        let l = c.at_closest_to_point(p).length_along();

        self.circles.retain(|cr| {
            let d = c.at_closest_to_point(&cr.center()).length_along();
            d < l
        });
    }

    pub fn get_end_curve(&self, distance: f64) -> Result<Curve2> {
        let mut total = 0.0;
        let mut points = Vec::new();
        if self.circles.len() < 2 {
            return Err(Box::from(
                "Cannot create a curve from less than two circles",
            ));
        }

        let mut i = if self.reversed {
            0
        } else {
            self.circles.len() - 1
        };

        while total < distance {
            let c = &self.circles[i];
            if let Some(last) = points.last() {
                total += dist(last, &c.center());
            }

            points.push(c.circle.center);

            if self.reversed {
                if i == self.circles.len() - 1 {
                    break;
                }
                i += 1;
            } else {
                if i == 0 {
                    break;
                }

                i -= 1;
            }
        }

        let points = points.into_iter().rev().collect::<Vec<_>>();
        Curve2::from_points(&points, 1e-4, false)
    }

    /// Get the inscribed circle at the end of the collection.
    pub fn last(&self) -> Option<&InscribedCircle> {
        if self.reversed {
            self.circles.first()
        } else {
            self.circles.last()
        }
    }

    /// Add a new inscribed circle to the end of the collection.  The circle will be reversed if
    /// it is not oriented in the same direction as the last circle in the collection.
    pub fn push(&mut self, circle: InscribedCircle) {
        let mut c = circle;

        if let Some(last) = self.last() {
            if last.spanning_ray.dir().dot(&c.spanning_ray.dir()) < 0.0 {
                c.reverse_in_place();
            }
        }

        if self.reversed {
            self.circles.insert(0, c);
        } else {
            self.circles.push(c);
        }
    }

    /// Using a camber curve generated at the end of the inscribed circles, this function will
    /// find the intersection point between a straight projection of the camber curve end and the
    /// airfoil section boundary.
    ///
    /// # Arguments
    ///
    /// * `section`: the airfoil section curve associated with the inscribed circles
    ///
    /// returns: OPoint<f64, Const<2>>
    pub fn intersect_from_end(&self, section: &Curve2) -> Point2 {
        let c = self.get_end_curve(self.last().unwrap().radius()).unwrap();
        let end = c.at_back().direction_point();
        let ts = section.intersection(&end);
        let t = ts.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        end.at_distance(*t)
    }

    /// Take the inscribed circles out of the container, leaving an empty container behind. The
    /// container cannot be used after this operation.
    pub fn take_circles(self) -> Vec<InscribedCircle> {
        self.circles
    }
}

/// Refines a stack of inscribed circles by checking the interpolation error between the circles
/// and adding new circles between them when the error is above a certain tolerance.
///
/// # Arguments
///
/// * `section`: the airfoil section curve
/// * `dest`: the destination vector which will receive the refined inscribed circles
/// * `stack`: the stack of inscribed circles to refine
/// * `tol`: the tolerance value which will determine when to add new circles between the existing
/// circles
///
/// returns: ()
pub fn refine_stations(
    section: &Curve2,
    dest: &mut OrientedCircles,
    stack: &mut Vec<InscribedCircle>,
    outer_tol: f64,
    inner_tol: f64,
) {
    while let Some(next) = stack.pop() {
        if let Some(last) = dest.last() {
            let n = if next.spanning_ray.dir().dot(&last.spanning_ray.dir()) < 0.0 {
                next.reversed()
            } else {
                next
            };

            let test_ray = n.spanning_ray.symmetry(&last.spanning_ray);

            if let Some(ray) = section.try_create_spanning_ray(&test_ray) {
                let mid = inscribed_from_spanning_ray(section, &ray, inner_tol);
                let error = mid.interpolation_error(&n, last);

                // TODO: check the distance between the centers to make sure we're not stuck
                if error > outer_tol {
                    // We are out of tolerance, we need to put next back on the stack and then put
                    // the mid-ray on top of it and try again
                    stack.push(n);
                    stack.push(mid);
                } else {
                    // We are within tolerance, we can put the next station in the destination. We
                    // will keep the mid-station since we've already gone through the trouble of
                    // creating it.
                    dest.push(mid);
                    dest.push(n);
                }
            }
        } else {
            dest.push(next);
        }
    }
}
