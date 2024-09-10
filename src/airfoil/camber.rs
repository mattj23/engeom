//! This module contains the implementation of airfoil camber line detection algorithms.

use super::helpers::reverse_inscribed_circles;
use crate::airfoil::inscribed_circle::InscribedCircle;
use crate::common::points::{dist, mid_point};
use crate::geom2::hull::farthest_pair_indices;
use crate::geom2::polyline2::SpanningRay;
use crate::geom2::{rot90, Line2};
use crate::AngleDir::{Ccw, Cw};
use crate::Result;
use crate::{Circle2, Curve2, Point2};
use parry2d_f64::query::Ray;
use parry2d_f64::shape::ConvexPolygon;

/// Attempts to calculate and extract the mean camber line from an airfoil section curve and its
/// convex hull using the inscribed circle method.
///
/// Functionally, this algorithm will work by first finding the farthest pair of points on the
/// convex hull, then bisecting the section at the halfway point of that longest span to create the
/// initial station. From there, the algorithm will advance along the camber direction, first in
/// one direction and then in the other, to extract the camber line. At each inscribed circle,
/// it will attempt to move forward by a fraction of the last circle's radius, incrementally
/// reducing the forward motion down to a minimum threshold when it encounters problems. As it
/// advances it looks forward for the farthest point in the section in the ever-updating camber
/// direction, and will terminate the search when the distance to that farthest point beyond the
/// last inscribed circle is less than a minimum threshold of the last circle's radius.
///
/// As it advances, it will look at the interpolation error between the inscribed circles and add
/// new circles between them when the error is above a certain tolerance, refining the camber line
/// to within the specified tolerance.
///
/// # Arguments
///
/// * `section`: the airfoil section curve
/// * `hull`: the convex hull of the airfoil section
/// * `tol`: an optional tolerance value which will determine when to add new circles between the
/// existing circles. This value will default to 1e-3 if not specified.
///
/// returns: Result<Vec<InscribedCircle, Global>, Box<dyn Error, Global>>
pub fn extract_camber_line(
    section: &Curve2,
    hull: &ConvexPolygon,
    tol: Option<f64>,
) -> Result<Vec<InscribedCircle>> {
    // The convex hull will be used to determine the starting spanning ray for the camber line
    let (i0, i1) = farthest_pair_indices(hull);
    let p0 = &hull.points()[i0];
    let p1 = &hull.points()[i1];
    let dir = p1 - p0;
    let normal = rot90(Cw) * dir;
    let mid_ray = Ray::new(mid_point(p0, p1), normal);

    let spanning = section
        .try_create_spanning_ray(&mid_ray)
        .ok_or("Failed to create first spanning ray")?;

    // There was an orientation step in here, but I'm not sure what it did
    let mut stations0 = extract_half_camber_line(section, &spanning, tol)?;
    let stations1 = extract_half_camber_line(section, &spanning.reversed(), tol)?;

    reverse_inscribed_circles(&mut stations0);
    stations0.extend(stations1);

    Ok(stations0)
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
/// limit in one direction and the limit in the other direction is less than this value
///
/// returns: InscribedCircle
fn inscribed_from_spanning_ray(curve: &Curve2, ray: &SpanningRay, tol: f64) -> InscribedCircle {
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

/// Attempt to find the next spanning ray (for the next inscribed circle search) by advancing along
/// the known camber direction based on the previous inscribed circle.  This function will return
/// the next spanning ray if it is valid, or it will return an end condition if the search should
/// get too close to the edge of the airfoil.  If the search fails to find a valid spanning ray, it
/// will return a failure condition.
///
/// The search will reduce the distance that it attempts to jump forward from 25% of the last
/// inscribed circle's radius down to 5% of the last circle's radius.  This allows the search to
/// tolerate some failures and still make progress.
///
/// # Arguments
///
/// * `section`: the airfoil section curve
/// * `last_station`: the last inscribed circle found
///
/// returns: RayAdvance
fn advance_search_along_ray(section: &Curve2, last_station: &InscribedCircle) -> RayAdvance {
    // We will begin by finding the camber point/direction of the last station, which will be used
    // to jump forward and create a new spanning ray.  However, we'll first check the distance from
    // the camber point to the farthest point on the section in the camber direction.  As we get
    // closer to the edge of the airfoil, we will want to terminate the search.
    let camber_point = last_station.camber_point();

    // We unwrap this because the only way it would fail is if the section is empty, which
    // would have prevented us from getting here in the first place.
    let farthest = section
        .max_point_in_direction(&camber_point.normal)
        .unwrap();
    let distance = camber_point.scalar_projection(&farthest);

    // When the distance beyond the last inscribed circle is less than 25% of the circle's radius,
    // we will consider ourselves close enough to the edge of the airfoil to terminate the search.
    // Getting closer to the edge will increase the probability that the assumptions of no local
    // maxima along the ray are violated.
    if distance - last_station.radius() < last_station.radius() * 0.25 {
        return RayAdvance::End;
    }

    // Now we will create a new spanning ray which will be used to find the next inscribed circle.
    // We will start by jumping forward 25% of the last circle's radius, and we will adjust this
    // value down as we have failures.  So long as we move forward at least 5% of the last circle's
    // radius, we will consider the search to have advanced.
    let mut frac = 0.25;
    while frac > 0.05 {
        let next_center = camber_point.at_distance(frac * last_station.radius());
        let test_dir = rot90(Ccw) * camber_point.normal;
        let test_ray = Ray::new(next_center, test_dir.into_inner());

        if let Some(ray) = section.try_create_spanning_ray(&test_ray) {
            // First, we want to test if the new ray spans at least 50% of the last station's
            // distance between the upper and lower contact points.  This is a heuristic to ensure
            // we haven't taken a step where the section thickness is dropping off too quickly.
            let last_dist = dist(&last_station.upper, &last_station.lower);
            let new_dist = ray.dir().norm();

            if new_dist < 0.5 * last_dist {
                frac *= 0.75;
                continue;
            }

            return RayAdvance::Valid(ray);
        } else {
            frac *= 0.75;
        }
    }

    RayAdvance::Failed
}

/// This enum represents the result of trying to advance the camber search along the airfoil by
/// computing the next spanning ray.
enum RayAdvance {
    /// A valid spanning ray was computed, and the search can continue.
    Valid(SpanningRay),

    /// The search is close to the edge of the airfoil and should terminate.
    End,

    /// The search has failed to find a valid spanning ray
    Failed,
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
fn refine_stations(
    section: &Curve2,
    dest: &mut Vec<InscribedCircle>,
    stack: &mut Vec<InscribedCircle>,
    outer_tol: f64,
    inner_tol: f64,
) {
    while let Some(next) = stack.pop() {
        if let Some(last) = dest.last() {
            let test_ray = next.spanning_ray.symmetry(&last.spanning_ray);

            if let Some(ray) = section.try_create_spanning_ray(&test_ray) {
                let mid = inscribed_from_spanning_ray(section, &ray, inner_tol);
                let error = mid.interpolation_error(&next, last);

                // TODO: check the distance between the centers to make sure we're not stuck
                if error > outer_tol {
                    // We are out of tolerance, we need to put next back on the stack and then put
                    // the mid-ray on top of it and try again
                    stack.push(next);
                    stack.push(mid);
                } else {
                    // We are within tolerance, we can put the next station in the destination
                    dest.push(next);
                }
            }
        } else {
            dest.push(next);
        }
    }
}

/// Extracts the unambiguous portion of a mean camber line in the orthogonal direction to a
/// starting spanning ray. This function will terminate when it gets close to the farthest point in
/// the camber line direction.
///
/// # Arguments
///
/// * `curve`: the airfoil section curve
/// * `starting_ray`: the starting spanning ray for the camber line, determines the direction the
/// algorithm will advance.
/// * `tol`: an optional tolerance value which will determine when to add new circles between the
/// existing circles. This value will default to 1e-3 if not specified.
///
/// returns: Result<Vec<InscribedCircle, Global>, Box<dyn Error, Global>>
fn extract_half_camber_line(
    curve: &Curve2,
    starting_ray: &SpanningRay,
    tol: Option<f64>,
) -> Result<Vec<InscribedCircle>> {
    let outer_tol = tol.unwrap_or(1e-3);
    let inner_tol = outer_tol * 1e-2;

    let mut stations: Vec<InscribedCircle> = Vec::new();
    let mut refine_stack: Vec<InscribedCircle> = Vec::new();
    let mut ray = starting_ray.clone();

    loop {
        let circle = inscribed_from_spanning_ray(curve, &ray, inner_tol);
        refine_stack.push(circle);

        refine_stations(
            curve,
            &mut stations,
            &mut refine_stack,
            outer_tol,
            inner_tol,
        );

        let station = &stations.last().expect("Station was not transferred");

        match advance_search_along_ray(curve, station) {
            RayAdvance::Valid(r) => {
                ray = r;
            }
            RayAdvance::End => {
                break;
            }
            RayAdvance::Failed => {
                return Err(Box::from("Failed to advance search along ray"));
            }
        };
    }

    Ok(stations)
}
