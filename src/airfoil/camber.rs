//! This module contains the implementation of airfoil camber line detection algorithms.

use super::helpers::{
    inscribed_from_spanning_ray, refine_stations, reverse_inscribed_circles, OrientedCircles,
};
use crate::airfoil::inscribed_circle::InscribedCircle;
use crate::common::points::{dist, mid_point};
use crate::geom2::hull::farthest_pair_indices;
use crate::geom2::polyline2::SpanningRay;
use crate::geom2::{rot90, Line2};
use crate::AngleDir::{Ccw, Cw};
use crate::Curve2;
use crate::Result;
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

    let mut stations0 = extract_half_camber_line(section, &spanning, tol)?;
    let stations1 = extract_half_camber_line(section, &spanning.reversed(), tol)?;

    reverse_inscribed_circles(&mut stations0);
    stations0.extend(stations1);

    Ok(stations0)
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
    let (_, farthest) = section
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

    let mut stations = OrientedCircles::create(false);
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

    Ok(stations.take_circles())
}
