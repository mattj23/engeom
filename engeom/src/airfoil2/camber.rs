use crate::AngleDir::Ccw;
use crate::common::dist;
use crate::geom2::hull::{convex_hull_2d, farthest_pair_on_hull};
use crate::geom2::{LineOps2, rot90};
use crate::{AngleDir, Circle2, Curve2, Line2, Point2, Result, SurfacePoint2};
use serde::{Deserialize, Serialize};

const CIRCLE_TOL_FACTOR: f64 = 0.1;

/// Represents an inscribed circle inside an airfoil section. Contains the circle itself (center
/// and radius) and the two contact point with the perimeter of the section. The circle center is
/// definitionally on the mean camber line. The two points are not guaranteed to be in any order.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Inscribed {
    c: Circle2,
    p0: Point2,
    p1: Point2,
}

impl Inscribed {
    fn new(c: Circle2, p0: Point2, p1: Point2) -> Self {
        Self { c, p0, p1 }
    }

    fn camber_point(&self) -> SurfacePoint2 {
        let d = rot90(AngleDir::Cw) * (self.p1 - self.p0);
        SurfacePoint2::new_normalize(self.c.center, d)
    }
}

/// Extract the unambiguous inscribed circles of an airfoil section.
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a curve. May be open at the leading or trailing
///   edge but not both.  Should not be open in the middle unless the gap is small compared to the
///   thickness of the airfoil at the open portion.
/// * `tol`: A tolerance for extracting the circles. More circles will be added until interpolation
///   between neighboring circles has an error less than this amount. The circle centers will be
///   located with a tolerance 1/10th of this value.
///
/// returns: Result<Vec<Inscribed, Global>, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn extract_inscribed_circles(section: &Curve2, tol: f64) -> Result<Vec<Inscribed>> {
    let hull = convex_hull_2d(section.points());
    let (i0, i1) = farthest_pair_on_hull(section.points(), &hull);
    let naive_chord = Line2::from_points(&section.points()[i0], &section.points()[i1]);

    // The starting line will be halfway down, but ultimately we're going to want to be able to
    // sample at multiple places along the chord and look for valid crossing lines, with some
    // mechanism to guess at their quality
    //
    // Alternately, perhaps something with the zhang suen line can be a good starting point
    let ref_line = Line2::new(naive_chord.at(0.5), naive_chord.orthogonal());
    let start_line = crossing_line(section, &ref_line).ok_or_else(|| {
        "Failed to find a crossing line for the initial reference line".to_string()
    })?;

    let mut working = extract_half_circles(section, &start_line, tol)?;
    working.reverse();

    let start_line = crossing_line(section, &ref_line.new_reversed()).ok_or_else(|| {
        "Failed to find a crossing line for the reversed reference line".to_string()
    })?;
    working.extend(extract_half_circles(section, &start_line, tol)?);

    Ok(working)
}

/// Extracts the inscribed circles starting at the initial line and proceeding in the direction of
/// the line orthogonal direction until reaching the edge neighborhood or throwing an error.
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a `Curve2`
/// * `line`: A valid crossing line for the airfoil (t=0 and t=1 must be on the section perimeter)
/// * `tol`: A tolerance for extracting the circles. More circles will be added until interpolation
///   between neighboring circles has an error less than this amount. The circle centers will be
///   located with a tolerance 1/10th of this value.
///
/// returns: Result<Vec<Inscribed, Global>, Box<dyn Error, Global>>
fn extract_half_circles(section: &Curve2, line: &Line2, tol: f64) -> Result<Vec<Inscribed>> {
    let mut results = Vec::new();

    let mut working_line = line.clone();

    loop {
        let circle = get_inscribed(section, &working_line, tol * CIRCLE_TOL_FACTOR);
        add_refined_up_to(&mut results, section, circle, tol);

        let last = results
            .last()
            .ok_or("Failed to get last inscribed circle".to_string())?;

        match advance_search(section, last, tol)? {
            Some(next_line) => {
                working_line = next_line;
            }
            None => {
                break;
            }
        };
    }

    Ok(results)
}

/// Given a `Vec` of inscribed circles, this function will append more inscribed circles to it such
/// that no pair of adjacent circles, starting with the current last element, violates the
/// interpolation tolerance. The last circle appended will be the one supplied in the argument.
fn add_refined_up_to(items: &mut Vec<Inscribed>, section: &Curve2, next: Inscribed, tol: f64) {
    if items.is_empty() {
        items.push(next);
        return;
    }

    let mut stack = vec![next];
    while let Some(top) = stack.pop() {
        let last = items.last().unwrap().clone();
        if let Some(to_push) = refine_check(section, &last, &top, tol) {
            stack.push(top);
            stack.push(to_push);
        } else {
            items.push(top);
        }
    }
}

/// Checks if a refinement step is needed between inscribed circles `a` and `b` to bring the linear
/// interpolation of radius and circle center to within the supplied tolerance. If refinement is
/// necessary, the inscribed circle between the two stations is returned. If it isn't necessary, a
/// `None` is returned instead.
fn refine_check(section: &Curve2, a: &Inscribed, b: &Inscribed, tol: f64) -> Option<Inscribed> {
    let ca = a.camber_point();
    let cb = b.camber_point();

    // Get the interpolation fraction (may need to be varied)
    let f = 0.5;

    // Get the interpolated point and radius
    let cp = ca.slerp(&cb, f);
    let cr = (b.c.r() - a.c.r()) * f + a.c.r();

    // Now find the measured circle
    let line = Line2::from(&cp.rot_normal_90(Ccw));
    let line = crossing_line(section, &line).expect("Implement failed crossings checks");
    let circle = get_inscribed(section, &line, tol * CIRCLE_TOL_FACTOR);

    if (circle.c.r() - cr).abs() < tol && dist(&cp, &circle.c.center) < tol {
        None
    } else {
        Some(circle)
    }
}

/// Attempts to advance the inscribed circle search by jumping forward from the last inscribed
/// circle by a fraction of its radius and seeking a valid crossing line. If it fails it will
/// attempt again at reduced distances. If it is close to the edge it will return `None` to indicate
/// that we are getting close to the ambiguous portion of the MCL and the search should terminate.
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a `Curve2`
/// * `last`: The last valid inscribed circle
/// * `tol`: The circle finding tolerance, which is used as a lower limit of when to terminate the
///   advancement attempt when all else fails.
///
/// returns: Result<Option<Line2>, Box<dyn Error, Global>>
///
/// # Warnings
///
/// If the airfoil is curved enough that part of the other side lies in front of the current
/// station's search direction, this function will produce unknown results and will need to be
/// modified.
fn advance_search(section: &Curve2, last: &Inscribed, tol: f64) -> Result<Option<Line2>> {
    // We will begin by finding the camber point/direction of the last station, which will be used
    // to jump forward and create a new spanning ray.  However, we'll first check the distance from
    // the camber point to the farthest point on the section in the camber direction.  As we get
    // closer to the edge of the airfoil, we will want to terminate the search.
    let cp = last.camber_point();

    // We unwrap this because the only way it would fail is if the section is empty, which
    // would have prevented us from getting here in the first place.
    // let (_, farthest) = section.max_point_in_direction(&cp.normal).unwrap();
    // let distance = cp.scalar_projection(&farthest);
    let distance = remaining(section, &cp);

    // When the distance beyond the last inscribed circle is less than 3/8 of the circle's radius,
    // we will consider ourselves close enough to the edge of the airfoil to terminate the search.
    // Getting closer to the edge will increase the probability that the assumptions of no local
    // maxima along the ray are violated.
    if distance - last.c.r() < last.c.r() * 0.375 {
        return Ok(None);
    }

    // Now we will create a new spanning ray which will be used to find the next inscribed circle.
    // We will start by jumping forward 25% of the last circle's radius, and we will adjust this
    // value down as we have failures.  So long as we move forward at least 5% of the last circle's
    // radius, we will consider the search to have advanced.
    let mut frac = 0.125;
    while frac > 0.05 {
        let advance_dist = frac * last.c.r();
        if advance_dist < tol {
            return Ok(None);
        }
        let next_center = cp.at_distance(advance_dist);
        let test_dir = rot90(Ccw) * cp.normal;
        let test_line = Line2::new(next_center, test_dir.into_inner());

        if let Some(line) = crossing_line(section, &test_line) {
            // First, we want to test if the new ray spans at least 50% of the last station's
            // distance between the upper and lower contact points.  This is a heuristic to ensure
            // we haven't taken a step where the section thickness is dropping off too quickly.
            let last_dist = dist(&last.p0, &last.p1);
            let new_dist = line.direction.norm();

            if new_dist < 0.5 * last_dist {
                frac *= 0.75;
                continue;
            }

            return Ok(Some(line));
        } else {
            frac *= 0.75;
        }
    }

    Err("Failed to advance search".into())
}

/// This is a local helper function which tries to find the maximum remaining distance in front of
/// a surface point. It's used to see how close to the end of the section we are.
fn remaining(section: &Curve2, cp: &SurfacePoint2) -> f64 {
    let mut result = f64::NEG_INFINITY;
    for p in section.points() {
        if cp.scalar_projection(p).is_sign_positive() {
            result = result.max(dist(cp, p));
        }
    }

    result
}

/// Given a valid crossing line, find the inscribed circle center along the line with a tolerance
/// of `tol`.  Internally, this uses a binary search.
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a `Curve2`
/// * `line`: A valid crossing line for the airfoil (t=0 and t=1 must be on the section perimeter)
/// * `tol`: The search will end when the located position is within this amount of the real value
///
/// returns: Inscribed
fn get_inscribed(section: &Curve2, line: &Line2, tol: f64) -> Inscribed {
    struct InscribedEnd {
        f: f64,
        d: f64,
        p: Point2,
    }

    impl InscribedEnd {
        pub fn new(f: f64, d: f64, p: Point2) -> Self {
            InscribedEnd { f, d, p }
        }
    }

    let mut pos = InscribedEnd::new(1.0, 0.0, line.at(1.0));
    let mut neg = InscribedEnd::new(0.0, 0.0, line.at(0.0));
    let mut working;

    // While the distance between the positive and negative search bounds is greater than the
    // tolerance, continue to search for the inscribed circle center.
    while (pos.f - neg.f) * line.direction.norm() > tol {
        // TODO: Once we are close, there may be a geometric way to skip a bunch of iterations
        // We will update the working point to be right in the middle of the positive and negative
        // direction limits.
        let fraction = (pos.f + neg.f) * 0.5;
        working = line.at(fraction);

        // Now we find the closest position on the curve to the working point, and calculate the
        // distance and direction to that point. The direction will be used to determine which side
        // of the limits we will adjust.
        let closest = section.at_closest_to_point(&working);
        let to_closest = closest.point() - working; // The direction vector to the closest point
        let distance = dist(&working, &closest.point());
        let update = InscribedEnd::new(fraction, distance, closest.point());

        // If the direction vector to the closest point is in the positive direction of the ray,
        // then we will adjust the positive limit.  Otherwise, we will adjust the negative limit.
        if to_closest.dot(&line.direction) > 0.0 {
            pos = update;
        } else {
            neg = update;
        }
    }

    // Finally, we will put the center of the inscribed circle at the midpoint of the positive and
    // negative limits, splitting the difference one last time, and we will set the radius to be
    // the average of the positive and negative distances. By this point the difference will be
    // below the tolerance value.
    let c = Circle2::from_point(line.at((pos.f + neg.f) * 0.5), (pos.d + neg.d) * 0.5);

    Inscribed::new(c, neg.p, pos.p)
}

/// From a donor line, try to find a valid line which crosses the airfoil section. That means that
/// it has exactly two intersections with the section, and is returned such that t=0 and t=1 are
/// both exactly on the section perimeter.
///
/// If the donor line does not result in exactly two intersections this function returns a `None`
///
/// # Arguments
///
/// * `section`: An airfoil section represented by a `Curve2`
/// * `donor`: Any parametric line coincident with the desired crossing line, does not need to
///   start or end on the section perimeter itself.
///
/// returns: Option<Line2>
fn crossing_line(section: &Curve2, donor: &Line2) -> Option<Line2> {
    let ts = section
        .intersections_with_line(&donor)
        .iter()
        .map(|(t, _)| *t)
        .collect::<Vec<_>>();
    if ts.len() != 2 {
        None
    } else {
        Some(Line2::from_points(&donor.at(ts[0]), &donor.at(ts[1])))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::airfoil_curve;
    use std::env::temp_dir;

    #[test]
    fn simple_extract() -> Result<()> {
        use std::io::Write;
        let section = airfoil_curve();
        let circles = extract_inscribed_circles(&section, 1e-3)?;
        let output = temp_dir().join("circles.txt");

        let mut file = std::fs::File::create(&output).unwrap();
        for circle in &circles {
            writeln!(
                file,
                "{} {} {} {} {} {} {}",
                circle.c.center.x,
                circle.c.center.y,
                circle.c.r(),
                circle.p0.x,
                circle.p0.y,
                circle.p1.x,
                circle.p1.y
            )
            .unwrap();
        }

        Ok(())
    }
}
