use crate::AngleDir::Ccw;
use crate::airfoil2::SectionInput;
use crate::airfoil2::inscribed::{Inscribed, InscribedVec};
use crate::common::dist;
use crate::geom2::hull::{convex_hull_2d, farthest_pair_on_hull};
use crate::geom2::{LineOps2, rot90};
use crate::{Curve2, Line2, Result, SurfacePoint2, Vector2};

const EDGE_STOP_FRACTION: f64 = 0.375;

/// Extract the unambiguous inscribed circles of an airfoil section.
///
/// - The circles will be returned in a consistent order, but the order may be _either_ from leading
///   to trailing edge _or_ the reverse.
/// - The `p0` and `p1` points will all be oriented in a consistent way, but the order may be
///   _either_ from upper to lower surface _or_ the reverse.
///
/// The camber line extraction algorithm will terminate when the distance left to the farthest
/// point on the edge is more than 3/8's of the radius of the last inscribed circle beyond its
/// perimeter.  From there, edge-aware algorithms should be used to extend the camber line and any
/// remaining inscribed circles.
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
pub fn extract_inscribed_circles(input: &SectionInput) -> Result<Vec<Inscribed>> {
    let hull = convex_hull_2d(input.section.points());
    let (i0, i1) = farthest_pair_on_hull(input.section.points(), &hull);
    let naive_chord = Line2::from_points(&input.section.points()[i0], &input.section.points()[i1]);

    // The starting line will be halfway down, but ultimately we're going to want to be able to
    // sample at multiple places along the chord and look for valid crossing lines, with some
    // mechanism to guess at their quality
    //
    // Alternately, perhaps something with the zhang suen line can be a good starting point
    let ref_line = Line2::new(naive_chord.at(0.5), naive_chord.orthogonal());

    // First we'll gather the circles in the front half
    let start_line = input.crossing_line(&ref_line).ok_or_else(|| {
        "Failed to find a crossing line for the initial reference line".to_string()
    })?;
    let mut working = extract_half_circles(input, &start_line)?;

    // We'll reverse the vector so that we can add the elements from the back half of the section,
    // removing the last element since it will be the same as the first element of the back half.
    working.reverse_order();
    working.pop();

    // Now we gather the circles in the back half
    let start_line = input
        .crossing_line(&ref_line.new_reversed())
        .ok_or_else(|| {
            "Failed to find a crossing line for the reversed reference line".to_string()
        })?;
    working.extend(extract_half_circles(input, &start_line)?);

    Ok(working.take_vec())
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
fn extract_half_circles(input: &SectionInput, line: &Line2) -> Result<InscribedVec> {
    let mut results = InscribedVec::empty();
    let mut working_line = line.clone();

    loop {
        let circle = input.inscribed_from_crossing(&working_line);
        results.refine_and_push(circle, input);
        let last = results
            .last()
            .ok_or("Failed to get last inscribed circle".to_string())?;

        match advance_search(input, last)? {
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
fn advance_search(input: &SectionInput, last: &Inscribed) -> Result<Option<Line2>> {
    // We will begin by finding the camber point/direction of the last station, which will be used
    // to jump forward and create a new spanning ray.  However, we'll first check the distance from
    // the camber point to the farthest point on the section in the camber direction.  As we get
    // closer to the edge of the airfoil, we will want to terminate the search.
    let cp = last.camber_point();

    // We unwrap this because the only way it would fail is if the section is empty, which
    // would have prevented us from getting here in the first place.
    // let (_, farthest) = section.max_point_in_direction(&cp.normal).unwrap();
    // let distance = cp.scalar_projection(&farthest);
    let distance = remaining(input.section, &cp);

    // When the distance beyond the last inscribed circle is less than 3/8 of the circle's radius,
    // we will consider ourselves close enough to the edge of the airfoil to terminate the search.
    // Getting closer to the edge will increase the probability that the assumptions of no local
    // maxima along the ray are violated.
    if distance - last.c.r() < last.c.r() * EDGE_STOP_FRACTION {
        return Ok(None);
    }

    // Now we will create a new spanning ray which will be used to find the next inscribed circle.
    // We will start by jumping forward 1/8 of the last circle's radius, and we will adjust this
    // value down as we have failures.  So long as we move forward at least 5% of the last circle's
    // radius, we will consider the search to have advanced.
    let mut frac = 0.125;
    while frac > 0.05 {
        let advance_dist = frac * last.c.r();
        if advance_dist < input.general_tol {
            return Ok(None);
        }
        let next_center = cp.at_distance(advance_dist);
        let test_dir = rot90(Ccw) * cp.normal;
        let test_line = Line2::new(next_center, test_dir.into_inner());

        if let Some(line) = input.crossing_line(&test_line) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::airfoil2::{OrientFwdAft, OrientUpperLower};
    use crate::tests::airfoil_curve;
    use std::env::temp_dir;

    #[test]
    fn simple_extract() -> Result<()> {
        use std::io::Write;
        let section = airfoil_curve();
        let input = SectionInput::new(&section, 1e-3);
        let circles = extract_inscribed_circles(&input)?;
        let circles = OrientFwdAft::TmaxFwd.apply(circles)?;
        let circles = OrientUpperLower::Curvature.apply(circles)?;

        let output = temp_dir().join("circles.txt");
        let mut file = std::fs::File::create(&output)?;
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
            )?
        }

        Ok(())
    }
}
