use crate::AngleDir::Ccw;
use crate::airfoil::InscribedCircle;
use crate::common::dist;
use crate::geom2::hull::{convex_hull_2d, farthest_pair_on_hull};
use crate::geom2::{LineOps2, Segment2, rot90};
use crate::{AngleDir, Circle2, Curve2, Line2, Point2, Result, SurfacePoint2};
use parry2d_f64::query::Ray;

pub struct Inscribed {
    c: Circle2,
    p0: Point2,
    p1: Point2,
}

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

fn extract_half_circles(section: &Curve2, line: &Line2, tol: f64) -> Result<Vec<Inscribed>> {
    let mut results = Vec::new();

    let mut working_line = line.clone();

    loop {
        let circle = get_inscribed(section, &working_line, tol);
        results.push(circle);

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

fn get_inscribed(section: &Curve2, search_line: &Line2, tol: f64) -> Inscribed {
    let mut pos = InscribedSearchState::new(1.0, 0.0, search_line.at(1.0));
    let mut neg = InscribedSearchState::new(0.0, 0.0, search_line.at(0.0));
    let mut working;

    // While the distance between the positive and negative search bounds is greater than the
    // tolerance, continue to search for the inscribed circle center.
    while (pos.f - neg.f) * search_line.direction.norm() > tol {
        // We will update the working point to be right in the middle of the positive and negative
        // direction limits.
        let fraction = (pos.f + neg.f) * 0.5;
        working = search_line.at(fraction);

        // Now we find the closest position on the curve to the working point, and calculate the
        // distance and direction to that point. The direction will be used to determine which side
        // of the limits we will adjust.
        let closest = section.at_closest_to_point(&working);
        let to_closest = closest.point() - working; // The direction vector to the closest point
        let distance = dist(&working, &closest.point());
        let update = InscribedSearchState::new(fraction, distance, closest.point());

        // If the direction vector to the closest point is in the positive direction of the ray,
        // then we will adjust the positive limit.  Otherwise, we will adjust the negative limit.
        if to_closest.dot(&search_line.direction) > 0.0 {
            pos = update;
        } else {
            neg = update;
        }
    }

    // Finally, we will put the center of the inscribed circle at the midpoint of the positive and
    // negative limits, splitting the difference one last time, and we will set the radius to be
    // the average of the positive and negative distances. By this point the difference will be
    // below the tolerance value.
    let c = Circle2::from_point(search_line.at((pos.f + neg.f) * 0.5), (pos.d + neg.d) * 0.5);

    Inscribed {
        c,
        p0: neg.p,
        p1: pos.p,
    }
}

fn remaining(section: &Curve2, cp: &SurfacePoint2) -> f64 {
    let mut result = f64::NEG_INFINITY;
    for p in section.points() {
        if cp.scalar_projection(p).is_sign_positive() {
            result = result.max(dist(cp, p));
        }
    }

    result
}

fn advance_search(section: &Curve2, last: &Inscribed, tol: f64) -> Result<Option<Line2>> {
    // We will begin by finding the camber point/direction of the last station, which will be used
    // to jump forward and create a new spanning ray.  However, we'll first check the distance from
    // the camber point to the farthest point on the section in the camber direction.  As we get
    // closer to the edge of the airfoil, we will want to terminate the search.
    let d = rot90(AngleDir::Cw) * (last.p1 - last.p0);
    let cp = SurfacePoint2::new_normalize(last.c.center, d);

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
        // let test_line = Ray::new(next_center, test_dir.into_inner());
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

struct InscribedSearchState {
    f: f64,
    d: f64,
    p: Point2,
}

impl InscribedSearchState {
    pub fn new(f: f64, d: f64, p: Point2) -> Self {
        InscribedSearchState { f, d, p }
    }

    pub fn update(&mut self, fraction: f64, distance: f64, point: Point2) {
        self.f = fraction;
        self.d = distance;
        self.p = point;
    }
}

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

    #[test]
    fn simple_extract() -> Result<()> {
        use std::io::Write;
        let section = airfoil_curve();
        let circles = extract_inscribed_circles(&section, 1e-4)?;

        let mut file = std::fs::File::create("/tmp/circles.txt").unwrap();
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
