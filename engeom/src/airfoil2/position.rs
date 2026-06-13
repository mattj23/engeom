//! This module has the implementations for locating positions on the
//! airfoil sides based on the different methods

use crate::common::Intersection;
use crate::{Circle2, Curve2, CurveStation2, SurfacePoint2};
use std::cmp::Ordering;

/// Locate a station on `target` by intersecting it with a circle centered at one of the camber
/// endpoints. A positive `value` centers the circle at the first camber vertex (leading edge); a
/// negative `value` centers it at the last camber vertex (trailing edge). The radius of the test
/// circle is `value.abs()`.
///
/// When more than one intersection exists, the one whose projection onto `target` is closest in
/// arc length to the projection of the edge point is returned, which keeps the result on the side
/// of the edge nearest the edge feature.
pub fn pos_radius<'a>(
    value: f64,
    camber: &Curve2,
    target: &'a Curve2,
) -> Option<CurveStation2<'a>> {
    let center = if value >= 0.0 {
        camber.at_front().point()
    } else {
        camber.at_back().point()
    };
    let circle = Circle2::from_point(center, value.abs());
    let hits = target.intersection(&circle);
    if hits.is_empty() {
        return None;
    }

    let l_edge = target.at_closest_to_point(&center).length_along();
    let best = hits.into_iter().min_by(|a, b| {
        let la = (target.at_closest_to_point(a).length_along() - l_edge).abs();
        let lb = (target.at_closest_to_point(b).length_along() - l_edge).abs();
        la.partial_cmp(&lb).unwrap_or(Ordering::Equal)
    })?;
    Some(target.at_closest_to_point(&best))
}

/// Locate a station on `target` by walking an arc distance along the camber curve and casting an
/// orthogonal line at that station. A positive `value` is measured from the leading edge (first
/// camber vertex) along the camber; a negative `value` is measured from the trailing edge.
///
/// Returns `None` if `value.abs()` exceeds the camber length or if the orthogonal line does not
/// intersect `target`. The intersection nearest the camber station is returned.
pub fn pos_camber<'a>(
    value: f64,
    camber: &Curve2,
    target: &'a Curve2,
) -> Option<CurveStation2<'a>> {
    let l = if value >= 0.0 {
        value
    } else {
        camber.length() + value
    };
    let station = camber.at_length(l)?;
    closest_orthogonal_hit(&station.surface_point(), target)
}

/// Locate a station on `target` by stepping along a ray that starts at the leading or trailing
/// edge point and runs in the direction of the camber tangency at that edge, then casting an
/// orthogonal line at the stepped position. A positive `value` originates at the leading edge; a
/// negative `value` originates at the trailing edge.
///
/// Returns `None` if the orthogonal line at the offset position does not intersect `target`.
pub fn pos_offset<'a>(
    value: f64,
    camber: &Curve2,
    target: &'a Curve2,
) -> Option<CurveStation2<'a>> {
    let edge = if value >= 0.0 {
        camber.at_front()
    } else {
        camber.at_back()
    };
    let origin = edge.point() + edge.direction().into_inner() * value;
    let sp = SurfacePoint2::new(origin, edge.normal());
    closest_orthogonal_hit(&sp, target)
}

fn closest_orthogonal_hit<'a>(sp: &SurfacePoint2, target: &'a Curve2) -> Option<CurveStation2<'a>> {
    let hits = target.intersections_with_line(sp);
    let (t, _) = hits
        .into_iter()
        .min_by(|a, b| a.0.abs().partial_cmp(&b.0.abs()).unwrap_or(Ordering::Equal))?;
    Some(target.at_closest_to_point(&sp.at_distance(t)))
}
