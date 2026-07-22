//! This module contains tools for taking measurements on established [`AfGeometry`] entities.

use crate::common::dist;
use crate::geom2::Segment2;
use crate::geom2::hull::convex_hull_2d;
use crate::{Curve2, Point2, Result, SurfacePoint2};

/// This function calculates the chord line of an airfoil section using the "caliper method". The
/// caliper method is a simple method that works on highly curved airfoils, but is an artifact of
/// legacy airfoil analysis methods and is not recommended for use with modern airfoil sections.
/// Don't use this method unless you know that you need it. Depending on the use case, it is
/// unlikely that the aerodynamic properties of the airfoil will be well represented by the chord
/// length calculated by this method.
///
/// The "caliper method" gained prominence when trying to measure highly cambered turbine airfoils,
/// and replicates a physical method that would consist of the following:
///
/// 1. The pressure side of the airfoil is rested against a straight-edge or a flat surface, such
///    that it makes contact with the surface at the leading and trailing edges, while its center
///    bows up away from the surface.
///
/// 2. A pair of calipers is used to measure the span of the leading to trailing edge by putting
///    tips of the jaws of the calipers in contact with the straight-edge, and then closing them
///    until the flats of the jaw touch the airfoil somewhere near the leading and trailing edges.
///    The jaws and the straight-edge form a rectangle with right angles that closes on the airfoil.
///
/// Computationally, this method involves calculating the convex hull of the airfoil points and then
/// finding the longest straight line that can be drawn between two points on the hull. This line
/// represents the flat surface that the airfoil would be resting against in the physical method,
/// and is also a line of tangency sometimes used to measure airfoil twist.
///
/// Once the line of tangency from leading to trailing edge is found, all points in the airfoil
/// section are projected onto the line and the two extremes are found.  These points would
/// represent the location of the tips of the calipers in the physical method, and the distance
/// between them is the chord length found by this technique.
///
/// # Arguments
///
/// * `section`: the airfoil section to analyze
/// * `camber`: the mean camber line associated with the airfoil section
///
/// returns: Result<CaliperChord, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn caliper_chord_line(section: &Curve2, camber: &Curve2) -> Result<(Segment2, Segment2)> {
    // The tangent chord line is found through the caliper method.  We look at the convex hull and
    // find the longest straight line that can be drawn between two points on the hull.  This line
    // is the line of tangency for the section.  Next we find the furthest forward and furthest
    // backwards projections of the airfoil outer boundary onto this line.  These points are the
    // leading and trailing edges of the chord, and the distance between them is equivalent to the
    // result of the caliper chord method, problematic as it is.

    let hull_indices = convex_hull_2d(section.points());

    // First find the longest leg of the hull
    let mut max_dist = 0.0;
    let mut max_p1 = Point2::origin();
    let mut max_p2 = Point2::origin();

    for i in 0..hull_indices.len() {
        let i1 = hull_indices[i];
        let i2 = hull_indices[(i + 1) % hull_indices.len()];
        let p1 = section.points()[i1];
        let p2 = section.points()[i2];
        let d = dist(&p1, &p2);
        if d > max_dist {
            max_dist = d;
            max_p1 = p1;
            max_p2 = p2;
        }
    }

    // Now orient it from the leading edge to the trailing edge
    let camber_le = camber.at_front().point();
    let chord = if dist(&max_p1, &camber_le) < dist(&max_p2, &camber_le) {
        SurfacePoint2::new_normalize(max_p1, max_p2 - max_p1)
    } else {
        SurfacePoint2::new_normalize(max_p2, max_p1 - max_p2)
    };

    // Now find the highest and lowest projection parameters on the chord line
    let te = section
        .max_point_in_direction(&chord.normal)
        .ok_or("Failed to find trailing edge")?;
    let le = section
        .max_point_in_direction(&-chord.normal)
        .ok_or("Failed to find leading edge")?;

    let _chord_line = Segment2::new(&le.1, &te.1)?;
    let _tangent_line = Segment2::new(&chord.projection(&le.1), &chord.projection(&te.1))?;

    // Ok((chord_line, tangent_line))
    todo!()
}
