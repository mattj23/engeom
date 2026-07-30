//! This module contains tools for taking measurements on established [`AfGeometry`] entities.

use crate::airfoil2::{AfGeometry, AfPos, AfSide};
use crate::common::dist;
use crate::geom2::Segment2;
use crate::geom2::hull::convex_hull_2d;
use crate::metrology::Distance2;
use crate::{Curve2, Point2, Result, SurfacePoint2};

/// Measure the thickness of an airfoil section at a location specified by one of the gage point
/// location methods.
///
/// A corresponding pair of points is located on the lower and upper surfaces using
/// [`AfGeometry::af_point`] with the same `method` and `value` on each side, and the thickness is
/// the Euclidean distance between them.
///
/// # Arguments
///
/// * `geom`: the established airfoil geometry to measure.
/// * `method`: how `value` is interpreted, see [`AfPos`] for the three schemes.
/// * `value`: the position value, whose meaning and sign convention depend on `method`. A positive
///   value is measured from the leading edge and a negative one from the trailing edge.
///
/// returns: `Option<Distance2>` going from the lower surface point to the upper surface point, with
/// no fixed direction so that the value is the full unsigned distance between them. `None` when the
/// position does not land on both surfaces.
pub fn thickness(geom: &AfGeometry, method: AfPos, value: f64) -> Option<Distance2> {
    let lower = geom.af_point(AfSide::Lower, method, value)?;
    let upper = geom.af_point(AfSide::Upper, method, value)?;
    Some(Distance2::new(lower.point(), upper.point(), None))
}

/// Measure the maximum thickness of an airfoil section, taken from the largest inscribed circle.
///
/// The maximum thickness of an airfoil is conventionally defined as the diameter of its largest
/// inscribed circle, and the inscribed circle stack already contains it: the camber refinement in
/// [`SectionInput::inscribed_refined`](crate::airfoil2::SectionInput::inscribed_refined) inserts a
/// bisected circle whenever it is larger than both of its neighbors, so the search actively hunts
/// the radius peak rather than merely sampling through it.
///
/// The returned measurement runs between the two contact points of that circle, which are points
/// measured on the section itself rather than constructed on the circle. The value therefore agrees
/// with `geom.max_thickness_circle().radius() * 2.0` only to within the `resolve_tol` of the
/// analysis, since the contact points come from a binary search resolved to that tolerance while
/// the radius is the average of the two search bounds.
///
/// # This is not the maximum of [`thickness`]
///
/// Sweeping [`thickness`] with [`AfPos::OnCamber`] measures the chord cast orthogonal to the camber
/// line, which is not required to fit inside the section. On a cambered airfoil its maximum lands
/// at a different station and is larger: on the `airfoil-0` test section it exceeds the inscribed
/// diameter by roughly 0.8%. The two are equivalent only when the camber is straight or lightly
/// curved. This function reports the inscribed circle, which is the conventional definition of
/// maximum airfoil thickness, not the orthogonal-chord maximum.
///
/// # Arguments
///
/// * `geom`: the established airfoil geometry to measure.
///
/// returns: `Distance2` going from the lower surface contact point to the upper surface contact
/// point. This is infallible because an [`AfGeometry`] always holds at least one inscribed circle.
pub fn max_thickness(geom: &AfGeometry) -> Distance2 {
    let c = geom.max_thickness_circle();
    Distance2::new(c.p0, c.p1, None)
}

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
#[allow(dead_code)] // public helper not yet wired into the airfoil pipeline
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::airfoil2::{AfEdgeSearch, OrientFwdAft, OrientUpperLower};
    use crate::common::mid_point;
    use crate::metrology::Measurement;
    use crate::tests::airfoil_curve;
    use approx::assert_relative_eq;

    const TOL: f64 = 1e-3;

    fn analyzed() -> Result<AfGeometry> {
        AfGeometry::try_from_geometric_analysis(
            &airfoil_curve(),
            TOL,
            OrientFwdAft::TmaxFwd,
            OrientUpperLower::Curvature,
            AfEdgeSearch::Auto,
            AfEdgeSearch::Auto,
        )
    }

    #[test]
    fn thickness_endpoints_are_on_the_correct_surfaces() -> Result<()> {
        let geom = analyzed()?;
        let mid = geom.camber.length() * 0.5;
        let d = thickness(&geom, AfPos::OnCamber, mid).ok_or("Failed to measure thickness")?;

        assert!(d.value() > 0.0, "Thickness should be positive");

        // `a` is the lower surface point and `b` is the upper one
        assert_relative_eq!(geom.lower.dist_to_point(&d.a), 0.0, epsilon = 1e-9);
        assert_relative_eq!(geom.upper.dist_to_point(&d.b), 0.0, epsilon = 1e-9);

        // With no direction supplied the value is the full unsigned distance between the points
        assert_relative_eq!(d.value(), dist(&d.a, &d.b), epsilon = 1e-12);
        Ok(())
    }

    #[test]
    fn thickness_works_for_every_position_method() -> Result<()> {
        let geom = analyzed()?;
        let tmax = geom.max_thickness_circle().radius() * 2.0;

        // Each method is given a value that is meaningful for it. `EdgeOffset` steps along the
        // camber tangent at the edge rather than along the camber itself, so on a cambered section
        // it only reaches the surface near the edge, which is what it is intended for.
        let cases = [
            (AfPos::OnCamber, geom.camber.length() * 0.5),
            (AfPos::Radius, geom.camber.length() * 0.25),
            (AfPos::EdgeOffset, geom.camber.length() * 0.02),
        ];

        for (method, value) in cases {
            let d = thickness(&geom, method, value)
                .ok_or_else(|| format!("Failed to measure thickness with {method:?}"))?;
            assert!(
                d.value() > 0.0 && d.value() <= tmax,
                "Thickness of {} with {:?} is not in (0, {}]",
                d.value(),
                method,
                tmax
            );
        }
        Ok(())
    }

    #[test]
    fn thickness_beyond_the_camber_is_none() -> Result<()> {
        let geom = analyzed()?;
        let past_end = geom.camber.length() * 1.5;
        assert!(thickness(&geom, AfPos::OnCamber, past_end).is_none());
        Ok(())
    }

    #[test]
    fn max_thickness_matches_the_largest_inscribed_circle() -> Result<()> {
        let geom = analyzed()?;
        let d = max_thickness(&geom);

        // The contact points are located by a binary search resolved to `resolve_tol`, and the
        // circle radius is the average of the two search bounds, so the chord between the contact
        // points can differ from the diameter by about that much in either direction.
        let diameter = geom.max_thickness_circle().radius() * 2.0;
        assert_relative_eq!(d.value(), diameter, epsilon = TOL * 0.1);

        assert_relative_eq!(geom.lower.dist_to_point(&d.a), 0.0, epsilon = 1e-9);
        assert_relative_eq!(geom.upper.dist_to_point(&d.b), 0.0, epsilon = 1e-9);
        Ok(())
    }

    /// Pins the relationship between the two competing definitions of maximum thickness, because
    /// it is not the one the reference literature asserts.
    ///
    /// `max_thickness` reports the largest inscribed circle. Sweeping `thickness` along the camber
    /// instead reports the longest chord cast orthogonal to the camber. On a cambered section these
    /// are *not* the same: the orthogonal chord is not constrained to fit inside the section, so
    /// its maximum lands at a different station and comes out larger. On the test airfoil the sweep
    /// exceeds the inscribed diameter by roughly 0.8%.
    ///
    /// The check that the sweep peak is genuinely not an inscribed diameter is what distinguishes
    /// this from a sampling artifact.
    #[test]
    fn max_thickness_is_the_inscribed_circle_not_the_orthogonal_maximum() -> Result<()> {
        let geom = analyzed()?;
        let inscribed = max_thickness(&geom).value();

        let steps = 2000;
        let length = geom.camber.length();
        let mut best: Option<Distance2> = None;
        for i in 1..steps {
            let l = length * (i as f64) / (steps as f64);
            if let Some(d) = thickness(&geom, AfPos::OnCamber, l)
                && best.as_ref().is_none_or(|b| d.value() > b.value())
            {
                best = Some(d);
            }
        }
        let best = best.ok_or("The camber sweep found no valid positions")?;

        assert!(
            best.value() > inscribed,
            "Expected the orthogonal sweep ({}) to exceed the inscribed diameter ({})",
            best.value(),
            inscribed
        );

        // The sweep peak is a real chord across the section, but a circle of that diameter would
        // not fit inside it, which is why it is not the maximum inscribed circle
        let section = airfoil_curve();
        let clearance = section.dist_to_point(&mid_point(&best.a, &best.b));
        assert!(
            clearance < best.value() / 2.0,
            "The sweep peak chord of {} has clearance {} at its midpoint, so it would fit as an \
             inscribed circle and the two definitions should not have disagreed",
            best.value(),
            clearance
        );

        // The tmax circle, by contrast, really is inscribed
        let c = geom.max_thickness_circle();
        assert_relative_eq!(
            section.dist_to_point(&mid_point(&c.p0, &c.p1)),
            c.radius(),
            epsilon = TOL * 0.1
        );
        Ok(())
    }
}
