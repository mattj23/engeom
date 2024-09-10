//! Helper functions for implementing airfoil geometry algorithms.

use crate::airfoil::InscribedCircle;
use crate::{Curve2, Result};

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
