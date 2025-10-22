use crate::AngleDir::Ccw;
use crate::airfoil::helpers::{
    OrientedCircles, extract_curve_beyond_station, extract_edge_sub_curve,
    inscribed_from_spanning_ray, refine_stations,
};
use crate::airfoil::{AirfoilEdge, EdgeGeometry, EdgeLocate, InscribedCircle, OpenIntersectGap};
use crate::common::points::{dist, mid_point};
use crate::common::{BestFit, linear_space};
use crate::geom2::{Ray2, Segment2, rot90};
use crate::{Circle2, Curve2, Result};
use parry2d_f64::query::Ray;

/// This struct implements the `EdgeLocation` trait and attempts to locate the edge of the airfoil
/// using a sequence of simple detection methods.  It will check for an open edge, a sharp corner,
/// a square edge, a full radius, and finally a squared edge with rounded corners. This detection
/// method should only be used on nominal airfoil sections and well-defined edges.  It is most
/// useful on trailing edges, which tend to have simpler geometries than leading edges.
pub struct EdgeAutoDetect {
    tol: f64,
}

impl EdgeAutoDetect {
    pub fn new(tol: f64) -> Self {
        Self { tol }
    }

    /// Create a new boxed instance of the `Detect` struct.
    pub fn make(tol: f64) -> Box<dyn EdgeLocate> {
        Box::new(Self::new(tol))
    }
}

impl EdgeLocate for EdgeAutoDetect {
    fn find_edge(
        &self,
        section: &Curve2,
        stations: Vec<InscribedCircle>,
        front: bool,
        af_tol: f64,
    ) -> Result<(Option<AirfoilEdge>, Vec<InscribedCircle>)> {
        let mut working_stations = OrientedCircles::new(stations, front);

        let station = working_stations
            .last()
            .ok_or("Empty inscribed circles container.")?;

        // Open edge detection
        // ----------------------------------------------------------------------------------------
        // Check if we get an intersection at the end. If not, this is probably an open edge, and
        // we should return an open intersection edge
        let Ok(te_intersect) = working_stations.intersect_from_end(&section) else {
            let open = OpenIntersectGap::new(10);
            return open.find_edge(section, working_stations.take_circles(), front, af_tol);
        };

        // If this isn't an open edge, we can try to extract the edge curve beyond the last station
        let end_sp = working_stations.end_sp()?.normal;
        let edge_curve = extract_curve_beyond_station(section, station, &end_sp)
            .ok_or("Failed to extract edge curve on edge detection algorithm.")?;

        // Corner-based edge detection (sharp or square)
        // ----------------------------------------------------------------------------------------
        // First we'll check for sharp corners in order to try to identify a sharp or square edge.
        let corners = corner_vertices(&edge_curve, 110.0f64.to_radians());

        // A square edge is defined as having the first two corners have a combined angle of less
        // than 180 degrees
        if corners.len() >= 2 && (corners[0].1 + corners[1].1) < 210.0f64.to_radians() {
            let p0 = edge_curve.points()[corners[0].0];
            let p1 = edge_curve.points()[corners[1].0];
            let edge = AirfoilEdge::square(te_intersect, (p0, p1));
            return Ok((Some(edge), working_stations.take_circles()));
        }

        // A sharp corner edge is defined as having a single corner with an angle less than 90
        // degrees
        if corners.len() == 1 && corners[0].1 < 90.0f64.to_radians() {
            // TODO: Refine stations and extract edge curve up to corner
            let p0 = edge_curve.points()[corners[0].0];
            let edge = AirfoilEdge::sharp(p0);
            return Ok((Some(edge), working_stations.take_circles()));
        }

        // Full radius edge detection
        // ----------------------------------------------------------------------------------------

        // Rounded square edge detection
        // ----------------------------------------------------------------------------------------

        todo!()
    }
}

/// Identify the indices of vertices that represent corners in the curve, where the angle between
/// the preceding and following segments is _smaller_ than `angle_max`.  The function returns a
/// vector of tuples, each containing the index of the corner vertex and the angle at that vertex
/// in radians. The returned vector is sorted in ascending order based on the angle values.
///
/// # Arguments
///
/// * `curve`: The input curve to analyze for corner vertices.
/// * `angle_max`: The maximum angle (in radians) to consider a vertex as a corner.
///
/// returns: Vec<(usize, f64), Global>
fn corner_vertices(curve: &Curve2, angle_max: f64) -> Vec<(usize, f64)> {
    let mut corners = Vec::new();
    if curve.points().len() < 3 {
        return corners;
    }

    for i in 1..(curve.points().len() - 1) {
        let p_prev = curve.points()[i - 1];
        let p_curr = curve.points()[i];
        let p_next = curve.points()[i + 1];

        let v1 = (p_prev - p_curr).normalize();
        let v2 = (p_next - p_curr).normalize();

        let angle = v1.angle(&v2).abs();

        if angle < angle_max {
            corners.push((i, angle));
        }
    }

    corners.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    corners
}
