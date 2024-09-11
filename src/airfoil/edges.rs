//! This module contains tools to work with the leading and trailing edges of the airfoil section.

use crate::airfoil::helpers::{
    curve_from_inscribed_circles, extract_edge_sub_curve, inscribed_from_spanning_ray,
    refine_stations, OrientedCircles,
};
use crate::airfoil::{AfParams, InscribedCircle};
use crate::common::points::{dist, mid_point};
use crate::common::{linear_space, Intersection};
use crate::geom2::{rot90, UnitVec2};
use crate::AngleDir::Ccw;
use crate::{Curve2, Point2, Result};
use parry2d_f64::query::Ray;

pub trait EdgeLocation {
    /// Given the airfoil section, the oriented collection of inscribed circles, and a flag
    /// indicating if we're searching for the edge at the front or back of the camber line, this
    /// implementation should return a `Point2` with the actual edge point (where the camber line
    /// meets the section boundary) and the collection of inscribed circles. If the edge point
    /// can't be found, or if the edge point doesn't exist (in the case of an open section), the
    /// function should return `None` for the edge point and return the collection of inscribed
    /// circles without modification.
    ///
    /// This function will be given ownership of the existing vector of circles with the intention
    /// that it *may* make modifications to it, depending on the method. The method will return a
    /// vector of circles that may end up being the same as the input vector, may be different,
    /// or may be the original vector modified.
    ///
    /// # Arguments
    ///
    /// * `section`: the airfoil section curve
    /// * `stations`: the inscribed circles in the airfoil section, already oriented from
    /// leading to trailing edge
    /// * `front`: a flag indicating if we're searching for the edge at the front or the back of
    /// the camber line
    /// * `af_tol`: the tolerance value specified in the airfoil analysis parameters
    ///
    /// returns: Result<(Option<OPoint<f64, Const<2>>>, Vec<InscribedCircle, Global>), Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    fn find_edge(
        &self,
        section: &Curve2,
        stations: Vec<InscribedCircle>,
        front: bool,
        af_tol: f64,
    ) -> Result<(Option<Point2>, Vec<InscribedCircle>)>;
}

/// This struct implements the `EdgeLocation` trait and does not attempt to locate the edge of
/// the airfoil section. It will return `None` for the edge point and the original collection of
/// inscribed circles. Use it in cases where you have a partial airfoil section, and you know that
/// this edge of the section is open.
pub struct OpenEdge {}

impl OpenEdge {
    pub fn new() -> Self {
        OpenEdge {}
    }

    /// Create a new boxed instance of the `OpenEdge` struct.
    pub fn make() -> Box<dyn EdgeLocation> {
        Box::new(OpenEdge::new())
    }
}

impl EdgeLocation for OpenEdge {
    fn find_edge(
        &self,
        section: &Curve2,
        stations: Vec<InscribedCircle>,
        front: bool,
        af_tol: f64,
    ) -> Result<(Option<Point2>, Vec<InscribedCircle>)> {
        Ok((None, stations))
    }
}

/// This struct implements the `EdgeLocation` trait and attempts to locate the edge of the airfoil
/// section by projecting the ray at the end of the camber line until it intersects the section
/// boundary. If the ray does not intersect the boundary, the edge point will be `None`.
pub struct IntersectEdge {}

impl IntersectEdge {
    pub fn new() -> Self {
        IntersectEdge {}
    }

    /// Create a new boxed instance of the `IntersectEdge` struct.
    pub fn make() -> Box<dyn EdgeLocation> {
        Box::new(IntersectEdge::new())
    }
}

impl EdgeLocation for IntersectEdge {
    fn find_edge(
        &self,
        section: &Curve2,
        mut stations: Vec<InscribedCircle>,
        front: bool,
        af_tol: f64,
    ) -> Result<(Option<Point2>, Vec<InscribedCircle>)> {
        let circles = OrientedCircles::new(stations, front);
        let edge_point = circles.intersect_from_end(section);

        Ok((Some(edge_point), circles.take_circles()))
    }
}

/// This struct implements the `EdgeLocation` trait and attempts to locate the edge of the airfoil
/// by first identifying the region of maximum curvature on the portion of the surface beyond the
/// last found station, then tracing its way between the last station and that maximum curvature
/// region.  This will work OK for sections with very smooth curvature and rounded edges, but will
/// struggle on unfiltered/smoothed data and/or sections with sharp edges.
pub struct TraceToMaxCurvature {
    max_tol: Option<f64>,
}

impl TraceToMaxCurvature {
    pub fn new(max_tol: Option<f64>) -> Self {
        TraceToMaxCurvature { max_tol }
    }

    /// Create a new boxed instance of the `TraceToMaxCurvature` struct.  The `max_tol` value is
    /// the fractional tolerance used to determine the acceptable window of maximum tolerance. For
    /// instance, if `max_tol` is 0.01 (1%), then the algorithm will find the point of maximum
    /// curvature and accept the point on the end curve that is closest to the camber line
    /// intersection with the section boundary within 1% of the maximum curvature value.
    ///
    /// If `max_tol` is `None`, then the algorithm will use a default value of 0.005.
    ///
    /// # Arguments
    ///
    /// * `max_tol`: The fractional tolerance used to determine the acceptable window of maximum
    /// curvature.  If `None`, the default value of 0.005 will be used.
    ///
    /// returns: Box<dyn EdgeLocation, Global>
    pub fn make(max_tol: Option<f64>) -> Box<dyn EdgeLocation> {
        Box::new(TraceToMaxCurvature::new(max_tol))
    }
}

impl EdgeLocation for TraceToMaxCurvature {
    fn find_edge(
        &self,
        section: &Curve2,
        mut stations: Vec<InscribedCircle>,
        front: bool,
        af_tol: f64,
    ) -> Result<(Option<Point2>, Vec<InscribedCircle>)> {
        // The last station will vary depending on whether we're looking for the leading or
        // trailing edge
        let station = if front {
            stations
                .first()
                .ok_or("Empty inscribed circles container.")?
        } else {
            stations
                .last()
                .ok_or("Empty inscribed circles container.")?
        };

        let last_radius = station.radius();

        // Now we'll extract a curve that has just the data past the contact points of the last
        // station.  This much reduced dataset will allow for faster search operations. If this
        // fails, it means that the section was open, and we are working on the open portion.
        let edge_curve = extract_edge_sub_curve(&section, &station)
            .ok_or("Failed to extract edge curve on trace to max curvature algorithm.")?;

        // We'll re-assign the stations to the working variable, as we're going to be modifying
        // them in place.
        let mut working_stations = OrientedCircles::new(stations, front);

        // Now we generate the curvature and calculate the maximal plateau
        let curvature = edge_curve.get_curvature_series().abs();
        let (max_x, max_y) = curvature.abs().global_maxima_xy();
        let max = curvature
            .plateau_at_maxima(max_x, max_y * self.max_tol.unwrap_or(0.005))
            .ok_or("Failed to find maximum curvature plateau.")?;

        // Find the point of intersection with the existing camber line and the section boundary.
        let intersection = working_stations.intersect_from_end(&edge_curve);
        let length_i = edge_curve.at_closest_to_point(&intersection).length_along();

        // Now we'll find the point in the interval which is closest to length_i
        let edge_point = edge_curve
            .at_length(max.clamp(length_i))
            .ok_or("Failed to find edge point.")?
            .point();

        // Now that we have the edge point, we'll try to back-fill the camber line between the
        // last inscribed station and the edge point.  We'll use the refining method.
        let mut stack = Vec::new();

        for i in 0..3 {
            let camber_end = working_stations.get_end_curve(last_radius)?;
            let end_point = camber_end.at_back().direction_point();
            let mid = mid_point(&end_point.point, &edge_point);
            let dir = rot90(Ccw) * (edge_point - end_point.point).normalize();
            let test_ray = Ray::new(mid, dir);
            if let Some(spanning_ray) = edge_curve.try_create_spanning_ray(&test_ray) {
                let circle = inscribed_from_spanning_ray(&edge_curve, &spanning_ray, af_tol * 1e-2);
                stack.push(circle);
                refine_stations(
                    &edge_curve,
                    &mut working_stations,
                    &mut stack,
                    af_tol,
                    af_tol * 1e-2,
                );
            } else {
                break;
            }
        }

        Ok((Some(edge_point), working_stations.take_circles()))
    }
}

/// This struct implements the `EdgeLocation` trait and attempts to locate the rounded edge of an
/// airfoil by extending or contracting the end of the camber line until it intersects the section
/// boundary at a location where the section tangent is perpendicular to end of the camber line.
/// This method is useful if you know that the camber line should meet the section boundary at a
/// right angle, which is true of most rounded edge airfoils.
pub struct ConvergeTangentEdge {
    tol: Option<f64>,
}

impl ConvergeTangentEdge {
    pub fn new(tol: Option<f64>) -> Self {
        ConvergeTangentEdge { tol }
    }

    /// Create a new boxed instance of the `ConvergeTangentEdge` struct.  The tolerance value is
    /// the allowable lateral distance between the direction of the end of the camber line and the
    /// tangent point on the section boundary perpendicular to the end of the camber line.
    ///
    /// The algorithm will end the camber line at the closest point to the original last station
    /// that is within the tolerance value of the tangent point.  This accounts for the typical
    /// increase in camber line curvature which happens as one approaches the section boundary. The
    /// closer to the original last inscribed circle, the less overall curvature will be included
    /// in the final camber line.
    ///
    /// If no tolerance value is specified, the default value of 1% of the last inscribed circle's
    /// radius will be used.
    pub fn make(tol: Option<f64>) -> Box<dyn EdgeLocation> {
        Box::new(ConvergeTangentEdge::new(tol))
    }
}

impl EdgeLocation for ConvergeTangentEdge {
    fn find_edge(
        &self,
        section: &Curve2,
        stations: Vec<InscribedCircle>,
        front: bool,
        af_tol: f64,
    ) -> Result<(Option<Point2>, Vec<InscribedCircle>)> {
        // Capture the last inscribed circle
        let temp_circles = OrientedCircles::new(stations, front);
        let c0 = temp_circles
            .last()
            .ok_or("Empty inscribed circles container.")?
            .circle
            .clone();

        let check_tol = self.tol.unwrap_or(c0.r() * 0.01);

        // Now we'll use the maximum curvature method to update the camber line
        let max_curvature = TraceToMaxCurvature::new(Some(0.01));
        let (_, stations) =
            max_curvature.find_edge(section, temp_circles.take_circles(), front, af_tol)?;

        let mut working_stations = OrientedCircles::new(stations, front);
        let edge_curve = extract_edge_sub_curve(&section, &working_stations.last().unwrap())
            .ok_or("Failed to extract edge curve on converging tangent edge algorithm.")?;

        // We're going to sweep the edge curve from back about 2 c0 diameters to the front of the
        // camber line, projecting a ray from the edge curve's interpolated direction and finding
        // the max section point in that direction.  The error will be the out-of-direction
        // distance from the ray to the located point.
        //
        // Ultimately, we're going to look for the point closest to c0 which is within the
        // allowed tolerance of the ray.

        let oriented_camber = working_stations.get_full_curve()?;
        let x0 = oriented_camber
            .at_closest_to_point(&c0.center)
            .length_along();

        let steps = 1000;
        let x_start = x0 - 4.0 * c0.r();
        let space = linear_space(x_start, oriented_camber.length(), steps);
        let mut measured = Vec::new();

        for x in space.values().iter() {
            let camber_point = oriented_camber.at_length(*x).unwrap();
            let camber_dir = camber_point.interpolated_direction_point();

            let max_point = edge_curve
                .max_point_in_direction(&camber_dir.normal)
                .ok_or("Failed to find max point in direction, did the curve end up empty?.")?;

            let error = camber_dir.planar_distance(&max_point);
            if error < check_tol {
                measured.push(((x - x0).abs(), camber_dir, max_point));
            }
        }

        // Now we'll find the closest point to c0 that is within the tolerance
        if let Some((_, dir, point)) = measured
            .iter()
            .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap())
        {
            // We now need to remove all the stations beyond the point we found
            working_stations.discard_sections_beyond_point(&dir.point);
            Ok((Some(*point), working_stations.take_circles()))
        } else {
            Ok((None, working_stations.take_circles()))
        }
    }
}
