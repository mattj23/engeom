//! This module contains tools to work with the leading and trailing edges of the airfoil section.

use crate::airfoil::helpers::{
    curve_from_inscribed_circles, extract_edge_sub_curve, inscribed_from_spanning_ray,
    refine_stations, OrientedCircles,
};
use crate::airfoil::{AfParams, InscribedCircle};
use crate::common::points::{dist, mid_point};
use crate::common::Intersection;
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
    /// the optional angular tolerance value used to determine if the tangent of the section
    /// boundary is perpendicular to the end of the camber line.
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
        // We'll orient the camber curve so that we're working at the back end in either the
        // leading or trailing edge cases
        let camber = if front {
            curve_from_inscribed_circles(&stations, section.tol())?.reversed()
        } else {
            curve_from_inscribed_circles(&stations, section.tol())?
        };

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

        // Now we'll extract a curve that has just the data past the contact points of the last
        // station.  This much reduced dataset will allow for faster search operations. If this
        // fails, it means that the section was open, and we are working on the open portion.
        let edge_curve = extract_edge_sub_curve(&section, &station)
            .ok_or("Failed to extract edge curve on converging tangent edge algorithm.")?;

        // We'll re-assign the stations to the working variable, as we're going to be modifying
        // them in place.
        let mut working_stations = OrientedCircles::new(stations, front);
        let mut working_points = camber.points()[camber.points().len() - 5..].to_vec();
        let mut working_camber = Curve2::from_points(&working_points, 1e-4, false)?;

        // We're going to need to try to advance the end of the camber line closer into the edge
        // than the last station.  We have some advantages at this point compared to the general
        // camber line extraction.
        // - We know that the edge we're approaching should be generally rounded
        // - We know that the curvature between here and the end of the airfoil is relatively low

        // First, we'll try to create up to eight additional inscribed circles
        let mut count = 0;
        while count < 30 {
            let working_point = working_camber.at_back().direction_point();

            let d = *edge_curve
                .intersection(&working_point)
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .ok_or("Failed to find intersection with edge curve.")?;

            let test_ray = Ray::new(
                working_point.at_distance(d * 0.1),
                (rot90(Ccw) * working_point.normal).into_inner(),
            );

            if let Some(spanning_ray) = edge_curve.try_create_spanning_ray(&test_ray) {
                let circle = inscribed_from_spanning_ray(&edge_curve, &spanning_ray, af_tol * 1e-2);

                // We'll look at the angle formed between the contact points and the center of
                // the circle. When that angle drops to below 90 degrees, we'll stop.
                let v0 = circle.upper - circle.center();
                let v1 = circle.lower - circle.center();
                let angle = v0.angle(&v1);
                if angle < std::f64::consts::FRAC_PI_2 * 0.25 {
                    break;
                }

                working_points.push(circle.center());
                working_stations.push(circle);
                working_camber = Curve2::from_points(&working_points, 1e-4, false)?;
            } else {
                break;
            }

            count += 1;
        }

        Ok((None, working_stations.take_circles()))
    }
}
