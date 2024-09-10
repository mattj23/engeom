//! This module contains tools to work with the leading and trailing edges of the airfoil section.

use crate::airfoil::helpers::curve_from_inscribed_circles;
use crate::airfoil::{AfParams, InscribedCircle};
use crate::common::Intersection;
use crate::{Curve2, Point2, Result};

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
        // We construct the camber curve and get the direction point at either the front or the
        // back of the curve, depending on the value of the `front` flag.  If we take the front
        // of the curve we need to reverse the direction point to get the ray that extends forward
        // from the camber line.
        let camber = curve_from_inscribed_circles(&stations, 1e-4)?;
        let ray = if front {
            camber.at_front().direction_point().reversed()
        } else {
            camber.at_back().direction_point()
        };

        // Find all the intersections of the ray with the airfoil section, which will be returned
        // as distances (parameters) along the ray.
        let ts = section.intersection(&ray);

        // If there is no intersection, we return `None` for the edge point.  If there is one or
        // more we will take the intersection with the largest parameter value, which will be the
        // one furthest along the ray.
        if ts.is_empty() {
            Ok((None, stations))
        } else {
            let t = ts
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .ok_or("Failed to find maximum intersection parameter.")?;
            Ok((Some(ray.at_distance(*t)), stations))
        }
    }
}
