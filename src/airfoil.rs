//! This module contains structures and algorithms for performing dimensional analysis of
//! airfoil sections, such as calculating the camber line, identifying the leading and trailing
//! edges, computing angles, thicknesses, and other properties.

mod camber;
mod edges;
pub mod helpers;
mod inscribed_circle;
mod orientation;

use crate::{Arc2, Curve2, Point2, Result, SurfacePoint2};

use crate::common::points::dist;
use crate::geom2::hull::convex_hull_2d;
pub use camber::extract_camber_line;
pub use edges::{
    ConstRadiusEdge, ConvergeTangentEdge, EdgeLocation, FitRadiusEdge, IntersectEdge, OpenEdge,
    OpenIntersectGap, TraceToMaxCurvature,
};
pub use inscribed_circle::InscribedCircle;
pub use orientation::{CamberOrientation, DirectionFwd, TMaxFwd};
use serde::{Deserialize, Serialize};

/// This structure contains the parameters used in the airfoil analysis algorithms.  It specifies
/// the minimum tolerance value used in many parts of the analysis, as well as the methods for
/// detecting the orientation of the leading edge, and the leading and trailing edges themselves.
pub struct AfParams {
    /// The minimum tolerance value, used in many parts of the analysis.  Generally speaking, the
    /// various algorithms will attempt to iteratively refine results until the error/difference
    /// falls below this value.
    pub tol: f64,

    /// The method for trying to detect the orientation of the leading edge on the airfoil.
    pub orient: Box<dyn CamberOrientation>,

    /// The method for trying to detect the leading edge on the airfoil.
    pub leading: Box<dyn EdgeLocation>,

    /// The method for trying to detect the trailing edge on the airfoil.
    pub trailing: Box<dyn EdgeLocation>,
}

impl AfParams {
    /// Create a new set of airfoil analysis parameters with the specified tolerance value and
    /// other algorithm selections.
    ///
    /// # Arguments
    ///
    /// * `tol`: the minimum tolerance value used in many parts of the analysis, generally used to
    /// refine results until the error/difference falls below this value.
    /// * `orient`: the method for trying to detect the orientation of the leading edge
    /// * `leading`: the method for trying to detect the leading edge
    /// * `trailing`: the method for trying to detect the trailing edge
    ///
    /// returns: AfParams
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn new(
        tol: f64,
        orient: Box<dyn CamberOrientation>,
        leading: Box<dyn EdgeLocation>,
        trailing: Box<dyn EdgeLocation>,
    ) -> Self {
        AfParams {
            tol,
            orient,
            leading,
            trailing,
        }
    }
}

/// This enumeration represents the possible edge geometries that can be detected on an airfoil
/// by the analysis methods, and is used to return information located by the edge detection
/// methods.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum EdgeGeometry {
    /// No special edge geometry was found
    None,

    /// The edge has a constant radius region represented by an arc
    Arc(Arc2),
}

/// An airfoil edge is a generic construct used to represent the leading and trailing edges of an
/// airfoil. When an edge is detected, it consists of a point and a geometry.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AirfoilEdge {
    pub point: Point2,
    pub geometry: EdgeGeometry,
}

impl AirfoilEdge {
    fn new(point: Point2, geometry: EdgeGeometry) -> Self {
        AirfoilEdge { point, geometry }
    }

    pub fn point_only(point: Point2) -> Self {
        AirfoilEdge::new(point, EdgeGeometry::None)
    }
}

/// This struct contains the results of a geometric analysis of an airfoil section.  It includes
/// the camber line, optional leading and trailing edge information, and other properties.
#[derive(Clone)]
pub struct AirfoilGeometry {
    /// The leading edge point of the airfoil section, if it was detected.
    pub leading_edge: Option<AirfoilEdge>,

    /// The trailing edge point of the airfoil section, if it was detected.
    pub trailing_edge: Option<AirfoilEdge>,

    /// A vector of inscribed circles in order from leading edge to trailing edge.
    pub stations: Vec<InscribedCircle>,

    /// The known portion of the airfoil section camber line, represented as a curve oriented from
    /// the leading edge to the trailing edge. If the leading/trailing edges are known, the
    /// first/last points of the curve will be the leading/trailing edge points, respectively.
    /// Otherwise, the curve will stop at the first/last inscribed circle.
    pub camber: Curve2,
}

impl AirfoilGeometry {
    fn new(
        leading_edge: Option<AirfoilEdge>,
        trailing_edge: Option<AirfoilEdge>,
        stations: Vec<InscribedCircle>,
        camber: Curve2,
    ) -> Self {
        AirfoilGeometry {
            leading_edge,
            trailing_edge,
            stations,
            camber,
        }
    }

    /// Find the inscribed circle with the maximum radius, which is typically a circle near the
    /// center of the airfoil section.
    pub fn find_tmax(&self) -> &InscribedCircle {
        self.stations
            .iter()
            .max_by(|a, b| a.radius().partial_cmp(&b.radius()).unwrap())
            .unwrap()
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct ChordLine {
    pub le: Point2,
    pub te: Point2,
}

impl ChordLine {
    pub fn new(le: Point2, te: Point2) -> Self {
        ChordLine { le, te }
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct CaliperChord {
    pub chord: ChordLine,
    pub tangent: ChordLine,
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
pub fn caliper_chord_line(section: &Curve2, camber: &Curve2) -> Result<CaliperChord> {
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

    let chord_line = ChordLine::new(le.1, te.1);
    let tangent_line = ChordLine::new(chord.projection(&le.1), chord.projection(&te.1));

    Ok(CaliperChord {
        chord: chord_line,
        tangent: tangent_line,
    })
}

/// Perform a geometric analysis of an airfoil section, extracting the camber line, leading and
/// trailing edges, and other properties. Geometric airfoil section analysis is centered around the
/// MCL (mean camber line) extraction through the inscribed circle method, and detects features of
/// the airfoil based solely on the geometry of the section.  It is suitable for use with very
/// clean airfoil section data, especially nominal geometry such as that from CAD or sections
/// generated mathematically.
///
/// It is less suitable for use with measured data, which has noise that can "poison" the geometry
/// enough that features will not be detected as expected. For measured data, especially measured
/// data which can have noise or large deviations from ideal geometry (such as damage, wear, or
/// significant warping), an analysis using a nominal reference airfoil is recommended.
///
/// # Arguments
///
/// * `section`: a `Curve2` representing the airfoil section geometry. This curve should be closed
/// if the section is intended to be closed. No specific orientation is required.
/// * `params`: the `AfParams` structure containing the parameters used in the analysis. Select the
/// appropriate values for the tolerance, orientation, and edge detection methods with care.
///
/// returns: Result<AnalyzedAirfoil, Box<dyn Error, Global>>
///
/// # Examples
///
/// ```
///
/// ```
pub fn analyze_airfoil_geometry(section: &Curve2, params: &AfParams) -> Result<AirfoilGeometry> {
    // Calculate the hull, we will need this for the inscribed circle method and the tangency
    // line.
    let hull = section
        .make_hull()
        .ok_or("Failed to calculate the hull of the airfoil section")?;

    // Compute the mean camber line using the inscribed circle method
    let stations = extract_camber_line(section, &hull, Some(params.tol))
        .map_err(|e| format!("Error during initial camber line extraction: {e}"))?;

    // Orient the camber line
    let stations = params
        .orient
        .orient_camber_line(section, stations)
        .map_err(|e| format!("Error orienting the initial camber line: {e}"))?;

    // Find the leading and trailing edges
    let (leading_edge, stations) = params
        .leading
        .find_edge(section, stations, true, params.tol)
        .map_err(|e| format!("Error finding the leading edge: {e}"))?;

    let (trailing_edge, stations) = params
        .trailing
        .find_edge(section, stations, false, params.tol)
        .map_err(|e| format!("Error finding the trailing edge: {e}"))?;

    // Create the camber curve
    let mut camber_points = stations.iter().map(|c| c.circle.center).collect::<Vec<_>>();
    if let Some(leading) = &leading_edge {
        camber_points.insert(0, leading.point);
    }
    if let Some(trailing) = &trailing_edge {
        camber_points.push(trailing.point);
    }
    let camber = Curve2::from_points(&camber_points, params.tol, false)
        .map_err(|e| format!("Error creating the final camber curve: {e}"))?;

    Ok(AirfoilGeometry::new(
        leading_edge,
        trailing_edge,
        stations,
        camber,
    ))
}
