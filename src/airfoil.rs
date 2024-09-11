//! This module contains structures and algorithms for performing dimensional analysis of
//! airfoil sections, such as calculating the camber line, identifying the leading and trailing
//! edges, computing angles, thicknesses, and other properties.

mod camber;
mod edges;
pub mod helpers;
mod inscribed_circle;
mod orientation;

use crate::{Curve2, Point2, Result};

use crate::airfoil::camber::extract_camber_line;
pub use edges::{ConvergeTangentEdge, EdgeLocation, IntersectEdge, OpenEdge, TraceToMaxCurvature};
pub use inscribed_circle::InscribedCircle;
pub use orientation::{CamberOrientation, TMaxFwd};

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

/// This struct contains the results of a geometric analysis of an airfoil section.  It includes
/// the camber line, optional leading and trailing edge information, and other properties.
pub struct AirfoilGeometry {
    /// The leading edge point of the airfoil section, if it was detected.
    pub leading_edge: Option<Point2>,

    /// The trailing edge point of the airfoil section, if it was detected.
    pub trailing_edge: Option<Point2>,

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
        leading_edge: Option<Point2>,
        trailing_edge: Option<Point2>,
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
    let stations = extract_camber_line(section, &hull, Some(params.tol))?;

    // Orient the camber line
    let stations = params.orient.orient_camber_line(section, stations)?;

    // Find the leading and trailing edges
    let (leading_edge, stations) = params
        .leading
        .find_edge(section, stations, true, params.tol)?;
    let (trailing_edge, stations) = params
        .trailing
        .find_edge(section, stations, false, params.tol)?;

    // Create the camber curve
    let mut camber_points = stations.iter().map(|c| c.circle.center).collect::<Vec<_>>();
    if let Some(leading) = leading_edge {
        camber_points.insert(0, leading);
    }
    if let Some(trailing) = trailing_edge {
        camber_points.push(trailing);
    }
    let camber = Curve2::from_points(&camber_points, params.tol, false)?;

    Ok(AirfoilGeometry::new(
        leading_edge,
        trailing_edge,
        stations,
        camber,
    ))
}
