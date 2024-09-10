//! This module contains structures and algorithms for performing dimensional analysis of
//! airfoil sections, such as calculating the camber line, identifying the leading and trailing
//! edges, computing angles, thicknesses, and other properties.

mod camber;
mod inscribed_circle;

use crate::{Curve2, Point2, Result, Vector2};

use crate::airfoil::camber::extract_camber_line;
pub use inscribed_circle::InscribedCircle;

/// Enum specifying the method for trying to detect the orientation of the leading edge on the
/// airfoil, with the default being to use the point of max thickness.
#[derive(Debug, Clone, Copy)]
pub enum Orientation {
    /// This method is based on the principle that for most subsonic airfoils the point of max
    /// thickness is closer to the leading edge than to the trailing edge along the camber line.
    TmaxFwd,

    /// This method takes a direction vector and orients the leading edge such that the vector from
    /// the back to the front of the camber line has a positive dot product with that vector. Use
    /// this method when you have a specific direction in mind for the leading edge.
    Direction(Vector2),
}

impl Default for Orientation {
    fn default() -> Self {
        Orientation::TmaxFwd
    }
}

/// Enum specifying the method for trying to extract an edge (either leading or trailing) from the
/// airfoil section data.  The default is to use the auto-detection method, but it is always better
/// to specify the edge detection method explicitly if you know what the edge should look like.
#[derive(Debug, Clone, Copy)]
pub enum Edge {
    /// Attempts to automatically detect the edge based on the airfoil section data.  Use this
    /// method only if you are unsure of the edge type.
    Auto,

    /// Attempts to detect a rounded leading edge with or without a specified radius.  It can take
    /// two optional parameters: a boolean flag indicating whether to adjust the edge point based
    /// on the point of furthest extension, and a floating point value indicating the expected
    /// radius of the leading edge.
    Round((Option<bool>, Option<f64>)),

    /// Attempts to detect a square edge (one having two points of maximum curvature). This is
    /// sometimes common on trailing edges.
    Square,

    /// Attempts to detect a sharp edge (one having a single point of maximum curvature). This is
    /// uncommon on subsonic leading edges, but is common on some trailing edges.
    Sharp,

    /// Specify this method if only a partial airfoil section has been provided and this side is
    /// not a closed section. For instance, if only leading edge data is provided, select `Open`
    /// for the trailing edge to indicate that it is not present.  This will cause the edge to
    /// be skipped.
    Open,
}

impl Default for Edge {
    fn default() -> Self {
        Edge::Auto
    }
}

/// This structure contains the parameters used in the airfoil analysis algorithms.  It specifies
/// the minimum tolerance value used in many parts of the analysis, as well as the methods for
/// detecting the orientation of the leading edge, and the leading and trailing edges themselves.
#[derive(Debug, Clone)]
pub struct AfParams {
    /// The minimum tolerance value, used in many parts of the analysis.  Generally speaking, the
    /// various algorithms will attempt to iteratively refine results until the error/difference
    /// falls below this value.
    pub tol: f64,

    /// The method for trying to detect the orientation of the leading edge on the airfoil.
    pub orient: Orientation,

    /// The method for trying to detect the leading edge on the airfoil.
    pub leading: Edge,

    /// The method for trying to detect the trailing edge on the airfoil.
    pub trailing: Edge,
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
    pub fn new(tol: f64, orient: Orientation, leading: Edge, trailing: Edge) -> Self {
        AfParams {
            tol,
            orient,
            leading,
            trailing,
        }
    }
}

impl Default for AfParams {
    fn default() -> Self {
        AfParams::new(
            1e-4,
            Orientation::default(),
            Edge::default(),
            Edge::default(),
        )
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

    /// The parameters used in the analysis.
    pub af_params: AfParams,
}

impl AirfoilGeometry {
    fn new(
        leading_edge: Option<Point2>,
        trailing_edge: Option<Point2>,
        stations: Vec<InscribedCircle>,
        camber: Curve2,
        af_params: AfParams,
    ) -> Self {
        AirfoilGeometry {
            leading_edge,
            trailing_edge,
            stations,
            camber,
            af_params,
        }
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

    // Compute the mean camber line using the inscribed circle method, create the camber curve
    let stations = extract_camber_line(section, &hull, Some(params.tol))?;
    let camber_points = stations.iter().map(|c| c.circle.center).collect::<Vec<_>>();
    let camber = Curve2::from_points(&camber_points, params.tol, false)?;

    // Now we'll attempt to detect the leading and trailing edges.

    Ok(AirfoilGeometry::new(
        None,
        None,
        stations,
        camber,
        params.clone(),
    ))
}
