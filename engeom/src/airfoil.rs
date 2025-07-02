//! This module contains structures and algorithms for performing dimensional analysis of
//! airfoil sections, such as calculating the camber line, identifying the leading and trailing
//! edges, computing angles, thicknesses, and other properties.

mod camber;
mod edges;
pub mod helpers;
mod inscribed_circle;
mod orientation;

use crate::{Arc2, Circle2, Curve2, Point2, Result, SurfacePoint2, UnitVec2, Vector2};

use crate::airfoil::camber::camber_detect_upper_dir;
use crate::common::Intersection;
use crate::common::points::dist;
use crate::geom2::hull::convex_hull_2d;
use crate::metrology::Distance2;
use crate::stats::compute_mean;
pub use camber::extract_camber_line;
pub use edges::{
    ConstRadiusEdge, ConvergeTangentEdge, FitRadiusEdge, IntersectEdge, OpenEdge, OpenIntersectGap,
    RansacRadiusEdge, TraceToMaxCurvature,
};
pub use inscribed_circle::InscribedCircle;
pub use orientation::{DirectionFwd, TMaxFwd};
use serde::{Deserialize, Serialize};
//=================================================================================================
// Airfoil analysis parameters and interfaces
//-------------------------------------------------------------------------------------------------
// This section has enumerations and traits which specify how algorithms for different analysis
// methods will be identified and/or given to downstream code.
//=================================================================================================

#[derive(Serialize, Deserialize, Clone, Copy)]
pub enum AfGage {
    OnCamber(f64),
    Radius(f64),
}

/// This enumeration represents the possible methods for detecting which of the two surfaces of
/// an airfoil is the upper (suction/convex) surface and which is the lower (pressure/concave)
/// surface.
#[derive(Serialize, Deserialize, Clone)]
pub enum FaceOrient {
    Detect,
    UpperDir(Vector2),
}

/// This trait defines an interface to perform orientation of the camber line of an airfoil
/// section, specifically referring to its order and its relationship to the leading edge. A
/// camber line begins at the leading edge and ends at the trailing edge.
pub trait CamberOrient {
    /// Orient the camber line of an airfoil section based on the method specified by the
    /// implementation. The camber line at this stage is a series of inscribed circles whose
    /// adjacency in space is coupled to their ordering in the container. However, the orientation
    /// (whether the first circle is at the leading edge or the trailing edge) is not yet known.
    ///
    /// This method will return a new container of inscribed circles with the camber line oriented
    /// so that the first circle is closest to the leading edge and the last circle is closest to
    /// the trailing edge. The order of the circles will be preserved.
    ///
    /// This method will take ownership of the input container and return a new container with the
    /// circles in the correct order.
    ///
    /// # Arguments
    ///
    /// * `section`: the airfoil section curve used to generate the inscribed circles
    /// * `stations`: the inscribed circles in the airfoil section
    ///
    /// returns: Vec<InscribedCircle, Global>
    fn orient_camber_line(
        &self,
        section: &Curve2,
        stations: Vec<InscribedCircle>,
    ) -> Result<Vec<InscribedCircle>>;
}

/// This trait defines an interface to locate the edge (leading or trailing) of an airfoil
/// cross-section to identify the specific point where the camber line meets the section boundary
/// and any additional geometry information that may be present.
pub trait EdgeLocate {
    /// Given the airfoil section, the oriented collection of inscribed circles, and a flag
    /// indicating if we're searching for the edge at the front or back of the camber line, this
    /// implementation should return an `AirfoilEdge` with the actual edge point (where the camber
    /// line meets the section boundary), optional geometry information where it exists, and the
    /// collection of inscribed circles. If the edge point can't be found, or if the edge point
    /// doesn't exist (in the case of an open section), the function should return `None` for the
    /// edge point and return the collection of inscribed circles without modification.
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
    ///   leading to trailing edge
    /// * `front`: a flag indicating if we're searching for the edge at the front or the back of
    ///   the camber line
    /// * `af_tol`: the tolerance value specified in the airfoil analysis parameters
    ///
    /// returns: Result<(Option<AirfoilEdge>, Vec<InscribedCircle, Global>), Box<dyn Error, Global>>
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
    ) -> Result<(Option<AirfoilEdge>, Vec<InscribedCircle>)>;
}

//=================================================================================================
// Airfoil result structures
//-------------------------------------------------------------------------------------------------
// This section has structures and enumerations which represent a common set of outputs for
// different analysis methods that perform the same function.
//=================================================================================================

/// This enumeration represents the possible edge geometries that can be detected on an airfoil
/// by the analysis methods, and is used to return information located by the edge detection
/// methods.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum EdgeGeometry {
    /// No special geometry was found, but the section is known to be closed at this edge.
    Closed,

    /// The section is known to be open at this edge
    Open,

    /// The edge is closed and has a constant radius region represented by an arc
    Arc(Arc2),
}

/// An airfoil edge is a generic construct used to represent the leading and trailing edges of an
/// airfoil. When an edge is detected, it consists of a point and a geometry.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct AirfoilEdge {
    /// The leading or trailing edge point of the airfoil section. This point should lie on both
    /// the camber line and the airfoil section boundary, unless the edge is open.
    pub point: Point2,

    /// The geometry of the edge, if any special geometry was detected.  This can be used to provide
    /// additional information about the edge, such as a constant radius region.
    pub geometry: EdgeGeometry,
}

impl AirfoilEdge {
    fn new(point: Point2, geometry: EdgeGeometry) -> Self {
        AirfoilEdge { point, geometry }
    }

    fn open(point: Point2) -> Self {
        AirfoilEdge::new(point, EdgeGeometry::Open)
    }

    fn closed_only(point: Point2) -> Self {
        AirfoilEdge::new(point, EdgeGeometry::Closed)
    }

    fn arc(point: Point2, arc: Arc2) -> Self {
        AirfoilEdge::new(point, EdgeGeometry::Arc(arc))
    }
}

/// This struct contains the results of a geometric analysis of an airfoil section.  It includes
/// the camber line, optional leading and trailing edge information, and other properties.
#[derive(Clone, Serialize, Deserialize)]
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

    /// The upper (suction) surface of the airfoil section, represented as a curve oriented in the
    /// same winding order as the original section. The first point may be at either the leading or
    /// trailing edge based on the coordinate system of the original section.
    pub upper: Option<Curve2>,

    /// The lower (pressure) surface of the airfoil section, represented as a curve oriented in the
    /// same winding order as the original section. The first point may be at either the leading or
    /// trailing edge based on the coordinate system of the original section.
    pub lower: Option<Curve2>,

    /// The tolerance value used in the analysis
    pub core_tol: f64,

    /// The tolerance value of curves used in this analysis
    pub curve_tol: f64,
}

impl AirfoilGeometry {
    /// Attempt to extract the airfoil geometry from the given section using only the geometric data
    /// of the section.
    ///
    /// The mean camber line will be extracted iteratively using the inscribed circle method with
    /// no prior knowledge about the section shape. The orientation of the leading edge direction,
    /// the pressure/suction surfaces, and the locations of the leading and trailing edges must be
    /// found with methods that don't need additional information.
    ///
    /// This geometry-only analysis is suitable for use with very clean airfoil section data, such
    /// as nominal geometry from CAD or sections generated mathematically. It will work best if
    /// you can specify good, specific methods for feature identification based on knowledge you
    /// have of the section.
    ///
    /// It is less suitable for use with measured data, which can have noise that can "poison" the
    /// geometry enough that features may not be detected as expected. For measured data, especially
    /// data with noise or large deviations from ideal geometry (such as damage, wear, or
    /// significant warping), an analysis using a nominal reference airfoil is recommended instead.
    ///
    /// # Arguments
    ///
    /// * `section`: The `Curve2` representing the airfoil section geometry. This curve should be
    ///   closed if the section is intended to be closed. No specific orientation is required.
    /// * `core_tol`: a tolerance value used in many parts of the analysis, generally used to refine
    ///   results until the error/difference falls below this value.
    /// * `camber_orient`: the method for trying to detect which side of the mean camber line is
    ///   at the leading edge and which is at the trailing edge.
    /// * `leading`: the algorithm for trying to detect the leading edge of the airfoil section
    /// * `trailing`: the algorithm for trying to detect the trailing edge of the airfoil section
    /// * `face_orient`: the method to use to try to detect which surface of the airfoil is the
    ///   upper (suction/convex) surface and which is the lower (pressure/concave) surface.
    ///
    /// returns: Result<AirfoilGeometry, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn try_analyze(
        section: &Curve2,
        core_tol: f64,
        camber_orient: Box<dyn CamberOrient>,
        leading: Box<dyn EdgeLocate>,
        trailing: Box<dyn EdgeLocate>,
        face_orient: FaceOrient,
    ) -> Result<Self> {
        // Calculate the hull, we will need this for the inscribed circle method and the tangency
        // line.
        let hull = section
            .make_hull()
            .ok_or("Failed to calculate the hull of the airfoil section")?;

        // Compute the mean camber line using the inscribed circle method
        let stations = extract_camber_line(section, &hull, Some(core_tol))
            .map_err(|e| format!("Error during initial camber line extraction: {e}"))?;

        // Orient the camber line
        let stations = camber_orient
            .orient_camber_line(section, stations)
            .map_err(|e| format!("Error orienting the initial camber line: {e}"))?;

        // Find the leading and trailing edges
        let (leading_edge, stations) = leading
            .find_edge(section, stations, true, core_tol)
            .map_err(|e| format!("Error finding the leading edge: {e}"))?;

        let (trailing_edge, stations) = trailing
            .find_edge(section, stations, false, core_tol)
            .map_err(|e| format!("Error finding the trailing edge: {e}"))?;

        // Create the camber curve
        let mut camber_points = stations.iter().map(|c| c.circle.center).collect::<Vec<_>>();
        if let Some(leading) = &leading_edge {
            camber_points.insert(0, leading.point);
        }
        if let Some(trailing) = &trailing_edge {
            camber_points.push(trailing.point);
        }

        let camber = Curve2::from_points(&camber_points, section.tol(), false)
            .map_err(|e| format!("Error creating the final camber curve: {e}"))?;

        let upper_dir = match face_orient {
            FaceOrient::Detect => camber_detect_upper_dir(&camber)?,
            FaceOrient::UpperDir(dir) => UnitVec2::new_normalize(dir),
        };

        // Split the airfoil section into upper and lower curves. For this to work we need both
        // the leading and the trailing edge to have been located.
        let pieces = if let (Some(le), Some(te)) = (&leading_edge, &trailing_edge) {
            // One of the two may be open
            match (&le.geometry, &te.geometry) {
                (EdgeGeometry::Open, _) => {
                    // Leading edge is open, we have a trailing edge
                    let l = section.at_closest_to_point(&te.point).length_along();
                    Some(section.split_open_at_length(l)?)
                }
                (_, EdgeGeometry::Open) => {
                    // Trailing edge is open, we have a leading edge
                    let l = section.at_closest_to_point(&le.point).length_along();
                    Some(section.split_open_at_length(l)?)
                }
                (_, _) => {
                    let l0 = section.at_closest_to_point(&le.point).length_along();
                    let l1 = section.at_closest_to_point(&te.point).length_along();
                    Some(section.split_closed_at_lengths(l0, l1)?)
                }
            }
        } else {
            // In this case, we don't have enough information to perform the split
            None
        };

        let (upper, lower) = if let Some((a, b)) = pieces {
            let camber_mid = camber
                .at_fraction(0.5)
                .ok_or("Failed to find camber midpoint, this shouldn't happen")?
                .point();
            let test_point = SurfacePoint2::new(camber_mid, upper_dir);
            order_faces(a, b, test_point)?
        } else {
            (None, None)
        };

        Ok(AirfoilGeometry {
            leading_edge,
            trailing_edge,
            stations,
            camber,
            upper,
            lower,
            core_tol,
            curve_tol: section.tol(),
        })
    }

    /// Find the inscribed circle with the maximum radius, which is typically a circle near the
    /// center of the airfoil section.
    pub fn find_tmax(&self) -> &InscribedCircle {
        self.stations
            .iter()
            .max_by(|a, b| a.radius().partial_cmp(&b.radius()).unwrap())
            .unwrap()
    }

    /// Take a thickness measurement of the airfoil at the given gage location.
    ///
    /// # Arguments
    ///
    /// * `gage`:
    ///
    /// returns: Result<Distance2, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn get_thickness(&self, gage: AfGage) -> Result<Distance2> {
        let (upper, lower) = match gage {
            AfGage::OnCamber(x) => {
                let l = if x < 0.0 { self.camber.length() + x } else { x };
                let sp = self
                    .camber
                    .at_length(l)
                    .ok_or("Invalid camber length")?
                    .surface_point();
                let ut = self
                    .upper
                    .as_ref()
                    .ok_or("Upper surface not found")?
                    .intersection(&sp);
                if ut.is_empty() {
                    return Err("Failed to find upper surface intersection".into());
                }

                let lt = self
                    .lower
                    .as_ref()
                    .ok_or("Lower surface not found")?
                    .intersection(&sp);
                if lt.is_empty() {
                    return Err("Failed to find lower surface intersection".into());
                }

                // Take the minimum absolute value of each set
                let u = ut
                    .iter()
                    .min_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap())
                    .ok_or("Failed to find upper surface intersection")?;

                let l = lt
                    .iter()
                    .min_by(|a, b| a.abs().partial_cmp(&b.abs()).unwrap())
                    .ok_or("Failed to find lower surface intersection")?;

                (sp.at_distance(*u), sp.at_distance(*l))
            }
            AfGage::Radius(r) => {
                let c = if r < 0.0 {
                    Circle2::from_point(
                        self.trailing_edge
                            .as_ref()
                            .ok_or("Trailing edge not found")?
                            .point,
                        -r,
                    )
                } else {
                    Circle2::from_point(
                        self.leading_edge
                            .as_ref()
                            .ok_or("Leading edge not found")?
                            .point,
                        r,
                    )
                };

                let u = self
                    .upper
                    .as_ref()
                    .ok_or("Upper surface not found")?
                    .intersection(&c);
                let l = self
                    .lower
                    .as_ref()
                    .ok_or("Lower surface not found")?
                    .intersection(&c);

                if u.is_empty() || l.is_empty() {
                    return Err("Failed to find upper or lower surface intersection".into());
                }

                // Find the point which is the closest to the camber line
                let u = u
                    .iter()
                    .min_by(|a, b| {
                        self.camber
                            .dist_to_point(a)
                            .partial_cmp(&self.camber.dist_to_point(b))
                            .unwrap()
                    })
                    .ok_or("Failed to find upper surface intersection")?;
                let l = l
                    .iter()
                    .min_by(|a, b| {
                        self.camber
                            .dist_to_point(a)
                            .partial_cmp(&self.camber.dist_to_point(b))
                            .unwrap()
                    })
                    .ok_or("Failed to find lower surface intersection")?;

                (*u, *l)
            }
        };

        Ok(Distance2::new(lower, upper, None))
    }

    pub fn get_thickness_max(&self) -> Result<Distance2> {
        let tmax = self.find_tmax();
        let (upper, lower) = self.order_points(&tmax.contact_neg, &tmax.contact_pos)?;
        Ok(Distance2::new(lower, upper, None))
    }

    /// Order the points so that the first point returned is on the upper surface and the second
    /// point returned is on the lower surface. If the upper and lower surfaces are not known,
    /// this will return an error.
    ///
    /// # Arguments
    ///
    /// * `a`: the first point to order
    /// * `b`: the second point to order
    ///
    /// returns: Result<(OPoint<f64, Const<2>>, OPoint<f64, Const<2>>), Box<dyn Error, Global>>
    fn order_points(&self, a: &Point2, b: &Point2) -> Result<(Point2, Point2)> {
        let upper_dist = self
            .upper
            .as_ref()
            .ok_or("Upper surface not found")?
            .dist_to_point(a);

        let lower_dist = self
            .lower
            .as_ref()
            .ok_or("Lower surface not found")?
            .dist_to_point(a);

        if upper_dist < lower_dist {
            Ok((*a, *b))
        } else {
            Ok((*b, *a))
        }
    }
}

/// Order the curves based on their direction from the test point.  The curve that is in the
/// direction of the upper surface is returned as the first element of the tuple, and the curve
/// that is in the direction of the lower surface is returned as the second element of the tuple.
///
/// # Arguments
///
/// * `a`:
/// * `b`:
/// * `test_point`:
///
/// returns: Result<(Option<Curve2>, Option<Curve2>), Box<dyn Error, Global>>
fn order_faces(
    a: Curve2,
    b: Curve2,
    test_point: SurfacePoint2,
) -> Result<(Option<Curve2>, Option<Curve2>)> {
    let a_t = a.intersection(&test_point);
    let b_t = b.intersection(&test_point);

    // We should have intersections with both curves. If the outline is clean, we will have exactly
    // one intersection with each, but if not we might have more than one at a similar distance.
    if a_t.is_empty() || b_t.is_empty() {
        Err("Failed to find intersections with the test point while ordering faces".into())
    } else {
        let a_m = compute_mean(&a_t)?;
        let b_m = compute_mean(&b_t)?;

        if a_m > b_m {
            Ok((Some(a), Some(b)))
        } else {
            Ok((Some(b), Some(a)))
        }
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
