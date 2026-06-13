//! This module and its submodules have tools for performing dimensional analysis on 2D airfoil
//! cross-sections.

/// Camber line extraction: finding the chain of inscribed circles between the leading and trailing
/// edges of an airfoil section.
pub mod camber;

/// Edge geometry fitting: routines that take an oriented inscribed circle stack and a section
/// curve and fit a specific edge geometry (sharp, square, rounded-square, full round, blended
/// round) at the leading or trailing edge.
pub mod edges;

mod geometry;

/// The `Inscribed` circle type and the `InscribedVec` collection used to manipulate an oriented
/// or partially-oriented stack of inscribed circles along the camber line.
pub mod inscribed;

/// Orientation strategies for resolving the forward/aft (leading/trailing edge) direction and
/// the upper/lower (suction/pressure) surfaces of an airfoil section.
pub mod orient;
mod position;

use crate::airfoil2::geometry::geometry_only_analysis;
use crate::airfoil2::inscribed::Inscribed;
use crate::airfoil2::position::{pos_camber, pos_offset, pos_radius};
use crate::{Curve2, CurveStation2, Point2, Result};
pub use orient::{OrientFwdAft, OrientUpperLower};
use serde::{Deserialize, Serialize};

/// Enum to specify between the upper and lower side of the airfoil
pub enum AfSide {
    /// The upper/suction/convex side of the airfoil
    Upper,

    /// The lower/pressure/concave side of the airfoil
    Lower,
}

/// Enum representing a method for specifying a location on an airfoil surface.
#[derive(Debug, Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum AfPos {
    /// A position specified by an arc distance along the mean camber line. Positive distances
    /// are measured from the leading edge point, negative ones are from the trailing edge point.
    OnCamber,

    /// A position specified by an intersection with a circle of a specified radius centered at
    /// the edge point. Positive radii indicate that the circle is centered at the leading edge
    /// point, and negative ones are centered at the trailing edge point.
    Radius,

    /// A position specified along a ray at the leading or trailing edge point in the direction of
    /// the camber-line tangency.  Positive values are measured from the leading edge point,
    /// negative ones from the trailing edge point. At the specified distance, an intersection
    /// is taken orthogonal to the tangent direction.
    EdgeOffset,
}

/// General wrapper around an airfoil section input for common airfoil algorithms implemented in
/// this module and its submodules. It bundles a reference to the section curve together with the
/// general tolerance used by downstream algorithms and the derived sub-tolerance used when an
/// algorithm needs to refine a result more tightly than the general tolerance.
pub struct SectionInput<'a> {
    /// The 2D airfoil cross-section curve being analyzed. May be open at one edge but not both.
    pub section: &'a Curve2,

    /// The general fitting tolerance used as the primary control parameter throughout the
    /// analysis (camber spacing, edge fit convergence, etc.).
    pub general_tol: f64,

    /// A tighter sub-tolerance used by algorithms that need to refine internal values more
    /// precisely than the user-facing general tolerance. Defaults to one tenth of `general_tol`.
    pub resolve_tol: f64,
}

impl<'a> SectionInput<'a> {
    /// Build a new `SectionInput` from a section reference and a general tolerance. The
    /// `resolve_tol` is set to one tenth of `general_tol`.
    pub fn new(section: &'a Curve2, general_tol: f64) -> Self {
        Self {
            section,
            general_tol,
            resolve_tol: general_tol * 0.1,
        }
    }
}

/// Selects which edge geometry to fit at the leading or trailing edge of an airfoil section
/// when running a geometric analysis.
#[derive(Debug, Clone, Copy)]
pub enum AfEdgeSearch {
    /// Try every fittable variant and return the one with the lowest average residual.
    Auto,

    /// Treat the edge as open and skip edge fitting. (Not yet implemented.)
    Open,

    /// Fit a sharp apex edge: a single point where the two surfaces meet.
    Sharp,

    /// Fit a flat-faced edge with two sharp corners.
    Square,

    /// Fit a flat-faced edge whose corners are blended by circular arcs of a single equal radius.
    RoundedSquare,

    /// Fit a full circular round joined to the surrounding surfaces by two short straight
    /// segments tangent to the round.
    FullRound,

    /// Fit a full circular round joined to the surrounding surfaces by two tangent blending arcs.
    BlendedRound,
}

/// The geometric description of a fitted airfoil edge, returned alongside a canonical edge
/// location in [`AfEdge`]. Each variant carries the parameters of the corresponding edge shape.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum AfEdgeGeometry {
    /// The section is known to be open at this edge
    Open,

    /// The edge comes to a sharp point
    Sharp(Point2),

    /// The edge consists of a square face with two sharp corners
    Square(Point2, Point2),

    /// The edge consists of a square face with two rounded corners that have the same radius
    RoundedSquare(Point2, Point2, f32),

    /// The edge is a full round joined to the rest of the airfoil by two short linear segments
    /// tangent to the round
    FullRound(Point2, f32),

    /// The edge is a full round joined to the rest of the airfoil by two tangent arcs
    BlendedRound(Point2, f32),
}

/// A fitted airfoil edge: a canonical edge location point together with a description of the
/// edge geometry. The meaning of `point` depends on the geometry variant; see [`AfEdgeGeometry`].
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AfEdge {
    /// The canonical edge location. For a sharp edge this is the apex; for a square or
    /// rounded-square edge it is the midpoint of the flat face; for the round variants it is the
    /// outermost point on the camber axis of the section.
    pub point: Point2,
    /// The geometric description of the edge shape.
    pub geometry: AfEdgeGeometry,
}

impl AfEdge {
    /// Build a new `AfEdge` from a canonical edge point and its associated geometry.
    pub fn new(point: Point2, geometry: AfEdgeGeometry) -> Self {
        Self { point, geometry }
    }
}

/// The result of a geometric analysis of an airfoil section.
///
/// Contains the fitted leading and trailing edges, the mean camber line, the segregated upper
/// (suction) and lower (pressure) surfaces, and the oriented stack of inscribed circles used
/// during the analysis.
#[derive(Clone)]
pub struct AfGeometry {
    /// The fitted leading edge.
    pub leading: AfEdge,

    /// The fitted trailing edge.
    pub trailing: AfEdge,

    /// The mean camber line, oriented so the first point is the leading edge point and the last
    /// point is the trailing edge point.
    pub camber: Curve2,

    /// The upper (suction) surface curve, split from the original section. The points are ordered
    /// so that the surface normals point outward from the airfoil skin.
    pub upper: Curve2,

    /// The lower (pressure) surface curve, split from the original section. The points are ordered
    /// so that the surface normals point outward from the airfoil skin.
    pub lower: Curve2,

    /// The inscribed circle stack used during analysis, ordered leading-to-trailing with each
    /// circle's `p0` on the lower surface and `p1` on the upper surface.
    pub circles: Vec<Inscribed>,
}

impl AfGeometry {
    /// Get the point on the specified side of the airfoil corresponding to the given position
    /// method and value.
    ///
    /// # Arguments
    ///
    /// * `side`:
    /// * `method`:
    /// * `value`:
    ///
    /// returns: Option<CurveStation2>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn af_point(&self, side: AfSide, method: AfPos, value: f64) -> Option<CurveStation2<'_>> {
        let target = match side {
            AfSide::Upper => &self.upper,
            AfSide::Lower => &self.lower,
        };

        match method {
            AfPos::OnCamber => pos_camber(value, &self.camber, target),
            AfPos::Radius => pos_radius(value, &self.camber, target),
            AfPos::EdgeOffset => pos_offset(value, &self.camber, target),
        }
    }

    /// Run a purely geometric analysis of an airfoil section, attempting to extract the mean
    /// camber line (MCL), identify the leading and trailing edge features, and orient the upper
    /// and lower surfaces.
    ///
    /// # Arguments
    ///
    /// * `section`: the airfoil cross-section curve, optionally open at one edge.
    /// * `general_tol`: general fitting tolerance used throughout the analysis (camber spacing,
    ///   edge fit convergence, etc.).
    /// * `fwd_aft`: strategy for deciding which end of the camber line is the leading edge.
    /// * `upper_lower`: strategy for deciding which surface is the upper (suction) side.
    /// * `le_search`: edge-fitting strategy to run at the leading edge.
    /// * `te_search`: edge-fitting strategy to run at the trailing edge.
    ///
    /// returns: `Result<AfGeometry>`. Errors if any stage of the analysis fails (e.g. too few
    /// inscribed circles to orient, edge fit failure, or both edges resolving as open).
    pub fn try_from_geometric_analysis(
        section: &Curve2,
        general_tol: f64,
        fwd_aft: OrientFwdAft,
        upper_lower: OrientUpperLower,
        le_search: AfEdgeSearch,
        te_search: AfEdgeSearch,
    ) -> Result<Self> {
        geometry_only_analysis(
            section,
            general_tol,
            fwd_aft,
            upper_lower,
            le_search,
            te_search,
        )
    }

    /// Returns the inscribed circle with the largest radius, which corresponds to the maximum
    /// thickness location along the camber line.
    pub fn tmax_circle(&self) -> &Inscribed {
        self.circles
            .iter()
            .max_by(|a, b| {
                a.c.r()
                    .partial_cmp(&b.c.r())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .expect("AfGeometry is guaranteed to contain at least one inscribed circle")
    }
}
