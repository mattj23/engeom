//! This module and its submodules have tools for performing dimensional analysis on 2D airfoil
//! cross-sections.

pub mod camber;
pub mod edges;
mod geometry;
pub mod inscribed;
pub mod orient;

use crate::airfoil2::geometry::geometry_only_analysis;
use crate::airfoil2::inscribed::Inscribed;
use crate::{Curve2, Point2, Result};
pub use orient::{OrientFwdAft, OrientUpperLower};
use serde::{Deserialize, Serialize};

/// This struct is a general wrapper around an airfoil section input for common airfoil algorithms
/// implemented in this module and its submodules. It holds a reference to the airfoil section as
/// well as common tolerances that will be used by downstream algorithms.
///
/// Also, part of the reason for making this a separate struct is to handle
pub struct SectionInput<'a> {
    pub section: &'a Curve2,
    pub general_tol: f64,
    pub resolve_tol: f64,
}

impl<'a> SectionInput<'a> {
    pub fn new(section: &'a Curve2, general_tol: f64) -> Self {
        Self {
            section,
            general_tol,
            resolve_tol: general_tol * 0.1,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum AfEdgeSearch {
    Auto,
    Open,
    Sharp,
    Square,
    RoundedSquare,
    FullRound,
    BlendedRound,
}

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

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AfEdge {
    pub point: Point2,
    pub geometry: AfEdgeGeometry,
}

impl AfEdge {
    pub fn new(point: Point2, geometry: AfEdgeGeometry) -> Self {
        Self { point, geometry }
    }
}

#[derive(Clone)]
pub struct AfGeometry {
    pub leading: AfEdge,
    pub trailing: AfEdge,
    pub camber: Curve2,
    pub upper: Curve2,
    pub lower: Curve2,
    pub circles: Vec<Inscribed>,
}

impl AfGeometry {
    /// Conducts a purely geometric analysis of an airfoil section, attempting to extract the main
    /// camber line (MCL), identify the leading and trailing edge features and directions, and
    /// orient the upper and lower surfaces.
    ///
    /// # Arguments
    ///
    /// * `section`:
    /// * `general_tol`:
    /// * `fwd_aft`:
    /// * `upper_lower`:
    /// * `le_search`:
    /// * `te_search`:
    ///
    /// returns: Result<AfGeometry, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
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
