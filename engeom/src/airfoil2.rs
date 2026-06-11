//! This module and its submodules have tools for performing dimensional analysis on 2D airfoil
//! cross-sections.
//!
//! In the aerospace industry, common measurements taken on airfoils include:
//!
//! - Thickness measurements through the airfoil, especially around the leading/trailing edges and
//!   the position of maximum thickness, usually located by position along the mean camber line or
//!   in local directions at the edges.
//! - Chord lengths, often specified by a variety of different methods.
//! - Form and profile measurements, often with multiple tolerance regions, sometimes requiring
//!   special rules for partially constrained floating zones
//! - Section position and angle, measured from nominal references
//! - Leading and trailing edge position and shape, such as leading edge radius and trailing edge
//!   trim position, etc.
//!
//! It's important to start with the understanding that there's a huge difference between
//! running airfoil shape analysis tools on nominal section data, such as that exported from a CAD
//! system or a mathematical representation of design geometry, and running the same tools on
//! measured section data that came from a system like a 3D scanner or CMM.  Actual data will have
//! both noise from the measurement system and actual defects/roughness from the manufacturing
//! process. This brings a whole collection of problems that range from the practical to the
//! philosophical, all which need to be addressed at some level.
//!
//!

pub mod camber;
mod edges;
pub mod inscribed;
pub mod orient;

use crate::{Curve2, Point2};
use serde::{Deserialize, Serialize};

pub use orient::{OrientFwdAft, OrientUpperLower};

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

pub struct AfGeometry {
    pub leading: AfEdge,
    pub trailing: AfEdge,
    pub camber: Curve2,
}
