pub mod camber;
mod inscribed;

use crate::{Curve2, Point2};

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

pub struct AfGeometry {
    pub leading_edge: AfEdgeGeometry,
    pub trailing_edge: AfEdgeGeometry,

    pub camber: Curve2,
}
