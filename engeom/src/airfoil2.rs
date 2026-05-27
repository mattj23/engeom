pub mod camber;

use crate::{Curve2, Point2};

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
