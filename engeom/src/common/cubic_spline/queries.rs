//! This module is to group spatial query related actions on the cubic Bezier spline struct.
//!
//! The main objective of this module is to get an accurate, reliable, and reasonably fast
//! implementation of the closest point query algorithm.  The common approaches of lookup table
//! pre-search can produce incorrect results when a query is done near spots where the curve
//! loops around or crosses itself, and you end up pruning the interval with the true minimum value
//! because an endpoint happened to be closer.
//!
//! Besides the incorrect pruning problem, the other danger is not handling searches correctly in
//! areas of the curve with extreme changes that aren't represented well by the discrete points at
//! interval boundaries.
//!
//! Initial thoughts:
//! - I want to avoid having to supply a tolerance for an acceleration structure
//! - Generalizing over two and three dimensions may be making this more complicated
//! - If there's some way I could break the spline up into simpler segments with some known maximum
//!   complexity, that would probably let me accomplish what I want
//! - It would be nice to have a single acceleration structure that could also be used as a
//!   foundation for intersection queries
//! - Bounding boxes feel heavy, but may be necessary
//!
//! Some stuff I think I know about cubic Béziers:
//! - Derivative roots for the individual components can be found via quadratic formula
//! - In two or three dimensions, it's only possible to have one cusp, and that cusp can be found
//!   (I think?) where the derivatives of all dimensions go to zero at the same place.
//! - In 2D, splitting at a cusp leaves you with two extremely simple arc-like curves
//! - What if we split so that we get segments that are convex and non-intersecting?
//!     - If there's a cusp, split there
//!     - There should be at most two inflection points on both 2D or 3D(?) splines
//!     - Curves that are loops can have no inflection points, but they aren't convex
//! - The goal should be to split into segments where extreme points are easily found in arbitrary
//!   directions?
//!
//! What types of queries and so what types of extremes?
//! - Closest point query: minimal distance to coordinates
//! - Closest point to halfspace: minimum distance to halfspace
//! - Halfspace intersection: minimum and maximum distance to halfspace
//! - Closest point to circle/sphere: minimum and maximum distance to coordinates
//! - Circle/sphere intersection: minimum and maximum distance to coordinates
//! - Intersections with other spline?
//! - Closest distance to other spline?
//!
//! What condition would create more than one local minima? How would you find both?
//! - Under any(?) circumstances, a local minima _or_ maxima will be at one of these:
//!     - A point where the derivative is orthogonal to the query direction
//!     - One of the interval endpoints
