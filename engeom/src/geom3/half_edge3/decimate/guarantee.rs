//! This module has the tolerance-bounded mesh decimation implementation based the 1997 IBM research
//! paper **Guéziec, _Surface Simplification Inside a Tolerance Volume_, IBM RC 20440, 1997.**
//!
//! # Overview
//!
//! Guéziec's paper has two main parts of interest to us.
//!
//! 1. It has a scheme for representing a volume around a mesh by per-vertex radii and a
//!    corresponding proof that using those volumes to track error during the mesh collapse process
//!    yields a check that can be used to guarantee a symmetric Hausdorff distance limit is never
//!    violated.
//!
//! 2. It offers several different methods to update the volumes after a collapse to ensure the
//!    conditions required to get the guarantee in #1 are satisfied.
//!
//! The scheme and the proof are well reasoned and presented, and also extremely clever. The methods
//! were well thought out but of the two most complete versions, I was unable to get his preferred
//! one (the subdivision method) to meet the guarantee.  The other (the projection method) works as
//! expected.  This module is based on the projection method.
//!
//! # How the bound works: two volumes
//!
//! Imagine a triangle with a sphere of an arbitrary radius at each vertex.  Now take the convex
//! hull of those three spheres[^1] and look at its enclosing volume. This volume can be used to
//! represent a zone of varying thickness (within the resolution of the vertices) over the surface
//! of a mesh.
//!
//! The volume can also be thought of as the union of all balls (spheres with the centerpoint on the
//! surface) for every one of the infinite number of points lying on the face of the triangle.  
//! It's worth mentioning that along with the volume itself, the concept of the surface's "balls"
//! are important...there is one ball at every point in every triangle, and its radius can be found
//! by the barycentric interpolation of the radii at the vertices.
//!
//! Guéziec's method requires two such simultaneous volumes:
//!
//! **The tolerance volume** is the volume around the original mesh which the simplified mesh must
//! not exceed. Guéziec's method can work on per-vertex tolerances, which requires this volume to be
//! stored.  However, for a spatially-uniform tolerance, this volume never actually needs to be
//! calculated.
//!
//! **The error volume** is a volume over the *decimated* mesh which can be thought of as an upper
//! bound of simplification error to the original surface accumulated over the course of the
//! algorithm.  It starts at zero size on every vertex before the first edge collapse, and is
//! updated after every collapse in a way that meets the following two criteria:
//!
//! - (𝔸): The updated volume completely contains the pre-collapse volume.
//! - (𝔹): Every ball of the error volume contains at least one ball of the previous error
//!   volume.
//!
//! If the error volume started at zero size and coincident with the original mesh, and if updates
//! are always performed in a way that satisfies (𝔸) and (𝔹) then it is provable by induction that
//! at any step:
//!
//! - The error volume completely contains the original surface.
//! - Every ball of the error volume contains at least one point (because the original radii were
//!   zeros) of the original surface.
//!
//! That means that, at any step, the error volume represents the upper bound of the symmetrical
//! Hausdorff distance to the original mesh.  As long as the error volume never breaks out of the
//! tolerance volume, you can claim with certainty that the symmetrical Hausdorff distance between
//! the original and the decimated mesh is below the tolerance value.
//!
//! Notice what this does _not_ do: it doesn't tell you anything about the actual Hausdorff
//! distance besides the largest number it could possibly be. You will never actually know through
//! this method _when_ the Hausdorff threshold is violated, only when it _definitely isn't_. The
//! method's compactness and relative simplicity come from only choosing to keep track of the
//! upper bound.
//!
//! [^1]: Guéziec describes this as a ball swept along the mesh with its radius interpolated
//!       barycentrically, yielding a conical tube on the edges, spheres at the vertices, and
//!       triangluar prisms on the sides. This is, more simply, the convex hull of three spheres.
//!
//! # Theory behind the error bound calculation
//!
//! Like other mesh decimation methods, Guéziec works with edge collapse. He calculates what the
//! updated error volume would be for the post-collapsed mesh state, then looks to see if the error
//! radius at any of the affected vertices would be larger than tolerance, and if so he vetoes the
//! collapse and moves on.
//!
//! It's the problem of _how_ to do the updated error volume that ends up being the obstacle. I'm
//! going to explain Guéziec's method first, and then explain why it works, which is what took me
//! the longest time to get my head around.
//!
//! <div class="warning">
//!
//! This will be hard to understand if you don't have some idea of the notation and vocabulary. In
//! the paper, Guéziec refers to the edge being collapsed as having endpoints v₁ and v₂ and
//! collapsing to vertex v₀.
//!
//! You will need to understand the concept of a _star_ in a simplicial complex, which in this
//! case is the ring of triangles attached to a structure. The star of the vertex v₀, which he
//! denotes as (v₀)★, is the set of triangles which have v₀ as a vertex. The star of edge v₁v₂
//! he writes as (v₁,v₂)★, and it's the union of the stars of v₁ and v₂.
//!
//! Lastly, the _link_ of a vertex or edge is the closed polygonal curve created by the edges of
//! its star that don't touch the entity itself. He writes the link of v₀ as ℓ(v₀), and the link
//! of edge v₁v₂ as ℓ(v₁,v₂).  See Figure 4 in the paper for a visual.
//!
//! The code uses Guéziec's numbering, so `v1` and `v2` are the endpoints of the collapsed edge and
//! anything named `0` belongs to the merged vertex.
//!
//! | paper | code | what it is |
//! |---|---|---|
//! | v₁ | `v1` | `h.tail()`, the endpoint the collapse **deletes** |
//! | v₂ | `v2` | `h.head()`, the endpoint that **survives** |
//! | v₀ | `p0` | where the merged vertex goes, as a position |
//! | p(v₁), p(v₂) | `p1`, `p2` | where the two endpoints are now |
//!
//! </div>
//!
//! The first thing to note is that during an edge collapse, when edge v₁v₂ gets collapsed to vertex
//! v₀, the only affected triangles are the ones in the star (v₁,v₂)★, and all that happens to them
//! is that they get hot-swapped for (v₀)★, which contains two less triangles, three less edges, and
//! one less vertex.  The link of the two stars is the same, in that none of those vertices are
//! getting deleted or moved by the operation.
//!
//! That means that to satisfy (𝔸) and (𝔹) we only have to satisfy them for the volume attached to
//! the vertex star (v₀)★ by absorbing the volume attached to the pre-collapse edge star (v₁,v₂)★.
//!
//! To do this Guéziec creates two mappings (or one bijective map that does the job of both), one
//! that takes any point on (v₁,v₂)★ and maps it to a point on (v₀)★, and the other which takes any
//! point on (v₀)★ and maps it to a point on (v₁,v₂)★.  Points don't have to round-trip on the
//! two-map version, though they will by definition on the bijection.  The mapping will be used to
//! assemble a set of inequality constraints, which can then be "solved" through any method of
//! choice.
//!
//! "But _how_ is this mapping supposed to decide which points get mapped to which?", you might
//! think to ask. "Surely that matters, and the paper never mentions it!"
//!
//! Actually, it doesn't matter as long as it's affine...at least it doesn't matter for correctness.
//! You can map the points however you want and the tolerance guarantee will still hold.  The method
//! will not perform many edge collapses, but you will not violate the tolerance.
//!
//! The purpose of the two maps (or the two mapping directions) is to calculate the error radii
//! needed to satisfy the two conditions (𝔸) and (𝔹). The mapping from (v₁,v₂)★ to (v₀)★ is used to
//! enforce the nesting criteria (𝔸), and the reverse mapping from (v₀)★ to (v₁,v₂)★ is used to
//! enforce the ball containment criteria (𝔹).
//!
//! How? Let's consider them individually.
//!
//! For a sphere to be enclosed in another sphere, the radius of the larger sphere must be _at
//! least_ the radius of the smaller sphere plus the distance between their centers. For a ball on
//! the new surface (v₀)★ located at center cₙ with error radius εₙ to enclose a ball on the old
//! surface (v₁,v₂)★ located at center cₒ with error radius εₒ:
//!
//! > εₙ ≥ |cₙ - cₒ| + εₒ
//!
//! Condition (𝔸) says that the entire pre-collapse volume must be contained within the
//! post-collapse volume.  That's equivalent to saying that _every_ ball in (v₁,v₂)★ must be
//! contained within at least one ball in (v₀)★.  It doesn't matter which one, we just need to know
//! that _at least_ one ball is big enough and situated somewhere that it fully encloses each point.
//!
//! Condition (𝔹) says that _every_ ball in (v₀)★ must contain at least one ball in (v₁,v₂)★.
//! That's basically the same relationship, except that if the mapping _isn't_ a bijection, you
//! could have some points in (v₀)★ that the mapping from (v₁,v₂)★ misses.
//!
//! We can look up the error radius εₒ of any point on (v₁,v₂)★ by using the barycentric
//! coordinates. We can rearrange the terms on so that we can express the error raidus εₙ for
//! any point on (v₀)★ as a linear combination of the error radius of the vertices and its own
//! barycentric coordinates.  We can build the inequality where we specify what vertex error radii
//! have to be to satisfy the inequality constraint...if only we just had center coordinates to plug
//! into the terms.
//!
//! Remember that the error volumes are supposed to track the _upper bound_ of the Hausdorff
//! distance, not to compute it. The distance term |cₙ - cₒ| will be its smallest at the closest
//! point projection from one surface to the other, but the inequality is valid at any pair of
//! points...the error radius εₙ just grows larger than necessary to absorb the slack.  So, while
//! your error volume will grow the _slowest_ the better your mapping approximates the closest point
//! pairs from (v₁,v₂)★ to (v₀)★, any mapping at all will be safe.  It just may not be useful.
//!
//! The only constraint on the maps is that they have to be made of pieces that are locally
//! _affine_, and the only way to make affine mappings between two sets of triangles with different
//! topologies is to do it in pieces that lie on a single triangle in the domain and a single
//! triangle in the co-domain.
//!
//! Why this strange requirement?  Guéziec doesn't explain _why_ the maps have to be
//! piecewise-affine (he says piecewise-linear), instead he just describes how the constraints are
//! collected and mentions that they only have to be gathered from the corners of the pieces.  He
//! gives no explaination of why.
//!
//! It turns out there's an unwritten convexity shortcut he's using. In the inequality constraint,
//! the distance term |cₙ - cₒ| will be convex if the mapping between the two surfaces is
//! affine[^2].  The error radius εₒ is a barycentric interpolation, so it's also affine.  A convex
//! function plus an affine function makes another convex function, so the error radius inequality
//! is convex over each piece of the mapping.  Convexity means that...while the _minimum_ could be
//! literally anywhere in the piece, the _maximum_ is absolutely, positively going to be at one of
//! the corners.  The fact that we're computing the upper bound of the Hausdorff distance means that
//! the maximum was all we ever really cared about.  So, as long as we set up our constraints at the
//! corners, we're guaranteed to capture the worst case.
//!
//! [^2]: To think through why, think about two points on two surfaces under any affine mapping
//!     between them.  The affine property of any arbitrary mapping between the surfaces means that
//!     moving the point in a straight line at a constant speed in the domain results in the
//!     associated point in the co-domain also moving on some arbitrary straight line at a constant
//!     speed.  Regardless of what those lines are and what the speeds might be, there's always
//!     going to be some closest approach between the two points, forming a distance minimum. On
//!     either side of the closest approach the distance between them grows out to inifinty.  That's
//!     a textbook example of a convex function.
//!
//! # Error computation methods
//!
//! Alongside the theoretical method, Guéziec publishes two practical methods for building the
//! arbitrary piecewise-affine maps.
//!
//! The first, which he calls the "subdivision method", constructs a bijective map by first
//! projecting v₀ onto (v₁,v₂)★, then v₁ and v₂ onto (v₀)★, and then building a bunch of nodes and
//! edges linking vertices until he has two topologically identical graphs.  From there he can
//! plug the nodes into both, and since the map is bijective (𝔸) and (𝔹) yield the same constraint
//! inequality and he only has to evaluate it once.
//!
//! The other method he calls the "projection method", and it involves projecting the triangles of
//! each star onto the other, and then identifying the points where edges, corners, and vertices end
//! up. This gives him two mappings, where the correspondence is geometric rather than topological,
//! and while it's twice as much work to assemble the constraints it, in theory, could be a little
//! more represenatative of the Hausdorff distance.
//!
//! ## The subdivision method
//!
//! The best results that the paper reported were with this method, which involves creating a
//! piecewise-linear homeomorphism between the original edge star (v₁,v₂)★ and the edge-collapsed
//! vertex star (v₀)★.
//!
//! This works by starting with two separate graphs, one he calls _G_ which lives in (v₁,v₂)★, and
//! one called _H_ which lives in (v₀)★. They both contain the exact same link vertices, but _G_ has
//! v₁ and v₂ and the associated edges and _H_ has v₀ and its associated edges. The goal is to make
//! them topologically identical while also satisfying the need to produce the individual affine
//! pieces.  To start he projects v₀ into (v₁,v₂)★ and inserts it as a node into _G_, then projects
//! v₁ and v₂ into (v₀)★ and inserts them as nodes into _H_. He then iterates through all the
//! edge pairs, building nodes at the closest approaches between them, and sampling constraints at
//! the corners.
//!
//! At least, that's my understanding.  Despite two concerted attempts, this method never worked for
//! me. It held the guarantee at the mesh vertices, but the dense sampling tests showed it failed
//! elsewhere.
//!
//! ## The projection method
//!
//! I won't say much about Guéziec's projection method, beyond that it attempts to create the
//! mapping pieces through geometric projection and ends up creating two separate maps rather than
//! a bijection. Closest distance across two surfaces is not a bijection except under special
//! conditions, so in theory the dual mapping could produce more effective results, although the
//! lack of a bijection means that constraints have to be evaluated twice.
//!
//! I read through the method's implementation notes in the paper's appendix, and could not build
//! a strong spatial intuition for how it worked.  I set it aside to see if there was another
//! option to try first.
//!
//! ## The projected overlay
//!
//! This was an attempt at taking a shortcut to building something closer to the substitution
//! method's bijection, with the caveat that it doesn't work under all cases.
//!
//! The idea is that, if you consider that both stars are disks with the same boundary, there's
//! a space where they line up exactly: projected to a plane.  If you can find a plane where, in
//! a projection, none of the triangles are flipped or squashed to a line, you will be starting
//! with most of the work done already.  The two planar polygons have the exact same boundary and
//! the exact same internal area; the only thing different about them is where the interior
//! vertices are.
//!
//! By taking the projected triangles from (v₁,v₂)★ and the projected triangles from (v₀)★ and
//! finding all of their edge crossings, plus the projected locations of v₀, v₁, and v₂, you find
//! the corners of reasonable affine maps between the two stars, and thus the positions to pull
//! constraints from.
//!
//! It's conceptually simple, though it does require a lot of naive segment intersection tests to
//! find the edge crossings without any sort of acceleration structure.  That said, the subdivision
//! method has that same issue.
//!
//! In the end, this does something very similar to the substitution method with a little less
//! computation, but the price paid for that is incompleteness...any star that folds doesn't get
//! considered. That tends to not be a lot of triangles, but it's some finite percentage in every
//! practical mesh.
//!
//! When that happens, we use the affine face map, mentioned below. The `fallback_rate` harness
//! measures how often that is; see that harness for the table. On a locked-boundary run the overlay
//! judges every collapse the run actually performs on an open mesh just as well as on a closed one,
//! and what is left over is folding: under half a percent on the blade even at a tolerance five
//! times the edge length, and around five percent on the bunny, which is curved in every direction.
//!
//! Also, a collapse which moves the star's outline does make the two stars span different regions,
//! and the bijection stops existing there. That case is by the specialized boundary mechanics,
//! refer to the boundary section further down in this document.
//!
//! ## The affine face map
//!
//! This was actually a test component made to be plugged in just long enough to make sure that the
//! rest of the machinery worked, and then thrown away.  In the end, it didn't get to retire
//! because it happens to "work" in places where the projected overlay fails.  If the projected
//! overlay ends up getting replaced with something more complete in the future, the affine face
//! map can probably go with it.
//!
//! The affine face map does _not_ do a good job of approximating the closest point mappings between
//! the stars.  It's sound only because _any_ full affine mapping between the stars is sound. But
//! it's super cheap to assemble, and it doesn't care about the actual geometry (hence why it's
//! not great at its job).
//!
//! It works by looking for a replacement for each triangle in the other star and then just
//! mapping it whole, which ends up meaning only that it takes the constraint at any non-link
//! vertices the triangle has.
//!
//! - For faces in (v₀)★, it looks for the face in (v₁,v₂)★ that has the same two link vertices.
//! - For faces in (v₁,v₂)★ that have two link vertices, it looks for the face in (v₀)★ that has
//!   the same two link vertices. For faces that have only one link vertex (because the other two
//!   vertices are v₁ and v₂) it looks for every triangle in (v₀)★ that has that link vertex, and
//!   it uses the smallest constraint of them.
//!
//! # Algorithm mechanics
//!
//! ## Structural (not error based) vetoes
//!
//! There are also two structural conditions which will veto a collapse.
//!
//! - A collapse which turns a face through more than `max_normal_deviation` is refused, catching
//!   flips and folds where a triangle passes through the surface and lands within tolerance on the
//!   far side.
//!
//! - A collapse which drives a face below `min_aspect` is refused, since a sliver is numerically
//!   useless even sitting exactly on the surface.
//!
//! ## Repeated runs
//!
//! The radii live on the mesh as `alum` vertex properties rather than in the decimator, so they
//! survive the run that produced them. A second decimation starts from the error the first one
//! accumulated and checks it against a budget that has not moved, so three runs at a given
//! tolerance land inside that tolerance and not inside three of them.
//!
//! ## Ordering and staleness
//!
//! Collapses are queued per vertex by the quadric residual at the position that was accepted, and
//! the queue is refreshed for the one-ring after each collapse. Ordering and deciding are separate
//! questions, so the quadric acts as a ranking, and the gate(s) allow or disallow an actual
//! collapse.
//!
//! A queued cost can go stale because a collapse two rings away can move geometry or grow a radius
//! it depended on, so the full test is run again when the entry is popped and the collapse is
//! skipped if it no longer passes. `alum`'s own decimation driver cannot express this, since its
//! `Decimater::before_collapse` hook has no way to refuse, which is why the collapse loop here is
//! driven directly against `alum`'s primitives.
//!
//! # Boundary simplification
//!
//! Guéziec's Section 10 sorts boundary edges into three kinds and treats them separately.
//!
//! - A **Type I** edge is on the outline itself, carrying one triangle and two boundary vertices.
//! - A **Type II** edge has exactly one boundary vertex and two triangles.
//! - A **Type III** edge has two boundary vertices but two triangles, joining two different parts
//!   of the outline, and he refuses to collapse it at all because doing so merges two boundaries
//!   and leaves a singular vertex.
//!
//! Type II he handles in a paragraph: slide the interior vertex along the edge until it sits on
//! the boundary vertex, which becomes the simplified vertex, and the rest of the tests are as
//! before. That is what `lock_boundary` already does here, since the boundary vertex is pinned and
//! the interior one collapses into it, and the overlay handles it unaided because the outline does
//! not move. See the `fallback_rate` harness for the measurement.
//!
//! Type I is the one that is genuinely undone. Guéziec gives a positioning rule; placing the
//! simplified vertex on an ellipsoid whose focal points are the two outer neighbors so that the
//! boundary's length is preserved, which is the outline's analogue of his volume preservation
//! that I didn't implement. For the part that matters here, the error bound, he writes only that
//! the test "must be adapted to the particular topology of Types I and II boundary edge
//! neighborhoods", that "it would be tedious to specify completely these tests here", and pulls a
//! classic "_this is left as an exercise for the reader."_
//!
//! So, to that end, everything below got worked out separately from Guéziec's paper.
//!
//! ## The Type I correspondence
//!
//! The two stars' outlines differ in one place only. Both are disks; both are bounded partly by the
//! link, which the collapse neither moves nor deletes, and partly by a **chain** of mesh boundary
//! edges running through the endpoints. The collapse turns that chain, `(a, v₁, v₂, b)`, into
//! `(a, v₀, b)`. Everything else about the two regions is identical.
//!
//! What makes the difference tractable is that the two invariants are pointwise existentials.
//! (𝔸) asks that each old ball lie inside _some_ new ball, and (𝔹) that each new ball contain
//! _some_ old one. Neither asks for a single global map, which we wouldn't be able to easily give.
//! So they need a _cover_ by affine pieces rather than a bijection, and the pieces are free to
//! come from different constructions:
//!
//! - over the region the two stars share, the overlay's cells cover both, and every constraint it
//!   already collects stays valid untouched. A crossing lies on an old edge and a new edge at once,
//!   so it is always in the shared region;
//!
//! - what is left is a **sliver** of each star with no counterpart on the other, and those get
//!   pieces of their own.
//!
//! A sliver is never constructed as a polygon. The slack under any affine map is concave, and a
//! concave function on a polygon is minimized at a vertex _even when the polygon is not convex_,
//! since the polygon sits inside its own hull whose extreme points are all polygon vertices. So
//! only the corner set is needed, and it is three cheap enumerations: face corners the other star
//! does not cover, where the face's edges meet the chain, and where the chain bends inside the
//! face.
//!
//! The map a sliver is carried by is the _nearest point of the chain_, which is piecewise affine
//! with pieces given by the chain's own Voronoi regions: a slab either side of each segment, where
//! a perpendicular projection is admissible by construction, and a wedge at each vertex, where that
//! vertex is the nearest chain point and a constant map charges exactly the right distance.
//!
//! Finally, on a moving outline the overlay and the face-to-face map stop being ranked. Each
//! returns a radius which is valid on its own, so the smaller is taken, which picks a method rather
//! than blending two. The overlay wins by a mile on a straight run, where the slivers are slivers
//! indeed; the face map still wins on a badly bent one.
//!
//! ```text
//! boundary unlocked, faces remaining      declined   Type I correspondence
//! flat-grid-20, tol=1e-6                    80            7
//! flat-grid-60, tol=1e-6                   266            6
//! stanford-bunny, tol=0.005                346          345
//! stanford-bunny, tol=0.02                 147          142
//! stanford-bunny, tol=0.1                   59           59
//! ```
//!
//! The flat grids are the case this was for: a straight outline can be collapsed away for nothing,
//! and now is. Both are checked against the invariants directly, one collapse at a time,
//! by the `invariants` submodule.

/// Which part of a star's outline a collapse can move, which decides whether the two stars cover
/// the same region and so whether a correspondence between them exists at all.
mod boundary;

/// What an error method is told about a collapse, and the shape it has to implement.
mod collapse;

/// One constraint on the merged vertex's radius, and the fold that turns a set of them into one.
mod constraint;

/// Gueziec's Contraction Method, the loosest of the three and the reference the others are
/// checked against.
mod contraction;

/// Gueziec's face-to-face affine map, which is sound on a moving outline where the overlay is not.
mod face_map;

/// A direct check of the two invariants, on one hand-built collapse at a time, together with the
/// fixtures it runs on. Test-only, and the thing which decides whether a correspondence is right.
#[cfg(test)]
mod invariants;

/// The projected overlay method and the machinery only it needs.
mod overlay;

use super::{
    AlumPolyMesh, CollapseGate, DecimateReport, DecimateStats, DecimateTarget, Placement, Quadric,
    QuadricKind, StarFace, Structural, candidate_positions, quadric_of_vertex, run_decimation,
    set_boundary_locks,
};
use crate::geom3::half_edge3::HalfEdgeMesh3;
use crate::{Result, Vector3};
use alum::{HH, Handle, HasIterators, HasTopology, VH};
use boundary::boundary_chain;
use collapse::{Collapse, ErrorRule};
use constraint::ErrorBound;
use contraction::Contraction;
use face_map::AffineFaceMap;
use overlay::ProjectedOverlay;
use std::f64::consts::PI;

/// How the error volume's growth is worked out on each collapse.
///
/// Both keep the surfaces inside the tolerance. They differ in how much slack they leave, and the
/// slack is what decides how small a tolerance is usable at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ErrorMethod {
    /// Map both endpoints onto the merged vertex and widen its ball until it swallows both.
    ///
    /// One line, and valid without any argument beyond that sentence, which makes it the reference
    /// the others are checked against. As a working method it is very loose: the radius grows by
    /// the distance the vertex travelled on every collapse, so an edge shorter than about twice the
    /// tolerance can never be collapsed at all.
    Contraction,

    /// Map each face of the star affinely onto the face which replaces it, and solve for the
    /// smallest radius satisfying the resulting constraints.
    ///
    /// Sound, and tighter per face than contracting, because an endpoint only has to reach the face
    /// replacing its own rather than travel all the way to the merged vertex. On a flat region that
    /// distance is zero.
    ///
    /// It nonetheless decimates no better, because the other invariant reintroduces the same
    /// barrier: the merged vertex must keep reaching back to *every* face of the old star, and the
    /// far side of a star is an edge length away. One affine map per face pair cannot avoid that.
    /// Escaping it needs the correspondence subdivided so that a face maps onto the part of the old
    /// star nearest it, which is Gueziec's Steps IV and V. See the module documentation.
    ///
    /// Its other job is as the overlay's fallback, and the reason it can hold a bound the overlay
    /// once could not is that it never assumed the two stars cover the same region: its pieces are
    /// one triangle onto one triangle, so a collapse which moves the mesh outline costs it nothing
    /// in validity. It still earns its keep there on a badly bent outline, where a sliver of the
    /// overlay reaches further from the chain than a face does from its own replacement.
    AffineFaceMap,

    /// Overlay the two stars in a shared projection and constrain at every vertex of the overlay.
    ///
    /// The only one of the three which decimates a flat region, because it is the only one whose
    /// map between the stars is a bijection rather than a face-to-face correspondence.
    ///
    /// Falls back to [`ErrorMethod::AffineFaceMap`] on a star too curved for either projection to
    /// stay injective. Where the mesh has a boundary and it is not locked, a collapse can move the
    /// star's outline so that the bijection stops existing; that case is covered rather than
    /// declined, and the two methods are then run against each other and the smaller radius kept.
    /// See the module documentation.
    #[default]
    ProjectedOverlay,
}

/// How decimation should behave.
///
/// Marked `#[non_exhaustive]`, so build one with [`DecimateOpts::to_tolerance`] or
/// [`DecimateOpts::to_face_count`] and the chained setters.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DecimateOpts {
    /// Farthest the surface may move from the original, in the mesh's own units.
    ///
    /// A genuine bound in both directions and over the surfaces, not at sampled points: no part of
    /// the result is further than this from the original, and no part of the original is further
    /// than this from the result. It holds across a whole run and across repeated runs rather than
    /// per collapse.
    ///
    /// The bound is conservative, since it is enforced through a union of balls swept over the
    /// mesh and that is a coarse shape to describe a thin shell with, but not by much: the
    /// deviation actually achieved lands at 60% to 100% of what was asked for. See the module
    /// documentation.
    pub tol: f64,

    /// The largest angle, in radians, a face may be turned through by a collapse.
    ///
    /// Catches flips and folds, where a triangle passes through the surface and lands within
    /// tolerance on the far side. The default of `PI / 4` is permissive enough not to obstruct
    /// ordinary decimation of curved regions.
    pub max_normal_deviation: f64,

    /// The lowest triangle quality a collapse may produce, on a scale where an equilateral triangle
    /// is 1 and a degenerate one is 0.
    ///
    /// A sliver is numerically useless even when it lies exactly on the surface, so this is a
    /// separate concern from deviation. The default of `0.01` rejects only genuinely bad triangles.
    pub min_aspect: f64,

    /// Whether vertices on a boundary are held fixed.
    ///
    /// An open mesh loses its outline if the boundary is free to collapse. Locking costs nothing
    /// per collapse, since it is applied once as a status flag before the run.
    ///
    /// Each run decides this afresh: a run with it `false` frees a boundary an earlier run on the
    /// same mesh locked.
    pub lock_boundary: bool,

    /// A tighter tolerance applied to vertices on a boundary, when they are not locked outright.
    ///
    /// `None` uses `tol` everywhere. This is simply a different starting budget on those vertices,
    /// applied once when the tolerance volume is seeded rather than branched on per collapse.
    pub boundary_tol: Option<f64>,

    /// Where the surviving vertex of a collapse is placed.
    pub placement: Placement,

    /// When to stop.
    pub target: DecimateTarget,

    /// Which quadric formulation to accumulate.
    pub quadric: QuadricKind,

    /// How the error volume's growth is worked out on each collapse.
    pub error_method: ErrorMethod,
}

impl Default for DecimateOpts {
    fn default() -> Self {
        Self {
            tol: 0.0,
            max_normal_deviation: PI / 4.0,
            min_aspect: 0.01,
            lock_boundary: true,
            boundary_tol: None,
            placement: Placement::Optimal,
            target: DecimateTarget::ToTolerance,
            quadric: QuadricKind::Triangle,
            error_method: ErrorMethod::ProjectedOverlay,
        }
    }
}

impl DecimateOpts {
    /// Decimate as far as the given deviation tolerance allows.
    pub fn to_tolerance(tol: f64) -> Self {
        Self {
            tol,
            target: DecimateTarget::ToTolerance,
            ..Default::default()
        }
    }

    /// Decimate towards a face count, still bounded by the given deviation tolerance.
    pub fn to_face_count(tol: f64, faces: usize) -> Self {
        Self {
            tol,
            target: DecimateTarget::FaceCount(faces),
            ..Default::default()
        }
    }

    /// Decimate towards a fraction of the starting face count, still bounded by the tolerance.
    pub fn to_ratio(tol: f64, ratio: f64) -> Self {
        Self {
            tol,
            target: DecimateTarget::Ratio(ratio),
            ..Default::default()
        }
    }

    /// Set the maximum angle a face may be turned through, returning the modified options.
    pub fn with_max_normal_deviation(mut self, value: f64) -> Self {
        self.max_normal_deviation = value;
        self
    }

    /// Set the minimum triangle quality, returning the modified options.
    pub fn with_min_aspect(mut self, value: f64) -> Self {
        self.min_aspect = value;
        self
    }

    /// Set whether boundary vertices are locked, returning the modified options.
    pub fn with_lock_boundary(mut self, value: bool) -> Self {
        self.lock_boundary = value;
        self
    }

    /// Set a tighter tolerance for boundary vertices, returning the modified options.
    pub fn with_boundary_tol(mut self, value: f64) -> Self {
        self.boundary_tol = Some(value);
        self
    }

    /// Set where the surviving vertex is placed, returning the modified options.
    pub fn with_placement(mut self, value: Placement) -> Self {
        self.placement = value;
        self
    }

    /// Set which quadric formulation to use, returning the modified options.
    pub fn with_quadric(mut self, value: QuadricKind) -> Self {
        self.quadric = value;
        self
    }

    /// Set how the error volume grows on each collapse, returning the modified options.
    pub fn with_error_method(mut self, value: ErrorMethod) -> Self {
        self.error_method = value;
        self
    }

    fn validate(&self) -> Result<()> {
        if self.tol.is_nan() || self.tol <= 0.0 {
            return Err(
                format!("A decimation tolerance must be positive, got {}", self.tol).into(),
            );
        }

        if self.max_normal_deviation.is_nan() || self.max_normal_deviation <= 0.0 {
            return Err("The maximum normal deviation must be positive".into());
        }

        if let Some(t) = self.boundary_tol
            && (t.is_nan() || t <= 0.0)
        {
            return Err("A boundary tolerance must be positive".into());
        }

        self.target.validate()
    }
}

/// What a collapse would do, once it has been found acceptable.
pub(crate) struct Accepted {
    /// Where the merged vertex goes.
    position: Vector3,

    /// The quadric the merged vertex inherits, used for ordering the next round.
    quadric: Quadric,

    /// The quadric residual at `position`, which is the queue key.
    cost: f64,

    /// The error radius the merged vertex takes on.
    error: f64,

    /// The tolerance the merged vertex takes on, which is the minimum over the edge star.
    tolerance: f64,
}

/// A decimator which orders collapses by quadric error but decides them by whether the collapse
/// keeps the surfaces inside Gueziec's error and tolerance volumes.
///
/// See the module documentation for why ordering and deciding are two different questions, and for
/// what the two volumes each bound.
pub struct ToleranceVolumeDecimator {
    opts: DecimateOpts,
    quadrics: Vec<Quadric>,

    /// The error volume, as one radius per vertex.
    ///
    /// The union of balls swept over the mesh with these radii, interpolated barycentrically
    /// between vertices, contains the original surface, and every ball intersects it. So the radius
    /// at a vertex bounds how far that part of the surface has moved, in both directions at once,
    /// and it does so over the surface rather than at a handful of sampled sites.
    ///
    /// Radii start at zero, since a mesh is its own reference until something is collapsed, and
    /// grow monotonically. That growth is what carries the bound across a whole run without ever
    /// consulting the original mesh again: the volume after a collapse contains the volume before
    /// it, so by induction it contains the original surface.
    ///
    /// Held as an `alum` vertex property rather than a side vector so that it survives collapses
    /// and garbage collection, and so a second decimation run continues from the error already
    /// accumulated rather than resetting and measuring against its own previous output.
    error: alum::VProperty<f64>,

    /// The tolerance volume, as one budget per vertex.
    ///
    /// Gueziec's implemented form needs no geometry: a collapse takes the minimum tolerance over
    /// its edge star, and is refused if any error in the resulting vertex star would exceed it.
    /// Because the minimum propagates, a point of the original surface is always mapped onto a
    /// triangle at least one of whose vertices carries a tolerance no larger than its own, which is
    /// what closes the argument.
    ///
    /// Per-vertex rather than a single scalar because that costs nothing here and buys a tighter
    /// budget on features that matter, a looser one elsewhere, and a zero at a contact point that
    /// has to survive untouched.
    tolerance: alum::VProperty<f64>,

    /// Which vertices may not move, indexed by vertex handle.
    ///
    /// Locking a boundary stops those vertices being *deleted*, because `is_collapse_legal` refuses
    /// a halfedge whose tail carries the feature flag. That is not enough on its own: an interior
    /// vertex may still collapse *into* a locked one, and the surviving vertex is then placed
    /// wherever the placement rule says, which drags the outline off the boundary while leaving the
    /// boundary vertex count untouched. Measured on the bunny, that moved the outline by 0.025 at a
    /// tolerance of 0.02, further than not locking it at all.
    ///
    /// So a locked vertex additionally pins the merged position to itself. See
    /// [`ToleranceVolumeDecimator::candidates`].
    pinned: Vec<bool>,

    /// The structural vetoes and the work counters, shared with the best-effort path.
    structural: Structural,
}

/// The `alum` property key under which the error volume is stored on a mesh.
///
/// Named so that a second decimation run finds the radii the first one left behind rather than
/// starting over from zero, which is what keeps the bound anchored to the true original.
const ERROR_PROP: &str = "engeom_decimate_error";

/// The `alum` property key under which the tolerance volume is stored on a mesh.
const TOLERANCE_PROP: &str = "engeom_decimate_tolerance";

impl ToleranceVolumeDecimator {
    /// Build a decimator for a half-edge mesh.
    ///
    /// Nothing outside the mesh is needed. The error volume already on the mesh, if any, is what
    /// the run measures against; on a mesh that has never been decimated those radii are all zero,
    /// which says the mesh is its own reference.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the half-edge mesh which will be decimated
    /// * `opts`: how decimation should behave
    ///
    /// returns: `Result<ToleranceVolumeDecimator>`
    pub fn new(mesh: &mut HalfEdgeMesh3, opts: DecimateOpts) -> Result<Self> {
        opts.validate()?;

        // A fresh property starts at zero everywhere; an existing one is returned as it stands,
        // carrying whatever a previous run accumulated.
        let error = mesh.vertex_prop_f64(ERROR_PROP, 0.0)?;
        let mut tolerance = mesh.vertex_prop_f64(TOLERANCE_PROP, opts.tol)?;

        {
            let inner = mesh.as_alum();
            let mut tol = tolerance
                .try_borrow_mut()
                .map_err(|_| "Failed to borrow the tolerance volume")?;

            // The budget is set fresh each run, because it is this run's parameter rather than an
            // accumulated quantity. A boundary tolerance is simply a different budget on the
            // boundary vertices, which is cheaper and clearer than branching per collapse.
            for v in inner.vertices() {
                tol[v] = match opts.boundary_tol {
                    Some(t) if v.is_boundary(inner) => t,
                    _ => opts.tol,
                };
            }
        }

        let inner = mesh.as_alum();

        // `alum`'s own `QuadricDecimater` keeps its per-vertex quadrics private, so they are
        // accumulated here the same way it does.
        let pinned = {
            let vstatus = inner.vertex_status_prop();
            let vstatus = vstatus
                .try_borrow()
                .map_err(|_| "Failed to borrow vertex status")?;
            inner.vertices().map(|v| vstatus[v].feature()).collect()
        };

        let mut quadrics = vec![Quadric::default(); inner.num_vertices()];
        for v in inner.vertices() {
            quadrics[v.index() as usize] = quadric_of_vertex(inner, v, opts.quadric)?;
        }

        Ok(Self {
            opts,
            quadrics,
            error,
            tolerance,
            pinned,
            structural: Structural::new(opts.max_normal_deviation, opts.min_aspect),
        })
    }

    /// The candidate positions for the merged vertex, in the order they should be tried.
    fn candidates(&self, q: &Quadric, v2: VH, p1: Vector3, p2: Vector3) -> Vec<Vector3> {
        let pinned = self
            .pinned
            .get(v2.index() as usize)
            .copied()
            .unwrap_or(false);
        candidate_positions(self.opts.placement, pinned, q, p1, p2)
    }

    /// The tolerance budget a collapse of this edge has to work within, and whether the errors
    /// already standing in the star fit inside it.
    ///
    /// The budget is the tightest tolerance anywhere in the edge star, so it propagates downwards
    /// and never upwards. That monotonicity is the whole of the Goal II argument: a point of the
    /// original surface always maps onto a triangle at least one of whose vertices carries a
    /// tolerance no larger than the point's own.
    ///
    /// Radii already in the star have to fit it too, because the budget can only have tightened.
    /// This is what catches a link vertex carrying a large radius from an earlier collapse. Neither
    /// half depends on where the merged vertex ends up, so both are settled once per edge rather
    /// than once per candidate position.
    fn budget(&self, mesh: &AlumPolyMesh, v1: VH, v2: VH, e: &[f64], t: &[f64]) -> Option<f64> {
        let mut budget = t[v1.index() as usize].min(t[v2.index() as usize]);
        for v in mesh.vv_ccw_iter(v1).chain(mesh.vv_ccw_iter(v2)) {
            budget = budget.min(t[v.index() as usize]);
        }

        for v in mesh.vv_ccw_iter(v1).chain(mesh.vv_ccw_iter(v2)) {
            if v != v1 && v != v2 && e[v.index() as usize] > budget {
                return None;
            }
        }

        Some(budget)
    }

    /// Whether this collapse keeps the surfaces inside the tolerance volume.
    ///
    /// This is where an [`ErrorMethod`] is turned into the [`ErrorRule`] which implements it, and
    /// it is the only place that mapping lives. What the rules cannot decide for themselves is what
    /// to do when one of them has no answer, which is the whole of the rest of this function.
    fn accept(&self, c: &Collapse<'_>, budget: f64) -> Option<f64> {
        let bound = match self.opts.error_method {
            ErrorMethod::Contraction => Contraction.bound(c),
            ErrorMethod::AffineFaceMap => AffineFaceMap.bound(c),
            ErrorMethod::ProjectedOverlay => {
                let overlay = ProjectedOverlay.bound(c);

                // On a moving outline the two methods stop being ranked and become two answers to
                // the same question. The overlay is far ahead on a straight run, where the slivers
                // are thin and it charges almost nothing, but the face map can still win on a badly
                // bent one, where a sliver reaches further from the chain than a face does from its
                // own replacement.
                //
                // Both radii are valid on their own, so taking the smaller *picks a method* rather
                // than blending two, which would not be sound. `Unsatisfiable` joins the same
                // reckoning here: it says the overlay's own map cannot be made to work, not that the
                // collapse is impossible, and the face map's map is a different one. This costs a
                // second evaluation only on collapses which have an outline to move.
                if c.outline.parts_regions(c.p0) {
                    let mine = match overlay {
                        ErrorBound::Bound(v) => Some(v),
                        _ => {
                            self.structural.bump_not_applicable();
                            None
                        }
                    };

                    let theirs = match AffineFaceMap.bound(c) {
                        ErrorBound::Bound(v) => Some(v),
                        _ => None,
                    };

                    match mine.into_iter().chain(theirs).reduce(f64::min) {
                        Some(v) => ErrorBound::Bound(v),
                        None => ErrorBound::Unsatisfiable,
                    }
                } else {
                    match overlay {
                        // With the regions equal the overlay's map is a bijection covering both
                        // stars entirely, so a constraint it cannot meet is a fact about the
                        // collapse rather than about the map.
                        ErrorBound::NotApplicable => {
                            self.structural.bump_not_applicable();
                            AffineFaceMap.bound(c)
                        }
                        settled => settled,
                    }
                }
            }
        };

        let error = match bound {
            ErrorBound::Bound(v) => v,
            ErrorBound::Unsatisfiable | ErrorBound::NotApplicable => {
                self.structural.bump_error_veto();
                return None;
            }
        };

        if error > budget {
            self.structural.bump_error_veto();
            return None;
        }

        Some(error)
    }

    /// What this collapse would do, or `None` if no candidate position is acceptable.
    fn evaluate(&self, mesh: &AlumPolyMesh, h: HH) -> Option<Accepted> {
        self.structural.bump_evaluation();
        let v1 = h.tail(mesh);
        let v2 = h.head(mesh);

        let points = mesh.points();
        let points = points.try_borrow().ok()?;

        let e = self.error.try_borrow().ok()?;
        let t = self.tolerance.try_borrow().ok()?;

        let p1 = points[v1];
        let p2 = points[v2];
        let q = self.quadrics[v1.index() as usize] + self.quadrics[v2.index() as usize];

        let star = self.structural.stars(mesh, v1, v2, &points)?;
        if star.is_empty() {
            return None;
        }

        let budget = match self.budget(mesh, v1, v2, &e, &t) {
            Some(b) => b,
            None => {
                self.structural.bump_error_veto();
                return None;
            }
        };

        // The part of the star's outline this collapse can move, which is what decides whether the
        // two stars cover the same region and so whether the projected overlay has anything to say.
        let outline = boundary_chain(&star, v1, v2);

        for candidate in self.candidates(&q, v2, p1, p2) {
            self.structural.bump_acceptance();
            let Some(new) = self.structural.reshape(&star, v1, v2, candidate, &points) else {
                continue;
            };

            let c = Collapse {
                star: &star,
                new: &new,
                points: &points,
                e: &e,
                v1,
                v2,
                p0: candidate,
                outline,
            };

            if let Some(error) = self.accept(&c, budget) {
                return Some(Accepted {
                    position: candidate,
                    quadric: q,
                    cost: q.residual(candidate),
                    error,
                    tolerance: budget,
                });
            }
        }

        None
    }
}

impl CollapseGate for ToleranceVolumeDecimator {
    type Accepted = Accepted;

    fn evaluate(&self, mesh: &AlumPolyMesh, h: HH) -> Option<Accepted> {
        ToleranceVolumeDecimator::evaluate(self, mesh, h)
    }

    fn cost(accepted: &Accepted) -> f64 {
        accepted.cost
    }

    fn position(accepted: &Accepted) -> Vector3 {
        accepted.position
    }

    /// Grow the error volume and tighten the tolerance volume at the surviving vertex.
    ///
    /// The link is left alone. That is what the maps above are for: whichever one ran has already
    /// proved that one ball at the merged vertex absorbs everything the collapse displaced.
    fn commit(&mut self, v: VH, accepted: Accepted) -> Result<()> {
        self.quadrics[v.index() as usize] = accepted.quadric;

        let mut e = self
            .error
            .try_borrow_mut()
            .map_err(|_| "Failed to borrow the error volume")?;
        let mut t = self
            .tolerance
            .try_borrow_mut()
            .map_err(|_| "Failed to borrow the tolerance volume")?;

        e[v] = accepted.error;
        t[v] = accepted.tolerance;

        Ok(())
    }

    fn stats(&self) -> DecimateStats {
        self.structural.stats.get()
    }
}

impl HalfEdgeMesh3 {
    /// Decimate the mesh, keeping both surfaces within the tolerance of each other.
    ///
    /// The bound is two-directional and holds over the surfaces rather than at sampled points: no
    /// part of the result strays further than `opts.tol` from the shape this mesh started as, and
    /// no part of that shape is left further than `opts.tol` from the result.
    ///
    /// Nothing outside the mesh is consulted. The reference is the error volume carried on the
    /// vertices, which is zero on a mesh that has not been decimated and otherwise records
    /// everything that has happened to it, so calling this repeatedly does not walk the surface
    /// away a tolerance at a time. See the `decimate` module documentation for the mechanism.
    ///
    /// # Arguments
    ///
    /// * `opts`: how decimation should behave
    ///
    /// # Errors
    ///
    /// Refuses if [`HalfEdgeMesh3::decimate_best_effort`] has already run on this mesh. That path
    /// leaves the error volume describing less than what has happened to the surface, so a bound
    /// computed from it afterwards would not be one. Round-trip through `Mesh3` and start again if
    /// the best-effort result is deliberately the new reference.
    ///
    /// returns: `Result<DecimateReport>`
    pub fn decimate_guaranteed(&mut self, opts: &DecimateOpts) -> Result<DecimateReport> {
        if self.error_volume_stale() {
            return Err(
                "A guaranteed decimation cannot follow a best-effort one on the same \
                 mesh: the error volume it reasons from no longer accounts for what the \
                 best-effort pass moved, so its bound would not be a bound. Rebuild the half-edge \
                 mesh if that result is meant to be the new reference."
                    .into(),
            );
        }

        opts.validate()?;

        set_boundary_locks(self, opts.lock_boundary)?;

        let mut gate = ToleranceVolumeDecimator::new(self, *opts)?;
        run_decimation(self, &mut gate, opts.target)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::half_edge3::RepairOpts;
    use crate::geom3::half_edge3::difference::ShapeDifference;
    use crate::{Mesh3, Point3};
    use std::sync::OnceLock;

    /// The independent measurement every claim in this module is checked against.
    ///
    /// This consults nothing the decimator tracked and does not sample the sites the decimator
    /// enforces, which is the whole point: a check which samples a gate's own sites agrees with
    /// that gate whether or not the gate is right. See the `difference` module documentation.
    fn deviation(original: &Mesh3, decimated: &Mesh3) -> ShapeDifference {
        ShapeDifference::new(original, decimated).unwrap()
    }

    /// The tolerance the blade tests decimate at, in millimetres.
    ///
    /// Half the blade's own edge length, which the face-to-face methods cannot do anything at all
    /// with and the overlay takes down to under a third of its faces.
    const BLADE_TOL: f64 = 0.05;

    /// An `n` by `n` grid of unit squares in the z=0 plane, each split into two triangles.
    ///
    /// The one fixture whose right answer is known in advance: every collapse within the plane
    /// moves the surface by exactly nothing, so a correct method permits all of them at any
    /// tolerance however tight.
    fn flat_grid(n: usize) -> Mesh3 {
        let mut points = Vec::new();
        for i in 0..=n {
            for j in 0..=n {
                points.push(Point3::new(i as f64, j as f64, 0.0));
            }
        }

        let mut faces = Vec::new();
        let idx = |i: usize, j: usize| (i * (n + 1) + j) as u32;
        for i in 0..n {
            for j in 0..n {
                faces.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
                faces.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
            }
        }

        Mesh3::new(points, faces, false)
    }

    /// Decimate the blade, serving the default run at [`BLADE_TOL`] from a memo.
    ///
    /// Five tests in this module ask for that one run: `the_budget_is_nearly_spent`,
    /// `a_tighter_tolerance_decimates_less`, `the_overlay_is_what_makes_a_tight_tolerance_usable`,
    /// `both_quadric_kinds_respect_the_bound` and `every_method_holds_the_bound`. The last three
    /// reach it through a setter and so do not look like they are asking for the default, but
    /// `ProjectedOverlay` and `QuadricKind::Triangle` are the defaults, so the options compare
    /// equal and the runs are the same run. At roughly eight seconds each that was the largest
    /// single cost in the module.
    ///
    /// Keyed on the whole options struct rather than on the tolerance alone so that a test which
    /// changes any field gets its own run. `DecimateOpts` is `PartialEq`, so this cannot drift out
    /// of agreement with what the callers actually asked for.
    ///
    /// This is a memo of a pure function, not shared state. `decimate_guaranteed` reads nothing a
    /// test could have modified, so there is nothing here for one test to observe of another.
    fn decimate_blade(opts: &DecimateOpts) -> (Mesh3, Mesh3, DecimateReport) {
        static DEFAULT_RUN: OnceLock<(Mesh3, Mesh3, DecimateReport)> = OnceLock::new();

        if *opts == DecimateOpts::to_tolerance(BLADE_TOL) {
            return DEFAULT_RUN
                .get_or_init(|| decimate_blade_uncached(opts))
                .clone();
        }

        decimate_blade_uncached(opts)
    }

    fn decimate_blade_uncached(opts: &DecimateOpts) -> (Mesh3, Mesh3, DecimateReport) {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();
        let report = he.decimate_guaranteed(opts).unwrap();
        let out = Mesh3::try_from(&he).unwrap();
        (original, out, report)
    }

    /// The whole point of the error volume, and the thing no earlier version of this module could
    /// do: both directions of the surface Hausdorff distance inside the tolerance.
    ///
    /// Swept rather than checked at one tolerance because the failure this replaces was not uniform;
    /// the old gate was 8.6x over at 0.005 and 1.8x over at 0.1, so a single sample would have been
    /// easy to pick favourably.
    #[test]
    fn both_directions_are_within_tolerance() {
        for tol in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0] {
            let (original, out, _) = decimate_blade(&DecimateOpts::to_tolerance(tol));
            deviation(&original, &out).assert_within(tol);
        }
    }

    /// Every method holds the bound, and the overlay is the one which makes a tight tolerance
    /// usable.
    ///
    /// The two face-to-face maps cannot collapse an edge shorter than about twice the tolerance,
    /// because they force an old face onto the face in its own slots and an endpoint then has to
    /// reach a triangle which shrank away from it. The overlay maps the stars to each other through
    /// a shared projection instead, which is a bijection, and the restriction goes away.
    #[test]
    fn the_overlay_is_what_makes_a_tight_tolerance_usable() {
        let tol = BLADE_TOL;
        let opts = DecimateOpts::to_tolerance(tol);

        let mut faces = Vec::new();
        for method in [
            ErrorMethod::Contraction,
            ErrorMethod::AffineFaceMap,
            ErrorMethod::ProjectedOverlay,
        ] {
            let (original, out, report) = decimate_blade(&opts.with_error_method(method));
            let d = deviation(&original, &out);
            assert!(d.hausdorff() <= tol, "{:?} broke the bound: {}", method, d);
            faces.push((method, report.faces_after));
        }

        let overlay = faces[2].1;
        for (method, count) in faces.iter().take(2) {
            assert!(
                overlay * 2 < *count,
                "the overlay ({}) should be far ahead of {:?} ({})",
                overlay,
                method,
                count
            );
        }
    }

    /// Every method offered as a bound holds it, swept across tolerances.
    ///
    /// The point of running this per method rather than only on the default is that a bound is the
    /// one property none of them may trade away. Completeness and speed are what the bake-off in
    /// `decimate_profile` decides between; a method that is faster and more complete but lets the
    /// surface out of the tolerance is not in the running at all.
    ///
    /// Gueziec's Subdivision Method was such a method, exceeding its tolerance by up to 3.8x in
    /// the middle of this range while tying on faces removed and costing more time.
    #[test]
    fn every_method_holds_the_bound() {
        for method in [
            ErrorMethod::Contraction,
            ErrorMethod::AffineFaceMap,
            ErrorMethod::ProjectedOverlay,
        ] {
            for tol in [0.01, 0.05, 0.5] {
                let opts = DecimateOpts::to_tolerance(tol).with_error_method(method);
                let (original, out, _) = decimate_blade(&opts);
                let d = deviation(&original, &out);
                assert!(
                    d.hausdorff() <= tol,
                    "{:?} broke the bound at tol={}: {}",
                    method,
                    tol,
                    d
                );
            }
        }
    }

    /// A flat region has an exactly known answer: every in-plane collapse moves the surface by
    /// nothing, so any tolerance however tight should permit every one of them.
    ///
    /// This is the sharpest statement of what a face-to-face correspondence costs. Both of the
    /// other methods refuse every collapse here at a tolerance a tenth of the edge length, which is
    /// not conservative but wrong.
    #[test]
    fn a_flat_region_decimates_at_any_tolerance() {
        let plane = flat_grid(20);

        let mut he = HalfEdgeMesh3::try_from(&plane).unwrap();
        let report = he
            .decimate_guaranteed(&DecimateOpts::to_tolerance(1.0e-6))
            .unwrap();
        let out = Mesh3::try_from(&he).unwrap();

        assert!(
            report.ratio() < 0.2,
            "a plane should flatten at any tolerance, got {}",
            report
        );
        assert!(
            deviation(&plane, &out).hausdorff() < 1.0e-9,
            "collapsing within a plane must not move the surface"
        );
    }

    /// The overlay judges every collapse a locked-boundary run actually performs, with no help from
    /// the fallback.
    ///
    /// A collapse where the surviving vertex is a locked boundary vertex is Gueziec's Type II
    /// boundary edge, which his Section 10 dismisses as presenting "no particular difficulties":
    /// the interior endpoint slides onto the boundary one, so the outline does not move and the two
    /// stars still span the same region. That is precisely the overlay's premise, and this asserts
    /// it holds rather than leaving it as an argument.
    ///
    /// The flat grid is the fixture for it because a fifth of its vertices are on the outline and
    /// its zero curvature rules out folding, so any fallback here would have to be the boundary. The
    /// count was 3507 before locked vertices stopped being costed as collapse tails, and every one
    /// of those was on a collapse `is_collapse_legal` refuses anyway.
    #[test]
    fn a_locked_boundary_never_reaches_the_fallback() {
        let plane = flat_grid(20);
        let mut he = HalfEdgeMesh3::try_from(&plane).unwrap();
        let report = he
            .decimate_guaranteed(&DecimateOpts::to_tolerance(1.0e-6))
            .unwrap();

        assert_eq!(
            report.stats.method_not_applicable, 0,
            "the overlay should judge every one of these unaided: {}",
            report.stats
        );
        assert!(report.ratio() < 0.2, "expected real decimation: {}", report);
    }

    /// The budget should be nearly spent, and real decimation should happen while spending it.
    ///
    /// Both halves matter and neither is sufficient. A bound is trivially held by refusing
    /// everything, and the earlier face-to-face methods used only a third of what they were given;
    /// the overlay runs at 80% to 100% of the tolerance, which is what says the tolerance is the
    /// thing actually limiting decimation rather than an artifact of the correspondence.
    #[test]
    fn the_budget_is_nearly_spent() {
        let (original, out, report) = decimate_blade(&DecimateOpts::to_tolerance(BLADE_TOL));
        let d = deviation(&original, &out);

        let ratio = d.hausdorff() / BLADE_TOL;
        assert!(ratio <= 1.0, "over the bound: {}", d);
        assert!(ratio > 0.5, "leaving most of the budget unused: {}", d);
        assert!(
            report.ratio() < 0.85,
            "expected real decimation, got {}",
            report
        );
    }

    /// A tighter tolerance must give both a smaller deviation and a finer mesh, or the bound is
    /// not actually driving anything.
    #[test]
    fn a_tighter_tolerance_decimates_less() {
        let (original, coarse, coarse_report) =
            decimate_blade(&DecimateOpts::to_tolerance(BLADE_TOL));
        let (_, fine, fine_report) = decimate_blade(&DecimateOpts::to_tolerance(BLADE_TOL / 5.0));

        assert!(fine_report.faces_after > coarse_report.faces_after);
        assert!(deviation(&original, &fine).hausdorff() <= BLADE_TOL / 5.0);
        assert!(deviation(&original, &coarse).hausdorff() <= BLADE_TOL);
    }

    /// With endpoint placement every surviving point is one that was measured, which is the reason
    /// to choose it.
    #[test]
    fn endpoint_placement_keeps_the_points_a_subset() {
        let original = Mesh3::create_sphere(10.0, 0.06).unwrap();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        he.decimate_guaranteed(
            &DecimateOpts::to_tolerance(0.5).with_placement(Placement::Endpoint),
        )
        .unwrap();
        let out = Mesh3::try_from(&he).unwrap();

        assert!(out.points().len() < original.points().len());
        for p in out.points() {
            let d = original
                .points()
                .iter()
                .map(|o| (o - p).norm())
                .fold(f64::INFINITY, f64::min);
            assert!(d < 1e-12, "point {:?} was not one of the originals", p);
        }
    }

    /// Locking the boundary keeps the outline of an open mesh from moving.
    ///
    /// The vertex count is the weaker half of this and it is not sufficient. Keeping every boundary
    /// vertex while letting the outline move is a real failure mode and was the actual behaviour
    /// until the merged position was pinned: an interior vertex collapsing *into* a locked one left
    /// the count untouched and dragged the outline 0.025 off a tolerance of 0.02, further than not
    /// locking at all. The geometric assertion is what catches that.
    #[test]
    fn a_locked_boundary_is_preserved() {
        let original = crate::tests::stanford_bun_4();
        let (mut he, _) =
            HalfEdgeMesh3::from_mesh_repaired(&original, &RepairOpts::default()).unwrap();
        let before = Mesh3::try_from(&he).unwrap();
        let boundary_before = before.compute_nav().boundary_vertices(None).len();

        he.decimate_guaranteed(&DecimateOpts::to_tolerance(0.02))
            .unwrap();
        let after = Mesh3::try_from(&he).unwrap();
        let boundary_after = after.compute_nav().boundary_vertices(None).len();

        assert_eq!(
            boundary_before, boundary_after,
            "a locked boundary should not lose vertices"
        );

        let moved = deviation(&before, &after).boundary.unwrap();
        assert!(
            moved < 1.0e-12,
            "a locked boundary should not move at all, but the outline shifted by {}",
            moved
        );
    }

    /// The boundary lock is per-run state rather than a ratchet: a later run asking for the
    /// boundary free must actually free it, even though an earlier run locked it on this same mesh.
    ///
    /// This was real. The lock was applied by setting the alum feature flag and nothing ever
    /// cleared it, so the second run below kept the outline fully locked whatever its options
    /// said, and no test caught it because every test built a fresh mesh.
    #[test]
    fn a_later_run_can_unlock_a_previously_locked_boundary() {
        let plane = flat_grid(12);
        let mut he = HalfEdgeMesh3::try_from(&plane).unwrap();

        he.decimate_guaranteed(&DecimateOpts::to_tolerance(1.0e-6))
            .unwrap();
        let locked = Mesh3::try_from(&he).unwrap();
        let boundary_locked = locked.compute_nav().boundary_vertices(None).len();

        he.decimate_guaranteed(&DecimateOpts::to_tolerance(1.0e-6).with_lock_boundary(false))
            .unwrap();
        let unlocked = Mesh3::try_from(&he).unwrap();
        let boundary_unlocked = unlocked.compute_nav().boundary_vertices(None).len();

        assert!(
            boundary_unlocked < boundary_locked,
            "the second run asked for a free boundary and should have simplified the outline, \
             but it stayed at {} vertices",
            boundary_locked
        );
    }

    #[test]
    fn an_unlocked_boundary_can_collapse() {
        let original = crate::tests::stanford_bun_4();
        let (mut he, _) =
            HalfEdgeMesh3::from_mesh_repaired(&original, &RepairOpts::default()).unwrap();
        let before = Mesh3::try_from(&he).unwrap();
        let boundary_before = before.compute_nav().boundary_vertices(None).len();

        he.decimate_guaranteed(&DecimateOpts::to_tolerance(0.02).with_lock_boundary(false))
            .unwrap();
        let after = Mesh3::try_from(&he).unwrap();
        let boundary_after = after.compute_nav().boundary_vertices(None).len();

        assert!(boundary_after < boundary_before);
    }

    /// Unlocking the boundary must not cost the bound, and it did.
    ///
    /// This is the configuration where the projected overlay's premise fails: a collapse which
    /// deletes a boundary vertex, or which drags a surviving one, moves the star's outline, so the
    /// two stars cover different regions and the shared projection is not a bijection between them.
    /// The part of the new star with no old counterpart is then covered by no constraint and never
    /// charged for, and the radius comes out too small. Nothing throws; the bound is simply not a
    /// bound.
    ///
    /// It was 1.9x over at the tightest tolerance here before `outline_moves` sent those collapses
    /// to the face-to-face map, which never assumed a shared region. Both of the other methods were
    /// always inside it, which is what identified the overlay as the one at fault.
    ///
    /// Swept over all three methods because "holds its bound" is the property none of them may
    /// trade away, and over three tolerances because the failure was not uniform: the same run was
    /// 1.9x over at 0.005 and comfortably inside at 0.02 and 0.1, so a single sample would have
    /// missed it.
    ///
    /// `an_unlocked_boundary_can_collapse` only shows the outline loses vertices, which a method
    /// that moved the surface anywhere at all would also satisfy.
    #[test]
    fn an_unlocked_boundary_still_holds_the_bound() {
        let original = crate::tests::stanford_bun_4();

        for method in [
            ErrorMethod::Contraction,
            ErrorMethod::AffineFaceMap,
            ErrorMethod::ProjectedOverlay,
        ] {
            for tol in [0.005, 0.02, 0.1] {
                let (mut he, _) =
                    HalfEdgeMesh3::from_mesh_repaired(&original, &RepairOpts::default()).unwrap();
                let before = Mesh3::try_from(&he).unwrap();

                let report = he
                    .decimate_guaranteed(
                        &DecimateOpts::to_tolerance(tol)
                            .with_lock_boundary(false)
                            .with_error_method(method),
                    )
                    .unwrap();
                let after = Mesh3::try_from(&he).unwrap();

                assert!(
                    report.ratio() < 0.9,
                    "{:?} at tol={} did nothing: {}",
                    method,
                    tol,
                    report
                );
                let d = deviation(&before, &after);
                assert!(
                    d.hausdorff() <= tol,
                    "{:?} broke the bound at tol={} with the boundary unlocked: {}",
                    method,
                    tol,
                    d
                );
            }
        }
    }

    /// Running repeatedly must not let the error creep past the tolerance.
    ///
    /// The error volume is what makes this hold, and it holds structurally rather than by anything
    /// this test's arrangement arranges. The radii live on the mesh, so the second run starts from
    /// what the first accumulated and the budget it is checked against has not moved. Nothing here
    /// re-consults the original mesh, which is the point: a run has no way to measure itself
    /// against the previous run's output even if it wanted to.
    #[test]
    fn repeated_runs_do_not_compound_error() {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        for _ in 0..3 {
            he.decimate_guaranteed(&DecimateOpts::to_tolerance(BLADE_TOL))
                .unwrap();
        }

        let out = Mesh3::try_from(&he).unwrap();
        let d = deviation(&original, &out);
        assert!(d.hausdorff() <= BLADE_TOL, "three runs drifted to {}", d);
    }

    /// A second run must also not simply repeat the first, or the persistence above would be
    /// untested: a decimator which refused everything after the first pass would satisfy the bound.
    #[test]
    fn a_second_run_continues_from_the_first() {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        let first = he
            .decimate_guaranteed(&DecimateOpts::to_tolerance(BLADE_TOL))
            .unwrap();
        let second = he
            .decimate_guaranteed(&DecimateOpts::to_tolerance(BLADE_TOL))
            .unwrap();

        assert_eq!(second.faces_before, first.faces_after);
        assert!(
            second.faces_after <= second.faces_before,
            "a second run grew the mesh: {}",
            second
        );
    }

    #[test]
    fn invalid_options_are_rejected() {
        let original = crate::tests::engine_blade();
        let mut he = HalfEdgeMesh3::try_from(&original).unwrap();

        assert!(
            he.decimate_guaranteed(&DecimateOpts::to_tolerance(0.0))
                .is_err()
        );
        assert!(
            he.decimate_guaranteed(&DecimateOpts::to_tolerance(-1.0))
                .is_err()
        );
        assert!(
            he.decimate_guaranteed(&DecimateOpts::to_ratio(0.01, 0.0))
                .is_err()
        );
        assert!(
            he.decimate_guaranteed(&DecimateOpts::to_ratio(0.01, 1.5))
                .is_err()
        );
    }

    /// Both quadric formulations must respect the bound; they differ in ordering, not in what is
    /// allowed.
    #[test]
    fn both_quadric_kinds_respect_the_bound() {
        for kind in [QuadricKind::Triangle, QuadricKind::Probabilistic] {
            let (original, out, report) =
                decimate_blade(&DecimateOpts::to_tolerance(BLADE_TOL).with_quadric(kind));

            let d = deviation(&original, &out);
            assert!(d.hausdorff() <= BLADE_TOL, "{:?} drifted to {}", kind, d);
            assert!(report.ratio() < 0.95, "{:?} did nothing: {}", kind, report);
        }
    }

    /// Not a test but a measurement harness, so it is ignored by default. Prints the tolerance
    /// sweep with the independently measured deviation and the work counters, which is the
    /// baseline every performance change is judged against.
    ///
    /// `cargo test -r -p engeom --lib decimate_profile -- --ignored --nocapture`
    #[test]
    #[ignore = "measurement harness, not a correctness test"]
    fn decimate_profile() {
        let original = crate::tests::engine_blade();
        println!(
            "\nengine-blade: {} faces, {} points\n",
            original.faces().len(),
            original.points().len()
        );

        for tol in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5] {
            for method in [
                ErrorMethod::Contraction,
                ErrorMethod::AffineFaceMap,
                ErrorMethod::ProjectedOverlay,
            ] {
                println!("--- {:?}", method);
                let mut he = HalfEdgeMesh3::try_from(&original).unwrap();
                let start = std::time::Instant::now();
                let report = he
                    .decimate_guaranteed(&DecimateOpts::to_tolerance(tol).with_error_method(method))
                    .unwrap();
                let elapsed = start.elapsed();

                let out = Mesh3::try_from(&he).unwrap();
                let d = deviation(&original, &out);
                let verdict = if d.hausdorff() <= tol { "OK  " } else { "OVER" };

                let (rev, fwd) = (d.dense.reference_to_test, d.dense.test_to_reference);
                println!(
                    "{verdict} tol={tol:<7} {report}  rev={rev:.6} ({:.1}x) fwd={fwd:.6} ({:.1}x)  \
                 [{elapsed:?}]",
                    rev / tol,
                    fwd / tol
                );
                println!("       {}", report.stats);
                println!(
                    "       acceptance tests per collapse: {:.1}, evaluations per collapse: {:.1}",
                    report.stats.acceptance_tests as f64 / report.collapses.max(1) as f64,
                    report.stats.evaluations as f64 / report.collapses.max(1) as f64,
                );
            }
        }
    }

    /// How often the projected overlay has nothing to say, and the face-to-face map decides instead.
    ///
    /// When the overlay returns [`ErrorBound::NotApplicable`] that is not a refusal:
    /// [`ErrorMethod::AffineFaceMap`] takes over and judges the collapse by the much looser
    /// face-to-face correspondence. The entire three-way [`ErrorMethod`] split exists because the
    /// overlay can fail to have an answer, so this measures how much of that machinery is servicing
    /// a real case.
    ///
    /// # A correction, and what the rate is really about
    ///
    /// This harness previously reported the following, and the conclusion drawn from it was wrong:
    ///
    /// ```text
    /// mesh              boundary vertices    fell back
    /// engine-blade                  0.0%     0.06% to 0.46%
    /// stanford-bunny               14.7%     6.6% to 8.2%
    /// flat-grid-20                 18.1%     17.8%
    /// flat-grid-60                  6.4%     5.4%
    /// ```
    ///
    /// The flat grids have exactly zero curvature, so nothing can fold, yet they had the highest
    /// rates; the rate tracked the boundary fraction almost one for one. That looked like proof
    /// that the overlay's premise fails on a boundary, since it assumes both stars are disks across
    /// the same link polygon.
    ///
    /// It was an artifact. Breaking the count apart by cause showed every one of the flat grids'
    /// fallbacks, and a third of the bunny's, came from a collapse whose *tail* was a locked
    /// boundary vertex. `is_collapse_legal` refuses those outright, so they never happened; they
    /// were being costed and thrown away. On flat-grid-20 that was 53% of every acceptance test the
    /// run performed.
    ///
    /// `queue_vertex` now skips a locked vertex before costing anything, and what is left is the
    /// bare rate:
    ///
    /// ```text
    /// mesh              boundary vertices    fell back      cause
    /// engine-blade                  0.0%     0.06% to 0.46% folding
    /// stanford-bunny               14.7%     5.1% to 6.1%   folding
    /// flat-grid-20                 18.1%     0.0%           -
    /// flat-grid-60                  6.4%     0.0%           -
    /// ```
    ///
    /// **The boundary was never the problem, and folding was.** A locked-boundary run only ever
    /// collapses an interior vertex into its neighbour, and where that neighbour is on the outline
    /// it is pinned there, so the outline does not move and the two stars do span the same region
    /// after all. That is Gueziec's Type II boundary edge, which his Section 10 calls presenting
    /// "no particular difficulties", and the overlay handles it unaided; the flat grid columns above
    /// are the evidence and `a_locked_boundary_never_reaches_the_fallback` is the regression test.
    ///
    /// The residual is folding on a curved star, which is where the shared projection stops being
    /// injective. It is negligible on the ruled blade and around five percent on the bunny, which is
    /// curved in every direction. So a correspondence more robust to folding has a little to win on
    /// organic data and essentially nothing on the sort of surface this library is aimed at.
    ///
    /// A collapse which *deletes* a boundary vertex does move the outline, and there the premise
    /// really does fail. That needs `lock_boundary` off to arise at all, and is measured separately
    /// below.
    ///
    /// `cargo test -r -p engeom --lib fallback_rate -- --ignored --nocapture`
    #[test]
    #[ignore = "measurement harness, not a correctness test"]
    fn fallback_rate() {
        fn report_on(label: &str, build: impl Fn() -> HalfEdgeMesh3, tol: f64) {
            run(label, build, tol, true);
        }

        fn run(label: &str, build: impl Fn() -> HalfEdgeMesh3, tol: f64, lock: bool) {
            let mut he = build();

            // Printed alongside rather than inferred, since this is the quantity the rate was once
            // thought to track.
            let before = Mesh3::try_from(&he).unwrap();
            let open = before.compute_nav().boundary_vertices(None).len();
            let open = 100.0 * open as f64 / before.points().len().max(1) as f64;

            let report = he
                .decimate_guaranteed(&DecimateOpts::to_tolerance(tol).with_lock_boundary(lock))
                .unwrap();
            let s = report.stats;

            // Per acceptance test rather than per evaluation: the fallback is decided once for each
            // candidate position tried, not once per edge considered.
            let rate = 100.0 * s.method_not_applicable as f64 / s.acceptance_tests.max(1) as f64;
            println!(
                "{label:<16} tol={tol:<8} lock={lock:<5} boundary {open:>5.1}%   {report}\n  \
                 {} of {} acceptance tests fell back ({rate:.3}%)",
                s.method_not_applicable, s.acceptance_tests,
            );
        }

        println!();
        let blade = crate::tests::engine_blade();
        for tol in [0.001, 0.005, 0.01, 0.05, 0.1, 0.5] {
            report_on(
                "engine-blade",
                || HalfEdgeMesh3::try_from(&blade).unwrap(),
                tol,
            );
        }

        // An organic surface, curved in every direction rather than ruled like the blade, which is
        // where a shared projection has the most opportunity to fold. It needs repairing first.
        let bunny = crate::tests::stanford_bun_4();
        for tol in [0.005, 0.02, 0.1] {
            report_on(
                "stanford-bunny",
                || {
                    HalfEdgeMesh3::from_mesh_repaired(&bunny, &RepairOpts::default())
                        .unwrap()
                        .0
                },
                tol,
            );
        }

        // Zero curvature, so nothing can fold and the rate has to be zero unless the boundary is
        // costing something. Two sizes, because they differ in one thing only: the share of the
        // mesh that is boundary. These are what the correction above rests on.
        for n in [20usize, 60] {
            let plane = flat_grid(n);
            report_on(
                &format!("flat-grid-{n}"),
                || HalfEdgeMesh3::try_from(&plane).unwrap(),
                1.0e-6,
            );
        }

        // Unlocking the boundary lets a collapse move the outline, so the two stars no longer span
        // the same region and the shared projection is not a bijection between them. The overlay
        // covers the leftover slivers instead of declining, and runs against the face-to-face map
        // with the smaller radius kept, so the counter here means something different from the rows
        // above: it counts the collapses where the overlay had nothing to offer and the face map
        // decided alone, rather than collapses judged by a looser method.
        println!("\n--- boundary unlocked, which is Gueziec's Type I collapse");
        for n in [20usize, 60] {
            let plane = flat_grid(n);
            run(
                &format!("flat-grid-{n}"),
                || HalfEdgeMesh3::try_from(&plane).unwrap(),
                1.0e-6,
                false,
            );
        }
        for tol in [0.005, 0.02, 0.1] {
            run(
                "stanford-bunny",
                || {
                    HalfEdgeMesh3::from_mesh_repaired(&bunny, &RepairOpts::default())
                        .unwrap()
                        .0
                },
                tol,
                false,
            );
        }
    }

    /// How much of the tolerance budget is actually spent, per vertex.
    ///
    /// The error volume is an over-estimate, and this says by how much: it prints the distribution
    /// of the radii the run left behind against the deviation that was actually achieved. A method
    /// which improves the constraint solving should show the radii falling towards the measured
    /// deviation while the measured deviation itself barely moves.
    ///
    /// `cargo test -r -p engeom --lib error_volume_slack -- --ignored --nocapture`
    #[test]
    #[ignore = "measurement harness, not a correctness test"]
    fn error_volume_slack() {
        let original = crate::tests::engine_blade();

        for tol in [0.1, 0.5, 1.0] {
            let mut he = HalfEdgeMesh3::try_from(&original).unwrap();
            let report = he
                .decimate_guaranteed(&DecimateOpts::to_tolerance(tol))
                .unwrap();

            let mut radii = he.vertex_prop_f64(ERROR_PROP, 0.0).unwrap();
            let mut radii: Vec<f64> = {
                let r = radii.try_borrow_mut().unwrap();
                r.iter().copied().collect()
            };
            radii.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let at = |q: f64| radii[((radii.len() - 1) as f64 * q) as usize];
            let out = Mesh3::try_from(&he).unwrap();
            let d = deviation(&original, &out);

            println!(
                "tol={tol:<5} {report}\n       radii: median {:.6} p90 {:.6} max {:.6}  \
                 measured {:.6} ({:.0}% of the worst radius)",
                at(0.5),
                at(0.9),
                at(1.0),
                d.hausdorff(),
                100.0 * d.hausdorff() / at(1.0).max(f64::MIN_POSITIVE)
            );
        }
    }
}
