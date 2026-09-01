//! Which part of a star's outline a collapse can move.
//!
//! This is a fact about the collapse rather than about any error method, and every method needs it
//! in the same form. Both of the ones which build a correspondence between the two stars rest on
//! the premise that the stars cover the *same* region, and that premise is what a mesh boundary
//! running through the collapsed edge breaks. See the parent module's documentation for what goes
//! wrong when it isn't explicitly checked, and for Gueziec's edge types.

use super::StarFace;
use crate::{Point3, Vector3};
use alum::VH;

/// The run of mesh boundary edges passing through the collapsed edge's endpoints.
///
/// This is the only part of a star's outline a collapse can move, and knowing where it is turns the
/// boundary case from something to refuse into something to reason about. The rest of the outline is
/// the link, whose vertices the collapse neither moves nor deletes.
///
/// At most four vertices: the two link vertices anchoring the run, and the one or two collapsed
/// endpoints between them. Held inline rather than in a `Vec` because this is built once per
/// candidate edge on the hot path.
#[derive(Debug, Clone, Copy)]
pub(super) struct Chain {
    pub(super) v: [VH; 4],
    pub(super) p: [Point3; 4],
    pub(super) len: usize,
}

impl Chain {
    /// The vertices in order along the outline.
    pub(super) fn iter(&self) -> impl Iterator<Item = (VH, Point3)> + '_ {
        (0..self.len).map(|i| (self.v[i], self.p[i]))
    }

    /// Whether the collapse actually moves the run, which it does unless every vertex it deletes or
    /// relocates ends up exactly where it already was.
    ///
    /// Exact equality rather than a tolerance: the pinned case produces the position itself, and
    /// anything else has moved the outline by however little, which is the direction it is safe to
    /// be strict about.
    pub(super) fn moves_under(&self, p0: Vector3) -> bool {
        self.iter()
            .skip(1)
            .take(self.len - 2)
            .any(|(_, q)| q != Point3::from(p0))
    }
}

/// What a star's outline is made of, which decides whether the two stars cover the same region.
#[derive(Debug, Clone, Copy)]
pub(super) enum Outline {
    /// No mesh boundary edge touches either endpoint, so the outline is the link polygon alone and
    /// nothing the collapse does can move it. This is every collapse on a watertight mesh.
    Closed,

    /// The outline runs through the endpoints along this chain.
    Open(Chain),

    /// A configuration this module declines to judge.
    ///
    /// Two separate runs is Gueziec's Type III edge, where the collapsed edge joins two different
    /// parts of the outline; collapsing it merges two boundaries and leaves a singular vertex, and
    /// he refuses it outright. Anything else which is not a single simple run is refused for the
    /// same reason it cannot be reasoned about: the shape of the region change is not known.
    Unsupported,
}

impl Outline {
    /// Whether this placement parts the two stars' regions, so that the overlay has to cover the
    /// leftovers rather than resting on a bijection between them.
    pub(super) fn parts_regions(&self, p0: Vector3) -> bool {
        match self {
            Outline::Closed => false,
            Outline::Open(chain) => chain.moves_under(p0),
            Outline::Unsupported => true,
        }
    }
}

/// Find the star's boundary chain, using only the star itself.
///
/// The test needs no access to the mesh. An interior mesh edge touching one of the endpoints is
/// carried by two faces and *both* are in the star, since both touch that endpoint; a mesh boundary
/// edge is carried by one. So among the edges touching an endpoint, appearing in exactly one star
/// face is the same statement as being on the mesh boundary.
///
/// Edges of the link are left alone even where they are on the mesh boundary themselves, which they
/// can be: they do not touch either endpoint, so the collapse does not move them.
pub(super) fn boundary_chain(star: &[StarFace], v1: VH, v2: VH) -> Outline {
    // Every star edge touching an endpoint, with how many star faces carry it.
    let mut edges: Vec<(VH, VH, usize)> = Vec::with_capacity(8);

    for f in star.iter() {
        for k in 0..3 {
            let (hk, hj) = (f.v[k], f.v[(k + 1) % 3]);
            if !(hk == v1 || hk == v2 || hj == v1 || hj == v2) {
                continue;
            }
            match edges
                .iter_mut()
                .find(|(a, b, _)| (*a == hk && *b == hj) || (*a == hj && *b == hk))
            {
                Some((_, _, n)) => *n += 1,
                None => edges.push((hk, hj, 1)),
            }
        }
    }

    edges.retain(|(_, _, n)| *n == 1);
    if edges.is_empty() {
        return Outline::Closed;
    }

    // A run through one or two endpoints has at most three edges, so anything longer is a shape this
    // does not describe.
    if edges.len() > 3 {
        return Outline::Unsupported;
    }

    let mut all: Vec<VH> = Vec::with_capacity(4);
    for (a, b, _) in edges.iter() {
        for v in [*a, *b] {
            if !all.contains(&v) {
                all.push(v);
            }
        }
    }

    // A single simple run has exactly two vertices of degree one, and everything between them is an
    // endpoint being merged. A link vertex of degree two would mean the run doubles back through a
    // vertex the collapse does not touch, which is not a shape the chain after the collapse can
    // describe.
    let mut ends: Vec<VH> = Vec::with_capacity(2);
    for v in all.iter() {
        let d = edges.iter().filter(|(a, b, _)| a == v || b == v).count();
        match d {
            1 => ends.push(*v),
            2 if *v == v1 || *v == v2 => {}
            _ => return Outline::Unsupported,
        }
    }
    if ends.len() != 2 {
        return Outline::Unsupported;
    }

    let mut chain = Chain {
        v: [VH::from(0u32); 4],
        p: [Point3::origin(); 4],
        len: 0,
    };
    let mut used = vec![false; edges.len()];
    let mut at = ends[0];

    loop {
        let Some(p) = star_position(star, at) else {
            return Outline::Unsupported;
        };
        chain.v[chain.len] = at;
        chain.p[chain.len] = p;
        chain.len += 1;

        let step = edges.iter().enumerate().find_map(|(i, (a, b, _))| {
            match (used[i], *a == at, *b == at) {
                (false, true, _) => Some((i, *b)),
                (false, _, true) => Some((i, *a)),
                _ => None,
            }
        });

        match step {
            Some((i, next)) => {
                used[i] = true;
                at = next;
            }
            None => break,
        }
    }

    // A disjoint cycle alongside the run would leave edges unwalked, and the degree test above does
    // not rule one out on its own.
    if chain.len != all.len() || used.iter().any(|u| !u) {
        return Outline::Unsupported;
    }

    Outline::Open(chain)
}

/// Where a vertex of the star sits, found from any face carrying it.
fn star_position(star: &[StarFace], v: VH) -> Option<Point3> {
    star.iter()
        .find_map(|f| f.v.iter().position(|h| *h == v).map(|k| f.p[k]))
}
