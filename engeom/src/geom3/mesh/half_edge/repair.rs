//! This module makes a triangle soup manifold enough for a half-edge structure to accept it.
//!
//! `alum`'s `add_face` refuses non-manifold topology, returning `ComplexVertex` or
//! `ComplexHalfedge` and leaving the mesh half-built, so a conversion which feeds it raw scan data
//! fails partway through with no account of what was wrong. Real measured data is routinely
//! non-manifold: a laser line that grazes an edge leaves two sheets joined at a single vertex, a
//! registration seam leaves an edge shared by three faces, and any exporter can emit a triangle
//! twice.
//!
//! So the repair runs on the point and face buffers *before* the half-edge structure sees them,
//! where a bad face can simply be dropped rather than aborting a partial build. Every pass records
//! what it changed into a [`RepairReport`]. A metrology library which silently deletes measured
//! faces is worse than one which refuses to load them, so nothing here is allowed to happen
//! quietly.
//!
//! The passes run in a fixed order, because each assumes the previous has finished:
//!
//! 1. Degenerate and duplicate faces, since they corrupt the edge counts every later pass reads.
//! 2. Non-manifold edges, bringing every edge down to at most two faces.
//! 3. Orientation, so each connected component agrees on which side is out.
//! 4. Any directed edge still used twice, which orientation could not resolve. Two faces walking
//!    an edge the same way cannot both survive in an orientable manifold, and this is the pass
//!    that makes the buffers something a half-edge structure will actually accept.
//! 5. Bowtie vertices, *after* all face removal, because dropping a face can itself split a fan
//!    in two and create a bowtie that was not there before.
//! 6. Isolated vertices, once nothing else will drop a face.
//!
//! Passes 2 and 4 both choose which faces to sacrifice by triangle quality, keeping the better
//! surface and discarding slivers. Pass 2 additionally prefers to keep a pair of faces which
//! traverse the shared edge in opposite directions, since a pair which agrees on direction cannot
//! be oriented consistently and would only push the problem into pass 4.

use crate::common::IndexMask;
use crate::geom3::mesh::algorithms::normals::compute_face_normal;
use crate::geom3::mesh::edges::edge_key;
use crate::geom3::mesh::half_edge::HalfEdgeMesh;
use crate::{Mesh3, Point3, Result};
use faer::prelude::default;
use parry3d_f64::utils::hashmap::HashMap;
use parry3d_f64::utils::hashset::HashSet;

/// Which repairs to attempt when making a mesh manifold.
///
/// The default enables every pass, which is what feeding real scan data to a half-edge structure
/// needs. Turn one off when you would rather the conversion fail than have that particular thing
/// changed under you.
///
/// Marked `#[non_exhaustive]`, so build one with the chained setters rather than a struct literal:
///
/// ```
/// use engeom::geom3::mesh::half_edge::RepairOpts;
///
/// let opts = RepairOpts::default().with_orient_consistently(false);
/// ```
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepairOpts {
    /// Drop faces which repeat a vertex or enclose no area.
    pub drop_degenerate: bool,

    /// Drop faces which repeat a triangle already present, in either winding.
    pub drop_duplicate_faces: bool,

    /// Drop faces until no edge is shared by more than two of them.
    pub resolve_nonmanifold_edges: bool,

    /// Split a vertex whose incident faces form more than one fan into one vertex per fan.
    pub split_bowtie_vertices: bool,

    /// Flip faces so that every face in a connected component agrees on which side is out.
    pub orient_consistently: bool,

    /// Drop points which no surviving face references.
    pub drop_isolated_vertices: bool,
}

impl Default for RepairOpts {
    fn default() -> Self {
        Self {
            drop_degenerate: true,
            drop_duplicate_faces: true,
            resolve_nonmanifold_edges: true,
            split_bowtie_vertices: true,
            orient_consistently: true,
            drop_isolated_vertices: true,
        }
    }
}

impl RepairOpts {
    /// Every pass disabled, as a base for turning on only what you want.
    pub fn none() -> Self {
        Self {
            drop_degenerate: false,
            drop_duplicate_faces: false,
            resolve_nonmanifold_edges: false,
            split_bowtie_vertices: false,
            orient_consistently: false,
            drop_isolated_vertices: false,
        }
    }

    /// Set whether degenerate faces are dropped, returning the modified options.
    pub fn with_drop_degenerate(mut self, value: bool) -> Self {
        self.drop_degenerate = value;
        self
    }

    /// Set whether duplicate faces are dropped, returning the modified options.
    pub fn with_drop_duplicate_faces(mut self, value: bool) -> Self {
        self.drop_duplicate_faces = value;
        self
    }

    /// Set whether non-manifold edges are resolved, returning the modified options.
    pub fn with_resolve_nonmanifold_edges(mut self, value: bool) -> Self {
        self.resolve_nonmanifold_edges = value;
        self
    }

    /// Set whether bowtie vertices are split, returning the modified options.
    pub fn with_split_bowtie_vertices(mut self, value: bool) -> Self {
        self.split_bowtie_vertices = value;
        self
    }

    /// Set whether faces are oriented consistently, returning the modified options.
    pub fn with_orient_consistently(mut self, value: bool) -> Self {
        self.orient_consistently = value;
        self
    }

    /// Set whether isolated vertices are dropped, returning the modified options.
    pub fn with_drop_isolated_vertices(mut self, value: bool) -> Self {
        self.drop_isolated_vertices = value;
        self
    }
}

/// An account of everything a repair changed.
///
/// Every field is a count of measured data that was altered or discarded, so a caller can decide
/// whether the result is still worth trusting. `is_clean` is the quick check that nothing happened
/// at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RepairReport {
    /// Faces dropped for repeating a vertex or enclosing no area.
    pub degenerate_removed: usize,

    /// Faces dropped for repeating a triangle already present, in either winding.
    pub duplicate_faces_removed: usize,

    /// Edges found to be shared by more than two faces.
    pub nonmanifold_edges: usize,

    /// Faces dropped in order to bring those edges down to two faces each.
    pub faces_dropped_at_nonmanifold: usize,

    /// Faces dropped because they walked an edge in the same direction as a face which was kept.
    ///
    /// Two such faces cannot both survive in an orientable manifold. A non-zero count here means
    /// the input had a fold or a non-orientable region that reorienting could not resolve.
    pub faces_dropped_for_orientation: usize,

    /// Vertices whose incident faces formed more than one fan, and were split apart.
    pub bowtie_vertices_split: usize,

    /// Points added by that splitting, which is one fewer than the fans at each split vertex.
    pub points_added_by_splitting: usize,

    /// Faces flipped to agree with the rest of their connected component.
    pub faces_reoriented: usize,

    /// Connected components which could not be oriented consistently, because they are
    /// non-orientable in the way a Moebius strip is. Their faces are left as they were.
    pub nonorientable_components: usize,

    /// Points dropped because no surviving face referenced them.
    pub isolated_vertices_removed: usize,

    /// Faces the half-edge structure rejected even after every pass had run.
    ///
    /// This should be zero. A non-zero value means the repair missed a case, and is reported
    /// rather than raised so that a mostly-good mesh still loads.
    pub faces_rejected_by_ingest: usize,
}

impl RepairReport {
    /// Whether the repair changed nothing at all.
    pub fn is_clean(&self) -> bool {
        *self == Self::default()
    }

    /// The total number of faces the repair discarded.
    pub fn faces_removed(&self) -> usize {
        self.degenerate_removed
            + self.duplicate_faces_removed
            + self.faces_dropped_at_nonmanifold
            + self.faces_dropped_for_orientation
            + self.faces_rejected_by_ingest
    }
}

impl std::fmt::Display for RepairReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_clean() {
            return write!(f, "mesh was already manifold, nothing changed");
        }

        let mut parts = Vec::new();
        let mut push = |n: usize, label: &str| {
            if n > 0 {
                parts.push(format!("{} {}", n, label));
            }
        };

        push(self.degenerate_removed, "degenerate faces removed");
        push(self.duplicate_faces_removed, "duplicate faces removed");
        push(self.nonmanifold_edges, "non-manifold edges");
        push(self.faces_dropped_at_nonmanifold, "faces dropped there");
        push(
            self.faces_dropped_for_orientation,
            "faces dropped as unorientable",
        );
        push(self.bowtie_vertices_split, "bowtie vertices split");
        push(self.points_added_by_splitting, "points added");
        push(self.faces_reoriented, "faces reoriented");
        push(self.nonorientable_components, "non-orientable components");
        push(self.isolated_vertices_removed, "isolated vertices removed");
        push(self.faces_rejected_by_ingest, "faces rejected by ingest");

        write!(f, "{}", parts.join(", "))
    }
}

/// The point and face buffers of a repaired mesh, with the account of what was done to them.
#[derive(Debug, Clone)]
pub struct Repaired {
    /// The surviving points, renumbered compactly if isolated vertices were dropped.
    pub points: Vec<Point3>,

    /// The surviving faces, re-indexed to match.
    pub faces: Vec<[u32; 3]>,

    /// What the repair changed.
    pub report: RepairReport,
}

/// Make a triangle soup manifold, returning the repaired buffers and an account of the changes.
///
/// See the module documentation for why the passes run in this order.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`
/// * `opts`: which passes to run
///
/// returns: `Result<Repaired>`, failing only if a face indexes a point which does not exist
pub fn repair_buffers(
    points: &[Point3],
    faces: &[[u32; 3]],
    opts: &RepairOpts,
) -> Result<Repaired> {
    for (i, face) in faces.iter().enumerate() {
        for &v in face.iter() {
            if v as usize >= points.len() {
                return Err(format!(
                    "Face {} references point {}, but the mesh has only {} points",
                    i,
                    v,
                    points.len()
                )
                .into());
            }
        }
    }

    let mut report = RepairReport::default();
    let mut points = points.to_vec();
    let mut faces = faces.to_vec();

    if opts.drop_degenerate {
        let before = faces.len();
        faces.retain(|face| !is_degenerate(&points, face));
        report.degenerate_removed = before - faces.len();
    }

    if opts.drop_duplicate_faces {
        let before = faces.len();
        let mut seen = HashSet::with_hasher(default());
        faces.retain(|face| {
            // Sorting the triple collapses both windings of the same triangle onto one key, which
            // is what catches the mirrored-duplicate case rather than only the identical one.
            let mut key = *face;
            key.sort_unstable();
            seen.insert(key)
        });
        report.duplicate_faces_removed = before - faces.len();
    }

    if opts.resolve_nonmanifold_edges {
        resolve_nonmanifold_edges(&points, &mut faces, &mut report);
    }

    if opts.orient_consistently {
        orient_consistently(&mut faces, &mut report);
        resolve_directed_collisions(&points, &mut faces, &mut report);
    }

    // After every pass which can drop a face, since dropping one can split a fan and create a
    // bowtie which was not in the input.
    if opts.split_bowtie_vertices {
        split_bowtie_vertices(&mut points, &mut faces, &mut report);
    }

    if opts.drop_isolated_vertices {
        drop_isolated_vertices(&mut points, &mut faces, &mut report);
    }

    Ok(Repaired {
        points,
        faces,
        report,
    })
}

/// Whether a triangle repeats a vertex or encloses no area.
fn is_degenerate(points: &[Point3], face: &[u32; 3]) -> bool {
    if face[0] == face[1] || face[1] == face[2] || face[2] == face[0] {
        return true;
    }

    let p = [
        points[face[0] as usize],
        points[face[1] as usize],
        points[face[2] as usize],
    ];

    compute_face_normal(&p).is_none()
}

/// A scale-invariant triangle quality in `(0, 1]`, peaking at 1 for an equilateral triangle.
///
/// Used to choose which faces to keep at a non-manifold edge. A sliver is both the more likely
/// artifact and the less useful surface, so it is the one to discard.
fn triangle_quality(points: &[Point3], face: &[u32; 3]) -> f64 {
    let p = [
        points[face[0] as usize],
        points[face[1] as usize],
        points[face[2] as usize],
    ];

    let e = [p[1] - p[0], p[2] - p[1], p[0] - p[2]];
    let sum_sq: f64 = e.iter().map(|v| v.norm_squared()).sum();
    if sum_sq <= 0.0 {
        return 0.0;
    }

    let area = e[0].cross(&e[1]).norm() * 0.5;
    (4.0 * 3.0_f64.sqrt() * area) / sum_sq
}

/// Build the map from each undirected edge to the faces which use it.
fn edge_to_faces(faces: &[[u32; 3]]) -> HashMap<[u32; 2], Vec<u32>> {
    let mut map: HashMap<[u32; 2], Vec<u32>> = HashMap::with_hasher(default());

    for (i, face) in faces.iter().enumerate() {
        for k in 0..3 {
            let key = edge_key(&[face[k], face[(k + 1) % 3]]);
            map.entry(key).or_default().push(i as u32);
        }
    }

    map
}

/// Drop faces until no edge is shared by more than two of them.
///
/// At each offending edge the two highest-quality faces are kept and the rest marked. Because
/// dropping a face only ever lowers the face count of the edges it touched, marking every edge
/// independently and then deleting in one pass cannot leave an edge still non-manifold, so this
/// does not need to iterate.
fn resolve_nonmanifold_edges(
    points: &[Point3],
    faces: &mut Vec<[u32; 3]>,
    report: &mut RepairReport,
) {
    let map = edge_to_faces(faces);
    let mut drop = IndexMask::new(faces.len(), false);
    let mut nonmanifold = 0usize;

    for (key, users) in map.iter() {
        if users.len() <= 2 {
            continue;
        }
        nonmanifold += 1;

        // Rank by quality, breaking ties by face index so the choice is deterministic.
        let mut ranked = users.clone();
        ranked.sort_by(|&a, &b| {
            let qa = triangle_quality(points, &faces[a as usize]);
            let qb = triangle_quality(points, &faces[b as usize]);
            qb.partial_cmp(&qa)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.cmp(&b))
        });

        // Keep the best face, then the best remaining one which walks the edge the *other* way.
        // A pair which agrees on direction cannot be oriented consistently, so choosing the second
        // on quality alone would only defer the problem to `resolve_directed_collisions`, which
        // would then have to drop a face this pass had deliberately kept.
        let first = ranked[0];
        let first_forward = walks_forward(&faces[first as usize], key);

        let second = ranked
            .iter()
            .skip(1)
            .copied()
            .find(|&f| walks_forward(&faces[f as usize], key) != first_forward)
            .unwrap_or(ranked[1]);

        for &f in ranked.iter() {
            if f != first && f != second {
                drop.set(f as usize, true);
            }
        }
    }

    report.nonmanifold_edges = nonmanifold;
    report.faces_dropped_at_nonmanifold = drop.count_true();

    if report.faces_dropped_at_nonmanifold > 0 {
        let mut index = 0usize;
        faces.retain(|_| {
            let keep = !drop.get(index);
            index += 1;
            keep
        });
    }
}

/// Whether a face traverses an undirected edge in the same direction the key is written, that is
/// from the lower vertex index to the higher.
fn walks_forward(face: &[u32; 3], key: &[u32; 2]) -> bool {
    (0..3).any(|k| face[k] == key[0] && face[(k + 1) % 3] == key[1])
}

/// Drop faces until no directed edge is used by more than one face.
///
/// In an orientable manifold each directed edge belongs to exactly one face, so two faces walking
/// an edge the same way is precisely what a half-edge structure cannot represent. Orientation
/// resolves this wherever a consistent assignment exists; whatever is left is a genuine fold or a
/// non-orientable region, and the only repair is to sacrifice a face.
///
/// The lower-quality face goes, and because dropping a face can only reduce the use count of the
/// edges it touched, a single pass is enough.
fn resolve_directed_collisions(
    points: &[Point3],
    faces: &mut Vec<[u32; 3]>,
    report: &mut RepairReport,
) {
    let mut directed: HashMap<[u32; 2], Vec<u32>> = HashMap::with_hasher(default());
    for (i, face) in faces.iter().enumerate() {
        for k in 0..3 {
            directed
                .entry([face[k], face[(k + 1) % 3]])
                .or_default()
                .push(i as u32);
        }
    }

    let mut drop = IndexMask::new(faces.len(), false);

    for users in directed.values() {
        if users.len() <= 1 {
            continue;
        }

        let mut ranked = users.clone();
        ranked.sort_by(|&a, &b| {
            let qa = triangle_quality(points, &faces[a as usize]);
            let qb = triangle_quality(points, &faces[b as usize]);
            qb.partial_cmp(&qa)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.cmp(&b))
        });

        for &f in ranked.iter().skip(1) {
            drop.set(f as usize, true);
        }
    }

    let dropped = drop.count_true();
    if dropped == 0 {
        return;
    }

    let mut index = 0usize;
    faces.retain(|_| {
        let keep = !drop.get(index);
        index += 1;
        keep
    });

    report.faces_dropped_for_orientation = dropped;
}

/// Split a vertex whose incident faces form more than one fan into one vertex per fan.
///
/// A bowtie is two sheets of surface meeting at a single point. The half-edge structure cannot
/// represent it, because the outgoing halfedges around the vertex would form more than one cycle.
/// Duplicating the point once per fan separates the sheets without moving any geometry, so no
/// measurement is altered; only the sharing of a coincident point is given up.
fn split_bowtie_vertices(
    points: &mut Vec<Point3>,
    faces: &mut [[u32; 3]],
    report: &mut RepairReport,
) {
    // Faces incident on each vertex.
    let mut vertex_faces: Vec<Vec<u32>> = vec![Vec::new(); points.len()];
    for (i, face) in faces.iter().enumerate() {
        for &v in face.iter() {
            vertex_faces[v as usize].push(i as u32);
        }
    }

    let mut split_count = 0usize;
    let mut added = 0usize;

    for v in 0..vertex_faces.len() {
        let incident = &vertex_faces[v];
        if incident.len() < 2 {
            continue;
        }

        let fans = group_into_fans(faces, incident, v as u32);
        if fans.len() < 2 {
            continue;
        }

        split_count += 1;

        // The first fan keeps the original point; every other fan gets a copy of it.
        for fan in fans.iter().skip(1) {
            let new_index = points.len() as u32;
            points.push(points[v]);
            added += 1;

            for &f in fan.iter() {
                let face = &mut faces[f as usize];
                for slot in face.iter_mut() {
                    if *slot == v as u32 {
                        *slot = new_index;
                    }
                }
            }
        }
    }

    report.bowtie_vertices_split = split_count;
    report.points_added_by_splitting = added;
}

/// Group the faces incident on a vertex into fans, where two faces share a fan when they share an
/// edge which contains that vertex.
fn group_into_fans(faces: &[[u32; 3]], incident: &[u32], vertex: u32) -> Vec<Vec<u32>> {
    // The other two vertices of each incident face, which are the far ends of the two edges that
    // pass through `vertex`. Two faces are joined when they share one of those ends.
    let spokes: Vec<[u32; 2]> = incident
        .iter()
        .map(|&f| {
            let face = &faces[f as usize];
            let mut others = [u32::MAX; 2];
            let mut n = 0;
            for &v in face.iter() {
                if v != vertex && n < 2 {
                    others[n] = v;
                    n += 1;
                }
            }
            others
        })
        .collect();

    let mut spoke_owner: HashMap<u32, Vec<usize>> = HashMap::with_hasher(default());
    for (local, pair) in spokes.iter().enumerate() {
        for &end in pair.iter() {
            if end != u32::MAX {
                spoke_owner.entry(end).or_default().push(local);
            }
        }
    }

    let mut fan_of = vec![usize::MAX; incident.len()];
    let mut fans: Vec<Vec<u32>> = Vec::new();
    let mut stack = Vec::new();

    for start in 0..incident.len() {
        if fan_of[start] != usize::MAX {
            continue;
        }

        let fan_index = fans.len();
        let mut fan = Vec::new();

        fan_of[start] = fan_index;
        stack.push(start);

        while let Some(local) = stack.pop() {
            fan.push(incident[local]);

            for &end in spokes[local].iter() {
                let Some(neighbors) = spoke_owner.get(&end) else {
                    continue;
                };
                for &other in neighbors.iter() {
                    if fan_of[other] == usize::MAX {
                        fan_of[other] = fan_index;
                        stack.push(other);
                    }
                }
            }
        }

        fans.push(fan);
    }

    fans
}

/// Flip faces so every face in a connected component agrees on which side is out.
///
/// Two triangles sharing an edge agree when they traverse that edge in opposite directions. The
/// walk fixes the first face of each component and flips whatever disagrees with it, which decides
/// orientation only up to a whole-component flip; a component whose overall sense is wrong stays
/// wrong, and that is `Mesh3::flip_normals_in_place`'s business rather than this pass's.
fn orient_consistently(faces: &mut [[u32; 3]], report: &mut RepairReport) {
    let map = edge_to_faces(faces);
    let mut visited = IndexMask::new(faces.len(), false);
    let mut flipped = 0usize;
    let mut nonorientable = 0usize;

    let mut stack: Vec<u32> = Vec::new();

    for start in 0..faces.len() {
        if visited.get(start) {
            continue;
        }

        visited.set(start, true);
        stack.push(start as u32);
        let mut component_conflict = false;

        while let Some(current) = stack.pop() {
            let face = faces[current as usize];

            for k in 0..3 {
                let a = face[k];
                let b = face[(k + 1) % 3];
                let key = edge_key(&[a, b]);

                let Some(users) = map.get(&key) else {
                    continue;
                };

                for &other in users.iter() {
                    if other == current {
                        continue;
                    }

                    let neighbor = faces[other as usize];
                    // The neighbor agrees when it walks this edge from b to a.
                    let agrees = (0..3).any(|j| neighbor[j] == b && neighbor[(j + 1) % 3] == a);

                    if visited.get(other as usize) {
                        // Already fixed. If it disagrees now, no consistent assignment exists.
                        if !agrees {
                            component_conflict = true;
                        }
                        continue;
                    }

                    visited.set(other as usize, true);
                    if !agrees {
                        faces[other as usize].swap(1, 2);
                        flipped += 1;
                    }
                    stack.push(other);
                }
            }
        }

        if component_conflict {
            nonorientable += 1;
        }
    }

    report.faces_reoriented = flipped;
    report.nonorientable_components = nonorientable;
}

/// Drop points which no surviving face references, compacting the buffers.
fn drop_isolated_vertices(
    points: &mut Vec<Point3>,
    faces: &mut [[u32; 3]],
    report: &mut RepairReport,
) {
    let mut used = vec![false; points.len()];
    for face in faces.iter() {
        for &v in face.iter() {
            used[v as usize] = true;
        }
    }

    let removed = used.iter().filter(|u| !**u).count();
    if removed == 0 {
        return;
    }

    let mut remap = vec![u32::MAX; points.len()];
    let mut kept = Vec::with_capacity(points.len() - removed);
    for (old, is_used) in used.iter().enumerate() {
        if *is_used {
            remap[old] = kept.len() as u32;
            kept.push(points[old]);
        }
    }

    for face in faces.iter_mut() {
        for slot in face.iter_mut() {
            *slot = remap[*slot as usize];
        }
    }

    *points = kept;
    report.isolated_vertices_removed = removed;
}

impl HalfEdgeMesh {
    /// Build a half-edge mesh from a `Mesh3`, repairing whatever topology the structure cannot
    /// accept and reporting what had to change.
    ///
    /// Use this rather than `TryFrom<&Mesh3>` on measured data. The `TryFrom` path is strict and
    /// fails on the first element it cannot insert, which is the right behavior when the input is
    /// supposed to be clean and the wrong one when it came off a scanner.
    ///
    /// Attributes, `is_solid`, and the UV mapping are not carried across, because the half-edge
    /// structure has nowhere to put them and the repair may renumber or drop the elements they
    /// were indexed by.
    ///
    /// # Arguments
    ///
    /// * `mesh`: the mesh to convert
    /// * `opts`: which repairs to attempt
    ///
    /// returns: `Result<(HalfEdgeMesh, RepairReport)>`
    pub fn from_mesh_repaired(mesh: &Mesh3, opts: &RepairOpts) -> Result<(Self, RepairReport)> {
        let repaired = repair_buffers(mesh.points(), mesh.faces(), opts)?;
        Self::from_repaired(repaired)
    }

    /// Build a half-edge mesh from already-repaired buffers, carrying the report through.
    ///
    /// A face the structure still rejects is skipped and counted in
    /// `RepairReport::faces_rejected_by_ingest` rather than aborting a partial build.
    pub fn from_repaired(repaired: Repaired) -> Result<(Self, RepairReport)> {
        let Repaired {
            points,
            faces,
            mut report,
        } = repaired;

        let mut result = Self::new();
        let inner = result.as_alum_mut();
        let mut handles = Vec::with_capacity(points.len());

        for p in points.iter() {
            let handle = inner
                .add_vertex(p.coords)
                .map_err(|e| format!("Failed to add vertex: {:?}", e))?;
            handles.push(handle);
        }

        for face in faces.iter() {
            if inner
                .add_tri_face(
                    handles[face[0] as usize],
                    handles[face[1] as usize],
                    handles[face[2] as usize],
                )
                .is_err()
            {
                report.faces_rejected_by_ingest += 1;
            }
        }

        Ok((result, report))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tri(a: u32, b: u32, c: u32) -> [u32; 3] {
        [a, b, c]
    }

    /// A unit square as two triangles sharing the diagonal 1-2.
    fn square() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        let faces = vec![tri(0, 1, 2), tri(1, 3, 2)];
        (points, faces)
    }

    #[test]
    fn a_clean_mesh_is_left_alone() {
        let (points, faces) = square();
        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        assert!(out.report.is_clean(), "report was {}", out.report);
        assert_eq!(out.points, points);
        assert_eq!(out.faces, faces);
    }

    #[test]
    fn a_face_referencing_a_missing_point_is_an_error() {
        let (points, _) = square();
        let faces = vec![tri(0, 1, 9)];

        assert!(repair_buffers(&points, &faces, &RepairOpts::default()).is_err());
    }

    #[test]
    fn degenerate_faces_are_dropped() {
        let (mut points, mut faces) = square();

        // A repeated vertex, and three collinear points.
        points.push(Point3::new(2.0, 0.0, 0.0));
        points.push(Point3::new(3.0, 0.0, 0.0));
        points.push(Point3::new(4.0, 0.0, 0.0));
        faces.push(tri(0, 1, 1));
        faces.push(tri(4, 5, 6));

        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.degenerate_removed, 2);
        assert_eq!(out.faces.len(), 2);
    }

    /// Both windings of the same triangle collapse onto one key, which is the case the Stanford
    /// bunny fixture actually contains.
    #[test]
    fn duplicate_faces_are_dropped_in_either_winding() {
        let (points, mut faces) = square();
        faces.push(tri(0, 1, 2)); // identical
        faces.push(tri(0, 2, 1)); // reversed

        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.duplicate_faces_removed, 2);
        assert_eq!(out.faces.len(), 2);
    }

    /// Three faces on one edge is the classic registration-seam artifact. The two best-quality
    /// faces survive and the sliver is the one dropped.
    #[test]
    fn a_nonmanifold_edge_keeps_the_two_best_faces() {
        let mut points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(0.5, -1.0, 0.0),
        ];
        // A sliver hanging off the same edge 0-1.
        points.push(Point3::new(0.5, 0.001, 0.0));

        let faces = vec![tri(0, 1, 2), tri(1, 0, 3), tri(0, 1, 4)];

        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.nonmanifold_edges, 1);
        assert_eq!(out.report.faces_dropped_at_nonmanifold, 1);
        assert_eq!(out.faces.len(), 2);

        // Point 4 was only used by the sliver, so it goes too.
        assert_eq!(out.report.isolated_vertices_removed, 1);
        assert_eq!(out.points.len(), 4);
    }

    /// Two triangles meeting at a single vertex, which alum cannot represent. The point is
    /// duplicated so each sheet gets its own, and no geometry moves.
    #[test]
    fn a_bowtie_vertex_is_split() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0), // the shared apex
            Point3::new(-1.0, 1.0, 0.0),
            Point3::new(-1.0, -1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, -1.0, 0.0),
        ];
        let faces = vec![tri(0, 1, 2), tri(0, 4, 3)];

        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.bowtie_vertices_split, 1);
        assert_eq!(out.report.points_added_by_splitting, 1);
        assert_eq!(out.points.len(), 6);
        assert_eq!(out.faces.len(), 2, "no face should have been dropped");

        // The duplicate sits exactly on top of the original, so nothing was measured differently.
        let apexes: Vec<_> = out
            .points
            .iter()
            .filter(|p| (*p - Point3::new(0.0, 0.0, 0.0)).norm() < 1e-12)
            .collect();
        assert_eq!(apexes.len(), 2);

        // The two faces no longer share any vertex.
        let a = out.faces[0];
        let b = out.faces[1];
        assert!(a.iter().all(|v| !b.contains(v)));
    }

    /// A fan of triangles around a vertex is a single fan and must not be split.
    #[test]
    fn a_normal_fan_is_not_split() {
        let mut points = vec![Point3::new(0.0, 0.0, 0.0)];
        for i in 0..6 {
            let a = std::f64::consts::TAU * (i as f64) / 6.0;
            points.push(Point3::new(a.cos(), a.sin(), 0.0));
        }
        let faces: Vec<_> = (0..6).map(|i| tri(0, 1 + i, 1 + (i + 1) % 6)).collect();

        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.bowtie_vertices_split, 0);
        assert_eq!(out.points.len(), 7);
    }

    #[test]
    fn inconsistent_winding_is_corrected() {
        let (points, mut faces) = square();
        faces[1].swap(1, 2); // flip the second triangle

        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.faces_reoriented, 1);
        assert_eq!(out.report.nonorientable_components, 0);
        assert_eq!(out.faces, square().1);
    }

    /// Orientation is only fixed up to a whole-component flip. A component which is internally
    /// consistent but globally inside-out is left alone, because deciding which way is out is not
    /// this pass's job.
    #[test]
    fn a_wholly_reversed_component_is_left_alone() {
        let (points, faces) = square();
        let reversed: Vec<_> = faces
            .iter()
            .map(|f| {
                let mut g = *f;
                g.swap(1, 2);
                g
            })
            .collect();

        let out = repair_buffers(&points, &reversed, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.faces_reoriented, 0);
        assert_eq!(out.faces, reversed);
    }

    #[test]
    fn isolated_vertices_are_dropped_and_faces_reindexed() {
        let (mut points, faces) = square();
        points.insert(0, Point3::new(9.0, 9.0, 9.0));
        let shifted: Vec<_> = faces
            .iter()
            .map(|f| tri(f[0] + 1, f[1] + 1, f[2] + 1))
            .collect();

        let out = repair_buffers(&points, &shifted, &RepairOpts::default()).unwrap();

        assert_eq!(out.report.isolated_vertices_removed, 1);
        assert_eq!(out.points.len(), 4);
        assert_eq!(out.faces, faces);
    }

    #[test]
    fn opts_can_turn_a_pass_off() {
        let (points, mut faces) = square();
        faces.push(tri(0, 2, 1));

        let opts = RepairOpts::default().with_drop_duplicate_faces(false);
        let out = repair_buffers(&points, &faces, &opts).unwrap();

        assert_eq!(out.report.duplicate_faces_removed, 0);
    }

    #[test]
    fn none_opts_change_nothing() {
        let (points, mut faces) = square();
        faces.push(tri(0, 1, 1));

        let out = repair_buffers(&points, &faces, &RepairOpts::none()).unwrap();

        assert!(out.report.is_clean());
        assert_eq!(out.faces.len(), 3);
    }

    // ---- Ingest ----------------------------------------------------------------------------

    /// Three faces on one edge is a case the strict path genuinely refuses. The repairing path
    /// takes it.
    ///
    /// Note that alum does *not* reject every non-manifold input. A bowtie whose fans are both
    /// open is accepted and even passes alum's own `check_topology`, so the strict path is not a
    /// reliable manifoldness test, only a reliable "can this be inserted" test.
    #[test]
    fn ingest_accepts_what_try_from_rejects() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(0.5, -1.0, 0.0),
            Point3::new(0.5, 0.0, 1.0),
        ];
        let faces = vec![tri(0, 1, 2), tri(1, 0, 3), tri(0, 1, 4)];
        let mesh = Mesh3::new(points, faces, false);

        assert!(
            HalfEdgeMesh::try_from(&mesh).is_err(),
            "three faces on one edge should be refused by the strict path"
        );

        let (he, report) = HalfEdgeMesh::from_mesh_repaired(&mesh, &RepairOpts::default()).unwrap();

        assert_eq!(report.nonmanifold_edges, 1);
        assert_eq!(report.faces_rejected_by_ingest, 0);
        assert_eq!(he.clone_faces().unwrap().len(), 2);
    }

    /// Two faces walking the same edge in the same direction is a fold, which no orientation
    /// assignment can fix and which the strict path also refuses.
    #[test]
    fn a_fold_is_resolved_by_dropping_a_face() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(0.5, -1.0, 0.0),
        ];
        let faces = vec![tri(0, 1, 2), tri(0, 1, 3)];
        let mesh = Mesh3::new(points.clone(), faces.clone(), false);

        assert!(HalfEdgeMesh::try_from(&mesh).is_err());

        let out = repair_buffers(&points, &faces, &RepairOpts::default()).unwrap();

        // Orientation flips one of them, which resolves the fold without losing a face.
        assert_eq!(out.faces.len(), 2);
        assert_eq!(out.report.faces_reoriented, 1);
        assert_eq!(out.report.faces_dropped_for_orientation, 0);

        let (_, report) = HalfEdgeMesh::from_mesh_repaired(&mesh, &RepairOpts::default()).unwrap();
        assert_eq!(report.faces_rejected_by_ingest, 0);
    }

    /// The bunny fixture is considerably dirtier than the one documented defect at `mesh.rs:1131`
    /// suggests: 31 distinct duplicated vertex-triples and 47 edges shared by more than two faces.
    #[test]
    fn the_bunny_fixture_is_repaired() {
        let mesh = Mesh3::stanford_bunny_res4();

        let out = repair_buffers(mesh.points(), mesh.faces(), &RepairOpts::default()).unwrap();

        assert_eq!(out.report.duplicate_faces_removed, 40);
        assert_eq!(out.report.nonmanifold_edges, 47);
        assert_eq!(out.report.faces_dropped_at_nonmanifold, 40);
        assert_eq!(out.report.bowtie_vertices_split, 15);
        assert_eq!(out.report.isolated_vertices_removed, 5);
        assert_eq!(out.faces.len(), 868);

        // Whatever else it did, the result must be something a half-edge structure can hold.
        assert_no_nonmanifold_edges(&out.faces);
        assert_no_directed_collisions(&out.faces);
    }

    /// The round trip the whole editing layer depends on: every repaired face is accepted, and
    /// converting back gives exactly what went in.
    #[test]
    fn the_bunnies_round_trip_through_the_half_edge_structure() {
        for mesh in [
            Mesh3::stanford_bunny_res4(),
            Mesh3::stanford_bunny_res3(),
            Mesh3::stanford_bunny_res2(),
        ] {
            let repaired =
                repair_buffers(mesh.points(), mesh.faces(), &RepairOpts::default()).unwrap();
            let expected = repaired.faces.len();

            let (he, report) =
                HalfEdgeMesh::from_mesh_repaired(&mesh, &RepairOpts::default()).unwrap();

            assert_eq!(
                report.faces_rejected_by_ingest, 0,
                "every repaired face should have been accepted"
            );

            let back = Mesh3::try_from(&he).unwrap();
            assert_eq!(back.faces().len(), expected);
        }
    }

    /// An already-manifold mesh must come through untouched. A repair which quietly edits clean
    /// measured data would be worse than one which fails.
    #[test]
    fn a_clean_fixture_is_untouched() {
        let mesh = crate::tests::engine_blade();

        let out = repair_buffers(mesh.points(), mesh.faces(), &RepairOpts::default()).unwrap();

        assert!(out.report.is_clean(), "report was {}", out.report);
        assert_eq!(out.faces.len(), mesh.faces().len());
        assert_eq!(out.points.len(), mesh.points().len());
    }

    /// Splitting bowties is not optional on real data. alum tolerates a bowtie whose fans are both
    /// open, but rejects one whose fan has closed, and the bunny has those.
    #[test]
    fn the_bowtie_pass_is_load_bearing() {
        let mesh = Mesh3::stanford_bunny_res4();

        let without = RepairOpts::default().with_split_bowtie_vertices(false);
        let (_, report) = HalfEdgeMesh::from_mesh_repaired(&mesh, &without).unwrap();
        assert_eq!(
            report.faces_rejected_by_ingest, 13,
            "without the split, alum should refuse these faces"
        );

        let (_, report) = HalfEdgeMesh::from_mesh_repaired(&mesh, &RepairOpts::default()).unwrap();
        assert_eq!(report.faces_rejected_by_ingest, 0);
    }

    fn assert_no_nonmanifold_edges(faces: &[[u32; 3]]) {
        let map = edge_to_faces(faces);
        let worst = map.values().map(|v| v.len()).max().unwrap_or(0);
        assert!(worst <= 2, "an edge is still shared by {} faces", worst);
    }

    fn assert_no_directed_collisions(faces: &[[u32; 3]]) {
        let mut directed: HashMap<[u32; 2], usize> = HashMap::with_hasher(default());
        for face in faces.iter() {
            for k in 0..3 {
                *directed.entry([face[k], face[(k + 1) % 3]]).or_default() += 1;
            }
        }
        let worst = directed.values().copied().max().unwrap_or(0);
        assert!(worst <= 1, "a directed edge is still used {} times", worst);
    }

    #[test]
    fn a_clean_report_displays_as_such() {
        assert_eq!(
            RepairReport::default().to_string(),
            "mesh was already manifold, nothing changed"
        );

        let report = RepairReport {
            duplicate_faces_removed: 3,
            isolated_vertices_removed: 1,
            ..Default::default()
        };
        assert_eq!(
            report.to_string(),
            "3 duplicate faces removed, 1 isolated vertices removed"
        );
    }
}
