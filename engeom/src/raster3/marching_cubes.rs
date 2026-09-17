//! Extraction of a triangle surface from the zero set of a [`SparseGrid3`] of signed distances.
//!
//! # What the band boundary does to the result
//!
//! The extractor considers a cell only when all eight corners contain written values. At the end of
//! the band, surrounding cells have at least one corner that was never evaluated, so those cells
//! contribute no geometry and the surface stops. This behavior produces an open mesh from a scan
//! of an open surface, with a boundary near the end of the measured data. It does not close the
//! mesh with unmeasured geometry.
//!
//! The cost is a ragged edge. The boundary follows cell corners, so it steps by up to a voxel
//! instead of tracing a smooth curve, and isolated fringe cells can produce small slivers. Use the
//! patch filtering in `geom3::mesh` to remove these artifacts.
//!
//! # Why a loop is not simply fanned
//!
//! The case table returns closed loops of crossing edges. A triangle fan from the first vertex can
//! fill an isolated cell, but this approach can fail when neighboring cells are considered.
//!
//! A fan over a loop of `n` vertices draws `n - 3` chords between vertices that are not adjacent
//! around the loop. Those chords run through the inside of the cell, but their endpoints are
//! vertices on the cube's edges, and a pair of them can lie on the same face. The cell on the other
//! side of that face has its own loop through the same two vertices and can draw a chord between
//! the same pair. The two chords are different segments in space but use the same pair of mesh
//! indices, so four faces share the edge between them. Of the 256
//! configurations, 144 contain a loop with a chord positioned to do this.
//!
//! A loop longer than three is therefore filled around a new vertex at its centroid, with one
//! triangle per loop segment. Every edge of the result is then either a loop segment, shared with
//! exactly one neighboring cell, or a spoke to a centroid that belongs only to this cell and cannot
//! be referenced by another cell. The surface is manifold along every edge created by filling the
//! loops. This approach costs one extra
//! vertex and two extra triangles per long loop, and it avoids the slivers a fan produces on a loop
//! of six or seven, which also improves triangle quality.
//!
//! The centroid sits slightly off the true surface, by an amount that goes as the square of the
//! voxel size over the radius of curvature. This is the same order as the existing error from
//! linear interpolation of a curved field along an edge, so the centroid does not change the
//! asymptotic accuracy of the result.
//!
//! # How a vertex is shared
//!
//! Every vertex sits on a grid edge, and an edge is named by its lower voxel key and the axis it
//! runs along. The name is independent of the requesting cell, so the four or fewer cells around an
//! edge resolve to the same vertex and form a welded mesh. The block holding the lower voxel owns
//! the vertex, which partitions the work by block without coordination between blocks.
//!
//! # Determinism
//!
//! Blocks are processed in sorted key order and vertices are numbered by a prefix sum over that
//! order, so the same grid always produces the same buffers in the same order. See the note on
//! [`SparseGrid3::block_keys`] for why that matters.

use super::mc_tables::{CORNER_OFFSETS, EDGE_CORNERS, EDGE_LOWER, MAX_LOOPS, case_table};
use super::sparse_grid::{BLOCK_EDGE, BLOCK_VOXELS, Block3, BlockKey3, SdfVoxel, SparseGrid3};
use crate::geom3::mesh::MeshData3;
use crate::{Point3, Result};
use faer::prelude::default;
use parry3d_f64::utils::hashmap::HashMap;
use rayon::prelude::*;

/// Statistics about the cells considered and the mesh produced during extraction.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IsoSurfaceStats3 {
    /// Cells whose eight corners were all written and were therefore considered for triangulation.
    /// Uniform cells are included even though they produce no triangles.
    pub cells_visited: usize,

    /// Cells that held at least one written corner and at least one unwritten corner and were
    /// skipped. This is the size of the fringe at the edge of the band.
    pub cells_skipped_unknown: usize,

    /// Vertices in the resulting mesh.
    pub vertices: usize,

    /// Triangles in the resulting mesh.
    pub faces: usize,
}

/// Extract the zero-level surface of a signed distance grid as a triangle mesh.
///
/// Triangles are wound so that their normals point toward the positive side of the field, which
/// for a distance field built with outward normals means they face out of the solid.
///
/// The result is welded: cells that share a grid edge also share its vertex. The mesh contains no
/// unreferenced points. However, marching cubes can produce a nonmanifold vertex or edge where two
/// sheets meet only at that location. This is a property of the method. Run the repair pass in
/// `geom3::half_edge3` over the result if a manifold mesh is required.
///
/// # Arguments
///
/// * `grid`: The signed distance field to extract from. Voxels that were never written are
///   treated as unknown, and any cell touching one is left out.
///
/// returns: `Result<(MeshData3, IsoSurfaceStats3)>`
pub fn extract_isosurface(grid: &SparseGrid3<SdfVoxel>) -> Result<(MeshData3, IsoSurfaceStats3)> {
    let keys = grid.block_keys();

    let mut block_index: HashMap<BlockKey3, usize> = HashMap::with_hasher(default());
    for (i, key) in keys.iter().enumerate() {
        block_index.insert(*key, i);
    }

    // Phase one: every block computes the vertices on the edges it owns, independently.
    let per_block: Vec<BlockVertices> = keys
        .par_iter()
        .map(|key| compute_block_vertices(grid, key))
        .collect();

    // Phase two: number the vertices by a prefix sum in sorted block order.
    let mut starts: Vec<u32> = Vec::with_capacity(keys.len() + 1);
    let mut total: u32 = 0;
    for block in per_block.iter() {
        starts.push(total);
        total = total
            .checked_add(block.points.len() as u32)
            .ok_or("The extracted surface has more vertices than a u32 index can hold")?;
    }
    starts.push(total);

    // Phase three: every block emits faces for cells whose lowest corner it holds.
    let emitted: Vec<CellFaces> = keys
        .par_iter()
        .map(|key| compute_block_faces(grid, key, &block_index, &per_block, &starts))
        .collect();

    let mut points: Vec<Point3> = Vec::with_capacity(total as usize);
    for block in per_block.iter() {
        points.extend_from_slice(&block.points);
    }

    let mut faces: Vec<[u32; 3]> = Vec::new();
    let mut stats = IsoSurfaceStats3::default();
    for part in emitted.iter() {
        faces.extend_from_slice(&part.faces);
        stats.cells_visited += part.cells_visited;
        stats.cells_skipped_unknown += part.cells_skipped_unknown;
    }

    let (points, faces) = compact_unused_points(points, faces);

    stats.vertices = points.len();
    stats.faces = faces.len();

    Ok((MeshData3::new(points, faces)?, stats))
}

/// The vertices a single block owns, and where to find each of them.
///
/// A block owns two kinds of vertices. Edge vertices sit on grid edges whose lower voxel belongs to
/// the block and are shared with neighboring cells that use those edges. Centroid vertices belong
/// to a single cell rooted in this block and cannot be referenced by another cell.
struct BlockVertices {
    points: Vec<Point3>,

    /// For each owned edge, the index of its vertex within `points`, or `u32::MAX` for an edge
    /// that carries none. Indexed by `local_index * 3 + axis`.
    edges: Vec<u32>,

    /// For each cell rooted in this block, the centroid vertex of each of its loops, or `u32::MAX`
    /// where the loop was short enough to need none. Indexed by `local_index * MAX_LOOPS + loop`.
    centroids: Vec<u32>,
}

impl BlockVertices {
    fn empty() -> Self {
        Self {
            points: Vec::new(),
            edges: vec![u32::MAX; BLOCK_VOXELS * 3],
            centroids: vec![u32::MAX; BLOCK_VOXELS * MAX_LOOPS],
        }
    }
}

/// The faces one block contributed, with the cell counts that went with them.
struct CellFaces {
    faces: Vec<[u32; 3]>,
    cells_visited: usize,
    cells_skipped_unknown: usize,
}

/// A read-only view of a block and the seven neighbors reachable by stepping one block in the
/// positive direction on any axis.
///
/// A cell rooted in a block reaches one voxel past it on each axis, so those eight blocks hold
/// every value the block's cells can need. Gathering them once per block turns the eight hash
/// lookups a cell would otherwise do into eight per block.
struct BlockWindow<'a, T> {
    parts: [Option<&'a T>; 8],
}

impl<'a, T> BlockWindow<'a, T> {
    fn new(base: &BlockKey3, mut fetch: impl FnMut(&BlockKey3) -> Option<&'a T>) -> Self {
        let mut parts: [Option<&'a T>; 8] = [None; 8];
        for dx in 0..2usize {
            for dy in 0..2usize {
                for dz in 0..2usize {
                    let key = [
                        base[0] + dx as i32,
                        base[1] + dy as i32,
                        base[2] + dz as i32,
                    ];
                    parts[(dx * 2 + dy) * 2 + dz] = fetch(&key);
                }
            }
        }
        Self { parts }
    }

    /// The part holding a coordinate in `0..=BLOCK_EDGE`, with the coordinate rebased into it.
    fn locate(&self, x: usize, y: usize, z: usize) -> Option<(&'a T, usize)> {
        let (dx, lx) = (x / BLOCK_EDGE, x % BLOCK_EDGE);
        let (dy, ly) = (y / BLOCK_EDGE, y % BLOCK_EDGE);
        let (dz, lz) = (z / BLOCK_EDGE, z % BLOCK_EDGE);

        self.parts[(dx * 2 + dy) * 2 + dz]
            .map(|part| (part, Block3::<SdfVoxel>::local_index_of(lx, ly, lz)))
    }
}

/// The value at a coordinate in `0..=BLOCK_EDGE` relative to a block, if it was written.
fn sample(window: &BlockWindow<'_, Block3<SdfVoxel>>, x: usize, y: usize, z: usize) -> Option<f64> {
    let (block, local) = window.locate(x, y, z)?;
    block.at(local).value()
}

/// The corner values of the cell rooted at a local coordinate, packed into a case index, or `None`
/// if any of the eight was never written.
fn cell_case(window: &BlockWindow<'_, Block3<SdfVoxel>>, x: usize, y: usize, z: usize) -> CellCase {
    let mut case = 0u8;
    let mut known = 0;

    for (i, offset) in CORNER_OFFSETS.iter().enumerate() {
        let cx = x + offset[0] as usize;
        let cy = y + offset[1] as usize;
        let cz = z + offset[2] as usize;

        if let Some(v) = sample(window, cx, cy, cz) {
            known += 1;
            if is_inside(v) {
                case |= 1 << i;
            }
        }
    }

    if known == 8 {
        CellCase::Known(case)
    } else if known > 0 {
        CellCase::Partial
    } else {
        CellCase::Absent
    }
}

enum CellCase {
    /// All eight corners were written, giving this configuration.
    Known(u8),
    /// Some corners were written and some were not, so the cell sits on the edge of the band.
    Partial,
    /// Nothing here was ever written.
    Absent,
}

/// The world position of the lowest corner of the cell rooted at a local coordinate in a block.
fn cell_origin(
    grid: &SparseGrid3<SdfVoxel>,
    base: [i32; 3],
    x: usize,
    y: usize,
    z: usize,
) -> Point3 {
    Point3::new(
        grid.origin().x + (base[0] + x as i32) as f64 * grid.voxel_size(),
        grid.origin().y + (base[1] + y as i32) as f64 * grid.voxel_size(),
        grid.origin().z + (base[2] + z as i32) as f64 * grid.voxel_size(),
    )
}

/// The surface crossing on one cube edge of a cell, in world coordinates.
///
/// Both corners are known to be written and to straddle the surface, because the caller only asks
/// about edges the case table listed.
fn cell_edge_point(
    grid: &SparseGrid3<SdfVoxel>,
    window: &BlockWindow<'_, Block3<SdfVoxel>>,
    base: [i32; 3],
    x: usize,
    y: usize,
    z: usize,
    edge: u8,
) -> Option<Point3> {
    let (a, b) = EDGE_CORNERS[edge as usize];
    let (oa, ob) = (CORNER_OFFSETS[a as usize], CORNER_OFFSETS[b as usize]);

    let va = sample(
        window,
        x + oa[0] as usize,
        y + oa[1] as usize,
        z + oa[2] as usize,
    )?;
    let vb = sample(
        window,
        x + ob[0] as usize,
        y + ob[1] as usize,
        z + ob[2] as usize,
    )?;

    let origin = cell_origin(grid, base, x, y, z);
    let h = grid.voxel_size();
    let pa = origin + crate::Vector3::new(oa[0] as f64, oa[1] as f64, oa[2] as f64) * h;
    let pb = origin + crate::Vector3::new(ob[0] as f64, ob[1] as f64, ob[2] as f64) * h;

    Some(interpolate(&pa, va, &pb, vb))
}

fn block_base(key: &BlockKey3) -> [i32; 3] {
    [
        key[0] * BLOCK_EDGE as i32,
        key[1] * BLOCK_EDGE as i32,
        key[2] * BLOCK_EDGE as i32,
    ]
}

fn compute_block_vertices(grid: &SparseGrid3<SdfVoxel>, key: &BlockKey3) -> BlockVertices {
    let window = BlockWindow::new(key, |k| grid.block(k));
    let mut out = BlockVertices::empty();
    let base = block_base(key);
    let table = case_table();

    for x in 0..BLOCK_EDGE {
        for y in 0..BLOCK_EDGE {
            for z in 0..BLOCK_EDGE {
                let local = Block3::<SdfVoxel>::local_index_of(x, y, z);

                // The vertices on the three grid edges leaving this sample in the positive
                // direction, which are the ones this block owns.
                if let Some(here) = sample(&window, x, y, z) {
                    let coords = [x, y, z];

                    for axis in 0..3 {
                        let mut other = coords;
                        other[axis] += 1;

                        let Some(there) = sample(&window, other[0], other[1], other[2]) else {
                            continue;
                        };

                        if is_inside(here) == is_inside(there) {
                            continue;
                        }

                        let low = cell_origin(grid, base, x, y, z);
                        let mut high = low;
                        high[axis] += grid.voxel_size();

                        out.edges[local * 3 + axis] = out.points.len() as u32;
                        out.points.push(interpolate(&low, here, &high, there));
                    }
                }

                // The centroid of each loop long enough to need one, for the cell rooted here.
                let CellCase::Known(case) = cell_case(&window, x, y, z) else {
                    continue;
                };

                for l in 0..table.loop_count(case) {
                    let cycle = table.loop_edges(case, l);
                    if cycle.len() <= 3 {
                        continue;
                    }

                    let mut sum = crate::Vector3::zeros();
                    let mut count = 0.0;
                    for edge in cycle {
                        if let Some(p) = cell_edge_point(grid, &window, base, x, y, z, *edge) {
                            sum += p.coords;
                            count += 1.0;
                        }
                    }

                    debug_assert_eq!(
                        count as usize,
                        cycle.len(),
                        "a known cell was missing a crossing the table listed"
                    );

                    if count > 0.0 {
                        out.centroids[local * MAX_LOOPS + l] = out.points.len() as u32;
                        out.points.push(Point3::from(sum / count));
                    }
                }
            }
        }
    }

    out
}

fn compute_block_faces(
    grid: &SparseGrid3<SdfVoxel>,
    key: &BlockKey3,
    block_index: &HashMap<BlockKey3, usize>,
    per_block: &[BlockVertices],
    starts: &[u32],
) -> CellFaces {
    let values = BlockWindow::new(key, |k| grid.block(k));
    let vertices = BlockWindow::new(key, |k| block_index.get(k).map(|i| &per_block[*i]));

    // The vertex numbering offset of each block in the window, in the same slot order.
    let mut offsets = [0u32; 8];
    for dx in 0..2usize {
        for dy in 0..2usize {
            for dz in 0..2usize {
                let k = [key[0] + dx as i32, key[1] + dy as i32, key[2] + dz as i32];
                if let Some(i) = block_index.get(&k) {
                    offsets[(dx * 2 + dy) * 2 + dz] = starts[*i];
                }
            }
        }
    }

    // This block's own vertices, which is where every centroid of a cell rooted here lives.
    let own_index = block_index[key];
    let own = &per_block[own_index];
    let own_start = starts[own_index];

    let table = case_table();
    let mut out = CellFaces {
        faces: Vec::new(),
        cells_visited: 0,
        cells_skipped_unknown: 0,
    };

    // Scratch space for one loop's resolved vertex indices, reused across cells.
    let mut ring: Vec<u32> = Vec::with_capacity(12);

    for x in 0..BLOCK_EDGE {
        for y in 0..BLOCK_EDGE {
            for z in 0..BLOCK_EDGE {
                let case = match cell_case(&values, x, y, z) {
                    CellCase::Known(case) => case,
                    CellCase::Partial => {
                        out.cells_skipped_unknown += 1;
                        continue;
                    }
                    CellCase::Absent => continue,
                };

                out.cells_visited += 1;
                let local = Block3::<SdfVoxel>::local_index_of(x, y, z);

                for l in 0..table.loop_count(case) {
                    let cycle = table.loop_edges(case, l);

                    ring.clear();
                    let mut resolved = true;

                    for edge in cycle {
                        let (offset, axis) = EDGE_LOWER[*edge as usize];
                        let ex = x + offset[0] as usize;
                        let ey = y + offset[1] as usize;
                        let ez = z + offset[2] as usize;

                        let Some((owner, owner_local)) = vertices.locate(ex, ey, ez) else {
                            resolved = false;
                            break;
                        };

                        let index = owner.edges[owner_local * 3 + axis];
                        if index == u32::MAX {
                            resolved = false;
                            break;
                        }

                        let (dx, dy, dz) = (ex / BLOCK_EDGE, ey / BLOCK_EDGE, ez / BLOCK_EDGE);
                        ring.push(offsets[(dx * 2 + dy) * 2 + dz] + index);
                    }

                    // Every crossed edge of a fully known cell carries a vertex, so resolution
                    // should always succeed. If an implementation defect violates this invariant,
                    // omit the loop to avoid terminating the caller's process.
                    debug_assert!(resolved, "a fully known cell referenced a missing vertex");
                    if !resolved {
                        continue;
                    }

                    if ring.len() == 3 {
                        out.faces.push([ring[0], ring[1], ring[2]]);
                        continue;
                    }

                    let centroid = own.centroids[local * MAX_LOOPS + l];
                    debug_assert_ne!(centroid, u32::MAX, "a long loop has no centroid");
                    if centroid == u32::MAX {
                        continue;
                    }
                    let centroid = own_start + centroid;

                    for k in 0..ring.len() {
                        out.faces
                            .push([centroid, ring[k], ring[(k + 1) % ring.len()]]);
                    }
                }
            }
        }
    }

    out
}

/// Whether a field value is on the negative side of the surface.
///
/// A value of zero counts as outside. Either classification is valid, but all cells must use the
/// same one. Otherwise, two cells can disagree about a corner on the surface and tear the mesh
/// between them.
fn is_inside(value: f64) -> bool {
    value < 0.0
}

/// The point on a grid edge where the field passes through zero, by linear interpolation.
///
/// The caller has already established that the two values straddle zero, so their difference
/// cannot vanish.
fn interpolate(low: &Point3, low_value: f64, high: &Point3, high_value: f64) -> Point3 {
    let t = low_value / (low_value - high_value);
    low + (high - low) * t
}

/// Drop points that no face refers to and renumber the faces to match.
///
/// An edge can carry a vertex while every surrounding cell is missing a corner, leaving an unused
/// point in the buffer. These points occur at the band fringe. Removing them keeps the reported
/// mesh size accurate and prevents loose points from affecting measurements.
fn compact_unused_points(
    points: Vec<Point3>,
    mut faces: Vec<[u32; 3]>,
) -> (Vec<Point3>, Vec<[u32; 3]>) {
    let mut used = vec![false; points.len()];
    for face in faces.iter() {
        for i in face {
            used[*i as usize] = true;
        }
    }

    if used.iter().all(|u| *u) {
        return (points, faces);
    }

    let mut remap = vec![u32::MAX; points.len()];
    let mut kept = Vec::with_capacity(points.len());
    for (i, point) in points.into_iter().enumerate() {
        if used[i] {
            remap[i] = kept.len() as u32;
            kept.push(point);
        }
    }

    for face in faces.iter_mut() {
        for i in face.iter_mut() {
            *i = remap[*i as usize];
        }
    }

    (kept, faces)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// A deterministic generator, so that a failing seed can be rerun.
    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }

        /// A value in `-1.0..1.0`.
        fn next_signed(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
        }
    }

    /// Build a grid by evaluating a field at every sample in a key range, writing only where it
    /// returns a value.
    fn build_grid(
        voxel_size: f64,
        lo: i32,
        hi: i32,
        field: impl Fn(&Point3) -> Option<f64>,
    ) -> SparseGrid3<SdfVoxel> {
        let mut grid =
            SparseGrid3::new(voxel_size, Point3::origin()).expect("grid creation failed");

        for x in lo..=hi {
            for y in lo..=hi {
                for z in lo..=hi {
                    let key = [x, y, z];
                    let p = grid.corner_position(&key);
                    if let Some(v) = field(&p) {
                        *grid.get_mut_or_insert(&key) = SdfVoxel::new(v, 1.0);
                    }
                }
            }
        }

        grid
    }

    /// How many times each directed edge is used by the faces.
    fn directed_edge_counts(faces: &[[u32; 3]]) -> HashMap<(u32, u32), usize> {
        let mut counts = HashMap::new();
        for face in faces {
            for k in 0..3 {
                *counts.entry((face[k], face[(k + 1) % 3])).or_insert(0) += 1;
            }
        }
        counts
    }

    /// How many faces touch each undirected edge.
    fn undirected_edge_counts(faces: &[[u32; 3]]) -> HashMap<(u32, u32), usize> {
        let mut counts = HashMap::new();
        for face in faces {
            for k in 0..3 {
                let (a, b) = (face[k], face[(k + 1) % 3]);
                *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        counts
    }

    /// Assert that a mesh is closed, manifold along its edges, and consistently wound.
    fn assert_closed_and_oriented(faces: &[[u32; 3]], context: &str) {
        let directed = directed_edge_counts(faces);

        for (&(a, b), &count) in directed.iter() {
            assert_eq!(
                count, 1,
                "{context}: directed edge {a}->{b} was used {count} times"
            );
            assert_eq!(
                directed.get(&(b, a)).copied().unwrap_or(0),
                1,
                "{context}: edge {a}-{b} has no opposing face, so the surface is open or twisted"
            );
        }
    }

    // ===============================================================================================
    // Degenerate inputs
    // ===============================================================================================

    #[test]
    fn an_empty_grid_extracts_an_empty_mesh() {
        let grid = SparseGrid3::<SdfVoxel>::new(1.0, Point3::origin()).expect("grid failed");
        let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

        assert_eq!(mesh.points().len(), 0);
        assert_eq!(mesh.faces().len(), 0);
        assert_eq!(stats, IsoSurfaceStats3::default());
    }

    #[test]
    fn a_grid_with_no_sign_change_extracts_nothing() {
        let grid = build_grid(1.0, 0, 6, |_| Some(1.0));
        let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

        assert_eq!(mesh.points().len(), 0);
        assert_eq!(mesh.faces().len(), 0);

        // Samples were written for keys 0 through 6, so the cells with all eight corners written
        // are those rooted at 0 through 5, and the ones rooted at 6 reach a key that was never
        // evaluated. That fringe is the band boundary in miniature.
        assert_eq!(stats.cells_visited, 6 * 6 * 6);
        assert_eq!(stats.cells_skipped_unknown, 7 * 7 * 7 - 6 * 6 * 6);
    }

    /// A single negative corner in an otherwise positive cell gives one triangle, and it must face
    /// away from that corner. This ties the table's orientation to real geometry, where the earlier
    /// table test only checked the edge numbering.
    #[test]
    fn one_negative_corner_gives_an_outward_triangle() {
        let mut grid = SparseGrid3::<SdfVoxel>::new(1.0, Point3::origin()).expect("grid failed");

        for x in 0..2 {
            for y in 0..2 {
                for z in 0..2 {
                    let inside = x == 0 && y == 0 && z == 0;
                    *grid.get_mut_or_insert(&[x, y, z]) =
                        SdfVoxel::new(if inside { -1.0 } else { 1.0 }, 1.0);
                }
            }
        }

        let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

        assert_eq!(stats.cells_visited, 1);
        assert_eq!(mesh.faces().len(), 1);
        assert_eq!(mesh.points().len(), 3);

        let f = mesh.faces()[0];
        let p: Vec<Point3> = f.iter().map(|&i| mesh.points()[i as usize]).collect();
        let normal = (p[1] - p[0]).cross(&(p[2] - p[0]));

        // The negative corner is at the origin, so the surface must face away from it.
        let centroid = Point3::from((p[0].coords + p[1].coords + p[2].coords) / 3.0);
        assert!(
            normal.dot(&centroid.coords) > 0.0,
            "the triangle faces the negative corner"
        );

        // Each vertex sits at the midpoint of an edge leaving the origin.
        for point in p.iter() {
            let on_axis = point.coords.iter().filter(|c| **c > 0.0).count();
            assert_eq!(on_axis, 1, "vertex {point:?} is not on an axis edge");
            assert!((point.coords.sum() - 0.5).abs() < 1e-12);
        }
    }

    /// A loop longer than three gets a centroid vertex and one triangle per loop segment. This
    /// verifies the construction described in the module documentation, which avoids the four
    /// faces per edge that a fan can produce.
    #[test]
    fn a_long_loop_is_filled_around_a_centroid() {
        // Corners 0, 2, 3, 4 and 5 negative is a configuration whose single loop runs over seven
        // edges, which is the longest any case produces.
        let case: u8 = 0b0011_1101;

        let mut grid = SparseGrid3::<SdfVoxel>::new(1.0, Point3::origin()).expect("grid failed");
        for (i, offset) in CORNER_OFFSETS.iter().enumerate() {
            let inside = (case >> i) & 1 == 1;
            *grid.get_mut_or_insert(&[offset[0], offset[1], offset[2]]) =
                SdfVoxel::new(if inside { -1.0 } else { 1.0 }, 1.0);
        }

        let table = case_table();
        assert_eq!(
            table.loop_count(case),
            1,
            "the chosen case is not a single loop"
        );
        let length = table.loop_edges(case, 0).len();
        assert_eq!(length, 7, "the chosen case does not have a loop of seven");

        let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

        assert_eq!(stats.cells_visited, 1);
        assert_eq!(
            mesh.faces().len(),
            length,
            "a loop of {length} should give one triangle per segment"
        );
        assert_eq!(
            mesh.points().len(),
            length + 1,
            "a loop of {length} should add exactly one centroid vertex"
        );

        // The centroid is the one vertex not sitting on a cube edge, and every triangle uses it.
        let mut counts = vec![0usize; mesh.points().len()];
        for face in mesh.faces() {
            for i in face {
                counts[*i as usize] += 1;
            }
        }
        let apex = counts
            .iter()
            .position(|&c| c == length)
            .expect("no vertex is shared by every triangle");

        // Each loop vertex is used by two triangles, the apex by all of them.
        for (i, &c) in counts.iter().enumerate() {
            if i != apex {
                assert_eq!(c, 2, "loop vertex {i} is used {c} times");
            }
        }

        // Every spoke to the apex is traversed once each way, which is what makes the fill
        // manifold.
        let directed = directed_edge_counts(mesh.faces());
        for (&(a, b), &n) in directed.iter() {
            assert_eq!(n, 1, "directed edge {a}->{b} was used {n} times");
        }
    }

    // ===============================================================================================
    // The structural guarantee
    // ===============================================================================================

    /// The central claim: over an arbitrary field whose outer shell is positive, the extracted
    /// surface is closed, manifold along every edge, and consistently wound.
    ///
    /// Random corner values exercise all 256 cases many times, including the ambiguous cases. A
    /// mistake in the generated table or in the way vertices are shared
    /// between blocks would show up here as an unmatched edge. The region spans several blocks on
    /// purpose because sharing a vertex within a block and sharing one across a block boundary are
    /// different code paths.
    #[test]
    fn a_random_field_extracts_a_closed_oriented_surface() {
        const SEEDS: u64 = 400;
        const LO: i32 = -9;
        const HI: i32 = 9;

        let mut total_faces = 0usize;

        for seed in 1..=SEEDS {
            let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));

            // Draw every value up front so the field is a plain lookup, since `build_grid` visits
            // the samples in a fixed order.
            let mut values = HashMap::new();
            for x in LO..=HI {
                for y in LO..=HI {
                    for z in LO..=HI {
                        let shell = x == LO || x == HI || y == LO || y == HI || z == LO || z == HI;
                        let v = if shell { 1.0 } else { rng.next_signed() };
                        values.insert((x, y, z), v);
                    }
                }
            }

            let grid = SparseGrid3::new(1.0, Point3::origin()).expect("grid failed");
            let mut grid = grid;
            for (&(x, y, z), &v) in values.iter() {
                *grid.get_mut_or_insert(&[x, y, z]) = SdfVoxel::new(v, 1.0);
            }

            let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

            // Cells rooted at the last populated key reach one key further and are skipped, which
            // leaves the fringe below. Everything inside it was triangulated.
            let span = (HI - LO) as usize;
            assert_eq!(
                stats.cells_visited,
                span * span * span,
                "seed {seed}: not every fully populated cell was visited"
            );

            assert_closed_and_oriented(mesh.faces(), &format!("seed {seed}"));
            total_faces += mesh.faces().len();

            // Nothing may be left over: every point must be used, and no point may repeat.
            let mut used = vec![false; mesh.points().len()];
            for face in mesh.faces() {
                for i in face {
                    used[*i as usize] = true;
                }
            }
            assert!(
                used.iter().all(|u| *u),
                "seed {seed}: the mesh carries unreferenced points"
            );
        }

        assert!(
            total_faces > SEEDS as usize * 100,
            "the random fields produced too little surface to have proven much"
        );
    }

    // ===============================================================================================
    // A known surface
    // ===============================================================================================

    const SPHERE_R: f64 = 5.0;
    const SPHERE_H: f64 = 0.5;

    fn sphere_center() -> Point3 {
        // Offset from the grid so that no sample lands exactly on the surface.
        Point3::new(0.13, -0.07, 0.21)
    }

    fn sphere_grid() -> SparseGrid3<SdfVoxel> {
        let c = sphere_center();
        let band = 2.5 * SPHERE_H;

        build_grid(SPHERE_H, -16, 16, |p| {
            let d = (p - c).norm() - SPHERE_R;
            if d.abs() <= band { Some(d) } else { None }
        })
    }

    /// An exact sphere sample must produce a closed mesh with the expected topology and size. These
    /// properties exercise the sign convention, interpolation, and welding.
    #[test]
    fn a_sampled_sphere_extracts_correctly() {
        let grid = sphere_grid();
        let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

        assert!(mesh.faces().len() > 1000, "the sphere came out too coarse");
        assert_closed_and_oriented(mesh.faces(), "sphere");

        // A closed surface of genus zero.
        let v = mesh.points().len() as i64;
        let f = mesh.faces().len() as i64;
        let e = undirected_edge_counts(mesh.faces()).len() as i64;
        assert_eq!(v - e + f, 2, "the sphere is not topologically a sphere");

        // The band was made wide enough that no crossing cell could lose a corner.
        assert!(stats.cells_visited > 0);
        assert_eq!(stats.vertices, v as usize);
        assert_eq!(stats.faces, f as usize);

        // Every vertex must sit on the sphere, to the accuracy linear interpolation of a curved
        // field allows. The error is second order in the voxel size over the radius.
        let c = sphere_center();
        let mut worst = 0.0f64;
        for p in mesh.points() {
            worst = worst.max(((p - c).norm() - SPHERE_R).abs());
        }
        // Measured on this fixture at 0.0355 voxels worst and 0.0137 mean, for a sphere ten
        // voxels in radius. The worst case is a loop centroid, which sits inside the surface by
        // the sagitta of its loop; the edge vertices alone do better. There is nothing random
        // here, so the bound is slightly above the measurement.
        assert!(
            worst < 0.05 * SPHERE_H,
            "worst radial error was {worst}, which is {:.4} voxels",
            worst / SPHERE_H
        );

        // Faces must point away from the middle, which is where the positive side is.
        for face in mesh.faces() {
            let p: Vec<Point3> = face.iter().map(|&i| mesh.points()[i as usize]).collect();
            let normal = (p[1] - p[0]).cross(&(p[2] - p[0]));
            let centroid = (p[0].coords + p[1].coords + p[2].coords) / 3.0;
            assert!(
                normal.dot(&(centroid - c.coords)) > 0.0,
                "a face on the sphere points inward"
            );
        }
    }

    /// Flipping the field's sign must reverse the surface winding without changing its size.
    #[test]
    fn negating_the_field_reverses_the_winding() {
        let c = sphere_center();
        let band = 2.5 * SPHERE_H;
        let flipped = build_grid(SPHERE_H, -16, 16, |p| {
            let d = (p - c).norm() - SPHERE_R;
            if d.abs() <= band { Some(-d) } else { None }
        });

        let (normal_mesh, _) = extract_isosurface(&sphere_grid()).expect("extraction failed");
        let (flipped_mesh, _) = extract_isosurface(&flipped).expect("extraction failed");

        assert_eq!(normal_mesh.points().len(), flipped_mesh.points().len());
        assert_eq!(normal_mesh.faces().len(), flipped_mesh.faces().len());
        assert_closed_and_oriented(flipped_mesh.faces(), "flipped sphere");

        for face in flipped_mesh.faces() {
            let p: Vec<Point3> = face
                .iter()
                .map(|&i| flipped_mesh.points()[i as usize])
                .collect();
            let normal = (p[1] - p[0]).cross(&(p[2] - p[0]));
            let centroid = (p[0].coords + p[1].coords + p[2].coords) / 3.0;
            assert!(
                normal.dot(&(centroid - c.coords)) < 0.0,
                "a face on the negated sphere still points outward"
            );
        }
    }

    // ===============================================================================================
    // The edge of the band
    // ===============================================================================================

    /// The surface must stop with the band and leave an open boundary without adding unmeasured
    /// geometry.
    #[test]
    fn an_unknown_region_leaves_the_surface_open() {
        let c = sphere_center();
        let band = 2.5 * SPHERE_H;

        // Only the upper half of the shell is populated.
        let grid = build_grid(SPHERE_H, -16, 16, |p| {
            let d = (p - c).norm() - SPHERE_R;
            if d.abs() <= band && p.z >= c.z {
                Some(d)
            } else {
                None
            }
        });

        let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

        assert!(
            mesh.faces().len() > 100,
            "the half sphere came out too coarse"
        );
        assert!(
            stats.cells_skipped_unknown > 0,
            "cutting the band should have left a fringe of skipped cells"
        );

        let undirected = undirected_edge_counts(mesh.faces());
        let boundary = undirected.values().filter(|&&n| n == 1).count();

        assert!(boundary > 0, "the cut surface has no open boundary");
        assert!(
            undirected.values().all(|&n| n <= 2),
            "an edge was shared by more than two faces"
        );

        // Every face must lie in the populated half, proving that no triangle was built from a
        // corner that was never evaluated.
        for face in mesh.faces() {
            for &i in face {
                let p = mesh.points()[i as usize];
                assert!(
                    p.z >= c.z - SPHERE_H - 1e-9,
                    "a vertex at {p:?} came from outside the populated region"
                );
            }
        }
    }

    // ===============================================================================================
    // Determinism
    // ===============================================================================================

    /// Two extractions of the same grid must produce byte-identical buffers. Hash-map block order
    /// is unstable, and vertex numbering depends on that order unless the extractor sorts it.
    #[test]
    fn extraction_is_reproducible() {
        let grid = sphere_grid();

        let (first, first_stats) = extract_isosurface(&grid).expect("extraction failed");
        let (second, second_stats) = extract_isosurface(&grid).expect("extraction failed");

        assert_eq!(first_stats, second_stats);
        assert_eq!(first.faces(), second.faces());
        assert_eq!(first.points().len(), second.points().len());
        for (a, b) in first.points().iter().zip(second.points().iter()) {
            assert_eq!(a, b);
        }
    }

    /// The same sampled surface must have the same structure at every position relative to the
    /// block grid. Duplicating a vertex across a block boundary would violate this property.
    #[test]
    fn shifting_the_surface_across_block_boundaries_changes_nothing_structural() {
        let band = 2.5 * SPHERE_H;

        let counts: Vec<(usize, usize)> = [0.0, 1.0, 2.0, 3.0]
            .iter()
            .map(|shift| {
                let c = sphere_center() + crate::Vector3::new(*shift, *shift, *shift) * SPHERE_H;
                let grid = build_grid(SPHERE_H, -24, 24, |p| {
                    let d = (p - c).norm() - SPHERE_R;
                    if d.abs() <= band { Some(d) } else { None }
                });

                let (mesh, _) = extract_isosurface(&grid).expect("extraction failed");
                assert_closed_and_oriented(mesh.faces(), "shifted sphere");
                (mesh.points().len(), mesh.faces().len())
            })
            .collect();

        // Shifting by whole voxels moves the sphere through the block lattice without changing how
        // it is sampled, so the mesh size must not move either.
        for window in counts.windows(2) {
            assert_eq!(
                window[0], window[1],
                "the mesh size changed when the surface moved across block boundaries: {counts:?}"
            );
        }
    }

    /// The edge numbering the extractor uses to name a vertex has to match the corners the table
    /// expects, or triangles would be built from vertices on the wrong edges.
    #[test]
    fn edge_lower_matches_the_corners_the_table_uses() {
        for edge in 0..12usize {
            let (a, b) = EDGE_CORNERS[edge];
            let (offset, axis) = EDGE_LOWER[edge];

            let oa = CORNER_OFFSETS[a as usize];
            let ob = CORNER_OFFSETS[b as usize];
            let lower = if oa[axis] < ob[axis] { oa } else { ob };
            let upper = if oa[axis] < ob[axis] { ob } else { oa };

            assert_eq!(offset, lower);
            assert_eq!(upper[axis] - lower[axis], 1);
            for d in 0..3 {
                if d != axis {
                    assert_eq!(lower[d], upper[d], "edge {edge} is not axis aligned");
                }
            }
        }
    }
}
