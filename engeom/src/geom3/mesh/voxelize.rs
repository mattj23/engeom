//! Reducing a mesh to a point cloud by finding which cells of a regular voxel grid its surface
//! passes through.
//!
//! Each face is rasterized onto the grid by flooring the triangle's bounding box to an inclusive
//! range of integer cell keys, then keeping each candidate cell if it genuinely intersects the
//! triangle. This is a full separating-axis test, so the result should be exact, and a triangle
//! that only clips the corner of a cell is still included.
//!
//! Every occupied cell then emits one point, which lies **on** the surface: the cell center is
//! projected onto the triangles passing through that cell and the nearest projection is kept, along
//! with that triangle's normal. No BVH is involved, since the projection uses only the three
//! vertices of a triangle already in hand, which is why this is available on `MeshData3` as well as
//! `Mesh3`.
//!
//! # Why the cell center itself is not the output
//!
//! Outputting the cell centers directly was the obvious first step.  But it turns out that the
//! projection only ends up being 13% to 15% of the total runtime, and the voxel sampling is
//! inherently fast to begin with.  Using the centers rather than the projection means that the
//! points lie on the grid lattice and so can be off the surface by up to `v√3 / 2` and about `0.4v`
//! RMS.  For the small additional cost, projecting the cells back down to the mesh seemed like a
//! good tradeoff.
//!
//! # Which point a cell emits when several triangles pass through it
//!
//! When several triangles pass through the same voxel, the nearest projection wins with ties being
//! broken towards the lower face index.  I considered averaging the points, but that runs into the
//! issue of having the points come up off the surface.
//!
//! Interestingly, equidistant projections are common, which I hadn't originally anticipated. If the
//! cell center is in the voronoi region of an edge it will project onto it, and if the edge is
//! shared by two faces as most edges are, the projection will be exactly identical (after all, the
//! faces may be different but they index to the exact same two vertices for their shared edge, only
//! the order is different).
//!
//! # Cost
//!
//! Work scales with the number of (face, cell) pairs, which is roughly the face count while faces
//! are smaller than a cell and roughly `surface_area / v^2` once they are larger. Choosing a voxel
//! size far below the mesh resolution is therefore expensive without producing more information
//! than the mesh contains.

use super::{Mesh3, MeshData3};
use crate::common::IndexMask;
use crate::geom3::point_cloud::PointCloud3;
use crate::{Point3, Result, SurfacePoint3, UnitVec3};
use faer::prelude::default;
use parry3d_f64::bounding_volume::Aabb;
use parry3d_f64::query::PointQuery;
use parry3d_f64::query::details::intersection_test_aabb_triangle;
use parry3d_f64::shape::Triangle;
use parry3d_f64::utils::hashmap::{Entry, HashMap};
use rayon::prelude::*;

/// How many faces one rayon task rasterizes before its partial result is merged into the whole.
///
/// The merge is order-independent, so this only trades task overhead against the size of the
/// partial maps and has no effect on the output.
const FACE_CHUNK: usize = 4096;

// ===============================================================================================
// Public entry points
// ===============================================================================================

impl MeshData3 {
    /// Reduce the surface to one point per occupied cell of a regular voxel grid, placing each
    /// point on the surface itself.
    ///
    /// For each occupied cell, the cell center is projected onto every triangle passing through
    /// that cell and the nearest projection is kept, along with that triangle's normal. The result
    /// is a cloud whose points lie on the mesh to within floating point.
    ///
    /// Degenerate faces, meaning those with no computable normal, are skipped. They have zero area
    /// and contribute no surface, so a cell touched by nothing else is simply not occupied.
    ///
    /// The returned cloud carries positions and normals.
    ///
    /// # Arguments
    ///
    /// * `voxel_size`: the edge length of the grid cells, which must be finite and positive
    /// * `face_mask`: optional mask over the faces, restricting the sampling to those set true
    ///
    /// returns: `Result<PointCloud3>`
    pub fn sample_voxel_surface(
        &self,
        voxel_size: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<PointCloud3> {
        let sps = compute_voxel_surface_points(self.points(), self.faces(), voxel_size, face_mask)?;
        Ok(PointCloud3::from_surface_points(&sps))
    }
}

impl Mesh3 {
    /// Reduce the surface to one point per occupied cell of a regular voxel grid, placing each
    /// point on the surface itself.
    ///
    /// See [`MeshData3::sample_voxel_surface`], which this delegates to. Nothing here uses the
    /// mesh's BVH: the projection needs only the three vertices of the triangle it is projecting
    /// onto.
    pub fn sample_voxel_surface(
        &self,
        voxel_size: f64,
        face_mask: Option<&IndexMask>,
    ) -> Result<PointCloud3> {
        let sps = compute_voxel_surface_points(self.points(), self.faces(), voxel_size, face_mask)?;
        Ok(PointCloud3::from_surface_points(&sps))
    }
}

/// One point on the surface for every voxel of a regular grid which the given triangles pass
/// through, in ascending cell-key order, each carrying the normal of the face it lies on.
///
/// This is the primitive behind [`MeshData3::sample_voxel_surface`]; see there and the module
/// documentation for how a cell touched by several faces resolves.
///
/// # Arguments
///
/// * `points`: the vertex buffer
/// * `faces`: triangles as indices into `points`
/// * `voxel_size`: the edge length of the grid cells, which must be finite and positive
/// * `face_mask`: optional mask over `faces`, restricting the sampling to those set true
///
/// returns: `Result<Vec<SurfacePoint3>>`
pub fn compute_voxel_surface_points(
    points: &[Point3],
    faces: &[[u32; 3]],
    voxel_size: f64,
    face_mask: Option<&IndexMask>,
) -> Result<Vec<SurfacePoint3>> {
    validate(points, faces, voxel_size, face_mask)?;

    let best = faces
        .par_chunks(FACE_CHUNK)
        .enumerate()
        .map(|(chunk_i, chunk)| {
            let base = chunk_i * FACE_CHUNK;
            let mut local: HashMap<[i32; 3], CellBest> = HashMap::with_hasher(default());

            for (offset, face) in chunk.iter().enumerate() {
                let face_index = (base + offset) as u32;
                if skip_face(face_mask, base + offset) {
                    continue;
                }

                let tri = triangle(points, face);

                // A degenerate face has no normal to report and no area to occupy, so it is not
                // allowed to claim a cell. A cell touched by nothing else stays unoccupied.
                let Some(normal) = tri.normal() else {
                    continue;
                };

                rasterize_face(&tri, voxel_size, |key| {
                    let center = cell_center(key, voxel_size);
                    let point = tri.project_local_point(&center, false).point;
                    let candidate = CellBest {
                        d2: (point - center).norm_squared(),
                        face: face_index,
                        point,
                        normal,
                    };

                    match local.entry(key) {
                        Entry::Occupied(mut e) => {
                            if candidate.beats(e.get()) {
                                e.insert(candidate);
                            }
                        }
                        Entry::Vacant(e) => {
                            e.insert(candidate);
                        }
                    }
                });
            }

            local
        })
        .reduce(
            || HashMap::with_hasher(default()),
            |mut a, b| {
                for (key, candidate) in b {
                    match a.entry(key) {
                        Entry::Occupied(mut e) => {
                            if candidate.beats(e.get()) {
                                e.insert(candidate);
                            }
                        }
                        Entry::Vacant(e) => {
                            e.insert(candidate);
                        }
                    }
                }
                a
            },
        );

    let mut cells: Vec<([i32; 3], CellBest)> = best.into_iter().collect();
    cells.sort_unstable_by_key(|(key, _)| *key);

    Ok(cells
        .into_iter()
        .map(|(_, c)| SurfacePoint3::new(c.point, c.normal))
        .collect())
}

// ===============================================================================================
// Internals
// ===============================================================================================

/// The best claim on a cell seen so far: the nearest projection of the cell center onto any of the
/// triangles passing through it.
#[derive(Clone, Copy)]
struct CellBest {
    /// Squared distance from the cell center to `point`, which is what the claims are ranked by.
    d2: f64,

    /// The face `point` lies on, used only to break ties.
    face: u32,

    /// The projection of the cell center onto that face.
    point: Point3,

    /// That face's normal, carried here so the winning face does not have to be rebuilt later.
    normal: UnitVec3,
}

impl CellBest {
    /// Whether this claim displaces `other`.
    ///
    /// Ties in distance break toward the lower face index. Combined with the exact comparison on
    /// `d2` that makes this a total order, so merging two partial maps gives the same answer in
    /// either direction and the parallel result does not depend on how rayon scheduled the work.
    fn beats(&self, other: &CellBest) -> bool {
        self.d2 < other.d2 || (self.d2 == other.d2 && self.face < other.face)
    }
}

/// Whether a face is excluded by the caller's mask.
#[inline]
fn skip_face(face_mask: Option<&IndexMask>, face_index: usize) -> bool {
    face_mask.is_some_and(|m| !m.get(face_index))
}

/// Build the triangle for a face.
#[inline]
fn triangle(points: &[Point3], face: &[u32; 3]) -> Triangle {
    Triangle::new(
        points[face[0] as usize],
        points[face[1] as usize],
        points[face[2] as usize],
    )
}

/// The grid cell a coordinate falls into on one axis.
#[inline]
fn axis_key(coord: f64, voxel_size: f64) -> i32 {
    (coord / voxel_size).floor() as i32
}

/// The geometric center of a cell.
#[inline]
fn cell_center(key: [i32; 3], voxel_size: f64) -> Point3 {
    Point3::new(
        (key[0] as f64 + 0.5) * voxel_size,
        (key[1] as f64 + 0.5) * voxel_size,
        (key[2] as f64 + 0.5) * voxel_size,
    )
}

/// The bounds of a cell.
#[inline]
fn cell_aabb(key: [i32; 3], voxel_size: f64) -> Aabb {
    let mins = Point3::new(
        key[0] as f64 * voxel_size,
        key[1] as f64 * voxel_size,
        key[2] as f64 * voxel_size,
    );
    let maxs = Point3::new(
        mins.x + voxel_size,
        mins.y + voxel_size,
        mins.z + voxel_size,
    );
    Aabb::new(mins, maxs)
}

/// Call `emit` once for every cell of the grid which `tri` passes through.
///
/// The candidate range comes from the triangle's bounding box, and each candidate is confirmed with
/// a separating-axis test so that cells the box covers but the triangle misses are not reported.
/// The overwhelmingly common case on a mesh finer than the grid is a triangle inside a single cell,
/// which is detected from the key range alone and skips the test entirely.
fn rasterize_face(tri: &Triangle, voxel_size: f64, mut emit: impl FnMut([i32; 3])) {
    let mut lo = [0i32; 3];
    let mut hi = [0i32; 3];
    for d in 0..3 {
        let a = tri.a[d];
        let b = tri.b[d];
        let c = tri.c[d];
        lo[d] = axis_key(a.min(b).min(c), voxel_size);
        hi[d] = axis_key(a.max(b).max(c), voxel_size);
    }

    if lo == hi {
        emit(lo);
        return;
    }

    for i in lo[0]..=hi[0] {
        for j in lo[1]..=hi[1] {
            for k in lo[2]..=hi[2] {
                let key = [i, j, k];
                if intersection_test_aabb_triangle(&cell_aabb(key, voxel_size), tri) {
                    emit(key);
                }
            }
        }
    }
}

/// Reject the inputs a rasterization cannot be run on.
///
/// The interesting one is the coordinate range. Cell keys are `i32`, so a mesh far enough from the
/// origin relative to `voxel_size` would saturate the cast and collide unrelated cells; checking it
/// once here also bounds the per-face candidate loop, which is driven by the same keys and would
/// otherwise run for a very long time on a single wild vertex.
fn validate(
    points: &[Point3],
    faces: &[[u32; 3]],
    voxel_size: f64,
    face_mask: Option<&IndexMask>,
) -> Result<()> {
    if !voxel_size.is_finite() || voxel_size <= 0.0 {
        return Err(format!("Voxel size must be finite and positive, got {voxel_size}").into());
    }

    if let Some(mask) = face_mask
        && mask.len() != faces.len()
    {
        return Err(format!(
            "Face mask has {} entries but the mesh has {} faces",
            mask.len(),
            faces.len()
        )
        .into());
    }

    // A generous margin below the true limit, so that the `+/- 1` cells around a key are also
    // representable and nothing has to reason about the boundary.
    let limit = i32::MAX as f64 - 2.0;

    for p in points {
        for d in 0..3 {
            let scaled = p[d] / voxel_size;
            if !scaled.is_finite() || scaled.abs() > limit {
                return Err(format!(
                    "A vertex at {} is {:.3e} cells from the origin on axis {d}, which a voxel \
                     grid cannot index. Use a larger voxel size, or move the mesh nearer the \
                     origin.",
                    p, scaled
                )
                .into());
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Iso3;

    fn engine_blade() -> Mesh3 {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/engine-blade.tcmesh");
        let data = crate::io::read_tc_mesh_file(&path).expect("failed to load engine blade");
        Mesh3::from_data(data, false).expect("failed to accelerate engine blade")
    }

    /// Two parallel sheets a fraction of a cell apart, each spanning several cells laterally, so
    /// that every occupied cell is straddled by surfaces facing opposite ways.
    fn thin_wall(z_lower: f64, z_upper: f64) -> MeshData3 {
        let mut points = Vec::new();
        let mut faces = Vec::new();

        for z in [z_lower, z_upper] {
            let base = points.len() as u32;
            points.push(Point3::new(0.0, 0.0, z));
            points.push(Point3::new(4.0, 0.0, z));
            points.push(Point3::new(4.0, 4.0, z));
            points.push(Point3::new(0.0, 4.0, z));
            faces.push([base, base + 1, base + 2]);
            faces.push([base, base + 2, base + 3]);
        }

        MeshData3::new(points, faces).expect("the sheets are well formed")
    }

    /// Every cell the surface passes through, taken straight from the rasterizer.
    ///
    /// The sampler no longer reports cells, only the points inside them, and a projected point can
    /// land outside the cell that produced it. Tests about which cells are occupied therefore go to
    /// the rasterizer directly rather than trying to recover keys from the output.
    fn occupied_cells(points: &[Point3], faces: &[[u32; 3]], voxel_size: f64) -> Vec<[i32; 3]> {
        let mut keys = Vec::new();
        for face in faces {
            let tri = triangle(points, face);
            if tri.normal().is_none() {
                continue;
            }
            rasterize_face(&tri, voxel_size, |key| keys.push(key));
        }
        keys.sort_unstable();
        keys.dedup();
        keys
    }

    // ===========================================================================================
    // Where the points land
    // ===========================================================================================

    /// The whole reason this sampler projects rather than emitting cell centers, so it is pinned to
    /// floating point rather than to a tolerance. A projected point is on a triangle of the mesh by
    /// construction, so the mesh's own closest-point query has to agree.
    #[test]
    fn surface_sampling_lands_on_the_mesh() {
        for (name, mesh, voxel_size) in [
            ("sphere", Mesh3::create_sphere(10.0, 0.027).unwrap(), 1.0),
            ("bunny", Mesh3::stanford_bunny_res2(), 0.01),
            ("blade", engine_blade(), 2.0),
        ] {
            let cloud = mesh
                .sample_voxel_surface(voxel_size, None)
                .expect("sampling failed");

            assert!(!cloud.is_empty(), "{name} produced nothing");

            let worst = cloud
                .points()
                .iter()
                .map(|p| mesh.distance_closest_to(p))
                .fold(0.0, f64::max);

            assert!(
                worst < 1e-9,
                "{name}: a sampled point was {worst:.3e} from the surface"
            );
        }
    }

    /// The sampler emits one point per occupied cell and no more, which is the invariant that lets
    /// a caller reason about the output size from the grid alone.
    #[test]
    fn one_point_comes_out_per_occupied_cell() {
        let mesh = Mesh3::create_sphere(10.0, 0.027).unwrap();
        let voxel_size = 1.0;

        let cloud = mesh.sample_voxel_surface(voxel_size, None).unwrap();
        let cells = occupied_cells(mesh.points(), mesh.faces(), voxel_size);

        assert!(!cells.is_empty());
        assert_eq!(cloud.point_count(), cells.len());
    }

    // ===========================================================================================
    // Rasterization
    // ===========================================================================================

    /// Sweeping a region larger than the triangle checks the parts written here, which are the
    /// candidate loop bounds and the deduplication, against the same intersection test used
    /// inside. It does not independently verify parry's test, only that no cell that test accepts
    /// is missed and none it rejects is reported.
    #[test]
    fn rasterization_matches_a_brute_force_sweep() {
        let points = vec![
            Point3::new(0.3, 0.2, 0.1),
            Point3::new(7.4, 0.6, 2.2),
            Point3::new(1.1, 6.9, 5.3),
        ];
        let faces = vec![[0u32, 1, 2]];
        let voxel_size = 0.5;

        let got = occupied_cells(&points, &faces, voxel_size);

        let tri = Triangle::new(points[0], points[1], points[2]);
        let mut expected = Vec::new();
        for i in -3..=20 {
            for j in -3..=20 {
                for k in -3..=20 {
                    let key = [i, j, k];
                    if intersection_test_aabb_triangle(&cell_aabb(key, voxel_size), &tri) {
                        expected.push(key);
                    }
                }
            }
        }
        expected.sort_unstable();

        assert!(expected.len() > 100, "the sweep was too small to be a test");
        assert_eq!(got, expected);
    }

    /// A face has to occupy at least the cell its own centroid falls in. Run on a box, whose twelve
    /// triangles are each far larger than a cell, so this exercises the multi-cell path.
    #[test]
    fn every_face_centroid_lands_in_an_occupied_cell() {
        let mesh = Mesh3::create_box(3.0, 5.0, 7.0, true);
        let voxel_size = 0.25;

        let occupied = occupied_cells(mesh.points(), mesh.faces(), voxel_size);

        for centroid in mesh.compute_face_centers().unwrap() {
            let key = [
                axis_key(centroid.x, voxel_size),
                axis_key(centroid.y, voxel_size),
                axis_key(centroid.z, voxel_size),
            ];
            assert!(
                occupied.binary_search(&key).is_ok(),
                "cell {key:?} was not reported"
            );
        }
    }

    /// Rayon's scheduling and the hash map's iteration order must not reach the output, so the same
    /// input has to give a bitwise identical answer every time it is run.
    #[test]
    fn repeated_runs_are_bitwise_identical() {
        let mesh = Mesh3::stanford_bunny_res2();
        let voxel_size = 0.01;

        let surface =
            || compute_voxel_surface_points(mesh.points(), mesh.faces(), voxel_size, None);
        let a = surface().unwrap();
        let b = surface().unwrap();

        assert!(!a.is_empty());
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.point, y.point);
            assert_eq!(x.normal, y.normal);
        }
    }

    /// Which cells the surface occupies is a property of the geometry, so renumbering the faces
    /// cannot change it.
    ///
    /// The points inside those cells are a weaker guarantee, and the difference is worth pinning.
    /// Reordering changes face indices, which changes how an exact distance tie resolves, and ties
    /// are not rare: they happen wherever a cell center is equidistant from two faces. The winning
    /// position still agrees to floating point, because a tie means the two projections landed on
    /// the same place. The normal need not, since the tied faces can be opposite sides of a
    /// sub-voxel wall.
    #[test]
    fn the_occupied_cells_do_not_depend_on_face_order() {
        let mesh = Mesh3::stanford_bunny_res2();
        let voxel_size = 0.01;

        let mut reversed = mesh.faces().to_vec();
        reversed.reverse();

        let forward = occupied_cells(mesh.points(), mesh.faces(), voxel_size);
        let backward = occupied_cells(mesh.points(), &reversed, voxel_size);
        assert!(!forward.is_empty());
        assert_eq!(forward, backward);

        let forward =
            compute_voxel_surface_points(mesh.points(), mesh.faces(), voxel_size, None).unwrap();
        let backward =
            compute_voxel_surface_points(mesh.points(), &reversed, voxel_size, None).unwrap();

        assert_eq!(forward.len(), backward.len());
        for (a, b) in forward.iter().zip(backward.iter()) {
            assert!(
                (a.point - b.point).norm() < 1e-12,
                "a reordering moved a point from {} to {}",
                a.point,
                b.point
            );
        }
    }

    // ===========================================================================================
    // How a cell with several faces resolves
    // ===========================================================================================

    /// Every normal on a box must be an axis direction. An averaged normal at an edge cell would
    /// have two components near `1/sqrt(2)`, so this fails loudly if the selection rule is ever
    /// replaced by a blend.
    #[test]
    fn normals_come_from_one_facet_and_are_never_blended() {
        let mesh = Mesh3::create_box(10.0, 10.0, 10.0, true);
        let cloud = mesh.sample_voxel_surface(1.0, None).unwrap();

        let normals = cloud.point_normals().expect("normals were not written");
        assert!(!normals.is_empty());

        for n in normals {
            let longest = n.x.abs().max(n.y.abs()).max(n.z.abs());
            assert!(
                (longest - 1.0).abs() < 1e-12,
                "normal {n:?} is a blend of facets"
            );
        }
    }

    /// Across a wall thinner than a cell, every cell holds two sheets facing opposite ways.
    /// Selecting the nearer one keeps the point on a sheet; averaging would put it in the void
    /// between them and cancel the normals.
    #[test]
    fn a_thin_wall_puts_points_on_a_sheet_not_between_them() {
        let voxel_size = 1.0;
        let data = thin_wall(0.40, 0.55);

        let cloud = data.sample_voxel_surface(voxel_size, None).unwrap();
        assert!(!cloud.is_empty());

        for p in cloud.points() {
            // The cell centers sit at z = 0.5, so the upper sheet is nearer and must win outright.
            assert!(
                (p.z - 0.55).abs() < 1e-12,
                "point at z = {} is not on the nearer sheet (an average would be 0.475)",
                p.z
            );
        }

        for n in cloud.point_normals().unwrap() {
            assert!((n.z.abs() - 1.0).abs() < 1e-12, "normal {n:?} was blended");
        }
    }

    /// With the two sheets equidistant there is no nearer one, and the tie has to resolve the same
    /// way every run rather than by whichever thread got there first.
    #[test]
    fn an_exact_tie_breaks_toward_the_lower_face_index() {
        let voxel_size = 1.0;
        let data = thin_wall(0.45, 0.55);

        let cloud = data.sample_voxel_surface(voxel_size, None).unwrap();
        assert!(!cloud.is_empty());

        for p in cloud.points() {
            assert!(
                (p.z - 0.45).abs() < 1e-12,
                "point at z = {} did not come from the lower-indexed sheet",
                p.z
            );
        }
    }

    /// A zero-area face is not a surface, so it cannot claim a cell on its own.
    #[test]
    fn a_degenerate_face_claims_nothing() {
        let points = vec![
            Point3::new(0.1, 0.1, 0.1),
            Point3::new(0.2, 0.1, 0.1),
            Point3::new(0.3, 0.1, 0.1),
        ];
        let faces = vec![[0u32, 1, 2]];

        assert!(
            compute_voxel_surface_points(&points, &faces, 1.0, None)
                .unwrap()
                .is_empty()
        );
        assert!(occupied_cells(&points, &faces, 1.0).is_empty());
    }

    // ===========================================================================================
    // Arguments
    // ===========================================================================================

    #[test]
    fn a_face_mask_restricts_the_sampling() {
        let mesh = Mesh3::create_box(10.0, 10.0, 10.0, true);
        let voxel_size = 1.0;

        // The first two triangles are one side of the box, so masking them has to drop most of it.
        let mask = IndexMask::try_from_indices(&[0, 1], mesh.faces().len()).unwrap();

        let full = mesh.sample_voxel_surface(voxel_size, None).unwrap();
        let part = mesh.sample_voxel_surface(voxel_size, Some(&mask)).unwrap();

        assert!(!part.is_empty());
        assert!(part.point_count() < full.point_count() / 4);

        for p in part.points() {
            assert!(mesh.distance_closest_to(p) < 1e-9);
        }
    }

    #[test]
    fn an_empty_mesh_samples_to_an_empty_cloud() {
        let data = MeshData3::empty();
        assert!(data.sample_voxel_surface(1.0, None).unwrap().is_empty());
    }

    #[test]
    fn a_nonsense_voxel_size_is_rejected() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                mesh.sample_voxel_surface(bad, None).is_err(),
                "{bad} was accepted"
            );
        }
    }

    #[test]
    fn a_mask_of_the_wrong_length_is_rejected() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true);
        let wrong = IndexMask::new(mesh.faces().len() + 1, true);
        assert!(mesh.sample_voxel_surface(0.25, Some(&wrong)).is_err());
    }

    /// Cell keys are `i32`, so a mesh far enough from the origin relative to the voxel size has to
    /// be refused rather than silently colliding unrelated cells.
    #[test]
    fn a_mesh_too_far_from_the_origin_to_index_is_rejected() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true)
            .transform_copy(&Iso3::translation(1e6, 0.0, 0.0));

        // 1e12 cells out is past the i32 key space; 1e6 is comfortably inside it.
        assert!(mesh.sample_voxel_surface(1e-6, None).is_err());
        assert!(mesh.sample_voxel_surface(1.0, None).is_ok());
    }

    // ===========================================================================================
    // Reporting
    // ===========================================================================================

    /// Best-of-three wall time for `f`, because a single shot on this workload has a noise floor
    /// wide enough to invert a comparison.
    fn time_best<T>(mut f: impl FnMut() -> T) -> (T, f64) {
        let mut best = f64::INFINITY;
        let mut out = None;
        for _ in 0..3 {
            let t = std::time::Instant::now();
            let r = f();
            best = best.min(t.elapsed().as_secs_f64() * 1000.0);
            out = Some(r);
        }
        (out.expect("ran at least once"), best)
    }

    /// How far a cloud sits from the surface it came from.
    fn deviation(mesh: &Mesh3, cloud: &PointCloud3) -> (f64, f64) {
        let d: Vec<f64> = cloud
            .points()
            .iter()
            .map(|p| mesh.distance_closest_to(p))
            .collect();
        let rms = (d.iter().map(|x| x * x).sum::<f64>() / d.len() as f64).sqrt();
        (rms, d.iter().cloned().fold(0.0, f64::max))
    }

    /// The measurement that decides whether this module earns its place, since the route
    /// `sample_dense(v/2)` -> `reduce_by_voxel(v)` already produced a
    /// comparable cloud before it existed. Run with:
    ///
    /// ```text
    /// cargo test -r -p engeom --lib voxelize -- --ignored --nocapture
    /// ```
    #[test]
    #[ignore]
    fn report_against_the_existing_dense_sample_route() {
        let mesh = engine_blade();
        println!(
            "\nengine blade: {} points, {} faces\n",
            mesh.point_count(),
            mesh.face_count()
        );
        println!(
            "{:>6}  {:<22} {:>10} {:>10} {:>12} {:>12}",
            "v", "route", "points", "ms", "rms dev", "max dev"
        );

        for v in [4.0, 2.0, 1.0, 0.5] {
            let (dense, dense_ms) =
                time_best(|| mesh.sample_dense(v / 2.0, None).reduce_by_voxel(v).unwrap());
            let (surface, surface_ms) = time_best(|| mesh.sample_voxel_surface(v, None).unwrap());

            for (name, cloud, ms) in [
                ("dense sample + reduce", &dense, dense_ms),
                ("voxel surface", &surface, surface_ms),
            ] {
                let (rms, max) = deviation(&mesh, cloud);
                println!(
                    "{:>6}  {:<22} {:>10} {:>10.1} {:>12.3e} {:>12.3e}",
                    v,
                    name,
                    cloud.point_count(),
                    ms,
                    rms,
                    max
                );
            }
            println!();
        }
    }
}
