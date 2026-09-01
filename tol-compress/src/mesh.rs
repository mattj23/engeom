//! Triangle mesh containers, conventionally `.tcmesh`.
//!
//! An item is a points block and an index block, with an optional name. Faces are validated
//! against the item's own vertex count on both write and read, so a file can never carry a face
//! pointing past the end of its vertex array.
//!
//! # Collections
//!
//! Like every container in this crate, a mesh file is an ordered collection. Most hold one item,
//! which is what [`read_one_from`] and [`write_one_to`] are for, but that is an API convenience
//! rather than a different format. Reading one item out of a file holding several **errors**
//! rather than silently discarding the rest.
//!
//! # Tolerance
//!
//! [`write_to`] applies one tolerance to every item. The format itself does not require that,
//! since bit widths are recorded per block, so a future API could give each item its own tolerance
//! without a format change.
//!
//! # Vertex and face order
//!
//! <div class="warning">
//!
//! Unless you use the non-reordering version, **writing a mesh does not preserve the order of its
//! vertices or its faces.**
//!
//! </div>
//!
//! [`write_to`] hands the mesh to [`crate::reorder`] first, which renumbers vertices so that
//! indices compress, and that is worth -25% to -50% of the index block on real meshes. The mesh
//! that comes back out is the same geometry: the same triangles over the same positions, with
//! everything renumbered in a consistent way. Nothing about the reordering is stored, so decoding
//! costs nothing extra but if you want the original order persisted in some way you need to do it
//! yourself.
//!
//! To avoid the reordering version, use [`write_to_preserving_order`] instead. It writes the arrays
//! exactly as they were given, at the cost of a larger index block. The block records which coding
//! it used so the readers are unaffected.
//!
//! ## Keeping external data in step, without giving up the saving
//!
//! There is a third option, and it is usually the one you want if you have per-vertex data the
//! crate does not hold: **do the renumbering yourself**, move your own arrays through the same
//! permutation, then ask the encoder not to do it again. The file is byte-for-byte identical to
//! what [`write_to`] would have produced, because [`crate::indices`] picks the compact coding
//! whenever its input is already in first-use order.
//!
//! ```
//! use tol_compress::{Mesh3, mesh, reorder};
//!
//! let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]];
//! let faces = vec![[0u32, 1, 2], [1, 3, 2]];
//!
//! // Something the crate knows nothing about, one entry per vertex.
//! let scan_pass = vec![10u32, 11, 12, 13];
//!
//! let plan = reorder::optimize(&faces, points.len())?;
//! let points = reorder::permute(&points, &plan.vertex_order);
//! let scan_pass = reorder::permute(&scan_pass, &plan.vertex_order);
//!
//! let mut buf = Vec::new();
//! let mesh = Mesh3::new(points, plan.faces);
//! mesh::write_one_to_preserving_order(&mut buf, &mesh, 1e-6)?;
//!
//! // `scan_pass[i]` still belongs to vertex `i` of the mesh that was written.
//! assert_eq!(scan_pass.len(), mesh.points.len());
//! # Ok::<(), tol_compress::Error>(())
//! ```
//!
//! Face-domain data moves through `plan.face_order` the same way.

use crate::container::{self, Kind, Named, item};
use crate::effort::Effort;
use crate::error::{Error, Result};
use crate::indices::{read_indices, write_indices};
use crate::metadata::Metadata;
use crate::points::{read_points, write_points_with};
use crate::raw::MAX_PREALLOC;
use crate::reorder;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Whether the encoder may renumber a mesh's vertices and faces to make its indices smaller.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum VertexOrder {
    /// Renumber for size. The default, and much the smaller of the two.
    #[default]
    Optimize,
    /// Write the arrays exactly as given, so vertex and face indices mean the same thing after a
    /// round trip as they did before it.
    Preserve,
}

/// Everything [`write_to_with`] can be told beyond the geometry and the tolerance.
///
/// A struct rather than more arguments, because the settings are independent and none of them is
/// wanted often enough to deserve its own function on every entry point.
///
/// Non-exhaustive, so that later settings do not break callers. Build one with [`WriteOptions::new`]
/// and the `with_` methods, or from [`Default`], and set the fields on it directly from there.
#[non_exhaustive]
#[derive(Debug, Clone, Default)]
pub struct WriteOptions {
    /// File-level metadata. Stored, never interpreted.
    pub metadata: Metadata,
    /// Whether vertices and faces may be renumbered. Semantic, so it is never folded into
    /// [`Effort`]: it changes what the caller gets back rather than only what it costs.
    pub order: VertexOrder,
    /// How hard to search for a smaller file. Changes size and encoding time, nothing else.
    pub effort: Effort,
}

impl WriteOptions {
    /// The defaults: no metadata, vertices renumbered for size, [`Effort::Balanced`].
    pub fn new() -> Self {
        Self::default()
    }

    /// The same options carrying file-level metadata.
    pub fn with_metadata(mut self, metadata: Metadata) -> Self {
        self.metadata = metadata;
        self
    }

    /// The same options at a different search effort.
    pub fn with_effort(mut self, effort: Effort) -> Self {
        self.effort = effort;
        self
    }

    /// The same options with vertex and face order left as the caller gave it.
    pub fn preserving_order(mut self) -> Self {
        self.order = VertexOrder::Preserve;
        self
    }
}

/// A triangle mesh with an optional name.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Mesh3 {
    /// Optional identifier, preserved through a round trip. Order is meaningful regardless.
    pub name: Option<String>,
    /// Vertex positions.
    pub points: Vec<[f64; 3]>,
    /// Triangles, as indices into `points`.
    pub faces: Vec<[u32; 3]>,
    /// Caller metadata, empty when the mesh carries none. Stored, never interpreted.
    pub metadata: Metadata,
}

impl Mesh3 {
    /// A mesh with no name.
    pub fn new(points: Vec<[f64; 3]>, faces: Vec<[u32; 3]>) -> Self {
        Self {
            name: None,
            points,
            faces,
            metadata: Metadata::new(),
        }
    }

    /// The same mesh carrying a name.
    pub fn named(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    /// The same mesh carrying one metadata entry.
    pub fn with_meta(mut self, key: impl Into<String>, value: impl Into<crate::Value>) -> Self {
        self.metadata.insert(key.into(), value.into());
        self
    }
}

impl Named for Mesh3 {
    fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }
}

/// Write a collection of meshes, every item at the same tolerance.
///
/// Vertices and faces are renumbered for size; see the module documentation and
/// [`write_to_preserving_order`].
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if any axis is too wide to meet `tol`, and
/// [`Error::Malformed`] for a non-finite coordinate or a face pointing past the vertex array.
pub fn write_to<W: Write>(writer: &mut W, meshes: &[Mesh3], tol: f64) -> Result<()> {
    write_to_with(writer, meshes, tol, &WriteOptions::default())
}

/// Write a collection of meshes without renumbering anything. See [`write_to`].
pub fn write_to_preserving_order<W: Write>(
    writer: &mut W,
    meshes: &[Mesh3],
    tol: f64,
) -> Result<()> {
    write_to_with(writer, meshes, tol, &WriteOptions::new().preserving_order())
}

/// Write a collection of meshes with full control over metadata and ordering. See [`write_to`].
pub fn write_to_with<W: Write>(
    writer: &mut W,
    meshes: &[Mesh3],
    tol: f64,
    options: &WriteOptions,
) -> Result<()> {
    let count = u32::try_from(meshes.len())
        .map_err(|_| Error::Malformed("container holds more items than a u32 can count"))?;
    container::write_header(writer, Kind::Mesh3, count, &options.metadata)?;

    for mesh in meshes {
        write_item(writer, mesh, tol, options)?;
    }

    Ok(())
}

/// Write a single mesh as a one-item collection.
pub fn write_one_to<W: Write>(writer: &mut W, mesh: &Mesh3, tol: f64) -> Result<()> {
    write_to(writer, std::slice::from_ref(mesh), tol)
}

/// Write a single mesh without renumbering anything. See [`write_to`].
pub fn write_one_to_preserving_order<W: Write>(
    writer: &mut W,
    mesh: &Mesh3,
    tol: f64,
) -> Result<()> {
    write_to_preserving_order(writer, std::slice::from_ref(mesh), tol)
}

/// Read a collection of meshes.
pub fn read_from<R: Read>(reader: &mut R) -> Result<Vec<Mesh3>> {
    let header = container::read_header(reader, Kind::Mesh3)?;

    let mut out = Vec::with_capacity((header.count as usize).min(MAX_PREALLOC));
    for _ in 0..header.count {
        out.push(read_item(reader)?);
    }

    Ok(out)
}

/// Read a container that holds exactly one mesh.
///
/// # Errors
///
/// [`Error::NotASingleItem`] if the container holds any other number, because returning the first
/// of several would silently discard data.
pub fn read_one_from<R: Read>(reader: &mut R) -> Result<Mesh3> {
    let header = container::read_header(reader, Kind::Mesh3)?;
    if header.count != 1 {
        return Err(Error::NotASingleItem {
            found: header.count,
        });
    }
    read_item(reader)
}

/// Write a collection of meshes to a file. See [`write_to`].
pub fn write_file(path: &Path, meshes: &[Mesh3], tol: f64) -> Result<()> {
    write_file_with(path, meshes, tol, &WriteOptions::default())
}

/// Write a collection of meshes to a file with full control over metadata, ordering and effort.
///
/// The counterpart of [`write_to_with`] for callers who are writing a path rather than a writer,
/// which is otherwise the only way to reach those settings.
pub fn write_file_with(
    path: &Path,
    meshes: &[Mesh3],
    tol: f64,
    options: &WriteOptions,
) -> Result<()> {
    // Buffering is required, not tidiness: the bit reader and writer move a byte at a time, so an
    // unbuffered File would mean a syscall per byte.
    let mut writer = BufWriter::new(File::create(path)?);
    write_to_with(&mut writer, meshes, tol, options)?;
    writer.flush()?;
    Ok(())
}

/// Write a collection of meshes to a file without renumbering anything. See [`write_to`].
pub fn write_file_preserving_order(path: &Path, meshes: &[Mesh3], tol: f64) -> Result<()> {
    write_file_with(path, meshes, tol, &WriteOptions::new().preserving_order())
}

/// Write a single mesh to a file as a one-item collection.
pub fn write_one_file(path: &Path, mesh: &Mesh3, tol: f64) -> Result<()> {
    write_file(path, std::slice::from_ref(mesh), tol)
}

/// Write a single mesh to a file without renumbering anything. See [`write_to`].
pub fn write_one_file_preserving_order(path: &Path, mesh: &Mesh3, tol: f64) -> Result<()> {
    write_file_preserving_order(path, std::slice::from_ref(mesh), tol)
}

/// Read a collection of meshes from a file.
pub fn read_file(path: &Path) -> Result<Vec<Mesh3>> {
    let mut reader = BufReader::new(File::open(path)?);
    read_from(&mut reader)
}

/// Read a file holding exactly one mesh. See [`read_one_from`].
pub fn read_one_file(path: &Path) -> Result<Mesh3> {
    let mut reader = BufReader::new(File::open(path)?);
    read_one_from(&mut reader)
}

fn write_item<W: Write>(
    writer: &mut W,
    mesh: &Mesh3,
    tol: f64,
    options: &WriteOptions,
) -> Result<()> {
    item::write_preamble(writer, mesh.name.as_deref(), &mesh.metadata, false)?;

    // The limit is this item's own vertex count, so the encoder refuses a face that points past
    // the end of the array it will be decoded against.
    let limit = u32::try_from(mesh.points.len())
        .map_err(|_| Error::Malformed("mesh holds more vertices than a u32 can index"))?;

    // Renumbering has to move the points and the faces together, which is why it happens here,
    // where both are in hand, rather than inside either encoder. When attributes arrive they are
    // permuted at this same point, through `reorder::permute`.
    if options.order == VertexOrder::Optimize && !mesh.faces.is_empty() {
        let plan = reorder::optimize(&mesh.faces, mesh.points.len())?;
        let points = reorder::permute(&mesh.points, &plan.vertex_order);

        write_points_with(writer, &points, tol, options.effort)?;
        write_indices(writer, &plan.faces, limit)?;
    } else {
        write_points_with(writer, &mesh.points, tol, options.effort)?;
        write_indices(writer, &mesh.faces, limit)?;
    }

    Ok(())
}

fn read_item<R: Read>(reader: &mut R) -> Result<Mesh3> {
    let preamble = item::read_preamble(reader, false)?;

    let points: Vec<[f64; 3]> = read_points(reader)?;
    let limit = u32::try_from(points.len())
        .map_err(|_| Error::Malformed("mesh holds more vertices than a u32 can index"))?;
    let faces = read_indices(reader, limit)?;

    Ok(Mesh3 {
        name: preamble.name,
        points,
        faces,
        metadata: preamble.metadata,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::corpus;
    use crate::testgen::Rng;
    use std::io::Cursor;

    /// A scratch file that removes itself, so a failing assertion does not leave litter behind.
    struct TempPath(std::path::PathBuf);

    impl TempPath {
        fn new(tag: &str) -> Self {
            let mut rng = Rng::new(std::process::id() as u64);
            let name = format!("tol-compress-{tag}-{:x}.tcmesh", rng.next_u64());
            Self(std::env::temp_dir().join(name))
        }
    }

    impl Drop for TempPath {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    fn distance(a: &[f64; 3], b: &[f64; 3]) -> f64 {
        (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
    }

    /// A round trip through [`write_to_preserving_order`], where the arrays come back as they went
    /// in and index-for-index comparison is the right check.
    fn assert_matches(original: &Mesh3, recovered: &Mesh3, tol: f64, what: &str) {
        assert_eq!(recovered.name, original.name, "{what}: name");
        assert_eq!(
            recovered.points.len(),
            original.points.len(),
            "{what}: vertex count"
        );
        assert_eq!(
            recovered.faces, original.faces,
            "{what}: faces must be exact"
        );

        for (i, (o, r)) in original
            .points
            .iter()
            .zip(recovered.points.iter())
            .enumerate()
        {
            let d = distance(o, r);
            assert!(d <= tol, "{what}: vertex {i} recovered {d} away, tol {tol}");
        }
    }

    /// A round trip through [`write_to`], which renumbers. Array equality is the wrong question
    /// here; the right one is whether it is the same mesh.
    ///
    /// Checked as geometry rather than as indices: every triangle must join the same three
    /// positions it joined before. That is what catches the failure this whole design risks, which
    /// is a points array permuted out of step with the faces that index it.
    fn assert_same_mesh(original: &Mesh3, recovered: &Mesh3, tol: f64, what: &str) {
        assert_eq!(recovered.name, original.name, "{what}: name");
        assert_eq!(
            recovered.points.len(),
            original.points.len(),
            "{what}: vertex count"
        );
        assert_eq!(
            recovered.faces.len(),
            original.faces.len(),
            "{what}: face count"
        );

        let plan = reorder::optimize(&original.faces, original.points.len()).unwrap();

        for (new, &old) in plan.face_order.iter().enumerate() {
            for corner in 0..3 {
                let got = recovered.points[recovered.faces[new][corner] as usize];
                let want = original.points[original.faces[old as usize][corner] as usize];

                let d = distance(&got, &want);
                assert!(
                    d <= tol,
                    "{what}: face {new} corner {corner} landed {d} away, tol {tol}"
                );
            }
        }

        // Vertices no face references have to survive too, and nothing above would notice them.
        for (i, &old) in plan.vertex_order.iter().enumerate() {
            let d = distance(&recovered.points[i], &original.points[old as usize]);
            assert!(d <= tol, "{what}: vertex {i} recovered {d} away, tol {tol}");
        }
    }

    #[test]
    fn every_corpus_mesh_round_trips() {
        for case in corpus::all() {
            let mesh = Mesh3::new(case.points.clone(), case.faces.clone());

            let mut buf = Vec::new();
            write_one_to(&mut buf, &mesh, case.tol).unwrap();

            let mut cursor = Cursor::new(&buf);
            let back = read_one_from(&mut cursor).unwrap();

            assert_eq!(
                cursor.position() as usize,
                buf.len(),
                "{}: decoder left bytes unread",
                case.name
            );
            assert_same_mesh(&mesh, &back, case.tol, case.name);
        }
    }

    /// The same corpus through the opt-out, where the arrays must come back untouched.
    #[test]
    fn every_corpus_mesh_round_trips_preserving_order() {
        for case in corpus::all() {
            let mesh = Mesh3::new(case.points.clone(), case.faces.clone());

            let mut buf = Vec::new();
            write_one_to_preserving_order(&mut buf, &mesh, case.tol).unwrap();

            let back = read_one_from(&mut Cursor::new(&buf)).unwrap();
            assert_matches(&mesh, &back, case.tol, case.name);
        }
    }

    /// Renumbering is the whole point, so it has to actually pay on real geometry.
    ///
    /// It no longer pays on every mesh, and that is a real interaction rather than a slack test.
    /// Renumbering orders vertices for the index block, by breadth-first traversal of connectivity.
    /// The points block now wants them ordered for *space*, because contiguous runs with tight boxes
    /// are what [`crate::segment`] trades on. On a mesh whose connectivity is a poor guide to
    /// position, which `boundary_heavy` is by construction, the traversal scrambles the caller's
    /// spatially coherent order and the points block loses more than the index block gains.
    ///
    /// So this asserts that renumbering pays across the corpus and on the clear majority of meshes
    /// in it, not that it pays on all of them.
    #[test]
    fn renumbering_makes_meshes_smaller() {
        let mut wins = 0;
        let mut losses = Vec::new();
        let (mut total_optimized, mut total_preserved) = (0usize, 0usize);

        for case in corpus::all() {
            if case.faces.is_empty() || case.faces.len() < 100 {
                continue;
            }
            let mesh = Mesh3::new(case.points.clone(), case.faces.clone());

            let mut optimized = Vec::new();
            write_one_to(&mut optimized, &mesh, case.tol).unwrap();

            let mut preserved = Vec::new();
            write_one_to_preserving_order(&mut preserved, &mesh, case.tol).unwrap();

            total_optimized += optimized.len();
            total_preserved += preserved.len();

            if optimized.len() < preserved.len() {
                wins += 1;
            } else {
                losses.push((case.name, optimized.len(), preserved.len()));
            }
        }

        assert!(
            total_optimized < total_preserved,
            "renumbering must pay across the corpus: {total_optimized} against {total_preserved}"
        );
        assert!(
            wins > 2 * losses.len(),
            "renumbering should still win on most meshes, lost on {losses:?}"
        );
    }

    #[test]
    fn a_collection_preserves_order_and_names() {
        let smooth = corpus::smooth_surface();
        let planar = corpus::planar();

        let meshes = vec![
            Mesh3::new(smooth.points.clone(), smooth.faces.clone()).named("smooth"),
            Mesh3::new(planar.points.clone(), planar.faces.clone()),
            Mesh3::new(smooth.points.clone(), smooth.faces.clone()).named("smooth again"),
        ];

        let mut buf = Vec::new();
        write_to(&mut buf, &meshes, 1e-3).unwrap();

        let back = read_from(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(back.len(), 3);
        assert_eq!(back[0].name.as_deref(), Some("smooth"));
        assert_eq!(back[1].name, None);
        assert_eq!(back[2].name.as_deref(), Some("smooth again"));

        for (o, r) in meshes.iter().zip(back.iter()) {
            assert_same_mesh(o, r, 1e-3, "collection");
        }
    }

    #[test]
    fn an_empty_collection_round_trips() {
        let mut buf = Vec::new();
        write_to(&mut buf, &[], 1e-3).unwrap();

        assert_eq!(buf.len(), 12, "an empty collection is just a header");
        assert!(read_from(&mut Cursor::new(&buf)).unwrap().is_empty());
    }

    /// The recipe documented in the module header has to cost nothing, or callers with external
    /// per-vertex data are choosing between correctness and size.
    ///
    /// This is the load-bearing test for that documentation: without it, a later change to how
    /// [`crate::indices`] chooses a coding could quietly make the caller-driven path produce the
    /// larger file, and the docs would go on promising otherwise.
    #[test]
    fn the_caller_driven_recipe_is_byte_identical() {
        for case in corpus::all() {
            if case.faces.is_empty() {
                continue;
            }
            let mesh = Mesh3::new(case.points.clone(), case.faces.clone());

            let mut internal = Vec::new();
            write_one_to(&mut internal, &mesh, case.tol).unwrap();

            let plan = reorder::optimize(&mesh.faces, mesh.points.len()).unwrap();
            let by_hand = Mesh3::new(
                reorder::permute(&mesh.points, &plan.vertex_order),
                plan.faces,
            );
            let mut manual = Vec::new();
            write_one_to_preserving_order(&mut manual, &by_hand, case.tol).unwrap();

            assert_eq!(
                internal, manual,
                "{}: reordering by hand produced a different file",
                case.name
            );
        }
    }

    /// Taking the first of several would silently discard data, so it is refused.
    #[test]
    fn read_one_refuses_a_multi_item_file() {
        let m = Mesh3::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let mut buf = Vec::new();
        write_to(&mut buf, &[m.clone(), m.clone()], 1e-3).unwrap();

        assert!(matches!(
            read_one_from(&mut Cursor::new(&buf)),
            Err(Error::NotASingleItem { found: 2 })
        ));
    }

    #[test]
    fn read_one_refuses_an_empty_file() {
        let mut buf = Vec::new();
        write_to(&mut buf, &[], 1e-3).unwrap();

        assert!(matches!(
            read_one_from(&mut Cursor::new(&buf)),
            Err(Error::NotASingleItem { found: 0 })
        ));
    }

    /// Faces are checked against this item's own vertex count, not some ambient number, which is
    /// what keeps a collection's items independent.
    #[test]
    fn a_face_past_the_vertex_array_is_refused() {
        let mesh = Mesh3::new(
            vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 3]],
        );

        let mut buf = Vec::new();
        assert!(matches!(
            write_one_to(&mut buf, &mesh, 1e-3),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn a_mesh_with_no_faces_is_fine() {
        let mesh = Mesh3::new(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], Vec::new());

        let mut buf = Vec::new();
        write_one_to(&mut buf, &mesh, 1e-3).unwrap();
        let back = read_one_from(&mut Cursor::new(&buf)).unwrap();

        assert_matches(&mesh, &back, 1e-3, "no faces");
    }

    #[test]
    fn round_trips_through_a_real_file() {
        let case = corpus::smooth_surface();
        let mesh = Mesh3::new(case.points.clone(), case.faces.clone()).named("on disk");
        let path = TempPath::new("single");

        write_one_file(&path.0, &mesh, case.tol).unwrap();
        let back = read_one_file(&path.0).unwrap();

        assert_same_mesh(&mesh, &back, case.tol, "file");

        // And the collection path over the same file.
        write_file(&path.0, &[mesh.clone(), mesh.clone()], case.tol).unwrap();
        let all = read_file(&path.0).unwrap();
        assert_eq!(all.len(), 2);
        assert_same_mesh(&mesh, &all[1], case.tol, "file collection");

        // And the order-preserving file path, which must hand the arrays back untouched.
        write_one_file_preserving_order(&path.0, &mesh, case.tol).unwrap();
        assert_matches(
            &mesh,
            &read_one_file(&path.0).unwrap(),
            case.tol,
            "preserved",
        );
    }

    #[test]
    fn a_different_kind_is_rejected() {
        let mut buf = Vec::new();
        container::write_header(&mut buf, Kind::Cloud3, 0, &Metadata::new()).unwrap();

        assert!(matches!(
            read_from(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    /// The closed and chord-tolerance bits belong to polylines. A mesh item setting them is a
    /// malformed file, not something to quietly ignore.
    #[test]
    fn polyline_only_flags_are_rejected() {
        let mesh = Mesh3::new(vec![[0.0; 3]], Vec::new());
        let mut buf = Vec::new();
        write_one_to(&mut buf, &mesh, 1e-3).unwrap();

        // The flag byte is the first thing after the 12 byte header.
        buf[12] |= item::CLOSED;

        assert!(matches!(
            read_one_from(&mut Cursor::new(&buf)),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn truncated_input_is_an_error() {
        let case = corpus::planar();
        let mesh = Mesh3::new(case.points.clone(), case.faces.clone()).named("truncate me");

        let mut buf = Vec::new();
        write_one_to(&mut buf, &mesh, case.tol).unwrap();

        for cut in [0, 5, 11, 12, 20, 60, buf.len() / 2, buf.len() - 1] {
            assert!(
                read_one_from(&mut Cursor::new(&buf[..cut])).is_err(),
                "truncating to {cut} bytes should fail"
            );
        }
    }

    #[test]
    fn an_absurd_item_count_does_not_allocate() {
        let mut buf = Vec::new();
        container::write_header(&mut buf, Kind::Mesh3, u32::MAX, &Metadata::new()).unwrap();

        assert!(
            read_from(&mut Cursor::new(&buf)).is_err(),
            "should run out of input, not memory"
        );
    }
}
