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

use crate::container::{self, Kind, item};
use crate::error::{Error, Result};
use crate::indices::{read_indices, write_indices};
use crate::metadata::Metadata;
use crate::points::{read_points, write_points};
use crate::raw::MAX_PREALLOC;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

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

/// Write a collection of meshes, every item at the same tolerance.
///
/// # Errors
///
/// [`Error::ToleranceNotRepresentable`] if any axis is too wide to meet `tol`, and
/// [`Error::Malformed`] for a non-finite coordinate or a face pointing past the vertex array.
pub fn write_to<W: Write>(writer: &mut W, meshes: &[Mesh3], tol: f64) -> Result<()> {
    write_to_with_meta(writer, meshes, tol, &Metadata::new())
}

/// Write a collection of meshes with file-level metadata attached. See [`write_to`].
pub fn write_to_with_meta<W: Write>(
    writer: &mut W,
    meshes: &[Mesh3],
    tol: f64,
    file_metadata: &Metadata,
) -> Result<()> {
    let count = u32::try_from(meshes.len())
        .map_err(|_| Error::Malformed("container holds more items than a u32 can count"))?;
    container::write_header(writer, Kind::Mesh3, count, file_metadata)?;

    for mesh in meshes {
        write_item(writer, mesh, tol)?;
    }

    Ok(())
}

/// Write a single mesh as a one-item collection.
pub fn write_one_to<W: Write>(writer: &mut W, mesh: &Mesh3, tol: f64) -> Result<()> {
    write_to(writer, std::slice::from_ref(mesh), tol)
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
    // Buffering is required, not tidiness: the bit reader and writer move a byte at a time, so an
    // unbuffered File would mean a syscall per byte.
    let mut writer = BufWriter::new(File::create(path)?);
    write_to(&mut writer, meshes, tol)?;
    writer.flush()?;
    Ok(())
}

/// Write a single mesh to a file as a one-item collection.
pub fn write_one_file(path: &Path, mesh: &Mesh3, tol: f64) -> Result<()> {
    write_file(path, std::slice::from_ref(mesh), tol)
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

fn write_item<W: Write>(writer: &mut W, mesh: &Mesh3, tol: f64) -> Result<()> {
    item::write_preamble(writer, mesh.name.as_deref(), &mesh.metadata, false)?;

    write_points(writer, &mesh.points, tol)?;

    // The limit is this item's own vertex count, so the encoder refuses a face that points past
    // the end of the array it will be decoded against.
    let limit = u32::try_from(mesh.points.len())
        .map_err(|_| Error::Malformed("mesh holds more vertices than a u32 can index"))?;
    write_indices(writer, &mesh.faces, limit)?;

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
            assert_matches(&mesh, &back, case.tol, case.name);
        }
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
            assert_matches(o, r, 1e-3, "collection");
        }
    }

    #[test]
    fn an_empty_collection_round_trips() {
        let mut buf = Vec::new();
        write_to(&mut buf, &[], 1e-3).unwrap();

        assert_eq!(buf.len(), 12, "an empty collection is just a header");
        assert!(read_from(&mut Cursor::new(&buf)).unwrap().is_empty());
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

        assert_matches(&mesh, &back, case.tol, "file");

        // And the collection path over the same file.
        write_file(&path.0, &[mesh.clone(), mesh.clone()], case.tol).unwrap();
        let all = read_file(&path.0).unwrap();
        assert_eq!(all.len(), 2);
        assert_matches(&mesh, &all[1], case.tol, "file collection");
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
