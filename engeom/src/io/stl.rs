#![cfg(feature = "stl")]

//! STL reading and writing.
//!
//! STL is the lowest common denominator of triangle mesh formats: an unindexed list of triangles
//! with 32-bit float coordinates and no per-element data of any kind. It's historically been
//! supported just about everywhere, which is the only reason to use it. Even in binary format it's
//! wildly inefficient and all it can hold is an unstructured list of faces. Even in the best case
//! vertices are duplicated multiple times.
//!
//! Only use STL when you're importing from or exporting to systems that don't support other
//! options. Even a binary list of vertices and faces would be preferable, because then at least
//! points wouldn't be duplicated.
//!
//! Also note that in binary format the coordinates are 32-bit floating points, so precision is
//! lost with distance from the origin.  In real units, the quantization precision grows by about
//! 1 micron every 10 m from the origin (this ratio holds across units), meaning that you need to
//! make an effort to save STL files with the points transformed to reasonably close to the origin.
//!
//! Everything else...attributes, colors, etc...are lost.  Only use this format if you have to.

use crate::geom3::mesh::MeshData3;
use crate::geom3::mesh::algorithms::normals::compute_face_normal;
use crate::{Mesh3, Point3, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, Write};
use std::path::Path;

// ===============================================================================================
// Reading
// ===============================================================================================

/// Read a triangle mesh from an STL file, in either the ascii or binary encoding.
///
/// The encoding is detected from the contents rather than trusted from the extension. Points are
/// recovered by welding triangle corners with exactly equal coordinates, which needs no tolerance
/// to tune and cannot collapse a genuine thin feature. Note that `0.0` and `-0.0` have different
/// bit patterns and so do not weld with each other, and that the file's per-facet normals are
/// discarded, since they are redundant with the winding and are frequently wrong in files written
/// by other tools.
///
/// # Arguments
///
/// * `path`: the path to the STL file
///
/// returns: `Result<MeshData3>`
pub fn load_stl_mesh_data(path: &Path) -> Result<MeshData3> {
    let mut file = BufReader::new(File::open(path)?);
    read_stl_mesh_data(&mut file)
}

/// Read a triangle mesh from any seekable source of STL data. See `load_stl_mesh_data`.
pub fn read_stl_mesh_data<R: Read + Seek>(source: &mut R) -> Result<MeshData3> {
    // TODO: drop the `stl_io` dependency and read the format directly.
    //
    // The writer above is ours already, so half the format is implemented here and half is not.
    // This call is the only remaining use of the crate, and replacing it needs three small pieces:
    // an ascii/binary probe, a binary reader, and an ascii tokenizer. Roughly 120 lines. The weld
    // is a `HashMap` keyed on `f32::to_bits`, which the round trip tests in this module already
    // cover.
    //
    // Three reasons to do it, beyond the asymmetry:
    //
    // - The `Seek` bound comes from this call, and it stops `read_stl_mesh_data` accepting a plain
    //   stream the way `read_ply_mesh_data` does. A pipe or a network body has to be buffered into
    //   a `Cursor` first.
    // - The ascii probe demands the literal "solid " including the trailing space, so a file
    //   beginning `solid\n` is routed to the binary reader and fails with an EOF error that says
    //   nothing about what is actually wrong. It is why `DEFAULT_HEADER` exists.
    // - Errors would be ours, so a malformed facet could report its byte offset.
    //
    // Worth doing alongside the `Mesh3` migration, when the `Mesh3` bridges in this module get
    // deleted anyway. Not urgent: the crate is small, and a hand written weld is the kind of thing
    // that fails silently by changing topology rather than loudly, so the replacement wants the
    // same care the writer got.
    let raw = stl_io::read_stl(source)?;

    let points = raw
        .vertices
        .iter()
        .map(|v| Point3::new(v[0] as f64, v[1] as f64, v[2] as f64))
        .collect::<Vec<_>>();

    let faces = raw
        .faces
        .iter()
        .map(|f| {
            [
                f.vertices[0] as u32,
                f.vertices[1] as u32,
                f.vertices[2] as u32,
            ]
        })
        .collect::<Vec<_>>();

    MeshData3::new(points, faces)
}

/// Read an STL file into a `Mesh3`.
///
/// This is a temporary bridge for callers which need the accelerated type, and will be replaced by
/// a conversion from `MeshData3` once that exists. Prefer `load_stl_mesh_data`.
///
/// # Arguments
///
/// * `path`: the path to the STL file
/// * `merge_duplicates`: additionally merge points which compare equal as `f64`, which after the
///   exact weld the reader already performs means only merging `0.0` with `-0.0`
/// * `delete_degenerate`: drop triangles with zero area or bad topology
pub fn read_mesh_stl(
    path: &Path,
    merge_duplicates: bool,
    delete_degenerate: bool,
) -> Result<Mesh3> {
    let (points, faces, _) = load_stl_mesh_data(path)?.into_parts();
    Mesh3::new_with_options(
        points,
        faces,
        false,
        merge_duplicates,
        delete_degenerate,
        None,
    )
}

// ===============================================================================================
// Writing
// ===============================================================================================

/// Options for writing an STL file.
///
/// Marked `#[non_exhaustive]` and given a `Default`, so construct it as
/// `StlWriteOpts { allow_attribute_loss: true, ..Default::default() }` and new fields will not
/// break your code.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct StlWriteOpts {
    /// Accept that every attribute on the mesh will be discarded, since STL cannot carry any.
    ///
    /// Defaults to false, so a mesh carrying attributes is refused rather than silently stripped.
    /// A mesh with no attributes writes either way, so the safe setting costs nothing.
    pub allow_attribute_loss: bool,

    /// Write the binary encoding rather than ascii. Defaults to true: binary is about a fifth the
    /// size and does not round the coordinates a second time on the way to decimal.
    pub binary: bool,

    /// A comment to place in the file.
    ///
    /// In binary this fills the 80-byte header, truncated if longer and zero padded if shorter. In
    /// ascii it becomes the name on the `solid` line. Leading whitespace is trimmed and a leading
    /// `solid` is rejected, since a header beginning with that word makes naive format detectors
    /// read a binary file as ascii.
    pub header: Option<String>,
}

impl Default for StlWriteOpts {
    fn default() -> Self {
        Self {
            allow_attribute_loss: false,
            binary: true,
            header: None,
        }
    }
}

/// The fixed size of the binary STL header, which is a free-form comment field.
const BINARY_HEADER_LEN: usize = 80;

/// Write a mesh to an STL file, which carries geometry and nothing else.
///
/// Every attribute on the mesh is discarded, so this fails unless the mesh has none or the caller
/// set `allow_attribute_loss`. See the module documentation for the precision and point identity
/// the format also costs.
///
/// Unlike the previous implementation, a triangle whose normal cannot be computed is **not**
/// silently skipped: a degenerate face is written with a zero normal, which is what the format
/// says to do when the normal is not supplied, so the output always has as many facets as the mesh
/// has faces.
///
/// # Arguments
///
/// * `path`: the path to write to, which is overwritten if it already exists
/// * `mesh`: the mesh to write
/// * `opts`: encoding, header, and attribute loss options
///
/// returns: `Result<()>`
pub fn write_stl_mesh_data(path: &Path, mesh: &MeshData3, opts: &StlWriteOpts) -> Result<()> {
    let mut file = BufWriter::new(File::create(path)?);
    write_stl_to(&mut file, mesh, opts)?;
    file.flush()?;
    Ok(())
}

/// Write a mesh in the STL format to any sink. See `write_stl_mesh_data`.
pub fn write_stl_to<W: Write>(out: &mut W, mesh: &MeshData3, opts: &StlWriteOpts) -> Result<()> {
    mesh.check_attribute_loss("STL", opts.allow_attribute_loss)?;

    if mesh.face_count() > u32::MAX as usize {
        return Err(format!(
            "This mesh has {} faces, but the binary STL facet count is a 32 bit value",
            mesh.face_count()
        )
        .into());
    }

    let header = validated_header(opts.header.as_deref())?;

    if opts.binary {
        write_binary(out, mesh, header)
    } else {
        write_ascii(out, mesh, header)
    }
}

/// The name written when the caller supplies no header.
///
/// The ascii encoding cannot use an empty one: readers key on the literal `"solid "`, including
/// the trailing space, so a bare `solid` line is rejected as not being ascii at all. Rather than
/// have the two encodings disagree about what an absent header means, both use this.
const DEFAULT_HEADER: &str = "engeom";

/// Check a caller supplied header for the one thing that actually breaks files, and substitute the
/// default when there is nothing usable.
fn validated_header(header: Option<&str>) -> Result<&str> {
    let text = header.unwrap_or("").trim_start();

    // A binary file whose header starts with "solid" is read as ascii by detectors which sniff
    // only the first five bytes, and an ascii file gets a doubled keyword. Compare as bytes, since
    // the header may begin with a multi byte character and slicing it would panic.
    let head = &text.as_bytes()[..text.len().min(5)];
    if head.eq_ignore_ascii_case(b"solid") {
        return Err(
            "An STL header may not begin with the word 'solid', because it makes the \
                    file look like ascii to a naive format detector"
                .into(),
        );
    }

    if text.is_empty() {
        return Ok(DEFAULT_HEADER);
    }

    Ok(text)
}

/// Look up a face's three points, which construction guarantees are in range.
fn face_points(mesh: &MeshData3, face: &[u32; 3]) -> [Point3; 3] {
    [
        mesh.points()[face[0] as usize],
        mesh.points()[face[1] as usize],
        mesh.points()[face[2] as usize],
    ]
}

/// The facet normal to record, which is zero for a degenerate face. The specification allows a
/// zero normal to mean "derive it from the winding", which is the honest thing to write when the
/// face has no normal at all.
fn facet_normal(p: &[Point3; 3]) -> [f32; 3] {
    match compute_face_normal(p) {
        Some(n) => [n.x as f32, n.y as f32, n.z as f32],
        None => [0.0, 0.0, 0.0],
    }
}

fn write_binary<W: Write>(out: &mut W, mesh: &MeshData3, header: &str) -> Result<()> {
    // The header is a fixed size field, so truncate on a character boundary rather than mid
    // codepoint, and zero pad the rest.
    let mut header_bytes = [0u8; BINARY_HEADER_LEN];
    let source = header.as_bytes();
    let mut take = source.len().min(BINARY_HEADER_LEN);
    while take > 0 && !header.is_char_boundary(take) {
        take -= 1;
    }
    header_bytes[..take].copy_from_slice(&source[..take]);
    out.write_all(&header_bytes)?;

    out.write_all(&(mesh.face_count() as u32).to_le_bytes())?;

    for f in mesh.faces() {
        let p = face_points(mesh, f);
        let n = facet_normal(&p);

        for c in n.iter() {
            out.write_all(&c.to_le_bytes())?;
        }
        for v in p.iter() {
            out.write_all(&(v.x as f32).to_le_bytes())?;
            out.write_all(&(v.y as f32).to_le_bytes())?;
            out.write_all(&(v.z as f32).to_le_bytes())?;
        }

        // The attribute byte count, which the specification says is zero. See the module
        // documentation for why this is not used to carry a color.
        out.write_all(&0u16.to_le_bytes())?;
    }

    Ok(())
}

fn write_ascii<W: Write>(out: &mut W, mesh: &MeshData3, header: &str) -> Result<()> {
    // Values go out through `Display`, which for a float is the shortest decimal that reads back
    // as the same value. That is both compact and unambiguous, and it avoids the exponent
    // formatting that some readers of this format handle poorly.
    //
    // The name is never empty, so the first line always carries the "solid " prefix that readers
    // probe for.
    writeln!(out, "solid {}", header)?;

    for f in mesh.faces() {
        let p = face_points(mesh, f);
        let n = facet_normal(&p);

        writeln!(out, "facet normal {} {} {}", n[0], n[1], n[2])?;
        writeln!(out, "  outer loop")?;
        for v in p.iter() {
            // Narrowed to f32 to match what the binary encoding stores, so the two forms carry
            // exactly the same values rather than the ascii one carrying more.
            writeln!(
                out,
                "    vertex {} {} {}",
                v.x as f32, v.y as f32, v.z as f32
            )?;
        }
        writeln!(out, "  endloop")?;
        writeln!(out, "endfacet")?;
    }

    writeln!(out, "endsolid {}", header)?;
    Ok(())
}

/// Write a `Mesh3` to a binary STL file.
///
/// This is a temporary bridge for callers holding the accelerated type, and will be replaced by a
/// conversion to `MeshData3` once that exists. Prefer `write_stl_mesh_data`.
pub fn write_mesh_stl(path: &Path, mesh: &Mesh3) -> Result<()> {
    let data = MeshData3::new(mesh.vertices().to_vec(), mesh.faces().to_vec())?;
    write_stl_mesh_data(path, &data, &StlWriteOpts::default())
}

// ===============================================================================================
// Tests
// ===============================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::MeshAttr3;
    use crate::tests::get_test_file_path;
    use crate::{UnitVec3, Vector3};
    use std::io::Cursor;

    /// A mesh with distinct geometry in every axis, no attributes, and no degeneracy.
    fn simple() -> MeshData3 {
        MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 2.0, 0.0),
                Point3::new(0.0, 2.0, 0.5),
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
        .unwrap()
    }

    fn round_trip(mesh: &MeshData3, opts: &StlWriteOpts) -> Result<MeshData3> {
        let mut buf = Vec::new();
        write_stl_to(&mut buf, mesh, opts)?;
        read_stl_mesh_data(&mut Cursor::new(buf))
    }

    /// STL has no point identity, so a round trip renumbers the points. What has to survive is the
    /// set of triangles, each described by its three corner positions at f32 precision.
    fn triangle_set(mesh: &MeshData3) -> Vec<[[u32; 3]; 3]> {
        let mut out = mesh
            .faces()
            .iter()
            .map(|f| {
                let mut corners = face_points(mesh, f).map(|p| {
                    [
                        (p.x as f32).to_bits(),
                        (p.y as f32).to_bits(),
                        (p.z as f32).to_bits(),
                    ]
                });
                // Sorted so that a rotation of the winding does not read as a different triangle.
                corners.sort();
                corners
            })
            .collect::<Vec<_>>();
        out.sort();
        out
    }

    #[test]
    fn a_mesh_round_trips_through_the_binary_encoding() -> Result<()> {
        let mesh = simple();
        let back = round_trip(&mesh, &StlWriteOpts::default())?;

        assert_eq!(back.face_count(), mesh.face_count());
        assert_eq!(back.point_count(), mesh.point_count());
        assert_eq!(triangle_set(&back), triangle_set(&mesh));

        Ok(())
    }

    #[test]
    fn a_mesh_round_trips_through_the_ascii_encoding() -> Result<()> {
        let mesh = simple();
        let opts = StlWriteOpts {
            binary: false,
            ..Default::default()
        };
        let back = round_trip(&mesh, &opts)?;

        assert_eq!(back.face_count(), mesh.face_count());
        assert_eq!(triangle_set(&back), triangle_set(&mesh));

        Ok(())
    }

    /// The two encodings must describe exactly the same geometry, which is only true because the
    /// ascii writer narrows to f32 the same way the binary one does.
    #[test]
    fn the_two_encodings_agree() -> Result<()> {
        let mesh = load_ply_fixture()?;

        let binary = round_trip(
            &mesh,
            &StlWriteOpts {
                allow_attribute_loss: true,
                ..Default::default()
            },
        )?;
        let ascii = round_trip(
            &mesh,
            &StlWriteOpts {
                binary: false,
                allow_attribute_loss: true,
                ..Default::default()
            },
        )?;

        assert_eq!(binary.point_count(), ascii.point_count());
        assert_eq!(binary.face_count(), ascii.face_count());
        assert_eq!(triangle_set(&binary), triangle_set(&ascii));

        Ok(())
    }

    fn load_ply_fixture() -> Result<MeshData3> {
        crate::io::load_ply_mesh_data(&get_test_file_path("bun_zipper_res4.ply"))
    }

    /// The Stanford bunny at res4, which is a real mesh with 453 points and 948 faces, and carries
    /// two per-point scalars that STL cannot represent.
    #[test]
    fn a_real_mesh_round_trips_with_its_topology_intact() -> Result<()> {
        let mesh = load_ply_fixture()?;
        assert_eq!(mesh.point_count(), 453);
        assert_eq!(mesh.face_count(), 948);

        let back = round_trip(
            &mesh,
            &StlWriteOpts {
                allow_attribute_loss: true,
                ..Default::default()
            },
        )?;

        // The exact weld recovers the original point count, because every shared corner was
        // written from the same f64 value and so narrowed to the same f32.
        assert_eq!(back.point_count(), 453);
        assert_eq!(back.face_count(), 948);
        assert_eq!(triangle_set(&back), triangle_set(&mesh));

        // Nothing came back but geometry.
        assert!(back.attrs().is_empty());

        Ok(())
    }

    /// The guard exists so that measured data cannot be thrown away without the caller saying so.
    #[test]
    fn a_mesh_with_attributes_is_refused_by_default() {
        let mut mesh = simple();
        mesh.set_point_stdev(Some(vec![0.01, 0.02, 0.03, 0.04]))
            .unwrap();

        let mut buf = Vec::new();
        let err = write_stl_to(&mut buf, &mesh, &StlWriteOpts::default())
            .expect_err("writing a mesh with attributes must fail by default");

        let message = err.to_string();
        assert!(message.contains("point_stdev"), "message was: {}", message);
        assert!(message.contains("allow_attribute_loss"));

        // Nothing was written, so a refused write cannot leave a truncated file behind.
        assert!(buf.is_empty());
    }

    /// Every attribute in the way is named, so the caller can see the full cost of accepting.
    #[test]
    fn the_refusal_names_every_attribute_that_would_be_lost() {
        let mut mesh = simple();
        mesh.set_point_normals(Some(vec![UnitVec3::new_normalize(Vector3::z()); 4]))
            .unwrap();
        mesh.set_face_labels(Some(vec![7, 9])).unwrap();
        mesh.insert_point_attr("intensity", MeshAttr3::Scalar(vec![1.0, 2.0, 3.0, 4.0]))
            .unwrap();

        let mut buf = Vec::new();
        let message = write_stl_to(&mut buf, &mesh, &StlWriteOpts::default())
            .expect_err("must fail")
            .to_string();

        for name in ["point_normals", "face_labels", "intensity"] {
            assert!(message.contains(name), "{} missing from: {}", name, message);
        }
    }

    #[test]
    fn accepting_the_loss_writes_the_geometry() -> Result<()> {
        let mut mesh = simple();
        mesh.set_point_stdev(Some(vec![0.01, 0.02, 0.03, 0.04]))?;

        let back = round_trip(
            &mesh,
            &StlWriteOpts {
                allow_attribute_loss: true,
                ..Default::default()
            },
        )?;

        assert_eq!(triangle_set(&back), triangle_set(&mesh));
        assert!(back.point_stdev().is_none());

        Ok(())
    }

    /// A mesh with nothing to lose writes without the flag, so the safe default costs nothing in
    /// the common case.
    #[test]
    fn a_mesh_with_no_attributes_needs_no_flag() -> Result<()> {
        let mut buf = Vec::new();
        write_stl_to(&mut buf, &simple(), &StlWriteOpts::default())?;
        assert!(!buf.is_empty());

        Ok(())
    }

    /// The previous implementation skipped any face whose normal it could not compute, so the file
    /// could quietly have fewer facets than the mesh had faces.
    #[test]
    fn a_degenerate_face_is_written_rather_than_skipped() -> Result<()> {
        let mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(2.0, 0.0, 0.0), // collinear with the other two
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2], [0, 1, 3]],
        )?;

        let mut buf = Vec::new();
        write_stl_to(&mut buf, &mesh, &StlWriteOpts::default())?;

        // 80 byte header, a 32 bit facet count, and 50 bytes per facet.
        assert_eq!(buf.len(), 80 + 4 + 50 * 2);
        assert_eq!(u32::from_le_bytes(buf[80..84].try_into().unwrap()), 2);

        // The degenerate facet carries a zero normal, which the format reads as "derive it from
        // the winding" rather than as a claim about the geometry.
        let normal = &buf[84..96];
        assert!(normal.iter().all(|b| *b == 0));

        Ok(())
    }

    #[test]
    fn the_header_is_written_into_the_binary_file() -> Result<()> {
        let opts = StlWriteOpts {
            header: Some("engeom test part".to_string()),
            ..Default::default()
        };

        let mut buf = Vec::new();
        write_stl_to(&mut buf, &simple(), &opts)?;

        assert_eq!(&buf[..16], b"engeom test part");
        // The remainder of the fixed size field is zero padded.
        assert!(buf[16..80].iter().all(|b| *b == 0));

        Ok(())
    }

    /// With no header supplied the file still gets a name, because the ascii encoding cannot use
    /// an empty one: readers probe for the literal "solid " including its trailing space.
    #[test]
    fn a_default_header_is_written_when_none_is_given() -> Result<()> {
        let mut binary = Vec::new();
        write_stl_to(&mut binary, &simple(), &StlWriteOpts::default())?;
        assert_eq!(&binary[..6], b"engeom");

        let mut ascii = Vec::new();
        write_stl_to(
            &mut ascii,
            &simple(),
            &StlWriteOpts {
                binary: false,
                ..Default::default()
            },
        )?;
        assert!(String::from_utf8(ascii)?.starts_with("solid engeom\n"));

        Ok(())
    }

    #[test]
    fn the_header_becomes_the_solid_name_in_ascii() -> Result<()> {
        let opts = StlWriteOpts {
            binary: false,
            header: Some("engeom test part".to_string()),
            ..Default::default()
        };

        let mut buf = Vec::new();
        write_stl_to(&mut buf, &simple(), &opts)?;
        let text = String::from_utf8(buf)?;

        assert!(text.starts_with("solid engeom test part\n"));
        assert!(text.trim_end().ends_with("endsolid engeom test part"));

        Ok(())
    }

    /// A binary header beginning with "solid" makes format detectors that sniff the first five
    /// bytes read the file as ascii, so it is refused rather than written.
    #[test]
    fn a_header_beginning_with_solid_is_refused() {
        for header in ["solid part", "SOLID PART", "  solid part"] {
            let opts = StlWriteOpts {
                header: Some(header.to_string()),
                ..Default::default()
            };
            let mut buf = Vec::new();
            assert!(
                write_stl_to(&mut buf, &simple(), &opts).is_err(),
                "header {:?} should have been refused",
                header
            );
        }

        // A header which merely contains the word elsewhere is fine.
        let opts = StlWriteOpts {
            header: Some("a solid part".to_string()),
            ..Default::default()
        };
        let mut buf = Vec::new();
        assert!(write_stl_to(&mut buf, &simple(), &opts).is_ok());
    }

    /// A header longer than the fixed field is truncated, and must not be cut through the middle
    /// of a multi byte character.
    #[test]
    fn an_oversized_header_is_truncated_on_a_character_boundary() -> Result<()> {
        // One ascii character followed by 26 four byte ones is 105 bytes, and the 80 byte cut
        // lands inside the 20th emoji rather than between two of them.
        let opts = StlWriteOpts {
            header: Some(format!("x{}", "😀".repeat(26))),
            ..Default::default()
        };

        let mut buf = Vec::new();
        write_stl_to(&mut buf, &simple(), &opts)?;

        assert_eq!(buf.len(), 80 + 4 + 50 * 2);
        let header = std::str::from_utf8(&buf[..80]).expect("the header must remain valid utf8");
        let text = header.trim_end_matches('\0');

        // 1 + 19 whole emoji is 77 bytes; the 20th would have run past the field.
        assert_eq!(text.chars().count(), 20);
        assert_eq!(text.len(), 77);

        Ok(())
    }

    /// Reading rebuilds the point buffer by welding corners with identical coordinates, which is
    /// what recovers a topology the format does not store.
    #[test]
    fn reading_welds_the_triangle_soup() -> Result<()> {
        // Two triangles sharing an edge: six corners in the file, four points after welding.
        let mesh = simple();
        assert_eq!(mesh.point_count(), 4);

        let back = round_trip(&mesh, &StlWriteOpts::default())?;
        assert_eq!(back.point_count(), 4);
        assert_eq!(back.face_count(), 2);

        Ok(())
    }

    /// Coordinates are stored as 32 bit floats, so a value that needs more than 24 bits of mantissa
    /// does not survive. This is the format's cost, pinned so it is not mistaken for a bug.
    #[test]
    fn coordinates_are_quantized_to_f32() -> Result<()> {
        // The f32 step at 1000 mm is 2^-13, about 61 nanometres, so a value offset half a step
        // from a representable one is as wrong as this format can make it.
        let exact = 1000.000_03_f64;
        let mesh = MeshData3::new(
            vec![
                Point3::new(exact, 0.0, 0.0),
                Point3::new(exact + 1.0, 0.0, 0.0),
                Point3::new(exact, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )?;

        let back = round_trip(&mesh, &StlWriteOpts::default())?;

        let recovered = back.points().iter().map(|p| p.x).fold(f64::MAX, f64::min);
        assert_ne!(recovered, exact);
        assert_eq!(recovered, exact as f32 as f64);

        // Half of the f32 step at this distance, roughly 30 nanometres. Small in isolation, but
        // it grows linearly with distance from the origin: at 10 m it is about half a micron, and
        // a part held in machine coordinates can easily be that far out.
        let step = 1000.0_f32.to_bits() + 1;
        let step = f32::from_bits(step) as f64 - 1000.0;
        assert!((recovered - exact).abs() > step * 0.4);
        assert!((recovered - exact).abs() <= step * 0.5);

        Ok(())
    }
}
