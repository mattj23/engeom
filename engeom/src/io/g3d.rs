//! Reader for GOM's `.g3d` binary file format ("GOM 3D File Format", version 1.1b, July 2010),
//! used by the Atos structured-light scanner and its inspection software (GOM/Zeiss Inspect) to
//! carry 3D measurement data.
//!
//! The format can in principle carry several kinds of view: triangle meshes, several flavors of
//! point cloud, scan sections, and colored triangle meshes. In practice GOM's own software only
//! ever wrote triangle meshes, even at the time the spec was published ("The current software...
//! uses the triangle meshes only. The other data types will not be used anymore."), so this
//! reader only understands the triangle mesh view (type 0) and skips any other view it finds
//! while walking the view chain.
//!
//! It also only supports the "all views in one file" layout, which is both what real files use
//! and the only layout the spec itself recommends ("We advise you not to use every view in a
//! separate file... for creating new files.").
//!
//! Every data block in the file is reached through an absolute file offset stored in a
//! previously read block, and every block carries its own declared size. Real files use this to
//! add fields beyond what the 2010 specification documents: the sample files this reader was
//! developed against declare triangle-mesh view headers of 432 bytes, against the 168 the
//! specification describes. This reader follows the specification's own advice for that case:
//! blocks are located with their stored pointer rather than assumed to be adjacent, and any
//! bytes in a point or triangle record beyond the fields read here are skipped rather than
//! interpreted.

use crate::{MeshData3, Point3, Result};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::Path;

const GLOBAL_MAGIC: &[u8; 8] = b"%GOM-3DH";

/// A triangle mesh view, the only kind of view this reader interprets.
const VIEW_TYPE_TRIANGLE_MESH: u32 = 0;

/// The size in bytes of the x, y, z coordinates this reader takes from a point record. A point
/// record's declared size may be larger (the specification's own default appends a trailing
/// `f32` quality value), in which case the remainder is skipped.
const POINT_XYZ_SIZE: u32 = 24;

/// The size in bytes of the three point indices this reader takes from a triangle record.
const TRIANGLE_SIZE: u32 = 12;

/// The minimum triangle-mesh view header size that gives us the point and triangle block
/// pointers this reader depends on.
const MIN_TRIANGLE_VIEW_HEADER_SIZE: u32 = 168;

/// Load a triangle mesh from a GOM `.g3d` file.
///
/// Only the triangle mesh view is read; any other view (point cloud, section, colored mesh,
/// etc.) present in the file is skipped. See the module documentation for the format and what
/// this reader does and does not support.
///
/// # Arguments
///
/// * `path`: the path to the `.g3d` file
///
/// returns: `Result<MeshData3>`
pub fn load_g3d_mesh_data(path: &Path) -> Result<MeshData3> {
    let mut file = BufReader::new(File::open(path)?);
    read_g3d_mesh_data(&mut file)
}

/// Read a triangle mesh from any seekable source of `.g3d` data. See [`load_g3d_mesh_data`].
pub fn read_g3d_mesh_data<R: Read + Seek>(source: &mut R) -> Result<MeshData3> {
    source.seek(SeekFrom::Start(0))?;

    let mut magic = [0u8; 8];
    source.read_exact(&mut magic)?;
    if &magic != GLOBAL_MAGIC {
        return Err("not a GOM g3d file: the global header magic bytes do not match".into());
    }

    // The byte order marker is a literal two byte tag, not a value to be decoded with a chosen
    // endianness: 0x0001 (bytes 01 00) means the rest of the file is little-endian, 0x0100
    // (bytes 00 01) means big-endian.
    let mut order = [0u8; 2];
    source.read_exact(&mut order)?;
    let little_endian = match order {
        [0x01, 0x00] => true,
        [0x00, 0x01] => false,
        other => {
            return Err(format!("g3d file has an unrecognized byte order marker {other:?}").into());
        }
    };

    let mut r = G3dReader {
        source,
        little_endian,
    };

    // Version is explicitly documented as untrustworthy for classifying file contents, so it is
    // read only to advance the cursor.
    let _version = r.u16()?;
    let _header_size = r.u32()?;
    // The project flag distinguishes files that belong to a multi-file project from ones that
    // don't; it says nothing about how to parse the file we were handed directly, so it is
    // ignored.
    let _project_flag = r.u32()?;
    let all_in_one_file = r.u32()?;
    let view_count = r.u32()?;
    let first_view = r.u32()?;
    // The remaining 64 bytes of the global header are a free-form comment, unused here.

    if all_in_one_file == 0 {
        return Err(
            "g3d files that split their views across separate files are not supported".into(),
        );
    }
    if view_count == 0 || first_view == 0 {
        return Err("g3d file declares no views".into());
    }

    let mut ptr = first_view;
    loop {
        r.seek_to(ptr)?;
        let next_view = r.u32()?;
        let header_size = r.u32()?;
        let _id = r.u32()?;
        let view_type = r.u32()?;

        if view_type == VIEW_TYPE_TRIANGLE_MESH {
            return read_triangle_mesh_view(&mut r, header_size);
        }

        if next_view == 0 {
            return Err("g3d file has no triangle mesh view".into());
        }
        ptr = next_view;
    }
}

/// Read the body of a triangle mesh view, with the reader positioned just after the `type`
/// field. `header_size` is the view's own declared header size, checked only to make sure the
/// fields this reader depends on are actually present.
fn read_triangle_mesh_view<R: Read + Seek>(
    r: &mut G3dReader<R>,
    header_size: u32,
) -> Result<MeshData3> {
    if header_size < MIN_TRIANGLE_VIEW_HEADER_SIZE {
        return Err(format!(
            "g3d triangle mesh view header is {header_size} bytes, but at least \
             {MIN_TRIANGLE_VIEW_HEADER_SIZE} are needed to hold the point and triangle block \
             pointers"
        )
        .into());
    }

    // name (64 bytes) + comment (64 bytes), neither of which this reader uses.
    r.skip(128)?;

    let num_points = r.u32()?;
    let ptr_points = r.u32()?;
    let size_point = r.u32()?;
    let num_triangles = r.u32()?;
    let ptr_triangles = r.u32()?;
    let size_triangle = r.u32()?;

    if size_point < POINT_XYZ_SIZE {
        return Err(format!(
            "g3d point records are {size_point} bytes, too small to hold an x, y, z coordinate"
        )
        .into());
    }
    if size_triangle < TRIANGLE_SIZE {
        return Err(format!(
            "g3d triangle records are {size_triangle} bytes, too small to hold three point \
             indices"
        )
        .into());
    }

    let mut points = Vec::with_capacity(num_points as usize);
    if num_points > 0 {
        r.seek_to(ptr_points)?;
        let padding = (size_point - POINT_XYZ_SIZE) as i64;
        for _ in 0..num_points {
            let x = r.f64()?;
            let y = r.f64()?;
            let z = r.f64()?;
            points.push(Point3::new(x, y, z));
            r.skip(padding)?;
        }
    }

    let mut faces = Vec::with_capacity(num_triangles as usize);
    if num_triangles > 0 {
        r.seek_to(ptr_triangles)?;
        let padding = (size_triangle - TRIANGLE_SIZE) as i64;
        for _ in 0..num_triangles {
            let a = r.u32()?;
            let b = r.u32()?;
            let c = r.u32()?;
            faces.push([a, b, c]);
            r.skip(padding)?;
        }
    }

    MeshData3::new(points, faces)
}

// ================================================================================================
// Byte reader
// ================================================================================================

/// A cursor over a seekable source, tracking the endianness declared by the file's global
/// header byte order marker.
struct G3dReader<'a, R: Read + Seek> {
    source: &'a mut R,
    little_endian: bool,
}

impl<R: Read + Seek> G3dReader<'_, R> {
    fn seek_to(&mut self, pos: u32) -> Result<()> {
        self.source.seek(SeekFrom::Start(pos as u64))?;
        Ok(())
    }

    fn skip(&mut self, n: i64) -> Result<()> {
        if n != 0 {
            self.source.seek(SeekFrom::Current(n))?;
        }
        Ok(())
    }

    fn bytes<const N: usize>(&mut self) -> Result<[u8; N]> {
        let mut buf = [0u8; N];
        self.source.read_exact(&mut buf)?;
        Ok(buf)
    }

    fn u16(&mut self) -> Result<u16> {
        let b = self.bytes::<2>()?;
        Ok(if self.little_endian {
            u16::from_le_bytes(b)
        } else {
            u16::from_be_bytes(b)
        })
    }

    fn u32(&mut self) -> Result<u32> {
        let b = self.bytes::<4>()?;
        Ok(if self.little_endian {
            u32::from_le_bytes(b)
        } else {
            u32::from_be_bytes(b)
        })
    }

    fn f64(&mut self) -> Result<f64> {
        let b = self.bytes::<8>()?;
        Ok(if self.little_endian {
            f64::from_le_bytes(b)
        } else {
            f64::from_be_bytes(b)
        })
    }
}

// ================================================================================================
// Tests
// ================================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::get_test_file_path;
    use approx::assert_relative_eq;
    use std::io::Cursor;

    #[test]
    fn reads_the_stud_bolt_file() -> Result<()> {
        let mesh = load_g3d_mesh_data(&get_test_file_path("stud-bolt.g3d"))?;

        assert_eq!(mesh.point_count(), 8565);
        assert_eq!(mesh.face_count(), 16957);

        let first = mesh.points()[0];
        assert_relative_eq!(first.x, 3.8822855365478506, epsilon = 1.0e-9);
        assert_relative_eq!(first.y, -16.53775066421784, epsilon = 1.0e-9);
        assert_relative_eq!(first.z, -22.36773866300889, epsilon = 1.0e-9);

        let last = mesh.points()[mesh.point_count() - 1];
        assert_relative_eq!(last.x, 51.102225220201504, epsilon = 1.0e-9);
        assert_relative_eq!(last.y, -17.207949890147148, epsilon = 1.0e-9);
        assert_relative_eq!(last.z, -27.249870490302808, epsilon = 1.0e-9);

        assert_eq!(mesh.faces()[0], [0, 44, 23]);
        assert_eq!(mesh.faces()[mesh.face_count() - 1], [8564, 8561, 8563]);

        Ok(())
    }

    /// Build a minimal, single-triangle g3d file by hand: a global header pointing at one
    /// triangle mesh view, whose point and triangle records are exactly the documented size.
    ///
    /// Pointer fields (`ptr_points`, `ptr_triangles`, and view chain pointers) are written as
    /// placeholders and backpatched once the position they refer to is known, which is far less
    /// error prone than pre-computing block sizes by hand.
    ///
    /// `point_padding` and `triangle_padding` let a test declare `size_point`/`size_triangle`
    /// larger than the fields actually written, to exercise the forward-compatibility skip
    /// logic. `view_header_padding` inflates the declared view header size the same way, without
    /// writing any corresponding bytes: since this reader locates blocks purely through their
    /// pointers, a header size that overstates the physical layout must not throw parsing off.
    struct Builder {
        little_endian: bool,
        point_padding: u32,
        triangle_padding: u32,
        view_header_padding: u32,
        second_view_type: Option<u32>,
    }

    impl Default for Builder {
        fn default() -> Self {
            Self {
                little_endian: true,
                point_padding: 0,
                triangle_padding: 0,
                view_header_padding: 0,
                second_view_type: None,
            }
        }
    }

    impl Builder {
        fn push_u16(&self, buf: &mut Vec<u8>, v: u16) {
            if self.little_endian {
                buf.extend_from_slice(&v.to_le_bytes());
            } else {
                buf.extend_from_slice(&v.to_be_bytes());
            }
        }

        fn push_u32(&self, buf: &mut Vec<u8>, v: u32) {
            if self.little_endian {
                buf.extend_from_slice(&v.to_le_bytes());
            } else {
                buf.extend_from_slice(&v.to_be_bytes());
            }
        }

        fn push_f64(&self, buf: &mut Vec<u8>, v: f64) {
            if self.little_endian {
                buf.extend_from_slice(&v.to_le_bytes());
            } else {
                buf.extend_from_slice(&v.to_be_bytes());
            }
        }

        /// Overwrite a previously written `u32` placeholder, at the offset returned when it was
        /// pushed, with its final value.
        fn patch_u32(&self, buf: &mut [u8], offset: usize, v: u32) {
            let bytes = if self.little_endian {
                v.to_le_bytes()
            } else {
                v.to_be_bytes()
            };
            buf[offset..offset + 4].copy_from_slice(&bytes);
        }

        /// Build a file with a triangle covering points (0,0,0), (1,0,0), (0,1,0), optionally
        /// preceded by an extra, unrecognized view type.
        fn build(&self) -> Vec<u8> {
            let points = [
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ];

            let mut buf = Vec::new();

            // ---- Global header ----
            buf.extend_from_slice(GLOBAL_MAGIC);
            if self.little_endian {
                buf.extend_from_slice(&[0x01, 0x00]);
            } else {
                buf.extend_from_slice(&[0x00, 0x01]);
            }
            self.push_u16(&mut buf, 100); // version
            self.push_u32(&mut buf, 96); // header size
            self.push_u32(&mut buf, 0); // project flag
            self.push_u32(&mut buf, 1); // all views in one file
            self.push_u32(
                &mut buf,
                if self.second_view_type.is_some() {
                    2
                } else {
                    1
                },
            );
            self.push_u32(&mut buf, 96); // pointer to first view (always right after this header)
            buf.extend_from_slice(&[0u8; 64]); // comment
            assert_eq!(buf.len(), 96);

            // ---- Optional leading view of an unrecognized type, with no data of its own ----
            if let Some(view_type) = self.second_view_type {
                let next_view_ptr_offset = buf.len();
                self.push_u32(&mut buf, 0); // pointer to next view, patched below
                self.push_u32(&mut buf, 168); // header size
                self.push_u32(&mut buf, 1); // id
                self.push_u32(&mut buf, view_type); // type
                buf.extend_from_slice(&[0u8; 64]); // name
                buf.extend_from_slice(&[0u8; 64]); // comment
                self.push_u32(&mut buf, 0); // num points
                self.push_u32(&mut buf, 0); // ptr points
                self.push_u32(&mut buf, 28); // size point
                self.push_u32(&mut buf, 0); // num triangles
                self.push_u32(&mut buf, 0); // ptr triangles
                self.push_u32(&mut buf, 12); // size triangle

                let mesh_view_ptr = buf.len() as u32;
                self.patch_u32(&mut buf, next_view_ptr_offset, mesh_view_ptr);
            }

            // ---- Triangle mesh view header ----
            self.push_u32(&mut buf, 0); // pointer to next view: this is the last one
            self.push_u32(&mut buf, 168 + self.view_header_padding); // declared header size
            self.push_u32(&mut buf, 2); // id
            self.push_u32(&mut buf, VIEW_TYPE_TRIANGLE_MESH); // type
            buf.extend_from_slice(&[0u8; 64]); // name
            buf.extend_from_slice(&[0u8; 64]); // comment

            // Only x, y, z are actually written per point below (unlike a real GOM file, which
            // appends a trailing quality float), so the declared stride must match that, plus
            // whatever padding this test asked for.
            let size_point = POINT_XYZ_SIZE + self.point_padding;
            let size_triangle = TRIANGLE_SIZE + self.triangle_padding;

            self.push_u32(&mut buf, points.len() as u32); // num points
            let ptr_points_offset = buf.len();
            self.push_u32(&mut buf, 0); // ptr points, patched below
            self.push_u32(&mut buf, size_point);
            self.push_u32(&mut buf, 1); // num triangles
            let ptr_triangles_offset = buf.len();
            self.push_u32(&mut buf, 0); // ptr triangles, patched below
            self.push_u32(&mut buf, size_triangle);

            // ---- Point block ----
            let ptr_points = buf.len() as u32;
            self.patch_u32(&mut buf, ptr_points_offset, ptr_points);
            for p in &points {
                self.push_f64(&mut buf, p.x);
                self.push_f64(&mut buf, p.y);
                self.push_f64(&mut buf, p.z);
                buf.extend_from_slice(&vec![0u8; self.point_padding as usize]);
            }

            // ---- Triangle block ----
            let ptr_triangles = buf.len() as u32;
            self.patch_u32(&mut buf, ptr_triangles_offset, ptr_triangles);
            self.push_u32(&mut buf, 0);
            self.push_u32(&mut buf, 1);
            self.push_u32(&mut buf, 2);
            buf.extend_from_slice(&vec![0u8; self.triangle_padding as usize]);

            buf
        }
    }

    fn assert_is_the_test_triangle(mesh: &MeshData3) {
        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 1);
        assert_eq!(mesh.faces()[0], [0, 1, 2]);
        assert_relative_eq!(mesh.points()[1].x, 1.0, epsilon = 0.0);
        assert_relative_eq!(mesh.points()[2].y, 1.0, epsilon = 0.0);
    }

    #[test]
    fn a_minimal_file_round_trips() -> Result<()> {
        let bytes = Builder::default().build();
        let mesh = read_g3d_mesh_data(&mut Cursor::new(bytes))?;
        assert_is_the_test_triangle(&mesh);
        Ok(())
    }

    #[test]
    fn a_big_endian_file_is_read_correctly() -> Result<()> {
        let bytes = Builder {
            little_endian: false,
            ..Default::default()
        }
        .build();
        let mesh = read_g3d_mesh_data(&mut Cursor::new(bytes))?;
        assert_is_the_test_triangle(&mesh);
        Ok(())
    }

    /// Blocks declared larger than this reader's known fields (a newer file version, per the
    /// specification's own forward-compatibility rule) must have their extra bytes skipped
    /// rather than misread as the next field.
    #[test]
    fn oversized_blocks_are_tolerated() -> Result<()> {
        let bytes = Builder {
            point_padding: 8,
            triangle_padding: 4,
            view_header_padding: 16,
            ..Default::default()
        }
        .build();
        let mesh = read_g3d_mesh_data(&mut Cursor::new(bytes))?;
        assert_is_the_test_triangle(&mesh);
        Ok(())
    }

    /// A view type this reader does not recognize must be skipped in favor of a later triangle
    /// mesh view, rather than causing a failure.
    #[test]
    fn an_unrecognized_view_is_skipped() -> Result<()> {
        let bytes = Builder {
            second_view_type: Some(3), // unsorted point cloud
            ..Default::default()
        }
        .build();
        let mesh = read_g3d_mesh_data(&mut Cursor::new(bytes))?;
        assert_is_the_test_triangle(&mesh);
        Ok(())
    }

    #[test]
    fn an_invalid_magic_is_rejected() {
        let mut bytes = Builder::default().build();
        bytes[0] = b'X';
        assert!(read_g3d_mesh_data(&mut Cursor::new(bytes)).is_err());
    }

    #[test]
    fn a_file_with_no_triangle_mesh_view_is_rejected() {
        let mut buf = Vec::new();
        buf.extend_from_slice(GLOBAL_MAGIC);
        buf.extend_from_slice(&[0x01, 0x00]);
        buf.extend_from_slice(&100u16.to_le_bytes());
        buf.extend_from_slice(&96u32.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(&1u32.to_le_bytes()); // all in one file
        buf.extend_from_slice(&1u32.to_le_bytes()); // view count
        buf.extend_from_slice(&96u32.to_le_bytes()); // first view pointer
        buf.extend_from_slice(&[0u8; 64]);

        // A lone rastered point cloud view (type 1), which this reader has no interest in.
        buf.extend_from_slice(&0u32.to_le_bytes()); // next view: none
        buf.extend_from_slice(&212u32.to_le_bytes()); // header size
        buf.extend_from_slice(&1u32.to_le_bytes()); // id
        buf.extend_from_slice(&1u32.to_le_bytes()); // type: rastered point cloud

        assert!(read_g3d_mesh_data(&mut Cursor::new(buf)).is_err());
    }

    #[test]
    fn a_split_file_layout_is_rejected() {
        let mut buf = Vec::new();
        buf.extend_from_slice(GLOBAL_MAGIC);
        buf.extend_from_slice(&[0x01, 0x00]);
        buf.extend_from_slice(&100u16.to_le_bytes());
        buf.extend_from_slice(&96u32.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes()); // all in one file = 0 (split layout)
        buf.extend_from_slice(&1u32.to_le_bytes());
        buf.extend_from_slice(&96u32.to_le_bytes());
        buf.extend_from_slice(&[0u8; 64]);

        assert!(read_g3d_mesh_data(&mut Cursor::new(buf)).is_err());
    }
}
