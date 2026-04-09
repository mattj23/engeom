//! This module contains a no-nonsense, extremely simple binary mesh format which serializes only
//! the vertices and face indices of a mesh. The vertices are stored as 32-bit floats and the face
//! indices as 32-bit unsigned integers.
//!
//! There is no default extension for this format, but `.binmsh` is one suggestion.
//!
//! Binary layout (all values little-endian):
//! - 4 bytes magic: `b"BMSH"`
//! - 1 × `u32`: vertex count
//! - vertex count × 3 × `f32`: x, y, z per vertex
//! - 1 × `u32`: face count
//! - face count × 3 × `u32`: vertex indices per face

use crate::{Mesh, Point3, Result};
use std::fs::File;
use std::io::{BufWriter, Cursor, Read, Write};
use std::path::Path;

const MAGIC: &[u8; 4] = b"BMSH";

/// Save a [`Mesh`] to a file in the portable binary mesh format.
///
/// This is a no-nonsense, extremely simple binary mesh format that serializes only the vertices
/// and face indices of a mesh. The vertices are stored as 32-bit floats and the face indices as
/// 32-bit unsigned integers.
///
/// There is no default extension for this format, but `.binmsh` is one suggestion.
///
/// # Errors
///
/// Returns an error if the file cannot be created or written.
pub fn write_mesh_binary_file(path: &Path, mesh: &Mesh) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    write_mesh_binary(&mut w, mesh)
}

/// Load a [`Mesh`] from a file written by [`write_mesh_binary_file`].
///
/// This is a no-nonsense, extremely simple binary mesh format that serializes only the vertices
/// and face indices of a mesh. The vertices are stored as 32-bit floats and the face indices as
/// 32-bit unsigned integers.
///
/// There is no default extension for this format, but `.binmsh` is one suggestion.
///
/// # Errors
///
/// Returns an error if the file cannot be read, or if the magic bytes do not match.
pub fn read_mesh_binary_file(path: &Path) -> Result<Mesh> {
    let mut file = File::open(path)?;
    read_mesh_binary(&mut file)
}

/// Serialize a [`Mesh`] to a `Vec<u8>`.
///
/// This is a no-nonsense, extremely simple binary mesh format that serializes only the vertices
/// and face indices of a mesh. The vertices are stored as 32-bit floats and the face indices as
/// 32-bit unsigned integers.
///
/// # Errors
///
/// Returns an error if serialization fails (in practice this should never happen for in-memory
/// writes).
pub fn mesh_to_binary_bytes(mesh: &Mesh) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    write_mesh_binary(&mut buf, mesh)?;
    Ok(buf)
}

/// Deserialize a [`Mesh`] from a byte slice.
///
/// This is a no-nonsense, extremely simple binary mesh format that serializes only the vertices
/// and face indices of a mesh. The vertices are stored as 32-bit floats and the face indices as
/// 32-bit unsigned integers.
///
/// # Errors
///
/// Returns an error if the bytes are too short, malformed, or the magic bytes do not match.
pub fn mesh_from_binary_bytes(bytes: &[u8]) -> Result<Mesh> {
    let mut cursor = Cursor::new(bytes);
    read_mesh_binary(&mut cursor)
}

/// Serialize a [`Mesh`] into any [`Write`] destination.
///
/// This is a no-nonsense, extremely simple binary mesh format that serializes only the vertices
/// and face indices of a mesh. The vertices are stored as 32-bit floats and the face indices as
/// 32-bit unsigned integers.
///
/// # Errors
///
/// Returns an error if any write fails.
pub fn write_mesh_binary<W: Write>(writer: &mut W, mesh: &Mesh) -> Result<()> {
    writer.write_all(MAGIC)?;

    let vertices = mesh.vertices();
    writer.write_all(&(vertices.len() as u32).to_le_bytes())?;
    for v in vertices {
        writer.write_all(&(v.x as f32).to_le_bytes())?;
        writer.write_all(&(v.y as f32).to_le_bytes())?;
        writer.write_all(&(v.z as f32).to_le_bytes())?;
    }

    let faces = mesh.faces();
    writer.write_all(&(faces.len() as u32).to_le_bytes())?;
    for face in faces {
        writer.write_all(&face[0].to_le_bytes())?;
        writer.write_all(&face[1].to_le_bytes())?;
        writer.write_all(&face[2].to_le_bytes())?;
    }

    Ok(())
}

/// Deserialize a [`Mesh`] from any [`Read`] source.
///
/// This is a no-nonsense, extremely simple binary mesh format that serializes only the vertices
/// and face indices of a mesh. The vertices are stored as 32-bit floats and the face indices as
/// 32-bit unsigned integers.
///
/// # Errors
///
/// Returns an error if any read fails or if the magic bytes do not match.
pub fn read_mesh_binary<R: Read>(reader: &mut R) -> Result<Mesh> {
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes)?;

    let mut r = ByteReader::new(&bytes);

    let magic = r.read_bytes::<4>();
    if &magic != MAGIC {
        return Err("Not a binary mesh file: invalid magic bytes".into());
    }

    let vertex_count = r.read_u32() as usize;
    let mut vertices = Vec::with_capacity(vertex_count);
    for _ in 0..vertex_count {
        let x = f64::from(r.read_f32());
        let y = f64::from(r.read_f32());
        let z = f64::from(r.read_f32());
        vertices.push(Point3::new(x, y, z));
    }

    let face_count = r.read_u32() as usize;
    let mut faces = Vec::with_capacity(face_count);
    for _ in 0..face_count {
        faces.push([r.read_u32(), r.read_u32(), r.read_u32()]);
    }

    Ok(Mesh::new(vertices, faces, false))
}

// ------------------------------------------------------------------------------------------------
// Internal byte reader
// ------------------------------------------------------------------------------------------------

struct ByteReader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ByteReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn read_bytes<const N: usize>(&mut self) -> [u8; N] {
        let value: [u8; N] = self.bytes[self.offset..self.offset + N].try_into().unwrap();
        self.offset += N;
        value
    }

    fn read_f32(&mut self) -> f32 {
        f32::from_le_bytes(self.read_bytes::<4>())
    }

    fn read_u32(&mut self) -> u32 {
        u32::from_le_bytes(self.read_bytes::<4>())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::stanford_bun_2;
    use approx::assert_relative_eq;

    fn check_round_trip(mesh: &Mesh, recovered: &Mesh) {
        assert_eq!(mesh.vertices().len(), recovered.vertices().len());
        assert_eq!(mesh.faces().len(), recovered.faces().len());
        for (a, b) in mesh.vertices().iter().zip(recovered.vertices().iter()) {
            assert_relative_eq!(a, b, epsilon = 1e-8);
        }
        for (a, b) in mesh.faces().iter().zip(recovered.faces().iter()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn round_trip_bytes() {
        let mesh = stanford_bun_2();
        let bytes = mesh_to_binary_bytes(&mesh).unwrap();
        let recovered = mesh_from_binary_bytes(&bytes).unwrap();
        check_round_trip(&mesh, &recovered);
    }

    #[test]
    fn round_trip_file() {
        let mesh = stanford_bun_2();
        let path = std::env::temp_dir().join("stanford_bun_2_round_trip.binmsh");
        write_mesh_binary_file(&path, &mesh).unwrap();
        let recovered = read_mesh_binary_file(&path).unwrap();
        check_round_trip(&mesh, &recovered);
    }

    #[test]
    fn invalid_magic_rejected() {
        let bad = b"NOPE\x01\x00\x00\x00";
        assert!(mesh_from_binary_bytes(bad).is_err());
    }
}
