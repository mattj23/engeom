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
use std::io::{BufWriter, Read, Write};
use std::path::Path;

const MAGIC: &[u8; 4] = b"BMSH";

/// Save a [`Mesh`] to a file in the portable binary mesh format.
///
/// Vertices are written as 32-bit floats and face indices as 32-bit unsigned integers,
/// all in little-endian byte order.
///
/// # Errors
///
/// Returns an error if the file cannot be created or written.
pub fn write_mesh_binary_file(path: &Path, mesh: &Mesh) -> Result<()> {
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    w.write_all(MAGIC)?;

    let vertices = mesh.vertices();
    w.write_all(&(vertices.len() as u32).to_le_bytes())?;
    for v in vertices {
        w.write_all(&(v.x as f32).to_le_bytes())?;
        w.write_all(&(v.y as f32).to_le_bytes())?;
        w.write_all(&(v.z as f32).to_le_bytes())?;
    }

    let faces = mesh.faces();
    w.write_all(&(faces.len() as u32).to_le_bytes())?;
    for face in faces {
        w.write_all(&face[0].to_le_bytes())?;
        w.write_all(&face[1].to_le_bytes())?;
        w.write_all(&face[2].to_le_bytes())?;
    }

    Ok(())
}

/// Load a [`Mesh`] from a file written by [`write_mesh_binary_file`].
///
/// Vertex coordinates are read as 32-bit floats and widened to `f64`. Face indices are
/// read as 32-bit unsigned integers.
///
/// # Errors
///
/// Returns an error if the file cannot be read, or if the magic bytes do not match.
pub fn read_mesh_binary_file(path: &Path) -> Result<Mesh> {
    let mut bytes = Vec::new();
    File::open(path)?.read_to_end(&mut bytes)?;

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
}