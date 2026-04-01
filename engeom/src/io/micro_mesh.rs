//! Micro mesh format, used to store meshes in internal binary data.  Works for meshes that have
//! less than u16::MAX vertices and discretizes the positions to 1/u16::MAX increments in an
//! axis-aligned bounding box of a size specified.
use crate::geom3::Aabb3;
use crate::{Mesh, Point3, Result};

/// Serialize a [`Mesh`] into the micro mesh binary format.
///
/// This is a convenience wrapper around [`u_mesh_data_to_bytes`] that operates directly on a
/// [`Mesh`] object. The mesh must have fewer than [`u16::MAX`] vertices.
///
/// # Errors
///
/// Returns an error if the mesh has more than [`u16::MAX`] vertices.
pub fn u_mesh_to_bytes(mesh: &Mesh) -> Result<Vec<u8>> {
    u_mesh_data_to_bytes(mesh.vertices(), mesh.faces())
}

/// Deserialize a [`Mesh`] from a micro mesh binary buffer.
///
/// This is a convenience wrapper around [`u_bytes_to_mesh_data`] that constructs and returns a
/// [`Mesh`] directly.
///
/// # Errors
///
/// Returns an error if the byte buffer is malformed or too short to be read.
pub fn u_bytes_to_mesh(bytes: &[u8]) -> Result<Mesh> {
    let (vertices, triangles) = u_bytes_to_mesh_data(bytes)?;
    Ok(Mesh::new(vertices, triangles, false))
}

/// Deserialize raw vertex and triangle data from a micro mesh binary buffer.
///
/// Reads an axis-aligned bounding box from the buffer header, then decodes vertex positions stored
/// as [`u16`] triples back into [`f64`] coordinates by interpolating within that bounding box.
/// Triangle indices are stored as [`u16`] values and widened to [`u32`] on read.
///
/// Returns a tuple of `(vertices, triangles)` where each triangle is an array of three vertex
/// indices.
///
/// # Errors
///
/// Returns an error if the byte buffer is malformed or too short to be read.
pub fn u_bytes_to_mesh_data(bytes: &[u8]) -> Result<(Vec<Point3>, Vec<[u32; 3]>)> {
    let mut reader = ByteRead::new(bytes);

    // Read the bounding box
    let min = Point3::new(reader.read_f64(), reader.read_f64(), reader.read_f64());
    let max = Point3::new(reader.read_f64(), reader.read_f64(), reader.read_f64());

    // Read the number of vertices
    let vertex_count = reader.read_u16() as usize;

    // Read the vertices
    let mut vertices = Vec::with_capacity(vertex_count);
    for _ in 0..vertex_count {
        let x = min.x + (f64::from(reader.read_u16()) / u16::MAX as f64) * (max.x - min.x);
        let y = min.y + (f64::from(reader.read_u16()) / u16::MAX as f64) * (max.y - min.y);
        let z = min.z + (f64::from(reader.read_u16()) / u16::MAX as f64) * (max.z - min.z);
        vertices.push(Point3::new(x, y, z));
    }

    // Read the number of triangles
    let triangle_count = reader.read_u32() as usize;

    // Read the triangles
    let mut triangles = Vec::with_capacity(triangle_count);

    for _ in 0..triangle_count {
        triangles.push([
            reader.read_u16() as u32,
            reader.read_u16() as u32,
            reader.read_u16() as u32,
        ]);
    }

    Ok((vertices, triangles))
}

/// Serialize raw vertex and triangle data into the micro mesh binary format.
///
/// The format encodes the axis-aligned bounding box of `vertices` as six [`f64`] values, followed
/// by each vertex position discretized to [`u16`] resolution within that box. Triangle indices are
/// stored as [`u16`] values, so all vertex indices must also be within [`u16::MAX`].
///
/// The binary layout is:
/// - 6 × `f64` (little-endian): bounding box min (x, y, z) then max (x, y, z)
/// - 1 × `u16` (little-endian): vertex count
/// - vertex count × 3 × `u16` (little-endian): discretized x, y, z per vertex
/// - 1 × `u32` (little-endian): triangle count
/// - triangle count × 3 × `u16` (little-endian): vertex indices per triangle
///
/// # Errors
///
/// Returns an error if `vertices` contains more than [`u16::MAX`] entries.
pub fn u_mesh_data_to_bytes(vertices: &[Point3], triangles: &[[u32; 3]]) -> Result<Vec<u8>> {
    // Check if the number of vertices is less than u16::MAX
    if vertices.len() > u16::MAX as usize {
        return Err("Mesh has too many vertices for the small format".into());
    }

    let mut output = Vec::new();

    // Get the AABB of the mesh
    let bounds = Aabb3::from_points_ref(vertices);
    output.extend_from_slice(&bounds.mins.x.to_le_bytes());
    output.extend_from_slice(&bounds.mins.y.to_le_bytes());
    output.extend_from_slice(&bounds.mins.z.to_le_bytes());
    output.extend_from_slice(&bounds.maxs.x.to_le_bytes());
    output.extend_from_slice(&bounds.maxs.y.to_le_bytes());
    output.extend_from_slice(&bounds.maxs.z.to_le_bytes());

    // Write the number of vertices
    output.extend_from_slice(&(vertices.len() as u16).to_le_bytes());

    // Write the vertices
    for p in vertices {
        output.extend_from_slice(&to_u16(p.x, bounds.mins.x, bounds.maxs.x).to_le_bytes());
        output.extend_from_slice(&to_u16(p.y, bounds.mins.y, bounds.maxs.y).to_le_bytes());
        output.extend_from_slice(&to_u16(p.z, bounds.mins.z, bounds.maxs.z).to_le_bytes());
    }

    // Write the number of triangles
    output.extend_from_slice(&(triangles.len() as u32).to_le_bytes());

    // Write the triangles
    for triangle in triangles {
        output.extend_from_slice(&(triangle[0] as u16).to_le_bytes());
        output.extend_from_slice(&(triangle[1] as u16).to_le_bytes());
        output.extend_from_slice(&(triangle[2] as u16).to_le_bytes());
    }

    Ok(output)
}

fn to_u16(value: f64, min: f64, max: f64) -> u16 {
    let range = max - min;
    let scale = u16::MAX as f64 / range;
    ((value - min) * scale).round() as u16
}

struct ByteRead<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> ByteRead<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        ByteRead { bytes, offset: 0 }
    }

    fn read_f64(&mut self) -> f64 {
        let value =
            f64::from_le_bytes(self.bytes[self.offset..self.offset + 8].try_into().unwrap());
        self.offset += 8;
        value
    }

    fn read_u16(&mut self) -> u16 {
        let value =
            u16::from_le_bytes(self.bytes[self.offset..self.offset + 2].try_into().unwrap());
        self.offset += 2;
        value
    }

    fn read_u32(&mut self) -> u32 {
        let value =
            u32::from_le_bytes(self.bytes[self.offset..self.offset + 4].try_into().unwrap());
        self.offset += 4;
        value
    }
}
