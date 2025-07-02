//! Mesh edge structure

use super::Mesh;
use crate::{Point3, Result};
use std::collections::{HashMap, HashSet};

pub struct MeshEdges<'a> {
    /// The original mesh associated with the edge structure
    mesh: &'a Mesh,

    /// A list of edges in the mesh. Each edge consists of two indices into the `vertices` list,
    /// which are the two vertices of the edge. The order of the vertices in the edge is not
    /// the same as the order of the same two vertices in the `vertices` list.
    pub edges: Vec<[u32; 2]>,

    /// A list of the lengths of each edge in the mesh. The order of the lengths is the same as the
    /// order of the edges in the `edges` list.
    pub edge_lengths: Vec<f64>,

    /// A list of edges associated with each face. A face at index `i` in `face_edges` corresponds
    /// with the face at index `i` in `faces`. The `faces` list references the three vertices,
    /// while the `face_edges` list references the three edges of the face.
    pub face_edges: Vec<[u32; 3]>,

    /// A list of the different boundaries in the mesh. Each boundary is a list of indices into the
    /// vertices list that form a loop.
    pub boundary_loops: Vec<Vec<u32>>,
}

// impl MeshEdges {
//     pub fn mesh(&self) -> &Mesh {
//         &self.mesh
//     }
//
//     /// Get a reference to the AABB of the underlying mesh in the local coordinate system.
//     pub fn aabb(&self) -> &Aabb3 {
//         self.mesh.shape.local_aabb()
//     }
//
//     /// Gets a reference to the underlying `TriMesh` object to provide direct access to
//     /// the `parry3d` API.
//     pub fn tri_mesh(&self) -> &TriMesh {
//         &self.mesh.shape
//     }
//
//     /// Return a flag indicating whether the mesh is considered "solid" or not for the purposes of
//     /// distance queries. If a mesh is "solid", then distance queries for points on the inside of
//     /// the mesh will return a zero distance.
//     pub fn is_solid(&self) -> bool {
//         self.mesh.is_solid
//     }
//
//     /// Get a reference to the vertices of the mesh.
//     pub fn vertices(&self) -> &[Point3] {
//         self.mesh.shape.vertices()
//     }
//
//     /// Get a reference to the face indices of the mesh.
//     pub fn faces(&self) -> &[[u32; 3]] {
//         self.mesh.shape.indices()
//     }
// }

impl<'a> MeshEdges<'a> {
    pub fn mesh(&self) -> &Mesh {
        self.mesh
    }

    /// Get a reference to the vertices of the mesh.
    pub fn vertices(&self) -> &[Point3] {
        self.mesh.shape.vertices()
    }

    /// Get a reference to the face indices of the mesh.
    pub fn faces(&self) -> &[[u32; 3]] {
        self.mesh.shape.indices()
    }

    pub fn new(mesh: &'a Mesh) -> Result<Self> {
        let (edges, face_edges, boundary_loops) = identify_edges(mesh.faces())?;

        let edge_lengths = edges
            .iter()
            .map(|edge| {
                let v0 = mesh.vertices()[edge[0] as usize];
                let v1 = mesh.vertices()[edge[1] as usize];
                (v1 - v0).norm()
            })
            .collect();

        Ok(Self {
            mesh,
            edges,
            edge_lengths,
            face_edges,
            boundary_loops,
        })
    }
}

/// Given an edge, return a key that can be used to identify the vertices that are connected by the
/// edge without regard to the order of the vertices.
pub fn edge_key(edge: &[u32; 2]) -> [u32; 2] {
    let x = edge[0].min(edge[1]);
    let y = edge[0].max(edge[1]);
    [x, y]
}

/// Generates a complete list of edges from the faces of a mesh, naively, and in the order that
/// the faces are presented. The edges are not deduplicated nor their direction normalized. The
/// resulting list will be 3x the length of the original faces. Elements 0-2 will be from face 0,
/// elements 3-5 will be from face 1, and so on.
pub fn naive_edges(faces: &[[u32; 3]]) -> Vec<[u32; 2]> {
    let mut edges = Vec::new();
    for face in faces {
        edges.push([face[1], face[2]]);
        edges.push([face[2], face[0]]);
        edges.push([face[0], face[1]]);
    }
    edges
}

/// Given a reference list of all edges, return a sorted list of unique edges and the number of
/// times each edge appeared in the original list.
pub fn unique_edges(all_edges: &[[u32; 2]]) -> Vec<([u32; 2], usize)> {
    let mut unique = HashMap::new();
    for edge in all_edges {
        let key = edge_key(edge);
        let count = unique.entry(key).or_insert(0);
        *count += 1;
    }

    let mut unique_count: Vec<_> = unique.into_iter().collect();
    unique_count.sort();
    unique_count
}

fn boundary_loops(boundary_map: HashMap<u32, u32>) -> Vec<Vec<u32>> {
    let mut all_loops = Vec::new();
    let mut working = Vec::new();
    let mut queue: HashSet<u32> = boundary_map.keys().copied().collect();

    while !queue.is_empty() {
        if let Some(last_id) = working.last() {
            let next_id = boundary_map[last_id];
            queue.remove(&next_id);

            if *working.first().unwrap() == next_id {
                working.reverse();
                all_loops.push(working);
                working = Vec::new();
            } else {
                working.push(next_id);
            }
        } else {
            let start_id = *queue.iter().next().unwrap();
            working.push(start_id);
        }
    }

    all_loops
}

fn identify_edges(faces: &[[u32; 3]]) -> Result<(Vec<[u32; 2]>, Vec<[u32; 3]>, Vec<Vec<u32>>)> {
    // The direct edges are the edges that are directly defined by the faces, kept in the same
    // order as they are defined in the faces.
    let direct_edges = naive_edges(faces);

    // We need to identify the unique edges, and put them into a sorted order where we can still
    // map from an original face's edge to the unique edge
    let unique_edge_count = unique_edges(&direct_edges);

    // Boundary edges are edges that only appear once in the mesh. Non-manifold edges are ones
    // that appear more than twice. Any non-manifold edges will cause this function to return an
    // error
    if unique_edge_count.iter().any(|(_, count)| *count > 2) {
        return Err("Non-manifold edges detected".into());
    }

    // Now we can create a mapping from the original edge to the corresponding unique edge
    let to_unique_index: HashMap<[u32; 2], usize> = unique_edge_count
        .iter()
        .enumerate()
        .map(|(i, (edge, _))| (*edge, i))
        .collect();

    // Let's remap the face edges to the unique edges and build the boundary map at the same time
    let mut boundary_map = HashMap::new();
    let mut face_edges = Vec::new();
    for face_chunk in direct_edges.chunks(3) {
        let i0 = to_unique_index[&edge_key(&face_chunk[0])];
        let i1 = to_unique_index[&edge_key(&face_chunk[1])];
        let i2 = to_unique_index[&edge_key(&face_chunk[2])];
        face_edges.push([i0 as u32, i1 as u32, i2 as u32]);

        if unique_edge_count[i0].1 == 1 {
            boundary_map.insert(face_chunk[0][0], face_chunk[0][1]);
        }
        if unique_edge_count[i1].1 == 1 {
            boundary_map.insert(face_chunk[1][0], face_chunk[1][1]);
        }
        if unique_edge_count[i2].1 == 1 {
            boundary_map.insert(face_chunk[2][0], face_chunk[2][1]);
        }
    }

    let loops = boundary_loops(boundary_map);
    let edges = unique_edge_count.iter().map(|(edge, _)| *edge).collect();

    Ok((edges, face_edges, loops))
}
