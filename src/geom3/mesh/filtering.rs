//! This module has implementations of different ways of filtering/reducing a mesh

use crate::{Mesh, Point3, Result, Vector3};
use itertools::Itertools;
use parry2d_f64::utils::hashmap::HashMap;
use std::collections::HashSet;

pub struct MeshFilter<'a> {
    mesh: &'a Mesh,
}


impl<'a> MeshFilter<'a> {

    /// Filter a list of indices to only include triangles that have a positive dot product with
    /// the specified direction vector.
    ///
    /// # Arguments
    ///
    /// * `normal`: a test normal direction to check against the triangle normals
    ///
    /// returns: Vec<usize, Global>
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Mesh, Vector3};
    /// let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
    /// let indices = mesh.filters().triangles_facing(&Vector3::z());
    /// assert_eq!(indices.len(), 2);
    /// ```
    pub fn triangles_facing(&self, normal: &Vector3) -> Vec<usize> {
        self.mesh.shape.triangles().enumerate()
            .filter_map(|(i, face)| {
                if let Some(n) = face.normal() {
                    if n.dot(normal) > 0.0 {
                        Some(i)
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect()
    }
}


impl Mesh {
    /// Get a filter handle for the mesh
    pub fn filters(&self) -> MeshFilter {
        MeshFilter { mesh: self }
    }

    /// Create a new mesh from a list of triangle indices. The indices correspond with elements in
    /// the `triangles()` slice. This function will iterate through the triangle indices,
    /// taking the three vertices associated with each index and marking them for inclusion in the
    /// new mesh. Then it will recreate the triangles, remapping them to the new vertex indices.
    ///
    /// # Arguments
    ///
    /// * `indices`: A slice of usize values that correspond to the indices of the triangles in the
    ///   original mesh. There cannot be any duplicate indices, or the function will return a
    ///   non-manifold mesh.
    ///
    /// returns: Mesh
    ///
    /// # Examples
    ///
    /// ```
    /// use engeom::{Mesh, Vector3};
    /// let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
    /// let indices = mesh.filters().triangles_facing(&Vector3::z());
    /// let new_mesh = mesh.create_from_indices(&indices);
    ///
    /// assert_eq!(new_mesh.triangles().len(), 2);
    /// assert_eq!(new_mesh.vertices().len(), 4);
    /// ```
    pub fn create_from_indices(&self, indices: &[usize]) -> Self {
        let to_keep = self.unique_vertices(indices);
        // The map_back array will map the old vertex indices to the new ones
        let map_back: HashMap<u32, u32> = to_keep
            .iter()
            .enumerate()
            .map(|(i, v)| (*v, i as u32))
            .collect();

        let vertices: Vec<Point3> = to_keep.iter().map(|i| self.vertices()[*i as usize]).collect();

        let triangles = indices
            .iter()
            .map(|i| {
                let t = self.triangles()[*i];
                [map_back[&t[0]], map_back[&t[1]], map_back[&t[2]]]
            })
            .collect_vec();

        Self::new(vertices, triangles, false)
    }

    fn unique_vertices(&self, triangle_indices: &[usize]) -> Vec<u32> {
        let mut to_save = HashSet::new();
        for i in triangle_indices {
            let t = self.triangles()[*i];
            to_save.insert(t[0]);
            to_save.insert(t[1]);
            to_save.insert(t[2]);
        }

        // Now we can sort them in order
        let mut keep_order = to_save.iter().map(|i| *i).collect_vec();
        keep_order.sort_unstable();

        keep_order
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_triangles_facing() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let indices = mesh.filters().triangles_facing(&Vector3::new(0.0, 0.0, 1.0));

        assert_eq!(indices.len(), 2);

        let new_mesh = mesh.create_from_indices(&indices);
        assert_eq!(new_mesh.triangles().len(), 2);

        for t in new_mesh.tri_mesh().triangles() {
            let n = t.normal().unwrap();
            assert!(n.dot(&Vector3::new(0.0, 0.0, 1.0)) > 0.0);
        }
    }

}