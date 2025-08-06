//! This module has a struct that provides quick lookups of associations between faces and
//! edges in a triangular mesh.

use crate::common::IndexMask;
use crate::geom3::mesh::edges::edge_key;
use crate::Mesh;
use crate::Result;
use parry3d_f64::utils::hashmap::HashMap;
use parry3d_f64::utils::hashset::HashSet;
use std::io::Write;

pub struct MeshNav<'a> {
    pub mesh: &'a Mesh,
    pub face_to_edges: Vec<[[u32; 2]; 3]>,
    pub edge_to_faces: HashMap<[u32; 2], Vec<u32>>,
}

impl<'a> MeshNav<'a> {
    pub fn new(mesh: &'a Mesh) -> Self {
        let mut face_to_edges = Vec::new();
        let mut edge_to_faces: HashMap<[u32; 2], Vec<u32>> = HashMap::new();

        for (i, face) in mesh.faces().iter().enumerate() {
            let e0 = edge_key(&[face[0], face[1]]);
            let e1 = edge_key(&[face[1], face[2]]);
            let e2 = edge_key(&[face[2], face[0]]);
            face_to_edges.push([e0, e1, e2]);

            edge_to_faces.entry(e0).or_default().push(i as u32);
            edge_to_faces.entry(e1).or_default().push(i as u32);
            edge_to_faces.entry(e2).or_default().push(i as u32);
        }

        Self {
            mesh,
            face_to_edges,
            edge_to_faces,
        }
    }

    pub fn patches(&self, mask: Option<&IndexMask>) -> Result<Vec<IndexMask>> {
        let mut results = Vec::new();
        let mut remaining = HashSet::new();
        if let Some(m) = mask {
            for i in m.iter_true() {
                remaining.insert(i as u32);
            }
        } else {
            for i in 0..self.mesh.faces().len() {
                remaining.insert(i as u32);
            }
        }

        while !remaining.is_empty() {
            let mut patch = IndexMask::new(self.mesh.faces().len(), false);
            let start_face = *remaining.iter().next().unwrap();
            remaining.remove(&start_face);

            let mut working_queue = vec![start_face];
            patch.set(start_face as usize, true);

            while let Some(face_index) = working_queue.pop() {
                for &edge in &self.face_to_edges[face_index as usize] {
                    if let Some(face_list) = self.edge_to_faces.get(&edge) {
                        for &neighbor_face in face_list {
                            if remaining.contains(&neighbor_face) {
                                remaining.remove(&neighbor_face);
                                patch.set(neighbor_face as usize, true);
                                working_queue.push(neighbor_face);
                            }
                        }
                    }
                }
            }

            results.push(patch);
        }

        Ok(results)
    }

    /// Gets a list of boundary edges of the mesh. If a mask is provided, it will only consider
    /// edges from faces that are included in the mask, similar to if the mesh had been pruned
    /// with the mask.
    ///
    /// Boundary edges are edges which are only associated with a single face.
    ///
    /// # Arguments
    ///
    /// * `mask`: an optional mask that filters the faces to consider when finding boundary edges.
    ///
    /// returns: Vec<[u32; 2], Global>
    pub fn boundary_edges(&self, mask: Option<&IndexMask>) -> Vec<[u32; 2]> {
        let mut edges = Vec::new();

        for (key, faces) in self.edge_to_faces.iter() {
            let face_indices = if let Some(m) = mask {
                faces
                    .iter()
                    .filter(|&&f| m.get(f as usize))
                    .collect::<Vec<_>>()
            } else {
                faces.iter().collect::<Vec<_>>()
            };

            if face_indices.len() != 1 {
                continue; // Only consider edges with exactly one face
            }

            let face = &self.mesh.faces()[*face_indices[0] as usize];
            let e0 = edge_key(&[face[0], face[1]]);
            let e1 = edge_key(&[face[1], face[2]]);
            let e2 = edge_key(&[face[2], face[0]]);

            if e0 == *key {
                edges.push([face[0], face[1]]);
            } else if e1 == *key {
                edges.push([face[1], face[2]]);
            } else if e2 == *key {
                edges.push([face[2], face[0]]);
            } else {
                panic!(
                    "Edge key {}-{} not found in face {}",
                    key[0], key[1], face_indices[0]
                );
            }
        }

        edges
    }

    /// Returns a list of boundary loops in the mesh. A boundary loop is a closed path of vertices
    /// that form a loop on the boundary of the mesh. If a mask is provided, the function will
    /// only consider faces that are included in the mask, similar to if the mesh had been pruned
    /// with the mask.
    ///
    /// # Arguments
    ///
    /// * `mask`: an optional mask that filters the faces to consider when finding boundary loops.
    ///
    /// returns: Vec<Vec<u32, Global>, Global>
    pub fn boundary_loops(&self, mask: Option<&IndexMask>) -> Vec<Vec<u32>> {
        let edges = self.boundary_edges(mask);
        let mut edge_map = HashMap::new();
        for edge in edges {
            edge_map.insert(edge[0], edge[1]);
        }

        let mut all_loops = Vec::new();
        let mut working = Vec::new();
        let mut queue: HashSet<u32> = edge_map.keys().copied().collect();

        while !queue.is_empty() {
            if let Some(last_id) = working.last() {
                let next_id = edge_map[last_id];
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
                queue.remove(&start_id);
                working.push(start_id);
            }
        }

        // If there's any remaining working loop, add it to the list
        if !working.is_empty() {
            working.reverse();
            all_loops.push(working);
        }

        all_loops
    }


    /// This is a morphological operation that flood-selects a patch of faces that are contained
    /// within a loop of vertices. It can be used to fill holes in a selection.
    ///
    /// The `vertex_loop` is a list of vertex indices that define the loop, and every index to the
    /// next (including the last to the first) must be a valid edge in the mesh. Faces which
    /// contain the edge going in the same order as the loop are considered "outside" the loop,
    /// and will form a boundary for the flood fill operation.  Faces which contain the edge
    /// going in the opposite order are considered "inside" the loop, and will be the seed faces
    /// for the flood select.
    ///
    /// # Arguments
    ///
    /// * `vertex_loop`: a slice of vertex indices that define the loop. The loop must be closed,
    ///   meaning the last vertex connects back to the first vertex, and each consecutive pair of
    ///   indices must be a valid edge that exists in one of the mesh's faces.
    ///
    /// returns: Result<IndexMask, Box<dyn Error, Global>>
    pub fn get_patch_inside_loop(&self, vertex_loop: &[u32]) -> Result<IndexMask> {
        let mut outside = IndexMask::new(self.mesh.faces().len(), false);
        let mut inside = IndexMask::new(self.mesh.faces().len(), false);
        let mut working = HashSet::new();

        // Prepare the inside and outside masks
        for i in 0..vertex_loop.len() {
            let i0 = vertex_loop[i];
            let i1 = vertex_loop[(i + 1) % vertex_loop.len()];
            let key = edge_key(&[i0, i1]);

            let face_indices = self
                .edge_to_faces
                .get(&key)
                .ok_or(format!("The edge {}-{} does not exist in the mesh", i0, i1))?;

            for fi in face_indices {
                let face = &self.mesh.faces()[*fi as usize];
                if i0 == face[0] && i1 == face[1]
                    || i0 == face[1] && i1 == face[2]
                    || i0 == face[2] && i1 == face[0]
                {
                    outside.set(*fi as usize, true);
                } else {
                    working.insert(*fi as usize);
                }
            }
        }

        // Now we'll expand the inside mask
        while !working.is_empty() {
            let fi = *working.iter().next().unwrap();
            working.remove(&fi);
            inside.set(fi, true);

            for edge in self.face_to_edges[fi].iter() {
                for f_op in self.edge_to_faces[edge].iter() {
                    if !outside.get(*f_op as usize) && !inside.get(*f_op as usize) {
                        // If the face is not already marked as outside or inside, add it to working
                        working.insert(*f_op as usize);
                    }
                }
            }
        }

        Ok(inside)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raster2::{Point2I, RasterMapping, RasterMask};
    use crate::{Point2, To2D, To3D};

    fn make_fixture() -> (RasterMask, RasterMapping, Mesh) {
        let mapping = RasterMapping::new(Point2::new(0.0, 0.0), (100, 100), 1.0, None);
        let mask = mapping.make_mask();
        let filled = mask.not();
        let (indices, faces) = filled.triangle_structure();
        let vertices = indices
            .iter()
            .map(|pi| mapping.point_of_image_point_i(*pi).to_3d())
            .collect::<Vec<_>>();

        let mesh = Mesh::new(vertices, faces, false);
        (mask, mapping, mesh)
    }

    fn faces_from_mask(mask: &RasterMask, mapping: &RasterMapping, mesh: &Mesh) -> IndexMask {
        let mut indices = IndexMask::new(mesh.faces().len(), false);
        for (i, face) in mesh.faces().iter().enumerate() {
            let v0 = mapping.image_index_of(&mesh.vertices()[face[0] as usize].to_2d());
            let v1 = mapping.image_index_of(&mesh.vertices()[face[1] as usize].to_2d());
            let v2 = mapping.image_index_of(&mesh.vertices()[face[2] as usize].to_2d());
            if mask.get_point(v0) || mask.get_point(v1) || mask.get_point(v2) {
                indices.set(i, true);
            }
        }
        indices
    }

    #[test]
    fn patch_split() {
        let (mut mask, mapping, mesh) = make_fixture();
        mask.draw_rect_mut(Point2I::new(10, 10), Point2I::new(40, 40), true, true);
        mask.draw_rect_mut(Point2I::new(10, 60), Point2I::new(90, 90), true, true);
        let indices = faces_from_mask(&mask, &mapping, &mesh);

        let nav = MeshNav::new(&mesh);
        let patches = nav.patches(Some(&indices)).unwrap();
        assert_eq!(patches.len(), 2, "Expected two patches from the mask");

        let f0 = patches[0].count_true();
        let f1 = patches[1].count_true();

        assert_eq!(f0.min(f1), 1920, "Smaller patch should have 1920 faces");
        assert_eq!(f0.max(f1), 5020, "Smaller patch should have 5020 faces");
    }

    #[test]
    fn boundary_loops() {
        let (mut mask, mapping, mesh) = make_fixture();
        mask.draw_circle_mut(Point2I::new(37, 50), 20, true, false);
        mask.draw_circle_mut(Point2I::new(63, 50), 20, true, false);
        mask.dilate_alternating_norms_mut(1);
        let indices = faces_from_mask(&mask, &mapping, &mesh);

        let nav = MeshNav::new(&mesh);

        let loops = nav.boundary_loops(Some(&indices));
        assert_eq!(loops.len(), 4, "Expected four boundary loops from the mask");

        let mut lengths = loops.iter().map(|loop_| loop_.len()).collect::<Vec<_>>();
        lengths.sort();

        assert_eq!(lengths[0], 62, "First loop should have 62 points");
        assert_eq!(lengths[1], 122, "Second loop should have 122 points");
        assert_eq!(lengths[2], 122, "Third loop should have 122 points");
        assert_eq!(lengths[3], 210, "Last loop should have 210 points");
    }

    #[test]
    fn extract_inner() {
        let (mut mask, mapping, mesh) = make_fixture();
        mask.draw_circle_mut(Point2I::new(50, 50), 30, true, false);
        mask.dilate_alternating_norms_mut(1);
        let indices = faces_from_mask(&mask, &mapping, &mesh);
        let nav = MeshNav::new(&mesh);

        let mut boundary_loops = nav.boundary_loops(Some(&indices));
        boundary_loops.sort_by(|a, b| a.len().cmp(&b.len()));
        let mut working = boundary_loops[0].clone();
        working.reverse();

        let inner = nav.get_patch_inside_loop(&working).unwrap();
        assert_eq!(inner.count_true(), 4960, "Inner loop filled should have 4960 faces");
    }
}
