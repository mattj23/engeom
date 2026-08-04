//! This module has a struct that provides quick lookups of associations between faces and
//! edges in a triangular mesh.

use crate::Mesh3;
use crate::Result;
use crate::common::IndexMask;
use crate::geom3::mesh::edges::edge_key;
use crate::geom3::mesh::patches::PatchLabels;
use faer::prelude::default;
use parry3d_f64::utils::hashmap::HashMap;
use parry3d_f64::utils::hashset::HashSet;

pub struct MeshNav<'a> {
    pub mesh: &'a Mesh3,
    pub face_to_edges: Vec<[[u32; 2]; 3]>,
    pub edge_to_faces: HashMap<[u32; 2], Vec<u32>>,
}

impl<'a> MeshNav<'a> {
    pub fn new(mesh: &'a Mesh3) -> Self {
        let mut face_to_edges = Vec::new();
        let mut edge_to_faces: HashMap<[u32; 2], Vec<u32>> = HashMap::with_hasher(default());

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

    pub fn faces_with_vertex(&self, vertex: u32) -> Vec<u32> {
        let mut indices = Vec::new();
        for (i, face) in self.mesh.faces().iter().enumerate() {
            if face.contains(&vertex) {
                indices.push(i as u32);
            }
        }
        indices
    }

    /// Label every face with the connected patch it belongs to.
    ///
    /// Two faces are in the same patch when a path of shared *edges* connects them. Faces which
    /// touch only at a single vertex are therefore separate patches, which is the useful
    /// definition when the point of the exercise is discarding junk: a piece of surface hanging
    /// off the body by one vertex is not attached in any meaningful sense.
    ///
    /// Patches are numbered in order of their lowest-indexed face, so the result is deterministic.
    ///
    /// This is the primitive behind `patches`, and is what to reach for on a mesh which may split
    /// into many pieces, because it costs one `u32` per face rather than one mask per patch.
    ///
    /// # Arguments
    ///
    /// * `mask`: an optional mask restricting which faces take part, as if the mesh had been
    ///   pruned to it. Excluded faces come back labeled `PatchLabels::NO_PATCH`.
    ///
    /// returns: `Result<PatchLabels>`
    pub fn patch_labels(&self, mask: Option<&IndexMask>) -> Result<PatchLabels> {
        let n_faces = self.mesh.faces().len();

        if let Some(m) = mask
            && m.len() != n_faces
        {
            return Err(format!(
                "A face mask of length {} does not match a mesh with {} faces",
                m.len(),
                n_faces
            )
            .into());
        }

        // Face indices are u32 throughout the adjacency maps, and NO_PATCH claims u32::MAX, so a
        // mesh at that scale could not be labeled unambiguously.
        if n_faces >= PatchLabels::NO_PATCH as usize {
            return Err(format!(
                "A mesh with {} faces is too large to label, the limit is {}",
                n_faces,
                PatchLabels::NO_PATCH - 1
            )
            .into());
        }

        let included = |f: usize| mask.is_none_or(|m| m.get(f));

        let mut labels = vec![PatchLabels::NO_PATCH; n_faces];
        let mut count = 0u32;
        let mut stack: Vec<u32> = Vec::new();

        for start in 0..n_faces {
            if labels[start] != PatchLabels::NO_PATCH || !included(start) {
                continue;
            }

            let label = count;
            count += 1;

            labels[start] = label;
            stack.push(start as u32);

            while let Some(face_index) = stack.pop() {
                for edge in &self.face_to_edges[face_index as usize] {
                    let Some(face_list) = self.edge_to_faces.get(edge) else {
                        continue;
                    };

                    for &neighbor in face_list {
                        let n = neighbor as usize;
                        if labels[n] == PatchLabels::NO_PATCH && included(n) {
                            labels[n] = label;
                            stack.push(neighbor);
                        }
                    }
                }
            }
        }

        PatchLabels::new(labels, count as usize)
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

    /// Returns a list of vertices that are part of non-manifold edges, meaning that the vertices
    /// are boundary vertices but are connected to more than two other vertices. These are the
    /// specific vertices which will result in an error when running the `boundary_loops`
    /// function.
    ///
    /// # Arguments
    ///
    /// * `mask`: an optional mask that filters the faces to consider when finding non-manifold
    ///   boundary vertices.
    ///
    /// returns: Vec<u32, Global>
    pub fn nonmanifold_boundary_vertices(&self, mask: Option<&IndexMask>) -> Vec<u32> {
        let mut duplicates = HashSet::with_hasher(default());
        let edges = self.boundary_edges(mask);
        let mut edge_map = HashMap::with_hasher(default());
        for edge in edges {
            if edge_map.insert(edge[0], edge[1]).is_some() {
                duplicates.insert(edge[0]);
            }
        }

        duplicates.into_iter().collect::<Vec<_>>()
    }

    pub fn boundary_vertices(&self, mask: Option<&IndexMask>) -> Vec<u32> {
        let edges = self.boundary_edges(mask);
        let mut vertices = HashSet::with_hasher(default());
        for edge in edges {
            vertices.insert(edge[0]);
            vertices.insert(edge[1]);
        }

        let mut result = vertices.into_iter().collect::<Vec<_>>();
        result.sort();
        result
    }

    /// Returns a list of boundary loops in the mesh. A boundary loop is a closed path of vertices
    /// that form a loop on the boundary of the mesh. If a mask is provided, the function will
    /// only consider faces that are included in the mask, similar to if the mesh had been pruned
    /// with the mask.
    ///
    /// If a non-manifold edge vertex is found, the function will return an error. To see which
    /// specific vertices are non-manifold, use the `nonmanifold_boundary_vertices` method.
    ///
    /// # Arguments
    ///
    /// * `mask`: an optional mask that filters the faces to consider when finding boundary loops.
    ///
    /// returns: Result<Vec<Vec<u32, Global>, Global>, Box<dyn Error, Global>>
    pub fn boundary_loops(&self, mask: Option<&IndexMask>) -> Result<Vec<Vec<u32>>> {
        let edges = self.boundary_edges(mask);
        let mut edge_map = HashMap::with_hasher(default());
        for edge in edges {
            if edge_map.insert(edge[0], edge[1]).is_some() {
                return Err(format!("Non-manifold edge found: {}-{}", edge[0], edge[1]).into());
            }
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

        Ok(all_loops)
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
    pub fn compute_patch_inside_loop(&self, vertex_loop: &[u32]) -> Result<IndexMask> {
        let mut outside = IndexMask::new(self.mesh.faces().len(), false);
        let mut inside = IndexMask::new(self.mesh.faces().len(), false);
        let mut working = HashSet::with_hasher(default());

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

    fn make_fixture() -> (RasterMask, RasterMapping, Mesh3) {
        let mapping = RasterMapping::new(Point2::new(0.0, 0.0), (100, 100), 1.0, None);
        let mask = mapping.make_mask();
        let filled = mask.not();
        let (indices, faces) = filled.triangle_structure();
        let vertices = indices
            .iter()
            .map(|pi| mapping.point_of_image_point_i(*pi).to_3d())
            .collect::<Vec<_>>();

        let mesh = Mesh3::new(vertices, faces, false);
        (mask, mapping, mesh)
    }

    fn faces_from_mask(mask: &RasterMask, mapping: &RasterMapping, mesh: &Mesh3) -> IndexMask {
        let mut indices = IndexMask::new(mesh.faces().len(), false);
        for (i, face) in mesh.faces().iter().enumerate() {
            let v0 = mapping.image_index_of(&mesh.points()[face[0] as usize].to_2d());
            let v1 = mapping.image_index_of(&mesh.points()[face[1] as usize].to_2d());
            let v2 = mapping.image_index_of(&mesh.points()[face[2] as usize].to_2d());
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
        let labels = nav.patch_labels(Some(&indices)).unwrap();
        assert_eq!(labels.count(), 2, "Expected two patches from the mask");

        let mut counts = labels.face_counts();
        counts.sort();

        assert_eq!(counts[0], 1920, "Smaller patch should have 1920 faces");
        assert_eq!(counts[1], 5020, "Larger patch should have 5020 faces");
    }

    /// The expanded mask form is a view onto the label buffer, so the two must agree face for
    /// face, including on which faces the mask excluded entirely.
    #[test]
    fn patch_labels_agree_with_expanded_masks() {
        let (mut mask, mapping, mesh) = make_fixture();
        mask.draw_rect_mut(Point2I::new(10, 10), Point2I::new(40, 40), true, true);
        mask.draw_rect_mut(Point2I::new(10, 60), Point2I::new(90, 90), true, true);
        let indices = faces_from_mask(&mask, &mapping, &mesh);

        let nav = MeshNav::new(&mesh);
        let labels = nav.patch_labels(Some(&indices)).unwrap();
        let patches = labels.to_masks();

        assert_eq!(labels.count(), patches.len());
        assert_eq!(labels.len(), mesh.faces().len());

        for face in 0..mesh.faces().len() {
            match labels.label_of(face) {
                Some(patch) => {
                    assert!(
                        indices.get(face),
                        "Face {} was labeled but the mask excluded it",
                        face
                    );
                    for (other, m) in patches.iter().enumerate() {
                        assert_eq!(
                            m.get(face),
                            other == patch,
                            "Face {} labeled {} disagrees with mask {}",
                            face,
                            patch,
                            other
                        );
                    }
                }
                None => {
                    assert!(
                        !indices.get(face),
                        "Face {} was left unlabeled but the mask included it",
                        face
                    );
                    for m in patches.iter() {
                        assert!(!m.get(face), "Unlabeled face {} appears in a patch", face);
                    }
                }
            }
        }
    }

    /// Without a mask every face belongs to some patch, and the fixture is a single sheet.
    #[test]
    fn patch_labels_without_a_mask_cover_every_face() {
        let (_, _, mesh) = make_fixture();
        let nav = MeshNav::new(&mesh);
        let labels = nav.patch_labels(None).unwrap();

        assert_eq!(labels.count(), 1, "The fixture is one connected sheet");
        assert!(labels.labels().iter().all(|&l| l == 0));
        assert_eq!(labels.face_counts(), vec![mesh.faces().len()]);
    }

    /// Patches are numbered by their lowest-indexed face, which the old hash-set traversal did not
    /// guarantee. Callers ranking patches depend on the numbering not moving between runs.
    #[test]
    fn patch_labels_are_deterministic() {
        let (mut mask, mapping, mesh) = make_fixture();
        mask.draw_rect_mut(Point2I::new(10, 10), Point2I::new(40, 40), true, true);
        mask.draw_rect_mut(Point2I::new(10, 60), Point2I::new(90, 90), true, true);
        let indices = faces_from_mask(&mask, &mapping, &mesh);

        let nav = MeshNav::new(&mesh);
        let first = nav.patch_labels(Some(&indices)).unwrap();

        for _ in 0..8 {
            assert_eq!(first, nav.patch_labels(Some(&indices)).unwrap());
        }

        // The first labeled face must belong to patch 0.
        let first_labeled = (0..mesh.faces().len())
            .find(|&f| first.label_of(f).is_some())
            .unwrap();
        assert_eq!(first.label_of(first_labeled), Some(0));
    }

    #[test]
    fn patch_labels_reject_a_mismatched_mask() {
        let (_, _, mesh) = make_fixture();
        let nav = MeshNav::new(&mesh);
        let wrong = IndexMask::new(mesh.faces().len() + 1, true);

        assert!(nav.patch_labels(Some(&wrong)).is_err());
    }

    #[test]
    fn boundary_loops() {
        let (mut mask, mapping, mesh) = make_fixture();
        mask.draw_circle_mut(Point2I::new(37, 50), 20, true, false);
        mask.draw_circle_mut(Point2I::new(63, 50), 20, true, false);
        mask.dilate_alternating_norms_mut(1);
        let indices = faces_from_mask(&mask, &mapping, &mesh);

        let nav = MeshNav::new(&mesh);

        let loops = nav.boundary_loops(Some(&indices)).unwrap();
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

        let mut boundary_loops = nav.boundary_loops(Some(&indices)).unwrap();
        boundary_loops.sort_by_key(|a| a.len());
        let mut working = boundary_loops[0].clone();
        working.reverse();

        let inner = nav.compute_patch_inside_loop(&working).unwrap();
        assert_eq!(
            inner.count_true(),
            4960,
            "Inner loop filled should have 4960 faces"
        );
    }
}
