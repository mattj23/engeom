//! This module has implementations of different ways of filtering/reducing a mesh

use crate::common::indices::index_vec;
use crate::common::{SelectOp, Selection};
use crate::{Mesh, Point3, SurfacePoint3, UnitVec3, Vector3};
use itertools::Itertools;
use std::collections::{HashMap, HashSet};

pub struct TriangleFilter<'a> {
    mesh: &'a Mesh,
    indices: HashSet<usize>,
}

impl TriangleFilter<'_> {
    /// Collect the indices of the triangles that have been filtered
    pub fn collect(self) -> Vec<usize> {
        self.indices.into_iter().collect()
    }

    /// Create a new mesh from the filtered indices
    pub fn create_mesh(self) -> Mesh {
        let i = self.indices.into_iter().collect_vec();
        self.mesh.create_from_indices(&i)
    }

    /// Get the indices of the triangles which would need to be checked for an operation of the
    /// specified type. If the operation is `SelectOp::Add`, then the triangles that are not in the
    /// current selection will be returned. If the operation is `SelectOp::Remove`, or
    /// `SelectOp::Keep` then the triangles that are in the current selection will be returned.
    fn to_check(&self, mode: SelectOp) -> Vec<usize> {
        match mode {
            SelectOp::Add => (0..self.mesh.faces().len())
                .filter(|i| !self.indices.contains(i))
                .collect(),
            SelectOp::Remove | SelectOp::Keep => self.indices.iter().copied().collect(),
        }
    }
    fn mutate_pass_list(mut self, mode: SelectOp, pass_list: Vec<usize>) -> Self {
        match mode {
            SelectOp::Add => {
                for i in pass_list {
                    self.indices.insert(i);
                }
            }
            SelectOp::Remove => {
                for i in pass_list {
                    self.indices.remove(&i);
                }
            }
            SelectOp::Keep => {
                let check_set: HashSet<usize> = pass_list.into_iter().collect();
                self.indices.retain(|i| check_set.contains(i));
            }
        };

        self
    }

    fn mutate(mut self, mode: SelectOp, predicate: &dyn Fn(usize, &Mesh) -> bool) -> Self {
        match mode {
            SelectOp::Add => {
                for i in 0..self.mesh.faces().len() {
                    if !self.indices.contains(&i) && predicate(i, self.mesh) {
                        self.indices.insert(i);
                    }
                }
            }
            SelectOp::Remove => {
                self.indices.retain(|&i| !predicate(i, self.mesh));
            }
            SelectOp::Keep => {
                self.indices.retain(|&i| predicate(i, self.mesh));
            }
        };

        self
    }

    pub fn facing(self, normal: &Vector3, angle: f64, mode: SelectOp) -> Self {
        self.mutate(mode, &|i, m| {
            let n = m.shape.triangle(i as u32).normal();
            if let Some(nv) = n {
                nv.angle(normal) < angle
            } else {
                false
            }
        })
    }

    /// Modify the list of indices to only include triangles that are within a certain distance of
    /// their closest projection onto another mesh. The distance can require that all points of the
    /// triangle are within the tolerance, or just one.
    ///
    /// There are two additional optional tolerances that can be applied.
    ///
    /// 1. A planar tolerance, which checks the distance of the vertex projected onto the plane of
    ///    the reference mesh triangle and looks at how far it is from the projection point. This
    ///    is useful to filter out triangles that go past the edge of the reference mesh.
    /// 2. An angle tolerance, which checks the angle between the normal of the current triangle
    ///    and the normal of the reference triangle. This is useful to filter out triangles that
    ///    are not facing the same direction as the reference mesh.
    ///
    /// # Arguments
    ///
    /// * `other`:
    /// * `all_points`:
    /// * `distance_tol`:
    /// * `planar_tol`:
    /// * `angle_tol`:
    ///
    /// returns: TriangleFilter
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn near_mesh(
        self,
        other: &Mesh,
        all_points: bool,
        distance_tol: f64,
        planar_tol: Option<f64>,
        angle_tol: Option<f64>,
        mode: SelectOp,
    ) -> Self {
        let mut check = MeshNearCheck::new(self.mesh, other, distance_tol, planar_tol, angle_tol);
        let to_check = self.to_check(mode);
        let passes = to_check
            .into_iter()
            .filter(|&i| {
                let tri = self.mesh.faces()[i];
                let face = self.mesh.shape.triangle(i as u32);

                if all_points {
                    check.near_check(tri[0], face.normal())
                        && check.near_check(tri[1], face.normal())
                        && check.near_check(tri[2], face.normal())
                } else {
                    check.near_check(tri[0], face.normal())
                        || check.near_check(tri[1], face.normal())
                        || check.near_check(tri[2], face.normal())
                }
            })
            .collect::<Vec<_>>();

        self.mutate_pass_list(mode, passes)
    }
}

impl Mesh {
    /// Start an operation to filter the faces of the mesh. This function will return a filter
    /// handle that can be used to add or remove faces from the selection while maintaining
    /// an immutable reference to the mesh.
    ///
    /// The filter can be started with no faces selected (`Selection::None`), all faces selected
    /// (`Selection::All`), or a specific set of faces selected (`Selection::Indices(Vec<usize>)`).
    /// Each successive filter operation will modify the selection the selected indices.
    ///
    /// # Arguments
    ///
    /// * `start`: The initial selection of faces to start with, either `Selection::None`,
    ///   `Selection::All`, or `Selection::Indices(Vec<usize>)`
    ///
    /// returns: TriangleFilter
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn face_select(&self, start: Selection) -> TriangleFilter {
        let indices = match start {
            Selection::None => vec![],
            Selection::All => index_vec(None, self.faces().len()),
            Selection::Indices(i) => i,
        }
        .into_iter()
        .collect();
        TriangleFilter {
            mesh: self,
            indices,
        }
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
    /// use std::f64::consts::PI;
    /// use engeom::{Mesh, Vector3, SelectOp, Selection};
    /// let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
    /// let indices = mesh.face_select(Selection::None)
    ///     .facing(&Vector3::z(), PI / 2.0, SelectOp::Add)
    ///     .collect();
    /// let new_mesh = mesh.create_from_indices(&indices);
    ///
    /// assert_eq!(new_mesh.faces().len(), 2);
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

        let vertices: Vec<Point3> = to_keep
            .iter()
            .map(|i| self.vertices()[*i as usize])
            .collect();

        let triangles = indices
            .iter()
            .map(|i| {
                let t = self.faces()[*i];
                [map_back[&t[0]], map_back[&t[1]], map_back[&t[2]]]
            })
            .collect_vec();

        Self::new(vertices, triangles, false)
    }

    fn unique_vertices(&self, triangle_indices: &[usize]) -> Vec<u32> {
        let mut to_save = HashSet::new();
        for i in triangle_indices {
            let t = self.faces()[*i];
            to_save.insert(t[0]);
            to_save.insert(t[1]);
            to_save.insert(t[2]);
        }

        // Now we can sort them in order
        let mut keep_order = to_save.iter().copied().collect_vec();
        keep_order.sort_unstable();

        keep_order
    }
}

struct MeshNearCheck<'a> {
    this_mesh: &'a Mesh,
    ref_mesh: &'a Mesh,
    checked: HashMap<u32, bool>,
    distance_tol: f64,
    planar_tol: Option<f64>,
    angle_tol: Option<f64>,
}

impl<'a> MeshNearCheck<'a> {
    fn new(
        this_mesh: &'a Mesh,
        ref_mesh: &'a Mesh,
        distance_tol: f64,
        planar_tol: Option<f64>,
        angle_tol: Option<f64>,
    ) -> Self {
        Self {
            this_mesh,
            ref_mesh,
            checked: HashMap::new(),
            distance_tol,
            planar_tol,
            angle_tol,
        }
    }

    fn store_and_return(&mut self, vertex_index: u32, result: bool) -> bool {
        self.checked.insert(vertex_index, result);
        result
    }

    fn near_check(&mut self, vertex_index: u32, face_normal: Option<UnitVec3>) -> bool {
        if let Some(&checked) = self.checked.get(&vertex_index) {
            checked
        } else {
            let p = self.this_mesh.vertices()[vertex_index as usize];

            let is_ok = if let Some((prj, ri, _loc)) =
                self.ref_mesh.project_with_max_dist(&p, self.distance_tol)
            {
                if self.planar_tol.is_none() && self.angle_tol.is_none() {
                    true
                } else if let Some(rn) = self.ref_mesh.shape.triangle(ri).normal() {
                    // We need to get the normal of the reference triangle
                    let rsp = SurfacePoint3::new(prj.point, rn);

                    let check_planar = if let Some(planar_tol) = self.planar_tol {
                        rsp.planar_distance(&p) <= planar_tol
                    } else {
                        true
                    };

                    let check_angle = if let Some(angle_tol) = self.angle_tol {
                        if let Some(face_normal) = face_normal {
                            face_normal.angle(&rn) <= angle_tol
                        } else {
                            // No face normal, so we can't check the angle, assume it's bad?
                            false
                        }
                    } else {
                        true
                    };

                    check_planar && check_angle
                } else {
                    false
                }
            } else {
                false
            };

            self.store_and_return(vertex_index, is_ok)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SelectOp::Add;
    use std::f64::consts::PI;

    #[test]
    fn test_triangles_facing() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let indices = mesh
            .face_select(Selection::None)
            .facing(&Vector3::z(), PI / 2.0, Add)
            .collect();

        assert_eq!(indices.len(), 2);

        let new_mesh = mesh.create_from_indices(&indices);
        assert_eq!(new_mesh.faces().len(), 2);

        for t in new_mesh.tri_mesh().triangles() {
            let n = t.normal().unwrap();
            assert!(n.dot(&Vector3::z()) > 0.0);
        }
    }
}
