//! This module has tools for computing a planar section with a mesh. I don't use the built-in
//! parry3d implementation because I've found it to get stuck in infinite loops under certain
//! conditions.

use super::Mesh;
use crate::geom3::Aabb3;
use crate::{Curve3, Line3, Plane3, Point3, Result};
use parry3d_f64::partitioning::TraversalAction;
use parry3d_f64::shape::TriMesh;

impl Mesh {
    pub fn section_with_plane(
        &self,
        plane: &Plane3,
        curve_tol: Option<f64>,
    ) -> Result<Vec<Curve3>> {
        let mut products = Vec::new();

        for face_i in candidate_faces(&self.shape, plane) {
            let ai = self.faces()[face_i as usize][0];
            let bi = self.faces()[face_i as usize][1];
            let ci = self.faces()[face_i as usize][2];
            let a = self.vertices()[ai as usize];
            let b = self.vertices()[bi as usize];
            let c = self.vertices()[ci as usize];

            let ab = edge_intersection(&a, &b, plane);
            let bc = edge_intersection(&b, &c, plane);
            let ca = edge_intersection(&c, &a, plane);

            let data = match (ab, bc, ca) {
                (Some(ab), Some(bc), None) => {},
                (None, Some(bc), Some(ca)) => {},
                (Some(ab), None, Some(ca)) => {},
                _ => Err("Something went wrong with the intersection calculation")?
            }

        }

        todo!()
    }
}

fn edge_intersection(a: &Point3, b: &Point3, plane: &Plane3) -> Option<Point3> {
    if !intersects_edge(a, b, plane) {
        return None;
    }
    // The points can't(?) be equal if they made it through the check above.
    let line = Line3::new(*a, *b - *a);
    let t = line.intersect_plane(plane)?;
    Some(line.at(t))
}

fn candidate_faces(shape: &TriMesh, plane: &Plane3) -> Vec<u32> {
    let mut candidates = Vec::new();
    shape.bvh().traverse(|node| {
        if !aabb_plane(&node.aabb(), plane) {
            return TraversalAction::Prune;
        }

        if let Some(index) = node.leaf_data() {
            let t = shape.triangle(index);
            if let Some(n) = t.normal() {
                if n.cross(&plane.normal).norm_squared() > 1e-10
                    && (intersects_edge(&t.a, &t.b, plane)
                        || intersects_edge(&t.b, &t.c, plane)
                        || intersects_edge(&t.c, &t.a, plane))
                {
                    candidates.push(index);
                }
            };
        }

        TraversalAction::Continue
    });

    candidates
}

fn intersects_edge(a: &Point3, b: &Point3, plane: &Plane3) -> bool {
    let ap = plane.signed_distance_to_point(a).is_sign_positive();
    let bp = plane.signed_distance_to_point(b).is_sign_positive();
    ap != bp
}

fn aabb_plane(aabb: &Aabb3, plane: &Plane3) -> bool {
    let mut pos = false;
    let mut neg = false;
    for v in aabb.vertices().iter() {
        let p = plane.signed_distance_to_point(v).is_sign_positive();
        pos = pos || p;
        neg = neg || !p;
    }

    neg && pos
}


#[cfg(test)]
mod tests {
    use crate::Vector3;
    use super::*;

    #[test]
    fn candidates_box_has_eight() {
        let mesh = Mesh::create_box(2.0, 2.0, 2.0, false);
        let plane = Plane3::new(Vector3::z_axis(), 0.0);
        let candidates = candidate_faces(&mesh.shape, &plane);
        assert_eq!(candidates.len(), 8);
    }

    #[test]
    fn candidates_parallel_face_empty() {
        let vertices = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)];
        let faces = vec![[0, 1, 2]];
        let mesh = Mesh::new(vertices, faces, false);
        let plane = Plane3::new(Vector3::z_axis(), 0.0);
        let candidates = candidate_faces(&mesh.shape, &plane);
        assert!(candidates.is_empty());
    }

}