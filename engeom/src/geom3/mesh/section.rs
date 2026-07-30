//! This module has tools for computing a planar section with a mesh. I don't use the built-in
//! parry3d implementation because I've found it to get stuck in infinite loops under certain
//! conditions.

use super::Mesh3;
use crate::geom3::Aabb3;
use crate::geom3::mesh::edges::edge_key;
use crate::{Curve3, IndexMask, Line3, Plane3, Point3, Result, Vector3};
use parry3d_f64::partitioning::TraversalAction;
use parry3d_f64::shape::TriMesh;
use std::collections::{HashMap, HashSet};

impl Mesh3 {
    pub fn section_with_plane(
        &self,
        plane: &Plane3,
        curve_tol: Option<f64>,
        faces: Option<&IndexMask>,
    ) -> Result<Vec<Curve3>> {
        let curve_tol = curve_tol.unwrap_or(1e-6);

        if let Some(mask) = faces
            && mask.len() != self.face_count()
        {
            return Err(format!(
                "A face mask of length {} does not match a mesh with {} faces",
                mask.len(),
                self.face_count()
            )
            .into());
        }

        // First, we'll find all triangles that intersect the plane and produce segments from them.
        // Each segment will have two endpoints. The order of the points will be such that when
        // transformed into the plane, the 2D segment normal will point in the same direction as
        // the triangle's original normal, preserving the inside/outside relationship from the
        // mesh.
        let mut segments = Vec::new();
        for face_i in candidate_faces(&self.shape, plane, faces) {
            let Some(tri_n) = self.shape.triangle(face_i).normal() else {
                continue;
            };

            let ai = self.faces()[face_i as usize][0];
            let bi = self.faces()[face_i as usize][1];
            let ci = self.faces()[face_i as usize][2];
            let a = self.points()[ai as usize];
            let b = self.points()[bi as usize];
            let c = self.points()[ci as usize];

            let ab = edge_intersection(&a, &b, plane);
            let bc = edge_intersection(&b, &c, plane);
            let ca = edge_intersection(&c, &a, plane);

            let data = match (ab, bc, ca) {
                (Some(ab), Some(bc), None) => {
                    TriIntr::new(ab, bc, edge_key(&[ai, bi]), edge_key(&[bi, ci]))
                }
                (None, Some(bc), Some(ca)) => {
                    TriIntr::new(bc, ca, edge_key(&[bi, ci]), edge_key(&[ci, ai]))
                }
                (Some(ab), None, Some(ca)) => {
                    TriIntr::new(ab, ca, edge_key(&[ai, bi]), edge_key(&[ci, ai]))
                }
                _ => Err("Something went wrong with the intersection calculation")?,
            };

            if tri_n.cross(&data.direction()).dot(&plane.normal) < 0.0 {
                segments.push(data.reversed());
            } else {
                segments.push(data);
            }
        }

        // We're going to find the count of the different edge keys.
        let mut edge_count = HashMap::<[u32; 2], usize>::new();
        for segment in segments.iter() {
            *edge_count.entry(segment.key_a).or_insert(0) += 1;
            *edge_count.entry(segment.key_b).or_insert(0) += 1;
        }

        // Terminations occur at keys that have a count of 1 or >2.
        let terminations = edge_count
            .iter()
            .filter(|(_, count)| **count == 1 || **count > 2)
            .map(|(&key, _)| key)
            .collect::<HashSet<[u32; 2]>>();

        // Now we'll start bunching the segments into groups of connected curves. We'll work until
        // every segment is accounted for, even if it's only a single segment long. Endpoints can
        // be connected if they have the same edge key, but point b can only be connected to point
        // a. Connections continue forward until `key_b` is in the terminations, and then reverse
        // until `key_a` is in the terminations.
        let mut connected = Vec::new();
        let mut work_bag = (0..segments.len()).collect::<Vec<_>>();

        while !work_bag.is_empty() {
            let mut current = vec![work_bag.pop().unwrap()];

            // Backward search
            while !terminations.contains(&segments[current[0]].key_a) {
                // Find the segment which has a key_b equal to the current key_a.
                let key_a = segments[current[0]].key_a;
                let mut did_something = false;
                for i in 0..work_bag.len() {
                    if segments[work_bag[i]].key_b == key_a {
                        current.insert(0, work_bag.remove(i));
                        did_something = true;
                        break;
                    }
                }

                if !did_something {
                    break;
                }
            }

            // Forward search
            while !terminations.contains(&segments[*current.last().unwrap()].key_b)
                && !is_loop(&segments, &current)
                && !work_bag.is_empty()
            {
                // Find the segment which has a key_a equal to the current key_b.
                let key_b = segments[*current.last().unwrap()].key_b;
                let mut did_something = false;
                for i in 0..work_bag.len() {
                    if segments[work_bag[i]].key_a == key_b {
                        current.push(work_bag.remove(i));
                        did_something = true;
                        break;
                    }
                }

                if !did_something {
                    break;
                }
            }

            connected.push(current);
        }

        let mut results = Vec::new();
        for group in connected.iter() {
            let mut curve_points = group.iter().map(|&i| segments[i].a).collect::<Vec<_>>();
            curve_points.push(segments[group[group.len() - 1]].b);
            results.push(Curve3::from_points(&curve_points, curve_tol)?);
        }

        Ok(results)
    }
}

fn is_loop(segments: &[TriIntr], working: &[usize]) -> bool {
    let start = segments[working[0]].a;
    let end = segments[working[working.len() - 1]].b;

    start == end
}

struct TriIntr {
    a: Point3,
    b: Point3,
    key_a: [u32; 2],
    key_b: [u32; 2],
}

impl TriIntr {
    fn new(a: Point3, b: Point3, key_a: [u32; 2], key_b: [u32; 2]) -> Self {
        Self { a, b, key_a, key_b }
    }

    fn reversed(&self) -> Self {
        Self::new(self.b, self.a, self.key_b, self.key_a)
    }

    fn direction(&self) -> Vector3 {
        (self.b - self.a).normalize()
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

fn candidate_faces(shape: &TriMesh, plane: &Plane3, faces: Option<&IndexMask>) -> Vec<u32> {
    let mut candidates = Vec::new();
    shape.bvh().traverse(|node| {
        if !aabb_plane(&node.aabb(), plane) {
            return TraversalAction::Prune;
        }

        if let Some(index) = node.leaf_data() {
            let t = shape.triangle(index);
            if let Some(n) = t.normal()
                && n.cross(&plane.normal).norm_squared() > 1e-10
                && (intersects_edge(&t.a, &t.b, plane)
                    || intersects_edge(&t.b, &t.c, plane)
                    || intersects_edge(&t.c, &t.a, plane))
            {
                if let Some(mask) = faces {
                    if mask.get(index as usize) {
                        candidates.push(index);
                    }
                } else {
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
    use super::*;
    use crate::common::points::dist;
    use crate::{Iso3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::TAU;

    #[test]
    fn candidates_box_has_eight() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, false);
        let plane = Plane3::new(Vector3::z_axis(), 0.0);
        let candidates = candidate_faces(&mesh.shape, &plane, None);
        assert_eq!(candidates.len(), 8);
    }

    #[test]
    fn candidates_parallel_face_empty() {
        let vertices = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let faces = vec![[0, 1, 2]];
        let mesh = Mesh3::new(vertices, faces, false);
        let plane = Plane3::new(Vector3::z_axis(), 0.0);
        let candidates = candidate_faces(&mesh.shape, &plane, None);
        assert!(candidates.is_empty());
    }

    #[test]
    fn single_loop() -> Result<()> {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, false);

        let curves = mesh.section_with_plane(&Plane3::xy(), None, None)?;

        assert_eq!(curves.len(), 1);

        Ok(())
    }

    #[test]
    fn unit_cylinder_in_xz_plane_creates_one_circle_curve() {
        let mesh = Mesh3::create_cylinder(1.0, 2.0, 256);
        let plane = Plane3::xz();

        let curves = mesh
            .section_with_plane(&plane, Some(1.0e-10), None)
            .unwrap();
        assert_eq!(curves.len(), 1);

        let curve = &curves[0];
        assert!(curve.count() >= 3);

        for vertex in curve.vertices() {
            assert_relative_eq!(vertex.y, 0.0, epsilon = 1.0e-12);

            let radius = (vertex.x * vertex.x + vertex.z * vertex.z).sqrt();
            // The tolerance has to be high enough to account for the fact that the cylinder
            // faces have a diagonal in them and where they pass through y=0 is halfway between
            // the arc endpoints formed by the vertices that were deliberately placed at the radius
            assert_relative_eq!(radius, 1.0, epsilon = 1.0e-4);
        }

        assert_relative_eq!(curve.length(), TAU, epsilon = 1.0e-2);
    }

    #[test]
    fn two_unit_cylinders_in_xz_plane_create_two_circle_curves() {
        let mut mesh = Mesh3::create_cylinder(1.0, 2.0, 256);
        mesh.transform_in_place(&Iso3::translation(0.0, 0.0, -2.0));
        let mut m1 = Mesh3::create_cylinder(1.0, 2.0, 256);
        m1.transform_in_place(&Iso3::translation(0.0, 0.0, 2.0));

        mesh.append_in_place(&m1).unwrap();

        let plane = Plane3::xz();

        let curves = mesh
            .section_with_plane(&plane, Some(1.0e-10), None)
            .unwrap();
        assert_eq!(curves.len(), 2);

        for curve in curves.iter() {
            assert!(curve.count() >= 3);
            assert_relative_eq!(curve.length(), TAU, epsilon = 1.0e-2);

            let expected_center = if curve.vertices()[0].z > 0.0 {
                Point3::new(0.0, 0.0, 2.0)
            } else {
                Point3::new(0.0, 0.0, -2.0)
            };
            for vertex in curve.vertices() {
                assert_relative_eq!(dist(&expected_center, vertex), 1.0, epsilon = 1.0e-4);
            }
        }
    }
}
