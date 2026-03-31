//! Distance queries and measurements on meshes

use super::{Mesh, MeshSurfPoint};
use crate::common::PCoords;
use crate::common::indices::chained_indices;
use crate::common::points::dist;
use crate::{Curve3, Iso3, Plane3, Point3, Result, SurfacePoint3};
use parry3d_f64::query::{IntersectResult, PointProjection, PointQueryWithLocation, SplitResult};
use parry3d_f64::shape::TrianglePointLocation;
use std::f64::consts::PI;

impl Mesh {
    /// This is an extremely simple closest distance query that returns only the scalar distance
    /// from the mesh to the point. It does not return any information about the face or which
    /// side of the corresponding face normal the point is on.  It will always return a single
    /// zero or positive value, which is the distance from the point to its closest projection on
    /// the mesh.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to seek the closest distance to
    ///
    /// returns: f64
    pub fn distance_closest_to(&self, point: &impl PCoords<3>) -> f64 {
        let point = Point3::from(point.coords());
        let result = self
            .shape
            .project_local_point_and_get_location(&point, self.is_solid);
        let (projection, _) = result;
        dist(&point, &projection.point)
    }

    /// Get the point and normal of a position on the mesh given a face ID and the barycentric
    /// coordinates of interest within the face.
    ///
    /// If the face ID is invalid or the normal is invalid, an error is returned.
    ///
    /// # Arguments
    ///
    /// * `face_id`: The ID of the face to query.
    /// * `bc`: An array of barycentric coordinates [u, v, w] where u + v + w = 1.0.
    ///
    /// returns: Result<SurfacePoint<3>, Box<dyn Error, Global>>
    pub fn at_barycentric(&self, face_id: u32, bc: [f64; 3]) -> Result<MeshSurfPoint> {
        if face_id >= self.faces().len() as u32 {
            return Err("Invalid face ID".into());
        }

        let face = self.shape.triangle(face_id);
        let coords = face.a.coords * bc[0] + face.b.coords * bc[1] + face.c.coords * bc[2];
        let normal = face.normal().ok_or("No face normal found")?;
        let sp = SurfacePoint3::new(coords.into(), normal);
        Ok(MeshSurfPoint {
            face_index: face_id,
            bc,
            sp,
        })
    }

    /// Find the index of the face that is closest to the given point in local coordinates.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to seek the closest face to
    ///
    /// returns: u32
    pub fn face_closest_to(&self, point: &impl PCoords<3>) -> u32 {
        let point = Point3::from(point.coords());
        let result = self
            .shape
            .project_local_point_and_get_location(&point, self.is_solid);
        let (_, (tri_id, _)) = result;
        tri_id
    }

    /// Find the closest point on the mesh surface to the specified test point. This method will
    /// return a descriptor which includes the face index, barycentric coordinates, and a
    /// point/normal combination.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to seek the closest surface point to
    ///
    /// returns: MeshSurfPoint
    pub fn surf_closest_to(&self, point: &impl PCoords<3>) -> MeshSurfPoint {
        let point = Point3::from(point.coords());
        let result = self
            .shape
            .project_local_point_and_get_location(&point, self.is_solid);
        let (projection, (tri_id, location)) = result;
        let triangle = self.shape.triangle(tri_id);
        let normal = triangle.normal().expect("Triangle doesn't have a normal");
        let sp = SurfacePoint3::new(projection.point, normal);
        let bc = location
            .barycentric_coordinates()
            .expect("Barycentric coordinates should be valid");

        MeshSurfPoint {
            face_index: tri_id,
            bc,
            sp,
        }
    }

    /// Find the closest point on the mesh surface to the specified point.
    ///
    /// This returns only the projected position, not the face index, barycentric
    /// coordinates, or surface normal. If the point projects onto an edge or vertex,
    /// the returned position is the nearest point on the mesh surface at that location.
    ///
    /// # Arguments
    ///
    /// * `point` - The test point to project onto the mesh.
    ///
    /// # Returns
    ///
    /// The closest point on the mesh surface to `point`.
    pub fn point_closest_to(&self, point: &impl PCoords<3>) -> Point3 {
        let query = Point3::from(point.coords());
        let (result, _) = self
            .shape
            .project_local_point_and_get_location(&query, self.is_solid);
        result.point
    }

    pub fn project_with_max_dist(
        &self,
        point: &impl PCoords<3>,
        max_dist: f64,
    ) -> Option<(PointProjection, u32, TrianglePointLocation)> {
        // TODO: Needs test coverage
        let point = Point3::from(point.coords());
        self.shape
            .project_local_point_and_get_location_with_max_dist(&point, self.is_solid, max_dist)
            .map(|(prj, (id, loc))| (prj, id, loc))
    }

    /// Given a test point, return its projection onto the mesh *if and only if* it is within the
    /// given distance tolerance from the mesh and the angle between the normal of the triangle and
    /// the +/- vector from the triangle to the point is less than the given angle tolerance.
    ///
    /// When a test point projects onto to the face of a triangle, the vector from the triangle
    /// point to the test point will be parallel to the triangle normal, by definition.  The angle
    /// tolerance will come into effect when the test point projects to an edge or vertex.  This
    /// will happen occasionally when the test point is near an edge with two triangles that reflex
    /// away from the point, and it will happen when the test point is beyond the edge of the mesh.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to project onto the mesh
    /// * `max_dist`: the maximum search distance from the test point to find a projection
    /// * `max_angle`: the max allowable angle deviation between the mesh normal at the projection
    ///   and the vector from the projection to the test point
    /// * `transform`: an optional transform to apply to the test point before projecting it onto
    ///   the mesh
    ///
    /// returns: Option<(PointProjection, u32, TrianglePointLocation)>
    pub fn project_with_tol(
        &self,
        point: &impl PCoords<3>,
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Option<(PointProjection, u32, TrianglePointLocation)> {
        // TODO: Needs test coverage
        let point = Point3::from(point.coords());
        let point = if let Some(transform) = transform {
            transform * point
        } else {
            point
        };

        let result = self
            .shape
            .project_local_point_and_get_location_with_max_dist(&point, self.is_solid, max_dist);
        if let Some((prj, (id, loc))) = result {
            let local = point - prj.point;
            let triangle = self.shape.triangle(id);
            if let Some(normal) = triangle.normal() {
                let angle = normal.angle(&local).abs();
                if angle < max_angle || angle > PI - max_angle {
                    Some((prj, id, loc))
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }

    /// Return the indices of the points in the given list that project onto the mesh within the
    /// given distance tolerance and angle tolerance.  An optional transform can be provided to
    /// transform the points before projecting them onto the mesh.
    ///
    /// # Arguments
    ///
    /// * `points`:
    /// * `max_dist`:
    /// * `max_angle`:
    /// * `transform`:
    ///
    /// returns: Vec<usize, Global>
    pub fn indices_in_tol(
        &self,
        points: &[Point3],
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Vec<usize> {
        // TODO: Needs test coverage
        let mut result = Vec::new();
        for (i, point) in points.iter().enumerate() {
            if self
                .project_with_tol(point, max_dist, max_angle, transform)
                .is_some()
            {
                result.push(i);
            }
        }
        result
    }

    /// Split the mesh by a plane into the portions on the negative and positive sides.
    ///
    /// If the plane intersects the mesh, the result contains two new meshes:
    /// - `Pair(negative, positive)` when the mesh is cut into two parts
    /// - `Negative` when the entire mesh lies on the negative side of the plane
    /// - `Positive` when the entire mesh lies on the positive side of the plane
    ///
    /// The returned meshes are newly constructed from the split triangle meshes and are not
    /// solid.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane used to split the mesh
    ///
    /// # Returns
    ///
    /// A `SplitResult<Mesh>` describing how the mesh relates to the plane.
    pub fn split(&self, plane: &Plane3) -> SplitResult<Mesh> {
        let result = self.shape.local_split(&plane.normal, plane.d, 1.0e-6);
        match result {
            SplitResult::Pair(a, b) => {
                let mesh_a = Mesh::new_take_trimesh(a, false);
                let mesh_b = Mesh::new_take_trimesh(b, false);
                SplitResult::Pair(mesh_a, mesh_b)
            }
            SplitResult::Negative => SplitResult::Negative,
            SplitResult::Positive => SplitResult::Positive,
        }
    }

    /// Perform a section of the mesh with a plane, returning a list of `Curve3` objects that
    /// trace the intersection of the mesh with the plane.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane used to section the mesh
    /// * `tol`: the tolerance used for the intersection and for the construction of the curves,
    ///   defaults to 1.0e-6
    ///
    /// returns: Result<Vec<Curve3, Global>, Box<dyn Error, Global>>
    pub fn section(&self, plane: &Plane3, tol: Option<f64>) -> Result<Vec<Curve3>> {
        let tol = tol.unwrap_or(1.0e-6);
        let mut collected = Vec::new();
        let result = self
            .shape
            .intersection_with_local_plane(&plane.normal, plane.d, tol);

        if let IntersectResult::Intersect(pline) = result {
            let chains = chained_indices(pline.indices());
            for chain in chains.iter() {
                let points = chain
                    .iter()
                    .map(|&i| pline.vertices()[i as usize])
                    .collect::<Vec<_>>();
                if let Ok(curve) = Curve3::from_points(&points, tol) {
                    collected.push(curve);
                }
            }
        }

        Ok(collected)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Mesh;
    use crate::tests::stanford_bun_4;
    use approx::assert_relative_eq;
    use parry3d_f64::na::Vector3;

    #[test]
    fn distance_closest_to_unit_box_face() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let point = Point3::new(1.5, 0.0, 0.0);

        let distance = mesh.distance_closest_to(&point);
        let closest = mesh.point_closest_to(&point);

        assert_relative_eq!(distance, 1.0, epsilon = 1.0e-12);
        assert_relative_eq!(closest, Point3::new(0.5, 0.0, 0.0), epsilon = 1.0e-12);
    }

    #[test]
    fn distance_closest_to_unit_box_edge() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let point = Point3::new(1.5, 1.5, 0.0);

        let distance = mesh.distance_closest_to(&point);
        let closest = mesh.point_closest_to(&point);

        assert_relative_eq!(distance, 2.0_f64.sqrt(), epsilon = 1.0e-12);
        assert_relative_eq!(closest, Point3::new(0.5, 0.5, 0.0), epsilon = 1.0e-12);
    }

    #[test]
    fn distance_closest_to_unit_box_corner_vertex() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let point = Point3::new(1.5, 1.5, 1.5);

        let distance = mesh.distance_closest_to(&point);
        let closest = mesh.point_closest_to(&point);

        assert_relative_eq!(distance, 3.0_f64.sqrt(), epsilon = 1.0e-12);
        assert_relative_eq!(closest, Point3::new(0.5, 0.5, 0.5), epsilon = 1.0e-12);
    }

    #[test]
    fn closest_to_matches_brute_force_on_unit_box() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let eps = 1.0e-12;

        for (face_index, face) in mesh.faces().iter().enumerate() {
            let face_index = face_index as u32;
            let a = mesh.vertices()[face[0] as usize];
            let b = mesh.vertices()[face[1] as usize];
            let c = mesh.vertices()[face[2] as usize];

            let centroid = Point3::from((a.coords + b.coords + c.coords) / 3.0);
            let normal = mesh.tri_mesh().triangle(face_index).normal().unwrap();

            let edge1 = b - a;
            let edge2 = c - a;
            let tangent = if edge1.norm_squared() > eps {
                edge1.normalize()
            } else {
                edge2.normalize()
            };
            let bitangent = normal.cross(&tangent).normalize();

            // Sample points across the face plane and slightly off the surface on both sides.
            for u in [-0.40, -0.20, 0.0, 0.20, 0.40] {
                for v in [-0.40, -0.20, 0.0, 0.20, 0.40] {
                    for n in [-0.25, 0.0, 0.25] {
                        let query =
                            centroid + tangent * u + bitangent * v + normal.into_inner() * n;

                        let brute = brute_force_closest(&query, &mesh);

                        // Check the closest face
                        let closest_face = mesh.face_closest_to(&query);
                        assert!(
                            brute.iter().any(|(id, _)| *id == closest_face),
                            "face_closest_to returned face {} for query {:?}, but brute force found {:?}",
                            closest_face,
                            query,
                            brute.iter().map(|(id, _)| *id).collect::<Vec<_>>()
                        );

                        // Check the closest surface
                        let closest_surf = mesh.surf_closest_to(&query);
                        assert!(
                            brute.iter().any(|(id, _)| *id == closest_surf.face_index),
                            "surf_closest_to returned face {} for query {:?}, but brute force found {:?}",
                            closest_surf.face_index,
                            query,
                            brute.iter().map(|(id, _)| *id).collect::<Vec<_>>()
                        );
                        let cp = brute
                            .iter()
                            .find(|(id, _)| *id == closest_surf.face_index)
                            .unwrap()
                            .1;
                        assert_relative_eq!(cp, closest_surf.sp.point, epsilon = 1.0e-12);
                        let abc = mesh
                            .at_barycentric(closest_surf.face_index, closest_surf.bc)
                            .unwrap();
                        assert_relative_eq!(cp, abc.sp.point, epsilon = 1.0e-12);
                        assert_relative_eq!(
                            abc.sp.normal,
                            closest_surf.sp.normal,
                            epsilon = 1.0e-12
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn at_barycentric_matches_face_vertices_and_center() {
        let mesh = stanford_bun_4();

        for (face_id, face) in mesh.faces().iter().enumerate() {
            let face_id = face_id as u32;
            let a = mesh.vertices()[face[0] as usize];
            let b = mesh.vertices()[face[1] as usize];
            let c = mesh.vertices()[face[2] as usize];

            let expected_normal = mesh.tri_mesh().triangle(face_id).normal().unwrap();

            let vertex_cases = [
                ([1.0, 0.0, 0.0], a),
                ([0.0, 1.0, 0.0], b),
                ([0.0, 0.0, 1.0], c),
                (
                    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                    Point3::from((a.coords + b.coords + c.coords) / 3.0),
                ),
            ];

            for (bc, expected_point) in vertex_cases {
                let sp = mesh.at_barycentric(face_id, bc).unwrap();

                assert_eq!(sp.face_index, face_id);
                assert_relative_eq!(sp.bc[0], bc[0], epsilon = 1.0e-12);
                assert_relative_eq!(sp.bc[1], bc[1], epsilon = 1.0e-12);
                assert_relative_eq!(sp.bc[2], bc[2], epsilon = 1.0e-12);
                assert_relative_eq!(sp.point(), expected_point, epsilon = 1.0e-12);
                assert_relative_eq!(
                    sp.normal().into_inner(),
                    expected_normal.into_inner(),
                    epsilon = 1.0e-12
                );
            }
        }
    }

    #[test]
    fn check_brute_force_closest_edge() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let point = Point3::new(1.5, 1.5, 0.0);

        let results = brute_force_closest(&point, &mesh);
        assert_eq!(results.len(), 2);
        for (_, p) in results.iter() {
            assert_relative_eq!(*p, Point3::new(0.5, 0.5, 0.0), epsilon = 1.0e-12);
        }
    }

    #[test]
    fn check_brute_force_closest_face() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let point = Point3::new(0.33, 0.0, 0.51);

        let results = brute_force_closest(&point, &mesh);
        assert_eq!(results.len(), 1);
        for (_, p) in results.iter() {
            assert_relative_eq!(*p, Point3::new(0.33, 0.0, 0.5), epsilon = 1.0e-12);
        }
    }

    fn brute_force_closest(query_point: &Point3, mesh: &Mesh) -> Vec<(u32, Point3)> {
        let mut closest: Vec<(u32, Point3)> = Vec::new();
        let mut best_dist = f64::INFINITY;

        for (face_index, face) in mesh.faces().iter().enumerate() {
            let a = mesh.vertices()[face[0] as usize];
            let b = mesh.vertices()[face[1] as usize];
            let c = mesh.vertices()[face[2] as usize];

            let closest_point = manual_closest_point_on_triangle(a, b, c, *query_point);
            let dist = dist(&closest_point, query_point);

            if dist + f64::EPSILON < best_dist {
                best_dist = dist;
                closest.clear();
                closest.push((face_index as u32, closest_point));
            } else if (dist - best_dist).abs() <= f64::EPSILON {
                closest.push((face_index as u32, closest_point));
            }
        }

        closest
    }
    #[test]
    fn split_unit_box_through_center_yields_two_nonempty_parts() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let plane = Plane3::x_axis();

        match mesh.split(&plane) {
            SplitResult::Pair(negative, positive) => {
                assert!(!negative.faces().is_empty());
                assert!(!positive.faces().is_empty());

                for vertex in negative.vertices() {
                    assert!(vertex.x <= 1.0e-12);
                }
                for vertex in positive.vertices() {
                    assert!(vertex.x >= -1.0e-12);
                }
            }
            _ => panic!("expected mesh to split into two parts"),
        }
    }

    #[test]
    fn split_unit_box_with_plane_outside_returns_positive_side() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let plane = Plane3::new(Vector3::x_axis(), -2.0);

        match mesh.split(&plane) {
            SplitResult::Positive => {}
            _ => panic!("expected whole mesh on positive side"),
        }
    }

    #[test]
    fn split_unit_box_with_plane_outside_returns_negative_side() {
        let mesh = Mesh::create_box(1.0, 1.0, 1.0, false);
        let plane = Plane3::new(Vector3::x_axis(), 2.0);

        match mesh.split(&plane) {
            SplitResult::Negative => {}
            _ => panic!("expected whole mesh on negative side"),
        }
    }

    #[test]
    fn section_unit_cylinder_in_xz_plane_creates_one_circle_curve() {
        use std::f64::consts::TAU;

        let mesh = Mesh::create_cylinder(1.0, 2.0, 256);
        let plane = Plane3::y_axis();

        let curves = mesh.section(&plane, Some(1.0e-10)).unwrap();
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
    fn section_two_unit_cylinder_in_xy_plane_creates_two_circles_curves() {
        use std::f64::consts::TAU;

        let mut mesh = Mesh::create_cylinder(1.0, 2.0, 256);
        mesh.transform_by(&Iso3::translation(0.0, 0.0, -2.0));
        let mut m1 = Mesh::create_cylinder(1.0, 2.0, 256);
        m1.transform_by(&Iso3::translation(0.0, 0.0, 2.0));

        mesh.append(&mut m1).unwrap();

        let plane = Plane3::y_axis();

        let curves = mesh.section(&plane, Some(1.0e-10)).unwrap();
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

    fn manual_closest_point_on_triangle(a: Point3, b: Point3, c: Point3, p: Point3) -> Point3 {
        let ab = b - a;
        let ac = c - a;
        let ap = p - a;

        let d1 = ab.dot(&ap);
        let d2 = ac.dot(&ap);
        if d1 <= 0.0 && d2 <= 0.0 {
            return a;
        }

        let bp = p - b;
        let d3 = ab.dot(&bp);
        let d4 = ac.dot(&bp);
        if d3 >= 0.0 && d4 <= d3 {
            return b;
        }

        let vc = d1 * d4 - d3 * d2;
        if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
            let v = d1 / (d1 - d3);
            return a + ab * v;
        }

        let cp = p - c;
        let d5 = ab.dot(&cp);
        let d6 = ac.dot(&cp);
        if d6 >= 0.0 && d5 <= d6 {
            return c;
        }

        let vb = d5 * d2 - d1 * d6;
        if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
            let w = d2 / (d2 - d6);
            return a + ac * w;
        }

        let va = d3 * d6 - d5 * d4;
        if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
            let bc = c - b;
            let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            return b + bc * w;
        }

        let denom = 1.0 / (va + vb + vc);
        let v = vb * denom;
        let w = vc * denom;
        a + ab * v + ac * w
    }
}
