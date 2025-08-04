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
    /// This is an extremely simple closest distance query which returns only the scalar distance
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

    pub fn point_closest_to(&self, point: &Point3) -> Point3 {
        let (result, _) = self
            .shape
            .project_local_point_and_get_location(point, self.is_solid);
        result.point
    }

    pub fn project_with_max_dist(
        &self,
        point: &Point3,
        max_dist: f64,
    ) -> Option<(PointProjection, u32, TrianglePointLocation)> {
        self.shape
            .project_local_point_and_get_location_with_max_dist(point, self.is_solid, max_dist)
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
        point: &Point3,
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Option<(PointProjection, u32, TrianglePointLocation)> {
        let point = if let Some(transform) = transform {
            transform * point
        } else {
            *point
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
    /// * `plane`:
    /// * `tol`:
    ///
    /// returns: Result<Vec<Curve3, Global>, Box<dyn Error, Global>>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn section(&self, plane: &Plane3, tol: Option<f64>) -> crate::Result<Vec<Curve3>> {
        let tol = tol.unwrap_or(1.0e-6);
        let mut collected = Vec::new();
        let result = self
            .shape
            .intersection_with_local_plane(&plane.normal, plane.d, 1.0e-6);

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
