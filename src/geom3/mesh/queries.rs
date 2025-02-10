//! Distance queries and measurements on meshes

use super::Mesh;
use crate::common::indices::chained_indices;
use crate::{Curve3, Iso3, Plane3, Point3, SurfacePoint3};
use parry3d_f64::query::{IntersectResult, PointProjection, PointQueryWithLocation, SplitResult};
use parry3d_f64::shape::TrianglePointLocation;
use std::f64::consts::PI;

impl Mesh {
    pub fn surf_closest_to(&self, point: &Point3) -> SurfacePoint3 {
        let result = self
            .shape
            .project_local_point_and_get_location(point, self.is_solid);
        let (projection, (tri_id, _location)) = result;
        let triangle = self.shape.triangle(tri_id);
        let normal = triangle.normal().unwrap(); // When could this fail? On a degenerate tri?
        SurfacePoint3::new(projection.point, normal)
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
