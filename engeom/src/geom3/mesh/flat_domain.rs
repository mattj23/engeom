//! Queries between a mesh surface and its flattened 2D chart.
//!
//! A mesh with the `point_flat` attribute has a planar position for every point. These positions
//! can come from a flattening such as boundary-first flattening
//! (`MeshEdges::boundary_first_flatten` in `conformal`) and lay out the mesh's faces as tiles in a
//! region of the plane. This module provides queries in both directions across that chart. A point
//! on or near the surface can be located in the flat domain together with its depth above the
//! surface. A point in the flat domain can be lifted back to its corresponding face, barycentric
//! coordinates, position, and normal on the surface. The flat domain supports measurements that
//! only make sense on a plane and provides the layout for rasters such as depth maps; the
//! [`RasterMapping`] produced here completes that mapping.
//!
//! This is not a texture mapping. The chart is in the mesh's own length units, it is expected to
//! have low distortion, and distances and directions measured in it are intended to be meaningful
//! on the surface.
//!
//! # The two halves
//!
//! The per-point coordinates are stored as an attribute
//! ([`PointAttrSet3::flat`](crate::PointAttrSet3::flat), accessed through `point_flat` and
//! `set_point_flat` on the containers), so they remain attached to the mesh through data round
//! trips, subsets, appends, scaling, and PLY. Queries require a 2D acceleration structure over the
//! flat triangles. It is built on demand as a [`FlatDomain`] that borrows the mesh, much as
//! `compute_index` builds a `CloudIndex3` over a cloud. Build it after finalizing the mesh, and
//! rebuild it if `point_flat` changes; the domain holds its own copy of the 2D positions and does
//! not observe later changes to the attribute.
//!
//! ```ignore
//! let flat = mesh.compute_edges()?.boundary_first_flatten()?;
//! mesh.set_point_flat(Some(flat))?;
//!
//! let domain = mesh.compute_flat_domain()?;
//! let (uv, depth) = domain.project_to_flat(&measured, 2.0, PI / 4.0, None).unwrap();
//! let back = domain.surface_at(&uv).unwrap();
//! let raster = domain.raster_mapping(0.1, 4);
//! ```
//!
//! # Self-overlapping charts
//!
//! Nothing here checks that the flat triangles do not overlap one another. A flattening of a
//! surface that is not a topological disk, or a poor flattening of one that is, can fold over
//! itself, and a flat point under a fold belongs to more than one face. `surface_at` returns
//! the first containing face reached by the traversal. On such a chart, the surface-to-flat
//! direction is well defined, but the reverse direction is ambiguous.
//!
//! # Converting from `UvMapping`
//!
//! Before 2026-08 this functionality was a `UvMapping` held as a special-cased field on `Mesh3`.
//! This section is for code written against that API.
//!
//! What changed:
//!
//! - The per-point 2D coordinates are now the `point_flat` attribute, set on a finished mesh with
//!   `set_point_flat` rather than passed into the constructor. They are validated against the
//!   point count at the time of the call, which the old path could not do because it built the
//!   mapping from the triangles *before* `merge_duplicates`/`delete_degenerate` renumbered them.
//! - The 2D acceleration structure and the queries are a borrowed [`FlatDomain`], built on demand
//!   with `mesh.compute_flat_domain()?`, instead of living inside the mesh.
//! - The coordinates now survive `to_data`/`from_data`, `scale_copy` (scaled by the magnitude of
//!   the factor), every subset, and PLY (`flat_x`/`flat_y` columns), and `append_in_place`
//!   concatenates them. Previously, each of these operations either dropped the mapping or refused
//!   to run.
//! - `MeshData3` and `PointCloud3` carry `point_flat` too. A cloud projected onto a flattened
//!   reference can hold each point's flat position alongside its depth.
//!
//! Construction:
//!
//! ```ignore
//! // before
//! let mesh = Mesh3::new_with_options(points, faces, solid, merge, delete, Some(uv))?;
//! // after
//! let mut mesh = Mesh3::new_with_options(points, faces, solid, merge, delete)?;
//! mesh.set_point_flat(Some(uv))?;
//! let domain = mesh.compute_flat_domain()?;
//! ```
//!
//! Every old call and its replacement:
//!
//! | old | new |
//! |---|---|
//! | `mesh.uv()`, `UvMapping::new(uv, faces)` | `let domain = mesh.compute_flat_domain()?` |
//! | `uv.point(face, bc)` | `domain.flat_at(face, bc)?` (a bad face index is now an `Err`, not a panic) |
//! | `uv.triangle(&p2)` giving `(face, bc)` | `domain.surface_at(&p2)` giving a `MeshSurfPoint`, which carries `face_index`, `bc`, and the 3D point and normal |
//! | `mesh.uv_to_3d(&p2)` | `domain.surface_at(&p2)` |
//! | `uv.at_closest_uv(&p2)` giving `(face, p2)` | `domain.surface_closest_to(&p2)?` giving a `MeshSurfPoint`; the snapped 2D point is `domain.flat_at(mp.face_index, mp.bc)?`. One difference: for a point *inside* the chart the old call snapped it to the nearest edge of its triangle (parry's non-solid projection); the new one returns the point itself, and only points outside the chart are moved |
//! | `mesh.project_to_uv(&p3)` | `domain.flat_closest_to(&p3)` |
//! | `mesh.uv_with_tol(&p3, max_dist, max_angle, transform)` giving `(uv, depth)` | `domain.project_to_flat(&p3, max_dist, max_angle, transform)`, same result and semantics |
//! | `uv.make_raster_mapping(px_size, padding)` | `domain.raster_mapping(px_size, padding)` |
//! | `uv.faces()` | `mesh.faces()`; the domain shares the mesh's faces |
//! | "Cannot append meshes with UV mappings" | appends work; both meshes must carry `point_flat` or neither, like any other attribute |
//!
//! The one behavioral difference to watch for: `FlatDomain` borrows the mesh, so it cannot be
//! stored in a struct alongside an owned `Mesh3` the way the old field was. Build it where the
//! queries happen.

use super::{Mesh3, MeshSurfPoint};
use crate::common::PCoords;
use crate::common::barycentric::barycentric;
use crate::geom2::Aabb2;
use crate::raster2::RasterMapping;
use crate::{Iso3, Point2, Point3, Result, SurfacePoint3, Vector2};
use parry2d_f64::partitioning::TraversalAction;
use parry2d_f64::query::{PointQuery, PointQueryWithLocation};
use parry2d_f64::shape::TriMesh;

/// A borrowed index over a mesh's flattened chart that provides queries between the surface and
/// the flat domain.
///
/// Build one with [`Mesh3::compute_flat_domain`]; see the module documentation for details.
pub struct FlatDomain<'a> {
    mesh: &'a Mesh3,
    flat: TriMesh,
}

impl Mesh3 {
    /// Build the flat-domain index over this mesh's `point_flat` attribute.
    ///
    /// This method copies the flat coordinates and faces into a 2D triangle mesh with its own
    /// acceleration structure. The result borrows the mesh and does not observe later changes to
    /// the attribute, so rebuild it after `set_point_flat`.
    ///
    /// returns: `Result<FlatDomain>`, failing if the mesh has no `point_flat` attribute
    pub fn compute_flat_domain(&self) -> Result<FlatDomain<'_>> {
        let flat = self.point_flat().ok_or(
            "The mesh has no point_flat attribute. Set one with set_point_flat, for example from \
             compute_edges()?.boundary_first_flatten()?, before building a flat domain.",
        )?;
        let flat = TriMesh::new(flat.to_vec(), self.faces().to_vec())?;
        Ok(FlatDomain { mesh: self, flat })
    }
}

impl<'a> FlatDomain<'a> {
    /// The mesh this domain was built over.
    pub fn mesh(&self) -> &'a Mesh3 {
        self.mesh
    }

    /// The bounding box of the flat domain.
    pub fn bounds(&self) -> Aabb2 {
        self.flat.local_aabb()
    }

    /// The flat position of a location on the surface given by a face and barycentric
    /// coordinates within it.
    ///
    /// # Arguments
    ///
    /// * `face_index`: the face the location lies on
    /// * `bc`: the barycentric coordinates of the location within the face, in vertex order
    ///
    /// returns: `Result<Point2>`, failing if the face index is out of range
    pub fn flat_at(&self, face_index: u32, bc: [f64; 3]) -> Result<Point2> {
        if face_index as usize >= self.mesh.face_count() {
            return Err(format!(
                "Face index {} is out of range for a mesh with {} faces",
                face_index,
                self.mesh.face_count()
            )
            .into());
        }

        let tri = self.flat.triangle(face_index);
        Ok(Point2::from(
            tri.a.coords * bc[0] + tri.b.coords * bc[1] + tri.c.coords * bc[2],
        ))
    }

    /// The flat position of the point on the surface closest to a 3D point, with no limit on how
    /// far away it may be. Use `project_to_flat` instead when points far from the surface or beyond
    /// its edge should be rejected.
    ///
    /// # Arguments
    ///
    /// * `point`: the 3D point to locate
    ///
    /// returns: `Point2`
    pub fn flat_closest_to(&self, point: &impl PCoords<3>) -> Point2 {
        let mp = self.mesh.surface_closest_to(point);
        self.flat_at(mp.face_index, mp.bc)
            .expect("the closest surface point of a mesh lies on one of its own faces")
    }

    /// Locate a 3D point in the flat domain, together with its signed depth above the surface,
    /// if and only if it projects onto the surface within the given tolerances.
    ///
    /// This is the query behind a depth map: each measured point lands at a flat position with a
    /// depth. The tolerances are those of [`Mesh3::project_with_tol`]. The distance bound rejects
    /// points too far from the surface to belong to it, and the angle bound rejects points whose
    /// projection lands on an edge or vertex at an angle that indicates they are beyond the edge of
    /// the mesh rather than above its face. The depth is the projection of the point onto the
    /// normal of the face it landed on: positive on the side the normal points to.
    ///
    /// # Arguments
    ///
    /// * `point`: the 3D point to locate
    /// * `max_dist`: the maximum distance from the point to its projection on the surface
    /// * `max_angle`: the maximum angle between the face normal at the projection and the vector
    ///   from the projection to the point
    /// * `transform`: an optional transform applied to the point before projecting it
    ///
    /// returns: `Option<(Point2, f64)>`, the flat position and the signed depth, or `None` if the
    /// point fails either tolerance or lands on a face with no computable normal
    pub fn project_to_flat(
        &self,
        point: &impl PCoords<3>,
        max_dist: f64,
        max_angle: f64,
        transform: Option<&Iso3>,
    ) -> Option<(Point2, f64)> {
        let point = Point3::from(point.coords());
        let point = transform.map_or(point, |t| t * point);

        let (prj, face_index, loc) = self
            .mesh
            .project_with_tol(&point, max_dist, max_angle, None)?;
        let bc = loc.barycentric_coordinates()?;
        let normal = self.mesh.shape.triangle(face_index).normal()?;

        let flat = self
            .flat_at(face_index, bc)
            .expect("project_with_tol returns a face index of this mesh");
        let depth = SurfacePoint3::new(prj.point, normal).scalar_projection(&point);
        Some((flat, depth))
    }

    /// Lift a point of the flat domain back to the surface, if it lies within the chart.
    ///
    /// This is a containment test over the flat triangles, not a nearest-point search: a point
    /// outside every triangle gives `None`. Use `surface_closest_to` to snap such points onto the
    /// edge of the chart.
    ///
    /// # Arguments
    ///
    /// * `flat`: the point in the flat domain
    ///
    /// returns: `Option<MeshSurfPoint>`, carrying the face, barycentric coordinates, and the 3D
    /// point and normal, or `None` if the point is outside the chart or the face it lands on has
    /// no computable normal
    pub fn surface_at(&self, flat: &impl PCoords<2>) -> Option<MeshSurfPoint> {
        let (face_index, bc) = self.containing_face(&Point2::from(flat.coords()))?;
        self.mesh.at_barycentric(face_index, bc).ok()
    }

    /// Lift a point of the flat domain to the surface location whose flat position is nearest to
    /// it. For a point inside the chart this is the same as `surface_at`; for a point outside, it
    /// is the nearest point on the chart's boundary.
    ///
    /// # Arguments
    ///
    /// * `flat`: the point in the flat domain
    ///
    /// returns: `Result<MeshSurfPoint>`, failing only if the face landed on has no computable
    /// normal
    pub fn surface_closest_to(&self, flat: &impl PCoords<2>) -> Result<MeshSurfPoint> {
        // Solid, so that a point inside a flat triangle projects onto itself rather than onto the
        // nearest edge of that triangle.
        let p = Point2::from(flat.coords());
        let (prj, (face_index, _)) = self.flat.project_local_point_and_get_location(&p, true);
        let tri = self.flat.triangle(face_index);
        let bc = barycentric(&tri.a, &tri.b, &tri.c, &prj.point);
        self.mesh.at_barycentric(face_index, bc)
    }

    /// A raster mapping covering the flat domain, with the given pixel size and a margin of
    /// `padding` pixels on every side.
    ///
    /// The mapping's origin is the minimum corner of the domain's bounding box moved out by the
    /// padding, so that the whole chart lands inside the raster with a border around it.
    ///
    /// # Arguments
    ///
    /// * `px_size`: the size of each pixel in the mesh's length units
    /// * `padding`: the number of pixels of margin on each side
    ///
    /// returns: `RasterMapping`
    pub fn raster_mapping(&self, px_size: f64, padding: usize) -> RasterMapping {
        let bounds = self.bounds();
        let size = bounds.maxs - bounds.mins;

        let width_px = (size.x / px_size).ceil() as usize + padding * 2;
        let height_px = (size.y / px_size).ceil() as usize + padding * 2;
        let pad = padding as f64 * px_size;
        let origin = bounds.mins - Vector2::new(pad, pad);
        RasterMapping::new(origin, (height_px, width_px), px_size, None)
    }

    /// The first flat triangle found to contain the point, with the point's barycentric
    /// coordinates in it.
    fn containing_face(&self, point: &Point2) -> Option<(u32, [f64; 3])> {
        let mut result = None;
        self.flat.bvh().traverse(|node| {
            if !node.aabb().contains_local_point(point) {
                return TraversalAction::Prune;
            }

            if let Some(index) = node.leaf_data() {
                let tri = self.flat.triangle(index);
                if tri.contains_local_point(point) {
                    result = Some((index, barycentric(&tri.a, &tri.b, &tri.c, point)));
                    return TraversalAction::EarlyExit;
                }
            }
            TraversalAction::Continue
        });

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use crate::common::points::dist;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    /// A unit square in the xy plane, split into two faces, whose chart is its own xy.
    fn square() -> Mesh3 {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let mut mesh = Mesh3::new(points, vec![[0, 1, 2], [0, 2, 3]], false);
        let flat = mesh
            .points()
            .iter()
            .map(|p| Point2::new(p.x, p.y))
            .collect();
        mesh.set_point_flat(Some(flat)).unwrap();
        mesh
    }

    /// Half of a cylinder of the given radius about the z axis, from theta = 0 to pi and z = 0 to
    /// height, with the isometric chart `(r * theta, z)`. Normals point outward.
    fn half_cylinder(radius: f64, height: f64, n_theta: usize, n_z: usize) -> Mesh3 {
        let mut points = Vec::new();
        let mut flat = Vec::new();
        for iz in 0..=n_z {
            let z = height * iz as f64 / n_z as f64;
            for it in 0..=n_theta {
                let theta = PI * it as f64 / n_theta as f64;
                points.push(Point3::new(radius * theta.cos(), radius * theta.sin(), z));
                flat.push(Point2::new(radius * theta, z));
            }
        }

        let row = (n_theta + 1) as u32;
        let mut faces = Vec::new();
        for iz in 0..n_z as u32 {
            for it in 0..n_theta as u32 {
                let a = iz * row + it;
                let b = a + 1;
                let c = a + row;
                let d = c + 1;
                faces.push([a, b, d]);
                faces.push([a, d, c]);
            }
        }

        let mut mesh = Mesh3::new(points, faces, false);
        mesh.set_point_flat(Some(flat)).unwrap();
        mesh
    }

    #[test]
    fn a_mesh_without_flat_coordinates_has_no_domain() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, false);
        let message = mesh.compute_flat_domain().err().unwrap().to_string();
        assert!(message.contains("point_flat"), "{message}");
    }

    #[test]
    fn flat_at_interpolates_the_chart_and_checks_the_face() -> Result<()> {
        let mesh = square();
        let domain = mesh.compute_flat_domain()?;

        let p = domain.flat_at(0, [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])?;
        assert_relative_eq!(p, Point2::new(2.0 / 3.0, 1.0 / 3.0), epsilon = 1e-12);

        let p = domain.flat_at(1, [0.0, 0.0, 1.0])?;
        assert_relative_eq!(p, Point2::new(0.0, 1.0), epsilon = 1e-12);

        assert!(domain.flat_at(2, [1.0, 0.0, 0.0]).is_err());
        Ok(())
    }

    #[test]
    fn surface_at_finds_the_containing_face_and_misses_outside() -> Result<()> {
        let mesh = square();
        let domain = mesh.compute_flat_domain()?;

        // The centroid of the first face.
        let mp = domain
            .surface_at(&Point2::new(2.0 / 3.0, 1.0 / 3.0))
            .expect("inside the chart");
        assert_eq!(mp.face_index, 0);
        for w in mp.bc {
            assert_relative_eq!(w, 1.0 / 3.0, epsilon = 1e-9);
        }
        assert_relative_eq!(
            mp.sp.point,
            Point3::new(2.0 / 3.0, 1.0 / 3.0, 0.0),
            epsilon = 1e-9
        );
        assert_relative_eq!(mp.sp.normal.into_inner(), Vector3::z(), epsilon = 1e-12);

        // A point in the second face.
        let mp = domain
            .surface_at(&Point2::new(0.25, 0.75))
            .expect("inside the chart");
        assert_eq!(mp.face_index, 1);
        assert_relative_eq!(mp.sp.point, Point3::new(0.25, 0.75, 0.0), epsilon = 1e-9);

        assert!(domain.surface_at(&Point2::new(2.0, 2.0)).is_none());
        assert!(domain.surface_at(&Point2::new(-0.01, 0.5)).is_none());
        Ok(())
    }

    #[test]
    fn surface_closest_to_snaps_outside_points_to_the_chart_boundary() -> Result<()> {
        let mesh = square();
        let domain = mesh.compute_flat_domain()?;

        let mp = domain.surface_closest_to(&Point2::new(1.5, 0.5))?;
        assert_relative_eq!(mp.sp.point, Point3::new(1.0, 0.5, 0.0), epsilon = 1e-9);
        let snapped = domain.flat_at(mp.face_index, mp.bc)?;
        assert_relative_eq!(snapped, Point2::new(1.0, 0.5), epsilon = 1e-9);

        // Inside the chart it agrees with surface_at.
        let inside = Point2::new(0.3, 0.2);
        let a = domain.surface_closest_to(&inside)?;
        let b = domain.surface_at(&inside).unwrap();
        assert_eq!(a.face_index, b.face_index);
        assert_relative_eq!(a.sp.point, b.sp.point, epsilon = 1e-9);
        Ok(())
    }

    #[test]
    fn project_to_flat_reports_depth_and_honors_the_tolerances() -> Result<()> {
        let mesh = square();
        let domain = mesh.compute_flat_domain()?;

        let above = Point3::new(0.6, 0.2, 0.05);
        let (flat, depth) = domain
            .project_to_flat(&above, 0.1, PI / 4.0, None)
            .expect("within tolerance");
        assert_relative_eq!(flat, Point2::new(0.6, 0.2), epsilon = 1e-9);
        assert_relative_eq!(depth, 0.05, epsilon = 1e-12);

        let below = Point3::new(0.6, 0.2, -0.05);
        let (_, depth) = domain.project_to_flat(&below, 0.1, PI / 4.0, None).unwrap();
        assert_relative_eq!(depth, -0.05, epsilon = 1e-12);

        // Too far away.
        assert!(
            domain
                .project_to_flat(&Point3::new(0.6, 0.2, 0.5), 0.1, PI / 4.0, None)
                .is_none()
        );

        // Beyond the edge of the mesh: the projection lands on an edge at a right angle to the
        // normal.
        assert!(
            domain
                .project_to_flat(&Point3::new(1.05, 0.5, 0.0), 0.1, PI / 4.0, None)
                .is_none()
        );

        // A transform is applied to the point before projecting.
        let shift = Iso3::translation(0.0, 0.0, -0.05);
        let (_, depth) = domain
            .project_to_flat(&above, 0.1, PI / 4.0, Some(&shift))
            .unwrap();
        assert_relative_eq!(depth, 0.0, epsilon = 1e-12);

        // The unbounded query does not care about any of this.
        let far = domain.flat_closest_to(&Point3::new(0.6, 0.2, 5.0));
        assert_relative_eq!(far, Point2::new(0.6, 0.2), epsilon = 1e-9);
        Ok(())
    }

    #[test]
    fn a_half_cylinder_round_trips_through_its_unrolled_chart() -> Result<()> {
        let radius = 10.0;
        let mesh = half_cylinder(radius, 20.0, 36, 8);
        let domain = mesh.compute_flat_domain()?;

        // Points a little outside the cylinder wall at assorted angles and heights.
        for (theta, z, offset) in [
            (0.3, 5.0, 0.2),
            (1.2, 12.5, -0.15),
            (2.0, 1.0, 0.05),
            (2.9, 18.0, 0.3),
        ] {
            let p = Point3::new(
                (radius + offset) * f64::cos(theta),
                (radius + offset) * f64::sin(theta),
                z,
            );
            let (flat, depth) = domain
                .project_to_flat(&p, 1.0, PI / 4.0, None)
                .expect("close to the wall");

            // The faceted wall sits inside the true cylinder, so the depth is the offset plus the
            // sagitta of the facet, and the chart position is the unrolled angle to within the
            // facet's chord error.
            assert!(
                (depth - offset).abs() < 0.05,
                "theta {theta}: depth {depth} vs offset {offset}"
            );
            assert!(
                (flat.x - radius * theta).abs() < 0.05 && (flat.y - z).abs() < 1e-9,
                "theta {theta}: flat {flat:?}"
            );

            // Lifting the flat point returns the foot of the projection.
            let back = domain.surface_at(&flat).expect("inside the chart");
            let foot = mesh.surface_closest_to(&p);
            assert!(dist(&back.sp.point, &foot.sp.point) < 1e-9);
        }

        let bounds = domain.bounds();
        assert_relative_eq!(bounds.mins, Point2::new(0.0, 0.0), epsilon = 1e-12);
        assert_relative_eq!(bounds.maxs, Point2::new(radius * PI, 20.0), epsilon = 1e-12);
        Ok(())
    }

    #[test]
    fn raster_mapping_covers_the_chart_with_padding() -> Result<()> {
        let mesh = square();
        let domain = mesh.compute_flat_domain()?;

        let mapping = domain.raster_mapping(0.25, 2);
        assert_eq!(mapping.shape(), (8, 8));
        assert_relative_eq!(mapping.origin(), Point2::new(-0.5, -0.5), epsilon = 1e-12);
        assert_relative_eq!(mapping.px_size(), 0.25, epsilon = 0.0);
        Ok(())
    }
}
