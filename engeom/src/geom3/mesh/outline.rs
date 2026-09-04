//! This module exists to help generate a visual outline of a mesh

use super::Mesh3;
use crate::common::points::{fill_gaps, mid_point};
use crate::geom3::mesh::algorithms;
use crate::geom3::mesh::edges::{edge_key, naive_edges, unique_edges};
use crate::{Point3, Result, UnitVec3, Vector3};
use parry3d_f64::query::{Ray, RayCast};
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
type CEdgeTypes = (HashSet<[u32; 2]>, HashMap<[u32; 2], [u32; 2]>);

/// The fraction of the mesh's bounding box diagonal by which a perspective outline displaces its
/// points off the surface.
///
/// Unlike the absolute constant used by the parallel path, this value is scale-free. A mesh modeled in
/// meters and the same mesh modeled in millimeters produce the same drawing. At this fraction a
/// part filling a two thousand pixel frame is displaced by about two hundredths of a pixel, which
/// no raster can show, while staying far above the rounding of a double at those coordinates.
const PERSPECTIVE_EPS_FRACTION: f64 = 1e-5;

/// Where an outline is being viewed from.
///
/// This type is private because the two public entry points name the two available view modes.
/// Add any future mode here.
#[derive(Debug, Clone, Copy)]
enum OutlineView {
    /// Every sight line is parallel to this direction, which points from the scene toward the
    /// viewer.
    Parallel(UnitVec3),

    /// Every sight line converges on this point.
    Perspective(Point3),
}

impl OutlineView {
    /// The unit direction from a point toward the viewer, or `None` when the point coincides with
    /// the eye and no direction is defined.
    fn toward_eye(&self, p: &Point3) -> Option<UnitVec3> {
        match self {
            Self::Parallel(d) => Some(*d),
            Self::Perspective(eye) => UnitVec3::try_new(eye - p, 1e-12),
        }
    }

    /// How far an occlusion ray cast from a point may travel before it passes the viewer.
    ///
    /// A parallel view has no viewer position to stop at, so it uses the same unbounded limit the
    /// original implementation did. A perspective view must stop at the eye: an unbounded ray
    /// continues past it and reports geometry behind the camera as an occluder, which on an
    /// interior view marks the whole outline hidden.
    fn ray_limit(&self, p: &Point3, eps: f64) -> f64 {
        match self {
            Self::Parallel(_) => f64::MAX,
            Self::Perspective(eye) => (eye - p).norm() - eps,
        }
    }
}

/// The resolved settings for one outline computation. The two public entry points differ only in
/// the values assigned here.
#[derive(Debug, Clone, Copy)]
struct OutlineParams {
    view: OutlineView,
    max_edge_length: f64,
    corner_angle: f64,

    /// The distance by which outline points are displaced off the surface, both along the point
    /// normals of the returned points and along the sight line of the ray origins.
    eps: f64,
}

impl Mesh3 {
    pub fn compute_visual_outline(
        &self,
        facing: UnitVec3,
        max_edge_length: f64,
        corner_angle: Option<f64>,
    ) -> Result<Vec<(Point3, Point3, u8)>> {
        self.compute_outline(OutlineParams {
            view: OutlineView::Parallel(facing),
            max_edge_length,
            corner_angle: corner_angle.unwrap_or(PI / 4.0 - 1e-2),
            // Kept as the literal it has always been so that every existing drawing is
            // unchanged by the refactor. The perspective path uses a scale free value instead.
            eps: 1e-2,
        })
    }

    /// Compute a technical line drawing of the mesh as seen from a single eye point.
    ///
    /// This is the perspective counterpart of `compute_visual_outline`. It returns the same kind
    /// of tagged world space segments, but every sight line converges on `eye` rather than
    /// running parallel, so a silhouette is found where the surface turns away from that point
    /// and a segment is hidden only when the mesh blocks the straight line between it and the
    /// eye.
    ///
    /// Segments longer than `max_edge_length` are subdivided before they are classified, so that
    /// a long edge which is only partly hidden is not classified as a whole. In a perspective
    /// view the useful length depends on the distance to the eye, and `Camera::project_outline`
    /// converts a tolerance in pixels into this length.
    ///
    /// Geometry behind the eye is returned rather than removed. The outline is world space
    /// geometry and is not clipped to any view; use `Camera::project_segment` or
    /// `Camera::project_polyline` to clip and project it.
    ///
    /// The returned points are displaced off the surface by a small fraction of the mesh's
    /// bounding box so that the occlusion test does not immediately strike the surface that
    /// contains the outline. The points are therefore near the mesh surface, rather than on it.
    ///
    /// # Arguments
    ///
    /// * `eye`: the point in world coordinates that every sight line converges on
    /// * `max_edge_length`: the greatest length of a returned segment, in model units, which must
    ///   be a finite number greater than zero
    /// * `corner_angle`: the smallest angle in radians between two adjacent faces for their
    ///   shared edge to be drawn as a corner. If `None`, a default of just under 45 degrees is
    ///   used.
    ///
    /// returns: Result<Vec<(Point3, Point3, u8)>>, a list of segments each tagged with 0 when the
    /// segment is visible from the eye and 1 when the mesh hides it
    pub fn compute_perspective_outline(
        &self,
        eye: &Point3,
        max_edge_length: f64,
        corner_angle: Option<f64>,
    ) -> Result<Vec<(Point3, Point3, u8)>> {
        self.compute_outline(OutlineParams {
            view: OutlineView::Perspective(*eye),
            max_edge_length,
            corner_angle: corner_angle.unwrap_or(PI / 4.0 - 1e-2),
            eps: self.aabb().extents().norm() * PERSPECTIVE_EPS_FRACTION,
        })
    }

    /// The shared implementation for both outline entry points. `params` supplies the three
    /// differences between parallel and perspective views: the sight line at a point, the maximum
    /// ray distance before it passes the viewer, and the displacement epsilon.
    fn compute_outline(&self, params: OutlineParams) -> Result<Vec<(Point3, Point3, u8)>> {
        // `fill_gaps` searches for a subdivision count satisfying `d / (n + 1) <= max_dist`,
        // which no finite count satisfies when the limit is zero or negative, so an unchecked
        // value here would hang rather than fail.
        if !params.max_edge_length.is_finite() || params.max_edge_length <= 0.0 {
            return Err("max_edge_length must be a finite number greater than zero".into());
        }
        let corner_angle = params.corner_angle;

        let (boundaries, mut corners) = self.classified_edge_types();
        // let mut working = KeyChainer::new();
        let mut working = Vec::new();

        for (i, indices) in self.shape.indices().iter().enumerate() {
            for (i0, i1) in [(0, 1), (1, 2), (2, 0)] {
                let k = edge_key(&[indices[i0], indices[i1]]);

                if boundaries.contains(&k) {
                    working.push(k);
                } else if let Some(corner) = corners.get_mut(&k) {
                    if corner[0] == u32::MAX {
                        corner[0] = i as u32;
                    } else {
                        corner[1] = i as u32;
                    }
                }
            }
        }

        // At this point, working contains boundary edges and corners contains corner face pairs
        // Now we need to process the corners
        for (key, corner) in corners.iter() {
            if corner[0] == u32::MAX || corner[1] == u32::MAX {
                continue;
            }

            let n0u = self.shape.triangle(corner[0]).normal();
            let n1u = self.shape.triangle(corner[1]).normal();

            if let (Some(n0), Some(n1)) = (n0u, n1u) {
                if n0.angle(&n1) > corner_angle {
                    // Is this a corner?
                    working.push(*key);
                } else {
                    // Both faces have to be tested against the same sight line, or the test is
                    // not a sign straddle at all. Taking a direction per face, at each face
                    // centroid, would let both faces read as front facing from their own vantage
                    // even where the silhouette passes between them, dropping edges out of the
                    // outline and flipping the classification as the eye moves slightly.
                    let m = mid_point(
                        &self.shape.vertices()[key[0] as usize],
                        &self.shape.vertices()[key[1] as usize],
                    );
                    let Some(toward_eye) = params.view.toward_eye(&m) else {
                        continue;
                    };

                    let f0 = toward_eye.dot(&n0);
                    let f1 = toward_eye.dot(&n1);
                    let f_max = f0.max(f1);
                    let f_min = f0.min(f1);

                    if f_max >= 0.0 && f_min < 0.0 {
                        // Is this a silhouette?
                        working.push(*key);
                    }
                }
            }
        }

        // The normal is only used to nudge the outline off the surface so it does not z-fight, so a
        // point whose normal is ambiguous gets no nudge rather than failing the whole outline.
        let point_normals =
            algorithms::compute_point_normals_where_defined(self.points(), self.faces())?;
        let nudge = |index: u32| {
            point_normals[index as usize]
                .map(|n| n.into_inner() * params.eps)
                .unwrap_or_else(Vector3::zeros)
        };

        let mut edges = Vec::new();
        for k in working {
            let k0 = k[0];
            let k1 = k[1];

            let p0: Point3 = self.shape.vertices()[k0 as usize] + nudge(k0);
            let p1: Point3 = self.shape.vertices()[k1 as usize] + nudge(k1);

            let points = fill_gaps(&[p0, p1], params.max_edge_length);

            for (p0, p1) in points.iter().zip(points.iter().skip(1)) {
                let m = mid_point(p0, p1);
                let Some(toward_eye) = params.view.toward_eye(&m) else {
                    edges.push((*p0, *p1, 0));
                    continue;
                };

                // Bounded at the viewer, so that nothing beyond it can be mistaken for an
                // occluder. For a parallel view there is no viewer to stop at and the limit is
                // the same unbounded one as before.
                let limit = params.view.ray_limit(&m, params.eps);
                if limit <= 0.0 {
                    // Nothing can fit between the point and the eye.
                    edges.push((*p0, *p1, 0));
                    continue;
                }

                let ray = Ray::new(m + toward_eye.into_inner() * params.eps, *toward_eye);

                if self.shape.intersects_local_ray(&ray, limit) {
                    edges.push((*p0, *p1, 1))
                } else {
                    edges.push((*p0, *p1, 0))
                }
            }
        }

        Ok(edges)
    }

    fn classified_edge_types(&self) -> CEdgeTypes {
        let naive = naive_edges(self.shape.indices());
        let unique = unique_edges(&naive);

        let mut boundaries = HashSet::new();
        let mut corners = HashMap::new();

        for (key, count) in unique {
            if count == 1 {
                boundaries.insert(key);
            } else if count == 2 {
                corners.insert(key, [u32::MAX, u32::MAX]);
            }
        }

        (boundaries, corners)
    }
}
