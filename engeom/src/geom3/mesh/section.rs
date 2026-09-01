//! Tools for computing a planar section through a mesh. This module does not use the built-in
//! parry3d implementation because that implementation can become stuck in an infinite loop under
//! certain conditions.
//!
//! The work is split into two crate-visible stages so that other mesh algorithms can use the raw
//! cut without paying for the curves: [`Mesh3::section_segments`] produces one [`TriIntr`] per
//! cut triangle, still tagged with the face it came from, and [`chain_segments`] joins those into
//! connected runs. [`Mesh3::section_with_plane`] is the public composition of the two.

use super::Mesh3;
use crate::geom2::CurveGroup2;
use crate::geom3::Aabb3;
use crate::geom3::mesh::edges::edge_key;
use crate::geom3::planar_map::vertex_length_weights;
use crate::{
    Curve3, CurveGroup3, IndexMask, Line3, PlanarMap, Plane3, PlaneFrame, Point3, Result, Vector3,
};
use parry3d_f64::partitioning::TraversalAction;
use parry3d_f64::shape::TriMesh;
use std::collections::{HashMap, HashSet};

/// The result of cutting a mesh with a plane: the section curves in world coordinates, together
/// with the [`PlanarMap`] that brings them into the plane's two-dimensional coordinates and maps
/// 2D results back into the world.
///
/// The map is chosen when the section is taken according to the requested [`PlaneFrame`], so
/// everything derived from the section in 2D shares one frame. Retaining the curves in 3D also
/// lets the caller query the section where it lies on the part.
#[derive(Clone)]
pub struct PlanarSection {
    /// The section curves in world coordinates, one member per loop or open strand.
    pub curves: CurveGroup3,

    /// The mapping between the world and the plane's 2D coordinate system.
    pub map: PlanarMap,
}

impl PlanarSection {
    /// Brings the section curves into the plane's two-dimensional coordinates. This convenience
    /// method calls `self.map.curve_group_to_2d(&self.curves)`, which is typically the next step
    /// after taking a section.
    ///
    /// returns: `Result<CurveGroup2, Box<dyn Error, Global>>`
    ///
    /// # Failure
    ///
    /// Fails if a member collapses under the projection. This cannot happen to curves produced by
    /// the cut itself, but may occur if the group is later modified or transformed.
    pub fn to_2d(&self) -> Result<CurveGroup2> {
        self.map.curve_group_to_2d(&self.curves)
    }

    /// Splits the section into its curves and its map.
    pub fn into_parts(self) -> (CurveGroup3, PlanarMap) {
        (self.curves, self.map)
    }
}

/// Builds the map carried by a section according to the requested frame.
fn frame_map(plane: &Plane3, frame: PlaneFrame, curves: &CurveGroup3) -> Result<PlanarMap> {
    match frame {
        PlaneFrame::Auto => Ok(PlanarMap::from_plane(plane)),
        PlaneFrame::Oriented { origin, x } => PlanarMap::from_plane_oriented(plane, &origin, &x),
        PlaneFrame::Svd => {
            let mut points = Vec::new();
            let mut weights = Vec::new();
            for curve in curves.curves() {
                points.extend_from_slice(curve.vertices());
                weights.extend(vertex_length_weights(curve));
            }
            PlanarMap::from_plane_svd(plane, &points, Some(&weights))
        }
    }
}

impl Mesh3 {
    /// Computes the intersection between this mesh and a plane, returning the resulting curves
    /// as a single [`CurveGroup3`] together with the [`PlanarMap`] that brings them into two
    /// dimensions.
    ///
    /// A planar section of a mesh naturally consists of several curves rather than one. A part with
    /// a hole produces an outer and an inner loop, while an open mesh produces unclosed strands.
    /// These curves move together as one rigid body, which is what the group represents.
    ///
    /// Closed loops come back with their first vertex repeated as the last, since a `Curve3` has
    /// no closed flag of its own. That is what lets [`PlanarMap::curve_group_to_2d`] and
    /// [`crate::common::To2D`] recover the closure, so do not strip it.
    ///
    /// Each curve is wound so that its 2D normal, after mapping into the plane, points in the same
    /// direction as the original triangle normals. This preserves the surface's inside/outside
    /// orientation.
    ///
    /// The plane fixes only the z axis of the 2D coordinate system; `frame` chooses where its
    /// origin sits and which way x points. See [`PlaneFrame`] for the choices. The map that
    /// results is returned with the curves so that everything derived from the section shares it.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane to cut with
    /// * `frame`: how the plane's 2D coordinate system is laid out
    /// * `curve_tol`: the chordal tolerance given to the resulting curves, which is also the
    ///   distance within which their vertices are de-duplicated. Defaults to `1e-6`.
    /// * `faces`: an optional mask limiting the cut to a subset of the faces. Must be the same
    ///   length as the mesh has faces.
    ///
    /// returns: `Result<PlanarSection, Box<dyn Error, Global>>`
    ///
    /// # Failure
    ///
    /// A plane that does not intersect the mesh, or that misses every selected face, is an
    /// error rather than an empty group. A group with no members has no bounding box and no
    /// closest point, so there is nothing meaningful to hand back; a caller who expects to miss
    /// should match on the error rather than counting members.
    ///
    /// A frame that cannot be built is also an error: an `Oriented` x direction along the plane
    /// normal, or an `Svd` whose direction of greatest extent is not in the plane.
    ///
    /// # Examples
    ///
    /// ```
    /// use approx::assert_relative_eq;
    /// use engeom::{Mesh3, Plane3, PlaneFrame, Point2};
    ///
    /// let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
    /// let section = mesh.section_with_plane(&Plane3::xy(), PlaneFrame::Auto, None, None).unwrap();
    ///
    /// // The 2 by 4 box sections into one closed loop, centered on the plane's origin.
    /// let flat = section.to_2d().unwrap();
    /// assert_eq!(flat.len(), 1);
    /// assert!(flat.curves()[0].is_closed());
    /// assert_relative_eq!(flat.curves()[0].length(), 12.0, epsilon = 1e-9);
    ///
    /// // The map takes 2D results back to where they lie on the part.
    /// let corner = section.map.point_to_3d(&Point2::new(1.0, 2.0));
    /// assert_relative_eq!(corner.z, 0.0, epsilon = 1e-12);
    /// ```
    pub fn section_with_plane(
        &self,
        plane: &Plane3,
        frame: PlaneFrame,
        curve_tol: Option<f64>,
        faces: Option<&IndexMask>,
    ) -> Result<PlanarSection> {
        let curve_tol = curve_tol.unwrap_or(1e-6);
        let segments = self.section_segments(plane, faces)?;

        let mut results = Vec::new();
        for group in chain_segments(&segments) {
            let mut curve_points = group.iter().map(|&i| segments[i].a).collect::<Vec<_>>();
            curve_points.push(segments[group[group.len() - 1]].b);
            results.push(Curve3::from_points(&curve_points, curve_tol)?);
        }

        // Report this here instead of leaving it to `CurveGroup3::new`, whose general message
        // about groups would not tell the caller that the plane missed the mesh.
        if results.is_empty() {
            return Err(if faces.is_some() {
                "the plane does not intersect any of the selected faces".into()
            } else {
                Box::<dyn std::error::Error>::from("the plane does not intersect the mesh")
            });
        }

        let curves = CurveGroup3::new(results)?;
        let map = frame_map(plane, frame, &curves)?;
        Ok(PlanarSection { curves, map })
    }

    /// Cuts the mesh with a plane and returns one segment per intersected triangle, in no
    /// particular order and without joining them into curves.
    ///
    /// Each segment is oriented so that, when transformed into the plane, its 2D normal points
    /// in the same direction as the triangle's normal. It also carries the index of the face it
    /// was cut from. Triangles lying in the plane produce nothing.
    ///
    /// # Arguments
    ///
    /// * `plane`: the plane to cut with
    /// * `faces`: an optional mask limiting the cut to a subset of the faces. Must be the same
    ///   length as the mesh has faces.
    ///
    /// returns: `Result<Vec<TriIntr, Global>, Box<dyn Error, Global>>`
    ///
    /// # Failure
    ///
    /// A mask whose length does not match the face count is an error. A plane that misses is
    /// not: the result is simply empty, and it is the caller's business whether that matters.
    pub(crate) fn section_segments(
        &self,
        plane: &Plane3,
        faces: Option<&IndexMask>,
    ) -> Result<Vec<TriIntr>> {
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
                    TriIntr::new(ab, bc, edge_key(&[ai, bi]), edge_key(&[bi, ci]), face_i)
                }
                (None, Some(bc), Some(ca)) => {
                    TriIntr::new(bc, ca, edge_key(&[bi, ci]), edge_key(&[ci, ai]), face_i)
                }
                (Some(ab), None, Some(ca)) => {
                    TriIntr::new(ab, ca, edge_key(&[ai, bi]), edge_key(&[ci, ai]), face_i)
                }
                _ => Err("Something went wrong with the intersection calculation")?,
            };

            if tri_n.cross(&data.direction()).dot(&plane.normal) < 0.0 {
                segments.push(data.reversed());
            } else {
                segments.push(data);
            }
        }

        Ok(segments)
    }
}

/// Joins section segments into connected runs by matching the mesh edges they cross, returning
/// the segment indices of each run in order from its start to its end.
///
/// Endpoints connect when they share an edge key, and only the `b` end of one segment may lead
/// into the `a` end of the next, which is what keeps a run consistently wound. A run ends at an
/// edge crossed by exactly one segment (an open boundary) or by more than two (a non-manifold
/// junction), or when it closes on itself. Every segment lands in exactly one run, even one that
/// connects to nothing.
///
/// # Arguments
///
/// * `segments`: the segments to chain, as produced by [`Mesh3::section_segments`]
///
/// returns: Vec<Vec<usize, Global>, Global>
pub(crate) fn chain_segments(segments: &[TriIntr]) -> Vec<Vec<usize>> {
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

    // Now we'll group the segments into connected curves. We'll work until
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
            && !is_loop(segments, &current)
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

    connected
}

fn is_loop(segments: &[TriIntr], working: &[usize]) -> bool {
    let start = segments[working[0]].a;
    let end = segments[working[working.len() - 1]].b;

    start == end
}

/// The intersection of a plane with one triangle: a segment from `a` to `b`, the keys of the
/// two mesh edges it crosses at each end, and the index of the triangle it came from.
#[derive(Debug, Clone, Copy)]
pub(crate) struct TriIntr {
    pub(crate) a: Point3,
    pub(crate) b: Point3,
    pub(crate) key_a: [u32; 2],
    pub(crate) key_b: [u32; 2],
    pub(crate) face: u32,
}

impl TriIntr {
    fn new(a: Point3, b: Point3, key_a: [u32; 2], key_b: [u32; 2], face: u32) -> Self {
        Self {
            a,
            b,
            key_a,
            key_b,
            face,
        }
    }

    fn reversed(&self) -> Self {
        Self::new(self.b, self.a, self.key_b, self.key_a, self.face)
    }

    fn direction(&self) -> Vector3 {
        (self.b - self.a).normalize()
    }

    /// The length of the segment.
    pub(crate) fn length(&self) -> f64 {
        (self.b - self.a).norm()
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

        let section = mesh.section_with_plane(&Plane3::xy(), PlaneFrame::Auto, None, None)?;

        assert_eq!(section.curves.len(), 1);

        Ok(())
    }

    /// A plane that misses is an error rather than an empty group, so a caller cannot mistake
    /// "no intersection" for a group it can query.
    #[test]
    fn a_plane_which_misses_the_mesh_is_an_error() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, false);
        let plane = Plane3::new(Vector3::z_axis(), 50.0);

        let err = match mesh.section_with_plane(&plane, PlaneFrame::Auto, None, None) {
            Err(e) => e.to_string(),
            Ok(s) => panic!("a plane 50 units away cut {} curves", s.curves.len()),
        };
        assert!(
            err.contains("does not intersect"),
            "unexpected message: {err}"
        );
    }

    /// A cut that would otherwise intersect the mesh is also an error when the face selection
    /// excludes every intersection.
    #[test]
    fn a_face_mask_which_selects_nothing_is_an_error() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, false);
        let mask = IndexMask::new(mesh.face_count(), false);

        let err = match mesh.section_with_plane(&Plane3::xy(), PlaneFrame::Auto, None, Some(&mask))
        {
            Err(e) => e.to_string(),
            Ok(s) => panic!("an empty face mask cut {} curves", s.curves.len()),
        };
        assert!(err.contains("selected faces"), "unexpected message: {err}");
    }

    #[test]
    fn unit_cylinder_in_xy_plane_creates_one_circle_curve() {
        let mesh = Mesh3::create_cylinder(1.0, 2.0, 1.0e-5).unwrap();
        let plane = Plane3::xy();

        let curves = mesh
            .section_with_plane(&plane, PlaneFrame::Auto, Some(1.0e-10), None)
            .unwrap()
            .curves;
        assert_eq!(curves.len(), 1);

        let curve = &curves.curves()[0];
        assert!(curve.count() >= 3);

        for vertex in curve.vertices() {
            assert_relative_eq!(vertex.z, 0.0, epsilon = 1.0e-12);

            let radius = (vertex.x * vertex.x + vertex.y * vertex.y).sqrt();
            // The tolerance must account for each cylinder face's diagonal. Its intersection with
            // z=0 lies halfway between the arc endpoints formed by vertices placed on the radius.
            assert_relative_eq!(radius, 1.0, epsilon = 1.0e-4);
        }

        assert_relative_eq!(curve.length(), TAU, epsilon = 1.0e-2);
    }

    #[test]
    fn two_unit_cylinders_in_xy_plane_create_two_circle_curves() {
        let mut mesh = Mesh3::create_cylinder(1.0, 2.0, 1.0e-5).unwrap();
        mesh.transform_in_place(&Iso3::translation(-2.0, 0.0, 0.0));
        let mut m1 = Mesh3::create_cylinder(1.0, 2.0, 1.0e-5).unwrap();
        m1.transform_in_place(&Iso3::translation(2.0, 0.0, 0.0));

        mesh.append_in_place(&m1).unwrap();

        let plane = Plane3::xy();

        let curves = mesh
            .section_with_plane(&plane, PlaneFrame::Auto, Some(1.0e-10), None)
            .unwrap()
            .curves;
        assert_eq!(curves.len(), 2);

        for curve in curves.curves().iter() {
            assert!(curve.count() >= 3);
            assert_relative_eq!(curve.length(), TAU, epsilon = 1.0e-2);

            let expected_center = if curve.vertices()[0].x > 0.0 {
                Point3::new(2.0, 0.0, 0.0)
            } else {
                Point3::new(-2.0, 0.0, 0.0)
            };
            for vertex in curve.vertices() {
                assert_relative_eq!(dist(&expected_center, vertex), 1.0, epsilon = 1.0e-4);
            }
        }
    }

    // ============================================================================================
    // Bringing the section into the plane
    // ============================================================================================

    /// Exercises the complete workflow: cut a mesh, bring the section into its plane, and recover
    /// a faithful 2D copy.
    #[test]
    fn a_mesh_section_projects_back_into_its_own_plane() {
        let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);

        let section = mesh
            .section_with_plane(&Plane3::xy(), PlaneFrame::Auto, Some(1e-9), None)
            .unwrap();
        let flat = section.to_2d().unwrap();

        assert_eq!(flat.len(), 1);

        let loop_2d = &flat.curves()[0];
        assert!(loop_2d.is_closed(), "a box sections into a closed loop");

        // The box is 2 by 4 in the cutting plane, so the loop is its perimeter.
        assert_relative_eq!(loop_2d.length(), 12.0, epsilon = 1e-6);

        let aabb = loop_2d.aabb();
        assert_relative_eq!(aabb.mins.x, -1.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.x, 1.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.mins.y, -2.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.y, 2.0, epsilon = 1e-9);
    }

    /// The same cut on a tilted plane must produce the same 2D shape because the frame
    /// follows the plane. This is what makes the projection independent of how the part is posed.
    #[test]
    fn a_tilted_section_has_the_same_shape_as_an_upright_one() {
        let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
        let upright = mesh
            .section_with_plane(&Plane3::xy(), PlaneFrame::Auto, Some(1e-9), None)
            .unwrap()
            .to_2d()
            .unwrap();

        // Rotate the part, and cut with the correspondingly rotated plane.
        let iso = Iso3::new(Vector3::new(10.0, -5.0, 2.0), Vector3::new(0.3, -0.5, 0.2));
        let moved = mesh.transform_copy(&iso);
        let plane = Plane3::xy().transformed_by(&iso);

        let tilted = moved
            .section_with_plane(&plane, PlaneFrame::Auto, Some(1e-9), None)
            .unwrap()
            .to_2d()
            .unwrap();

        assert_eq!(tilted.len(), upright.len());
        assert_relative_eq!(tilted.length(), upright.length(), epsilon = 1e-6);
    }

    /// The map goes both ways: a 2D result lifted back through it lands on the part.
    #[test]
    fn a_2d_result_lifts_back_onto_the_section() {
        let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
        let iso = Iso3::new(Vector3::new(1.0, 2.0, 3.0), Vector3::new(0.3, -0.5, 0.2));
        let moved = mesh.transform_copy(&iso);
        let plane = Plane3::xy().transformed_by(&iso);

        let section = moved
            .section_with_plane(&plane, PlaneFrame::Auto, Some(1e-9), None)
            .unwrap();
        let flat = section.to_2d().unwrap();
        let back = section.map.curve_group_to_3d(&flat).unwrap();

        assert_relative_eq!(back.length(), section.curves.length(), epsilon = 1e-9);
        let (_, station) = section
            .curves
            .at_closest_to_point(&back.curves()[0].vertices()[1]);
        assert_relative_eq!(
            dist(&station.point(), &back.curves()[0].vertices()[1]),
            0.0,
            epsilon = 1e-9
        );
    }

    /// An off-center box cut with the `Svd` frame returns centered with its long side along x,
    /// whichever way the part is posed.
    #[test]
    fn an_svd_frame_centers_the_section_on_its_long_axis() {
        let mesh = Mesh3::create_box(4.0, 2.0, 6.0, false);
        let iso = Iso3::new(Vector3::new(10.0, -5.0, 2.0), Vector3::new(0.3, -0.5, 0.2));
        let moved = mesh.transform_copy(&iso);
        // The plane is rotated with the part but its origin is not moved with it, so the Auto
        // frame would put the loop well off-center.
        let plane = Plane3::xy().transformed_by(&iso);

        let flat = moved
            .section_with_plane(&plane, PlaneFrame::Svd, Some(1e-9), None)
            .unwrap()
            .to_2d()
            .unwrap();

        let aabb = flat.curves()[0].aabb();
        assert_relative_eq!(aabb.center(), crate::Point2::origin(), epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.x - aabb.mins.x, 4.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.y - aabb.mins.y, 2.0, epsilon = 1e-9);
    }

    /// An `Oriented` frame places the origin and x axis at their requested locations after
    /// projection.
    #[test]
    fn an_oriented_frame_places_the_origin_and_x_axis() {
        let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
        let frame = PlaneFrame::Oriented {
            origin: Point3::new(-1.0, -2.0, 5.0),
            x: Vector3::new(0.0, 1.0, 0.5),
        };

        let flat = mesh
            .section_with_plane(&Plane3::xy(), frame, Some(1e-9), None)
            .unwrap()
            .to_2d()
            .unwrap();

        // The loop's corner at (-1, -2) is the origin, and the box's y extent is now along x.
        let aabb = flat.curves()[0].aabb();
        assert_relative_eq!(aabb.mins.x, 0.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.x, 4.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.mins.y, -2.0, epsilon = 1e-9);
        assert_relative_eq!(aabb.maxs.y, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn a_frame_which_cannot_be_built_is_an_error() {
        let mesh = Mesh3::create_box(2.0, 4.0, 6.0, false);
        let frame = PlaneFrame::Oriented {
            origin: Point3::origin(),
            x: Vector3::z(),
        };
        assert!(
            mesh.section_with_plane(&Plane3::xy(), frame, None, None)
                .is_err()
        );
    }
}
