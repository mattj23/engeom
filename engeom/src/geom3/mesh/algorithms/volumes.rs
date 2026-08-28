//! This module selects the points and faces of a mesh which fall inside a bounding volume, from the
//! point and face buffers.
//!
//! The motivating use is cropping raw scan data down to a region of interest *before* it is
//! committed to an accelerated representation, so nothing here builds or consults a BVH. A caller
//! which already has one is not worse off: the work is a linear pass over buffers it already owns.
//!
//! # The two face semantics
//!
//! `compute_face_mask_in_aabb` takes an `all_vertices` flag which chooses between the two things a
//! caller can mean by a face being "in" the volume:
//!
//! - `true` selects the faces which lie **entirely** inside, meaning all three vertices are in.
//! - `false` selects the faces which **touch** the volume at all, including one which spans it with
//!   no vertex inside.
//!
//! These are the same two meanings the `all_vertices` flag carries on `TriangleFilter::above_plane`.
//! They coincide for a half-space, because a triangle overlaps a half-space exactly when one of its
//! vertices does, but for a bounded volume they do not, which is why the second case needs a real
//! triangle-box intersection test rather than a vertex count.
//!
//! # Why the vertex mask comes first
//!
//! The face pass starts by classifying the points, then counts how many of a face's three vertices
//! landed inside. That count settles the answer outright in every case but one:
//!
//! - three in: the face is inside under either meaning
//! - one or two in: the face touches the volume, and since the volume is convex there is nothing
//!   further to check
//! - none in: the face is out under the strict meaning, and is the *only* case where the loose
//!   meaning has to run a separating-axis test to find out
//!
//! A crop keeps or drops nearly every face outright, so the expensive test only runs on the thin
//! shell of triangles straddling the boundary.

use super::faces::face_points;
use crate::common::IndexMask;
use crate::geom3::Aabb3;
use crate::{Point3, Result};
use parry3d_f64::query::details::intersection_test_aabb_triangle;
use parry3d_f64::shape::Triangle;

/// Compute the mask of points which lie inside an axis-aligned bounding box.
///
/// A point exactly on a face of the box counts as inside.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `aabb`: the box to test against, in the same frame as the points
///
/// returns: `IndexMask` of length equal to the point count
pub fn compute_point_mask_in_aabb(points: &[Point3], aabb: &Aabb3) -> IndexMask {
    let mut mask = IndexMask::new(points.len(), false);
    for (i, point) in points.iter().enumerate() {
        mask.set(i, aabb.contains_local_point(point));
    }

    mask
}

/// Compute the mask of faces which lie inside an axis-aligned bounding box, under either of the two
/// meanings described in the module documentation.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`, all of which must be in range
/// * `aabb`: the box to test against, in the same frame as the points
/// * `all_vertices`: if `true`, a face is selected only when all three of its vertices are inside
///   the box; if `false`, a face is selected when any part of it touches the box at all
///
/// returns: `Result<IndexMask>` of length equal to the face count
pub fn compute_face_mask_in_aabb(
    points: &[Point3],
    faces: &[[u32; 3]],
    aabb: &Aabb3,
    all_vertices: bool,
) -> Result<IndexMask> {
    let inside = compute_point_mask_in_aabb(points, aabb);
    let mut mask = IndexMask::new(faces.len(), false);

    for (i, face) in faces.iter().enumerate() {
        let corners = face_points(points, face, i)?;
        let count = face
            .iter()
            .filter(|index| inside.get(**index as usize))
            .count();

        let keep = match (count, all_vertices) {
            (3, _) => true,
            (_, true) => false,
            (0, false) => {
                let triangle = Triangle::new(corners[0], corners[1], corners[2]);
                intersection_test_aabb_triangle(aabb, &triangle)
            }
            (_, false) => true,
        };

        mask.set(i, keep);
    }

    Ok(mask)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The unit box at the origin, which every test here selects against.
    fn unit_box() -> Aabb3 {
        Aabb3::new(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0))
    }

    /// A single triangle, as its own point and face buffers.
    fn one_face(a: Point3, b: Point3, c: Point3) -> (Vec<Point3>, Vec<[u32; 3]>) {
        (vec![a, b, c], vec![[0, 1, 2]])
    }

    #[test]
    fn a_point_on_the_boundary_is_inside() {
        let points = vec![
            Point3::new(0.5, 0.5, 0.5),
            Point3::new(0.0, 0.5, 0.5),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(1.0 + 1e-9, 0.5, 0.5),
        ];

        let mask = compute_point_mask_in_aabb(&points, &unit_box());

        assert_eq!(mask.to_indices(), vec![0, 1, 2]);
    }

    #[test]
    fn a_face_entirely_inside_is_selected_either_way() -> Result<()> {
        let (points, faces) = one_face(
            Point3::new(0.1, 0.1, 0.5),
            Point3::new(0.9, 0.1, 0.5),
            Point3::new(0.5, 0.9, 0.5),
        );

        for all_vertices in [true, false] {
            let mask = compute_face_mask_in_aabb(&points, &faces, &unit_box(), all_vertices)?;
            assert_eq!(mask.to_indices(), vec![0]);
        }

        Ok(())
    }

    #[test]
    fn a_face_entirely_outside_is_selected_neither_way() -> Result<()> {
        let (points, faces) = one_face(
            Point3::new(5.0, 5.0, 5.0),
            Point3::new(6.0, 5.0, 5.0),
            Point3::new(5.0, 6.0, 5.0),
        );

        for all_vertices in [true, false] {
            let mask = compute_face_mask_in_aabb(&points, &faces, &unit_box(), all_vertices)?;
            assert!(mask.to_indices().is_empty());
        }

        Ok(())
    }

    /// A face with one vertex in the box is only whole-face-inside under the loose meaning.
    #[test]
    fn a_face_with_one_vertex_inside_splits_the_two_meanings() -> Result<()> {
        let (points, faces) = one_face(
            Point3::new(0.5, 0.5, 0.5),
            Point3::new(5.0, 0.5, 0.5),
            Point3::new(0.5, 5.0, 0.5),
        );

        let strict = compute_face_mask_in_aabb(&points, &faces, &unit_box(), true)?;
        let loose = compute_face_mask_in_aabb(&points, &faces, &unit_box(), false)?;

        assert!(strict.to_indices().is_empty());
        assert_eq!(loose.to_indices(), vec![0]);

        Ok(())
    }

    /// The case the separating-axis test exists for: a triangle which passes through the box
    /// without putting a single vertex in it. A vertex count alone would miss this face.
    #[test]
    fn a_face_straddling_the_box_with_no_vertex_inside_is_selected() -> Result<()> {
        let (points, faces) = one_face(
            Point3::new(-1.0, 0.5, 0.5),
            Point3::new(2.0, 0.5, 0.5),
            Point3::new(0.5, 0.5, 3.0),
        );

        let strict = compute_face_mask_in_aabb(&points, &faces, &unit_box(), true)?;
        let loose = compute_face_mask_in_aabb(&points, &faces, &unit_box(), false)?;

        assert!(strict.to_indices().is_empty());
        assert_eq!(loose.to_indices(), vec![0]);

        Ok(())
    }

    /// A triangle whose own bounding box overlaps, while the triangle itself does not. This is the
    /// converse of the case above, and is what makes the separating-axis test worth the call.
    #[test]
    fn a_face_which_only_its_bounding_box_overlaps_is_not_selected() -> Result<()> {
        let (points, faces) = one_face(
            Point3::new(0.5, 2.0, 0.5),
            Point3::new(2.0, 0.5, 0.5),
            Point3::new(2.0, 2.0, 0.5),
        );

        let loose = compute_face_mask_in_aabb(&points, &faces, &unit_box(), false)?;

        assert!(loose.to_indices().is_empty());

        Ok(())
    }

    #[test]
    fn an_out_of_range_face_index_is_an_error() {
        let (points, _) = one_face(
            Point3::new(0.1, 0.1, 0.5),
            Point3::new(0.9, 0.1, 0.5),
            Point3::new(0.5, 0.9, 0.5),
        );

        let faces = vec![[0, 1, 3]];

        assert!(compute_face_mask_in_aabb(&points, &faces, &unit_box(), false).is_err());
    }

    #[test]
    fn empty_buffers_produce_empty_masks() -> Result<()> {
        let mask = compute_point_mask_in_aabb(&[], &unit_box());
        assert_eq!(mask.len(), 0);

        let mask = compute_face_mask_in_aabb(&[], &[], &unit_box(), false)?;
        assert_eq!(mask.len(), 0);

        Ok(())
    }
}
