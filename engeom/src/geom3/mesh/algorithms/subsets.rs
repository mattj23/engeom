//! This module extracts a sub-mesh from the point and face buffers, given masks over both domains.
//!
//! Both mesh containers need this and neither should own it. `MeshData3` reaches it through
//! `extract_subset_faces`/`extract_subset_points`, and `Mesh3` through `extract_subset_faces`. The
//! attribute half is not here, because it is not a buffer operation: each container hands its
//! attribute set the same pair of masks and lets it select for itself.

use crate::common::IndexMask;
use crate::{Point3, Result};

/// Build the point and face buffers of a sub-mesh from a pair of masks.
///
/// The surviving points are renumbered into a compact buffer in their original order, and the
/// surviving faces are re-indexed to match, so the result is a valid mesh in its own right rather
/// than a selection over the original.
///
/// The two masks must be **consistent**: every point referenced by a selected face must itself be
/// selected. A caller derives one from the other with `unique_point_mask` or `unique_face_mask`
/// rather than supplying both independently, and passing an inconsistent pair is an error rather
/// than a silently corrupt mesh.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`
/// * `point_mask`: selects which points survive, and must match the point count
/// * `face_mask`: selects which faces survive, and must match the face count
///
/// returns: `Result<(Vec<Point3>, Vec<[u32; 3]>)>`
pub fn compact_by_masks(
    points: &[Point3],
    faces: &[[u32; 3]],
    point_mask: &IndexMask,
    face_mask: &IndexMask,
) -> Result<(Vec<Point3>, Vec<[u32; 3]>)> {
    if point_mask.len() != points.len() {
        return Err(format!(
            "A point mask of length {} does not match a mesh with {} points",
            point_mask.len(),
            points.len()
        )
        .into());
    }

    if face_mask.len() != faces.len() {
        return Err(format!(
            "A face mask of length {} does not match a mesh with {} faces",
            face_mask.len(),
            faces.len()
        )
        .into());
    }

    // Maps an old point index onto its position in the compacted buffer. Points which were not
    // selected keep the sentinel, and reaching one means the masks disagreed.
    let mut remap = vec![u32::MAX; points.len()];
    for (new, old) in point_mask.iter_true().enumerate() {
        remap[old] = new as u32;
    }

    let mut new_faces = Vec::with_capacity(face_mask.count_true());
    for f in face_mask.iter_true() {
        let face = faces[f];
        let mut mapped = [0u32; 3];
        for (slot, index) in mapped.iter_mut().zip(face.iter()) {
            let new = remap[*index as usize];
            if new == u32::MAX {
                return Err(format!(
                    "Face {f} refers to point {index}, which the point selection excludes"
                )
                .into());
            }
            *slot = new;
        }
        new_faces.push(mapped);
    }

    let new_points = point_mask.clone_indices_of(points)?;

    Ok((new_points, new_faces))
}

/// Given a mask over the faces, produce the mask over the points which those faces reference.
///
/// # Arguments
///
/// * `faces`: triangles given as triples of indices into a point buffer of `n_points`
/// * `face_mask`: a mask whose length must match the face count
/// * `n_points`: the point count, which sets the length of the returned mask
///
/// returns: `Result<IndexMask>`
pub fn unique_point_mask(
    faces: &[[u32; 3]],
    face_mask: &IndexMask,
    n_points: usize,
) -> Result<IndexMask> {
    if face_mask.len() != faces.len() {
        return Err(format!(
            "A face mask of length {} does not match a mesh with {} faces",
            face_mask.len(),
            faces.len()
        )
        .into());
    }

    let mut points = IndexMask::new(n_points, false);
    for f in face_mask.iter_true() {
        for index in faces[f] {
            points.set(index as usize, true);
        }
    }

    Ok(points)
}

/// Given a mask over the points, produce the mask over the faces whose three points all survive.
///
/// # Arguments
///
/// * `faces`: triangles given as triples of indices into a point buffer
/// * `point_mask`: a mask whose length must match the point count
///
/// returns: `Result<IndexMask>`
pub fn unique_face_mask(faces: &[[u32; 3]], point_mask: &IndexMask) -> Result<IndexMask> {
    let mut mask = IndexMask::new(faces.len(), false);
    for (i, face) in faces.iter().enumerate() {
        let keep = face
            .iter()
            .all(|index| (*index as usize) < point_mask.len() && point_mask.get(*index as usize));
        mask.set(i, keep);
    }

    Ok(mask)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two triangles sharing an edge, over four points.
    fn quad() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let faces = vec![[0, 1, 2], [0, 2, 3]];
        (points, faces)
    }

    #[test]
    fn keeping_one_face_drops_the_points_it_does_not_use() -> Result<()> {
        let (points, faces) = quad();
        let face_mask = IndexMask::try_from_indices(&[0], 2)?;
        let point_mask = unique_point_mask(&faces, &face_mask, points.len())?;

        assert_eq!(point_mask.to_indices(), vec![0, 1, 2]);

        let (new_points, new_faces) = compact_by_masks(&points, &faces, &point_mask, &face_mask)?;

        assert_eq!(new_points.len(), 3);
        assert_eq!(new_faces, vec![[0, 1, 2]]);

        Ok(())
    }

    /// The second face uses points 0, 2, and 3, which compact to 0, 1, and 2.
    #[test]
    fn the_surviving_faces_are_re_indexed() -> Result<()> {
        let (points, faces) = quad();
        let face_mask = IndexMask::try_from_indices(&[1], 2)?;
        let point_mask = unique_point_mask(&faces, &face_mask, points.len())?;

        let (new_points, new_faces) = compact_by_masks(&points, &faces, &point_mask, &face_mask)?;

        assert_eq!(new_points, vec![points[0], points[2], points[3]]);
        assert_eq!(new_faces, vec![[0, 1, 2]]);

        Ok(())
    }

    #[test]
    fn a_face_referring_to_an_excluded_point_is_an_error() -> Result<()> {
        let (points, faces) = quad();

        // Keep both faces but only the points of the first, which face 1 then dangles off.
        let face_mask = IndexMask::new(2, true);
        let point_mask = IndexMask::try_from_indices(&[0, 1, 2], 4)?;

        assert!(compact_by_masks(&points, &faces, &point_mask, &face_mask).is_err());

        Ok(())
    }

    #[test]
    fn a_mask_of_the_wrong_length_is_an_error() -> Result<()> {
        let (points, faces) = quad();

        assert!(
            compact_by_masks(
                &points,
                &faces,
                &IndexMask::new(3, true),
                &IndexMask::new(2, true)
            )
            .is_err()
        );
        assert!(
            compact_by_masks(
                &points,
                &faces,
                &IndexMask::new(4, true),
                &IndexMask::new(5, true)
            )
            .is_err()
        );
        assert!(unique_point_mask(&faces, &IndexMask::new(9, true), points.len()).is_err());

        Ok(())
    }

    #[test]
    fn a_face_survives_a_point_selection_only_if_all_three_points_do() -> Result<()> {
        let (_, faces) = quad();

        // Points 0, 1, 2 keep face 0 only; face 1 needs point 3.
        let point_mask = IndexMask::try_from_indices(&[0, 1, 2], 4)?;
        let face_mask = unique_face_mask(&faces, &point_mask)?;

        assert_eq!(face_mask.to_indices(), vec![0]);

        Ok(())
    }
}
