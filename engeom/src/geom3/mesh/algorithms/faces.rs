//! This module computes the per-face quantities of a triangle mesh from its point and face buffers:
//! normals, areas, and centroids.
//!
//! These are one-liners over a triangle, but they belong here rather than on a container because
//! both mesh types need them and because a caller who already holds buffers should not have to
//! build a mesh to ask.
//!
//! # Degenerate faces
//!
//! A triangle with a repeated point, or with all three points collinear, has no normal. Areas and
//! centroids are still well defined for it: the area is zero and the centroid is the mean of the
//! three points. So `compute_face_areas` and `compute_face_centers` are infallible while
//! `compute_face_normals` is not, and the error names the offending face rather than reporting
//! that something, somewhere, went wrong.

use super::normals::compute_face_normal;
use crate::{Point3, Result, UnitVec3};

/// Compute the unit normal of every face, from the winding of its three points.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`, all of which must be in range
///
/// returns: `Result<Vec<UnitVec3>>`, one normal per face, failing on a face which has no
/// well-defined normal
pub fn compute_face_normals(points: &[Point3], faces: &[[u32; 3]]) -> Result<Vec<UnitVec3>> {
    faces
        .iter()
        .enumerate()
        .map(|(i, face)| {
            let p = face_points(points, face, i)?;
            compute_face_normal(&p).ok_or_else(|| {
                format!("Face {i} has no normal, because its points are coincident or collinear")
                    .into()
            })
        })
        .collect()
}

/// Compute the area of every face.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`, all of which must be in range
///
/// returns: `Result<Vec<f64>>`, one area per face. A degenerate face has an area of zero rather
/// than being an error.
pub fn compute_face_areas(points: &[Point3], faces: &[[u32; 3]]) -> Result<Vec<f64>> {
    faces
        .iter()
        .enumerate()
        .map(|(i, face)| {
            let p = face_points(points, face, i)?;
            Ok((p[1] - p[0]).cross(&(p[2] - p[0])).norm() * 0.5)
        })
        .collect()
}

/// Compute the centroid of every face, which is the mean of its three points.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`, all of which must be in range
///
/// returns: `Result<Vec<Point3>>`, one centroid per face
pub fn compute_face_centers(points: &[Point3], faces: &[[u32; 3]]) -> Result<Vec<Point3>> {
    faces
        .iter()
        .enumerate()
        .map(|(i, face)| {
            let p = face_points(points, face, i)?;
            Ok(Point3::from(
                (p[0].coords + p[1].coords + p[2].coords) / 3.0,
            ))
        })
        .collect()
}

/// Look up the three points of a face, naming the face in the error if an index is out of range.
fn face_points(points: &[Point3], face: &[u32; 3], index: usize) -> Result<[Point3; 3]> {
    let mut out = [Point3::origin(); 3];
    for (slot, i) in out.iter_mut().zip(face.iter()) {
        *slot = *points.get(*i as usize).ok_or_else(|| {
            format!(
                "Face {index} refers to point {i}, but the mesh has only {} points",
                points.len()
            )
        })?;
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use approx::assert_relative_eq;

    /// A right triangle in the xy plane with legs of 3 and 4, wound so its normal is +z.
    fn triangle() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
            Point3::new(0.0, 4.0, 0.0),
        ];
        (points, vec![[0, 1, 2]])
    }

    #[test]
    fn a_right_triangle_has_the_expected_normal_area_and_center() -> Result<()> {
        let (points, faces) = triangle();

        assert_relative_eq!(
            compute_face_normals(&points, &faces)?[0].into_inner(),
            Vector3::z(),
            epsilon = 1.0e-12
        );
        assert_relative_eq!(
            compute_face_areas(&points, &faces)?[0],
            6.0,
            epsilon = 1.0e-12
        );
        assert_relative_eq!(
            compute_face_centers(&points, &faces)?[0],
            Point3::new(1.0, 4.0 / 3.0, 0.0),
            epsilon = 1.0e-12
        );

        Ok(())
    }

    #[test]
    fn reversing_the_winding_reverses_the_normal() -> Result<()> {
        let (points, _) = triangle();
        let normals = compute_face_normals(&points, &[[1, 0, 2]])?;

        assert_relative_eq!(normals[0].into_inner(), -Vector3::z(), epsilon = 1.0e-12);

        Ok(())
    }

    /// A degenerate face has no normal, but its area and centroid are still well defined.
    #[test]
    fn a_collinear_face_has_no_normal_but_does_have_an_area() -> Result<()> {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ];
        let faces = vec![[0, 1, 2]];

        assert!(compute_face_normals(&points, &faces).is_err());
        assert_relative_eq!(compute_face_areas(&points, &faces)?[0], 0.0);
        assert_relative_eq!(
            compute_face_centers(&points, &faces)?[0],
            Point3::new(1.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );

        Ok(())
    }

    #[test]
    fn an_out_of_range_index_is_an_error() {
        let (points, _) = triangle();
        let faces = vec![[0, 1, 7]];

        assert!(compute_face_normals(&points, &faces).is_err());
        assert!(compute_face_areas(&points, &faces).is_err());
        assert!(compute_face_centers(&points, &faces).is_err());
    }
}
