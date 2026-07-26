//! This module computes normals from the point and face buffers of a triangle mesh.
//!
//! # Why angle weighting
//!
//! In previous versions of this library I used area weighting for computing normals. That turned
//! out to be a mistake, which is obvious as soon as you try to compute the vertex normals of a
//! box mesh.  At the time I didn't have the time or attention for a fix so it was left in the
//! library until version 0.4, when I had to deal with it again and took the time to find an
//! actual solution.
//!
//! It turns out that weighting by internal angle is the standard method (Thürmer and Wüthrich,
//! 1998), and this is what libigl, CGAL, and MeshLab use by default.
//!
//! Interior angle is invariant to how a planar polygon is triangulated, so regardless of how a
//! square face divided up, the angles its triangles contribute at a given corner sum to that
//! corner's 90 degrees. On the cube every face contributes an equal weight and the corner normal
//! lands exactly on the diagonal.
//!
//! # What this is not
//!
//! This produces one normal per point, which is a smoothed quantity. It is not the axis of the cone
//! bounding the incident face normals, which is what you want when you need a direction guaranteed
//! to face every adjacent face, and it is not the right displacement direction for offsetting a
//! surface, which needs the intersection of the offset planes rather than any averaged normal.

use crate::{Point3, Result, UnitVec3, Vector3};

/// Relative tolerance used to decide that a triangle is degenerate.
///
/// The test compares the magnitude of the edge cross product against the product of the edge
/// lengths, so both sides scale with the square of the mesh's units and the threshold means the same
/// thing whether the mesh is in meters or in microns.
const DEGENERATE_SIN: f64 = 1.0e-12;

/// Compute a normal for every point by averaging the normals of the faces which touch it, weighting
/// each by the interior angle of that face at that point.
///
/// Angle weighting is used because it is invariant to how a flat region is triangulated, which area
/// weighting and plain averaging are not. See the module documentation for the details.
///
/// Triangles with no well-defined normal, meaning those with a repeated point or with all three
/// points collinear, contribute nothing and are skipped.
///
/// # Arguments
///
/// * `points`: the point buffer
/// * `faces`: triangles given as triples of indices into `points`, all of which must be in range
///
/// returns: `Result<Vec<UnitVec3>>`, one normal per point, failing if any point has no well-defined
/// normal, which happens when it belongs to no face, belongs only to degenerate faces, or sits where
/// the weighted contributions cancel exactly
pub fn compute_point_normals(points: &[Point3], faces: &[[u32; 3]]) -> Result<Vec<UnitVec3>> {
    let mut sums = vec![Vector3::zeros(); points.len()];

    for face in faces.iter() {
        let p = [
            *points
                .get(face[0] as usize)
                .ok_or("Face refers to a point which does not exist")?,
            *points
                .get(face[1] as usize)
                .ok_or("Face refers to a point which does not exist")?,
            *points
                .get(face[2] as usize)
                .ok_or("Face refers to a point which does not exist")?,
        ];

        let Some(normal) = compute_face_normal(&p) else {
            continue;
        };

        for k in 0..3 {
            let e1 = p[(k + 1) % 3] - p[k];
            let e2 = p[(k + 2) % 3] - p[k];

            // atan2 of the cross magnitude against the dot is well conditioned across the whole
            // range, unlike acos of a normalized dot product near 0 and pi.
            let angle = e1.cross(&e2).norm().atan2(e1.dot(&e2));
            sums[face[k] as usize] += normal.into_inner() * angle;
        }
    }

    let mut normals = Vec::with_capacity(sums.len());
    for (i, sum) in sums.iter().enumerate() {
        let unit = UnitVec3::try_new(*sum, 1.0e-12).ok_or_else(|| {
            format!(
                "Point {i} has no well-defined normal. It belongs to no face, belongs only to \
                 degenerate faces, or its incident face normals cancel exactly."
            )
        })?;
        normals.push(unit);
    }

    Ok(normals)
}

/// Compute the unit normal of a single triangle, returning `None` if it is degenerate.
///
/// A triangle is degenerate when two of its points coincide or all three are collinear, in which
/// case it encloses no area and has no meaningful direction.
///
/// # Arguments
///
/// * `p`: the three points of the triangle, in winding order
///
/// returns: `Option<UnitVec3>`
pub fn compute_face_normal(p: &[Point3; 3]) -> Option<UnitVec3> {
    let e1 = p[1] - p[0];
    let e2 = p[2] - p[0];
    let cross = e1.cross(&e2);

    // Scale-invariant degeneracy test: |e1 x e2| = |e1||e2|sin(theta), so dividing through by the
    // edge lengths leaves a pure sine which does not depend on the units of the mesh.
    if cross.norm() <= DEGENERATE_SIN * e1.norm() * e2.norm() {
        return None;
    }

    Some(UnitVec3::new_normalize(cross))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// A cube of side 2 centred on the origin, with each square face split into two triangles.
    ///
    /// The diagonals are laid out so that the failure mode being guarded against is present: at
    /// point 1 the +z face contributes two triangles while the -y and -x faces contribute one each,
    /// which is what drives area weighting to the wrong answer.
    fn cube() -> (Vec<Point3>, Vec<[u32; 3]>) {
        let mut points = Vec::new();
        for x in [-1.0, 1.0] {
            for y in [-1.0, 1.0] {
                for z in [-1.0, 1.0] {
                    points.push(Point3::new(x, y, z));
                }
            }
        }

        // Each quad is given as four points in outward winding order, split on its first diagonal.
        let quads = [
            [1, 5, 7, 3], // +z
            [0, 2, 6, 4], // -z
            [0, 4, 5, 1], // -y
            [2, 3, 7, 6], // +y
            [4, 6, 7, 5], // +x
            [0, 1, 3, 2], // -x
        ];

        let mut faces = Vec::new();
        for q in quads {
            faces.push([q[0], q[1], q[2]]);
            faces.push([q[0], q[2], q[3]]);
        }

        (points, faces)
    }

    /// A grid of points in the xy plane, split into triangles.
    fn flat_grid(n: u32) -> (Vec<Point3>, Vec<[u32; 3]>) {
        let mut points = Vec::new();
        for j in 0..n {
            for i in 0..n {
                points.push(Point3::new(i as f64, j as f64, 0.0));
            }
        }

        let mut faces = Vec::new();
        for j in 0..n - 1 {
            for i in 0..n - 1 {
                let a = j * n + i;
                faces.push([a, a + 1, a + n + 1]);
                faces.push([a, a + n + 1, a + n]);
            }
        }

        (points, faces)
    }

    #[test]
    fn cube_corners_point_along_the_diagonals() -> Result<()> {
        let (points, faces) = cube();
        let normals = compute_point_normals(&points, &faces)?;

        assert_eq!(normals.len(), 8);

        for (normal, point) in normals.iter().zip(points.iter()) {
            let expected = UnitVec3::new_normalize(point.coords);
            assert_relative_eq!(
                normal.into_inner(),
                expected.into_inner(),
                epsilon = 1.0e-12
            );
        }

        Ok(())
    }

    /// The specific corner that area weighting gets wrong, pinned on its own so that the regression
    /// is unmistakable if the weighting is ever changed.
    #[test]
    fn the_corner_that_area_weighting_gets_wrong_is_exact() -> Result<()> {
        let (points, faces) = cube();
        let normals = compute_point_normals(&points, &faces)?;

        // Point 1 is at (-1, -1, 1) and touches two triangles of the +z face but only one each of
        // the -y and -x faces. Area weighting yields (-0.408, -0.408, 0.816), which is 19.5 degrees
        // off the diagonal.
        assert_eq!(points[1], Point3::new(-1.0, -1.0, 1.0));

        let d = 1.0 / 3.0_f64.sqrt();
        assert_relative_eq!(
            normals[1].into_inner(),
            Vector3::new(-d, -d, d),
            epsilon = 1.0e-12
        );

        Ok(())
    }

    #[test]
    fn a_flat_sheet_is_normal_everywhere() -> Result<()> {
        let (points, faces) = flat_grid(4);

        for normal in compute_point_normals(&points, &faces)? {
            assert_relative_eq!(normal.into_inner(), Vector3::z(), epsilon = 1.0e-12);
        }

        Ok(())
    }

    #[test]
    fn a_smooth_surface_lands_close_to_the_analytic_normal() -> Result<()> {
        // A unit sphere, whose true normal at any point is the radial direction.
        let sphere = crate::Mesh3::create_sphere(1.0, 64, 64);
        let normals = compute_point_normals(sphere.points(), sphere.faces())?;

        let mut worst: f64 = 0.0;
        for (normal, point) in normals.iter().zip(sphere.points()) {
            let expected = UnitVec3::new_normalize(point.coords);
            worst = worst.max(normal.angle(&expected));
        }

        assert!(
            worst.to_degrees() < 3.0,
            "worst deviation from the radial normal was {} degrees",
            worst.to_degrees()
        );

        Ok(())
    }

    #[test]
    fn degenerate_faces_are_skipped() -> Result<()> {
        // A flat quad plus one face with a repeated point, which has no normal and must not poison
        // the result.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let faces = vec![[0, 1, 2], [0, 2, 3], [0, 1, 1]];

        for normal in compute_point_normals(&points, &faces)? {
            assert_relative_eq!(normal.into_inner(), Vector3::z(), epsilon = 1.0e-12);
        }

        Ok(())
    }

    #[test]
    fn collinear_faces_are_skipped() -> Result<()> {
        // A 3x2 grid fully covered by good triangles, plus one extra face whose three points lie on
        // a line. The collinear face has no normal and must be skipped without disturbing the
        // points it touches.
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ];
        let faces = vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4], [0, 1, 2]];

        for normal in compute_point_normals(&points, &faces)? {
            assert_relative_eq!(normal.into_inner(), Vector3::z(), epsilon = 1.0e-12);
        }

        Ok(())
    }

    #[test]
    fn a_point_with_no_faces_is_an_error() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(5.0, 5.0, 5.0), // orphan
        ];
        let faces = vec![[0, 1, 2]];

        let err = compute_point_normals(&points, &faces)
            .unwrap_err()
            .to_string();
        assert!(err.contains("Point 3"), "{err}");
    }

    #[test]
    fn an_out_of_range_face_index_is_an_error() {
        let (points, mut faces) = flat_grid(3);
        faces.push([0, 1, 99]);

        assert!(compute_point_normals(&points, &faces).is_err());
    }

    #[test]
    fn scale_does_not_change_the_result() -> Result<()> {
        // The degeneracy test is relative, so a mesh in microns and the same mesh in meters must
        // produce identical normals.
        let (big, faces) = cube();
        let small: Vec<Point3> = big.iter().map(|p| *p * 1.0e-6).collect();

        for (a, b) in compute_point_normals(&big, &faces)?
            .iter()
            .zip(compute_point_normals(&small, &faces)?.iter())
        {
            assert_relative_eq!(a.into_inner(), b.into_inner(), epsilon = 1.0e-12);
        }

        Ok(())
    }

    #[test]
    fn face_normal_rejects_degenerate_triangles() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);

        // Repeated point.
        assert!(compute_face_normal(&[a, b, b]).is_none());

        // Collinear points.
        assert!(compute_face_normal(&[a, b, Point3::new(2.0, 0.0, 0.0)]).is_none());

        // A valid triangle still works.
        let n = compute_face_normal(&[a, b, Point3::new(0.0, 1.0, 0.0)]).unwrap();
        assert_relative_eq!(n.into_inner(), Vector3::z(), epsilon = 1.0e-12);
    }
}
