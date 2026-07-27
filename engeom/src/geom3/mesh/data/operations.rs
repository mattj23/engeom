//! This module contains the whole-mesh operations on `MeshData3`: rigid transformation, uniform
//! scaling, winding reversal, and appending one mesh onto another.
//!
//! Every one of these carries the per-element attributes along correctly, which is most of the
//! reason they belong on this type rather than being left to callers.
//!
//! # Orientation
//!
//! Two of these operations reverse the orientation of the surface, and both handle it rather than
//! leaving the caller with an inside-out mesh.
//!
//! A **negative** uniform scale is a mirror. Mirroring is orientation-reversing, so while the face
//! normals implied by the winding keep pointing the same way in space, the solid has been turned
//! inside out around them and they now point inward. `scale_in_place` therefore reverses the
//! winding and negates any stored point normals when the factor is negative.
//!
//! `flip_faces_in_place` reverses the winding deliberately, and negates stored point normals for the
//! same reason.

use super::MeshData3;
use crate::{Iso3, Point3, Result};

// ===============================================================================================
// Transformation
// ===============================================================================================

impl MeshData3 {
    /// Transform the mesh in place by a rigid isometry.
    ///
    /// The points are moved, and any stored point normals and per-point or per-face `Vector`
    /// attributes are rotated with them. Because those hold directions rather than positions, only
    /// the rotation component of the isometry affects them.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    pub fn transform_in_place(&mut self, iso: &Iso3) {
        for p in self.points.iter_mut() {
            *p = iso * *p;
        }

        self.attrs.transform_in_place(iso);
    }

    /// Return a copy of the mesh transformed by a rigid isometry, leaving this one unchanged.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    ///
    /// returns: `MeshData3`
    pub fn transform_copy(&self, iso: &Iso3) -> Self {
        let mut result = self.clone();
        result.transform_in_place(iso);
        result
    }

    /// Scale the mesh in place about the origin by a uniform factor.
    ///
    /// Any stored point standard deviations are scaled with the geometry, since they are lengths in
    /// the mesh's own units and would otherwise silently come to mean something else. Point normals
    /// are unaffected by the magnitude of a uniform scale, being directions.
    ///
    /// A **negative** factor mirrors the mesh, which reverses its orientation. The face winding is
    /// reversed and any stored point normals are negated to keep the surface facing the way it did
    /// before, rather than leaving it inside out.
    ///
    /// # Arguments
    ///
    /// * `scale`: the factor to scale by, which must be finite and non-zero
    ///
    /// returns: `Result<()>`
    pub fn scale_in_place(&mut self, scale: f64) -> Result<()> {
        if !scale.is_finite() || scale == 0.0 {
            return Err(format!(
                "A scale factor must be finite and non-zero, but {scale} was given. Scaling by zero \
                 would collapse the mesh to a point irrecoverably."
            )
            .into());
        }

        for p in self.points.iter_mut() {
            *p = Point3::from(p.coords * scale);
        }

        self.attrs.scale_in_place(scale);

        if scale < 0.0 {
            self.reverse_orientation();
        }

        Ok(())
    }

    /// Return a copy of the mesh scaled about the origin by a uniform factor, leaving this one
    /// unchanged.
    ///
    /// # Arguments
    ///
    /// * `scale`: the factor to scale by, which must be finite and non-zero
    ///
    /// returns: `Result<MeshData3>`
    pub fn scale_copy(&self, scale: f64) -> Result<Self> {
        let mut result = self.clone();
        result.scale_in_place(scale)?;
        Ok(result)
    }

    /// Reverse the winding order of every face, turning the surface inside out.
    ///
    /// Any stored point normals are negated to match, since the direction the surface faces has
    /// changed.
    pub fn flip_faces_in_place(&mut self) {
        self.reverse_orientation();
    }

    /// Reverse the winding of every face and negate any stored point normals.
    ///
    /// Swapping the first two indices of a triangle reverses the sign of its cross product and so
    /// the direction of the normal implied by its winding.
    fn reverse_orientation(&mut self) {
        for face in self.faces.iter_mut() {
            face.swap(0, 1);
        }

        self.attrs.flip_in_place();
    }
}

// ===============================================================================================
// Combination
// ===============================================================================================

impl MeshData3 {
    /// Append another mesh onto the end of this one.
    ///
    /// The other mesh's points are added after this one's and its faces are re-indexed to match. No
    /// points are welded and no faces are merged, so a shared boundary between the two meshes stays
    /// as two coincident sets of points.
    ///
    /// Attributes are all-or-nothing: a typed field or an open-map key present on one side and
    /// absent on the other is an error, because there is no correct value to pad the missing side
    /// with. The whole append is validated before anything is modified, so a failure leaves this
    /// mesh untouched.
    ///
    /// # Arguments
    ///
    /// * `other`: the mesh to append
    ///
    /// returns: `Result<()>`
    pub fn append_in_place(&mut self, other: &MeshData3) -> Result<()> {
        let offset = u32::try_from(self.points.len()).map_err(|_| {
            "The mesh already holds the maximum number of points a u32 can index".to_string()
        })?;

        u32::try_from(self.points.len() + other.points.len())
            .map_err(|_| "Appending would produce more points than a u32 can index".to_string())?;

        // Attributes are checked and merged first, because that is the step which can fail. Once it
        // succeeds, extending the buffers cannot.
        self.attrs.extend_from(&other.attrs)?;

        self.points.extend_from_slice(&other.points);
        self.faces.extend(
            other
                .faces
                .iter()
                .map(|f| [f[0] + offset, f[1] + offset, f[2] + offset]),
        );

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::mesh::data::Attr3;
    use crate::{UnitVec3, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::FRAC_PI_2;

    /// A single triangle in the xy plane, wound so that its normal is +z, carrying one attribute of
    /// every kind that an operation might need to touch.
    fn triangle() -> MeshData3 {
        let mut mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )
        .unwrap();

        mesh.set_point_normals(Some(vec![UnitVec3::new_normalize(Vector3::z()); 3]))
            .unwrap();
        mesh.set_point_stdev(Some(vec![0.1, 0.2, 0.3])).unwrap();
        mesh.set_face_labels(Some(vec![7])).unwrap();
        mesh.insert_point_attr("confidence", Attr3::Scalar(vec![0.5, 0.6, 0.7]))
            .unwrap();
        mesh.insert_point_attr("principal_dir", Attr3::Vector(vec![Vector3::x(); 3]))
            .unwrap();

        mesh
    }

    /// The unit normal implied by the winding of a face.
    fn winding_normal(mesh: &MeshData3, face: usize) -> Vector3 {
        let f = mesh.faces()[face];
        let p = mesh.points();
        let e1 = p[f[1] as usize] - p[f[0] as usize];
        let e2 = p[f[2] as usize] - p[f[0] as usize];
        e1.cross(&e2).normalize()
    }

    #[test]
    fn transform_moves_points_and_rotates_directions() {
        let mut mesh = triangle();

        // A quarter turn about +z, which maps +x onto +y, plus a translation.
        let iso = Iso3::new(Vector3::new(10.0, 0.0, 0.0), Vector3::z() * FRAC_PI_2);
        mesh.transform_in_place(&iso);

        // Points pick up both the rotation and the translation.
        assert_relative_eq!(
            mesh.points()[1],
            Point3::new(10.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );

        // Normals are directions, so only the rotation applies. +z is on the axis and is unmoved.
        assert_relative_eq!(
            mesh.point_normals().unwrap()[0].into_inner(),
            Vector3::z(),
            epsilon = 1.0e-12
        );

        // A Vector attribute rotates like a direction too.
        assert_relative_eq!(
            mesh.point_attr("principal_dir")
                .unwrap()
                .as_vector()
                .unwrap()[0],
            Vector3::y(),
            epsilon = 1.0e-12
        );

        // Scalars are untouched.
        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3]);
        assert_eq!(
            mesh.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.5, 0.6, 0.7]
        );
    }

    #[test]
    fn transform_copy_leaves_the_original_alone() {
        let mesh = triangle();
        let iso = Iso3::translation(5.0, 0.0, 0.0);

        let moved = mesh.transform_copy(&iso);

        assert_relative_eq!(mesh.points()[0], Point3::origin(), epsilon = 1.0e-12);
        assert_relative_eq!(
            moved.points()[0],
            Point3::new(5.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );
    }

    #[test]
    fn scale_scales_points_and_standard_deviations() -> Result<()> {
        let mut mesh = triangle();
        mesh.scale_in_place(25.4)?;

        assert_relative_eq!(
            mesh.points()[1],
            Point3::new(25.4, 0.0, 0.0),
            epsilon = 1.0e-12
        );

        for (actual, expected) in mesh
            .point_stdev()
            .unwrap()
            .iter()
            .zip([0.1, 0.2, 0.3].iter())
        {
            assert_relative_eq!(*actual, expected * 25.4, epsilon = 1.0e-12);
        }

        // Normals are directions and a positive uniform scale does not move them.
        assert_relative_eq!(
            mesh.point_normals().unwrap()[0].into_inner(),
            Vector3::z(),
            epsilon = 1.0e-12
        );

        // An open-map scalar has no declared units, so it must not be scaled.
        assert_eq!(
            mesh.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.5, 0.6, 0.7]
        );

        Ok(())
    }

    #[test]
    fn a_negative_scale_mirrors_and_keeps_the_surface_facing_outward() -> Result<()> {
        let mesh = triangle();
        assert_relative_eq!(winding_normal(&mesh, 0), Vector3::z(), epsilon = 1.0e-12);

        let mirrored = mesh.scale_copy(-1.0)?;

        // Every point is negated.
        assert_relative_eq!(
            mirrored.points()[1],
            Point3::new(-1.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );

        // Reflection through the origin is orientation-reversing, so without the winding fix the
        // face normal would still read +z while the solid had turned inside out. The reversal makes
        // it -z, matching the mirrored geometry.
        assert_relative_eq!(
            winding_normal(&mirrored, 0),
            -Vector3::z(),
            epsilon = 1.0e-12
        );

        // Stored normals follow.
        assert_relative_eq!(
            mirrored.point_normals().unwrap()[0].into_inner(),
            -Vector3::z(),
            epsilon = 1.0e-12
        );

        // Standard deviations stay positive despite the negative factor.
        for s in mirrored.point_stdev().unwrap() {
            assert!(*s > 0.0, "a standard deviation must not go negative");
        }

        Ok(())
    }

    #[test]
    fn scale_rejects_zero_and_non_finite_factors() {
        let mut mesh = triangle();

        assert!(mesh.scale_in_place(0.0).is_err());
        assert!(mesh.scale_in_place(f64::NAN).is_err());
        assert!(mesh.scale_in_place(f64::INFINITY).is_err());

        // The rejected calls must have left the mesh alone.
        assert_relative_eq!(
            mesh.points()[1],
            Point3::new(1.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );
    }

    #[test]
    fn flip_faces_reverses_winding_and_normals() {
        let mut mesh = triangle();
        let before = mesh.faces()[0];

        mesh.flip_faces_in_place();

        assert_eq!(mesh.faces()[0], [before[1], before[0], before[2]]);
        assert_relative_eq!(winding_normal(&mesh, 0), -Vector3::z(), epsilon = 1.0e-12);
        assert_relative_eq!(
            mesh.point_normals().unwrap()[0].into_inner(),
            -Vector3::z(),
            epsilon = 1.0e-12
        );

        // Flipping twice returns to the original state.
        mesh.flip_faces_in_place();
        assert_eq!(mesh.faces()[0], before);
        assert_relative_eq!(winding_normal(&mesh, 0), Vector3::z(), epsilon = 1.0e-12);
    }

    #[test]
    fn append_offsets_the_face_indices_and_unions_attributes() -> Result<()> {
        let mut mesh = triangle();
        let other = triangle().transform_copy(&Iso3::translation(10.0, 0.0, 0.0));

        mesh.append_in_place(&other)?;

        assert_eq!(mesh.point_count(), 6);
        assert_eq!(mesh.face_count(), 2);

        // The appended face must point at the appended points, not the originals.
        assert_eq!(mesh.faces()[0], [0, 1, 2]);
        assert_eq!(mesh.faces()[1], [3, 4, 5]);
        assert_relative_eq!(
            mesh.points()[3],
            Point3::new(10.0, 0.0, 0.0),
            epsilon = 1.0e-12
        );

        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.1, 0.2, 0.3]);
        assert_eq!(mesh.face_labels().unwrap(), &[7, 7]);
        assert_eq!(mesh.point_attr("confidence").unwrap().len(), 6);

        mesh.attrs()
            .validate(mesh.point_count(), mesh.face_count())?;

        Ok(())
    }

    #[test]
    fn append_rejects_mismatched_attributes_without_modifying_the_target() -> Result<()> {
        let mut mesh = triangle();
        let mut other = triangle();
        other.set_point_stdev(None)?;

        assert!(mesh.append_in_place(&other).is_err());

        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 1);
        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3]);

        Ok(())
    }

    #[test]
    fn appending_a_bare_mesh_onto_a_bare_mesh_works() -> Result<()> {
        let mut mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )?;
        let other = mesh.clone();

        mesh.append_in_place(&other)?;

        assert_eq!(mesh.point_count(), 6);
        assert_eq!(mesh.faces()[1], [3, 4, 5]);

        Ok(())
    }

    #[test]
    fn appending_an_empty_mesh_changes_nothing() -> Result<()> {
        let mut mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )?;

        mesh.append_in_place(&MeshData3::empty())?;

        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 1);

        Ok(())
    }
}
