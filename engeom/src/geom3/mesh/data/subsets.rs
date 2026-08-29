//! This module extracts sub-meshes from a `MeshData3`, selecting either by face or by point.
//!
//! The two directions are deliberately asymmetric, because in each case one domain is chosen by the
//! caller and the other is derived from it:
//!
//! - `extract_subset_faces` takes a face mask. The points are derived, so only points referenced by
//!   a surviving face are kept and any orphan is dropped.
//! - `extract_subset_points` takes a point mask. The faces are derived, so only faces whose three
//!   points *all* survive are kept. The points are exactly what the caller asked for, including any
//!   which end up belonging to no face, since dropping a point the caller explicitly selected would
//!   be surprising.
//!
//! Either way the surviving points are renumbered into a compact buffer and the faces are
//! re-indexed to match, so the result is a valid mesh in its own right rather than a mask over the
//! original. Every attribute array is carried through in both domains.

use super::MeshData3;
use crate::Result;
use crate::common::IndexMask;
use crate::geom3::mesh::algorithms::subsets::{
    compact_by_masks, compute_unique_face_mask, compute_unique_point_mask,
};

// ===============================================================================================
// Deriving one domain's mask from the other
// ===============================================================================================

impl MeshData3 {
    /// Given a mask over the faces, produce the mask over the points which those faces reference.
    ///
    /// # Arguments
    ///
    /// * `face_mask`: a mask whose length must match the face count
    ///
    /// returns: `Result<IndexMask>` of length equal to the point count
    pub fn compute_unique_point_mask(&self, face_mask: &IndexMask) -> Result<IndexMask> {
        compute_unique_point_mask(&self.faces, face_mask, self.points.len())
    }

    /// Given a mask over the points, produce the mask over the faces which can survive it.
    ///
    /// A face survives only if **all three** of its points are selected. A face with only some of
    /// its points selected cannot be kept, because the ones which were not selected would no longer
    /// exist for it to refer to.
    ///
    /// # Arguments
    ///
    /// * `point_mask`: a mask whose length must match the point count
    ///
    /// returns: `Result<IndexMask>` of length equal to the face count
    pub fn compute_unique_face_mask(&self, point_mask: &IndexMask) -> Result<IndexMask> {
        self.check_point_mask(point_mask)?;
        compute_unique_face_mask(&self.faces, point_mask)
    }
}

// ===============================================================================================
// Extracting sub-meshes
// ===============================================================================================

impl MeshData3 {
    /// Create a new mesh from a selection of faces.
    ///
    /// Only points referenced by a surviving face are kept, so any point which the selection orphans
    /// is dropped. The surviving points are renumbered and the faces re-indexed to match. Attributes
    /// are carried through in both domains.
    ///
    /// # Arguments
    ///
    /// * `face_mask`: a mask whose length must match the face count
    ///
    /// returns: `Result<MeshData3>`
    pub fn extract_subset_faces(&self, face_mask: &IndexMask) -> Result<Self> {
        let point_mask = self.compute_unique_point_mask(face_mask)?;
        self.compact(&point_mask, face_mask)
    }

    /// Create a new mesh from a selection of points.
    ///
    /// Only faces whose three points all survive are kept. Every selected point is kept, including
    /// any which ends up belonging to no surviving face. The surviving points are renumbered and the
    /// faces re-indexed to match. Attributes are carried through in both domains.
    ///
    /// # Arguments
    ///
    /// * `point_mask`: a mask whose length must match the point count
    ///
    /// returns: `Result<MeshData3>`
    pub fn extract_subset_points(&self, point_mask: &IndexMask) -> Result<Self> {
        let face_mask = self.compute_unique_face_mask(point_mask)?;
        self.compact(point_mask, &face_mask)
    }

    /// Build a new mesh from a pair of masks which are already known to be consistent, meaning that
    /// every point referenced by a selected face is itself selected.
    fn compact(&self, point_mask: &IndexMask, face_mask: &IndexMask) -> Result<Self> {
        let (points, faces) = compact_by_masks(&self.points, &self.faces, point_mask, face_mask)?;
        let point_attrs = self.point_attrs.subset(point_mask)?;
        let face_attrs = self.face_attrs.subset(face_mask)?;

        Self::new_with_attrs(points, faces, point_attrs, face_attrs)
    }
}

// ===============================================================================================
// Mask validation
// ===============================================================================================

impl MeshData3 {
    /// Verify that a mask is the right length to select over the points.
    fn check_point_mask(&self, mask: &IndexMask) -> Result<()> {
        if mask.len() != self.points.len() {
            return Err(format!(
                "A point mask of length {} does not match a mesh with {} points",
                mask.len(),
                self.points.len()
            )
            .into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point3;
    use crate::geom3::mesh::data::Attr3;

    /// Two triangles sharing the edge 0-2, over four points, plus a fifth point belonging to no
    /// face at all. Every attribute domain is populated so that subsetting can be checked on both.
    fn fixture() -> MeshData3 {
        let mut mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(5.0, 5.0, 5.0), // belongs to no face
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
        .unwrap();

        mesh.set_point_stdev(Some(vec![0.0, 0.1, 0.2, 0.3, 0.4]))
            .unwrap();
        mesh.set_face_labels(Some(vec![10, 20])).unwrap();
        mesh.insert_point_attr("confidence", Attr3::Scalar(vec![0.0, 1.0, 2.0, 3.0, 4.0]))
            .unwrap();
        mesh.insert_face_attr("quality", Attr3::Scalar(vec![0.5, 0.75]))
            .unwrap();

        mesh
    }

    fn mask(len: usize, selected: &[usize]) -> IndexMask {
        IndexMask::try_from_indices(selected, len).unwrap()
    }

    #[test]
    fn point_mask_collects_the_referenced_points() -> Result<()> {
        let mesh = fixture();

        // Only the first face: points 0, 1, 2.
        let points = mesh.compute_unique_point_mask(&mask(2, &[0]))?;
        assert_eq!(points.to_indices(), vec![0, 1, 2]);

        // Both faces: everything except the orphan.
        let points = mesh.compute_unique_point_mask(&mask(2, &[0, 1]))?;
        assert_eq!(points.to_indices(), vec![0, 1, 2, 3]);

        // No faces: nothing.
        let points = mesh.compute_unique_point_mask(&mask(2, &[]))?;
        assert_eq!(points.count_true(), 0);

        Ok(())
    }

    #[test]
    fn face_mask_requires_every_point_of_the_face() -> Result<()> {
        let mesh = fixture();

        // Points 0, 1, 2 fully contain the first face but only part of the second.
        let faces = mesh.compute_unique_face_mask(&mask(5, &[0, 1, 2]))?;
        assert_eq!(faces.to_indices(), vec![0]);

        // Dropping a single shared point loses both faces.
        let faces = mesh.compute_unique_face_mask(&mask(5, &[1, 2, 3, 4]))?;
        assert_eq!(faces.count_true(), 0);

        // Everything selected keeps everything.
        let faces = mesh.compute_unique_face_mask(&mask(5, &[0, 1, 2, 3, 4]))?;
        assert_eq!(faces.to_indices(), vec![0, 1]);

        Ok(())
    }

    #[test]
    fn subset_by_faces_drops_orphans_and_renumbers() -> Result<()> {
        let mesh = fixture();
        let sub = mesh.extract_subset_faces(&mask(2, &[1]))?;

        // The second face is [0, 2, 3], so points 0, 2, 3 survive and become 0, 1, 2.
        assert_eq!(sub.point_count(), 3);
        assert_eq!(sub.face_count(), 1);
        assert_eq!(sub.faces(), &[[0, 1, 2]]);

        assert_eq!(sub.points()[0], Point3::new(0.0, 0.0, 0.0));
        assert_eq!(sub.points()[1], Point3::new(1.0, 1.0, 0.0));
        assert_eq!(sub.points()[2], Point3::new(0.0, 1.0, 0.0));

        // Attributes follow the same selection in both domains.
        assert_eq!(sub.point_stdev().unwrap(), &[0.0, 0.2, 0.3]);
        assert_eq!(sub.face_labels().unwrap(), &[20]);
        assert_eq!(
            sub.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.0, 2.0, 3.0]
        );
        assert_eq!(
            sub.face_attr("quality").unwrap().as_scalar().unwrap(),
            &[0.75]
        );

        Ok(())
    }

    #[test]
    fn subset_by_points_keeps_selected_orphans() -> Result<()> {
        let mesh = fixture();

        // Select the first face's points plus the orphan. The orphan belongs to no face but was
        // asked for, so it stays.
        let sub = mesh.extract_subset_points(&mask(5, &[0, 1, 2, 4]))?;

        assert_eq!(sub.point_count(), 4);
        assert_eq!(sub.face_count(), 1);
        assert_eq!(sub.faces(), &[[0, 1, 2]]);
        assert_eq!(sub.points()[3], Point3::new(5.0, 5.0, 5.0));

        assert_eq!(sub.point_stdev().unwrap(), &[0.0, 0.1, 0.2, 0.4]);
        assert_eq!(sub.face_labels().unwrap(), &[10]);

        Ok(())
    }

    #[test]
    fn subset_by_points_drops_partially_selected_faces() -> Result<()> {
        let mesh = fixture();

        // Point 2 is shared by both faces, so excluding it loses both.
        let sub = mesh.extract_subset_points(&mask(5, &[0, 1, 3]))?;

        assert_eq!(sub.point_count(), 3);
        assert_eq!(sub.face_count(), 0);
        assert!(sub.face_labels().unwrap().is_empty());

        sub.point_attrs().validate(sub.point_count())?;
        sub.face_attrs().validate(sub.face_count())?;

        Ok(())
    }

    #[test]
    fn selecting_everything_reproduces_the_mesh() -> Result<()> {
        let mesh = fixture();

        let by_faces = mesh.extract_subset_faces(&mask(2, &[0, 1]))?;
        // The orphan point is not referenced by any face, so a face-driven selection drops it.
        assert_eq!(by_faces.point_count(), 4);
        assert_eq!(by_faces.faces(), mesh.faces());

        let by_points = mesh.extract_subset_points(&mask(5, &[0, 1, 2, 3, 4]))?;
        assert_eq!(by_points.point_count(), 5);
        assert_eq!(by_points.faces(), mesh.faces());
        assert_eq!(by_points.points(), mesh.points());
        assert_eq!(by_points.point_stdev(), mesh.point_stdev());

        Ok(())
    }

    #[test]
    fn selecting_nothing_gives_an_empty_mesh() -> Result<()> {
        let mesh = fixture();

        let sub = mesh.extract_subset_faces(&mask(2, &[]))?;
        assert!(sub.is_empty());
        sub.point_attrs().validate(0)?;
        sub.face_attrs().validate(0)?;

        let sub = mesh.extract_subset_points(&mask(5, &[]))?;
        assert!(sub.is_empty());

        Ok(())
    }

    #[test]
    fn a_mask_of_the_wrong_length_is_rejected() {
        let mesh = fixture();

        // A point mask handed to a face-taking method, and vice versa.
        assert!(mesh.compute_unique_point_mask(&mask(5, &[0])).is_err());
        assert!(mesh.compute_unique_face_mask(&mask(2, &[0])).is_err());
        assert!(mesh.extract_subset_faces(&mask(5, &[0])).is_err());
        assert!(mesh.extract_subset_points(&mask(2, &[0])).is_err());
    }

    #[test]
    fn the_result_is_a_valid_mesh_in_its_own_right() -> Result<()> {
        let mesh = fixture();
        let sub = mesh.extract_subset_faces(&mask(2, &[0]))?;

        // Every face index must be in range for the compacted point buffer, which is exactly what
        // the constructor checks, so rebuilding from the parts must succeed.
        MeshData3::new(sub.points().to_vec(), sub.faces().to_vec())?;
        sub.point_attrs().validate(sub.point_count())?;
        sub.face_attrs().validate(sub.face_count())?;

        Ok(())
    }
}
