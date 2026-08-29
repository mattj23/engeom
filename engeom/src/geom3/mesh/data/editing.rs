//! This module contains the element-wise editing operations on `MeshData3`: adding, modifying, and
//! removing individual points and faces. They maintain the invariant type guarantees, so a mesh
//! can never be observed with a face pointing at a point which does not exist or with an attribute
//! array of the wrong length.
//!
//! # Insertion and attributes
//!
//! For now, adding to a domain is blocked if there are attributes on that domain. My longer term
//! plan is to allow it if attribute key/value pairs are supplied, but I don't yet know what that
//! would look like.
//!
//! Removal needs no new values, so it always works.
//!
//! # Removal shifts indices
//!
//! Removal is ordered rather than a swap-remove, so the elements after the one removed shift down
//! by one. This will obviously invalidate any indices above the current one, so be advised that if
//! some code is holding on to an index after you do this it will be pointing at something else.
//!
//! In the case of points, whose indices are referenced in the faces, removing a point causes the
//! faces to all be updated.

use super::MeshData3;
use crate::common::IndexMask;
use crate::{Point3, Result};

// ===============================================================================================
// Adding elements
// ===============================================================================================

impl MeshData3 {
    /// Add a point to the end of the point buffer and return its index.
    ///
    /// The index is returned so that it can be fed directly into a face while assembling a mesh
    /// point by point.
    ///
    /// # Arguments
    ///
    /// * `point`: the position to add
    ///
    /// returns: `Result<u32>`, failing if the mesh carries per-point attributes, since there would
    /// be no value to give them for the new point
    pub fn push_point(&mut self, point: Point3) -> Result<u32> {
        let labels = self.point_attrs.attr_labels();
        if !labels.is_empty() {
            return Err(format!(
                "Cannot add a point to a mesh carrying per-point attributes ({}), because there is \
                 no value to give them for the new point. Remove them first, or build the point \
                 buffer before attaching attributes.",
                labels.join(", ")
            )
            .into());
        }

        let index = u32::try_from(self.points.len())
            .map_err(|_| "The mesh already holds the maximum number of points a u32 can index")?;
        self.points.push(point);
        Ok(index)
    }

    /// Add a face to the end of the face buffer and return its index.
    ///
    /// # Arguments
    ///
    /// * `face`: three indices into the point buffer, each of which must refer to a point which
    ///   exists
    ///
    /// returns: `Result<u32>`, failing if any index is out of range or if the mesh carries per-face
    /// attributes
    pub fn push_face(&mut self, face: [u32; 3]) -> Result<u32> {
        let labels = self.face_attrs.attr_labels();
        if !labels.is_empty() {
            return Err(format!(
                "Cannot add a face to a mesh carrying per-face attributes ({}), because there is no \
                 value to give them for the new face. Remove them first, or build the face buffer \
                 before attaching attributes.",
                labels.join(", ")
            )
            .into());
        }

        self.check_face(&face)?;

        let index = u32::try_from(self.faces.len())
            .map_err(|_| "The mesh already holds the maximum number of faces a u32 can index")?;
        self.faces.push(face);
        Ok(index)
    }
}

// ===============================================================================================
// Modifying elements
// ===============================================================================================

impl MeshData3 {
    /// Move a point to a new position.
    ///
    /// Nothing else is touched. In particular, a stored normal at this point is not recomputed, and
    /// will be stale if the surrounding geometry changed meaningfully.
    ///
    /// # Arguments
    ///
    /// * `index`: the point to move
    /// * `point`: its new position
    ///
    /// returns: `Result<()>`, failing if the index is out of range
    pub fn set_point(&mut self, index: u32, point: Point3) -> Result<()> {
        let i = self.point_index(index)?;
        self.points[i] = point;
        Ok(())
    }

    /// Get mutable access to the whole point buffer.
    ///
    /// This is safe to hand out because no invariant of the type depends on where a point is: the
    /// face indices stay valid whatever the positions become, and the length cannot change through
    /// a slice. There is deliberately no equivalent for faces, whose indices have to be checked.
    ///
    /// As with `set_point`, stored normals are not recomputed and may end up stale.
    pub fn points_mut(&mut self) -> &mut [Point3] {
        &mut self.points
    }

    /// Replace a face with a different triple of point indices.
    ///
    /// # Arguments
    ///
    /// * `index`: the face to replace
    /// * `face`: three indices into the point buffer, each of which must refer to a point which
    ///   exists
    ///
    /// returns: `Result<()>`, failing if the face index or any point index is out of range
    pub fn set_face(&mut self, index: u32, face: [u32; 3]) -> Result<()> {
        let i = self.face_index(index)?;
        self.check_face(&face)?;
        self.faces[i] = face;
        Ok(())
    }
}

// ===============================================================================================
// Removing elements
// ===============================================================================================

impl MeshData3 {
    /// Remove a point, along with every face which references it.
    ///
    /// Attribute values are carried through the deletion for both domains, so a mesh with per-point
    /// or per-face attributes stays coherent.
    ///
    /// Because removal is ordered rather than a swap, **every point index above the one removed
    /// shifts down by one**, and any index a caller was holding above it is invalidated. Face
    /// indices shift too, wherever a face was dropped.
    ///
    /// # Arguments
    ///
    /// * `index`: the point to remove
    ///
    /// returns: `Result<()>`, failing if the index is out of range
    pub fn remove_point(&mut self, index: u32) -> Result<()> {
        let i = self.point_index(index)?;

        let mut point_mask = IndexMask::new(self.points.len(), true);
        point_mask.set(i, false);

        let mut face_mask = IndexMask::new(self.faces.len(), true);
        for (f, face) in self.faces.iter().enumerate() {
            if face.contains(&index) {
                face_mask.set(f, false);
            }
        }

        // Subset both attribute domains before modifying the mesh so that a failure in either
        // domain leaves it untouched.
        let point_attrs = self.point_attrs.subset(&point_mask)?;
        let face_attrs = self.face_attrs.subset(&face_mask)?;
        self.point_attrs = point_attrs;
        self.face_attrs = face_attrs;

        self.points.remove(i);
        self.faces.retain(|face| !face.contains(&index));

        // Everything above the removed point has shifted down by one.
        for face in self.faces.iter_mut() {
            for v in face.iter_mut() {
                if *v > index {
                    *v -= 1;
                }
            }
        }

        Ok(())
    }

    /// Remove a face, leaving the points alone.
    ///
    /// Points which the removed face was the last to reference become orphans, which is a legal
    /// state for this type. Nothing renumbers them.
    ///
    /// Face indices above the one removed shift down by one.
    ///
    /// # Arguments
    ///
    /// * `index`: the face to remove
    ///
    /// returns: `Result<()>`, failing if the index is out of range
    pub fn remove_face(&mut self, index: u32) -> Result<()> {
        let i = self.face_index(index)?;

        let mut face_mask = IndexMask::new(self.faces.len(), true);
        face_mask.set(i, false);

        self.face_attrs = self.face_attrs.subset(&face_mask)?;
        self.faces.remove(i);

        Ok(())
    }
}

// ===============================================================================================
// Bounds checking
// ===============================================================================================

impl MeshData3 {
    /// Convert a point index into a slice index, failing if it is out of range.
    fn point_index(&self, index: u32) -> Result<usize> {
        let i = index as usize;
        if i >= self.points.len() {
            return Err(format!(
                "Point index {index} is out of range for a mesh with {} points",
                self.points.len()
            )
            .into());
        }
        Ok(i)
    }

    /// Convert a face index into a slice index, failing if it is out of range.
    fn face_index(&self, index: u32) -> Result<usize> {
        let i = index as usize;
        if i >= self.faces.len() {
            return Err(format!(
                "Face index {index} is out of range for a mesh with {} faces",
                self.faces.len()
            )
            .into());
        }
        Ok(i)
    }

    /// Verify that every index of a face refers to a point which exists.
    fn check_face(&self, face: &[u32; 3]) -> Result<()> {
        for index in face {
            if *index as usize >= self.points.len() {
                return Err(format!(
                    "Face refers to point {index}, but the mesh has only {} points",
                    self.points.len()
                )
                .into());
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::mesh::data::Attr3;

    /// Two triangles sharing an edge, over four points.
    fn square_mesh() -> MeshData3 {
        MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
        .unwrap()
    }

    #[test]
    fn pushed_points_return_usable_indices() -> Result<()> {
        let mut mesh = MeshData3::empty();

        let a = mesh.push_point(Point3::new(0.0, 0.0, 0.0))?;
        let b = mesh.push_point(Point3::new(1.0, 0.0, 0.0))?;
        let c = mesh.push_point(Point3::new(0.0, 1.0, 0.0))?;

        assert_eq!([a, b, c], [0, 1, 2]);

        let f = mesh.push_face([a, b, c])?;
        assert_eq!(f, 0);
        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 1);

        Ok(())
    }

    #[test]
    fn push_face_rejects_an_out_of_range_index() {
        let mut mesh = square_mesh();
        assert!(mesh.push_face([0, 1, 4]).is_err());
        assert_eq!(mesh.face_count(), 2);
    }

    #[test]
    fn insertion_is_blocked_while_attributes_are_present() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.set_point_colors(Some(vec![[0, 0, 0]; 4]))?;

        let err = mesh.push_point(Point3::origin()).unwrap_err().to_string();
        assert!(err.contains("point_colors"), "{err}");
        assert_eq!(mesh.point_count(), 4);

        // A per-point attribute does not block adding a face.
        mesh.push_face([0, 1, 3])?;

        mesh.set_face_labels(Some(vec![1, 2, 3]))?;
        let err = mesh.push_face([0, 1, 2]).unwrap_err().to_string();
        assert!(err.contains("face_labels"), "{err}");

        Ok(())
    }

    #[test]
    fn insertion_names_an_open_map_attribute_too() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.insert_point_attr("confidence", Attr3::Scalar(vec![0.5; 4]))?;

        let err = mesh.push_point(Point3::origin()).unwrap_err().to_string();
        assert!(err.contains("confidence"), "{err}");

        Ok(())
    }

    #[test]
    fn points_can_be_moved() -> Result<()> {
        let mut mesh = square_mesh();

        mesh.set_point(2, Point3::new(9.0, 9.0, 9.0))?;
        assert_eq!(mesh.points()[2], Point3::new(9.0, 9.0, 9.0));

        mesh.points_mut()[0].x = -1.0;
        assert_eq!(mesh.points()[0].x, -1.0);

        assert!(mesh.set_point(4, Point3::origin()).is_err());

        Ok(())
    }

    #[test]
    fn faces_can_be_replaced() -> Result<()> {
        let mut mesh = square_mesh();

        mesh.set_face(0, [1, 2, 3])?;
        assert_eq!(mesh.faces()[0], [1, 2, 3]);

        assert!(mesh.set_face(0, [1, 2, 9]).is_err());
        assert!(mesh.set_face(7, [0, 1, 2]).is_err());

        // The rejected calls must have left the face as it was.
        assert_eq!(mesh.faces()[0], [1, 2, 3]);

        Ok(())
    }

    #[test]
    fn removing_a_point_drops_its_faces_and_renumbers_the_rest() -> Result<()> {
        let mut mesh = square_mesh();

        // Point 3 is used only by the second face.
        mesh.remove_point(3)?;

        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.faces(), &[[0, 1, 2]]);

        Ok(())
    }

    #[test]
    fn removal_shifts_the_indices_above_it() -> Result<()> {
        let mut mesh = square_mesh();

        // Point 0 is used by both faces, so both go, and points 1, 2, 3 become 0, 1, 2.
        mesh.remove_point(0)?;

        assert_eq!(mesh.point_count(), 3);
        assert_eq!(mesh.face_count(), 0);
        assert_eq!(mesh.points()[0], Point3::new(1.0, 0.0, 0.0));

        // Now build a face over the renumbered points and remove the middle one.
        let mut mesh = square_mesh();
        mesh.remove_point(1)?;

        assert_eq!(mesh.point_count(), 3);
        // Only the face which avoided point 1 survives, with 2 and 3 shifted down to 1 and 2.
        assert_eq!(mesh.faces(), &[[0, 1, 2]]);

        Ok(())
    }

    #[test]
    fn removing_a_point_carries_attributes_through_both_domains() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.set_point_stdev(Some(vec![0.1, 0.2, 0.3, 0.4]))?;
        mesh.set_face_labels(Some(vec![10, 20]))?;
        mesh.insert_point_attr("confidence", Attr3::Scalar(vec![1.0, 2.0, 3.0, 4.0]))?;

        // Point 3 belongs only to face 1.
        mesh.remove_point(3)?;

        assert_eq!(mesh.point_stdev().unwrap(), &[0.1, 0.2, 0.3]);
        assert_eq!(
            mesh.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[1.0, 2.0, 3.0]
        );
        assert_eq!(mesh.face_labels().unwrap(), &[10]);

        // And the result is still internally coherent.
        mesh.point_attrs().validate(mesh.point_count())?;
        mesh.face_attrs().validate(mesh.face_count())?;

        Ok(())
    }

    #[test]
    fn removing_a_face_leaves_orphan_points() -> Result<()> {
        let mut mesh = square_mesh();
        mesh.set_face_labels(Some(vec![10, 20]))?;

        mesh.remove_face(1)?;

        assert_eq!(mesh.face_count(), 1);
        assert_eq!(mesh.point_count(), 4);
        assert_eq!(mesh.face_labels().unwrap(), &[10]);

        Ok(())
    }

    #[test]
    fn removal_rejects_an_out_of_range_index() {
        let mut mesh = square_mesh();

        assert!(mesh.remove_point(4).is_err());
        assert!(mesh.remove_face(2).is_err());
        assert_eq!(mesh.point_count(), 4);
        assert_eq!(mesh.face_count(), 2);
    }
}
