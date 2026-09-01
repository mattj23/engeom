//! This module contains `MeshView3`, a borrowed, read-only view of a triangle mesh. It provides the
//! point and face buffers and the attribute sets for both domains, without taking ownership or
//! providing an acceleration structure.
//!
//! It allows code that only needs to *read* a mesh—especially file writers—to accept either mesh
//! container through one interface. Both `MeshData3` and `Mesh3` convert into it without copying.
//! This lets `Mesh3::save_*` write directly from an accelerated mesh instead of first cloning every
//! buffer into a `MeshData3`. Because `parry3d` provides no way to move the buffers out of a
//! `TriMesh`, that copy would otherwise be unavoidable.
//!
//! A function that takes `impl Into<MeshView3<'a>>` can be called with `&mesh_data`, `&mesh`, or
//! an existing view.

use super::{Mesh3, MeshData3};
use crate::geom3::attributes3::{FaceAttrSet3, PointAttrSet3};
use crate::geom3::mesh::data::check_face_indices;
use crate::{Point3, Result};

/// A borrowed, read-only view of a triangle mesh and its attributes.
///
/// A view can be built from `&MeshData3` or `&Mesh3` through `From`, preserving the source's
/// invariants: every face index is in range, and every attribute array matches its domain's count.
/// It can also be built from raw parts through `try_new`, which validates those invariants.
#[derive(Clone, Copy)]
pub struct MeshView3<'a> {
    points: &'a [Point3],
    faces: &'a [[u32; 3]],
    point_attrs: &'a PointAttrSet3,
    face_attrs: &'a FaceAttrSet3,
}

impl<'a> MeshView3<'a> {
    /// Assemble a view from raw parts, verifying that they describe a coherent mesh.
    ///
    /// # Arguments
    ///
    /// * `points`: the point positions
    /// * `faces`: triangles given as triples of indices into `points`
    /// * `point_attrs`: the per-point attributes, whose arrays must match the point count
    /// * `face_attrs`: the per-face attributes, whose arrays must match the face count
    ///
    /// returns: `Result<MeshView3>`, failing if any face refers to a point which does not exist or
    /// if any attribute array is the wrong length
    pub fn try_new(
        points: &'a [Point3],
        faces: &'a [[u32; 3]],
        point_attrs: &'a PointAttrSet3,
        face_attrs: &'a FaceAttrSet3,
    ) -> Result<Self> {
        check_face_indices(faces, points.len())?;
        point_attrs.validate(points.len())?;
        face_attrs.validate(faces.len())?;

        Ok(Self {
            points,
            faces,
            point_attrs,
            face_attrs,
        })
    }

    /// The point positions.
    pub fn points(&self) -> &'a [Point3] {
        self.points
    }

    /// The faces, as triples of indices into `points`.
    pub fn faces(&self) -> &'a [[u32; 3]] {
        self.faces
    }

    /// The per-point attributes.
    pub fn point_attrs(&self) -> &'a PointAttrSet3 {
        self.point_attrs
    }

    /// The per-face attributes.
    pub fn face_attrs(&self) -> &'a FaceAttrSet3 {
        self.face_attrs
    }

    /// The number of points.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// The number of faces.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Returns true if the mesh carries any attributes in either domain.
    pub fn has_attrs(&self) -> bool {
        !self.point_attrs.is_empty() || !self.face_attrs.is_empty()
    }

    /// List the names of all present attributes, including typed fields, with point-domain
    /// attributes followed by face-domain attributes.
    ///
    /// Operations that cannot preserve the attributes report this list so callers can see what
    /// would be lost.
    pub fn attr_labels(&self) -> Vec<&'a str> {
        let mut labels = self.point_attrs.attr_labels();
        labels.extend(self.face_attrs.attr_labels());
        labels
    }

    /// Verify that the caller has accepted the loss of this mesh's attributes when using a format
    /// that cannot represent them.
    ///
    /// Writers for geometry-only formats call this before doing any work. If the mesh has no
    /// attributes, there is nothing to lose and this always succeeds, so the flag matters only
    /// when data would actually be discarded.
    ///
    /// This check prevents silent data loss that may otherwise be discovered much later. For
    /// example, a user could save a mesh carrying measured uncertainty, close the session, and
    /// discover weeks later that the data is gone. Raising an error at the moment of loss makes
    /// the risk explicit instead of merely recording it in a value or log line that could be
    /// overlooked.
    ///
    /// # Arguments
    ///
    /// * `format`: the name of the target format, used in the error message
    /// * `allow_loss`: whether the caller has accepted the loss, which comes from the
    ///   `allow_attribute_loss` field of the format's options struct
    ///
    /// returns: `Result<()>`
    pub fn check_attribute_loss(&self, format: &str, allow_loss: bool) -> Result<()> {
        if allow_loss || !self.has_attrs() {
            return Ok(());
        }

        Err(format!(
            "Writing to {} would discard the attributes on this mesh ({}), because the format \
             cannot represent them. Set `allow_attribute_loss` to accept this.",
            format,
            self.attr_labels().join(", ")
        )
        .into())
    }
}

impl<'a> From<&'a MeshData3> for MeshView3<'a> {
    fn from(mesh: &'a MeshData3) -> Self {
        Self {
            points: mesh.points(),
            faces: mesh.faces(),
            point_attrs: mesh.point_attrs(),
            face_attrs: mesh.face_attrs(),
        }
    }
}

impl<'a> From<&'a Mesh3> for MeshView3<'a> {
    fn from(mesh: &'a Mesh3) -> Self {
        Self {
            points: mesh.points(),
            faces: mesh.faces(),
            point_attrs: mesh.point_attrs(),
            face_attrs: mesh.face_attrs(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom3::Attr3;

    fn square() -> MeshData3 {
        let mut mesh = MeshData3::new(
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )
        .unwrap();
        mesh.set_point_stdev(Some(vec![0.1; 4])).unwrap();
        mesh.set_face_labels(Some(vec![1, 2])).unwrap();
        mesh.insert_face_attr("material_index", Attr3::Label(vec![5, 6]))
            .unwrap();
        mesh
    }

    #[test]
    fn both_containers_produce_the_same_view() -> Result<()> {
        let data = square();
        let mesh = Mesh3::from_data(data.clone(), false)?;

        let from_data = MeshView3::from(&data);
        let from_mesh = MeshView3::from(&mesh);

        assert_eq!(from_data.points(), from_mesh.points());
        assert_eq!(from_data.faces(), from_mesh.faces());
        assert_eq!(from_data.point_attrs(), from_mesh.point_attrs());
        assert_eq!(from_data.face_attrs(), from_mesh.face_attrs());
        assert_eq!(from_data.point_count(), 4);
        assert_eq!(from_data.face_count(), 2);

        Ok(())
    }

    #[test]
    fn attr_labels_list_both_domains_point_first() {
        let data = square();
        let view = MeshView3::from(&data);

        assert!(view.has_attrs());
        assert_eq!(
            view.attr_labels(),
            vec!["point_stdev", "face_labels", "material_index"]
        );

        let bare = MeshData3::create_box(1.0, 1.0, 1.0);
        let view = MeshView3::from(&bare);
        assert!(!view.has_attrs());
        assert!(view.attr_labels().is_empty());
    }

    #[test]
    fn attribute_loss_is_refused_unless_accepted() {
        let data = square();
        let view = MeshView3::from(&data);

        let message = match view.check_attribute_loss("a test format", false) {
            Err(e) => e.to_string(),
            Ok(()) => panic!("an attributed mesh must be refused"),
        };
        assert!(message.contains("a test format"), "{message}");
        assert!(message.contains("point_stdev"), "{message}");
        assert!(message.contains("material_index"), "{message}");

        assert!(view.check_attribute_loss("a test format", true).is_ok());

        let bare = MeshData3::create_box(1.0, 1.0, 1.0);
        assert!(
            MeshView3::from(&bare)
                .check_attribute_loss("a test format", false)
                .is_ok()
        );
    }

    #[test]
    fn try_new_checks_the_parts() {
        let points = [Point3::origin(), Point3::new(1.0, 0.0, 0.0)];
        let faces = [[0, 1, 2]];
        let point_attrs = PointAttrSet3::empty();
        let face_attrs = FaceAttrSet3::empty();

        assert!(MeshView3::try_new(&points, &faces, &point_attrs, &face_attrs).is_err());

        let points = [
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        assert!(MeshView3::try_new(&points, &faces, &point_attrs, &face_attrs).is_ok());

        let mut bad = FaceAttrSet3::empty();
        bad.set_labels(Some(vec![1, 2]), 2).unwrap();
        assert!(MeshView3::try_new(&points, &faces, &point_attrs, &bad).is_err());
    }
}
