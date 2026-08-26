//! This module crops a `MeshData3` to a bounding volume.
//!
//! It exists on the raw container rather than on `Mesh3` because the motivating use is reducing scan
//! data *before* it is committed to an accelerated representation. Cropping through `Mesh3` would
//! mean paying to build a BVH over the whole mesh in order to throw most of it away.
//!
//! The two halves are kept separate, the same way the rest of `MeshData3` keeps them: the
//! `compute_*_mask_in_aabb` methods answer which elements are in the volume, and the existing
//! `extract_subset_faces`/`extract_subset_points` turn a mask into a mesh. `extract_subset_faces_in_aabb`
//! is the two composed, for the common case where the caller only wants the cropped mesh.
//!
//! # Which domain to crop by
//!
//! Both directions are available, and they are not equivalent:
//!
//! - Going through the **face** mask and `extract_subset_faces` keeps the selected triangles and
//!   only the points they use, so a point in the volume which belongs to no surviving face is
//!   dropped.
//! - Going through the **point** mask and `extract_subset_points` keeps every point in the volume,
//!   including any which ends up orphaned, and keeps only the faces whose three points all survive.
//!
//! For cropping ahead of an accelerated structure the face route is the one you want, since an
//! orphan point is dead weight to a mesh. The point route is there for callers who are selecting
//! points and mean it.

use super::MeshData3;
use crate::Result;
use crate::common::IndexMask;
use crate::geom3::Aabb3;
use crate::geom3::mesh::algorithms::volumes::{
    compute_face_mask_in_aabb, compute_point_mask_in_aabb,
};

// ===============================================================================================
// Selecting by volume
// ===============================================================================================

impl MeshData3 {
    /// Compute the mask of points which lie inside an axis-aligned bounding box.
    ///
    /// A point exactly on a face of the box counts as inside.
    ///
    /// # Arguments
    ///
    /// * `aabb`: the box to test against, in the same frame as the mesh
    ///
    /// returns: `IndexMask` of length equal to the point count
    pub fn compute_point_mask_in_aabb(&self, aabb: &Aabb3) -> IndexMask {
        compute_point_mask_in_aabb(&self.points, aabb)
    }

    /// Compute the mask of faces which lie inside an axis-aligned bounding box.
    ///
    /// The `all_vertices` flag chooses which of the two meanings of "inside" applies, and they
    /// differ for a bounded volume in a way they do not for a half-space: a triangle can span the
    /// box without a single vertex landing in it.
    ///
    /// # Arguments
    ///
    /// * `aabb`: the box to test against, in the same frame as the mesh
    /// * `all_vertices`: if `true`, a face is selected only when it lies entirely inside the box;
    ///   if `false`, a face is selected when any part of it touches the box at all
    ///
    /// returns: `Result<IndexMask>` of length equal to the face count
    pub fn compute_face_mask_in_aabb(&self, aabb: &Aabb3, all_vertices: bool) -> Result<IndexMask> {
        compute_face_mask_in_aabb(&self.points, &self.faces, aabb, all_vertices)
    }
}

// ===============================================================================================
// Cropping
// ===============================================================================================

impl MeshData3 {
    /// Crop the mesh to an axis-aligned bounding box, keeping the faces it selects.
    ///
    /// Only points used by a surviving face are kept, so a point inside the box which the crop
    /// orphans is dropped. Attributes are carried through in both domains. A crop which selects
    /// nothing yields an empty mesh rather than an error, since a box which misses the mesh is a
    /// normal thing for a caller to have asked about.
    ///
    /// Equivalent to [`MeshData3::extract_subset_faces`] with the mask from
    /// [`MeshData3::compute_face_mask_in_aabb`].
    ///
    /// # Arguments
    ///
    /// * `aabb`: the box to crop to, in the same frame as the mesh
    /// * `all_vertices`: if `true`, keep only the faces lying entirely inside the box; if `false`,
    ///   keep every face which touches it at all
    ///
    /// returns: `Result<MeshData3>`
    pub fn extract_subset_faces_in_aabb(&self, aabb: &Aabb3, all_vertices: bool) -> Result<Self> {
        let mask = self.compute_face_mask_in_aabb(aabb, all_vertices)?;
        self.extract_subset_faces(&mask)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point3;
    use crate::geom3::mesh::data::Attr3;

    /// Two triangles sharing the edge 0-2, plus a fifth point which belongs to no face and sits
    /// well outside any box these tests use. Both attribute domains are populated so the crop can
    /// be checked for carrying them.
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

    /// A box around the half of the quad with x below 0.5, which contains points 0 and 3 outright
    /// and cuts through both faces.
    fn left_half() -> Aabb3 {
        Aabb3::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(0.5, 2.0, 1.0))
    }

    #[test]
    fn the_point_mask_selects_the_points_in_the_box() {
        let mesh = fixture();
        let mask = mesh.compute_point_mask_in_aabb(&left_half());

        assert_eq!(mask.to_indices(), vec![0, 3]);
    }

    #[test]
    fn the_face_mask_splits_on_all_vertices() -> Result<()> {
        let mesh = fixture();

        // Neither face lies wholly in the left half, but both cross into it.
        assert_eq!(
            mesh.compute_face_mask_in_aabb(&left_half(), true)?
                .count_true(),
            0
        );
        assert_eq!(
            mesh.compute_face_mask_in_aabb(&left_half(), false)?
                .to_indices(),
            vec![0, 1]
        );

        Ok(())
    }

    #[test]
    fn a_box_around_everything_keeps_every_face() -> Result<()> {
        let mesh = fixture();
        let all = Aabb3::new(
            Point3::new(-10.0, -10.0, -10.0),
            Point3::new(10.0, 10.0, 10.0),
        );

        let sub = mesh.extract_subset_faces_in_aabb(&all, true)?;

        assert_eq!(sub.face_count(), 2);

        // The orphan point is not used by a surviving face, so the crop drops it even though the
        // box contains it. Cropping by points instead would have kept it.
        assert_eq!(sub.point_count(), 4);

        Ok(())
    }

    #[test]
    fn a_crop_carries_the_attributes_of_the_faces_it_keeps() -> Result<()> {
        let mesh = fixture();

        // A box holding only point 1, which face 0 uses and face 1 does not.
        let corner = Aabb3::new(Point3::new(0.75, -0.5, -0.5), Point3::new(1.5, 0.5, 0.5));
        let sub = mesh.extract_subset_faces_in_aabb(&corner, false)?;

        assert_eq!(sub.face_count(), 1);
        assert_eq!(sub.point_count(), 3);
        assert_eq!(sub.faces(), &[[0, 1, 2]]);

        assert_eq!(sub.face_labels().unwrap(), &[10]);
        assert_eq!(sub.point_stdev().unwrap(), &[0.0, 0.1, 0.2]);
        assert_eq!(
            sub.face_attr("quality").unwrap().as_scalar().unwrap(),
            &[0.5]
        );

        Ok(())
    }

    /// A box which misses the mesh is a normal answer, not an error. `Mesh3` cannot represent this,
    /// which is part of why the crop lives here.
    #[test]
    fn a_box_which_selects_nothing_yields_an_empty_mesh() -> Result<()> {
        let mesh = fixture();
        let far = Aabb3::new(
            Point3::new(100.0, 100.0, 100.0),
            Point3::new(101.0, 101.0, 101.0),
        );

        let sub = mesh.extract_subset_faces_in_aabb(&far, false)?;

        assert!(sub.is_empty());
        assert_eq!(sub.face_count(), 0);
        assert_eq!(sub.point_count(), 0);

        Ok(())
    }

    /// The straddle case at the container level: a box entirely inside one triangle, touching none
    /// of its vertices, still selects it.
    #[test]
    fn a_box_inside_a_single_face_selects_it() -> Result<()> {
        let mesh = MeshData3::new(
            vec![
                Point3::new(-10.0, -10.0, 0.0),
                Point3::new(10.0, -10.0, 0.0),
                Point3::new(0.0, 10.0, 0.0),
            ],
            vec![[0, 1, 2]],
        )?;

        let inner = Aabb3::new(Point3::new(-1.0, -1.0, -1.0), Point3::new(1.0, 1.0, 1.0));

        assert_eq!(mesh.compute_point_mask_in_aabb(&inner).count_true(), 0);
        assert_eq!(
            mesh.compute_face_mask_in_aabb(&inner, false)?.to_indices(),
            vec![0]
        );

        Ok(())
    }
}
