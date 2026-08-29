//! This module contains `FaceAttrSet3`, the bundle of per-face attributes carried by a container
//! which has a face domain, which today means `MeshData3` and `Mesh3`.
//!
//! It is the face-domain peer of `PointAttrSet3` and follows the same rules: a small set of
//! first-class attributes are stored as typed fields, and everything else lives in an open,
//! name-keyed map of `Attr3` values. An attribute earns a typed field if an algorithm reads it in a
//! hot loop, or if the container has to do something to it that it doesn't do to an arbitrary
//! scalar.
//!
//! # What is not here
//!
//! There is no typed field for face normals. Nothing in the library stores them: a face normal is
//! fully determined by the three vertices and their winding, so it is always computed on demand
//! (see `compute_face_normals`) rather than carried as data which could disagree with the geometry.
//!
//! There is also no `scale_in_place` or `flip_in_place`. The face domain holds nothing with length
//! units and nothing whose sign depends on orientation, so a uniform scale or a winding reversal
//! leaves every face attribute alone.
//!
//! # Error message naming
//!
//! The typed fields report themselves as `face_colors` and `face_labels` in error messages rather
//! than by their bare field names, so that a message from a mesh with both domains says which one
//! it is talking about.

use super::{
    Attr3, check_both_or_neither, check_keys_match, check_len, check_reserved, check_same_variant,
    clone_indexed, clone_masked, extend_option,
};
use crate::common::IndexMask;
use crate::{Iso3, Result};
use std::collections::HashMap;

/// The set of per-face attributes attached to a container, holding both the first-class typed
/// attributes and an open, name-keyed map of everything else.
///
/// This type does not know the face count of the container it belongs to and so it can't enforce
/// that its arrays are the right length. Validation is the job of the owner, which does know the
/// count; there is a `validate(n)` method for that, and the setters take the count as an argument.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FaceAttrSet3 {
    colors: Option<Vec<[u8; 3]>>,
    labels: Option<Vec<u32>>,
    open: HashMap<String, Attr3>,
}

// ===============================================================================================
// Construction and access
// ===============================================================================================

impl FaceAttrSet3 {
    /// Create an attribute set with no attributes of any kind.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Returns true if no attributes of any kind are present.
    pub fn is_empty(&self) -> bool {
        self.colors.is_none() && self.labels.is_none() && self.open.is_empty()
    }

    /// List the names of every per-face attribute which is present, typed fields included.
    ///
    /// This is what an operation which cannot supply values for a new face reports back, so the
    /// caller can see what is standing in the way.
    pub fn attr_labels(&self) -> Vec<&str> {
        let mut names = Vec::new();
        if self.colors.is_some() {
            names.push("face_colors");
        }
        if self.labels.is_some() {
            names.push("face_labels");
        }
        names.extend(self.open.keys().map(|k| k.as_str()));
        names
    }

    /// Get the per-face RGB colors, if present.
    pub fn colors(&self) -> Option<&[[u8; 3]]> {
        self.colors.as_deref()
    }

    /// Get the per-face labels, if present. These identify which region, patch, scan pass, or
    /// material each face belongs to.
    pub fn labels(&self) -> Option<&[u32]> {
        self.labels.as_deref()
    }

    /// Get the open-map attribute stored under the given name, if present.
    pub fn attr(&self, name: &str) -> Option<&Attr3> {
        self.open.get(name)
    }

    /// Iterate over the names of the open-map attributes, in no particular order.
    pub fn attr_names(&self) -> impl Iterator<Item = &str> {
        self.open.keys().map(|k| k.as_str())
    }
}

// ===============================================================================================
// Mutation
// ===============================================================================================

impl FaceAttrSet3 {
    /// Set or clear the per-face RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them
    /// * `n_faces`: the face count of the owning container, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_colors(&mut self, values: Option<Vec<[u8; 3]>>, n_faces: usize) -> Result<()> {
        check_len(values.as_deref().map(|v| v.len()), n_faces, "face_colors")?;
        self.colors = values;
        Ok(())
    }

    /// Set or clear the per-face labels.
    ///
    /// # Arguments
    ///
    /// * `values`: the labels to store, or `None` to clear them
    /// * `n_faces`: the face count of the owning container, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_labels(&mut self, values: Option<Vec<u32>>, n_faces: usize) -> Result<()> {
        check_len(values.as_deref().map(|v| v.len()), n_faces, "face_labels")?;
        self.labels = values;
        Ok(())
    }

    /// Insert an open-map attribute under the given name, replacing any attribute already stored
    /// there.
    ///
    /// # Arguments
    ///
    /// * `name`: the key to store under, which must not be one of the reserved names
    /// * `attr`: the attribute array to store
    /// * `n_faces`: the face count of the owning container, which `attr` must match
    ///
    /// returns: `Result<()>`
    pub fn insert_attr(&mut self, name: &str, attr: Attr3, n_faces: usize) -> Result<()> {
        check_reserved(name)?;
        check_len(Some(attr.len()), n_faces, name)?;
        self.open.insert(name.to_string(), attr);
        Ok(())
    }

    /// Remove and return the open-map attribute stored under the given name.
    pub fn remove_attr(&mut self, name: &str) -> Option<Attr3> {
        self.open.remove(name)
    }
}

// ===============================================================================================
// Validation, subsetting, and combination
// ===============================================================================================

impl FaceAttrSet3 {
    /// Verify that every attribute array present has a length matching the face count.
    ///
    /// # Arguments
    ///
    /// * `n_faces`: the face count of the owning container
    ///
    /// returns: `Result<()>`
    pub fn validate(&self, n_faces: usize) -> Result<()> {
        check_len(
            self.colors.as_ref().map(|v| v.len()),
            n_faces,
            "face_colors",
        )?;
        check_len(
            self.labels.as_ref().map(|v| v.len()),
            n_faces,
            "face_labels",
        )?;

        for (name, attr) in self.open.iter() {
            check_len(Some(attr.len()), n_faces, name)?;
        }

        Ok(())
    }

    /// Create a new attribute set containing only the faces selected by the given mask, preserving
    /// their original order.
    ///
    /// # Arguments
    ///
    /// * `mask`: selects which faces survive, and must match the owning container's face count
    ///
    /// returns: `Result<FaceAttrSet3>`
    pub fn subset(&self, mask: &IndexMask) -> Result<Self> {
        let mut open = HashMap::with_capacity(self.open.len());
        for (name, attr) in self.open.iter() {
            open.insert(name.clone(), attr.clone_indices_of(mask)?);
        }

        Ok(Self {
            colors: clone_masked(self.colors.as_deref(), mask)?,
            labels: clone_masked(self.labels.as_deref(), mask)?,
            open,
        })
    }

    /// Create a new attribute set containing the faces at the given indices, in the order the
    /// indices are given. Indices may repeat.
    ///
    /// # Arguments
    ///
    /// * `indices`: the faces to take, each of which must be less than the owning container's
    ///   face count
    ///
    /// returns: `Result<FaceAttrSet3>`
    pub fn subset_indices(&self, indices: &[usize]) -> Result<Self> {
        let mut open = HashMap::with_capacity(self.open.len());
        for (name, attr) in self.open.iter() {
            open.insert(name.clone(), attr.clone_indices(indices)?);
        }

        Ok(Self {
            colors: clone_indexed(self.colors.as_deref(), indices)?,
            labels: clone_indexed(self.labels.as_deref(), indices)?,
            open,
        })
    }

    /// Append the contents of another attribute set onto the end of this one.
    ///
    /// Attributes are all-or-nothing: a typed field or an open-map key which is present on one side
    /// and absent on the other is an error, because there is no correct value to pad the missing
    /// side with. The check runs over the whole set before anything is modified, so a failure
    /// leaves this attribute set untouched.
    ///
    /// # Arguments
    ///
    /// * `other`: the attribute set to append
    ///
    /// returns: `Result<()>`
    pub fn extend_from(&mut self, other: &Self) -> Result<()> {
        self.check_extend_from(other)?;
        self.apply_extend_from(other)
    }

    /// Run every check an append has to pass, without modifying anything.
    ///
    /// This is separate from the append itself so that a mesh can check both of its domains before
    /// modifying either. Without it, a face-domain mismatch discovered after the point domain had
    /// already been extended would leave the mesh half appended.
    pub(crate) fn check_extend_from(&self, other: &Self) -> Result<()> {
        check_both_or_neither(self.colors.is_some(), other.colors.is_some(), "face_colors")?;
        check_both_or_neither(self.labels.is_some(), other.labels.is_some(), "face_labels")?;

        check_keys_match(&self.open, &other.open, "face")?;

        // Variant mismatches between two attributes of the same name would fail partway through the
        // append, so they are checked up front as well.
        for (name, attr) in self.open.iter() {
            check_same_variant(attr, &other.open[name], name)?;
        }

        Ok(())
    }

    /// Perform the append. Only valid after `check_extend_from` has succeeded for the same pair.
    pub(crate) fn apply_extend_from(&mut self, other: &Self) -> Result<()> {
        extend_option(&mut self.colors, other.colors.as_deref());
        extend_option(&mut self.labels, other.labels.as_deref());

        for (name, attr) in self.open.iter_mut() {
            attr.extend_from(&other.open[name])?;
        }

        Ok(())
    }
}

// ===============================================================================================
// Geometric operations
// ===============================================================================================

impl FaceAttrSet3 {
    /// Transform the spatial attributes in place by the given isometry.
    ///
    /// This rotates every `Attr3::Vector` attribute. Because these hold directions rather than
    /// positions, only the rotation component of the isometry has any effect. The typed fields are
    /// untouched: neither a color nor a label has a direction.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    pub fn transform_in_place(&mut self, iso: &Iso3) {
        for attr in self.open.values_mut() {
            attr.transform_in_place(iso);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use crate::geom3::attributes3::RESERVED_ATTR_NAMES;
    use approx::assert_relative_eq;
    use std::f64::consts::FRAC_PI_2;

    const N: usize = 4;

    /// An attribute set carrying every typed field and one open attribute.
    fn full_attrs() -> FaceAttrSet3 {
        let mut attrs = FaceAttrSet3::empty();

        attrs
            .set_colors(Some(vec![[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]]), N)
            .unwrap();

        attrs.set_labels(Some(vec![10, 11, 12, 13]), N).unwrap();

        attrs
            .insert_attr("material_index", Attr3::Label(vec![1, 2, 3, 4]), N)
            .unwrap();

        attrs
    }

    #[test]
    fn empty_set_reports_empty() {
        let attrs = FaceAttrSet3::empty();
        assert!(attrs.is_empty());
        assert!(attrs.colors().is_none());
        assert!(attrs.labels().is_none());
        assert_eq!(attrs.attr_names().count(), 0);

        assert!(!full_attrs().is_empty());
    }

    #[test]
    fn attr_labels_names_every_present_attribute() {
        let attrs = full_attrs();
        let labels = attrs.attr_labels();

        for expected in ["face_colors", "face_labels", "material_index"] {
            assert!(
                labels.contains(&expected),
                "missing '{expected}' in {labels:?}"
            );
        }

        assert!(FaceAttrSet3::empty().attr_labels().is_empty());
    }

    #[test]
    fn setters_reject_a_length_mismatch() {
        let mut attrs = FaceAttrSet3::empty();

        assert!(attrs.set_colors(Some(vec![[0, 0, 0]]), N).is_err());
        assert!(attrs.set_labels(Some(Vec::new()), N).is_err());
        assert!(attrs.insert_attr("q", Attr3::Scalar(vec![1.0]), N).is_err());

        // A rejected set must not have stored anything.
        assert!(attrs.is_empty());
    }

    #[test]
    fn setters_accept_none_regardless_of_count() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.set_colors(None, N)?;
        assert!(attrs.colors().is_none());
        attrs.set_labels(None, 0)?;
        assert!(attrs.labels().is_none());
        Ok(())
    }

    #[test]
    fn open_map_rejects_reserved_names() {
        let mut attrs = FaceAttrSet3::empty();

        for name in RESERVED_ATTR_NAMES {
            assert!(
                attrs
                    .insert_attr(name, Attr3::Scalar(vec![0.0; N]), N)
                    .is_err(),
                "expected '{name}' to be rejected as an attribute name"
            );
        }

        // A name which merely contains a reserved word is fine.
        assert!(
            attrs
                .insert_attr("label_confidence", Attr3::Scalar(vec![0.0; N]), N)
                .is_ok()
        );
    }

    #[test]
    fn remove_attr_returns_what_was_stored() {
        let mut attrs = full_attrs();

        let taken = attrs.remove_attr("material_index").unwrap();
        assert_eq!(taken.as_label().unwrap(), &[1, 2, 3, 4]);
        assert!(attrs.attr("material_index").is_none());
        assert!(attrs.remove_attr("material_index").is_none());
    }

    #[test]
    fn validate_catches_a_mismatch() -> Result<()> {
        let attrs = full_attrs();
        attrs.validate(N)?;

        assert!(attrs.validate(N + 1).is_err());
        assert!(attrs.validate(0).is_err());

        // An empty set validates against any count.
        FaceAttrSet3::empty().validate(0)?;
        FaceAttrSet3::empty().validate(9999)?;

        Ok(())
    }

    #[test]
    fn subset_selects_the_masked_faces() -> Result<()> {
        let attrs = full_attrs();
        let mask = IndexMask::try_from_indices(&[1, 3], N)?;

        let sub = attrs.subset(&mask)?;
        sub.validate(2)?;

        assert_eq!(sub.colors().unwrap(), &[[1, 1, 1], [3, 3, 3]]);
        assert_eq!(sub.labels().unwrap(), &[11, 13]);
        assert_eq!(
            sub.attr("material_index").unwrap().as_label().unwrap(),
            &[2, 4]
        );

        Ok(())
    }

    #[test]
    fn subset_indices_preserves_order_and_allows_repeats() -> Result<()> {
        let attrs = full_attrs();

        let sub = attrs.subset_indices(&[3, 0, 3])?;
        sub.validate(3)?;

        assert_eq!(sub.labels().unwrap(), &[13, 10, 13]);
        assert_eq!(sub.colors().unwrap(), &[[3, 3, 3], [0, 0, 0], [3, 3, 3]]);
        assert_eq!(
            sub.attr("material_index").unwrap().as_label().unwrap(),
            &[4, 1, 4]
        );

        assert!(attrs.subset_indices(&[0, N]).is_err());

        Ok(())
    }

    #[test]
    fn extend_from_unions_matching_sets() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.extend_from(&full_attrs())?;

        attrs.validate(N * 2)?;
        assert_eq!(attrs.labels().unwrap(), &[10, 11, 12, 13, 10, 11, 12, 13]);
        assert_eq!(attrs.attr("material_index").unwrap().len(), N * 2);

        Ok(())
    }

    #[test]
    fn extend_from_rejects_a_typed_field_on_only_one_side() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.set_labels(None, N).unwrap();

        assert!(attrs.extend_from(&other).is_err());

        // The rejected append must have left the target untouched.
        assert_eq!(attrs, full_attrs());
    }

    #[test]
    fn extend_from_rejects_an_open_key_on_only_one_side() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.remove_attr("material_index");

        assert!(attrs.extend_from(&other).is_err());
        assert_eq!(attrs, full_attrs());

        // And the reverse direction, where the other side has the extra key.
        let mut attrs = full_attrs();
        attrs.remove_attr("material_index");
        assert!(attrs.extend_from(&full_attrs()).is_err());
    }

    #[test]
    fn extend_from_rejects_a_variant_mismatch_under_the_same_name() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.remove_attr("material_index");
        other
            .insert_attr("material_index", Attr3::Scalar(vec![1.0, 2.0, 3.0, 4.0]), N)
            .unwrap();

        assert!(attrs.extend_from(&other).is_err());
        assert_eq!(attrs, full_attrs());
    }

    #[test]
    fn transform_rotates_vector_attrs_only() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.insert_attr(
            "face_dir",
            Attr3::Vector(vec![Vector3::new(1.0, 0.0, 0.0); N]),
            N,
        )?;

        // A quarter turn about +z with a translation, which maps +x onto +y.
        let iso = Iso3::new(Vector3::new(5.0, 6.0, 7.0), Vector3::z() * FRAC_PI_2);
        attrs.transform_in_place(&iso);

        assert_relative_eq!(
            attrs.attr("face_dir").unwrap().as_vector().unwrap()[0],
            Vector3::new(0.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );

        // Colors and labels are untouched by a transform.
        assert_eq!(attrs.colors(), full_attrs().colors());
        assert_eq!(attrs.labels(), full_attrs().labels());
        assert_eq!(
            attrs.attr("material_index"),
            full_attrs().attr("material_index")
        );

        Ok(())
    }
}
