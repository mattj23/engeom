//! This module contains `MeshAttrSet3`, the bundle of per-element attributes carried by a
//! `MeshData3`.
//!
//! The bundle holds two kinds of data. A small set of first-class attributes are stored as typed
//! fields, and everything else lives in an open, name-keyed map of `MeshAttr3` values. An attribute
//! earns a typed field if either of the following is true:
//!
//! 1. An algorithm in the library reads it in a hot loop, so it shouldn't pay for a hash lookup and
//!    a variant match.
//! 2. The container has to do something to it that it doesn't do to an arbitrary scalar, such as
//!    rotating it under a transform, scaling it, or negating it.
//!
//! # Naming and distributions
//!
//! `point_stdev` is named for the distribution parameter it holds, not for the general idea of
//! uncertainty, because the container's behavior depends on which one it is. A standard deviation
//! scales by `k` under a uniform scale of `k`, while a variance scales by `k²`. There is no way to
//! scale the value correctly (for example, when converting a mesh from mm to inches) without
//! knowing which is stored, so the name has to say.

use super::MeshAttr3;
use crate::common::IndexMask;
use crate::{Iso3, Result, UnitVec3};
use std::collections::HashMap;

/// Attribute names which may not be used as keys in the open maps, because they name a quantity
/// that either already has a typed field or is computed on demand. As a precaution we're going to
/// reject them to prevent any quantity from having two homes that can silently disagree.
const RESERVED_ATTR_NAMES: [&str; 8] = [
    "normal", "normals", "color", "colors", "stdev", "std_dev", "label", "labels",
];

/// The set of per-element attributes attached to a mesh, holding both the first-class typed
/// attributes and an open, name-keyed map of everything else.
///
/// This type does not know the point or face count of the mesh it belongs to and so it can't
/// enforce that its arrays are the right length. Unfortunately that means that validation has to
/// be the job of the owner, which does know the counts. To ease the pain somewhat I have created
/// the `validate(n_points, n_faces)` method.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MeshAttrSet3 {
    point_normals: Option<Vec<UnitVec3>>,
    point_colors: Option<Vec<[u8; 3]>>,
    point_stdev: Option<Vec<f64>>,

    face_colors: Option<Vec<[u8; 3]>>,
    face_labels: Option<Vec<u32>>,

    point_attrs: HashMap<String, MeshAttr3>,
    face_attrs: HashMap<String, MeshAttr3>,
}

// ===============================================================================================
// Construction and access
// ===============================================================================================

impl MeshAttrSet3 {
    /// Create an attribute set with no attributes of any kind.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Returns true if no attributes of any kind are present.
    pub fn is_empty(&self) -> bool {
        self.point_normals.is_none()
            && self.point_colors.is_none()
            && self.point_stdev.is_none()
            && self.face_colors.is_none()
            && self.face_labels.is_none()
            && self.point_attrs.is_empty()
            && self.face_attrs.is_empty()
    }

    /// List the names of every per-point attribute which is present, typed fields included.
    ///
    /// This is what an operation which cannot supply values for a new point reports back, so the
    /// caller can see what is standing in the way.
    pub fn point_attr_labels(&self) -> Vec<&str> {
        let mut names = Vec::new();
        if self.point_normals.is_some() {
            names.push("point_normals");
        }
        if self.point_colors.is_some() {
            names.push("point_colors");
        }
        if self.point_stdev.is_some() {
            names.push("point_stdev");
        }
        names.extend(self.point_attrs.keys().map(|k| k.as_str()));
        names
    }

    /// List the names of every per-face attribute which is present, typed fields included.
    pub fn face_attr_labels(&self) -> Vec<&str> {
        let mut names = Vec::new();
        if self.face_colors.is_some() {
            names.push("face_colors");
        }
        if self.face_labels.is_some() {
            names.push("face_labels");
        }
        names.extend(self.face_attrs.keys().map(|k| k.as_str()));
        names
    }

    /// Get the per-point unit normals, if present.
    pub fn point_normals(&self) -> Option<&[UnitVec3]> {
        self.point_normals.as_deref()
    }

    /// Get the per-point RGB colors, if present.
    pub fn point_colors(&self) -> Option<&[[u8; 3]]> {
        self.point_colors.as_deref()
    }

    /// Get the per-point standard deviations, if present.
    ///
    /// These are 1-sigma values expressed in the mesh's own length units, representing the spread
    /// of positions that would be expected if the point were measured repeatedly. They scale with
    /// the mesh under a uniform scale.
    pub fn point_stdev(&self) -> Option<&[f64]> {
        self.point_stdev.as_deref()
    }

    /// Get the per-face RGB colors, if present.
    pub fn face_colors(&self) -> Option<&[[u8; 3]]> {
        self.face_colors.as_deref()
    }

    /// Get the per-face labels, if present. These identify which region, patch, scan pass, or
    /// material each face belongs to.
    pub fn face_labels(&self) -> Option<&[u32]> {
        self.face_labels.as_deref()
    }

    /// Get the open-map per-point attribute stored under the given name, if present.
    pub fn point_attr(&self, name: &str) -> Option<&MeshAttr3> {
        self.point_attrs.get(name)
    }

    /// Get the open-map per-face attribute stored under the given name, if present.
    pub fn face_attr(&self, name: &str) -> Option<&MeshAttr3> {
        self.face_attrs.get(name)
    }

    /// Iterate over the names of the open-map per-point attributes, in no particular order.
    pub fn point_attr_names(&self) -> impl Iterator<Item = &str> {
        self.point_attrs.keys().map(|k| k.as_str())
    }

    /// Iterate over the names of the open-map per-face attributes, in no particular order.
    pub fn face_attr_names(&self) -> impl Iterator<Item = &str> {
        self.face_attrs.keys().map(|k| k.as_str())
    }
}

// ===============================================================================================
// Mutation
// ===============================================================================================

impl MeshAttrSet3 {
    /// Set or clear the per-point unit normals.
    ///
    /// # Arguments
    ///
    /// * `values`: the normals to store, or `None` to clear them
    /// * `n_points`: the point count of the owning mesh, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_point_normals(
        &mut self,
        values: Option<Vec<UnitVec3>>,
        n_points: usize,
    ) -> Result<()> {
        check_len(
            values.as_deref().map(|v| v.len()),
            n_points,
            "point_normals",
        )?;
        self.point_normals = values;
        Ok(())
    }

    /// Set or clear the per-point RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them
    /// * `n_points`: the point count of the owning mesh, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_point_colors(
        &mut self,
        values: Option<Vec<[u8; 3]>>,
        n_points: usize,
    ) -> Result<()> {
        check_len(values.as_deref().map(|v| v.len()), n_points, "point_colors")?;
        self.point_colors = values;
        Ok(())
    }

    /// Set or clear the per-point standard deviations, which must be 1-sigma values in the mesh's
    /// own length units. Negative values are rejected.
    ///
    /// # Arguments
    ///
    /// * `values`: the standard deviations to store, or `None` to clear them
    /// * `n_points`: the point count of the owning mesh, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_point_stdev(&mut self, values: Option<Vec<f64>>, n_points: usize) -> Result<()> {
        check_len(values.as_deref().map(|v| v.len()), n_points, "point_stdev")?;

        if let Some(v) = &values
            && let Some(i) = v.iter().position(|s| *s < 0.0 || s.is_nan())
        {
            return Err(format!(
                "point_stdev[{}] is {}, but a standard deviation must be finite and non-negative",
                i, v[i]
            )
            .into());
        }

        self.point_stdev = values;
        Ok(())
    }

    /// Set or clear the per-face RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them
    /// * `n_faces`: the face count of the owning mesh, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_face_colors(&mut self, values: Option<Vec<[u8; 3]>>, n_faces: usize) -> Result<()> {
        check_len(values.as_deref().map(|v| v.len()), n_faces, "face_colors")?;
        self.face_colors = values;
        Ok(())
    }

    /// Set or clear the per-face labels.
    ///
    /// # Arguments
    ///
    /// * `values`: the labels to store, or `None` to clear them
    /// * `n_faces`: the face count of the owning mesh, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_face_labels(&mut self, values: Option<Vec<u32>>, n_faces: usize) -> Result<()> {
        check_len(values.as_deref().map(|v| v.len()), n_faces, "face_labels")?;
        self.face_labels = values;
        Ok(())
    }

    /// Insert an open-map per-point attribute under the given name, replacing any attribute
    /// already stored there.
    ///
    /// # Arguments
    ///
    /// * `name`: the key to store under, which must not be one of the reserved names
    /// * `attr`: the attribute array to store
    /// * `n_points`: the point count of the owning mesh, which `attr` must match
    ///
    /// returns: `Result<()>`
    pub fn insert_point_attr(
        &mut self,
        name: &str,
        attr: MeshAttr3,
        n_points: usize,
    ) -> Result<()> {
        check_reserved(name)?;
        check_len(Some(attr.len()), n_points, name)?;
        self.point_attrs.insert(name.to_string(), attr);
        Ok(())
    }

    /// Insert an open-map per-face attribute under the given name, replacing any attribute already
    /// stored there.
    ///
    /// # Arguments
    ///
    /// * `name`: the key to store under, which must not be one of the reserved names
    /// * `attr`: the attribute array to store
    /// * `n_faces`: the face count of the owning mesh, which `attr` must match
    ///
    /// returns: `Result<()>`
    pub fn insert_face_attr(&mut self, name: &str, attr: MeshAttr3, n_faces: usize) -> Result<()> {
        check_reserved(name)?;
        check_len(Some(attr.len()), n_faces, name)?;
        self.face_attrs.insert(name.to_string(), attr);
        Ok(())
    }

    /// Remove and return the open-map per-point attribute stored under the given name.
    pub fn remove_point_attr(&mut self, name: &str) -> Option<MeshAttr3> {
        self.point_attrs.remove(name)
    }

    /// Remove and return the open-map per-face attribute stored under the given name.
    pub fn remove_face_attr(&mut self, name: &str) -> Option<MeshAttr3> {
        self.face_attrs.remove(name)
    }
}

// ===============================================================================================
// Validation, subsetting, and combination
// ===============================================================================================

impl MeshAttrSet3 {
    /// Verify that every attribute array present has a length matching the element count of the
    /// domain it belongs to.
    ///
    /// # Arguments
    ///
    /// * `n_points`: the point count of the owning mesh
    /// * `n_faces`: the face count of the owning mesh
    ///
    /// returns: `Result<()>`
    pub fn validate(&self, n_points: usize, n_faces: usize) -> Result<()> {
        check_len(
            self.point_normals.as_ref().map(|v| v.len()),
            n_points,
            "point_normals",
        )?;
        check_len(
            self.point_colors.as_ref().map(|v| v.len()),
            n_points,
            "point_colors",
        )?;
        check_len(
            self.point_stdev.as_ref().map(|v| v.len()),
            n_points,
            "point_stdev",
        )?;
        check_len(
            self.face_colors.as_ref().map(|v| v.len()),
            n_faces,
            "face_colors",
        )?;
        check_len(
            self.face_labels.as_ref().map(|v| v.len()),
            n_faces,
            "face_labels",
        )?;

        for (name, attr) in self.point_attrs.iter() {
            check_len(Some(attr.len()), n_points, name)?;
        }

        for (name, attr) in self.face_attrs.iter() {
            check_len(Some(attr.len()), n_faces, name)?;
        }

        Ok(())
    }

    /// Create a new attribute set containing only the elements selected by the given masks,
    /// preserving their original order.
    ///
    /// # Arguments
    ///
    /// * `point_mask`: selects which points survive, and must match the owning mesh's point
    ///   count
    /// * `face_mask`: selects which faces survive, and must match the owning mesh's face count
    ///
    /// returns: `Result<MeshAttrSet3>`
    pub fn subset(&self, point_mask: &IndexMask, face_mask: &IndexMask) -> Result<Self> {
        let mut point_attrs = HashMap::with_capacity(self.point_attrs.len());
        for (name, attr) in self.point_attrs.iter() {
            point_attrs.insert(name.clone(), attr.clone_indices_of(point_mask)?);
        }

        let mut face_attrs = HashMap::with_capacity(self.face_attrs.len());
        for (name, attr) in self.face_attrs.iter() {
            face_attrs.insert(name.clone(), attr.clone_indices_of(face_mask)?);
        }

        Ok(Self {
            point_normals: clone_masked(self.point_normals.as_deref(), point_mask)?,
            point_colors: clone_masked(self.point_colors.as_deref(), point_mask)?,
            point_stdev: clone_masked(self.point_stdev.as_deref(), point_mask)?,
            face_colors: clone_masked(self.face_colors.as_deref(), face_mask)?,
            face_labels: clone_masked(self.face_labels.as_deref(), face_mask)?,
            point_attrs,
            face_attrs,
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
        check_both_or_neither(
            self.point_normals.is_some(),
            other.point_normals.is_some(),
            "point_normals",
        )?;
        check_both_or_neither(
            self.point_colors.is_some(),
            other.point_colors.is_some(),
            "point_colors",
        )?;
        check_both_or_neither(
            self.point_stdev.is_some(),
            other.point_stdev.is_some(),
            "point_stdev",
        )?;
        check_both_or_neither(
            self.face_colors.is_some(),
            other.face_colors.is_some(),
            "face_colors",
        )?;
        check_both_or_neither(
            self.face_labels.is_some(),
            other.face_labels.is_some(),
            "face_labels",
        )?;

        check_keys_match(&self.point_attrs, &other.point_attrs, "point")?;
        check_keys_match(&self.face_attrs, &other.face_attrs, "face")?;

        // Variant mismatches between two attributes of the same name would fail partway through
        // the append below, so they are checked up front as well.
        for (name, attr) in self.point_attrs.iter() {
            check_same_variant(attr, &other.point_attrs[name], name)?;
        }
        for (name, attr) in self.face_attrs.iter() {
            check_same_variant(attr, &other.face_attrs[name], name)?;
        }

        extend_option(&mut self.point_normals, other.point_normals.as_deref());
        extend_option(&mut self.point_colors, other.point_colors.as_deref());
        extend_option(&mut self.point_stdev, other.point_stdev.as_deref());
        extend_option(&mut self.face_colors, other.face_colors.as_deref());
        extend_option(&mut self.face_labels, other.face_labels.as_deref());

        for (name, attr) in self.point_attrs.iter_mut() {
            attr.extend_from(&other.point_attrs[name])?;
        }
        for (name, attr) in self.face_attrs.iter_mut() {
            attr.extend_from(&other.face_attrs[name])?;
        }

        Ok(())
    }
}

// ===============================================================================================
// Geometric operations
// ===============================================================================================

impl MeshAttrSet3 {
    /// Transform the spatial attributes in place by the given isometry.
    ///
    /// This rotates the point normals and every per-point `MeshAttr3::Vector` attribute. Because
    /// all of these hold directions rather than positions, only the rotation component of the
    /// isometry has any effect. Per-face vector attributes are rotated as well.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    pub fn transform_in_place(&mut self, iso: &Iso3) {
        if let Some(normals) = &mut self.point_normals {
            for n in normals.iter_mut() {
                *n = iso * *n;
            }
        }

        for attr in self.point_attrs.values_mut() {
            attr.transform_in_place(iso);
        }

        for attr in self.face_attrs.values_mut() {
            attr.transform_in_place(iso);
        }
    }

    /// Scale the length-dimensioned attributes in place by the given uniform scale factor.
    ///
    /// Only `point_stdev` is affected, and it is scaled by the absolute value of the factor, since
    /// a standard deviation is a magnitude and must stay non-negative even under a mirroring scale.
    /// Nothing else is touched: normals are directions and are unchanged by a uniform scale, and
    /// the container has no way to know whether an open-map scalar carries length units.
    ///
    /// # Arguments
    ///
    /// * `scale`: the uniform scale factor which was applied to the mesh points
    pub fn scale_in_place(&mut self, scale: f64) {
        if let Some(stdev) = &mut self.point_stdev {
            let k = scale.abs();
            for s in stdev.iter_mut() {
                *s *= k;
            }
        }
    }

    /// Flip the orientation-dependent attributes in place, to accompany a reversal of the mesh's
    /// face winding. Only the point normals are negated.
    pub fn flip_in_place(&mut self) {
        if let Some(normals) = &mut self.point_normals {
            for n in normals.iter_mut() {
                *n = -*n;
            }
        }
    }
}

// ===============================================================================================
// Helpers
// ===============================================================================================

/// Verify that an optional attribute length matches the expected element count.
fn check_len(actual: Option<usize>, expected: usize, name: &str) -> Result<()> {
    match actual {
        Some(n) if n != expected => Err(format!(
            "Attribute '{name}' has {n} values, but the mesh has {expected} elements in that domain"
        )
        .into()),
        _ => Ok(()),
    }
}

/// Verify that a name is not one of the reserved open-map keys.
fn check_reserved(name: &str) -> Result<()> {
    if RESERVED_ATTR_NAMES.contains(&name) {
        return Err(format!(
            "'{name}' is a reserved attribute name. Quantities which have a typed field or are \
             computed on demand (normals, colors, standard deviations, labels) must be set through \
             their own accessor, so that they cannot have two homes which disagree."
        )
        .into());
    }

    Ok(())
}

/// Verify that a typed attribute is either present on both sides of an append or absent on both.
fn check_both_or_neither(a: bool, b: bool, name: &str) -> Result<()> {
    if a != b {
        let (has, lacks) = if a {
            ("this", "the other")
        } else {
            ("the other", "this")
        };
        return Err(format!(
            "Cannot append: {has} attribute set has '{name}' but {lacks} one does not"
        )
        .into());
    }

    Ok(())
}

/// Verify that two open maps hold exactly the same set of keys.
fn check_keys_match(
    a: &HashMap<String, MeshAttr3>,
    b: &HashMap<String, MeshAttr3>,
    domain: &str,
) -> Result<()> {
    for name in a.keys() {
        if !b.contains_key(name) {
            return Err(format!(
                "Cannot append: this attribute set has a {domain} attribute '{name}' but the other \
                 one does not"
            )
            .into());
        }
    }

    for name in b.keys() {
        if !a.contains_key(name) {
            return Err(format!(
                "Cannot append: the other attribute set has a {domain} attribute '{name}' but this \
                 one does not"
            )
            .into());
        }
    }

    Ok(())
}

/// Verify that two attributes stored under the same name hold the same variant.
fn check_same_variant(a: &MeshAttr3, b: &MeshAttr3, name: &str) -> Result<()> {
    if a.kind() != b.kind() {
        return Err(format!(
            "Cannot append: attribute '{}' is a {} on this side and a {} on the other",
            name,
            a.kind(),
            b.kind()
        )
        .into());
    }

    Ok(())
}

/// Select the masked elements of an optional attribute array.
fn clone_masked<T: Clone>(values: Option<&[T]>, mask: &IndexMask) -> Result<Option<Vec<T>>> {
    values.map(|v| mask.clone_indices_of(v)).transpose()
}

/// Append the contents of an optional attribute array onto another, where both are known to be
/// present or both absent.
fn extend_option<T: Clone>(target: &mut Option<Vec<T>>, source: Option<&[T]>) {
    if let (Some(t), Some(s)) = (target.as_mut(), source) {
        t.extend_from_slice(s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vector3;
    use approx::assert_relative_eq;
    use std::f64::consts::FRAC_PI_2;

    const N_POINTS: usize = 4;
    const N_FACES: usize = 2;

    /// An attribute set carrying every typed field and one open attribute in each domain.
    fn full_attrs() -> MeshAttrSet3 {
        let mut attrs = MeshAttrSet3::empty();

        attrs
            .set_point_normals(
                Some(vec![
                    UnitVec3::new_normalize(Vector3::new(1.0, 0.0, 0.0)),
                    UnitVec3::new_normalize(Vector3::new(0.0, 1.0, 0.0)),
                    UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0)),
                    UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 1.0)),
                ]),
                N_POINTS,
            )
            .unwrap();

        attrs
            .set_point_colors(
                Some(vec![[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]]),
                N_POINTS,
            )
            .unwrap();

        attrs
            .set_point_stdev(Some(vec![0.1, 0.2, 0.3, 0.4]), N_POINTS)
            .unwrap();

        attrs
            .set_face_colors(Some(vec![[10, 10, 10], [20, 20, 20]]), N_FACES)
            .unwrap();

        attrs.set_face_labels(Some(vec![7, 8]), N_FACES).unwrap();

        attrs
            .insert_point_attr(
                "confidence",
                MeshAttr3::Scalar(vec![0.5, 0.6, 0.7, 0.8]),
                N_POINTS,
            )
            .unwrap();

        attrs
            .insert_face_attr("material_index", MeshAttr3::Label(vec![1, 2]), N_FACES)
            .unwrap();

        attrs
    }

    #[test]
    fn empty_set_reports_empty() {
        let attrs = MeshAttrSet3::empty();
        assert!(attrs.is_empty());
        assert!(attrs.point_normals().is_none());
        assert!(attrs.face_labels().is_none());
        assert_eq!(attrs.point_attr_names().count(), 0);

        assert!(!full_attrs().is_empty());
    }

    #[test]
    fn setters_reject_a_length_mismatch() {
        let mut attrs = MeshAttrSet3::empty();

        assert!(
            attrs
                .set_point_colors(Some(vec![[0, 0, 0]]), N_POINTS)
                .is_err()
        );
        assert!(attrs.set_face_labels(Some(vec![1, 2, 3]), N_FACES).is_err());
        assert!(
            attrs
                .insert_point_attr("q", MeshAttr3::Scalar(vec![1.0]), N_POINTS)
                .is_err()
        );

        // A rejected set must not have stored anything.
        assert!(attrs.is_empty());
    }

    #[test]
    fn setters_accept_none_regardless_of_count() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.set_point_colors(None, N_POINTS)?;
        assert!(attrs.point_colors().is_none());
        Ok(())
    }

    #[test]
    fn stdev_rejects_negative_and_nan() {
        let mut attrs = MeshAttrSet3::empty();

        assert!(
            attrs
                .set_point_stdev(Some(vec![0.1, -0.2, 0.3, 0.4]), N_POINTS)
                .is_err()
        );
        assert!(
            attrs
                .set_point_stdev(Some(vec![0.1, f64::NAN, 0.3, 0.4]), N_POINTS)
                .is_err()
        );
        assert!(
            attrs
                .set_point_stdev(Some(vec![0.0, 0.2, 0.3, 0.4]), N_POINTS)
                .is_ok()
        );
    }

    #[test]
    fn open_maps_reject_reserved_names() {
        let mut attrs = MeshAttrSet3::empty();

        for name in RESERVED_ATTR_NAMES {
            assert!(
                attrs
                    .insert_point_attr(name, MeshAttr3::Scalar(vec![0.0; N_POINTS]), N_POINTS)
                    .is_err(),
                "expected '{name}' to be rejected as a point attribute name"
            );
            assert!(
                attrs
                    .insert_face_attr(name, MeshAttr3::Scalar(vec![0.0; N_FACES]), N_FACES)
                    .is_err(),
                "expected '{name}' to be rejected as a face attribute name"
            );
        }

        // A name which merely contains a reserved word is fine.
        assert!(
            attrs
                .insert_point_attr(
                    "scanner_color_temp",
                    MeshAttr3::Scalar(vec![0.0; N_POINTS]),
                    N_POINTS
                )
                .is_ok()
        );
    }

    #[test]
    fn validate_catches_a_mismatch_in_either_domain() -> Result<()> {
        let attrs = full_attrs();
        attrs.validate(N_POINTS, N_FACES)?;

        assert!(attrs.validate(N_POINTS + 1, N_FACES).is_err());
        assert!(attrs.validate(N_POINTS, N_FACES + 1).is_err());

        Ok(())
    }

    #[test]
    fn subset_selects_from_the_correct_domain() -> Result<()> {
        let attrs = full_attrs();
        let point_mask = IndexMask::try_from_indices(&[1, 3], N_POINTS)?;
        let face_mask = IndexMask::try_from_indices(&[0], N_FACES)?;

        let sub = attrs.subset(&point_mask, &face_mask)?;
        sub.validate(2, 1)?;

        assert_eq!(sub.point_colors().unwrap(), &[[1, 1, 1], [3, 3, 3]]);
        assert_eq!(sub.point_stdev().unwrap(), &[0.2, 0.4]);
        assert_eq!(sub.face_colors().unwrap(), &[[10, 10, 10]]);
        assert_eq!(sub.face_labels().unwrap(), &[7]);

        assert_eq!(
            sub.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.6, 0.8]
        );
        assert_eq!(
            sub.face_attr("material_index").unwrap().as_label().unwrap(),
            &[1]
        );

        Ok(())
    }

    #[test]
    fn extend_from_unions_matching_sets() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.extend_from(&full_attrs())?;

        attrs.validate(N_POINTS * 2, N_FACES * 2)?;
        assert_eq!(attrs.point_stdev().unwrap().len(), N_POINTS * 2);
        assert_eq!(attrs.face_labels().unwrap(), &[7, 8, 7, 8]);
        assert_eq!(attrs.point_attr("confidence").unwrap().len(), N_POINTS * 2);

        Ok(())
    }

    #[test]
    fn extend_from_rejects_a_typed_field_on_only_one_side() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.set_point_stdev(None, N_POINTS).unwrap();

        assert!(attrs.extend_from(&other).is_err());

        // The rejected append must have left the target untouched.
        assert_eq!(attrs, full_attrs());
    }

    #[test]
    fn extend_from_rejects_an_open_key_on_only_one_side() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.remove_point_attr("confidence");

        assert!(attrs.extend_from(&other).is_err());
        assert_eq!(attrs, full_attrs());

        // And the reverse direction, where the other side has the extra key.
        let mut attrs = full_attrs();
        attrs.remove_point_attr("confidence");
        assert!(attrs.extend_from(&full_attrs()).is_err());
    }

    #[test]
    fn extend_from_rejects_a_variant_mismatch_under_the_same_name() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.remove_point_attr("confidence");
        other
            .insert_point_attr("confidence", MeshAttr3::Label(vec![1, 2, 3, 4]), N_POINTS)
            .unwrap();

        assert!(attrs.extend_from(&other).is_err());
        assert_eq!(attrs, full_attrs());
    }

    #[test]
    fn transform_rotates_normals_and_vector_attrs() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.insert_point_attr(
            "principal_dir",
            MeshAttr3::Vector(vec![Vector3::new(1.0, 0.0, 0.0); N_POINTS]),
            N_POINTS,
        )?;

        // A quarter turn about +z with a translation, which maps +x onto +y.
        let iso = Iso3::new(Vector3::new(5.0, 6.0, 7.0), Vector3::z() * FRAC_PI_2);
        attrs.transform_in_place(&iso);

        assert_relative_eq!(
            attrs.point_normals().unwrap()[0].into_inner(),
            Vector3::new(0.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );
        assert_relative_eq!(
            attrs
                .point_attr("principal_dir")
                .unwrap()
                .as_vector()
                .unwrap()[0],
            Vector3::new(0.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );

        // Scalars are untouched by a transform.
        assert_eq!(attrs.point_stdev().unwrap(), &[0.1, 0.2, 0.3, 0.4]);
        assert_eq!(
            attrs.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.5, 0.6, 0.7, 0.8]
        );

        Ok(())
    }

    #[test]
    fn scale_affects_only_the_standard_deviations() {
        let mut attrs = full_attrs();
        attrs.scale_in_place(25.4);

        for (actual, expected) in attrs
            .point_stdev()
            .unwrap()
            .iter()
            .zip([0.1, 0.2, 0.3, 0.4].iter())
        {
            assert_relative_eq!(*actual, expected * 25.4, epsilon = 1.0e-12);
        }

        // An open-map scalar has no declared units, so it must not be scaled.
        assert_eq!(
            attrs.point_attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.5, 0.6, 0.7, 0.8]
        );
    }

    #[test]
    fn scale_keeps_standard_deviations_positive_under_mirroring() {
        let mut attrs = full_attrs();
        attrs.scale_in_place(-2.0);

        for s in attrs.point_stdev().unwrap() {
            assert!(*s >= 0.0, "a standard deviation must not go negative");
        }
        assert_relative_eq!(attrs.point_stdev().unwrap()[1], 0.4, epsilon = 1.0e-12);
    }

    #[test]
    fn flip_negates_the_point_normals() {
        let original = full_attrs();
        let mut attrs = full_attrs();
        attrs.flip_in_place();

        for (flipped, source) in attrs
            .point_normals()
            .unwrap()
            .iter()
            .zip(original.point_normals().unwrap().iter())
        {
            assert_relative_eq!(
                flipped.into_inner(),
                -source.into_inner(),
                epsilon = 1.0e-12
            );
        }

        // Nothing else changes.
        assert_eq!(attrs.point_colors(), original.point_colors());
        assert_eq!(attrs.point_stdev(), original.point_stdev());
    }
}
