//! This module contains the per-element attribute arrays which can be attached to any of the 3D
//! geometry containers through an attribute set: `MeshData3` and `Mesh3` on either the point or the
//! face domain, and `PointCloudData3` on the point domain.
//!
//! It lives here, outside `mesh`, because nothing about the array itself is specific to a mesh. The
//! attribute set types which own these arrays and know the element counts to validate against are
//! the ones which belong to a particular container.
//!
//! `PointAttrSet3` also lives here. The point domain is common to every container, so it is written
//! once here and composed by `MeshAttrSet3`, which adds the face domain on top of it. The validation
//! helpers at the bottom of this file are shared by both for the same reason.

mod point_set;

pub use point_set::PointAttrSet3;

use crate::common::IndexMask;
use crate::{Iso3, Result, Vector3};
use std::collections::HashMap;

/// Attribute names which may not be used as keys in the open attribute maps, because they name a
/// quantity that either already has a typed field or is computed on demand. As a precaution we're
/// going to reject them to prevent any quantity from having two homes that can silently disagree.
pub(crate) const RESERVED_ATTR_NAMES: [&str; 8] = [
    "normal", "normals", "color", "colors", "stdev", "std_dev", "label", "labels",
];

/// A single per-element attribute array. The length of the underlying vector is expected to match
/// the element count of the domain it is attached to, but this type does not know which domain that
/// is, and does not enforce it. Validation is the responsibility of the owning attribute set, which
/// knows the counts.
///
/// # The `Vector` variant is spatial
///
/// `Attr3::Vector` is documented as holding *spatial* directions, which means it is rotated
/// when the geometry is transformed. If you need to store a non-spatial triple of numbers, store it
/// as three separate `Scalar` attributes instead, or it will be silently rotated out from under you.
#[derive(Debug, Clone, PartialEq)]
pub enum Attr3 {
    /// A floating point value per element, such as a curvature, a deviation, or a confidence.
    Scalar(Vec<f64>),

    /// An unsigned integer per element, such as a region/patch identifier, a scan pass number, or
    /// a material index.
    Label(Vec<u32>),

    /// A spatial direction per element. Rotated by `transform_in_place`.
    Vector(Vec<Vector3>),

    /// An 8-bit RGB color per element.
    Color(Vec<[u8; 3]>),
}

impl Attr3 {
    /// Get the number of elements in the attribute array.
    pub fn len(&self) -> usize {
        match self {
            Attr3::Scalar(v) => v.len(),
            Attr3::Label(v) => v.len(),
            Attr3::Vector(v) => v.len(),
            Attr3::Color(v) => v.len(),
        }
    }

    /// Returns true if the attribute array has no elements.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Get a short, lowercase name for this attribute's variant, suitable for use in error
    /// messages.
    pub fn kind(&self) -> &'static str {
        match self {
            Attr3::Scalar(_) => "scalar",
            Attr3::Label(_) => "label",
            Attr3::Vector(_) => "vector",
            Attr3::Color(_) => "color",
        }
    }

    /// Get the underlying values if this is a `Scalar` attribute, otherwise `None`.
    pub fn as_scalar(&self) -> Option<&[f64]> {
        match self {
            Attr3::Scalar(v) => Some(v),
            _ => None,
        }
    }

    /// Get the underlying values if this is a `Label` attribute, otherwise `None`.
    pub fn as_label(&self) -> Option<&[u32]> {
        match self {
            Attr3::Label(v) => Some(v),
            _ => None,
        }
    }

    /// Get the underlying values if this is a `Vector` attribute, otherwise `None`.
    pub fn as_vector(&self) -> Option<&[Vector3]> {
        match self {
            Attr3::Vector(v) => Some(v),
            _ => None,
        }
    }

    /// Get the underlying values if this is a `Color` attribute, otherwise `None`.
    pub fn as_color(&self) -> Option<&[[u8; 3]]> {
        match self {
            Attr3::Color(v) => Some(v),
            _ => None,
        }
    }

    /// Create a new attribute array containing only the elements whose indices are marked `true`
    /// in the given mask, preserving their original order.
    ///
    /// # Arguments
    ///
    /// * `mask`: a mask whose length must match the length of this attribute array
    ///
    /// returns: `Result<Attr3>`
    pub fn clone_indices_of(&self, mask: &IndexMask) -> Result<Self> {
        match self {
            Attr3::Scalar(v) => Ok(Attr3::Scalar(mask.clone_indices_of(v)?)),
            Attr3::Label(v) => Ok(Attr3::Label(mask.clone_indices_of(v)?)),
            Attr3::Vector(v) => Ok(Attr3::Vector(mask.clone_indices_of(v)?)),
            Attr3::Color(v) => Ok(Attr3::Color(mask.clone_indices_of(v)?)),
        }
    }

    /// Create a new attribute array containing the elements at the given indices, in the order the
    /// indices are given. Indices may repeat.
    ///
    /// # Arguments
    ///
    /// * `indices`: the indices to take, each of which must be less than the length of this
    ///   attribute array
    ///
    /// returns: `Result<Attr3>`
    pub fn clone_indices(&self, indices: &[usize]) -> Result<Self> {
        let n = self.len();
        if let Some(bad) = indices.iter().find(|&&i| i >= n) {
            return Err(
                format!("Index {bad} is out of bounds for an attribute of length {n}").into(),
            );
        }

        Ok(match self {
            Attr3::Scalar(v) => Attr3::Scalar(take_indices(v, indices)),
            Attr3::Label(v) => Attr3::Label(take_indices(v, indices)),
            Attr3::Vector(v) => Attr3::Vector(take_indices(v, indices)),
            Attr3::Color(v) => Attr3::Color(take_indices(v, indices)),
        })
    }

    /// Append the contents of another attribute array onto the end of this one. Both attributes
    /// must be the same variant, otherwise an error is returned and this attribute is left
    /// unmodified.
    ///
    /// # Arguments
    ///
    /// * `other`: the attribute array to append
    ///
    /// returns: `Result<()>`
    pub fn extend_from(&mut self, other: &Self) -> Result<()> {
        match (self, other) {
            (Attr3::Scalar(a), Attr3::Scalar(b)) => a.extend_from_slice(b),
            (Attr3::Label(a), Attr3::Label(b)) => a.extend_from_slice(b),
            (Attr3::Vector(a), Attr3::Vector(b)) => a.extend_from_slice(b),
            (Attr3::Color(a), Attr3::Color(b)) => a.extend_from_slice(b),
            (a, b) => {
                return Err(format!(
                    "Cannot append a {} attribute onto a {} attribute",
                    b.kind(),
                    a.kind()
                )
                .into());
            }
        }

        Ok(())
    }

    /// Transform this attribute array in place by the given isometry.
    ///
    /// Only the `Vector` variant is affected, and only by the rotation component of the isometry,
    /// because its values are directions rather than positions. All other variants are left
    /// unchanged.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    pub fn transform_in_place(&mut self, iso: &Iso3) {
        if let Attr3::Vector(v) = self {
            for value in v.iter_mut() {
                *value = iso * *value;
            }
        }
    }
}

/// Clone the elements of `items` at the given indices, in the order the indices are given. The
/// caller is responsible for having verified that the indices are in bounds.
fn take_indices<T: Clone>(items: &[T], indices: &[usize]) -> Vec<T> {
    indices.iter().map(|&i| items[i].clone()).collect()
}

// ===============================================================================================
// Shared attribute set helpers
// ===============================================================================================
//
// These are used by both `PointAttrSet3` and `MeshAttrSet3`. They live here rather than in either
// one so that the two cannot drift apart on what counts as a valid attribute or a legal append.

/// Verify that an optional attribute length matches the expected element count.
pub(crate) fn check_len(actual: Option<usize>, expected: usize, name: &str) -> Result<()> {
    match actual {
        Some(n) if n != expected => Err(format!(
            "Attribute '{name}' has {n} values, but there are {expected} elements in that domain"
        )
        .into()),
        _ => Ok(()),
    }
}

/// Verify that a name is not one of the reserved open-map keys.
pub(crate) fn check_reserved(name: &str) -> Result<()> {
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
pub(crate) fn check_both_or_neither(a: bool, b: bool, name: &str) -> Result<()> {
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
pub(crate) fn check_keys_match(
    a: &HashMap<String, Attr3>,
    b: &HashMap<String, Attr3>,
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
pub(crate) fn check_same_variant(a: &Attr3, b: &Attr3, name: &str) -> Result<()> {
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
pub(crate) fn clone_masked<T: Clone>(
    values: Option<&[T]>,
    mask: &IndexMask,
) -> Result<Option<Vec<T>>> {
    values.map(|v| mask.clone_indices_of(v)).transpose()
}

/// Append the contents of an optional attribute array onto another, where both are known to be
/// present or both absent.
pub(crate) fn extend_option<T: Clone>(target: &mut Option<Vec<T>>, source: Option<&[T]>) {
    if let (Some(t), Some(s)) = (target.as_mut(), source) {
        t.extend_from_slice(s);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::FRAC_PI_2;

    fn sample_scalar() -> Attr3 {
        Attr3::Scalar(vec![0.0, 1.0, 2.0, 3.0])
    }

    fn sample_label() -> Attr3 {
        Attr3::Label(vec![10, 11, 12, 13])
    }

    fn sample_vector() -> Attr3 {
        Attr3::Vector(vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
        ])
    }

    fn sample_color() -> Attr3 {
        Attr3::Color(vec![[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])
    }

    fn all_samples() -> Vec<Attr3> {
        vec![
            sample_scalar(),
            sample_label(),
            sample_vector(),
            sample_color(),
        ]
    }

    #[test]
    fn len_matches_underlying_vector() {
        for attr in all_samples() {
            assert_eq!(attr.len(), 4);
            assert!(!attr.is_empty());
        }

        assert!(Attr3::Scalar(Vec::new()).is_empty());
    }

    #[test]
    fn accessors_return_only_their_own_variant() {
        assert!(sample_scalar().as_scalar().is_some());
        assert!(sample_scalar().as_label().is_none());
        assert!(sample_scalar().as_vector().is_none());
        assert!(sample_scalar().as_color().is_none());

        assert!(sample_label().as_label().is_some());
        assert!(sample_label().as_scalar().is_none());

        assert!(sample_vector().as_vector().is_some());
        assert!(sample_vector().as_scalar().is_none());

        assert!(sample_color().as_color().is_some());
        assert!(sample_color().as_scalar().is_none());
    }

    #[test]
    fn clone_indices_of_selects_masked_elements() -> Result<()> {
        let mask = IndexMask::try_from_indices(&[1, 3], 4)?;

        assert_eq!(
            sample_scalar().clone_indices_of(&mask)?,
            Attr3::Scalar(vec![1.0, 3.0])
        );
        assert_eq!(
            sample_label().clone_indices_of(&mask)?,
            Attr3::Label(vec![11, 13])
        );
        assert_eq!(
            sample_color().clone_indices_of(&mask)?,
            Attr3::Color(vec![[3, 4, 5], [9, 10, 11]])
        );

        let vectors = sample_vector().clone_indices_of(&mask)?;
        assert_eq!(
            vectors.as_vector().unwrap(),
            &[Vector3::new(0.0, 1.0, 0.0), Vector3::new(1.0, 1.0, 1.0)]
        );

        Ok(())
    }

    #[test]
    fn clone_indices_of_rejects_a_mismatched_mask() {
        let mask = IndexMask::new(5, true);
        for attr in all_samples() {
            assert!(attr.clone_indices_of(&mask).is_err());
        }
    }

    #[test]
    fn clone_indices_preserves_order_and_allows_repeats() -> Result<()> {
        let taken = sample_scalar().clone_indices(&[3, 0, 3])?;
        assert_eq!(taken, Attr3::Scalar(vec![3.0, 0.0, 3.0]));

        let taken = sample_label().clone_indices(&[2, 1])?;
        assert_eq!(taken, Attr3::Label(vec![12, 11]));

        Ok(())
    }

    #[test]
    fn clone_indices_rejects_out_of_bounds() {
        for attr in all_samples() {
            assert!(attr.clone_indices(&[0, 4]).is_err());
        }
    }

    #[test]
    fn extend_from_appends_matching_variants() -> Result<()> {
        let mut attr = sample_scalar();
        attr.extend_from(&Attr3::Scalar(vec![4.0]))?;
        assert_eq!(attr, Attr3::Scalar(vec![0.0, 1.0, 2.0, 3.0, 4.0]));

        let mut attr = sample_color();
        attr.extend_from(&Attr3::Color(vec![[12, 13, 14]]))?;
        assert_eq!(attr.len(), 5);

        Ok(())
    }

    #[test]
    fn extend_from_rejects_a_variant_mismatch() {
        let mut attr = sample_scalar();
        let result = attr.extend_from(&sample_label());

        assert!(result.is_err());

        // The failed append must leave the original untouched.
        assert_eq!(attr, sample_scalar());
    }

    #[test]
    fn transform_rotates_only_the_vector_variant() {
        // A quarter turn about +z, which maps +x onto +y.
        let iso = Iso3::rotation(Vector3::z() * FRAC_PI_2);

        let mut vectors = sample_vector();
        vectors.transform_in_place(&iso);
        let values = vectors.as_vector().unwrap();
        assert_relative_eq!(values[0], Vector3::new(0.0, 1.0, 0.0), epsilon = 1.0e-12);
        assert_relative_eq!(values[1], Vector3::new(-1.0, 0.0, 0.0), epsilon = 1.0e-12);
        assert_relative_eq!(values[2], Vector3::new(0.0, 0.0, 1.0), epsilon = 1.0e-12);

        for mut attr in [sample_scalar(), sample_label(), sample_color()] {
            let original = attr.clone();
            attr.transform_in_place(&iso);
            assert_eq!(attr, original);
        }
    }

    #[test]
    fn transform_ignores_the_translation_component() {
        // Directions must not pick up the translation of the isometry.
        let iso = Iso3::translation(10.0, 20.0, 30.0);

        let mut vectors = sample_vector();
        let original = vectors.clone();
        vectors.transform_in_place(&iso);

        assert_eq!(vectors, original);
    }
}
