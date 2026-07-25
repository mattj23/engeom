//! This module contains the per-element attribute arrays which can be attached to a `MeshData3`
//! or `Mesh3` through an attribute set.

use crate::common::IndexMask;
use crate::{Iso3, Result, Vector3};

/// A single per-element attribute array attached to a mesh. The length of the underlying vector is
/// expected to match either the point count or the face count of the mesh it belongs to, but this
/// type does not know which, and does not enforce it. Validation is the responsibility of the
/// owner, which knows the counts.
///
/// # The `Vector` variant is spatial
///
/// `MeshAttr3::Vector` is documented as holding *spatial* directions, which means it is rotated
/// when the mesh is transformed. If you need to store a non-spatial triple of numbers, store it as
/// three separate `Scalar` attributes instead, or it will be silently rotated out from under you.
#[derive(Debug, Clone, PartialEq)]
pub enum MeshAttr3 {
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

impl MeshAttr3 {
    /// Get the number of elements in the attribute array.
    pub fn len(&self) -> usize {
        match self {
            MeshAttr3::Scalar(v) => v.len(),
            MeshAttr3::Label(v) => v.len(),
            MeshAttr3::Vector(v) => v.len(),
            MeshAttr3::Color(v) => v.len(),
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
            MeshAttr3::Scalar(_) => "scalar",
            MeshAttr3::Label(_) => "label",
            MeshAttr3::Vector(_) => "vector",
            MeshAttr3::Color(_) => "color",
        }
    }

    /// Get the underlying values if this is a `Scalar` attribute, otherwise `None`.
    pub fn as_scalar(&self) -> Option<&[f64]> {
        match self {
            MeshAttr3::Scalar(v) => Some(v),
            _ => None,
        }
    }

    /// Get the underlying values if this is a `Label` attribute, otherwise `None`.
    pub fn as_label(&self) -> Option<&[u32]> {
        match self {
            MeshAttr3::Label(v) => Some(v),
            _ => None,
        }
    }

    /// Get the underlying values if this is a `Vector` attribute, otherwise `None`.
    pub fn as_vector(&self) -> Option<&[Vector3]> {
        match self {
            MeshAttr3::Vector(v) => Some(v),
            _ => None,
        }
    }

    /// Get the underlying values if this is a `Color` attribute, otherwise `None`.
    pub fn as_color(&self) -> Option<&[[u8; 3]]> {
        match self {
            MeshAttr3::Color(v) => Some(v),
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
    /// returns: `Result<MeshAttr3>`
    pub fn clone_indices_of(&self, mask: &IndexMask) -> Result<Self> {
        match self {
            MeshAttr3::Scalar(v) => Ok(MeshAttr3::Scalar(mask.clone_indices_of(v)?)),
            MeshAttr3::Label(v) => Ok(MeshAttr3::Label(mask.clone_indices_of(v)?)),
            MeshAttr3::Vector(v) => Ok(MeshAttr3::Vector(mask.clone_indices_of(v)?)),
            MeshAttr3::Color(v) => Ok(MeshAttr3::Color(mask.clone_indices_of(v)?)),
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
    /// returns: `Result<MeshAttr3>`
    pub fn clone_indices(&self, indices: &[usize]) -> Result<Self> {
        let n = self.len();
        if let Some(bad) = indices.iter().find(|&&i| i >= n) {
            return Err(
                format!("Index {bad} is out of bounds for an attribute of length {n}").into(),
            );
        }

        Ok(match self {
            MeshAttr3::Scalar(v) => MeshAttr3::Scalar(take_indices(v, indices)),
            MeshAttr3::Label(v) => MeshAttr3::Label(take_indices(v, indices)),
            MeshAttr3::Vector(v) => MeshAttr3::Vector(take_indices(v, indices)),
            MeshAttr3::Color(v) => MeshAttr3::Color(take_indices(v, indices)),
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
            (MeshAttr3::Scalar(a), MeshAttr3::Scalar(b)) => a.extend_from_slice(b),
            (MeshAttr3::Label(a), MeshAttr3::Label(b)) => a.extend_from_slice(b),
            (MeshAttr3::Vector(a), MeshAttr3::Vector(b)) => a.extend_from_slice(b),
            (MeshAttr3::Color(a), MeshAttr3::Color(b)) => a.extend_from_slice(b),
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
        if let MeshAttr3::Vector(v) = self {
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::FRAC_PI_2;

    fn sample_scalar() -> MeshAttr3 {
        MeshAttr3::Scalar(vec![0.0, 1.0, 2.0, 3.0])
    }

    fn sample_label() -> MeshAttr3 {
        MeshAttr3::Label(vec![10, 11, 12, 13])
    }

    fn sample_vector() -> MeshAttr3 {
        MeshAttr3::Vector(vec![
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 1.0),
        ])
    }

    fn sample_color() -> MeshAttr3 {
        MeshAttr3::Color(vec![[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]])
    }

    fn all_samples() -> Vec<MeshAttr3> {
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

        assert!(MeshAttr3::Scalar(Vec::new()).is_empty());
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
            MeshAttr3::Scalar(vec![1.0, 3.0])
        );
        assert_eq!(
            sample_label().clone_indices_of(&mask)?,
            MeshAttr3::Label(vec![11, 13])
        );
        assert_eq!(
            sample_color().clone_indices_of(&mask)?,
            MeshAttr3::Color(vec![[3, 4, 5], [9, 10, 11]])
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
        assert_eq!(taken, MeshAttr3::Scalar(vec![3.0, 0.0, 3.0]));

        let taken = sample_label().clone_indices(&[2, 1])?;
        assert_eq!(taken, MeshAttr3::Label(vec![12, 11]));

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
        attr.extend_from(&MeshAttr3::Scalar(vec![4.0]))?;
        assert_eq!(attr, MeshAttr3::Scalar(vec![0.0, 1.0, 2.0, 3.0, 4.0]));

        let mut attr = sample_color();
        attr.extend_from(&MeshAttr3::Color(vec![[12, 13, 14]]))?;
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
