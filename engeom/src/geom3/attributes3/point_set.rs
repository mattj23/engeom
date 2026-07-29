//! This module contains `PointAttrSet3`, the bundle of per-point attributes carried by any 3D
//! container which has a point domain.
//!
//! The bundle holds two kinds of data. A small set of first-class attributes are stored as typed
//! fields, and everything else lives in an open, name-keyed map of `Attr3` values. An attribute
//! earns a typed field if either of the following is true:
//!
//! 1. An algorithm in the library reads it in a hot loop, so it shouldn't pay for a hash lookup and
//!    a variant match.
//! 2. The container has to do something to it that it doesn't do to an arbitrary scalar, such as
//!    rotating it under a transform, scaling it, or negating it.
//!
//! # Naming and distributions
//!
//! `stdev` is named for the distribution parameter it holds, not for the general idea of
//! uncertainty, because the container's behavior depends on which one it is. A standard deviation
//! scales by `k` under a uniform scale of `k`, while a variance scales by `k²`. There is no way to
//! scale the value correctly (for example, when converting from mm to inches) without knowing which
//! is stored, so the name has to say.
//!
//! # Error message naming
//!
//! The typed fields report themselves as `point_normals`, `point_colors`, and `point_stdev` in
//! error messages rather than by their bare field names. `MeshAttrSet3` composes this type and has
//! a face domain to disambiguate against, and a message which names the domain is useful to a point
//! cloud caller too.

use super::{
    Attr3, check_both_or_neither, check_keys_match, check_len, check_reserved, check_same_variant,
    clone_masked, extend_option,
};
use crate::common::IndexMask;
use crate::{Iso3, Result, UnitVec3};
use std::collections::HashMap;

/// The set of per-point attributes attached to a container, holding both the first-class typed
/// attributes and an open, name-keyed map of everything else.
///
/// This type does not know the point count of the container it belongs to and so it can't enforce
/// that its arrays are the right length. Unfortunately that means that validation has to be the job
/// of the owner, which does know the count. To ease the pain somewhat there is a `validate(n)`
/// method, and the setters take the count as an argument.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct PointAttrSet3 {
    normals: Option<Vec<UnitVec3>>,
    colors: Option<Vec<[u8; 3]>>,
    stdev: Option<Vec<f64>>,
    open: HashMap<String, Attr3>,
}

// ===============================================================================================
// Construction and access
// ===============================================================================================

impl PointAttrSet3 {
    /// Create an attribute set with no attributes of any kind.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Returns true if no attributes of any kind are present.
    pub fn is_empty(&self) -> bool {
        self.normals.is_none()
            && self.colors.is_none()
            && self.stdev.is_none()
            && self.open.is_empty()
    }

    /// List the names of every per-point attribute which is present, typed fields included.
    ///
    /// This is what an operation which cannot supply values for a new point reports back, so the
    /// caller can see what is standing in the way.
    pub fn attr_labels(&self) -> Vec<&str> {
        let mut names = Vec::new();
        if self.normals.is_some() {
            names.push("point_normals");
        }
        if self.colors.is_some() {
            names.push("point_colors");
        }
        if self.stdev.is_some() {
            names.push("point_stdev");
        }
        names.extend(self.open.keys().map(|k| k.as_str()));
        names
    }

    /// Get the per-point unit normals, if present.
    pub fn normals(&self) -> Option<&[UnitVec3]> {
        self.normals.as_deref()
    }

    /// Get the per-point RGB colors, if present.
    pub fn colors(&self) -> Option<&[[u8; 3]]> {
        self.colors.as_deref()
    }

    /// Get the per-point standard deviations, if present.
    ///
    /// These are 1-sigma values expressed in the container's own length units, representing the
    /// spread of positions that would be expected if the point were measured repeatedly. They scale
    /// with the geometry under a uniform scale.
    pub fn stdev(&self) -> Option<&[f64]> {
        self.stdev.as_deref()
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

impl PointAttrSet3 {
    /// Set or clear the per-point unit normals.
    ///
    /// # Arguments
    ///
    /// * `values`: the normals to store, or `None` to clear them
    /// * `n_points`: the point count of the owning container, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_normals(&mut self, values: Option<Vec<UnitVec3>>, n_points: usize) -> Result<()> {
        check_len(
            values.as_deref().map(|v| v.len()),
            n_points,
            "point_normals",
        )?;
        self.normals = values;
        Ok(())
    }

    /// Set or clear the per-point RGB colors.
    ///
    /// # Arguments
    ///
    /// * `values`: the colors to store, or `None` to clear them
    /// * `n_points`: the point count of the owning container, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_colors(&mut self, values: Option<Vec<[u8; 3]>>, n_points: usize) -> Result<()> {
        check_len(values.as_deref().map(|v| v.len()), n_points, "point_colors")?;
        self.colors = values;
        Ok(())
    }

    /// Set or clear the per-point standard deviations, which must be 1-sigma values in the
    /// container's own length units. Negative values are rejected.
    ///
    /// # Arguments
    ///
    /// * `values`: the standard deviations to store, or `None` to clear them
    /// * `n_points`: the point count of the owning container, which `values` must match
    ///
    /// returns: `Result<()>`
    pub fn set_stdev(&mut self, values: Option<Vec<f64>>, n_points: usize) -> Result<()> {
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

        self.stdev = values;
        Ok(())
    }

    /// Insert an open-map attribute under the given name, replacing any attribute already stored
    /// there.
    ///
    /// # Arguments
    ///
    /// * `name`: the key to store under, which must not be one of the reserved names
    /// * `attr`: the attribute array to store
    /// * `n_points`: the point count of the owning container, which `attr` must match
    ///
    /// returns: `Result<()>`
    pub fn insert_attr(&mut self, name: &str, attr: Attr3, n_points: usize) -> Result<()> {
        check_reserved(name)?;
        check_len(Some(attr.len()), n_points, name)?;
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

impl PointAttrSet3 {
    /// Verify that every attribute array present has a length matching the point count.
    ///
    /// # Arguments
    ///
    /// * `n_points`: the point count of the owning container
    ///
    /// returns: `Result<()>`
    pub fn validate(&self, n_points: usize) -> Result<()> {
        check_len(
            self.normals.as_ref().map(|v| v.len()),
            n_points,
            "point_normals",
        )?;
        check_len(
            self.colors.as_ref().map(|v| v.len()),
            n_points,
            "point_colors",
        )?;
        check_len(
            self.stdev.as_ref().map(|v| v.len()),
            n_points,
            "point_stdev",
        )?;

        for (name, attr) in self.open.iter() {
            check_len(Some(attr.len()), n_points, name)?;
        }

        Ok(())
    }

    /// Create a new attribute set containing only the points selected by the given mask, preserving
    /// their original order.
    ///
    /// # Arguments
    ///
    /// * `mask`: selects which points survive, and must match the owning container's point count
    ///
    /// returns: `Result<PointAttrSet3>`
    pub fn subset(&self, mask: &IndexMask) -> Result<Self> {
        let mut open = HashMap::with_capacity(self.open.len());
        for (name, attr) in self.open.iter() {
            open.insert(name.clone(), attr.clone_indices_of(mask)?);
        }

        Ok(Self {
            normals: clone_masked(self.normals.as_deref(), mask)?,
            colors: clone_masked(self.colors.as_deref(), mask)?,
            stdev: clone_masked(self.stdev.as_deref(), mask)?,
            open,
        })
    }

    /// Create a new attribute set containing the points at the given indices, in the order the
    /// indices are given. Indices may repeat.
    ///
    /// # Arguments
    ///
    /// * `indices`: the points to take, each of which must be less than the owning container's
    ///   point count
    ///
    /// returns: `Result<PointAttrSet3>`
    pub fn subset_indices(&self, indices: &[usize]) -> Result<Self> {
        let mut open = HashMap::with_capacity(self.open.len());
        for (name, attr) in self.open.iter() {
            open.insert(name.clone(), attr.clone_indices(indices)?);
        }

        Ok(Self {
            normals: clone_indexed(self.normals.as_deref(), indices)?,
            colors: clone_indexed(self.colors.as_deref(), indices)?,
            stdev: clone_indexed(self.stdev.as_deref(), indices)?,
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
    /// This is separate from the append itself so that `MeshAttrSet3` can check both of its domains
    /// before modifying either. Without it, a face-domain mismatch discovered after the point
    /// domain had already been extended would leave the set half appended.
    pub(crate) fn check_extend_from(&self, other: &Self) -> Result<()> {
        check_both_or_neither(
            self.normals.is_some(),
            other.normals.is_some(),
            "point_normals",
        )?;
        check_both_or_neither(
            self.colors.is_some(),
            other.colors.is_some(),
            "point_colors",
        )?;
        check_both_or_neither(self.stdev.is_some(), other.stdev.is_some(), "point_stdev")?;

        check_keys_match(&self.open, &other.open, "point")?;

        // Variant mismatches between two attributes of the same name would fail partway through the
        // append, so they are checked up front as well.
        for (name, attr) in self.open.iter() {
            check_same_variant(attr, &other.open[name], name)?;
        }

        Ok(())
    }

    /// Perform the append. Only valid after `check_extend_from` has succeeded for the same pair.
    pub(crate) fn apply_extend_from(&mut self, other: &Self) -> Result<()> {
        extend_option(&mut self.normals, other.normals.as_deref());
        extend_option(&mut self.colors, other.colors.as_deref());
        extend_option(&mut self.stdev, other.stdev.as_deref());

        for (name, attr) in self.open.iter_mut() {
            attr.extend_from(&other.open[name])?;
        }

        Ok(())
    }
}

// ===============================================================================================
// Geometric operations
// ===============================================================================================

impl PointAttrSet3 {
    /// Transform the spatial attributes in place by the given isometry.
    ///
    /// This rotates the normals and every `Attr3::Vector` attribute. Because all of these hold
    /// directions rather than positions, only the rotation component of the isometry has any
    /// effect.
    ///
    /// # Arguments
    ///
    /// * `iso`: the isometry to apply
    pub fn transform_in_place(&mut self, iso: &Iso3) {
        if let Some(normals) = &mut self.normals {
            for n in normals.iter_mut() {
                *n = iso * *n;
            }
        }

        for attr in self.open.values_mut() {
            attr.transform_in_place(iso);
        }
    }

    /// Scale the length-dimensioned attributes in place by the given uniform scale factor.
    ///
    /// Only `stdev` is affected, and it is scaled by the absolute value of the factor, since a
    /// standard deviation is a magnitude and must stay non-negative even under a mirroring scale.
    /// Nothing else is touched: normals are directions and are unchanged by a uniform scale, and
    /// the container has no way to know whether an open-map scalar carries length units.
    ///
    /// # Arguments
    ///
    /// * `scale`: the uniform scale factor which was applied to the points
    pub fn scale_in_place(&mut self, scale: f64) {
        if let Some(stdev) = &mut self.stdev {
            let k = scale.abs();
            for s in stdev.iter_mut() {
                *s *= k;
            }
        }
    }

    /// Flip the orientation-dependent attributes in place, to accompany a reversal of the surface
    /// the points belong to. Only the normals are negated.
    pub fn flip_in_place(&mut self) {
        if let Some(normals) = &mut self.normals {
            for n in normals.iter_mut() {
                *n = -*n;
            }
        }
    }
}

/// Select the elements of an optional attribute array at the given indices, in the order given.
fn clone_indexed<T: Clone>(values: Option<&[T]>, indices: &[usize]) -> Result<Option<Vec<T>>> {
    let Some(values) = values else {
        return Ok(None);
    };

    let n = values.len();
    if let Some(bad) = indices.iter().find(|&&i| i >= n) {
        return Err(format!("Index {bad} is out of bounds for an attribute of length {n}").into());
    }

    Ok(Some(indices.iter().map(|&i| values[i].clone()).collect()))
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
    fn full_attrs() -> PointAttrSet3 {
        let mut attrs = PointAttrSet3::empty();

        attrs
            .set_normals(
                Some(vec![
                    UnitVec3::new_normalize(Vector3::new(1.0, 0.0, 0.0)),
                    UnitVec3::new_normalize(Vector3::new(0.0, 1.0, 0.0)),
                    UnitVec3::new_normalize(Vector3::new(0.0, 0.0, 1.0)),
                    UnitVec3::new_normalize(Vector3::new(1.0, 1.0, 1.0)),
                ]),
                N,
            )
            .unwrap();

        attrs
            .set_colors(Some(vec![[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]]), N)
            .unwrap();

        attrs.set_stdev(Some(vec![0.1, 0.2, 0.3, 0.4]), N).unwrap();

        attrs
            .insert_attr("confidence", Attr3::Scalar(vec![0.5, 0.6, 0.7, 0.8]), N)
            .unwrap();

        attrs
    }

    #[test]
    fn empty_set_reports_empty() {
        let attrs = PointAttrSet3::empty();
        assert!(attrs.is_empty());
        assert!(attrs.normals().is_none());
        assert!(attrs.colors().is_none());
        assert!(attrs.stdev().is_none());
        assert_eq!(attrs.attr_names().count(), 0);

        assert!(!full_attrs().is_empty());
    }

    #[test]
    fn attr_labels_names_every_present_attribute() {
        let attrs = full_attrs();
        let labels = attrs.attr_labels();

        for expected in ["point_normals", "point_colors", "point_stdev", "confidence"] {
            assert!(
                labels.contains(&expected),
                "missing '{expected}' in {labels:?}"
            );
        }

        assert!(PointAttrSet3::empty().attr_labels().is_empty());
    }

    #[test]
    fn setters_reject_a_length_mismatch() {
        let mut attrs = PointAttrSet3::empty();

        assert!(attrs.set_colors(Some(vec![[0, 0, 0]]), N).is_err());
        assert!(attrs.set_normals(Some(Vec::new()), N).is_err());
        assert!(attrs.insert_attr("q", Attr3::Scalar(vec![1.0]), N).is_err());

        // A rejected set must not have stored anything.
        assert!(attrs.is_empty());
    }

    #[test]
    fn setters_accept_none_regardless_of_count() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.set_colors(None, N)?;
        assert!(attrs.colors().is_none());
        Ok(())
    }

    #[test]
    fn stdev_rejects_negative_and_nan() {
        let mut attrs = PointAttrSet3::empty();

        assert!(attrs.set_stdev(Some(vec![0.1, -0.2, 0.3, 0.4]), N).is_err());
        assert!(
            attrs
                .set_stdev(Some(vec![0.1, f64::NAN, 0.3, 0.4]), N)
                .is_err()
        );
        assert!(attrs.set_stdev(Some(vec![0.0, 0.2, 0.3, 0.4]), N).is_ok());
    }

    #[test]
    fn open_map_rejects_reserved_names() {
        let mut attrs = PointAttrSet3::empty();

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
                .insert_attr("scanner_color_temp", Attr3::Scalar(vec![0.0; N]), N)
                .is_ok()
        );
    }

    #[test]
    fn remove_attr_returns_what_was_stored() {
        let mut attrs = full_attrs();

        let taken = attrs.remove_attr("confidence").unwrap();
        assert_eq!(taken.as_scalar().unwrap(), &[0.5, 0.6, 0.7, 0.8]);
        assert!(attrs.attr("confidence").is_none());
        assert!(attrs.remove_attr("confidence").is_none());
    }

    #[test]
    fn validate_catches_a_mismatch() -> Result<()> {
        let attrs = full_attrs();
        attrs.validate(N)?;

        assert!(attrs.validate(N + 1).is_err());
        assert!(attrs.validate(0).is_err());

        // An empty set validates against any count.
        PointAttrSet3::empty().validate(0)?;
        PointAttrSet3::empty().validate(9999)?;

        Ok(())
    }

    #[test]
    fn subset_selects_the_masked_points() -> Result<()> {
        let attrs = full_attrs();
        let mask = IndexMask::try_from_indices(&[1, 3], N)?;

        let sub = attrs.subset(&mask)?;
        sub.validate(2)?;

        assert_eq!(sub.colors().unwrap(), &[[1, 1, 1], [3, 3, 3]]);
        assert_eq!(sub.stdev().unwrap(), &[0.2, 0.4]);
        assert_eq!(
            sub.attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.6, 0.8]
        );

        Ok(())
    }

    #[test]
    fn extend_from_unions_matching_sets() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.extend_from(&full_attrs())?;

        attrs.validate(N * 2)?;
        assert_eq!(attrs.stdev().unwrap().len(), N * 2);
        assert_eq!(attrs.attr("confidence").unwrap().len(), N * 2);

        Ok(())
    }

    #[test]
    fn extend_from_rejects_a_typed_field_on_only_one_side() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.set_stdev(None, N).unwrap();

        assert!(attrs.extend_from(&other).is_err());

        // The rejected append must have left the target untouched.
        assert_eq!(attrs, full_attrs());
    }

    #[test]
    fn extend_from_rejects_an_open_key_on_only_one_side() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.remove_attr("confidence");

        assert!(attrs.extend_from(&other).is_err());
        assert_eq!(attrs, full_attrs());

        // And the reverse direction, where the other side has the extra key.
        let mut attrs = full_attrs();
        attrs.remove_attr("confidence");
        assert!(attrs.extend_from(&full_attrs()).is_err());
    }

    #[test]
    fn extend_from_rejects_a_variant_mismatch_under_the_same_name() {
        let mut attrs = full_attrs();
        let mut other = full_attrs();
        other.remove_attr("confidence");
        other
            .insert_attr("confidence", Attr3::Label(vec![1, 2, 3, 4]), N)
            .unwrap();

        assert!(attrs.extend_from(&other).is_err());
        assert_eq!(attrs, full_attrs());
    }

    #[test]
    fn transform_rotates_normals_and_vector_attrs() -> Result<()> {
        let mut attrs = full_attrs();
        attrs.insert_attr(
            "principal_dir",
            Attr3::Vector(vec![Vector3::new(1.0, 0.0, 0.0); N]),
            N,
        )?;

        // A quarter turn about +z with a translation, which maps +x onto +y.
        let iso = Iso3::new(Vector3::new(5.0, 6.0, 7.0), Vector3::z() * FRAC_PI_2);
        attrs.transform_in_place(&iso);

        assert_relative_eq!(
            attrs.normals().unwrap()[0].into_inner(),
            Vector3::new(0.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );
        assert_relative_eq!(
            attrs.attr("principal_dir").unwrap().as_vector().unwrap()[0],
            Vector3::new(0.0, 1.0, 0.0),
            epsilon = 1.0e-12
        );

        // Scalars are untouched by a transform.
        assert_eq!(attrs.stdev().unwrap(), &[0.1, 0.2, 0.3, 0.4]);
        assert_eq!(
            attrs.attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.5, 0.6, 0.7, 0.8]
        );

        Ok(())
    }

    #[test]
    fn scale_affects_only_the_standard_deviations() {
        let mut attrs = full_attrs();
        attrs.scale_in_place(25.4);

        for (actual, expected) in attrs
            .stdev()
            .unwrap()
            .iter()
            .zip([0.1, 0.2, 0.3, 0.4].iter())
        {
            assert_relative_eq!(*actual, expected * 25.4, epsilon = 1.0e-12);
        }

        // An open-map scalar has no declared units, so it must not be scaled.
        assert_eq!(
            attrs.attr("confidence").unwrap().as_scalar().unwrap(),
            &[0.5, 0.6, 0.7, 0.8]
        );
    }

    #[test]
    fn scale_keeps_standard_deviations_positive_under_mirroring() {
        let mut attrs = full_attrs();
        attrs.scale_in_place(-2.0);

        for s in attrs.stdev().unwrap() {
            assert!(*s >= 0.0, "a standard deviation must not go negative");
        }
        assert_relative_eq!(attrs.stdev().unwrap()[1], 0.4, epsilon = 1.0e-12);
    }

    #[test]
    fn flip_negates_the_normals() {
        let original = full_attrs();
        let mut attrs = full_attrs();
        attrs.flip_in_place();

        for (flipped, source) in attrs
            .normals()
            .unwrap()
            .iter()
            .zip(original.normals().unwrap().iter())
        {
            assert_relative_eq!(
                flipped.into_inner(),
                -source.into_inner(),
                epsilon = 1.0e-12
            );
        }

        // Nothing else changes.
        assert_eq!(attrs.colors(), original.colors());
        assert_eq!(attrs.stdev(), original.stdev());
    }
}
