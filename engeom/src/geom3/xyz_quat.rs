use crate::na::{Quaternion, Translation3, UnitQuaternion};
use crate::{Iso3, Vector3};
use serde::{Deserialize, Serialize};

/// A struct representing a 6D pose as a translation plus a rotation quaternion, in a flat,
/// self-describing form suitable for serialization.
///
/// This exists because `Iso3` is a type alias over a `nalgebra` type, and so inherits `nalgebra`'s
/// `serde` derive. That derive emits bare arrays with no field names, ordered rotation first:
///
/// ```json
/// {"rotation":[0.30,-0.12,-0.74,0.58],"translation":[-80.6,14.4,23.6]}
/// ```
///
/// Nothing in that output says whether the rotation is stored `ijkw` or `wijk`, and both the
/// ordering and the transparent representation of `Unit<Quaternion>` are `nalgebra` implementation
/// details rather than a stable format. Anything written to disk and expected to be readable later
/// should go through this type instead:
///
/// ```json
/// {"tx":-80.6,"ty":14.4,"tz":23.6,"i":0.30,"j":-0.12,"k":-0.74,"w":0.58}
/// ```
///
/// The field order matches the arguments of the Python binding's `Iso3.from_quaternion`, so the
/// Rust and Python surfaces describe a pose the same way.
///
/// # Normalization
///
/// Converting to an `Iso3` normalizes the quaternion rather than requiring it to be unit length,
/// which is what makes `From` infallible. A quaternion of zero length has no rotation to recover
/// and will produce a non-finite result; use [`XyzQuat::is_valid`] when the values come from an
/// untrusted source.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct XyzQuat {
    /// The x component of the translation.
    pub tx: f64,

    /// The y component of the translation.
    pub ty: f64,

    /// The z component of the translation.
    pub tz: f64,

    /// The coefficient of the quaternion's first imaginary basis element.
    pub i: f64,

    /// The coefficient of the quaternion's second imaginary basis element.
    pub j: f64,

    /// The coefficient of the quaternion's third imaginary basis element.
    pub k: f64,

    /// The real (scalar) part of the quaternion.
    pub w: f64,
}

impl XyzQuat {
    pub fn new(tx: f64, ty: f64, tz: f64, i: f64, j: f64, k: f64, w: f64) -> Self {
        XyzQuat {
            tx,
            ty,
            tz,
            i,
            j,
            k,
            w,
        }
    }

    /// Returns true if every component is finite and the quaternion has enough magnitude to
    /// normalize into a rotation.
    pub fn is_valid(&self) -> bool {
        self.to_array().iter().all(|v| v.is_finite())
            && (self.i * self.i + self.j * self.j + self.k * self.k + self.w * self.w) > 1e-12
    }

    /// Returns true if this pose and `other` describe the same rigid transform to within
    /// `epsilon`.
    ///
    /// The comparison is made on the transform matrices rather than component by component,
    /// because a quaternion and its negation describe the same rotation.
    pub fn approx_eq(&self, other: &XyzQuat, epsilon: f64) -> bool {
        let m0 = Iso3::from(self).to_matrix();
        let m1 = Iso3::from(other).to_matrix();
        let diff = m0 - m1;
        diff.amax() < epsilon
    }

    pub fn to_array(&self) -> [f64; 7] {
        [self.tx, self.ty, self.tz, self.i, self.j, self.k, self.w]
    }
}

impl std::fmt::Display for XyzQuat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "XyzQuat {{ tx: {}, ty: {}, tz: {}, i: {}, j: {}, k: {}, w: {} }}",
            self.tx, self.ty, self.tz, self.i, self.j, self.k, self.w
        )
    }
}

impl From<[f64; 7]> for XyzQuat {
    fn from(value: [f64; 7]) -> Self {
        XyzQuat::new(
            value[0], value[1], value[2], value[3], value[4], value[5], value[6],
        )
    }
}

impl From<&Iso3> for XyzQuat {
    fn from(isometry: &Iso3) -> Self {
        let t = isometry.translation.vector;
        let q = isometry.rotation.quaternion();
        XyzQuat::new(t.x, t.y, t.z, q.i, q.j, q.k, q.w)
    }
}

impl From<&XyzQuat> for Iso3 {
    fn from(value: &XyzQuat) -> Self {
        let translation = Vector3::new(value.tx, value.ty, value.tz);

        // `Quaternion::new` takes the real part first, unlike the field order of `XyzQuat`, which
        // follows the convention of putting the translation and the imaginary parts first.
        let quat = Quaternion::new(value.w, value.i, value.j, value.k);
        let rotation = UnitQuaternion::from_quaternion(quat);

        Iso3::from_parts(Translation3::from(translation), rotation)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::common::random_geometry::RandomGeometry3;
    use approx::assert_relative_eq;

    /// Round trips a set of random isometries through `XyzQuat` and back, comparing the transform
    /// matrices.
    ///
    /// Comparing matrices rather than components sidesteps the quaternion double cover, where `q`
    /// and `-q` describe the same rotation.
    #[test]
    fn iso3_round_trip_preserves_the_transform() {
        let mut rg = RandomGeometry3::from_seed(0xb0a7_51de);

        for _ in 0..2000 {
            let iso = rg.iso3(100.0);
            let record = XyzQuat::from(&iso);
            let back = Iso3::from(&record);

            assert_relative_eq!(iso.to_matrix(), back.to_matrix(), epsilon = 1e-12);
        }
    }

    /// A record built from a unit quaternion comes back component for component, not merely up to
    /// sign. This is what makes a serialized pose diffable against the value that produced it.
    #[test]
    fn a_unit_record_round_trips_component_wise() {
        let mut rg = RandomGeometry3::from_seed(0x51de_b0a7);

        for _ in 0..2000 {
            let record = XyzQuat::from(&rg.iso3(100.0));
            let back = XyzQuat::from(&Iso3::from(&record));

            for (a, b) in record.to_array().iter().zip(back.to_array().iter()) {
                assert_relative_eq!(a, b, epsilon = 1e-12);
            }
        }
    }

    /// A non-unit quaternion is normalized on the way in, so scaling all four components leaves
    /// the rotation unchanged.
    #[test]
    fn a_non_unit_quaternion_is_normalized() {
        let unit = XyzQuat::new(1.0, 2.0, 3.0, 0.5, 0.5, 0.5, 0.5);
        let scaled = XyzQuat::new(1.0, 2.0, 3.0, 1.5, 1.5, 1.5, 1.5);

        assert!(unit.approx_eq(&scaled, 1e-12));
    }

    /// The serialized form must be a labeled object with these exact seven keys. The private
    /// alignment data on disk is written this way, so a change here silently breaks reading it.
    #[test]
    fn serialization_uses_labeled_fields() {
        let record = XyzQuat::new(-80.6, 14.4, 23.6, 0.30, -0.12, -0.74, 0.58);
        let text = serde_json::to_string(&record).unwrap();

        assert_eq!(
            text,
            r#"{"tx":-80.6,"ty":14.4,"tz":23.6,"i":0.3,"j":-0.12,"k":-0.74,"w":0.58}"#
        );
    }

    #[test]
    fn deserialization_reads_labeled_fields() {
        let text = r#"{"tx":-80.6,"ty":14.4,"tz":23.6,"i":0.3,"j":-0.12,"k":-0.74,"w":0.58}"#;
        let record: XyzQuat = serde_json::from_str(text).unwrap();

        assert_relative_eq!(record.tx, -80.6, epsilon = 1e-12);
        assert_relative_eq!(record.ty, 14.4, epsilon = 1e-12);
        assert_relative_eq!(record.tz, 23.6, epsilon = 1e-12);
        assert_relative_eq!(record.i, 0.3, epsilon = 1e-12);
        assert_relative_eq!(record.j, -0.12, epsilon = 1e-12);
        assert_relative_eq!(record.k, -0.74, epsilon = 1e-12);
        assert_relative_eq!(record.w, 0.58, epsilon = 1e-12);
    }

    #[test]
    fn a_degenerate_quaternion_is_rejected() {
        assert!(!XyzQuat::new(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0).is_valid());
        assert!(!XyzQuat::new(f64::NAN, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0).is_valid());
        assert!(XyzQuat::new(1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 1.0).is_valid());
    }
}
