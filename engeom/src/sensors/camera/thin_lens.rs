//! This module provides a thin-lens model of an imaging lens. This simple optical model predicts
//! the quantities commonly used to design an imaging system:
//! magnification, field of view, depth of field, and defocus blur.

use crate::Result;

/// The wavelength of green light, 550 nm, expressed in millimeters. This is the conventional
/// wavelength to use when one value represents visible light. Pass it to
/// `ThinLens::airy_disk_diameter` when the rest of the model uses millimeters.
pub const WAVELENGTH_550NM_MM: f64 = 0.55e-3;

/// A thin lens with a circular aperture, focused so that objects at `focus_distance` are imaged
/// sharply onto the sensor plane.
///
/// The thin-lens model collapses the whole optical assembly to a single plane, so all distances
/// are measured from that plane along the optical axis: `focus_distance` is on the object side
/// and the image distance returned by `image_distance` is on the sensor side. Real lenses have
/// two separated principal planes, so distances measured against the physical body of a lens
/// can differ from this model. That offset does not affect the predicted magnification, blur, or
/// depth of field.
///
/// The caller chooses the length unit, and every derived length uses that same unit. A wavelength
/// passed to `airy_disk_diameter` must also use that unit. For millimeters, use
/// `WAVELENGTH_550NM_MM`.
///
/// A lens focused at infinity is represented by a `focus_distance` of `f64::INFINITY`, and every
/// method handles that case without special treatment by the caller.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThinLens {
    /// The focal length of the lens.
    pub focal_length: f64,

    /// The f-number of the lens at infinity focus, which is the focal length divided by the
    /// diameter of the entrance pupil.
    pub f_number: f64,

    /// The object-side distance from the lens plane to the plane of best focus. May be
    /// `f64::INFINITY`.
    pub focus_distance: f64,
}

impl ThinLens {
    /// Create a new thin lens from its focal length, f-number, and focus distance.
    ///
    /// # Arguments
    ///
    /// * `focal_length`: the focal length of the lens, which must be a finite number greater
    ///   than zero
    /// * `f_number`: the f-number of the lens at infinity focus, which must be a finite number
    ///   greater than zero
    /// * `focus_distance`: the object-side distance to the plane of best focus, which must be
    ///   greater than the focal length. `f64::INFINITY` is permitted and means the lens is
    ///   focused at infinity; distances at or below the focal length are rejected because they
    ///   do not form a real image on the sensor side.
    ///
    /// returns: Result<ThinLens>
    pub fn new(focal_length: f64, f_number: f64, focus_distance: f64) -> Result<Self> {
        if !focal_length.is_finite() || focal_length <= 0.0 {
            return Err("ThinLens focal_length must be a finite number greater than zero".into());
        }
        if !f_number.is_finite() || f_number <= 0.0 {
            return Err("ThinLens f_number must be a finite number greater than zero".into());
        }
        if focus_distance.is_nan() {
            return Err("ThinLens focus_distance must not be NaN".into());
        }
        if focus_distance <= focal_length {
            return Err(
                "ThinLens focus_distance must be greater than the focal length; use \
                 f64::INFINITY for a lens focused at infinity"
                    .into(),
            );
        }
        Ok(Self {
            focal_length,
            f_number,
            focus_distance,
        })
    }

    /// Create a thin lens focused to achieve a specified magnification, as commonly used to
    /// specify a macro imaging system.
    ///
    /// # Arguments
    ///
    /// * `focal_length`: the focal length of the lens, which must be a finite number greater
    ///   than zero
    /// * `f_number`: the f-number of the lens at infinity focus, which must be a finite number
    ///   greater than zero
    /// * `magnification`: the ratio of image size to object size, which must be a finite number
    ///   greater than zero
    ///
    /// returns: Result<ThinLens>
    pub fn from_magnification(
        focal_length: f64,
        f_number: f64,
        magnification: f64,
    ) -> Result<Self> {
        if !magnification.is_finite() || magnification <= 0.0 {
            return Err("ThinLens magnification must be a finite number greater than zero".into());
        }
        let focus_distance = focal_length * (1.0 + magnification) / magnification;
        Self::new(focal_length, f_number, focus_distance)
    }

    /// Create a thin lens whose focal length fills a sensor of a given width with a specified
    /// width of the object plane at a given working distance. Use this constructor to calculate
    /// the focal length from the part width and working distance.
    ///
    /// # Arguments
    ///
    /// * `f_number`: the f-number of the lens at infinity focus, which must be a finite number
    ///   greater than zero
    /// * `focus_distance`: the object-side working distance, which must be a finite number
    ///   greater than zero
    /// * `sensor_width`: the physical width of the sensor, which must be a finite number
    ///   greater than zero
    /// * `field_width`: the width of the object plane to be imaged across that sensor width,
    ///   which must be a finite number greater than zero
    ///
    /// returns: Result<ThinLens>
    pub fn from_field_width(
        f_number: f64,
        focus_distance: f64,
        sensor_width: f64,
        field_width: f64,
    ) -> Result<Self> {
        if !focus_distance.is_finite() || focus_distance <= 0.0 {
            return Err(
                "ThinLens focus_distance must be a finite number greater than zero when solving \
                 for a focal length"
                    .into(),
            );
        }
        if !sensor_width.is_finite() || sensor_width <= 0.0 {
            return Err("ThinLens sensor_width must be a finite number greater than zero".into());
        }
        if !field_width.is_finite() || field_width <= 0.0 {
            return Err("ThinLens field_width must be a finite number greater than zero".into());
        }
        let focal_length = focus_distance * sensor_width / (field_width + sensor_width);
        Self::new(focal_length, f_number, focus_distance)
    }

    /// Create a copy of this lens refocused to a different distance, leaving the focal length
    /// and f-number unchanged.
    ///
    /// # Arguments
    ///
    /// * `focus_distance`: the new object-side distance to the plane of best focus, subject to
    ///   the same rules as in `ThinLens::new`
    ///
    /// returns: Result<ThinLens>
    pub fn focused_at(&self, focus_distance: f64) -> Result<Self> {
        Self::new(self.focal_length, self.f_number, focus_distance)
    }

    /// Whether this lens is focused at infinity.
    pub fn is_focused_at_infinity(&self) -> bool {
        self.focus_distance.is_infinite()
    }

    /// The image-side distance from the lens plane to the sensor plane, calculated from the
    /// thin-lens equation `1 / f = 1 / s_o + 1 / s_i`. A lens focused at infinity returns its
    /// focal length; any closer focus returns something longer, which is why a lens focused close
    /// behaves as though it were longer than its nameplate focal length.
    pub fn image_distance(&self) -> f64 {
        // Computed through reciprocals so that an infinite focus distance falls out as 1/f
        // instead of producing an indeterminate form.
        1.0 / (1.0 / self.focal_length - 1.0 / self.focus_distance)
    }

    /// The magnification of the lens at its focus distance, which is the ratio of image size to
    /// object size. A lens focused at infinity has a magnification of zero.
    pub fn magnification(&self) -> f64 {
        self.image_distance() / self.focus_distance
    }

    /// The effective f-number of the lens at its focus distance, `N * (1 + m)`. The f-number
    /// marked on a lens applies at infinity focus; as the lens is focused closer the cone of
    /// light reaching the sensor narrows, and at 1:1 magnification the effective f-number is
    /// twice the marked one. This is the value that governs both diffraction and exposure.
    pub fn effective_f_number(&self) -> f64 {
        self.f_number * (1.0 + self.magnification())
    }

    /// The diameter of the entrance pupil, which is the focal length divided by the f-number.
    pub fn aperture_diameter(&self) -> f64 {
        self.focal_length / self.f_number
    }

    /// The hyperfocal distance for an acceptable circle of confusion: the focus distance beyond
    /// which the far limit of the depth of field reaches infinity. This does not depend on where
    /// the lens is currently focused.
    ///
    /// # Arguments
    ///
    /// * `coc`: the largest circle of confusion diameter on the sensor plane that is still
    ///   considered acceptably sharp. A value of zero or less returns `f64::INFINITY`, since no
    ///   finite focus distance can hold infinity in focus under a zero blur tolerance.
    ///
    /// returns: f64
    pub fn hyperfocal_distance(&self, coc: f64) -> f64 {
        if coc.is_nan() || coc <= 0.0 {
            return f64::INFINITY;
        }
        self.focal_length * self.focal_length / (self.f_number * coc) + self.focal_length
    }

    /// The diameter of the Airy disk on the sensor plane, measured to the first zero of the
    /// diffraction pattern. The calculation uses the effective f-number to account for high
    /// magnification. This diameter is the lower limit on the image size of a point, regardless
    /// of focus. Comparing it with the pixel pitch indicates whether diffraction or sampling
    /// limits the system.
    ///
    /// # Arguments
    ///
    /// * `wavelength`: the wavelength of the light, in the same length unit as the rest of the
    ///   lens. `WAVELENGTH_550NM_MM` is the usual choice for visible light in millimeters.
    ///
    /// returns: f64
    pub fn airy_disk_diameter(&self, wavelength: f64) -> f64 {
        2.44 * wavelength * self.effective_f_number()
    }

    /// The diameter of the circle of confusion on the sensor plane produced by a point at a
    /// given object distance. A point at the focus distance produces zero, and the blur grows
    /// in both directions from there.
    ///
    /// # Arguments
    ///
    /// * `z`: the object-side distance from the lens plane to the point, which must be greater
    ///   than zero. `f64::INFINITY` is permitted.
    ///
    /// returns: Option<f64>, or `None` if `z` is not greater than zero
    pub fn coc_diameter(&self, z: f64) -> Option<f64> {
        if z.is_nan() || z <= 0.0 {
            return None;
        }
        // Written through reciprocals so that both an infinite focus distance and an infinite
        // object distance stay finite. This is the usual closed form
        // `(f^2 / N) * |z - s_o| / (z * (s_o - f))` rearranged, and it remains valid for object
        // distances shorter than the focal length, where the light leaving the lens diverges.
        let scale = self.aperture_diameter() * self.image_distance();
        Some(scale * (1.0 / self.focus_distance - 1.0 / z).abs())
    }

    /// The near and far limits of the depth of field for an acceptable circle of confusion,
    /// which are the object distances at which the blur grows to the given tolerance. The far
    /// limit is `f64::INFINITY` when the lens is focused at or beyond its hyperfocal distance.
    ///
    /// # Arguments
    ///
    /// * `coc`: the largest circle of confusion diameter on the sensor plane that is still
    ///   considered acceptably sharp, which must be zero or greater
    ///
    /// returns: Option<(f64, f64)> of the near and far limits, or `None` if `coc` is negative
    ///   or NaN
    pub fn depth_of_field(&self, coc: f64) -> Option<(f64, f64)> {
        if coc.is_nan() || coc < 0.0 {
            return None;
        }
        // The reciprocal of the object distance moves linearly away from the reciprocal of the
        // focus distance as the blur grows, which makes both limits a single expression and puts
        // the transition to an infinite far limit on the sign of one subtraction.
        let spread = coc / (self.aperture_diameter() * self.image_distance());
        let near = 1.0 / (1.0 / self.focus_distance + spread);
        let far_recip = 1.0 / self.focus_distance - spread;
        let far = if far_recip > 0.0 {
            1.0 / far_recip
        } else {
            f64::INFINITY
        };
        Some((near, far))
    }

    /// The total depth of field for an acceptable circle of confusion, which is the distance
    /// between the near and far limits. This is `f64::INFINITY` when the far limit is infinite.
    ///
    /// # Arguments
    ///
    /// * `coc`: the largest circle of confusion diameter on the sensor plane that is still
    ///   considered acceptably sharp, which must be zero or greater
    ///
    /// returns: Option<f64>, or `None` if `coc` is negative or NaN
    pub fn total_depth_of_field(&self, coc: f64) -> Option<f64> {
        self.depth_of_field(coc).map(|(near, far)| far - near)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// A 50mm f/2.8 lens focused at 1m, in millimeters.
    fn lens_50mm() -> ThinLens {
        ThinLens::new(50.0, 2.8, 1000.0).unwrap()
    }

    #[test]
    fn conjugates_of_a_lens_focused_close() {
        let lens = lens_50mm();
        // 1 / (1/50 - 1/1000) = 1 / 0.019
        assert_relative_eq!(lens.image_distance(), 52.631578947368425, epsilon = 1e-9);
        // Also f / (s_o - f) = 50 / 950
        assert_relative_eq!(lens.magnification(), 50.0 / 950.0, epsilon = 1e-12);
        assert_relative_eq!(
            lens.effective_f_number(),
            2.8 * (1.0 + 50.0 / 950.0),
            epsilon = 1e-12
        );
        assert_relative_eq!(lens.aperture_diameter(), 50.0 / 2.8, epsilon = 1e-12);
        assert!(!lens.is_focused_at_infinity());
    }

    #[test]
    fn coc_matches_the_textbook_closed_form() {
        let lens = lens_50mm();
        let z = 1200.0;
        let coc = lens.coc_diameter(z).unwrap();

        // (f^2 / N) * |z - s_o| / (z * (s_o - f))
        let expected = (2500.0 / 2.8) * (z - 1000.0f64).abs() / (z * (1000.0 - 50.0));
        assert_relative_eq!(coc, expected, epsilon = 1e-12);
        assert_relative_eq!(coc, 0.15664160401002505, epsilon = 1e-9);
    }

    #[test]
    fn coc_is_zero_at_focus_and_bounded_at_infinity() {
        let lens = lens_50mm();
        assert_relative_eq!(lens.coc_diameter(1000.0).unwrap(), 0.0, epsilon = 1e-12);

        // At infinite object distance the blur converges to the aperture diameter times the
        // magnification.
        let expected = lens.aperture_diameter() * lens.magnification();
        assert_relative_eq!(
            lens.coc_diameter(f64::INFINITY).unwrap(),
            expected,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            lens.coc_diameter(f64::INFINITY).unwrap(),
            0.9398496240601504,
            epsilon = 1e-9
        );
    }

    #[test]
    fn coc_rejects_non_positive_distances() {
        let lens = lens_50mm();
        assert!(lens.coc_diameter(0.0).is_none());
        assert!(lens.coc_diameter(-1.0).is_none());
        assert!(lens.coc_diameter(f64::NAN).is_none());
    }

    #[test]
    fn depth_of_field_matches_the_textbook_closed_form() {
        let lens = lens_50mm();
        let coc = 0.00345;
        let (near, far) = lens.depth_of_field(coc).unwrap();

        // The textbook limits are s (H - f) / (H + s - 2f) and s (H - f) / (H - s)
        let h = lens.hyperfocal_distance(coc);
        let s = 1000.0;
        let f = 50.0;
        assert_relative_eq!(near, s * (h - f) / (h + s - 2.0 * f), epsilon = 1e-9);
        assert_relative_eq!(far, s * (h - f) / (h - s), epsilon = 1e-9);

        assert_relative_eq!(near, 996.3426254903501, epsilon = 1e-9);
        assert_relative_eq!(far, 1003.6843244180737, epsilon = 1e-9);
        assert_relative_eq!(
            lens.total_depth_of_field(coc).unwrap(),
            far - near,
            epsilon = 1e-12
        );
    }

    #[test]
    fn depth_of_field_limits_produce_the_stated_blur() {
        // The definition of the limits: a point at either one blurs to the tolerance.
        let lens = lens_50mm();
        let coc = 0.00345;
        let (near, far) = lens.depth_of_field(coc).unwrap();
        assert_relative_eq!(lens.coc_diameter(near).unwrap(), coc, epsilon = 1e-12);
        assert_relative_eq!(lens.coc_diameter(far).unwrap(), coc, epsilon = 1e-12);
    }

    #[test]
    fn focused_at_hyperfocal_reaches_infinity_and_half_of_it() {
        let coc = 0.00345;
        let h = lens_50mm().hyperfocal_distance(coc);
        assert_relative_eq!(h, 2500.0 / (2.8 * coc) + 50.0, epsilon = 1e-9);

        let lens = lens_50mm().focused_at(h).unwrap();
        let (near, far) = lens.depth_of_field(coc).unwrap();
        assert!(far.is_infinite());
        assert_relative_eq!(near, h / 2.0, epsilon = 1e-9);
    }

    #[test]
    fn zero_blur_tolerance_collapses_the_depth_of_field() {
        let lens = lens_50mm();
        let (near, far) = lens.depth_of_field(0.0).unwrap();
        assert_relative_eq!(near, 1000.0, epsilon = 1e-9);
        assert_relative_eq!(far, 1000.0, epsilon = 1e-9);
        assert!(lens.hyperfocal_distance(0.0).is_infinite());
        assert!(lens.depth_of_field(-1.0).is_none());
        assert!(lens.depth_of_field(f64::NAN).is_none());
    }

    #[test]
    fn airy_disk_uses_the_effective_f_number() {
        let lens = lens_50mm();
        let d = lens.airy_disk_diameter(WAVELENGTH_550NM_MM);
        assert_relative_eq!(
            d,
            2.44 * 0.00055 * lens.effective_f_number(),
            epsilon = 1e-15
        );
        assert_relative_eq!(d, 0.003955368421052632, epsilon = 1e-12);
    }

    #[test]
    fn macro_lens_at_one_to_one() {
        // At 1:1 the object and image distances are both twice the focal length, and the
        // effective f-number is twice the marked one.
        let lens = ThinLens::from_magnification(100.0, 4.0, 1.0).unwrap();
        assert_relative_eq!(lens.focus_distance, 200.0, epsilon = 1e-12);
        assert_relative_eq!(lens.image_distance(), 200.0, epsilon = 1e-9);
        assert_relative_eq!(lens.magnification(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(lens.effective_f_number(), 8.0, epsilon = 1e-12);
        assert_relative_eq!(lens.aperture_diameter(), 25.0, epsilon = 1e-12);
    }

    #[test]
    fn focused_at_infinity() {
        let lens = ThinLens::new(50.0, 2.8, f64::INFINITY).unwrap();
        assert!(lens.is_focused_at_infinity());
        assert_relative_eq!(lens.image_distance(), 50.0, epsilon = 1e-12);
        assert_relative_eq!(lens.magnification(), 0.0, epsilon = 1e-15);
        assert_relative_eq!(lens.effective_f_number(), 2.8, epsilon = 1e-12);

        // The blur of a near object reduces to f^2 / (N z)
        assert_relative_eq!(
            lens.coc_diameter(1000.0).unwrap(),
            2500.0 / (2.8 * 1000.0),
            epsilon = 1e-12
        );

        let (near, far) = lens.depth_of_field(0.00345).unwrap();
        assert!(far.is_infinite());
        assert_relative_eq!(near, 2500.0 / (2.8 * 0.00345), epsilon = 1e-9);

        // A lens focused at infinity holds focus from one focal length short of its own
        // hyperfocal distance
        assert_relative_eq!(
            near,
            lens.hyperfocal_distance(0.00345) - 50.0,
            epsilon = 1e-9
        );
    }

    #[test]
    fn field_width_solves_for_a_focal_length() {
        // A 50mm lens at 1m across an 8.4456mm sensor sees 160.4664mm of the object plane
        let lens = ThinLens::from_field_width(2.8, 1000.0, 8.4456, 160.4664).unwrap();
        assert_relative_eq!(lens.focal_length, 50.0, epsilon = 1e-9);
        assert_relative_eq!(lens.focus_distance, 1000.0, epsilon = 1e-12);
    }

    #[test]
    fn refocusing_leaves_the_other_parameters_alone() {
        let lens = lens_50mm().focused_at(500.0).unwrap();
        assert_relative_eq!(lens.focus_distance, 500.0, epsilon = 1e-12);
        assert_relative_eq!(lens.focal_length, 50.0, epsilon = 1e-12);
        assert_relative_eq!(lens.f_number, 2.8, epsilon = 1e-12);
        assert!(lens_50mm().focused_at(40.0).is_err());
    }

    #[test]
    fn invalid_lenses_rejected() {
        // Focus at or inside the focal length forms no real image
        assert!(ThinLens::new(50.0, 2.8, 40.0).is_err());
        assert!(ThinLens::new(50.0, 2.8, 50.0).is_err());
        assert!(ThinLens::new(50.0, 2.8, f64::NAN).is_err());

        assert!(ThinLens::new(0.0, 2.8, 1000.0).is_err());
        assert!(ThinLens::new(-50.0, 2.8, 1000.0).is_err());
        assert!(ThinLens::new(f64::INFINITY, 2.8, 1000.0).is_err());
        assert!(ThinLens::new(50.0, 0.0, 1000.0).is_err());
        assert!(ThinLens::new(50.0, -2.8, 1000.0).is_err());

        assert!(ThinLens::from_magnification(50.0, 2.8, 0.0).is_err());
        assert!(ThinLens::from_magnification(50.0, 2.8, -1.0).is_err());
        assert!(ThinLens::from_magnification(50.0, 2.8, f64::INFINITY).is_err());

        assert!(ThinLens::from_field_width(2.8, f64::INFINITY, 8.4456, 160.0).is_err());
        assert!(ThinLens::from_field_width(2.8, 1000.0, 0.0, 160.0).is_err());
        assert!(ThinLens::from_field_width(2.8, 1000.0, 8.4456, 0.0).is_err());
    }
}
