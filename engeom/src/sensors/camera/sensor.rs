//! This module models a digital imaging sensor as a rectangular grid of discrete
//! square pixels with a known physical size.

use crate::geom2::Aabb2;
use crate::{Point2, Result};

/// A digital imaging sensor, modeled as a rectangular grid of square pixels with a known
/// physical pitch.
///
/// The sensor does not model the lens in front of it or its position in the world. It converts
/// between the discrete pixel grid and physical lengths on the sensor plane. The caller chooses
/// the length unit for `pitch`, and every derived length uses that same unit. The unit is usually
/// millimeters; in that unit, a typical machine-vision pixel pitch is approximately 0.0035.
///
/// Pixel coordinates follow the same convention as `PinholeCamera`: they are continuous rather
/// than integer. The pixel with integer index `i` covers the interval `[i, i + 1)`, so the
/// center of the image lies at `(width_px / 2, height_px / 2)`.
///
/// The principal point is assumed to be at the center of the sensor, and pixels are assumed to
/// be square. A lens-distortion model should represent sensors for which either assumption is
/// false.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sensor {
    /// The number of pixel columns in the sensor.
    pub width_px: u32,

    /// The number of pixel rows in the sensor.
    pub height_px: u32,

    /// The physical center-to-center distance between adjacent pixels.
    pub pitch: f64,
}

impl Sensor {
    /// Create a new sensor from its pixel dimensions and physical pixel pitch.
    ///
    /// # Arguments
    ///
    /// * `width_px`: the number of pixel columns, which must not be zero
    /// * `height_px`: the number of pixel rows, which must not be zero
    /// * `pitch`: the physical center-to-center distance between adjacent pixels, which must be
    ///   a finite number greater than zero
    ///
    /// returns: Result<Sensor>
    pub fn new(width_px: u32, height_px: u32, pitch: f64) -> Result<Self> {
        if width_px == 0 {
            return Err("Sensor width_px must be greater than zero".into());
        }
        if height_px == 0 {
            return Err("Sensor height_px must be greater than zero".into());
        }
        if !pitch.is_finite() || pitch <= 0.0 {
            return Err("Sensor pitch must be a finite number greater than zero".into());
        }
        Ok(Self {
            width_px,
            height_px,
            pitch,
        })
    }

    /// Create a copy of this sensor with its pixel dimensions scaled by a factor. The pitch is
    /// inversely scaled to preserve the sensor's physical size. A camera built with the copy
    /// therefore has the same field of view with finer or coarser sampling. This property makes
    /// a reduced-resolution sensor useful for previewing a full-resolution image.
    ///
    /// The scaled pixel dimensions are rounded to the nearest whole pixel and never fall below
    /// one, so the physical size is preserved only to within that rounding.
    ///
    /// # Arguments
    ///
    /// * `factor`: the factor to scale the pixel dimensions by, which must be a finite number
    ///   greater than zero
    ///
    /// returns: Result<Sensor>
    pub fn rescaled(&self, factor: f64) -> Result<Self> {
        if !factor.is_finite() || factor <= 0.0 {
            return Err("Sensor rescale factor must be a finite number greater than zero".into());
        }
        let width_px = ((self.width_px as f64 * factor).round() as u32).max(1);
        let height_px = ((self.height_px as f64 * factor).round() as u32).max(1);
        Self::new(width_px, height_px, self.pitch / factor)
    }

    /// The physical width of the active sensor area.
    pub fn width(&self) -> f64 {
        self.width_px as f64 * self.pitch
    }

    /// The physical height of the active sensor area.
    pub fn height(&self) -> f64 {
        self.height_px as f64 * self.pitch
    }

    /// The physical diagonal of the active sensor area.
    pub fn diagonal(&self) -> f64 {
        self.width().hypot(self.height())
    }

    /// The ratio of the sensor's width to its height. Because the pixels are square this is the
    /// same whether it is computed in pixels or in physical units.
    pub fn aspect_ratio(&self) -> f64 {
        self.width_px as f64 / self.height_px as f64
    }

    /// The total number of pixels in the sensor.
    pub fn pixel_count(&self) -> u64 {
        self.width_px as u64 * self.height_px as u64
    }

    /// The center of the image in continuous pixel coordinates, which is where the principal
    /// point is assumed to lie.
    pub fn center_px(&self) -> Point2 {
        Point2::new(self.width_px as f64 / 2.0, self.height_px as f64 / 2.0)
    }

    /// The region of continuous pixel coordinates covered by the sensor, represented as a closed box running
    /// from the outside corner of the first pixel to the outside corner of the last.
    ///
    /// Because pixel `i` covers `[i, i + 1)`, the image spans `[0, width_px]` horizontally and
    /// `[0, height_px]` vertically. This is the region to clip or intersect against. It is closed
    /// on all four sides, so it includes its far edges, while `contains_pixel` excludes them.
    ///
    /// returns: Aabb2
    pub fn image_aabb(&self) -> Aabb2 {
        Aabb2::new(
            Point2::origin(),
            Point2::new(self.width_px as f64, self.height_px as f64),
        )
    }

    /// Whether a point in continuous pixel coordinates is recorded by a pixel of this sensor.
    ///
    /// The test is half-open on the far edges, matching the convention that pixel `i` covers
    /// `[i, i + 1)`: a point is contained when the floor of each coordinate is a valid column and
    /// row. A point at `x = width_px` is therefore not contained, even though it lies on the
    /// boundary of `image_aabb`. A point with a NaN coordinate is not contained.
    ///
    /// # Arguments
    ///
    /// * `pixel`: a point in continuous pixel coordinates
    ///
    /// returns: bool
    pub fn contains_pixel(&self, pixel: &Point2) -> bool {
        pixel.x >= 0.0
            && pixel.y >= 0.0
            && pixel.x < self.width_px as f64
            && pixel.y < self.height_px as f64
    }

    /// Convert a point in continuous pixel coordinates to a physical position on the sensor
    /// plane, measured from the center of the sensor. The physical axes match the pixel axes,
    /// so +x is to the right in the image and +y is down.
    ///
    /// # Arguments
    ///
    /// * `pixel`: a point in continuous pixel coordinates
    ///
    /// returns: Point2
    pub fn pixel_to_plane(&self, pixel: &Point2) -> Point2 {
        let c = self.center_px();
        Point2::new((pixel.x - c.x) * self.pitch, (pixel.y - c.y) * self.pitch)
    }

    /// Convert a physical position on the sensor plane, measured from the center of the sensor,
    /// to a point in continuous pixel coordinates. This is the inverse of `pixel_to_plane`.
    ///
    /// # Arguments
    ///
    /// * `plane`: a physical position on the sensor plane relative to the sensor center
    ///
    /// returns: Point2
    pub fn plane_to_pixel(&self, plane: &Point2) -> Point2 {
        let c = self.center_px();
        Point2::new(plane.x / self.pitch + c.x, plane.y / self.pitch + c.y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// A common machine vision sensor: 2448 x 2048 pixels on a 3.45 micron pitch, in mm.
    fn test_sensor() -> Sensor {
        Sensor::new(2448, 2048, 0.00345).unwrap()
    }

    #[test]
    fn physical_dimensions() {
        let s = test_sensor();
        assert_relative_eq!(s.width(), 8.4456, epsilon = 1e-10);
        assert_relative_eq!(s.height(), 7.0656, epsilon = 1e-10);
        assert_relative_eq!(s.diagonal(), 11.011397_f64, epsilon = 1e-6);
        assert_relative_eq!(s.aspect_ratio(), 2448.0 / 2048.0, epsilon = 1e-12);
        assert_eq!(s.pixel_count(), 2448 * 2048);
    }

    #[test]
    fn center_is_half_the_pixel_count() {
        let s = test_sensor();
        assert_relative_eq!(s.center_px().x, 1224.0, epsilon = 1e-12);
        assert_relative_eq!(s.center_px().y, 1024.0, epsilon = 1e-12);
    }

    #[test]
    fn plane_coords_measured_from_center() {
        let s = test_sensor();

        // The center pixel is the origin of the sensor plane
        let c = s.pixel_to_plane(&Point2::new(1224.0, 1024.0));
        assert_relative_eq!(c.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(c.y, 0.0, epsilon = 1e-12);

        // The upper left corner is half the physical size away in -x and -y
        let corner = s.pixel_to_plane(&Point2::new(0.0, 0.0));
        assert_relative_eq!(corner.x, -4.2228, epsilon = 1e-10);
        assert_relative_eq!(corner.y, -3.5328, epsilon = 1e-10);
    }

    #[test]
    fn plane_to_pixel_round_trip() {
        let s = test_sensor();
        let start = Point2::new(137.5, 1900.25);
        let back = s.plane_to_pixel(&s.pixel_to_plane(&start));
        assert_relative_eq!(back.x, start.x, epsilon = 1e-9);
        assert_relative_eq!(back.y, start.y, epsilon = 1e-9);
    }

    #[test]
    fn image_aabb_spans_the_whole_pixel_grid() {
        let s = Sensor::new(200, 100, 0.01).unwrap();
        let box2 = s.image_aabb();
        assert_relative_eq!(box2.mins.x, 0.0, epsilon = 1e-12);
        assert_relative_eq!(box2.mins.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(box2.maxs.x, 200.0, epsilon = 1e-12);
        assert_relative_eq!(box2.maxs.y, 100.0, epsilon = 1e-12);
    }

    #[test]
    fn contains_pixel_is_half_open_on_the_far_edges() {
        let s = Sensor::new(200, 100, 0.01).unwrap();
        assert!(s.contains_pixel(&Point2::new(0.0, 0.0)));
        assert!(s.contains_pixel(&Point2::new(199.999, 99.999)));
        assert!(!s.contains_pixel(&Point2::new(200.0, 50.0)));
        assert!(!s.contains_pixel(&Point2::new(100.0, 100.0)));
        assert!(!s.contains_pixel(&Point2::new(-1e-9, 50.0)));
        assert!(!s.contains_pixel(&Point2::new(f64::NAN, 50.0)));
        assert!(!s.contains_pixel(&Point2::new(50.0, f64::NAN)));
    }

    #[test]
    fn contains_pixel_agrees_with_flooring_to_a_valid_index() {
        let s = Sensor::new(7, 5, 1.0).unwrap();
        for i in -2..12 {
            for j in -2..9 {
                let p = Point2::new(i as f64 * 0.75, j as f64 * 0.75);
                let valid = (0..7).contains(&(p.x.floor() as i64))
                    && (0..5).contains(&(p.y.floor() as i64));
                assert_eq!(s.contains_pixel(&p), valid, "at {p:?}");
            }
        }
    }

    #[test]
    fn the_aabb_is_closed_where_contains_pixel_is_not() {
        // The two disagree on the far edges by design: the box is a region to clip against and
        // includes its boundary, while the predicate asks whether a pixel records the point.
        let s = Sensor::new(200, 100, 0.01).unwrap();
        let corner = Point2::new(200.0, 100.0);
        assert!(s.image_aabb().contains_local_point(&corner));
        assert!(!s.contains_pixel(&corner));
    }

    #[test]
    fn a_rescaled_sensor_reports_its_own_bounds() {
        let s = Sensor::new(200, 100, 0.01).unwrap().rescaled(0.5).unwrap();
        assert_relative_eq!(s.image_aabb().maxs.x, 100.0, epsilon = 1e-12);
        assert_relative_eq!(s.image_aabb().maxs.y, 50.0, epsilon = 1e-12);
        assert!(s.contains_pixel(&Point2::new(99.5, 49.5)));
        assert!(!s.contains_pixel(&Point2::new(100.0, 25.0)));
    }

    #[test]
    fn invalid_sensors_rejected() {
        assert!(Sensor::new(0, 10, 1.0).is_err());
        assert!(Sensor::new(10, 0, 1.0).is_err());
        assert!(Sensor::new(10, 10, 0.0).is_err());
        assert!(Sensor::new(10, 10, -1.0).is_err());
        assert!(Sensor::new(10, 10, f64::NAN).is_err());
        assert!(Sensor::new(10, 10, f64::INFINITY).is_err());
    }
}
