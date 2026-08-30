//! This module contains features for taking measurements on meshes

use crate::common::DistMode;
use crate::metrology::Distance3;
use crate::{Mesh3, Point3, Result, SurfacePoint3, UnitVec3};
use rayon::prelude::*;

/// The most samples either direction of a surface deviation measurement will take before refusing.
///
/// This exists because the spacing is easy to get catastrophically wrong. Asking for a spacing at
/// the scale of a tight tolerance on a coarse mesh requests a barycentric grid thousands of points
/// wide on every face, which does not fail gracefully: it exhausts memory and takes the process
/// with it. Refusing with a number the caller can act on is better than dying.
const MAX_DEVIATION_SAMPLES: f64 = 5.0e7;

/// A two-sided deviation measurement between the surfaces of two meshes.
///
/// The two one-sided values answer genuinely different questions and neither implies the other. If
/// the test surface drops a feature the reference had, `reference_to_test` grows while
/// `test_to_reference` may not move at all; if the test surface invents geometry the reference
/// never had, the reverse happens. Read both.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceDeviation {
    /// The furthest any sampled point of the reference surface sits from the test surface.
    ///
    /// This is the "is the reference still represented" direction. It grows when the test surface
    /// loses something the reference had.
    pub reference_to_test: f64,

    /// The furthest any sampled point of the test surface sits from the reference surface.
    ///
    /// This is the "does the test surface claim anything the reference did not" direction. It grows
    /// when the test surface bridges a hole, spans a concavity, or bows away between its vertices.
    pub test_to_reference: f64,

    /// How many samples were taken on the reference surface.
    pub reference_samples: usize,

    /// How many samples were taken on the test surface.
    pub test_samples: usize,
}

impl SurfaceDeviation {
    /// The symmetric Hausdorff distance, which is the larger of the two one-sided values.
    pub fn hausdorff(&self) -> f64 {
        self.reference_to_test.max(self.test_to_reference)
    }
}

impl std::fmt::Display for SurfaceDeviation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ref->test {:.6} ({} samples), test->ref {:.6} ({} samples), hausdorff {:.6}",
            self.reference_to_test,
            self.reference_samples,
            self.test_to_reference,
            self.test_samples,
            self.hausdorff()
        )
    }
}

/// A sample spacing derived from a mesh's own face sizes, giving a bounded number of samples per
/// face regardless of how coarse or fine the mesh is.
///
/// This is the edge length of an equilateral triangle of the mean face area, quartered, so a
/// typical face receives on the order of ten samples.
pub(crate) fn derived_spacing(mesh: &Mesh3) -> Result<f64> {
    let areas = mesh.compute_face_areas()?;
    if areas.is_empty() {
        return Err("Cannot derive a sample spacing for a mesh with no faces".into());
    }

    let mean = areas.iter().sum::<f64>() / areas.len() as f64;
    if mean <= 0.0 {
        return Err("Cannot derive a sample spacing for a mesh of zero area".into());
    }

    Ok((mean * 4.0 / 3.0_f64.sqrt()).sqrt() / 4.0)
}

/// Refuse a spacing which would demand an unreasonable number of samples on `from`.
fn check_sample_budget(from: &Mesh3, spacing: f64, label: &str) -> Result<()> {
    let areas = from.compute_face_areas()?;
    let total: f64 = areas.iter().sum();
    let estimate = total / (spacing * spacing);
    if estimate > MAX_DEVIATION_SAMPLES {
        return Err(format!(
            "A spacing of {} would take roughly {:.0} samples on the {} surface, over the limit of \
             {:.0}. Use a larger spacing.",
            spacing, estimate, label, MAX_DEVIATION_SAMPLES
        )
        .into());
    }

    Ok(())
}

/// The worst distance from any dense sample of `from` to the surface of `to`.
fn one_sided(from: &Mesh3, to: &Mesh3, spacing: f64, label: &str) -> Result<(f64, usize)> {
    check_sample_budget(from, spacing, label)?;

    let samples = from.sample_dense(spacing, None)?;
    let worst = samples
        .points()
        .par_iter()
        .map(|p| to.distance_closest_to(p))
        .reduce(|| 0.0f64, f64::max);

    Ok((worst, samples.point_count()))
}

/// Every distance from a dense sample of `from` to the surface of `to`, rather than only the worst.
///
/// This exists for characterizing a method which has no bound, where the shape of the distribution
/// is the whole question and a maximum says nothing: a decimator that is usually excellent and
/// occasionally terrible and one that is uniformly mediocre report the same maximum.
///
/// [`one_sided`] deliberately does not do this. It reduces with `max` so that nothing is retained,
/// which is what lets it run at spacings producing tens of millions of samples. Holding those would
/// be hundreds of megabytes, so this is for the derived spacing and diagnostics, not for a sweep.
///
/// Test-gated because the only caller is the shape-difference harness. Promote it if a shipped feature
/// ever wants a distribution rather than a worst case.
#[cfg(test)]
pub(crate) fn sample_distances(
    from: &Mesh3,
    to: &Mesh3,
    spacing: f64,
    label: &str,
) -> Result<Vec<f64>> {
    check_sample_budget(from, spacing, label)?;

    Ok(from
        .sample_dense(spacing, None)?
        .points()
        .par_iter()
        .map(|p| to.distance_closest_to(p))
        .collect())
}

impl Mesh3 {
    /// Measure the two-sided deviation between this mesh's surface and another's.
    ///
    /// Both surfaces are sampled densely and every sample is measured against the other mesh, so
    /// this is an approximation of the symmetric Hausdorff distance which converges from below as
    /// the spacing tightens. It is deliberately independent of whatever produced the test mesh:
    /// nothing here consults a decimator's bookkeeping, an alignment's residuals, or any other
    /// record of what was supposed to have happened.
    ///
    /// Sampling both surfaces is the point. Measuring only the vertices of one mesh against the
    /// other is a much cheaper test and a much weaker one, because the maximum over a triangle's
    /// three corners does not bound the triangle: a surface can agree exactly at every vertex and
    /// still bow well past tolerance in between. Any check that samples the same sites an algorithm
    /// enforces will agree with that algorithm whether or not it is right.
    ///
    /// This is a measurement rather than a bound. It is the reference a bound gets checked against,
    /// and it is slow enough to belong in tests and diagnostics rather than in a hot path.
    ///
    /// # Arguments
    ///
    /// * `test`: the mesh to compare against, treated as the test surface
    /// * `max_spacing`: the sample spacing on both surfaces. `None` derives a spacing from each
    ///   mesh's own mean face area, which keeps the sample count proportional to the face count
    ///   rather than to the ratio between the face size and some unrelated tolerance.
    ///
    /// returns: `Result<SurfaceDeviation>`, failing if either mesh has no faces or if the requested
    /// spacing would demand an unreasonable number of samples
    pub fn measure_surface_deviation(
        &self,
        test: &Mesh3,
        max_spacing: Option<f64>,
    ) -> Result<SurfaceDeviation> {
        if let Some(s) = max_spacing
            && (s.is_nan() || s <= 0.0)
        {
            return Err(format!("A sample spacing must be positive, got {}", s).into());
        }

        let reference_spacing = match max_spacing {
            Some(s) => s,
            None => derived_spacing(self)?,
        };
        let test_spacing = match max_spacing {
            Some(s) => s,
            None => derived_spacing(test)?,
        };

        let (reference_to_test, reference_samples) =
            one_sided(self, test, reference_spacing, "reference")?;
        let (test_to_reference, test_samples) = one_sided(test, self, test_spacing, "test")?;

        Ok(SurfaceDeviation {
            reference_to_test,
            test_to_reference,
            reference_samples,
            test_samples,
        })
    }
}

impl Mesh3 {
    /// Compute the deviation of a point from this mesh (this mesh is considered the reference) and
    /// return it as a Length Measurement object.
    ///
    /// The deviation is the distance from the point to its closest projection onto the mesh using
    /// the specified distance mode.  The direction of the measurement is the direction between the
    /// point and the projection, flipped into the positive half-space of the mesh surface at the
    /// projection point.
    ///
    /// If the distance is less than a very small floating point epsilon, the direction will be
    /// taken directly from the mesh surface normal.
    ///
    /// The first point `.a` of the measurement is the reference point, and the second point `.b`
    /// is the test point.
    ///
    /// # Arguments
    ///
    /// * `point`: the test point to measure the deviation from
    /// * `dist_mode`: whether to use the point-to-point distance or the scalar projection distance
    ///   when computing the deviation. This will have an effect near the edges of the mesh, in
    ///   which the `ToPlane` mode will not penalize a point for being off the mesh surface.
    ///
    /// returns: Length<3>
    ///
    /// # Examples
    ///
    /// ```
    ///
    /// ```
    pub fn measure_point_deviation(&self, point: &Point3, dist_mode: DistMode) -> Distance3 {
        let closest = self.surface_closest_to(point).sp;

        // In both cases, the measurement point `b` will remain the test point and `a` will be the
        // where the reference point, what will change is the direction of the measurement

        let d = deviation_direction(&closest, point, dist_mode);

        Distance3::new(closest.point, *point, Some(d))
    }

    /// Compute the signed deviation of many points from this mesh at once, returning the scalar
    /// values rather than full measurements.
    ///
    /// Each value is what `measure_point_deviation` would report for the same point and mode, and
    /// the two share their sign convention. This exists because the scalar is what a bulk
    /// comparison actually wants, and building a `Distance3` per point to read one number off it
    /// costs two points and a direction each time.
    ///
    /// # Arguments
    ///
    /// * `points`: the test points to measure
    /// * `dist_mode`: whether to use the point-to-point distance or the scalar projection distance,
    ///   which differ near the edges of the mesh
    ///
    /// returns: `Vec<f64>`, one signed deviation per test point, positive on the outside
    pub fn measure_deviations(&self, points: &[Point3], dist_mode: DistMode) -> Vec<f64> {
        points
            .iter()
            .map(|point| {
                let closest = self.surface_closest_to(point).sp;
                let d = deviation_direction(&closest, point, dist_mode);
                d.dot(&(point - closest.point))
            })
            .collect()
    }
}

/// The direction a deviation measurement is taken along, given a test point and its projection.
///
/// A test point sitting on the surface has no meaningful direction to its own projection, so the
/// surface normal is used rather than normalizing what would be floating point residue.
fn deviation_direction(closest: &SurfacePoint3, point: &Point3, dist_mode: DistMode) -> UnitVec3 {
    match dist_mode {
        DistMode::ToPoint => {
            let v = point - closest.point;
            if v.norm() < 1e-6 {
                closest.normal
            } else if closest.normal.dot(&v) > 0.0 {
                UnitVec3::new_normalize(v)
            } else {
                -UnitVec3::new_normalize(v)
            }
        }
        DistMode::ToPlane => closest.normal,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metrology::Measurement;
    use approx::assert_relative_eq;

    /// The bulk and single-point paths share a direction helper, so they must agree exactly rather
    /// than merely closely.
    #[test]
    fn measure_deviations_matches_the_single_point_measurement() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, true);
        let points = vec![
            Point3::new(0.0, 0.0, 3.0),   // clearly outside
            Point3::new(0.0, 0.0, -3.0),  // clearly outside, other side
            Point3::new(0.0, 0.0, 1.0),   // exactly on the surface
            Point3::new(5.0, 5.0, 5.0),   // off the corner, where the modes differ
            Point3::new(0.1, -0.2, 0.05), // inside
        ];

        for mode in [DistMode::ToPoint, DistMode::ToPlane] {
            let bulk = mesh.measure_deviations(&points, mode);
            assert_eq!(bulk.len(), points.len());

            for (i, point) in points.iter().enumerate() {
                let single = mesh.measure_point_deviation(point, mode).value();
                assert_relative_eq!(bulk[i], single, epsilon = 1.0e-15);
            }
        }
    }

    #[test]
    fn measure_deviations_signs_the_outside_positive() {
        let mesh = Mesh3::create_box(2.0, 2.0, 2.0, false);
        let values = mesh.measure_deviations(
            &[Point3::new(0.0, 0.0, 2.0), Point3::new(0.0, 0.0, 0.5)],
            DistMode::ToPoint,
        );

        assert_relative_eq!(values[0], 1.0, epsilon = 1.0e-12);
        assert!(values[1] < 0.0);
    }

    #[test]
    fn measure_deviations_of_nothing_is_nothing() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true);
        assert!(mesh.measure_deviations(&[], DistMode::ToPoint).is_empty());
    }

    #[test]
    fn a_mesh_does_not_deviate_from_itself() {
        let mesh = Mesh3::create_sphere(10.0, 0.06).unwrap();
        let d = mesh.measure_surface_deviation(&mesh, None).unwrap();

        assert!(d.hausdorff() < 1.0e-9, "got {}", d);
        assert!(d.reference_samples > 0);
        assert!(d.test_samples > 0);
    }

    /// Two parallel squares of the same footprint, where the answer is known in closed form. This
    /// pins the measurement to a number rather than to a sanity check.
    fn square_at(z: f64) -> Mesh3 {
        Mesh3::new(
            vec![
                Point3::new(-1.0, -1.0, z),
                Point3::new(1.0, -1.0, z),
                Point3::new(1.0, 1.0, z),
                Point3::new(-1.0, 1.0, z),
            ],
            vec![[0u32, 1, 2], [0, 2, 3]],
            false,
        )
    }

    #[test]
    fn a_uniform_offset_is_measured_in_both_directions() {
        let d = square_at(0.0)
            .measure_surface_deviation(&square_at(0.25), None)
            .unwrap();

        // The footprints coincide, so every point of either surface projects straight onto the
        // other and the separation is exactly the offset in both directions.
        assert_relative_eq!(d.reference_to_test, 0.25, epsilon = 1.0e-9);
        assert_relative_eq!(d.test_to_reference, 0.25, epsilon = 1.0e-9);
        assert_relative_eq!(d.hausdorff(), 0.25, epsilon = 1.0e-9);
    }

    /// The two directions are not interchangeable, and a measurement which quietly reported one of
    /// them would pass every symmetric test above.
    #[test]
    fn the_two_directions_are_reported_separately() {
        // A small square against a large one sharing a plane. Every point of the small surface sits
        // on the large one, but most of the large surface is far from the small one.
        let small = square_at(0.0);
        let large = Mesh3::new(
            vec![
                Point3::new(-5.0, -5.0, 0.0),
                Point3::new(5.0, -5.0, 0.0),
                Point3::new(5.0, 5.0, 0.0),
                Point3::new(-5.0, 5.0, 0.0),
            ],
            vec![[0u32, 1, 2], [0, 2, 3]],
            false,
        );

        let d = small.measure_surface_deviation(&large, None).unwrap();
        assert!(d.reference_to_test < 1.0e-9, "got {}", d);
        assert!(d.test_to_reference > 3.0, "got {}", d);
        assert_relative_eq!(d.hausdorff(), d.test_to_reference, epsilon = 1.0e-12);
    }

    /// The whole reason this measurement samples surfaces rather than vertices.
    ///
    /// Two surfaces which share every vertex can still disagree between them. A vertex-only check
    /// reports zero here; this one must not.
    #[test]
    fn surfaces_agreeing_at_every_vertex_can_still_deviate() {
        // The two triangulations of a saddle-shaped quad. They share all four corners exactly, so
        // no vertex-based comparison can tell them apart, but one bulges up along its diagonal and
        // the other bulges down along the opposite one.
        let corners = vec![
            Point3::new(-1.0, -1.0, 1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(-1.0, 1.0, -1.0),
        ];

        let up = Mesh3::new(corners.clone(), vec![[0u32, 1, 2], [0, 2, 3]], false);
        let down = Mesh3::new(corners, vec![[0u32, 1, 3], [1, 2, 3]], false);

        // The vertex-only view: every vertex of each mesh sits on the other.
        for p in up.points() {
            assert!(down.distance_closest_to(p) < 1.0e-9);
        }
        for p in down.points() {
            assert!(up.distance_closest_to(p) < 1.0e-9);
        }

        // The surface view disagrees, which is the point.
        let d = up.measure_surface_deviation(&down, None).unwrap();
        assert!(
            d.hausdorff() > 0.5,
            "a vertex-only check would report zero here, got {}",
            d
        );
    }

    #[test]
    fn an_absurd_spacing_is_refused_rather_than_attempted() {
        let mesh = Mesh3::create_sphere(100.0, 0.1).unwrap();
        let err = mesh
            .measure_surface_deviation(&mesh, Some(1.0e-4))
            .unwrap_err();

        assert!(
            err.to_string().contains("larger spacing"),
            "unhelpful message: {}",
            err
        );
    }

    #[test]
    fn a_nonsense_spacing_is_rejected() {
        let mesh = Mesh3::create_box(1.0, 1.0, 1.0, true);
        assert!(mesh.measure_surface_deviation(&mesh, Some(0.0)).is_err());
        assert!(mesh.measure_surface_deviation(&mesh, Some(-1.0)).is_err());
    }
}
