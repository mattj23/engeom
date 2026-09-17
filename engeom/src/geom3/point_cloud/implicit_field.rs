//! Turning an oriented point cloud into a signed distance field, and sampling that field onto a
//! narrow band of voxels.
//!
//! # What the field is
//!
//! Each cloud point and its normal define a tangent plane. The signed distance from a query position
//! to one tangent plane approximates the surface distance near that point. [`ImlsField3`] averages
//! the estimates from all points within a radius and weights them by proximity. This is the implicit
//! moving least squares surface of Kolluri, based on the earlier moving least squares point-set
//! surfaces of Levin and Alexa.
//!
//! Hoppe's 1992 method instead uses the plane of the single nearest point. That approach is less
//! expensive, but its field is discontinuous where the nearest point changes, creating small steps
//! at point boundaries. Each value also depends on one measurement and therefore retains its noise.
//! Averaging produces a smooth field and reduces individual measurement noise across a
//! neighborhood.
//!
//! # The bias this method has
//!
//! Averaging tangent planes over a neighborhood pulls the zero level set away from a curved
//! surface. Every tangent plane around a convex patch lies outside the surface it touches, so their
//! average reads slightly negative at a point on the surface. The zero crossing lies outside the
//! true surface by approximately
//!
//! ```text
//!     bias ~ 0.46 * h^2 / R
//! ```
//!
//! where `h` is the weight scale and `R` is the radius of curvature. The bias is second order and
//! systematic, and it vanishes on a flat surface. Most measured features of machined parts are
//! flat. At a weight scale of one-tenth the radius of curvature, the bias is approximately five
//! percent of a voxel, the same order as the existing error from linear interpolation along a voxel
//! edge.
//!
//! The test module measures this coefficient across several weight scales using points on an
//! analytic sphere and obtains values from 0.45 to 0.47. A tessellated sphere is unsuitable for
//! this measurement because its chordal facets lie inside the approximated surface and introduce
//! an additional offset.
//!
//! A smaller weight scale reduces bias but also reduces noise averaging. Eliminating this tradeoff
//! requires fitting a curved surface, such as a sphere or quadric, instead of a plane. Such a method
//! can implement the same trait without changing downstream processing.
//!
//! # Sharp edges
//!
//! Averaging across a sharp edge mixes planes with different directions and rounds the edge over the
//! neighborhood width. Robust weighting can recover the edge by discounting planes that disagree
//! with the local consensus. This is the RIMLS method of Oztireli and others, and it can be provided
//! as another implementation of the same trait.

use crate::common::kd_tree::KdTreeSearch;
use crate::raster3::{BLOCK_VOXELS, Block3, BlockKey3, SdfVoxel, SparseGrid3, VoxelValue};
use crate::{Point3, Result, UnitVec3};
use rayon::prelude::*;

/// One evaluation of an implicit field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldSample3 {
    /// The signed distance to the surface, positive on the side the normals point toward.
    pub value: f64,

    /// The evidence supporting the value, in field-specific units. It is meaningful only relative
    /// to other samples of the same field and permits several fields to be averaged later.
    pub weight: f64,
}

/// A scalar field over space whose zero level set is a surface.
///
/// This trait separates field construction from surface extraction. Downstream code uses the trait,
/// so another field implementation or a field assembled from several scans can use the same voxel
/// grid and extraction path.
pub trait ImplicitField3: Sync {
    /// The field at a position, or `None` where no data supports a value.
    ///
    /// Returning `None` prevents reconstruction in unsupported regions. A voxel whose evaluation
    /// returns `None` remains unwritten, and extraction omits every cell that touches it.
    fn evaluate(&self, p: &Point3) -> Option<FieldSample3>;

    /// The distance from the data beyond which [`ImplicitField3::evaluate`] always returns `None`.
    ///
    /// This value determines which voxels are evaluated. Overestimating it causes extra
    /// evaluations; underestimating it creates holes. An implementation should round up when the
    /// support limit is uncertain.
    fn support_radius(&self) -> f64;
}

/// A signed distance field built by averaging the tangent planes of an oriented point cloud.
///
/// See the module documentation for what the method is, the bias it has, and what it does to sharp
/// edges. Build one with [`ImlsField3::try_new`].
pub struct ImlsField3<'a, T> {
    points: &'a [Point3],
    normals: &'a [UnitVec3],
    confidence: Option<&'a [f64]>,
    tree: &'a T,
    radius: f64,
    sigma: f64,
}

impl<'a, T: KdTreeSearch<3> + Sync> ImlsField3<'a, T> {
    /// Build a field over an oriented point cloud.
    ///
    /// The normals must be oriented. A plane fit gives an axis without a direction, and arbitrary
    /// axis signs cause the field's distance sign to vary between points. Use
    /// [`super::orientation`] to select consistent directions.
    ///
    /// # Arguments
    ///
    /// * `points`: the cloud positions, which `tree` must have been built over
    /// * `normals`: the oriented normal at each point, one per point
    /// * `tree`: a k-d tree over `points`
    /// * `radius`: the maximum distance from a position at which points are gathered. The field
    ///   reports no value beyond this support radius. Must be finite and positive.
    /// * `sigma`: the weight scale, or the distance at which a neighbor's contribution falls to
    ///   `1/e`. Must be finite and positive. It should be at least the point spacing so that most
    ///   positions have enough contributing neighbors for averaging. A radius of twice `sigma`
    ///   captures all but approximately two percent of the weight.
    ///
    /// returns: `Result<ImlsField3<T>>`
    pub fn try_new(
        points: &'a [Point3],
        normals: &'a [UnitVec3],
        tree: &'a T,
        radius: f64,
        sigma: f64,
    ) -> Result<Self> {
        if points.len() != normals.len() {
            return Err(format!(
                "there are {} points but {} normals",
                points.len(),
                normals.len()
            )
            .into());
        }

        if !radius.is_finite() || radius <= 0.0 {
            return Err(format!("Radius must be finite and positive, got {radius}").into());
        }

        if !sigma.is_finite() || sigma <= 0.0 {
            return Err(format!("Sigma must be finite and positive, got {sigma}").into());
        }

        Ok(Self {
            points,
            normals,
            confidence: None,
            tree,
            radius,
            sigma,
        })
    }

    /// Weight each point by a per-point confidence as well as by distance.
    ///
    /// The confidence from [`super::NormalEstimates`] can reduce the contribution of a normal fitted
    /// from a sparse or nonplanar neighborhood relative to one fitted from a clearly planar
    /// neighborhood. A confidence of zero removes a point from the field, providing a weighted
    /// alternative to filtering the cloud first.
    ///
    /// # Arguments
    ///
    /// * `confidence`: one non-negative value per point
    ///
    /// returns: `Result<ImlsField3<T>>`
    pub fn with_confidence(mut self, confidence: &'a [f64]) -> Result<Self> {
        if confidence.len() != self.points.len() {
            return Err(format!(
                "there are {} points but {} confidence values",
                self.points.len(),
                confidence.len()
            )
            .into());
        }

        if confidence.iter().any(|c| !c.is_finite() || *c < 0.0) {
            return Err("Confidence values must be finite and non-negative".into());
        }

        self.confidence = Some(confidence);
        Ok(self)
    }

    /// The radius within which points are gathered.
    pub fn radius(&self) -> f64 {
        self.radius
    }

    /// The weight scale.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }
}

impl<T: KdTreeSearch<3> + Sync> ImplicitField3 for ImlsField3<'_, T> {
    fn evaluate(&self, p: &Point3) -> Option<FieldSample3> {
        let found = self.tree.within(p, self.radius);
        if found.is_empty() {
            return None;
        }

        let inv_sigma_sq = 1.0 / (self.sigma * self.sigma);
        let mut numerator = 0.0;
        let mut denominator = 0.0;

        for (i, distance) in found {
            let mut weight = (-distance * distance * inv_sigma_sq).exp();
            if let Some(confidence) = self.confidence {
                weight *= confidence[i];
            }

            if weight <= 0.0 {
                continue;
            }

            numerator += weight * self.normals[i].dot(&(p - self.points[i]));
            denominator += weight;
        }

        // Zero confidence can remove every neighbor, and all weights can underflow. Both cases
        // leave this position without supporting data.
        if denominator <= 0.0 {
            return None;
        }

        Some(FieldSample3 {
            value: numerator / denominator,
            weight: denominator,
        })
    }

    fn support_radius(&self) -> f64 {
        self.radius
    }
}

/// Sample a field onto a narrow band of voxels around a set of positions.
///
/// Only voxels near `seeds` are evaluated, and only evaluations supported by the field are stored.
/// The result is a shell of written voxels that follows the data, as expected by
/// [`crate::raster3::extract_isosurface`]. Non-finite field values and non-positive weights are
/// omitted.
///
/// This function is independent of field construction, so different field implementations and
/// fields accumulated from several scans can use the same sampling path.
///
/// # Arguments
///
/// * `field`: the field to sample
/// * `seeds`: the positions the band is built around, normally the cloud the field came from
/// * `voxel_size`: the spacing of the grid, in the units of the points
/// * `origin`: the world position of the sample with key `[0, 0, 0]`, normally the world origin
///
/// returns: `Result<SparseGrid3<SdfVoxel>>`, failing if `voxel_size` is not finite and positive or
/// the field reports a support radius that is negative or non-finite
///
/// # Choosing the voxel size
///
/// A cell can be triangulated only when all eight corners were evaluated, and the corners of a cell
/// containing a surface can be up to `sqrt(3)` voxels from it. The field's support radius should
/// therefore be at least approximately `1.75` voxels. A thinner band can create holes. The support
/// radius must also reach the data: if the voxel size is much smaller than the point spacing, most
/// grid samples fall outside every neighborhood.
pub fn compute_sdf_grid(
    field: &impl ImplicitField3,
    seeds: &[Point3],
    voxel_size: f64,
    origin: Point3,
) -> Result<SparseGrid3<SdfVoxel>> {
    let mut grid = SparseGrid3::<SdfVoxel>::new(voxel_size, origin)?;

    let support = field.support_radius();
    if !support.is_finite() || support < 0.0 {
        return Err(format!("The field reports a support radius of {support}").into());
    }

    let keys = grid.activate_blocks_near(seeds, support)?;

    // Fill each block independently. `Block3` stores its voxels in a box, so collecting the blocks
    // and moving them into the grid moves pointers without copying voxel data.
    let filled: Vec<(BlockKey3, Block3<SdfVoxel>)> = keys
        .par_iter()
        .filter_map(|key| {
            let mut block = Block3::<SdfVoxel>::new();
            let mut any = false;

            for local in 0..BLOCK_VOXELS {
                let voxel_key = SparseGrid3::<SdfVoxel>::voxel_key(key, local);
                let position = grid.corner_position(&voxel_key);

                let Some(sample) = field.evaluate(&position) else {
                    continue;
                };

                if !sample.value.is_finite() || sample.weight <= 0.0 {
                    continue;
                }

                // A voxel stores its weight as `f32`. Without this lower bound, a positive weight
                // that underflows during conversion would make the stored voxel appear unwritten.
                // Preserve supported samples by using the smallest representable positive weight.
                let weight = sample.weight.max(f32::MIN_POSITIVE as f64);

                *block.at_mut(local) = SdfVoxel::new(sample.value, weight);
                any |= block.at(local).is_known();
            }

            if any { Some((*key, block)) } else { None }
        })
        .collect();

    for (key, block) in filled {
        grid.insert_block(key, block);
    }

    Ok(grid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::raster3::extract_isosurface;
    use crate::{KdTree3, Mesh3, PointCloud3, Vector3};

    /// A flat patch of points on the z = 0 plane, all facing up.
    fn plane_cloud(half_width: f64, spacing: f64) -> (Vec<Point3>, Vec<UnitVec3>) {
        let n = (half_width / spacing).ceil() as i32;
        let mut points = Vec::new();

        for i in -n..=n {
            for j in -n..=n {
                points.push(Point3::new(i as f64 * spacing, j as f64 * spacing, 0.0));
            }
        }

        let normals = vec![UnitVec3::new_normalize(Vector3::z()); points.len()];
        (points, normals)
    }

    /// Points spread evenly over a true sphere, with their exact normals.
    ///
    /// A tessellated sphere is unsuitable for measuring method bias. Its chordal facets lie inside
    /// the true surface and add a fixed offset; fitting measurements from one shows an additional
    /// constant of approximately 0.006 beside the curvature term. These points lie on the analytic
    /// sphere, isolating the method's bias.
    fn analytic_sphere(n: usize, radius: f64) -> (Vec<Point3>, Vec<UnitVec3>) {
        let golden = std::f64::consts::PI * (3.0 - 5.0_f64.sqrt());

        (0..n)
            .map(|i| {
                let y = 1.0 - (i as f64 / (n as f64 - 1.0)) * 2.0;
                let r = (1.0 - y * y).max(0.0).sqrt();
                let theta = golden * i as f64;
                let u = Vector3::new(theta.cos() * r, y, theta.sin() * r);
                (Point3::from(u * radius), UnitVec3::new_unchecked(u))
            })
            .unzip()
    }

    fn sphere_cloud(radius: f64, spacing: f64) -> PointCloud3 {
        Mesh3::create_sphere(radius, radius * 0.002)
            .expect("sphere creation failed")
            .sample_poisson(spacing, None)
            .expect("sampling failed")
    }

    // ===========================================================================================
    // Construction
    // ===========================================================================================

    #[test]
    fn a_field_rejects_nonsense_parameters() {
        let (points, normals) = plane_cloud(1.0, 0.5);
        let tree = KdTree3::try_new(&points).expect("tree failed");

        assert!(ImlsField3::try_new(&points, &normals, &tree, 0.0, 0.5).is_err());
        assert!(ImlsField3::try_new(&points, &normals, &tree, -1.0, 0.5).is_err());
        assert!(ImlsField3::try_new(&points, &normals, &tree, f64::NAN, 0.5).is_err());
        assert!(ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.0).is_err());
        assert!(ImlsField3::try_new(&points, &normals, &tree, 1.0, f64::INFINITY).is_err());

        let short = normals[..normals.len() - 1].to_vec();
        assert!(ImlsField3::try_new(&points, &short, &tree, 1.0, 0.5).is_err());
    }

    #[test]
    fn confidence_must_match_the_cloud_and_be_sane() {
        let (points, normals) = plane_cloud(1.0, 0.5);
        let tree = KdTree3::try_new(&points).expect("tree failed");
        let build =
            || ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5).expect("field failed");

        let wrong = vec![1.0; points.len() + 1];
        assert!(build().with_confidence(&wrong).is_err());

        let negative = vec![-1.0; points.len()];
        assert!(build().with_confidence(&negative).is_err());

        let nan = vec![f64::NAN; points.len()];
        assert!(build().with_confidence(&nan).is_err());

        let fine = vec![0.5; points.len()];
        assert!(build().with_confidence(&fine).is_ok());
    }

    // ===========================================================================================
    // The field on a plane, where the answer is known exactly
    // ===========================================================================================

    /// Every tangent plane of a flat surface is identical, so their weighted average is independent
    /// of the selected neighbors and weights. The field must return the height above the plane
    /// exactly; any difference is an arithmetic error rather than method bias.
    #[test]
    fn the_field_over_a_plane_is_the_height_above_it() {
        let (points, normals) = plane_cloud(3.0, 0.25);
        let tree = KdTree3::try_new(&points).expect("tree failed");
        let field = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5).expect("field failed");

        for z in [-0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9] {
            for (x, y) in [(0.0, 0.0), (0.37, -1.2), (-2.0, 2.0)] {
                let p = Point3::new(x, y, z);
                let sample = field.evaluate(&p).expect("the field should reach here");

                assert!(
                    (sample.value - z).abs() < 1e-12,
                    "at {p:?} the field read {} rather than {z}",
                    sample.value
                );
                assert!(sample.weight > 0.0);
            }
        }
    }

    #[test]
    fn the_field_reports_nothing_beyond_its_support() {
        let (points, normals) = plane_cloud(3.0, 0.25);
        let tree = KdTree3::try_new(&points).expect("tree failed");
        let field = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5).expect("field failed");

        // Straight up, past the radius.
        assert!(field.evaluate(&Point3::new(0.0, 0.0, 1.5)).is_none());
        assert!(field.evaluate(&Point3::new(0.0, 0.0, -1.5)).is_none());

        // Off the side of the patch, past the radius from its edge.
        assert!(field.evaluate(&Point3::new(6.0, 0.0, 0.0)).is_none());

        assert_eq!(field.support_radius(), 1.0);
    }

    #[test]
    fn the_field_changes_sign_across_the_surface() {
        let (points, normals) = plane_cloud(3.0, 0.25);
        let tree = KdTree3::try_new(&points).expect("tree failed");
        let field = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5).expect("field failed");

        let above = field.evaluate(&Point3::new(0.0, 0.0, 0.3)).expect("above");
        let below = field.evaluate(&Point3::new(0.0, 0.0, -0.3)).expect("below");

        assert!(above.value > 0.0, "the normal side should read positive");
        assert!(below.value < 0.0, "the far side should read negative");
    }

    /// Reversing every normal must reverse the field value while preserving its weight because
    /// normal direction determines only the sign convention.
    #[test]
    fn reversing_the_normals_reverses_the_field() {
        let (points, normals) = plane_cloud(3.0, 0.25);
        let reversed: Vec<UnitVec3> = normals.iter().map(|n| -*n).collect();
        let tree = KdTree3::try_new(&points).expect("tree failed");

        let forward = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5).expect("field");
        let backward = ImlsField3::try_new(&points, &reversed, &tree, 1.0, 0.5).expect("field");

        for z in [-0.7, -0.2, 0.2, 0.7] {
            let p = Point3::new(0.1, -0.3, z);
            let a = forward.evaluate(&p).expect("forward");
            let b = backward.evaluate(&p).expect("backward");

            assert!((a.value + b.value).abs() < 1e-12);
            assert!((a.weight - b.weight).abs() < 1e-12);
        }
    }

    /// A point with zero confidence must not influence the field, which makes confidence weighting
    /// an alternative to filtering the cloud first.
    #[test]
    fn zero_confidence_removes_a_point_from_the_field() {
        // Two patches at different heights. Silencing the upper one must leave the field reading
        // as though only the lower one were there.
        let (lower, lower_normals) = plane_cloud(2.0, 0.25);
        let mut points = lower.clone();
        let mut normals = lower_normals.clone();

        let upper: Vec<Point3> = lower
            .iter()
            .map(|p| p + Vector3::new(0.0, 0.0, 0.4))
            .collect();
        points.extend(upper);
        normals.extend(lower_normals.iter().copied());

        let mut confidence = vec![1.0; points.len()];
        for c in confidence.iter_mut().skip(lower.len()) {
            *c = 0.0;
        }

        let tree = KdTree3::try_new(&points).expect("tree failed");
        let field = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5)
            .expect("field failed")
            .with_confidence(&confidence)
            .expect("confidence failed");

        // Only the lower patch counts, so this reads as height above z = 0.
        for z in [-0.3, 0.0, 0.2] {
            let sample = field.evaluate(&Point3::new(0.0, 0.0, z)).expect("sample");
            assert!(
                (sample.value - z).abs() < 1e-12,
                "at z = {z} the field read {}",
                sample.value
            );
        }

        // With every point assigned zero confidence, no data supports a field value.
        let silent = vec![0.0; points.len()];
        let field = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5)
            .expect("field failed")
            .with_confidence(&silent)
            .expect("confidence failed");
        assert!(field.evaluate(&Point3::origin()).is_none());
    }

    // ===========================================================================================
    // The field on a sphere, where the method's bias shows up
    // ===========================================================================================

    /// Find where the field crosses zero along a ray out from the middle of a sphere.
    fn zero_crossing(field: &impl ImplicitField3, direction: &Vector3, radius: f64) -> f64 {
        let at = |t: f64| {
            field
                .evaluate(&Point3::from(direction * (radius + t)))
                .map(|s| s.value)
        };

        // Stay well inside the field's support, since outside it there is no value to bracket
        // with.
        let reach = 0.4 * field.support_radius();
        let (mut lo, mut hi) = (-reach, reach);
        assert!(at(lo).expect("inside") < 0.0);
        assert!(at(hi).expect("outside") > 0.0);

        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if at(mid).expect("bracketed") < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// The zero level set lies outside a convex surface by an amount proportional to the square of
    /// the weight scale divided by the radius of curvature. This test verifies the coefficient from
    /// the module documentation across several weight scales.
    #[test]
    fn the_curvature_bias_follows_the_documented_formula() {
        let radius = 5.0;
        let (points, normals) = analytic_sphere(60_000, radius);
        let tree = KdTree3::try_new(&points).expect("tree failed");

        let direction = Vector3::new(0.31, -0.62, 0.72).normalize();

        // Measured coefficients on this fixture run from 0.45 to 0.47 across these weight scales,
        // against the 0.46 the module documentation quotes.
        for sigma in [0.2, 0.3, 0.4, 0.5] {
            let field =
                ImlsField3::try_new(&points, &normals, &tree, sigma * 2.0, sigma).expect("field");

            let measured = zero_crossing(&field, &direction, radius);
            let predicted = 0.46 * sigma * sigma / radius;

            assert!(
                measured > 0.0,
                "sigma {sigma}: the bias should push the surface outward, got {measured}"
            );
            assert!(
                (measured - predicted).abs() < 0.15 * predicted,
                "sigma {sigma}: bias was {measured:.6}, formula says {predicted:.6}"
            );
        }
    }

    /// Away from the zero crossing the field should still track the true distance to the sphere,
    /// to within the same curvature bias.
    #[test]
    fn the_field_tracks_the_distance_to_a_sphere() {
        let radius = 5.0;
        let sigma = 0.2;
        let (points, normals) = analytic_sphere(60_000, radius);
        let tree = KdTree3::try_new(&points).expect("tree failed");

        let field = ImlsField3::try_new(&points, &normals, &tree, sigma * 2.0, sigma)
            .expect("field failed");

        let direction = Vector3::new(0.0, 0.0, 1.0);
        let tolerance = 2.0 * 0.46 * sigma * sigma / radius;

        for t in [-0.3, -0.15, 0.15, 0.3] {
            let p = Point3::from(direction * (radius + t));
            let sample = field.evaluate(&p).expect("within the band");

            assert!(
                (sample.value - t).abs() < tolerance,
                "at an offset of {t} the field read {}",
                sample.value
            );
        }
    }

    // ===========================================================================================
    // Sampling the field onto a grid
    // ===========================================================================================

    #[test]
    fn a_grid_over_no_seeds_is_empty() {
        let (points, normals) = plane_cloud(1.0, 0.5);
        let tree = KdTree3::try_new(&points).expect("tree failed");
        let field = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5).expect("field failed");

        let grid = compute_sdf_grid(&field, &[], 0.25, Point3::origin()).expect("grid failed");
        assert_eq!(grid.block_count(), 0);
        assert_eq!(grid.known_count(), 0);
    }

    #[test]
    fn a_grid_rejects_a_nonsense_voxel_size() {
        let (points, normals) = plane_cloud(1.0, 0.5);
        let tree = KdTree3::try_new(&points).expect("tree failed");
        let field = ImlsField3::try_new(&points, &normals, &tree, 1.0, 0.5).expect("field failed");

        assert!(compute_sdf_grid(&field, &points, 0.0, Point3::origin()).is_err());
        assert!(compute_sdf_grid(&field, &points, -1.0, Point3::origin()).is_err());
    }

    /// Every written voxel must lie within the field's support, and every supported voxel must be
    /// written. Together, these conditions define the band. Missing supported voxels usually cause
    /// holes in a reconstruction.
    #[test]
    fn the_band_is_exactly_where_the_field_reaches() {
        let (points, normals) = plane_cloud(1.5, 0.25);
        let tree = KdTree3::try_new(&points).expect("tree failed");
        let support = 0.8;
        let field =
            ImlsField3::try_new(&points, &normals, &tree, support, 0.4).expect("field failed");

        let h = 0.2;
        let grid = compute_sdf_grid(&field, &points, h, Point3::origin()).expect("grid failed");
        assert!(grid.known_count() > 0);

        let point_tree = &tree;
        let mut checked = 0;

        for key in grid.block_keys() {
            for local in 0..BLOCK_VOXELS {
                let voxel_key = SparseGrid3::<SdfVoxel>::voxel_key(&key, local);
                let position = grid.corner_position(&voxel_key);
                let reachable = !point_tree.within(&position, support).is_empty();
                let written = grid.get(&voxel_key).map(|v| v.is_known()).unwrap_or(false);

                assert_eq!(
                    written, reachable,
                    "voxel {voxel_key:?} at {position:?} was written = {written} but reachable = {reachable}"
                );
                checked += 1;
            }
        }

        assert!(checked > 0);
    }

    /// The written voxels around a sphere should approximately fill the shell reached by the field.
    /// This checks that the band is neither empty nor the entire bounding box.
    #[test]
    fn the_band_around_a_sphere_is_about_the_size_it_should_be() {
        let radius = 3.0;
        let cloud = sphere_cloud(radius, 0.1);
        let normals = cloud.point_normals().expect("normals").to_vec();
        let tree = KdTree3::try_new(cloud.points()).expect("tree failed");

        let h = 0.15;
        let support = 2.0 * h;
        let field =
            ImlsField3::try_new(cloud.points(), &normals, &tree, support, h).expect("field failed");

        let grid =
            compute_sdf_grid(&field, cloud.points(), h, Point3::origin()).expect("grid failed");

        // A shell of thickness twice the support wrapped around the sphere.
        let shell = 4.0 * std::f64::consts::PI * radius * radius * 2.0 * support;
        let expected = shell / (h * h * h);
        let actual = grid.known_count() as f64;

        assert!(
            actual > expected * 0.5 && actual < expected * 2.0,
            "the band holds {actual} voxels where the shell volume suggests about {expected:.0}"
        );
    }

    #[test]
    fn sampling_a_field_onto_a_grid_is_reproducible() {
        let cloud = sphere_cloud(3.0, 0.15);
        let normals = cloud.point_normals().expect("normals").to_vec();
        let tree = KdTree3::try_new(cloud.points()).expect("tree failed");
        let field =
            ImlsField3::try_new(cloud.points(), &normals, &tree, 0.4, 0.2).expect("field failed");

        let build = || {
            compute_sdf_grid(&field, cloud.points(), 0.2, Point3::origin()).expect("grid failed")
        };

        let first = build();
        let second = build();

        assert_eq!(first.block_keys(), second.block_keys());
        assert_eq!(first.known_count(), second.known_count());

        for key in first.block_keys() {
            let a = first.block(&key).expect("block");
            let b = second.block(&key).expect("block");
            assert_eq!(a.as_slice(), b.as_slice(), "block {key:?} differs");
        }
    }

    /// Exercises the complete path from an oriented cloud to a mesh, which increment five will
    /// expose through a single call. The sampled sphere must produce a closed mesh of the expected
    /// size.
    #[test]
    fn a_sphere_goes_all_the_way_to_a_closed_mesh() {
        let radius = 5.0;
        let cloud = sphere_cloud(radius, 0.12);
        let normals = cloud.point_normals().expect("normals").to_vec();
        let tree = KdTree3::try_new(cloud.points()).expect("tree failed");

        let h = 0.25;
        let field =
            ImlsField3::try_new(cloud.points(), &normals, &tree, 2.0 * h, h).expect("field failed");

        let grid =
            compute_sdf_grid(&field, cloud.points(), h, Point3::origin()).expect("grid failed");
        let (mesh, stats) = extract_isosurface(&grid).expect("extraction failed");

        assert!(mesh.faces().len() > 1000, "the mesh came out too coarse");
        assert!(stats.cells_visited > 0);

        // Closed and manifold along every edge.
        let mut directed = std::collections::HashMap::new();
        for face in mesh.faces() {
            for k in 0..3 {
                *directed
                    .entry((face[k], face[(k + 1) % 3]))
                    .or_insert(0usize) += 1;
            }
        }
        for (&(a, b), &n) in directed.iter() {
            assert_eq!(n, 1, "directed edge {a}->{b} used {n} times");
            assert_eq!(
                directed.get(&(b, a)).copied().unwrap_or(0),
                1,
                "edge {a}-{b} has no opposing face, so the reconstruction is not closed"
            );
        }

        // Every vertex on the sphere, allowing for the curvature bias and the interpolation error.
        let tolerance = 0.46 * h * h / radius + 0.05 * h;
        let mut worst = 0.0f64;
        for p in mesh.points() {
            worst = worst.max((p.coords.norm() - radius).abs());
        }
        assert!(
            worst < tolerance,
            "worst radial error was {worst:.5}, over a tolerance of {tolerance:.5}"
        );
    }
}
