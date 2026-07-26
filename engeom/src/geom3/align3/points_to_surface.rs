use crate::common::dist;
use crate::common::kd_tree::{KdTree, KdTreeSearch};
use crate::common::ransac_tools::ransac_indices;
use crate::geom3::Alignment3;
use crate::geom3::align3::jacobian::{copy_jacobian, point_surf_jacobian};
use crate::geom3::align3::{AlignParams3, AlignSurfMatch3, SurfaceTarget3};
use crate::na::{Dyn, Matrix, Owned, U1, U6, Vector};
use crate::{Iso3, Point3, Result};
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
use rayon::prelude::*;

/// Performs a limited form of RANSAC on a set of points to find a rough alignment to the surface
/// target.
///
/// You provide a number of sample points. It will work best if they are roughly evenly spaced, so
/// consider generating or filtering them with Poisson disk sampling. The RANSAC algorithm will
/// randomly pick three of these points, then find the three nearest neighbors of each, for a total
/// of up to nine points (duplicates are removed) in the subsample. The `points_to_surface3`
/// operation is performed on these subsample points, and if the alignment converges, it iterates
/// through the entire set of original points, transforming them by the result, and counting the
/// number of inliers.  Inliers are points that are less than `inlier_threshold` from the target
/// after transformation.
///
/// The transformation with the highest number of inliers is returned.
///
/// # Arguments
///
/// * `points`: the points to be aligned, in their own local coordinate system
/// * `target`: the target surface entity which the points will be aligned to
/// * `params`: the alignment parameters, see [`AlignParams3`] for details
/// * `inlier_threshold`: the maximum distance between a point and the target that is still
///   considered an inlier
/// * `iterations`: the number of iterations of RANSAC to perform
/// * `ignore_off`: if the surface target can tell if a point does not project directly onto the
///   surface (such as if it projects onto the ends of a boundary), this flag allows such points
///   to be weighted 0.0 to prevent their influence on the alignment.
///
/// returns: Result<Isometry<f64, Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
pub fn ransac_points_to_surface3(
    points: &[Point3],
    target: &impl SurfaceTarget3,
    params: &AlignParams3,
    inlier_threshold: f64,
    iterations: usize,
    ignore_off: bool,
) -> Result<Iso3> {
    let tree = KdTree::try_new(points)?;
    let results = ransac_indices::<3>(iterations, points.len())?
        .par_iter()
        .map(|i| {
            let mut local_indices = vec![i[0], i[1], i[2]];
            for (j, _) in tree.nearest(&points[i[0]], 3) {
                local_indices.push(j);
            }
            for (j, _) in tree.nearest(&points[i[1]], 3) {
                local_indices.push(j);
            }
            for (j, _) in tree.nearest(&points[i[2]], 3) {
                local_indices.push(j);
            }
            // Sort and deduplicate
            local_indices.sort();
            local_indices.dedup();

            let local_points = local_indices.iter().map(|&i| points[i]).collect::<Vec<_>>();

            if let Ok(align) =
                points_to_surface3(&local_points, target, params.clone(), ignore_off, false)
            {
                let mut inliers = 0;
                for test_point in points.iter() {
                    let m = target.align_surf_closest_to(&(align.full() * test_point));
                    if dist(&m, test_point) <= inlier_threshold {
                        inliers += 1;
                    }
                }

                Some((*align.full(), inliers))
            } else {
                None
            }
        })
        .filter_map(|x| x)
        .collect::<Vec<_>>();
    let (best_transform, _) = results.iter().max_by_key(|(_, inliers)| *inliers).unwrap();
    Ok(*best_transform)
}

/// Performs a Levenberg-Marquardt minimization to align a set of points to a surface target.
///
/// # Arguments
///
/// * `points`: the 3D points to be aligned, in their own local coordinate system
/// * `target`: an entity that implements the [`SurfaceTarget3`] trait, which will be the target
///   of the points.
/// * `params`: the alignment parameters, see [`AlignParams3`] for details
/// * `ignore_off`: if the surface target can tell if a point does not project directly onto the
///   surface (such as if it projects onto the ends of a boundary), this flag allows such points to
///   be weighted 0.0 to prevent their influence on the alignment.
/// * `parallel`: if `false` the rayon parallel iteration will _not_ be used.
///
/// returns: Result<Alignment<Unit<Quaternion<f64>>, 3>, Box<dyn Error, Global>>
pub fn points_to_surface3(
    points: &[Point3],
    target: &impl SurfaceTarget3,
    params: AlignParams3,
    ignore_off: bool,
    parallel: bool,
) -> Result<Alignment3> {
    let problem = PointsToSurface::new(points, target, params, ignore_off, parallel);

    let (result, report) = LevenbergMarquardt::new().minimize(problem);

    if report.termination.was_successful() {
        let residuals = result.residuals().unwrap().as_slice().to_vec();
        let c = result.params.current_values();
        let align = Alignment3::new(
            c.transform,
            c.align,
            result.params.local,
            result.params.offset,
            residuals,
        );
        Ok(align)
    } else {
        Err(format!(
            "Failed to align points to surface: {:?}",
            report.termination
        )
        .into())
    }
}

struct PointsToSurface<'a, T: SurfaceTarget3> {
    points: &'a [Point3],
    target: &'a T,
    params: AlignParams3,
    moved: Vec<Point3>,
    closest: Vec<AlignSurfMatch3>,
    residuals: Vec<f64>,
    weights: Vec<f64>,
    ignore_off: bool,
    parallel: bool,
}

impl<'a, T: SurfaceTarget3> PointsToSurface<'a, T> {
    fn new(
        points: &'a [Point3],
        target: &'a T,
        params: AlignParams3,
        ignore_off: bool,
        parallel: bool,
    ) -> Self {
        let mut x = Self {
            points,
            target,
            params,
            moved: vec![Point3::default(); points.len()],
            closest: vec![AlignSurfMatch3::default(); points.len()],
            residuals: vec![0.0; points.len()],
            weights: vec![1.0; points.len()],
            ignore_off,
            parallel,
        };

        x.move_points();
        x
    }

    /// Internally, this moves the points and computes the closest surface point on the mesh to
    /// each.
    fn move_points(&mut self) {
        let current = self.params.current_values();
        let indices = (0..self.points.len()).collect::<Vec<_>>();

        if self.parallel {
            let collected = indices
                .par_iter()
                .map(|&i| {
                    let m = current.transform * self.points[i];
                    let c = self.target.align_surf_closest_to(&m);
                    (i, m, c)
                })
                .collect::<Vec<_>>();
            for (i, m, c) in collected {
                self.moved[i] = m;
                self.closest[i] = c;
            }
        } else {
            for (i, &j) in indices.iter().enumerate() {
                let m = current.transform * self.points[j];
                let c = self.target.align_surf_closest_to(&m);
                self.moved[i] = m;
                self.closest[i] = c;
            }
        }

        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            // The residual is the distance between the test point and the closest point on the
            // mesh surface, adjusted for the direction of the scalar projection.
            self.residuals[i] = dist(p, &c.point) * c.dn(p).signum();

            self.weights[i] = c.weight as f64;

            if self.ignore_off {
                self.weights[i] *= f64::from(c.is_on);
            }
        }
    }
}

impl<'a, T: SurfaceTarget3> LeastSquaresProblem<f64, Dyn, U6> for PointsToSurface<'a, T> {
    type ResidualStorage = Owned<f64, Dyn, U1>;
    type JacobianStorage = Owned<f64, Dyn, U6>;
    type ParameterStorage = Owned<f64, U6>;

    fn set_params(&mut self, x: &Vector<f64, U6, Self::ParameterStorage>) {
        self.params.set_storage(*x);
        self.move_points();
    }

    fn params(&self) -> Vector<f64, U6, Self::ParameterStorage> {
        self.params.get_storage()
    }

    fn residuals(&self) -> Option<Vector<f64, Dyn, Self::ResidualStorage>> {
        let mut res = Matrix::<f64, Dyn, U1, Self::ResidualStorage>::zeros(self.points.len());
        for i in 0..self.points.len() {
            res[i] = self.residuals[i] * self.weights[i];
        }

        Some(res)
    }

    fn jacobian(&self) -> Option<Matrix<f64, Dyn, U6, Self::JacobianStorage>> {
        let current = self.params.current_values();
        let mut jac = Matrix::<f64, Dyn, U6, Self::JacobianStorage>::zeros(self.points.len());
        for (i, (p, c)) in self.moved.iter().zip(self.closest.iter()).enumerate() {
            let values = point_surf_jacobian(p, c, &current) * self.weights[i];
            copy_jacobian(&values, &mut jac, i);
        }

        Some(jac)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::points::{clone_points, mean_point, transform_points};
    use crate::na::{Translation3, UnitQuaternion};
    use crate::tests::engine_blade;
    use crate::{Mesh3, SelectOp, Selection, Vector3};
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn simple_box_disturbed() -> Result<()> {
        // This test is to verify that a simple test against a box that doesn't have large rotations
        // produces a result that is roughly the inverse of the disturbance
        let mesh = Mesh3::create_box(10.0, 5.0, 2.0, false);
        let points = clone_points(&mesh.sample_poisson(0.1, None));
        let disturb = Iso3::from_parts(
            Translation3::new(3.0, 2.0, 1.0),
            UnitQuaternion::from_euler_angles(PI / 8.0, PI / 12.0, PI / 16.0),
        );

        let params = AlignParams3::new_at_origin(None);
        let to_align = transform_points(&points, &disturb);
        let result = points_to_surface3(&to_align, &mesh, params, false, false)?;

        assert_relative_eq!(disturb.inverse(), result.full(), epsilon = 1e-8);
        Ok(())
    }

    #[test]
    fn blade_example() -> Result<()> {
        let mesh = engine_blade();
        let mask = mesh
            .face_select(Selection::None)
            .facing(&Vector3::y(), PI / 4.0, SelectOp::Add)
            .take_mask();
        let expected_points = clone_points(&mesh.sample_poisson(2.0, Some(&mask)));

        let disturb = Iso3::from_parts(
            Translation3::new(-100.0, 150.0, 0.0),
            UnitQuaternion::new(Vector3::new(1.0, 1.0, 1.0).normalize() * PI / 6.0),
        );

        let to_align = transform_points(&expected_points, &disturb);

        let params = AlignParams3::new_at_center(mean_point(&to_align), None);
        let result = points_to_surface3(&to_align, &mesh, params, false, false)?;

        let aligned = transform_points(&to_align, result.full());

        let max_deviation = aligned
            .iter()
            .zip(expected_points.iter())
            .map(|(a, e)| (a - e).norm())
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        assert!(
            max_deviation < 1e-6,
            "Max deviation is too high: {}",
            max_deviation
        );

        Ok(())
    }
}
