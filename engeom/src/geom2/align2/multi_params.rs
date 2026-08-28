//! The parameterization for a simultaneous alignment of several rigid bodies.
//!
//! This is the multi-body counterpart to [`AlignParams2`]. Where that struct holds the three
//! parameters of one entity moving against a stationary target, this one holds a whole collection
//! of them and presents their concatenation as the single flat parameter vector a
//! Levenberg-Marquardt solve wants.
//!
//! The bookkeeping is dimension-independent and lives in
//! [`crate::common::align::multi_params::MultiAlignParams`]; this module supplies the 2D per-body
//! parameterization through the [`BodyParams`] implementation on [`AlignParams2`], and the
//! `geom3::align3::multi_params` module does the same for 3D.
//!
//! # Layout
//!
//! One body is held fixed and contributes no parameters (see the generic module's documentation
//! for why), so a set of `n` bodies has `3 * (n - 1)` parameters, laid out by body in index order
//! three at a time. With four bodies and body 1 static, the vector is
//!
//! ```text
//! [ body 0: tx ty rz | body 2: tx ty rz | body 3: tx ty rz ]
//! ```

use crate::common::align::multi_params::{BodyParams, MultiAlignParams};
use crate::geom2::align2::{AlignOrigin2, AlignParams2, AlignStorage2, AlignValues2, Dof3};
use crate::geom2::{Iso2, Point2};

/// The parameters of a simultaneous alignment of several 2D rigid bodies, one of which is held
/// fixed. See the module documentation for the layout of the flat parameter vector.
pub type MultiAlignParams2 = MultiAlignParams<AlignParams2, 3>;

impl BodyParams<3> for AlignParams2 {
    type Point = Point2;
    type Iso = Iso2;
    type Dof = Dof3;
    type Values = AlignValues2;

    fn from_posed_center(center: &Point2, start: Option<Iso2>, dof: Option<Dof3>) -> Self {
        let local = Iso2::translation(center.x, center.y);
        let start = start.unwrap_or_else(Iso2::identity);

        // With `transform = offset * align * local^-1` and the parameters at zero, the transform
        // is `offset * local^-1`. Setting `offset = start * local` therefore starts the body at
        // exactly `start`.
        AlignParams2::new(AlignOrigin2::Local(local), Some(start * local), dof)
    }

    fn set_storage(&mut self, storage: AlignStorage2) {
        AlignParams2::set_storage(self, storage);
    }

    fn compute_transform(&self) -> Iso2 {
        AlignParams2::compute_transform(self)
    }

    fn compute_values(&self) -> AlignValues2 {
        AlignParams2::compute_values(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::RandomGeometry2;
    use crate::geom2::Vector2;
    use crate::na::{DMatrix, DVector};
    use approx::assert_relative_eq;

    type Jacobian = DMatrix<f64>;

    fn centers() -> Vec<Point2> {
        vec![
            Point2::new(1.0, 2.0),
            Point2::new(-4.0, 5.0),
            Point2::new(7.0, -8.0),
            Point2::origin(),
        ]
    }

    #[test]
    fn the_static_body_contributes_no_parameters() {
        let params = MultiAlignParams2::from_centers(1, &centers(), None, None).unwrap();

        assert_eq!(params.body_count(), 4);
        assert_eq!(params.param_count(), 9);
        assert_eq!(params.column_offset(1), None);
    }

    #[test]
    fn columns_are_laid_out_in_body_order_skipping_the_static_one() {
        let params = MultiAlignParams2::from_centers(1, &centers(), None, None).unwrap();

        assert_eq!(params.column_offset(0), Some(0));
        assert_eq!(params.column_offset(1), None);
        assert_eq!(params.column_offset(2), Some(3));
        assert_eq!(params.column_offset(3), Some(6));
    }

    #[test]
    fn a_static_body_at_index_zero_shifts_everything_down() {
        let params = MultiAlignParams2::from_centers(0, &centers(), None, None).unwrap();

        assert_eq!(params.column_offset(0), None);
        assert_eq!(params.column_offset(1), Some(0));
        assert_eq!(params.column_offset(2), Some(3));
        assert_eq!(params.column_offset(3), Some(6));
    }

    #[test]
    fn initial_poses_are_reproduced_exactly_at_zero_parameters() {
        // The initial pose lives in the working offset, so a freshly built parameterization must
        // put every body where it was told to, with all parameters still at zero.
        let mut rg = RandomGeometry2::from_seed(0xb0d1e5);
        let centers = centers();

        for _ in 0..500 {
            let initial: Vec<Iso2> = (0..centers.len()).map(|_| rg.iso2(10.0)).collect();
            let params =
                MultiAlignParams2::from_centers(2, &centers, Some(&initial), None).unwrap();

            assert!(params.storage().iter().all(|v| *v == 0.0));
            for (i, expected) in initial.iter().enumerate() {
                assert_relative_eq!(
                    params.transform(i).to_matrix(),
                    expected.to_matrix(),
                    epsilon = 1e-10
                );
            }
        }
    }

    #[test]
    fn the_static_body_never_moves() {
        let centers = centers();
        let initial: Vec<Iso2> = vec![
            Iso2::translation(1.0, 0.0),
            Iso2::translation(0.0, 2.0),
            Iso2::translation(3.0, 0.0),
            Iso2::identity(),
        ];
        let mut params =
            MultiAlignParams2::from_centers(1, &centers, Some(&initial), None).unwrap();

        let before = params.transform(1);
        params.set_storage(&DVector::from_element(params.param_count(), 0.35));
        let after = params.transform(1);

        assert_relative_eq!(before.to_matrix(), after.to_matrix(), epsilon = 1e-12);
    }

    #[test]
    fn parameters_reach_the_body_they_belong_to() {
        let centers = centers();
        let mut params = MultiAlignParams2::from_centers(1, &centers, None, None).unwrap();

        // Translate body 2 by one unit in x, and nothing else.
        let mut x = DVector::zeros(params.param_count());
        x[params.column_offset(2).unwrap()] = 1.0;
        params.set_storage(&x);

        assert_relative_eq!(
            params.transform(2).translation.vector,
            Vector2::new(1.0, 0.0),
            epsilon = 1e-12
        );
        for other in [0usize, 1, 3] {
            assert_relative_eq!(
                params.transform(other).to_matrix(),
                Iso2::identity().to_matrix(),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn bodies_rotate_about_their_own_centers() {
        // A body's rotation center is the one point a pure rotation parameter leaves in place.
        let centers = centers();
        let mut params = MultiAlignParams2::from_centers(1, &centers, None, None).unwrap();

        let mut x = DVector::zeros(params.param_count());
        // rz on body 0, which is the third of its three parameters.
        x[params.column_offset(0).unwrap() + 2] = 0.7;
        params.set_storage(&x);

        let moved = params.transform(0) * centers[0];
        assert_relative_eq!(moved, centers[0], epsilon = 1e-12);
    }

    #[test]
    fn locked_dof_are_enforced_on_every_body() {
        let centers = centers();
        let dof = Dof3::new(false, true, false);
        let mut params = MultiAlignParams2::from_centers(1, &centers, None, Some(dof)).unwrap();

        params.set_storage(&DVector::from_element(params.param_count(), 0.5));

        for body in [0usize, 2, 3] {
            let storage = params.body(body).storage();
            assert_eq!(storage[0], 0.0, "tx should be locked on body {body}");
            assert_eq!(storage[2], 0.0, "rz should be locked on body {body}");
            assert_eq!(storage[1], 0.5, "ty should be free on body {body}");
        }
    }

    #[test]
    fn jacobian_blocks_land_in_the_right_columns() {
        let params = MultiAlignParams2::from_centers(1, &centers(), None, None).unwrap();
        let mut jac = Jacobian::zeros(3, params.param_count());

        let values = AlignStorage2::new(1.0, 2.0, 3.0);
        params.add_jacobian_block(&mut jac, 0, 2, &values);

        let start = params.column_offset(2).unwrap();
        for k in 0..3 {
            assert_eq!(jac[(0, start + k)], values[k]);
        }
        // Nothing outside that block was touched.
        assert_eq!(jac.row(0).iter().filter(|v| **v != 0.0).count(), 3);
        assert_eq!(jac.row(1).iter().filter(|v| **v != 0.0).count(), 0);
    }

    #[test]
    fn jacobian_blocks_for_the_static_body_are_dropped() {
        let params = MultiAlignParams2::from_centers(1, &centers(), None, None).unwrap();
        let mut jac = Jacobian::zeros(1, params.param_count());

        params.add_jacobian_block(&mut jac, 0, 1, &AlignStorage2::from_element(9.0));

        assert!(jac.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn adding_jacobian_blocks_accumulates_rather_than_overwrites() {
        // A correspondence between two bodies contributes to both, and a body matched against
        // itself would contribute twice to the same columns.
        let params = MultiAlignParams2::from_centers(1, &centers(), None, None).unwrap();
        let mut jac = Jacobian::zeros(1, params.param_count());

        let values = AlignStorage2::new(1.0, 2.0, 3.0);
        params.add_jacobian_block(&mut jac, 0, 0, &values);
        params.add_jacobian_block(&mut jac, 0, 0, &values);

        let start = params.column_offset(0).unwrap();
        for k in 0..3 {
            assert_eq!(jac[(0, start + k)], 2.0 * values[k]);
        }
    }

    #[test]
    fn compute_all_values_covers_every_body_in_order() {
        let params = MultiAlignParams2::from_centers(1, &centers(), None, None).unwrap();
        let values = params.compute_all_values();

        assert_eq!(values.len(), params.body_count());
        for (i, v) in values.iter().enumerate() {
            assert_relative_eq!(
                v.transform.to_matrix(),
                params.transform(i).to_matrix(),
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn degenerate_inputs_are_rejected() {
        let one = vec![Point2::origin()];
        assert!(MultiAlignParams2::from_centers(0, &one, None, None).is_err());

        let centers = centers();
        assert!(MultiAlignParams2::from_centers(9, &centers, None, None).is_err());

        let short = vec![Iso2::identity(); 2];
        assert!(MultiAlignParams2::from_centers(0, &centers, Some(&short), None).is_err());
    }
}
