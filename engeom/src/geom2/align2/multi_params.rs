//! The parameterization for a simultaneous alignment of several rigid bodies.
//!
//! This is the multi-body counterpart to [`AlignParams2`]. Where that struct holds the three
//! parameters of one entity moving against a stationary target, this one holds a whole collection
//! of them and presents their concatenation as the single flat parameter vector a
//! Levenberg-Marquardt solve wants.
//!
//! # The static body
//!
//! One body is held fixed. Without it the problem is singular in three directions, because rigidly
//! moving every body together changes no residual at all: only the *relative* poses are
//! observable. Fixing one body in place removes that freedom and makes its frame the one every
//! result is expressed in. It contributes no parameters, so a set of `n` bodies has `3 * (n - 1)`
//! of them.
//!
//! # Layout
//!
//! Parameters are laid out by body in index order, skipping the static one, three at a time. So
//! with four bodies and body 1 static, the vector is
//!
//! ```text
//! [ body 0: tx ty rz | body 2: tx ty rz | body 3: tx ty rz ]
//! ```
//!
//! [`MultiAlignParams2::column_offset`] is the only thing that needs to know this, and it is what
//! [`MultiAlignParams2::add_jacobian_block`] uses to place a body's three jacobian columns.
//!
//! This is the 2D counterpart of `geom3::align3::multi_params`, and is deliberately structured
//! identically to it.

use crate::Result;
use crate::geom2::align2::{AlignOrigin2, AlignParams2, AlignStorage2, AlignValues2, Dof3};
use crate::geom2::{Iso2, Point2};
use crate::na::{DVector, Dyn, Matrix, Owned};

/// The number of parameters contributed by each non-static body.
const PER_BODY: usize = 3;

type Jacobian = Matrix<f64, Dyn, Dyn, Owned<f64, Dyn, Dyn>>;

/// The parameters of a simultaneous alignment of several rigid bodies, one of which is held fixed.
///
/// See the module documentation for the layout of the flat parameter vector and for why a static
/// body is required.
#[derive(Clone, Debug)]
pub struct MultiAlignParams2 {
    /// The index of the body which is held fixed and contributes no parameters.
    static_i: usize,

    /// One parameterization per body, in body order. The static body's entry is present so that
    /// indexing is uniform, but its storage is never written.
    bodies: Vec<AlignParams2>,

    /// The concatenated parameters of every non-static body.
    storage: DVector<f64>,
}

impl MultiAlignParams2 {
    /// Builds the parameterization from one [`AlignParams2`] per body.
    ///
    /// Use this when the bodies need individually chosen local origins, working offsets, or
    /// degree-of-freedom locks. [`MultiAlignParams2::from_centers`] covers the common case.
    ///
    /// Returns an error if there are fewer than two bodies, or if `static_i` is out of range.
    pub fn new(static_i: usize, bodies: Vec<AlignParams2>) -> Result<Self> {
        if bodies.len() < 2 {
            return Err(format!(
                "a multi-body alignment needs at least two bodies, but {} were given",
                bodies.len()
            )
            .into());
        }
        if static_i >= bodies.len() {
            return Err(format!(
                "the static body index {} is out of range for {} bodies",
                static_i,
                bodies.len()
            )
            .into());
        }

        let storage = DVector::zeros((bodies.len() - 1) * PER_BODY);
        let mut item = Self {
            static_i,
            bodies,
            storage,
        };
        item.distribute();
        Ok(item)
    }

    /// Builds the parameterization from a rotation center per body, which is the usual case.
    ///
    /// Each body rotates about its own center and translates along the world axes. Any initial
    /// pose is placed in the body's working offset, so the parameters start at zero and describe
    /// motion *away from* the initial pose rather than motion from the origin. Putting the centers
    /// near the middle of each body keeps the rotation and translation parameters comparably
    /// scaled, which is what makes the solve well conditioned.
    ///
    /// # Arguments
    ///
    /// * `static_i`: the index of the body to hold fixed
    /// * `centers`: one rotation center per body, in that body's own coordinates
    /// * `initial`: an optional initial pose per body. `None` starts every body at the identity.
    /// * `dof`: an optional degree-of-freedom constraint applied to every non-static body
    pub fn from_centers(
        static_i: usize,
        centers: &[Point2],
        initial: Option<&[Iso2]>,
        dof: Option<Dof3>,
    ) -> Result<Self> {
        if let Some(initial) = initial
            && initial.len() != centers.len()
        {
            return Err(format!(
                "there are {} rotation centers but {} initial transforms",
                centers.len(),
                initial.len()
            )
            .into());
        }

        let bodies = centers
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let local = Iso2::translation(c.x, c.y);
                let start = initial.map_or_else(Iso2::identity, |t| t[i]);

                // With `transform = offset * align * local^-1` and the parameters at zero, the
                // transform is `offset * local^-1`. Setting `offset = start * local` therefore
                // starts the body at exactly `start`.
                AlignParams2::new(AlignOrigin2::Local(local), Some(start * local), dof)
            })
            .collect();

        Self::new(static_i, bodies)
    }

    /// The number of bodies, including the static one.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// The total number of free parameters, `3 * (body_count - 1)`.
    pub fn param_count(&self) -> usize {
        self.storage.len()
    }

    /// The concatenated parameter vector.
    pub fn storage(&self) -> &DVector<f64> {
        &self.storage
    }

    /// Replaces the parameter vector and pushes the new values out to the individual bodies.
    ///
    /// Each body's own degree-of-freedom locks are enforced on the way through, so a locked
    /// parameter stays at zero however the solver sets it.
    pub fn set_storage(&mut self, x: &DVector<f64>) {
        self.storage.copy_from(x);
        self.distribute();
    }

    /// The column in the flat parameter vector where a body's three parameters begin, or `None`
    /// for the static body, which has none.
    pub fn column_offset(&self, body: usize) -> Option<usize> {
        if body == self.static_i {
            None
        } else if body > self.static_i {
            Some((body - 1) * PER_BODY)
        } else {
            Some(body * PER_BODY)
        }
    }

    /// The parameterization of a single body.
    pub fn body(&self, body: usize) -> &AlignParams2 {
        &self.bodies[body]
    }

    /// The current world transform of a body.
    pub fn transform(&self, body: usize) -> Iso2 {
        self.bodies[body].compute_transform()
    }

    /// The precomputed alignment values of every body, in body order.
    ///
    /// A solve needs these once per parameter change and then once per correspondence, so they are
    /// worked out in a batch rather than recomputed inside the residual and jacobian loops.
    pub fn compute_all_values(&self) -> Vec<AlignValues2> {
        self.bodies.iter().map(|b| b.compute_values()).collect()
    }

    /// Adds a body's three jacobian values to the columns it owns.
    ///
    /// This is the form a multi-body residual needs, because a single correspondence touches two
    /// bodies and a body can appear on both sides of it. Overwriting would silently drop one of
    /// the two contributions when a body is matched against itself.
    pub fn add_jacobian_block(
        &self,
        matrix: &mut Jacobian,
        row: usize,
        body: usize,
        values: &AlignStorage2,
    ) {
        if let Some(start) = self.column_offset(body) {
            for (k, v) in values.iter().enumerate() {
                matrix[(row, start + k)] += *v;
            }
        }
    }

    /// Pushes the flat parameter vector out to the individual bodies.
    fn distribute(&mut self) {
        for i in 0..self.bodies.len() {
            if let Some(start) = self.column_offset(i) {
                let slice = self.storage.fixed_rows::<PER_BODY>(start).into_owned();
                self.bodies[i].set_storage(slice);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::random_geometry::RandomGeometry2;
    use crate::geom2::Vector2;
    use approx::assert_relative_eq;

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
