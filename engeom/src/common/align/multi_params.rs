//! The dimension-generic core of the multi-body alignment parameterization.
//!
//! A simultaneous alignment of several rigid bodies presents the concatenated parameters of every
//! non-static body as the single flat vector a Levenberg-Marquardt solve wants, whether the bodies
//! are 2D curve groups or 3D meshes. Everything about that bookkeeping is dimension-independent:
//! the layout of the flat vector, the static body's exclusion from it, and the placement of
//! per-body jacobian blocks all work the same way whether each body contributes three parameters
//! or six.
//!
//! [`MultiParams`] holds that shared logic once, generic over the per-body parameterization
//! `P` and the number of parameters `N` it contributes. The dimension-specific parts, which amount
//! to constructing a body about a rotation center and pushing parameter values through to a
//! concrete isometry, are supplied by the [`BodyParams`] implementation on `AlignParams2` and
//! `AlignParams3`.
//!
//! The public faces are the aliases `geom2::align2::MultiParams2` and
//! `geom3::align3::MultiParams3`, whose module docs show the concrete parameter layouts.

use crate::Result;
use crate::na::{DMatrix, DVector, SVector};

/// The per-body parameterization a dimension supplies to [`MultiParams`].
///
/// `N` is the number of parameters one body contributes to the flat vector: three in 2D
/// (`tx`, `ty`, `rz`) and six in 3D (`tx`, `ty`, `tz`, `rx`, `ry`, `rz`).
pub trait BodyParams<const N: usize>: Clone {
    /// The point type of the dimension, used to give each body its own rotation center.
    type Point;

    /// The isometry type of the dimension.
    type Iso: Copy;

    /// The per-body degree-of-freedom locks.
    type Dof: Copy;

    /// The precomputed values a solver's residual and jacobian loops consume.
    type Values;

    /// One body rotating about `center`, posed at `start` (the identity when `None`) with all of
    /// its parameters at zero, so that the parameters describe motion *away from* the starting
    /// pose rather than motion from the origin.
    fn from_posed_center(
        center: &Self::Point,
        start: Option<Self::Iso>,
        dof: Option<Self::Dof>,
    ) -> Self;

    /// Replaces the body's parameters, enforcing its own degree-of-freedom locks on the way
    /// through.
    fn set_storage(&mut self, storage: SVector<f64, N>);

    /// The body's current world transform.
    fn compute_transform(&self) -> Self::Iso;

    /// The body's precomputed alignment values.
    fn compute_values(&self) -> Self::Values;
}

/// The parameters of a simultaneous alignment of several rigid bodies, one of which is held fixed.
///
/// # The static body
///
/// One body is held fixed. Without it the problem is singular in `N` directions, because rigidly
/// moving every body together changes no residual at all: only the *relative* poses are
/// observable. Fixing one body in place removes that freedom and makes its frame the one every
/// result is expressed in. It contributes no parameters, so a set of `n` bodies has `N * (n - 1)`
/// of them.
///
/// # Layout
///
/// Parameters are laid out by body in index order, skipping the static one, `N` at a time.
/// [`MultiParams::column_offset`] is the only thing that needs to know this, and it is what
/// [`MultiParams::add_jacobian_block`] uses to place a body's jacobian columns.
#[derive(Clone, Debug)]
pub struct MultiParams<P: BodyParams<N>, const N: usize> {
    /// The index of the body which is held fixed and contributes no parameters.
    static_i: usize,

    /// One parameterization per body, in body order. The static body's entry is present so that
    /// indexing is uniform, but its storage is never written.
    bodies: Vec<P>,

    /// The concatenated parameters of every non-static body.
    storage: DVector<f64>,
}

impl<P: BodyParams<N>, const N: usize> MultiParams<P, N> {
    /// Builds the parameterization from one per-body parameterization per body.
    ///
    /// Use this when the bodies need individually chosen local origins, working offsets, or
    /// degree-of-freedom locks. [`MultiParams::from_centers`] covers the common case.
    ///
    /// Returns an error if there are fewer than two bodies, or if `static_i` is out of range.
    pub fn new(static_i: usize, bodies: Vec<P>) -> Result<Self> {
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

        let storage = DVector::zeros((bodies.len() - 1) * N);
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
        centers: &[P::Point],
        initial: Option<&[P::Iso]>,
        dof: Option<P::Dof>,
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
            .map(|(i, c)| P::from_posed_center(c, initial.map(|t| t[i]), dof))
            .collect();

        Self::new(static_i, bodies)
    }

    /// The number of bodies, including the static one.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// The total number of free parameters, `N * (body_count - 1)`.
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

    /// The column in the flat parameter vector where a body's parameters begin, or `None` for the
    /// static body, which has none.
    pub fn column_offset(&self, body: usize) -> Option<usize> {
        if body == self.static_i {
            None
        } else if body > self.static_i {
            Some((body - 1) * N)
        } else {
            Some(body * N)
        }
    }

    /// The parameterization of a single body.
    pub fn body(&self, body: usize) -> &P {
        &self.bodies[body]
    }

    /// The current world transform of a body.
    pub fn transform(&self, body: usize) -> P::Iso {
        self.bodies[body].compute_transform()
    }

    /// The precomputed alignment values of every body, in body order.
    ///
    /// A solve needs these once per parameter change and then once per correspondence, so they are
    /// worked out in a batch rather than recomputed inside the residual and jacobian loops.
    pub fn compute_all_values(&self) -> Vec<P::Values> {
        self.bodies.iter().map(|b| b.compute_values()).collect()
    }

    /// Adds a body's jacobian values to the columns it owns.
    ///
    /// This is the form a multi-body residual needs, because a single correspondence touches two
    /// bodies and a body can appear on both sides of it. Overwriting would silently drop one of
    /// the two contributions when a body is matched against itself.
    pub fn add_jacobian_block(
        &self,
        matrix: &mut DMatrix<f64>,
        row: usize,
        body: usize,
        values: &SVector<f64, N>,
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
                let slice = self.storage.fixed_rows::<N>(start).into_owned();
                self.bodies[i].set_storage(slice);
            }
        }
    }
}
