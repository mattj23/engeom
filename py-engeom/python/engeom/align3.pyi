from __future__ import annotations
import numpy
from .geom3 import Mesh3, Iso3, Point3, PointCloud3


class Dof6:
    """
    Specifies which of the 6 rigid-body degrees of freedom are active during an alignment
    optimization. Each field is a bool: `True` means the degree of freedom is free to be
    optimized; `False` locks it.

    The three translation DOFs are `tx`, `ty`, `tz` and the three rotation DOFs are
    `rx`, `ry``, `rz`` (rotations about the X, Y, and Z axes respectively).
    """

    tx: bool
    ty: bool
    tz: bool
    rx: bool
    ry: bool
    rz: bool

    def __init__(
        self,
        tx: bool = True,
        ty: bool = True,
        tz: bool = True,
        rx: bool = True,
        ry: bool = True,
        rz: bool = True,
    ) -> None:
        """
        Create a `Dof6` with explicit control over each degree of freedom.

        All arguments default to `True` (unconstrained), so `Dof6()` is equivalent to
        `Dof6.all()`. Pass `False` for any axis to lock that DOF during alignment.

        :param tx: allow translation along X.
        :param ty: allow translation along Y.
        :param tz: allow translation along Z.
        :param rx: allow rotation about X.
        :param ry: allow rotation about Y.
        :param rz: allow rotation about Z.
        """
        ...

    @staticmethod
    def all() -> Dof6:
        """Return a `Dof6` with all six degrees of freedom active."""
        ...


class AlignParams3:
    """
    Parameters for a 3-D rigid-body alignment problem.

    An `AlignParams3` encodes:
    - A *local origin* $L$ - the point/orientation around which rotation happens and relative
      to which translation directions are defined.
    - A *working offset* $O$ - an additional transformation applied after the alignment step,
      typically used to encode an initial guess.
    - A `Dof6` constraint that locks selected degrees of freedom.

    The full transformation applied to the test geometry is $O * A * L^{-1}$, where $A$ is
    the alignment transformation produced by the six optimized parameters.
    """

    @property
    def dof(self) -> Dof6:
        """The degrees-of-freedom constraint currently active on this alignment."""
        ...

    @property
    def local(self) -> Iso3:
        """The local origin transformation $L$."""
        ...

    @property
    def offset(self) -> Iso3:
        """The working offset transformation $O$."""
        ...

    def __init__(
        self,
        center: Point3 | None = None,
        local: Iso3 | None = None,
        offset: Iso3 | None = None,
        dof: Dof6 | None = None,
    ) -> None:
        """
        Create an `AlignParams3` describing how a 3-D alignment is parameterized.

        The local origin $L$ is selected by supplying at most one of `center` or `local`:

        - `center`: rotations happen about this point, and translations act along the world axes.
        - `local`: rotations happen about, and translations act along, the axes of this full
          `Iso3` frame. Use this for full control over the rotation center and translation
          directions, for example when applying DOF constraints along an arbitrary direction.
        - neither: the world origin is used, for geometry already near the origin where the
          numerical stability of rotations is not a concern.

        If `offset` is not given, it defaults to the local origin frame, so the test geometry
        starts in place and the alignment happens about that origin. Only pass an explicit
        `offset` (including the identity) if you specifically need the raw $O * A * L^{-1}$
        behavior where the geometry is displaced by $L^{-1}$ before alignment.

        :param center: Optional `Point3` rotation center. Mutually exclusive with `local`.
        :param local: Optional `Iso3` local origin frame. Mutually exclusive with `center`.
        :param offset: Optional `Iso3` working offset $O$. Defaults to the local origin frame.
        :param dof: Optional `Dof6` constraint. If `None`, all six degrees of freedom are active.
        :raises ValueError: if both `center` and `local` are supplied.
        """
        ...


class Align3:
    """
    The result of a completed 3-D alignment operation.

    Stores the component transformations and the per-sample residuals produced by the solver.
    """

    @property
    def full_transform(self) -> Iso3:
        """
        The full transformation from the test entity's space to the target's space.

        This is the composite $O * A * L^{-1}$ and is the transformation you typically apply
        to the test geometry once alignment has converged.
        """
        ...

    @property
    def local_transform(self) -> Iso3:
        """
        The alignment transformation $A$ produced by the six optimized parameters
        (`tx`, `ty`, `tz`, `rx`, `ry`, `rz`), expressed in the frame of the local origin.

        This is not the transformation to apply to the test geometry; use `full_transform` for
        that. Reading $O * A * L^{-1}$ right to left, $L^{-1}$ puts a point into the local
        origin's frame, $A$ moves it while it is there, and $O$ maps it back out, so $A$ is only
        meaningful applied to local-frame coordinates.
        """
        ...

    @property
    def local_origin(self) -> Iso3:
        """The local origin transformation $L$ that was active during alignment."""
        ...

    @property
    def offset(self) -> Iso3:
        """The working offset transformation $O$ that was active during alignment."""
        ...

    @property
    def residuals(self) -> numpy.ndarray:
        """
        Per-sample residuals from the alignment as a 1-D ``float64`` numpy array.

        Each value is the signed distance between a sampled point and the target surface after
        the alignment transformation is applied.
        """
        ...

    @property
    def residual_mean(self) -> float:
        """The mean of the alignment residuals."""
        ...

    @property
    def residual_mean_std_dev(self) -> tuple[float, float]:
        """The mean and standard deviation of the alignment residuals as a ``(mean, std_dev)`` tuple."""
        ...


class AlignOutcome3:
    """
    The full outcome of a 3-D alignment: the `Align3` itself, plus a record of how the solves
    which produced it terminated.

    An alignment is only ever returned when it is usable, so receiving one of these already means
    there is a real answer to work with. What it adds is the ability to tell a proven convergence
    from a merely plausible one, and to see whether robust refinement ran to completion.
    """

    @property
    def alignment(self) -> Align3:
        """The alignment which was produced."""
        ...

    @property
    def quality(self) -> str:
        """
        The quality of the weakest solve that contributed to the result.

        ``"converged"`` means every contributing solve met a convergence criterion.
        ``"unconverged"`` means at least one ran out of its evaluation budget, so the parameters
        are the best the solver found but convergence was never demonstrated. The alignment is
        valid geometry either way.

        An unconverged result is routine rather than alarming here: the correspondences are
        re-established every time the parameters change, so near the solution a point close to an
        edge can flip between two matches indefinitely without the criteria ever being satisfied.
        Raising `patience` will not help in that case.
        """
        ...

    @property
    def converged(self) -> bool:
        """Whether every solve that contributed to the result met a convergence criterion."""
        ...

    @property
    def refinement_rounds(self) -> int:
        """The number of robust refinement rounds which completed and contributed to the result."""
        ...

    @property
    def solves(self) -> list[str]:
        """
        How each contributing solve terminated, beginning with the initial solve and followed by
        one entry per completed refinement round.

        Values are ``"converged(ftol)"``, ``"converged(xtol)"``, ``"converged(ftol,xtol)"``,
        ``"residuals_zero"``, ``"orthogonal"``, or ``"lost_patience"``. A solve that broke down
        never appears here, because its result is rolled back rather than kept; see `halt`.
        """
        ...

    @property
    def halt(self) -> str | None:
        """
        Why robust refinement stopped before completing every requested round, or ``None`` if it
        ran them all.

        ``"no_noise_estimate"`` means the residual spread collapsed and there was nothing left to
        reweight. ``"underdetermined(...)"`` means reweighting would have left fewer weighted
        points than free parameters. ``"solve_failed(...)"`` means a refinement round broke down
        and was rolled back to the previous round's result.
        """
        ...


def points_to_mesh(
    points: numpy.ndarray,
    mesh: Mesh3,
    params: AlignParams3,
    ignore_off_target: bool = False,
    refinement_steps: int = 4,
    sigma_max: float | None = None,
    point_sigma: numpy.ndarray | list[float] | None = None,
    patience: int = 100,
) -> AlignOutcome3:
    """
    Align a set of 3-D points to a mesh surface, repeatedly projecting them onto their closest
    position on the surface as the solver moves them.

    By default this is a robust alignment: an initial unweighted solve followed by
    ``refinement_steps`` rounds of iteratively reweighted least squares using MAGSAC++ weights.
    Pass ``refinement_steps=0`` for a plain unweighted least-squares alignment.

    :param points: An ``(N, 3)`` ``float64`` numpy array of test points.
    :param mesh: The target `Mesh3` surface. If it carries per-vertex ``point_stdev``, that
        uncertainty is used automatically, interpolated to each match position.
    :param params: An `AlignParams3` that controls the local origin, working offset, and DOF
        constraints for the alignment.
    :param ignore_off_target: Weight points at 0.0 when they do not project onto the surface.
    :param refinement_steps: Rounds of robust reweighting after the initial solve. Zero disables
        robust weighting entirely.
    :param sigma_max: The MAGSAC++ upper noise bound. When ``None`` it is estimated from the
        residuals of the initial solve via the median absolute deviation. In units of sigma when
        any uncertainty is supplied, otherwise in the units of the geometry.
    :param point_sigma: Optional per-point standard deviations, one per input point. Combines in
        quadrature with any uncertainty the mesh reports. Every entry must be finite and positive.
    :param patience: The Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
        count. Must be greater than zero.
    :returns: An `AlignOutcome3` carrying the alignment and how the solves terminated.
    :raises ValueError: only when there is no answer at all, meaning the arguments were rejected
        or the initial solve broke down. A solve that merely exhausts its evaluation budget
        returns normally with ``quality == "unconverged"``.
    """
    ...


def points_to_cloud(
    points: numpy.ndarray,
    cloud: PointCloud3,
    params: AlignParams3,
    max_extrapolation: float,
    ignore_off_target: bool = False,
    refinement_steps: int = 4,
    sigma_max: float | None = None,
    point_sigma: numpy.ndarray | list[float] | None = None,
    patience: int = 100,
) -> AlignOutcome3:
    """
    Align a set of 3-D points to a point cloud, projecting them onto the tangent plane at their
    nearest neighbour as the solver moves them.

    This is the point-cloud counterpart of `points_to_mesh` and takes the same options, but a cloud
    is only samples of a surface rather than a surface itself, which brings two differences.

    The match is the query projected onto the tangent plane at the nearest cloud point, not that
    point itself. Matching to the nearest point would leave a residual floor of roughly half the
    sample spacing even at a perfect pose, since a test point between samples can get no closer
    than the gap between them. On an engine blade sampled at 2 mm, that is the difference between
    recovering a pose to 4.7e-6 mm and to 2.1e-2 mm.

    Because of that, the cloud **must carry per-point normals**: a normal supplies both the tangent
    plane and the sign of the residual, and neither can be recovered from positions alone. Use
    `PointCloud3.estimate_normals` if the cloud does not already have them.

    :param points: an ``(N, 3)`` ``float64`` numpy array of test points.
    :param cloud: the stationary `PointCloud3` to align to, which must carry normals. If it carries
        ``point_stdev`` that uncertainty is used automatically, and if it carries
        ``voxel_coherence`` from `reduce_by_voxel` that becomes the per-match weight, so voxels
        which straddled an edge count for less.
    :param params: an `AlignParams3` that controls the local origin, working offset, and DOF
        constraints for the alignment.
    :param max_extrapolation: how far *laterally* a point may sit from the nearest cloud point and
        still count as on-surface. This bounds the in-plane distance only, never the distance along
        the normal, because that component is the misalignment being solved for. It exists to catch
        points which have wandered past the edge of the cloud, where the tangent plane is an
        extrapolation into space nothing was measured in. Set it at a small multiple of the sample
        spacing. Points beyond it are still matched, and are only discarded if
        ``ignore_off_target`` is set.
    :param ignore_off_target: weight points at 0.0 when they fall beyond ``max_extrapolation``.
    :param refinement_steps: rounds of robust reweighting after the initial solve. Zero disables
        robust weighting entirely.
    :param sigma_max: the MAGSAC++ upper noise bound. When ``None`` it is estimated from the
        residuals of the initial solve via the median absolute deviation.
    :param point_sigma: optional per-point standard deviations, one per input point. Combines in
        quadrature with any uncertainty the cloud reports.
    :param patience: the Levenberg-Marquardt evaluation budget, as a multiplier on the parameter
        count. Must be greater than zero.
    :returns: an `AlignOutcome3` carrying the alignment and how the solves terminated.
    :raises ValueError: only when there is no answer at all, meaning the arguments were rejected,
        the cloud has no normals, or the initial solve broke down. A solve that merely exhausts its
        evaluation budget returns normally with ``quality == "unconverged"``.
    """
    ...
