from __future__ import annotations
import numpy
from .geom2 import Curve2, CurveGroup2, Iso2, Point2


class Dof3:
    """
    Specifies which of the 3 rigid-body degrees of freedom are active during an alignment
    optimization. Each field is a bool: `True` means the degree of freedom is free to be
    optimized; `False` locks it.

    The two translation DOFs are `tx` and `ty`, and `rz` is the single rotation in the plane.
    """

    tx: bool
    ty: bool
    rz: bool

    def __init__(self, tx: bool = True, ty: bool = True, rz: bool = True) -> None:
        """
        Create a `Dof3` with explicit control over each degree of freedom.

        All arguments default to `True` (unconstrained), so `Dof3()` is equivalent to
        `Dof3.all()`. Pass `False` for any axis to lock that DOF during alignment.

        :param tx: allow translation along X.
        :param ty: allow translation along Y.
        :param rz: allow rotation in the plane.
        """
        ...

    @staticmethod
    def all() -> Dof3:
        """Return a `Dof3` with all three degrees of freedom active."""
        ...


class AlignParams2:
    """
    Parameters for a 2-D rigid-body alignment problem.

    An `AlignParams2` encodes a *local origin* $L$, the point and orientation about which
    rotation happens and relative to which translation directions are defined, and a *working
    offset* $O$. The transformation produced is the composite $O * A * L^{-1}$, where $A$ is the
    motion the three optimized parameters describe.
    """

    def __init__(
        self,
        center: Point2 | None = None,
        local: Iso2 | None = None,
        offset: Iso2 | None = None,
        dof: Dof3 | None = None,
    ) -> None:
        """
        Create an `AlignParams2` describing how a 2D alignment is parameterized.

        The local origin $L$ is selected by supplying at most one of `center` or `local`:

        - `center`: rotation happens about this point, and translations act along the world axes.
        - `local`: rotation happens about, and translations act along, the axes of this full
          `Iso2` frame.
        - neither: the world origin is used.

        Supplying both `center` and `local` raises a `ValueError`.

        :param center: Optional `Point2` rotation center. Mutually exclusive with `local`.
        :param local: Optional `Iso2` local origin frame. Mutually exclusive with `center`.
        :param offset: Optional `Iso2` working offset $O$. Defaults to the local origin frame.
        :param dof: Optional `Dof3` constraint. If `None`, all three DOF are active.
        """
        ...

    @property
    def dof(self) -> Dof3:
        """The degrees-of-freedom constraint currently active on this alignment."""
        ...

    @property
    def local(self) -> Iso2:
        """The local origin transformation $L$."""
        ...

    @property
    def offset(self) -> Iso2:
        """The working offset transformation $O$."""
        ...


class Alignment2:
    """The transformation produced by a 2-D alignment, together with its residuals."""

    @property
    def full_transform(self) -> Iso2:
        """
        The full transformation from the test entity's coordinate system to the target's. This is
        the composite $O * A * L^{-1}$ and is what you apply to the test geometry.
        """
        ...

    @property
    def local_transform(self) -> Iso2:
        """
        The alignment transformation $A$, the motion produced by the three optimized parameters
        expressed in the frame of the local origin. Not the transformation to apply to the test
        geometry; use `full_transform` for that.
        """
        ...

    @property
    def local_origin(self) -> Iso2:
        """The local origin transformation $L$ that was used during alignment."""
        ...

    @property
    def offset(self) -> Iso2:
        """The working offset transformation $O$ that was used during alignment."""
        ...

    def residuals(self) -> numpy.ndarray[tuple[int], numpy.dtype[numpy.float64]]:
        """
        The per-sample residuals, as a 1-D array of `float64`. These are signed distances between
        each sampled point and the target after the alignment transformation is applied.
        """
        ...

    def residual_mean(self) -> float:
        """The mean of the residuals."""
        ...

    def residual_mean_std_dev(self) -> tuple[float, float]:
        """The mean and standard deviation of the residuals."""
        ...


class AlignOutcome2:
    """
    The full outcome of a 2-D alignment: the `Alignment2` itself, plus a record of how the solves
    which produced it terminated.
    """

    @property
    def alignment(self) -> Alignment2:
        """The alignment which was produced."""
        ...

    @property
    def quality(self) -> str:
        """
        The quality of the weakest contributing solve, as `"converged"` or `"unconverged"`. An
        unconverged result is still valid geometry.
        """
        ...

    @property
    def converged(self) -> bool:
        """Whether every contributing solve met a convergence criterion."""
        ...

    @property
    def refinement_rounds(self) -> int:
        """The number of robust refinement rounds which completed."""
        ...

    @property
    def solves(self) -> list[str]:
        """How each contributing solve terminated, the initial solve first."""
        ...

    @property
    def halt(self) -> str | None:
        """Why robust refinement stopped early, or `None` if it ran every round."""
        ...


class MultiAlignOutcome2:
    """
    The full outcome of a simultaneous alignment of several bodies: one `Alignment2` per body,
    plus a shared record of how the solves which produced them terminated.
    """

    @property
    def alignments(self) -> list[Alignment2]:
        """The alignment produced for each body, in the order the bodies were given."""
        ...

    def alignment(self, body: int) -> Alignment2:
        """
        The alignment produced for one body.

        :raises IndexError: if `body` is out of range for the number of bodies.
        """
        ...

    @property
    def quality(self) -> str:
        """The quality of the weakest contributing solve."""
        ...

    @property
    def converged(self) -> bool:
        """Whether every contributing solve met a convergence criterion."""
        ...

    @property
    def refinement_rounds(self) -> int:
        """The number of robust refinement rounds which completed."""
        ...

    @property
    def solves(self) -> list[str]:
        """How each contributing solve terminated."""
        ...

    @property
    def halt(self) -> str | None:
        """Why robust refinement stopped early, or `None` if it ran every round."""
        ...


def points_to_curve(
    points: numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]],
    curve: Curve2,
    params: AlignParams2,
    ignore_off_target: bool = False,
    refinement_steps: int = 4,
    sigma_max: float | None = None,
    point_sigma: list[float] | None = None,
    patience: int = 100,
) -> AlignOutcome2:
    """
    Align a set of 2-D points to a curve by repeatedly projecting them onto their closest position
    on it as the solver moves them.

    By default this is a robust alignment: an initial unweighted solve followed by
    `refinement_steps` rounds of iteratively reweighted least squares using MAGSAC++ weights. Pass
    `refinement_steps=0` for a plain unweighted least-squares alignment.

    :param points: an `(n, 2)` array of the points to align, in their own coordinate system.
    :param curve: the stationary `Curve2` to align to.
    :param params: an `AlignParams2` describing how the alignment is parameterized.
    :param ignore_off_target: weight points at 0.0 when they project past an open end.
    :param refinement_steps: rounds of robust reweighting after the initial solve.
    :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
    :param point_sigma: optional per-point standard deviations, one per input point.
    :param patience: the Levenberg-Marquardt evaluation budget.
    """
    ...


def points_to_group(
    points: numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]],
    group: CurveGroup2,
    params: AlignParams2,
    ignore_off_target: bool = False,
    refinement_steps: int = 4,
    sigma_max: float | None = None,
    point_sigma: list[float] | None = None,
    patience: int = 100,
) -> AlignOutcome2:
    """
    Align a set of 2-D points to a `CurveGroup2`, a collection of disjoint curves treated as one
    rigid entity, such as the loops and open segments a planar section of a mesh produces.

    Each point matches whichever member of the group is closest, and whether it landed past an
    open end is judged against that member.

    :param points: an `(n, 2)` array of the points to align, in their own coordinate system.
    :param group: the stationary `CurveGroup2` to align to.
    :param params: an `AlignParams2` describing how the alignment is parameterized.
    :param ignore_off_target: weight points at 0.0 when they project past an open end of a member.
    :param refinement_steps: rounds of robust reweighting after the initial solve.
    :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
    :param point_sigma: optional per-point standard deviations, one per input point.
    :param patience: the Levenberg-Marquardt evaluation budget.
    """
    ...


def points_to_cloud(
    points: numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]],
    target_points: numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]],
    target_normals: numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]],
    params: AlignParams2,
    max_extrapolation: float,
    target_sigma: list[float] | None = None,
    ignore_off_target: bool = False,
    refinement_steps: int = 4,
    sigma_max: float | None = None,
    point_sigma: list[float] | None = None,
    patience: int = 100,
) -> AlignOutcome2:
    """
    Align a set of 2-D points to an unordered set of measured points which carry normals.

    Each match is the query projected onto the tangent line at its nearest neighbor rather than the
    neighbor itself, which is what lets the alignment resolve below the spacing of the target
    samples.

    :param points: an `(n, 2)` array of the points to align, in their own coordinate system.
    :param target_points: an `(m, 2)` array of the stationary measured positions.
    :param target_normals: an `(m, 2)` array of normals, one per target point.
    :param params: an `AlignParams2` describing how the alignment is parameterized.
    :param max_extrapolation: how far *along the surface* a query may sit from its nearest target
        point and still count as on-surface. The normal component is excluded, because that
        component is the residual the alignment exists to remove. Set it at a small multiple of
        the target's sample spacing.
    :param target_sigma: optional per-target-point standard deviations.
    :param ignore_off_target: weight points at 0.0 when they sit beyond `max_extrapolation`.
    :param refinement_steps: rounds of robust reweighting after the initial solve.
    :param sigma_max: the MAGSAC++ upper noise bound. Estimated from the data when `None`.
    :param point_sigma: optional per-point standard deviations, one per input point.
    :param patience: the Levenberg-Marquardt evaluation budget.
    """
    ...


def multi_curve_adjustment(
    groups: list[CurveGroup2],
    max_distance: float,
    initial: list[Iso2] | None = None,
    max_normal_angle: float | None = None,
    refinement_steps: int = 4,
    sigma_max: float | None = None,
    patience: int = 100,
    dof: Dof3 | None = None,
    sample_spacing: float = 1.0,
    max_corner_angle: float | None = None,
    end_margin: float | None = None,
    filter_distances: float | None = 3.0,
) -> MultiAlignOutcome2:
    """
    Simultaneously align several `CurveGroup2` bodies to each other in one combined solve.

    This is a bundle adjustment rather than a pose graph optimization: one transformation is
    carried for each body except one, which is held fixed, and all of them are solved for at once
    against the raw correspondences. The static body is chosen automatically as the one the others
    reference most broadly.

    :param groups: the bodies taking part. At least two are required.
    :param max_distance: a hard cutoff on correspondence distance, in the units of the geometry.
        Required rather than optional: bodies overlap only partially, so a point in a region
        another body never saw has no meaningful match at any distance, and the opening unweighted
        solve has no defense against one. Choose it from the geometry, and err generous.
    :param initial: an optional starting pose per body. `None` starts every body at the identity.
    :param max_normal_angle: an optional cutoff, in radians, on the angle between a test point's
        normal and its match's.
    :param refinement_steps: rounds of robust reweighting after the initial solve.
    :param sigma_max: the MAGSAC++ upper noise bound, in the units of the geometry.
    :param patience: the Levenberg-Marquardt evaluation budget.
    :param dof: an optional `Dof3` constraint applied to every non-static body.
    :param sample_spacing: the arc length spacing between correspondence samples on each member.
    :param max_corner_angle: the largest turn, in radians, tolerated within one sample spacing of
        a sample, which keeps samples off corners. Defaults to 60 degrees.
    :param end_margin: the distance from either end of an open member within which samples are
        discarded. Defaults to twice `sample_spacing`.
    :param filter_distances: discard candidates further than this many standard deviations above
        the mean candidate distance.
    """
    ...
