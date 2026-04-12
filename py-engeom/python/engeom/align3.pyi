from __future__ import annotations
import numpy
from .geom3 import Mesh, Iso3, Point3


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
    - A *local origin* $L$ — the point/orientation around which rotation happens and relative
      to which translation directions are defined.
    - A *working offset* $O$ — an additional transformation applied after the alignment step,
      typically used to encode an initial guess.
    - A `Dof6` constraint that locks selected degrees of freedom.

    The full transformation applied to the test geometry is $O * A * L^{-1}$, where $A$ is
    the alignment transformation produced by the six optimized parameters.

    Use the factory class methods (`at_origin`, `at_center`, `at_local`, `new`) instead of
    constructing this class directly.
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

    @staticmethod
    def at_origin(dof: Dof6 | None = None) -> AlignParams3:
        """
        Create an `AlignParams3` with the local origin and working offset at the world origin.

        Use this when the test geometry is already near the origin and a good starting position,
        and you are not concerned about the numerical stability of rotations.

        :param dof: Optional `Dof6` constraint. If `None`, all six degrees of freedom are active.
        """
        ...

    @staticmethod
    def at_center(x: float, y: float, z: float, dof: Dof6 | None = None) -> AlignParams3:
        """
        Create an `AlignParams3` that rotates the test entity around a given center point.

        The local origin $L$ is placed at `center` with cardinal directions aligned to the world
        axes, and the working offset is set to match. Translations (`tx`, `ty`, `tz`) still act
        along the world axes; rotations happen around `center`.

        Use this when the test geometry is already in a good starting position but you want to
        provide an explicit rotation center for numerical stability (e.g. the geometry is far
        from the world origin).

        :param x: The X coordinate of the center point.
        :param y: The Y coordinate of the center point.
        :param z: The Z coordinate of the center point.
        :param dof: Optional `Dof6` constraint. If `None`, all six degrees of freedom are active.
        """
        ...

    @staticmethod
    def at_local(local: Iso3, dof: Dof6 | None = None) -> AlignParams3:
        """
        Create an `AlignParams3` whose transformation is expressed at a given local origin.

        Both the local origin $L$ and the working offset $O$ are set to `local`. The physical
        interpretation is that `tx`, `ty`, `tz` translate along the local origin's axes, and
        `rx`, `ry`, `rz` rotate around the local origin's center point and axes. Any DOF
        constraints refer to those same local axes.

        Use this when the test geometry is already in a good starting position and you need full
        control over translation directions and rotation axes — for example when applying DOF
        constraints along an arbitrary direction.

        :param local: The `Iso3` defining the local origin.
        :param dof: Optional `Dof6` constraint. If `None`, all six degrees of freedom are active.
        """
        ...

    @staticmethod
    def new(
        local: Iso3 | Point3 | None = None,
        offset: Iso3 | None = None,
        dof: Dof6 | None = None,
    ) -> AlignParams3:
        """
        Create an `AlignParams3` with full, explicit control over all three components.

        This is the most flexible constructor. The simpler factory methods (`at_origin`,
        `at_center`, `at_local`) cover the common cases and should be preferred.

        :param local: `Iso3` defining the local origin $L$, a `Point3` defining the desired center of rotation, or
            `None` to use the world origin.
        :param offset: `Iso3` working offset $O$ applied after the alignment step. Pass `None`
            to use the identity.
        :param dof: Optional `Dof6` constraint. If `None`, all six degrees of freedom are active.
        """
        ...


class Alignment3:
    """
    The result of a completed 3-D alignment operation.

    Stores the component transformations and the per-sample residuals produced by the solver.
    """

    @property
    def full(self) -> Iso3:
        """
        The full transformation from the test entity's space to the target's space.

        This is the composite $O * A * L^{-1}$ and is the transformation you typically apply
        to the test geometry once alignment has converged.
        """
        ...

    @property
    def alignment(self) -> Iso3:
        """
        The alignment transformation $A$ produced by the six optimized parameters
        (`tx`, `ty`, `tz`, `rx`, `ry`, `rz`) expressed about the local origin.
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

    def residuals(self) -> numpy.ndarray:
        """
        Per-sample residuals from the alignment as a 1-D ``float64`` numpy array.

        Each value is the signed distance between a sampled point and the target surface after
        the alignment transformation is applied.
        """
        ...

    def residual_mean(self) -> float:
        """The mean of the alignment residuals."""
        ...

    def residual_mean_std_dev(self) -> tuple[float, float]:
        """The mean and standard deviation of the alignment residuals as a ``(mean, std_dev)`` tuple."""
        ...


def points_to_mesh(
    points: numpy.ndarray,
    mesh: Mesh,
    params: AlignParams3,
) -> Alignment3:
    """
    Align a set of 3-D points to a mesh surface using an iterative least-squares solver.

    :param points: An ``(N, 3)`` ``float64`` numpy array of test points.
    :param mesh: The target `Mesh` surface.
    :param params: An `AlignParams3` that controls the local origin, working offset, and DOF
        constraints for the alignment.
    :returns: An `Alignment3` containing the resulting transformation and residuals.
    :raises ValueError: if the solver fails to converge or the inputs are invalid.
    """
    ...
