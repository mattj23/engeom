from __future__ import annotations
import numpy
from .engeom import DeviationMode
from .geom3 import Mesh, Iso3, PointCloud


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
