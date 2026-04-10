from __future__ import annotations
import numpy
from .engeom import DeviationMode
from .geom3 import Mesh, Iso3, PointCloud


class Dof6:
    """
    Specifies which of the 6 rigid-body degrees of freedom are active during an alignment
    optimization. Each field is a bool: ``True`` means the degree of freedom is free to be
    optimized; ``False`` locks it.

    The three translation DOFs are ``tx``, ``ty``, ``tz`` and the three rotation DOFs are
    ``rx``, ``ry``, ``rz`` (rotations about the X, Y, and Z axes respectively).
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
        Create a ``Dof6`` with explicit control over each degree of freedom.

        All arguments default to ``True`` (unconstrained), so ``Dof6()`` is equivalent to
        ``Dof6.all()``. Pass ``False`` for any axis to lock that DOF during alignment.

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
        """Return a ``Dof6`` with all six degrees of freedom active."""
        ...


def points_to_mesh(
        points: numpy.ndarray[float],
        mesh: Mesh,
        working: Iso3 | None = None,
        mode: DeviationMode | None = None,
        center: object | None = None,
        dof: Dof6 | None = None,
) -> Iso3:
    """
    Perform a Levenberg-Marquardt least-squares optimization to align a set of points to a mesh.

    :param points: array of shape ``(n, 3)`` containing the points to align.
    :param mesh: the target mesh.
    :param working: initial guess for the isometry; defaults to the identity transform.
    :param mode: deviation metric — point-to-point or point-to-plane; defaults to point-to-plane.
    :param center: optional center of rotation hint.
    :param dof: degrees of freedom to optimize; defaults to all six.
    :return: the isometry that best aligns the points to the mesh.
    :raises ValueError: if the optimization fails to converge.
    """
    ...

def points_to_cloud(
        points: numpy.ndarray[float],
        cloud: PointCloud,
        search_radius: float,
        initial: Iso3,
) -> Iso3:
    """
    Align a set of points to a point cloud using nearest-neighbor correspondences.

    Only correspondences within ``search_radius`` of each query point are considered.

    :param points: array of shape ``(n, 3)`` containing the points to align.
    :param cloud: the target point cloud.
    :param search_radius: maximum distance for a correspondence to be accepted.
    :param initial: initial guess for the isometry.
    :return: the isometry that best aligns the points to the cloud.
    :raises ValueError: if the optimization fails to converge.
    """
    ...

def mesh_to_mesh_iterative(
        mesh: Mesh,
        reference: Mesh,
        sample_spacing: float,
        initial: Iso3,
        mode: DeviationMode,
        max_iter: int,
) -> Iso3:
    """
    Iteratively align a mesh to a reference mesh using sampled point correspondences.

    At each iteration the moving mesh is sampled at ``sample_spacing`` intervals, nearest
    points on the reference are found, and a least-squares step is applied until convergence
    or ``max_iter`` iterations are reached.

    :param mesh: the mesh to align.
    :param reference: the reference mesh to align to.
    :param sample_spacing: distance between sample points on the moving mesh.
    :param initial: initial guess for the isometry.
    :param mode: deviation metric — point-to-point or point-to-plane.
    :param max_iter: maximum number of iterations to perform.
    :return: the isometry that best aligns the mesh to the reference mesh.
    :raises ValueError: if the optimization fails to converge.
    """
    ...