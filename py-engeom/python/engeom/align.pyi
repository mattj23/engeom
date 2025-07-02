from __future__ import annotations
import numpy
from .engeom import DeviationMode
from .geom3 import Mesh, Iso3, PointCloud


def points_to_mesh(
        points: numpy.ndarray[float],
        mesh: Mesh,
        initial: Iso3,
        mode: DeviationMode
) -> Iso3:
    """
    Perform a Levenberg-Marquardt, least squares optimization to align a set of points to a mesh. This will return the
    isometry that best aligns the points to the mesh, or will throw an exception if the optimization fails.

    :param points: a numpy array of shape (n, 3) containing the points to align.
    :param mesh: the mesh to align the points to.
    :param initial: the initial guess for the isometry. This will be used as the starting point for the optimization.
    :param mode: the mode to use for the deviation calculation. This will determine how the deviation between the points
    and the mesh is calculated.
    :return: the isometry that best aligns the points to the mesh.
    """
    ...

def points_to_cloud(
        points: numpy.ndarray[float],
        cloud: PointCloud,
        search_radius: float,
        initial: Iso3,
) -> Iso3:
    ...