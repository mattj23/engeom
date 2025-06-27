from __future__ import annotations
from .geom3 import Point3, Mesh, Iso3, Vector3
from numpy.typing import NDArray


class LaserLine:
    def __init__(
            self,
            ray_origin: Point3,
            detect_origin: Point3,
            line_start: Point3,
            line_end: Point3,
            min_range: float,
            max_range: float,
            rays: int,
            angle_limit: float | None = None,
    ):
        """

        :param ray_origin:
        :param detect_origin:
        :param line_start:
        :param line_end:
        :param min_range:
        :param max_range:
        :param rays:
        :param angle_limit:
        """
        ...

    def get_points(self, target: Mesh, obstruction: Mesh | None, iso: Iso3) -> NDArray[float]:
        """

        :param target:
        :param obstruction:
        :param iso:
        :return:
        """
        ...


class PanningLaserLine:
    def __init__(self, laser_line: LaserLine, pan_vector: Vector3, steps: int):
        """
        :param laser_line:
        :param pan_vector:
        :param steps:
        """
        ...

    def get_points(self, target: Mesh, obstruction: Mesh | None, iso: Iso3) -> NDArray[float]:
        """
        :param target:
        :param obstruction:
        :param iso:
        :return:
        """
        ...
