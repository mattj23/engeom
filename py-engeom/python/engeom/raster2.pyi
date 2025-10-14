from __future__ import annotations
from typing import List
from numpy.typing import NDArray
from engeom.geom3 import Mesh


class ScalarRaster:
    """
    A class representing a 2D scalar raster grid.
    """

    @property
    def px_size(self) -> float:
        """
        Get the pixel size of the raster grid.
        :return: the pixel size as a float.
        """
        ...

    @property
    def min_z(self) -> float:
        """
        Get the minimum Z value in the raster grid.
        :return: the minimum Z value as a float.
        """
        ...

    @property
    def max_z(self) -> float:
        """
        Get the maximum Z value in the raster grid.
        :return: the maximum Z value as a float.
        """
        ...

    def f_at(self, x: int, y: int) -> float:
        """
        Get the scalar value at the specified grid coordinates.
        :param x: The x-coordinate (column index).
        :param y: The y-coordinate (row index).
        :return: The scalar value at the specified coordinates. If out of bounds, returns NaN.
        """
        ...

    @staticmethod
    def from_serialized_bytes(data: bytes) -> ScalarRaster:
        """
        Create a ScalarRaster instance from serialized byte data.
        :param data: The byte data representing the serialized raster.
        :return: A ScalarRaster instance.
        """
        ...

    def build_depth_mesh(self) -> Mesh:
        """
        Build a 3D mesh representation of the raster grid.
        :return: A Mesh object representing the raster grid in 3D space.
        """
        ...