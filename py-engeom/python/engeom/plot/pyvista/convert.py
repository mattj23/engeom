"""
Conversions between `engeom` entities and the PyVista datasets that represent them.

The draw methods are a thin layer over `Plotter.add_mesh`, so nearly all of them come down to
building the right `PolyData` first. Keeping that here means the same conversion serves both the
helper and any caller who wants to hand an `engeom` entity to a PyVista filter or widget directly,
then bring the result back.
"""

from __future__ import annotations

import numpy
import pyvista

from engeom.geom3 import Curve3, CurveGroup3, Mesh3, MeshData3, PointCloud3, Segment3

from .._coerce import Coords3, to_tuple3

__all__ = ["FACE_ID", "POINT_ID", "LineBuilder", "to_mesh3", "to_polydata"]

# The names of the arrays that `to_polydata` stamps onto the datasets it builds, holding each
# face's and each point's index in the `engeom` entity it came from.
#
# PyVista's own filters carry cell and point data through, so these survive operations that
# renumber or discard elements, such as clipping, thresholding, or the extraction a rubber-band
# selection performs. That is what makes it possible to take a selection made in the render window
# and turn it back into indices into the original `engeom` entity.
FACE_ID = "engeom_face_id"
POINT_ID = "engeom_point_id"


def _stamp_ids(data: pyvista.PolyData) -> pyvista.PolyData:
    """
    Record each face's and point's own index on a dataset, so that it can be recovered after a
    filter has renumbered them.

    The arrays are written with `set_array` rather than by assignment, which stores them without
    making them the dataset's active scalars. An active scalar array is what `Plotter.add_mesh`
    reaches for when it is given no color of its own, so assigning these would silently color every
    mesh drawn by its face index.
    """
    data.cell_data.set_array(numpy.arange(data.n_cells, dtype=numpy.uint32), FACE_ID)
    data.point_data.set_array(numpy.arange(data.n_points, dtype=numpy.uint32), POINT_ID)
    return data


def to_polydata(entity) -> pyvista.PolyData:
    """
    Convert an `engeom` entity into the PyVista `PolyData` that represents it.

    This is what the draw methods build on, and it is public so that an `engeom` entity can be
    handed to any PyVista filter or widget that the helper does not wrap, with `to_mesh3` to bring
    a surface back afterwards.

    Faces and points carry their original indices in the `FACE_ID` and `POINT_ID` arrays, which
    PyVista's filters preserve, so elements can still be identified after a filter has renumbered
    or discarded some of them.

    :param entity: the entity to convert. A `Mesh3` or `MeshData3` becomes a triangle surface, a
        `PointCloud3` or an `(n, 3)` array of coordinates becomes a vertex cloud, a `Curve3` or
        `Segment3` becomes a single polyline, and a `CurveGroup3` becomes one polyline per curve.
    :return: the `PolyData` representing the entity.
    :raises TypeError: if the entity is not one of the types listed above.
    """
    if isinstance(entity, (Mesh3, MeshData3)):
        # VTK's face array carries each face's vertex count in front of that face's indices, where
        # an `engeom` mesh keeps a plain `(n, 3)` array of triangles, so the count column is added.
        prefix = numpy.ones((entity.faces.shape[0], 1), dtype=entity.faces.dtype)
        faces = numpy.hstack((prefix * 3, entity.faces))
        return _stamp_ids(pyvista.PolyData(entity.points, faces))

    if isinstance(entity, PointCloud3):
        return _stamp_ids(pyvista.PolyData(entity.points))

    if isinstance(entity, CurveGroup3):
        return _stamp_ids(_lines_polydata(*[c.points for c in entity.curves]))

    if isinstance(entity, Curve3):
        return _stamp_ids(_lines_polydata(entity.points))

    if isinstance(entity, Segment3):
        return _stamp_ids(_lines_polydata(numpy.array([to_tuple3(entity.a), to_tuple3(entity.b)])))

    if numpy.ndim(entity) == 2:
        return _stamp_ids(pyvista.PolyData(numpy.asarray(entity, dtype=numpy.float64)))

    raise TypeError(f"Cannot convert {type(entity).__name__} into a PyVista PolyData")


def to_mesh3(data, is_solid: bool = False) -> Mesh3:
    """
    Convert a PyVista dataset into the `engeom` mesh holding the same surface.

    Anything that is not already a surface has one extracted, and anything whose faces are not
    already triangles is triangulated, since a `Mesh3` is a triangle mesh. Neither step preserves
    the original faces, so the `FACE_ID` array that `to_polydata` stamps is the way to relate a
    filtered result back to the mesh it came from, rather than the face order here.

    :param data: the PyVista dataset to convert.
    :param is_solid: whether the resulting mesh should be treated as a closed solid, which governs
        how the signed side of a surface is determined.
    :return: the `engeom` mesh holding the dataset's surface.
    """
    surface = data if isinstance(data, pyvista.PolyData) else data.extract_surface()
    triangles = surface.triangulate()
    return Mesh3(
        numpy.asarray(triangles.points, dtype=numpy.float64),
        numpy.asarray(triangles.regular_faces, dtype=numpy.uint32),
        is_solid=is_solid,
    )


def _lines_polydata(*runs) -> pyvista.PolyData:
    """
    Build a PyVista `PolyData` holding one connected polyline per run of points.

    Every line the helper draws is built here and handed to `Plotter.add_mesh`, rather than going
    to `Plotter.add_lines`. `add_lines` accepts neither an opacity nor an open `**kwargs`, so
    anything drawn through it cannot take the styling arguments the rest of the helper offers. A
    `PolyData` of line cells renders identically, pixel for pixel, and accepts everything
    `add_mesh` does.

    :param runs: one or more `(n, 3)` arrays of points, each drawn as a single connected polyline.
        A run with fewer than two points has nothing to draw and is skipped.
    :return: the `PolyData` holding one polyline cell per run.
    :raises ValueError: if no run has at least two points.
    """
    arrays = [numpy.asarray(run, dtype=numpy.float64).reshape(-1, 3) for run in runs]
    arrays = [a for a in arrays if len(a) >= 2]
    if not arrays:
        raise ValueError("No run of points with at least two vertices, so there is nothing to draw")

    # VTK's cell array format for a polyline is the point count followed by that many point
    # indices, with the runs packed end to end into a single flat array.
    cells = numpy.empty(sum(len(a) + 1 for a in arrays), dtype=numpy.int64)
    at = 0
    start = 0
    for a in arrays:
        n = len(a)
        cells[at] = n
        cells[at + 1:at + 1 + n] = numpy.arange(start, start + n)
        at += n + 1
        start += n

    return pyvista.PolyData(numpy.concatenate(arrays), lines=cells)


class LineBuilder:
    """
    Accumulates points into runs, to be drawn as one set of connected polylines.

    Points added in sequence join into a single run. `skip` closes the current run and starts a new
    one, so that several disjoint polylines can go into a single actor and restyle as one thing.
    """

    def __init__(self):
        self._runs = [[]]

    def add(self, points: Coords3):
        """
        Append a point to the run currently being built.

        :param points: the point to append, in any coordinate form the helpers accept.
        """
        self._runs[-1].append(to_tuple3(points))

    def skip(self):
        """
        Close the current run, so that the next point added starts a new, disconnected polyline.
        """
        if self._runs[-1]:
            self._runs.append([])

    def build(self) -> pyvista.PolyData:
        """
        Build the `PolyData` holding every run added so far.

        :return: a `PolyData` with one polyline cell per run.
        :raises ValueError: if no run has at least two points.
        """
        return _lines_polydata(*self._runs)
