import numpy
import pytest
from numpy import linalg
from engeom.geom3 import Mesh3


def test_mesh_offset_points_copy():
    m = Mesh3.create_sphere(1.0, 100, 100)
    n = m.offset_points_copy(0.1)

    for v in n.vertices:
        assert pytest.approx(1.1, abs=1e-5) == linalg.norm(v)


def test_mesh_compute_point_normals_of_a_box_point_along_the_diagonals():
    """
    Angle weighting is invariant to how a flat face is triangulated, so every corner normal of a box
    lands on the diagonal. Area weighting, which this used to do, pulls them off it.
    """
    m = Mesh3.create_box(2.0, 2.0, 2.0)
    normals = m.compute_point_normals()
    diagonal = 1.0 / numpy.sqrt(3.0)

    assert normals.shape == (m.vertices.shape[0], 3)
    for point, normal in zip(m.vertices, normals):
        expected = numpy.sign(point) * diagonal
        assert numpy.allclose(normal, expected, atol=1e-12)


def test_mesh_compute_point_normals_is_cached():
    m = Mesh3.create_sphere(1.0, 20, 20)
    assert m.compute_point_normals() is m.compute_point_normals()
