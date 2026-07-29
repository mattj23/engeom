import pytest
from numpy import linalg
from engeom.geom3 import Mesh3


def test_mesh_offset_vertices_copy():
    m = Mesh3.create_sphere(1.0, 100, 100)
    n = m.offset_vertices_copy(0.1)

    for v in n.vertices:
        assert pytest.approx(1.1, abs=1e-5) == linalg.norm(v)
