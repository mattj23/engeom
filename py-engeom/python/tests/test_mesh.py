import pytest
from numpy import linalg
from engeom.geom3 import Mesh


def test_mesh_new_offset_vertices():
    m = Mesh.create_sphere(1.0, 100, 100)
    n = m.new_offset_vertices(0.1)

    for v in n.vertices:
        assert pytest.approx(1.1, abs=1e-5) == linalg.norm(v)
