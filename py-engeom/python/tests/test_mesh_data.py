"""
Tests for the unaccelerated data containers, `MeshData3` and `PointCloudData3`.

These focus on the parts the bindings are actually responsible for: the numpy shapes and dtypes
crossing the boundary, the attribute setters accepting and rejecting arrays, the serialization
round trips, and the bridges to the accelerated types.
"""

import numpy
import pytest

from engeom.geom3 import Iso3, Mesh, MeshData3, PointCloud, PointCloudData3


def triangle_points() -> numpy.ndarray:
    return numpy.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=numpy.float64
    )


def triangle_faces() -> numpy.ndarray:
    return numpy.array([[0, 1, 2]], dtype=numpy.uint32)


def loaded_mesh_data() -> MeshData3:
    """Mesh data carrying one attribute of every typed kind."""
    data = MeshData3(triangle_points(), triangle_faces())
    data.set_point_normals(numpy.tile([0.0, 0.0, 1.0], (3, 1)))
    data.set_point_colors(numpy.array([[255, 0, 0], [0, 255, 0], [0, 0, 255]], dtype=numpy.uint8))
    data.set_point_stdev(numpy.array([0.001, 0.002, 0.003]))
    data.set_face_colors(numpy.array([[10, 20, 30]], dtype=numpy.uint8))
    data.set_face_labels(numpy.array([42], dtype=numpy.uint32))
    return data


def loaded_cloud_data() -> PointCloudData3:
    cloud = PointCloudData3(triangle_points())
    cloud.set_point_normals(numpy.tile([0.0, 0.0, 1.0], (3, 1)))
    cloud.set_point_colors(numpy.array([[255, 0, 0], [0, 255, 0], [0, 0, 255]], dtype=numpy.uint8))
    cloud.set_point_stdev(numpy.array([0.001, 0.002, 0.003]))
    return cloud


# ================================================================================================
# Buffers and dtypes
# ================================================================================================


def test_mesh_data_buffers_have_the_documented_shapes_and_dtypes():
    data = MeshData3(triangle_points(), triangle_faces())

    assert data.points.shape == (3, 3)
    assert data.points.dtype == numpy.float64
    assert data.faces.shape == (1, 3)
    assert data.faces.dtype == numpy.uint32
    assert len(data) == 3


def test_cloud_data_buffers_have_the_documented_shapes_and_dtypes():
    cloud = PointCloudData3(triangle_points())

    assert cloud.points.shape == (3, 3)
    assert cloud.points.dtype == numpy.float64
    assert len(cloud) == 3


def test_attributes_are_none_on_a_bare_container():
    data = MeshData3(triangle_points(), triangle_faces())

    assert data.point_normals is None
    assert data.point_colors is None
    assert data.point_stdev is None
    assert data.face_colors is None
    assert data.face_labels is None

    assert PointCloudData3(triangle_points()).point_stdev is None


def test_attribute_shapes_and_dtypes():
    data = loaded_mesh_data()

    assert data.point_normals.shape == (3, 3)
    assert data.point_normals.dtype == numpy.float64
    assert data.point_colors.shape == (3, 3)
    assert data.point_colors.dtype == numpy.uint8
    assert data.point_stdev.shape == (3,)
    assert data.point_stdev.dtype == numpy.float64
    assert data.face_colors.shape == (1, 3)
    assert data.face_colors.dtype == numpy.uint8
    assert data.face_labels.shape == (1,)
    assert data.face_labels.dtype == numpy.uint32


def test_setting_none_clears_an_attribute():
    data = loaded_mesh_data()
    assert data.point_stdev is not None

    data.set_point_stdev(None)
    assert data.point_stdev is None

    data.set_point_stdev()
    assert data.point_stdev is None


def test_a_wrong_length_attribute_is_rejected():
    data = MeshData3(triangle_points(), triangle_faces())

    with pytest.raises(ValueError):
        data.set_point_stdev(numpy.array([0.1, 0.2]))

    with pytest.raises(ValueError):
        data.set_face_labels(numpy.array([1, 2], dtype=numpy.uint32))

    # The rejected calls must not have stored anything.
    assert data.point_stdev is None
    assert data.face_labels is None


def test_a_zero_length_normal_is_rejected():
    data = MeshData3(triangle_points(), triangle_faces())

    with pytest.raises(ValueError):
        data.set_point_normals(numpy.zeros((3, 3)))


def test_normals_are_normalized_on_the_way_in():
    data = MeshData3(triangle_points(), triangle_faces())
    data.set_point_normals(numpy.tile([0.0, 0.0, 5.0], (3, 1)))

    assert data.point_normals == pytest.approx(numpy.tile([0.0, 0.0, 1.0], (3, 1)))


def test_the_cached_array_is_invalidated_by_a_setter():
    data = loaded_mesh_data()
    before = data.point_stdev.copy()

    data.set_point_stdev(numpy.array([0.9, 0.9, 0.9]))

    assert before == pytest.approx([0.001, 0.002, 0.003])
    assert data.point_stdev == pytest.approx([0.9, 0.9, 0.9])


def test_transform_invalidates_the_cached_points():
    data = MeshData3(triangle_points(), triangle_faces())
    _ = data.points

    data.transform_by(Iso3.from_translation(10.0, 0.0, 0.0))

    assert data.points[0] == pytest.approx([10.0, 0.0, 0.0])


def test_transform_rotates_the_stored_normals():
    cloud = loaded_cloud_data()

    # A quarter turn about +x maps +z onto -y. The angle is in radians.
    cloud.transform_by(Iso3.from_rotation(numpy.pi / 2.0, 1.0, 0.0, 0.0))

    assert cloud.point_normals[0] == pytest.approx([0.0, -1.0, 0.0], abs=1e-12)


# ================================================================================================
# Serialization
# ================================================================================================


@pytest.mark.parametrize("binary", [True, False])
def test_mesh_data_ply_round_trip(tmp_path, binary):
    before = loaded_mesh_data()
    path = tmp_path / "mesh.ply"

    before.save_ply(path, binary=binary)
    after = MeshData3.load_ply(path)

    assert after.points == pytest.approx(before.points)
    assert numpy.array_equal(after.faces, before.faces)
    assert after.point_normals == pytest.approx(before.point_normals)
    assert numpy.array_equal(after.point_colors, before.point_colors)
    assert after.point_stdev == pytest.approx(before.point_stdev)
    assert numpy.array_equal(after.face_colors, before.face_colors)
    assert numpy.array_equal(after.face_labels, before.face_labels)


@pytest.mark.parametrize("binary", [True, False])
def test_cloud_data_ply_round_trip(tmp_path, binary):
    before = loaded_cloud_data()
    path = tmp_path / "cloud.ply"

    before.save_ply(path, binary=binary)
    after = PointCloudData3.load_ply(path)

    assert after.points == pytest.approx(before.points)
    assert after.point_normals == pytest.approx(before.point_normals)
    assert numpy.array_equal(after.point_colors, before.point_colors)
    assert after.point_stdev == pytest.approx(before.point_stdev)


def test_loading_a_mesh_as_a_point_cloud_is_refused(tmp_path):
    path = tmp_path / "actually_a_mesh.ply"
    loaded_mesh_data().save_ply(path)

    with pytest.raises(IOError) as info:
        PointCloudData3.load_ply(path)

    assert "MeshData3" in str(info.value)


def test_stl_refuses_to_silently_drop_attributes(tmp_path):
    path = tmp_path / "mesh.stl"

    with pytest.raises(IOError):
        loaded_mesh_data().save_stl(path)

    # Accepting the loss lets it through, and the geometry survives.
    loaded_mesh_data().save_stl(path, allow_attribute_loss=True)
    after = MeshData3.load_stl(path)

    assert len(after) == 3
    assert after.point_stdev is None


def test_a_bare_mesh_writes_stl_without_the_flag(tmp_path):
    path = tmp_path / "bare.stl"
    MeshData3(triangle_points(), triangle_faces()).save_stl(path)

    assert MeshData3.load_stl(path).faces.shape == (1, 3)


# ================================================================================================
# Bridges to the accelerated types
# ================================================================================================


def test_mesh_data_round_trips_through_mesh():
    before = loaded_mesh_data()
    mesh = before.to_mesh()

    assert isinstance(mesh, Mesh)
    assert mesh.vertices.shape == (3, 3)

    after = MeshData3.from_mesh(mesh)

    assert after.points == pytest.approx(before.points)
    assert after.point_stdev == pytest.approx(before.point_stdev)
    assert numpy.array_equal(after.face_labels, before.face_labels)


def test_a_faceless_mesh_cannot_become_a_mesh():
    data = MeshData3(triangle_points(), numpy.zeros((0, 3), dtype=numpy.uint32))

    with pytest.raises(ValueError):
        data.to_mesh()


def test_cloud_data_round_trips_through_point_cloud():
    before = loaded_cloud_data()
    cloud = before.to_cloud()

    assert isinstance(cloud, PointCloud)
    assert cloud.points.shape == (3, 3)

    after = PointCloudData3.from_cloud(cloud)

    assert after.points == pytest.approx(before.points)
    assert after.point_stdev == pytest.approx(before.point_stdev)
    assert numpy.array_equal(after.point_colors, before.point_colors)


# ================================================================================================
# Cloud operations
# ================================================================================================


def test_cloud_append_unions_the_attributes():
    cloud = loaded_cloud_data()
    cloud.append(loaded_cloud_data())

    assert len(cloud) == 6
    assert cloud.point_stdev.shape == (6,)
    assert cloud.point_colors.shape == (6, 3)


def test_cloud_append_rejects_a_mismatch_without_modifying_the_target():
    cloud = loaded_cloud_data()
    other = loaded_cloud_data()
    other.set_point_stdev(None)

    with pytest.raises(ValueError):
        cloud.append(other)

    assert len(cloud) == 3
    assert cloud.point_stdev == pytest.approx([0.001, 0.002, 0.003])


def test_cloud_subset_by_indices_carries_the_attributes():
    sub = loaded_cloud_data().create_from_indices([2, 0])

    assert len(sub) == 2
    assert sub.points[0] == pytest.approx([0.0, 1.0, 0.0])
    assert sub.point_stdev == pytest.approx([0.003, 0.001])
    assert numpy.array_equal(sub.point_colors, [[0, 0, 255], [255, 0, 0]])


def test_cloned_is_independent():
    original = loaded_cloud_data()
    copy = original.cloned()

    copy.set_point_stdev(numpy.array([9.0, 9.0, 9.0]))

    assert original.point_stdev == pytest.approx([0.001, 0.002, 0.003])
    assert copy.point_stdev == pytest.approx([9.0, 9.0, 9.0])


def test_repr_says_what_it_holds():
    assert "MeshData3" in repr(loaded_mesh_data())
    assert "PointCloudData3" in repr(loaded_cloud_data())
