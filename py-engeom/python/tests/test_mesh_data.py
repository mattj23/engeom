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


def test_transform_in_place_invalidates_the_cached_points():
    data = MeshData3(triangle_points(), triangle_faces())
    _ = data.points

    data.transform_in_place(Iso3.from_translation(10.0, 0.0, 0.0))

    assert data.points[0] == pytest.approx([10.0, 0.0, 0.0])


def test_transform_in_place_rotates_the_stored_normals():
    cloud = loaded_cloud_data()

    # A quarter turn about +x maps +z onto -y. The angle is in radians.
    cloud.transform_in_place(Iso3.from_rotation(numpy.pi / 2.0, 1.0, 0.0, 0.0))

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


def test_cloud_append_in_place_unions_the_attributes():
    cloud = loaded_cloud_data()
    cloud.append_in_place(loaded_cloud_data())

    assert len(cloud) == 6
    assert cloud.point_stdev.shape == (6,)
    assert cloud.point_colors.shape == (6, 3)


def test_cloud_append_in_place_rejects_a_mismatch_without_modifying_the_target():
    cloud = loaded_cloud_data()
    other = loaded_cloud_data()
    other.set_point_stdev(None)

    with pytest.raises(ValueError):
        cloud.append_in_place(other)

    assert len(cloud) == 3
    assert cloud.point_stdev == pytest.approx([0.001, 0.002, 0.003])


def test_cloud_subset_indices_carries_the_attributes():
    sub = loaded_cloud_data().create_subset_indices([2, 0])

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


# ================================================================================================
# The _in_place / _copy pairs
#
# A `_in_place` name only means something next to a `_copy` sibling, so both halves of every pair
# are checked here: the in-place one mutates the receiver, the copy one does not.
# ================================================================================================


@pytest.mark.parametrize("factory", [loaded_mesh_data, loaded_cloud_data])
def test_transform_copy_leaves_the_receiver_alone(factory):
    original = factory()
    moved = original.transform_copy(Iso3.from_translation(10.0, 0.0, 0.0))

    assert original.points[0] == pytest.approx([0.0, 0.0, 0.0])
    assert moved.points[0] == pytest.approx([10.0, 0.0, 0.0])


@pytest.mark.parametrize("factory", [loaded_mesh_data, loaded_cloud_data])
def test_scale_in_place_scales_points_and_standard_deviations(factory):
    item = factory()
    item.scale_in_place(10.0)

    assert item.points[1] == pytest.approx([10.0, 0.0, 0.0])
    assert item.point_stdev == pytest.approx([0.01, 0.02, 0.03])


@pytest.mark.parametrize("factory", [loaded_mesh_data, loaded_cloud_data])
def test_scale_copy_leaves_the_receiver_alone(factory):
    original = factory()
    scaled = original.scale_copy(10.0)

    assert original.points[1] == pytest.approx([1.0, 0.0, 0.0])
    assert scaled.points[1] == pytest.approx([10.0, 0.0, 0.0])


@pytest.mark.parametrize("factory", [loaded_mesh_data, loaded_cloud_data])
def test_scale_rejects_zero_and_non_finite_factors(factory):
    item = factory()

    for bad in (0.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            item.scale_in_place(bad)

    # The rejected calls must have left it alone.
    assert item.points[1] == pytest.approx([1.0, 0.0, 0.0])


@pytest.mark.parametrize("factory", [loaded_mesh_data, loaded_cloud_data])
def test_a_negative_scale_flips_the_stored_normals(factory):
    mirrored = factory().scale_copy(-1.0)

    assert mirrored.point_normals[0] == pytest.approx([0.0, 0.0, -1.0])
    assert all(s > 0.0 for s in mirrored.point_stdev)


def test_mesh_flip_faces_in_place_reverses_the_winding():
    data = loaded_mesh_data()
    before = data.faces.copy()

    data.flip_faces_in_place()

    assert numpy.array_equal(data.faces[0], [before[0][1], before[0][0], before[0][2]])
    assert data.point_normals[0] == pytest.approx([0.0, 0.0, -1.0])


# ================================================================================================
# Attributes surviving (or refusing) the derived-mesh operations on Mesh
# ================================================================================================


def attributed_mesh() -> Mesh:
    """A box carrying a per-point and a per-face attribute, as the accelerated type."""
    data = MeshData3.from_mesh(Mesh.create_box(1.0, 1.0, 1.0))
    data.set_point_stdev(numpy.full(len(data), 0.01))
    data.set_face_labels(numpy.arange(data.faces.shape[0], dtype=numpy.uint32))
    return data.to_mesh()


def test_split_refuses_an_attributed_mesh_unless_the_loss_is_accepted():
    from engeom.geom3 import Plane3

    mesh = attributed_mesh()
    plane = Plane3(1.0, 0.0, 0.0, 0.0)

    with pytest.raises(ValueError) as info:
        mesh.split(plane)

    assert "point_stdev" in str(info.value)

    # Accepting the loss lets it through.
    negative, positive = mesh.split(plane, allow_attribute_loss=True)
    assert negative is not None and positive is not None

    # A bare mesh needs no flag.
    bare = Mesh.create_box(1.0, 1.0, 1.0)
    assert bare.split(plane) is not None


def test_convex_hull_refuses_an_attributed_mesh_unless_the_loss_is_accepted():
    mesh = attributed_mesh()

    with pytest.raises(ValueError):
        mesh.convex_hull()

    assert mesh.convex_hull(allow_attribute_loss=True) is not None
    assert Mesh.create_box(1.0, 1.0, 1.0).convex_hull() is not None


def test_a_mesh_subset_carries_its_attributes():
    mesh = attributed_mesh()
    sub = mesh.create_from_indices([0, 1])

    data = MeshData3.from_mesh(sub)

    assert data.faces.shape == (2, 3)
    assert numpy.array_equal(data.face_labels, [0, 1])
    assert data.point_stdev.shape == (len(data),)


def test_mesh_ply_round_trip_keeps_the_attributes(tmp_path):
    """The old `Mesh.load_ply` discarded every property the file carried; this one does not."""
    before = attributed_mesh()
    path = tmp_path / "mesh.ply"

    before.save_ply(path)
    after = Mesh.load_ply(path)

    data = MeshData3.from_mesh(after)
    assert data.point_stdev == pytest.approx(0.01)
    assert numpy.array_equal(
        data.face_labels, numpy.arange(data.faces.shape[0], dtype=numpy.uint32)
    )


def test_mesh_write_stl_refuses_to_drop_attributes_silently(tmp_path):
    mesh = attributed_mesh()
    path = tmp_path / "mesh.stl"

    with pytest.raises(IOError):
        mesh.write_stl(path)

    mesh.write_stl(path, allow_attribute_loss=True)
    assert Mesh.load_stl(path) is not None


# ================================================================================================
# Primitive constructors
# ================================================================================================


def test_create_cone_uses_the_radius_and_full_height_it_was_given():
    """The two arguments used to reach parry swapped, so a cone came out with them exchanged."""
    cone = Mesh.create_cone(radius=2.0, height=10.0, steps=32)
    aabb = cone.aabb

    assert aabb.min.x == pytest.approx(-2.0)
    assert aabb.max.x == pytest.approx(2.0)
    assert aabb.min.y == pytest.approx(-5.0)
    assert aabb.max.y == pytest.approx(5.0)


def test_cone_and_cylinder_agree_on_what_height_means():
    cone = Mesh.create_cone(radius=1.0, height=6.0, steps=16)
    cyl = Mesh.create_cylinder(radius=1.0, height=6.0, steps=16)

    assert cone.aabb.max.y == pytest.approx(cyl.aabb.max.y)
    assert cone.aabb.min.y == pytest.approx(cyl.aabb.min.y)


def test_primitives_carry_no_attributes():
    for mesh in [
        Mesh.create_box(1.0, 2.0, 3.0),
        Mesh.create_sphere(1.0, 8, 8),
        Mesh.create_cylinder(1.0, 2.0, 8),
        Mesh.create_cone(1.0, 2.0, 8),
        Mesh.create_circle(1.0, 8),
    ]:
        data = MeshData3.from_mesh(mesh)
        assert data.point_normals is None
        assert data.point_stdev is None
        assert data.face_labels is None


def test_mesh_data_has_the_primitives_too():
    data = MeshData3.create_box(2.0, 4.0, 6.0)

    assert data.points.shape[1] == 3
    assert data.faces.shape[1] == 3
    assert data.to_mesh().aabb.max.z == pytest.approx(3.0)


def test_mesh_scale_copy_rejects_a_zero_factor():
    mesh = Mesh.create_box(1.0, 1.0, 1.0)

    assert mesh.scale_copy(2.0) is not None
    with pytest.raises(ValueError):
        mesh.scale_copy(0.0)


def test_mesh_append_in_place_offsets_the_face_indices():
    data = loaded_mesh_data()
    data.append_in_place(loaded_mesh_data())

    assert len(data) == 6
    assert numpy.array_equal(data.faces, [[0, 1, 2], [3, 4, 5]])
    assert data.point_stdev.shape == (6,)
    assert data.face_labels.shape == (2,)
