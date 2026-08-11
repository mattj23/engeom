"""
Tests for the unaccelerated data containers, `MeshData3` and `PointCloud3`.

These focus on the parts the bindings are actually responsible for: the numpy shapes and dtypes
crossing the boundary, the attribute setters accepting and rejecting arrays, the serialization
round trips, and the bridges to the accelerated types.
"""

from pathlib import Path

import numpy
import pytest

from engeom.geom3 import Iso3, Mesh3, MeshData3, Point3, PointCloud3, Vector3
from engeom.common import IndexMask


def triangle_points() -> numpy.ndarray:
    return numpy.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=numpy.float64
    )


def triangle_faces() -> numpy.ndarray:
    return numpy.array([[0, 1, 2]], dtype=numpy.uint32)


def loaded_mesh_data() -> MeshData3:
    """Mesh3 data carrying one attribute of every typed kind."""
    data = MeshData3(triangle_points(), triangle_faces())
    data.set_point_normals(numpy.tile([0.0, 0.0, 1.0], (3, 1)))
    data.set_point_colors(numpy.array([[255, 0, 0], [0, 255, 0], [0, 0, 255]], dtype=numpy.uint8))
    data.set_point_stdev(numpy.array([0.001, 0.002, 0.003]))
    data.set_face_colors(numpy.array([[10, 20, 30]], dtype=numpy.uint8))
    data.set_face_labels(numpy.array([42], dtype=numpy.uint32))
    return data


def loaded_cloud_data() -> PointCloud3:
    cloud = PointCloud3(triangle_points())
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
    cloud = PointCloud3(triangle_points())

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

    assert PointCloud3(triangle_points()).point_stdev is None


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
    after = PointCloud3.load_ply(path)

    assert after.points == pytest.approx(before.points)
    assert after.point_normals == pytest.approx(before.point_normals)
    assert numpy.array_equal(after.point_colors, before.point_colors)
    assert after.point_stdev == pytest.approx(before.point_stdev)


def test_loading_a_mesh_as_a_point_cloud_is_refused(tmp_path):
    path = tmp_path / "actually_a_mesh.ply"
    loaded_mesh_data().save_ply(path)

    with pytest.raises(IOError) as info:
        PointCloud3.load_ply(path)

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


STUD_BOLT_G3D = Path(__file__).parents[3] / "engeom" / "tests" / "data" / "stud-bolt.g3d"


def test_mesh_data_loads_a_g3d_file():
    data = MeshData3.load_g3d(STUD_BOLT_G3D)

    assert data.points.shape == (8565, 3)
    assert data.faces.shape == (16957, 3)
    assert data.points[0] == pytest.approx(
        [3.8822855365478506, -16.53775066421784, -22.36773866300889]
    )
    assert numpy.array_equal(data.faces[0], [0, 44, 23])


def test_mesh_loads_a_g3d_file():
    mesh = Mesh3.load_g3d(STUD_BOLT_G3D)

    assert mesh.points.shape == (8565, 3)
    assert mesh.faces.shape == (16957, 3)


# ================================================================================================
# Bridges to the accelerated types
# ================================================================================================


def test_mesh_data_round_trips_through_mesh():
    before = loaded_mesh_data()
    mesh = before.to_mesh()

    assert isinstance(mesh, Mesh3)
    assert mesh.points.shape == (3, 3)

    after = MeshData3.from_mesh(mesh)

    assert after.points == pytest.approx(before.points)
    assert after.point_stdev == pytest.approx(before.point_stdev)
    assert numpy.array_equal(after.face_labels, before.face_labels)


def test_a_faceless_mesh_cannot_become_a_mesh():
    data = MeshData3(triangle_points(), numpy.zeros((0, 3), dtype=numpy.uint32))

    with pytest.raises(ValueError):
        data.to_mesh()


def test_cloud_spatial_queries_do_not_disturb_attributes():
    """Indexing no longer copies the cloud into a second type, so nothing can be lost on the way to
    a spatial query. This replaced a lossy `to_cloud`/`from_cloud` round trip."""
    cloud = loaded_cloud_data()

    kept = cloud.sample_poisson_disk(0.0001)
    assert len(kept) == 3

    sampled = cloud.extract_poisson_sample(0.0001)
    assert sampled.points == pytest.approx(cloud.points)
    assert sampled.point_stdev == pytest.approx(cloud.point_stdev)
    assert numpy.array_equal(sampled.point_colors, cloud.point_colors)


def test_cloud_spatial_queries_see_mutations():
    """The class caches its k-d tree between queries, so what can go wrong is a query answered from
    a tree built before the points changed. Every mutator drops the cache; this checks the drop
    actually happens rather than trusting the discipline.

    It has to use an overlap query, because that is the only binding which touches the cached tree.
    Poisson sampling builds its own tree over a voxel-downsampled subset and would pass whether or
    not the cache were invalidated, which makes it useless as a check here."""
    other = PointCloud3(numpy.array([[5.0, 0.0, 0.0]]))

    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0]]))
    assert cloud.overlap_points_by_reciprocity(other, 0.001) == [0]

    # Add the point that `other` actually sits on. With a fresh tree the round trip from `other`
    # now lands on index 1, so index 0 stops overlapping and index 1 starts. A stale tree still
    # holds only the original point and would keep answering [0].
    cloud.append_in_place(PointCloud3(numpy.array([[5.0, 0.0, 0.0]])))
    assert cloud.overlap_points_by_reciprocity(other, 0.001) == [1]


def test_cloud_estimate_normals_recovers_a_plane():
    """A plane fit recovers an axis, not a direction, so `must_match` decides the sign. On a flat
    grid the answer should be +/-Z depending only on which way `must_match` points, and the
    confidence should be near 1 because the neighborhood really is planar."""
    xs, ys = numpy.meshgrid(numpy.arange(9.0), numpy.arange(9.0))
    pts = numpy.column_stack([xs.ravel(), ys.ravel(), numpy.zeros(xs.size)])
    cloud = PointCloud3(pts)

    up = numpy.tile([0.0, 0.0, 1.0], (len(pts), 1))
    normals, confidence = cloud.estimate_normals(up, 2.5)

    assert normals.shape == (len(pts), 3)
    assert confidence.shape == (len(pts),)
    assert normals[:, 2] == pytest.approx(1.0, abs=1e-9)
    assert numpy.all(confidence > 0.5)

    # Flip the reference and every normal flips with it.
    down = -up
    flipped, _ = cloud.estimate_normals(down, 2.5)
    assert flipped[:, 2] == pytest.approx(-1.0, abs=1e-9)


def test_cloud_estimate_normals_rejects_a_mismatched_reference():
    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    with pytest.raises(ValueError):
        cloud.estimate_normals(numpy.tile([0.0, 0.0, 1.0], (3, 1)), 2.0)


def test_cloud_extract_subset_points_takes_a_mask():
    cloud = loaded_cloud_data()
    mask = IndexMask.from_indices([0, 2], 3)

    sub = cloud.extract_subset_points(mask)

    assert len(sub) == 2
    assert sub.points == pytest.approx(cloud.points[[0, 2]])
    assert sub.point_stdev == pytest.approx([0.001, 0.003])


def test_cloud_reports_its_size_and_bounds():
    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]]))

    assert cloud.point_count() == 2
    assert len(cloud) == 2
    assert not cloud.is_empty()

    aabb = cloud.compute_aabb()
    assert aabb.min.x == pytest.approx(0.0)
    assert aabb.max.z == pytest.approx(3.0)

    blank = PointCloud3.empty()
    assert blank.is_empty()
    assert len(blank) == 0


def test_cloud_reduce_by_voxel_averages_and_reports_its_work():
    """The reduction creates new points, so the checks are that positions are centroids and that
    the two derived attributes come back through the bindings. Without those getters the coherence
    signal would exist in Rust and be invisible from Python."""
    cloud = PointCloud3(numpy.array([
        [0.0, 0.0, 0.0],
        [0.2, 0.0, 0.0],
        [0.0, 0.4, 0.0],
        [0.2, 0.4, 0.0],
        [5.5, 0.0, 0.0],
    ]))

    out = cloud.reduce_by_voxel(1.0)

    assert len(out) == 2
    assert out.points[0] == pytest.approx([0.1, 0.2, 0.0])
    assert out.points[1] == pytest.approx([5.5, 0.0, 0.0])

    assert list(out.voxel_count) == [4, 1]
    # No normals went in, so there is nothing to report coherence about.
    assert out.voxel_coherence is None


def test_cloud_reduce_by_voxel_reports_coherence_when_normals_are_present():
    cloud = PointCloud3(numpy.array([
        [0.1, 0.0, 0.0],
        [0.2, 0.0, 0.0],
        [5.1, 0.0, 0.0],
        [5.2, 0.0, 0.0],
    ]))
    cloud.set_point_normals(numpy.array([
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]))

    out = cloud.reduce_by_voxel(1.0)

    coherence = out.voxel_coherence
    assert coherence is not None
    assert coherence[0] == pytest.approx(1.0)   # agreeing normals
    assert coherence[1] == pytest.approx(0.0)   # opposed normals cancel


def test_cloud_reduce_by_voxel_rejects_a_nonsense_size():
    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0]]))
    with pytest.raises(ValueError):
        cloud.reduce_by_voxel(0.0)


def test_cloud_overlap_accepts_the_same_cloud_twice():
    """The cached tree is built through a shared borrow specifically so that passing one cloud as
    both arguments works rather than raising a borrow error."""
    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    assert cloud.overlap_points_by_reciprocity(cloud, 0.001) == [0, 1]


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
    sub = loaded_cloud_data().extract_subset_indices([2, 0])

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
    assert "PointCloud3" in repr(loaded_cloud_data())


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
# Attributes surviving (or refusing) the derived-mesh operations on Mesh3
# ================================================================================================


def attributed_mesh() -> Mesh3:
    """A box carrying a per-point and a per-face attribute, as the accelerated type."""
    data = MeshData3.from_mesh(Mesh3.create_box(1.0, 1.0, 1.0))
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
    bare = Mesh3.create_box(1.0, 1.0, 1.0)
    assert bare.split(plane) is not None


def test_convex_hull_refuses_an_attributed_mesh_unless_the_loss_is_accepted():
    mesh = attributed_mesh()

    with pytest.raises(ValueError):
        mesh.convex_hull()

    assert mesh.convex_hull(allow_attribute_loss=True) is not None
    assert Mesh3.create_box(1.0, 1.0, 1.0).convex_hull() is not None


def test_a_mesh_subset_carries_its_attributes():
    mesh = attributed_mesh()
    sub = mesh.extract_subset_faces_from_indices([0, 1])

    data = MeshData3.from_mesh(sub)

    assert data.faces.shape == (2, 3)
    assert numpy.array_equal(data.face_labels, [0, 1])
    assert data.point_stdev.shape == (len(data),)


def test_mesh_ply_round_trip_keeps_the_attributes(tmp_path):
    """The old `Mesh3.load_ply` discarded every property the file carried; this one does not."""
    before = attributed_mesh()
    path = tmp_path / "mesh.ply"

    before.save_ply(path)
    after = Mesh3.load_ply(path)

    data = MeshData3.from_mesh(after)
    assert data.point_stdev == pytest.approx(0.01)
    assert numpy.array_equal(
        data.face_labels, numpy.arange(data.faces.shape[0], dtype=numpy.uint32)
    )


def test_mesh_save_stl_refuses_to_drop_attributes_silently(tmp_path):
    mesh = attributed_mesh()
    path = tmp_path / "mesh.stl"

    with pytest.raises(IOError):
        mesh.save_stl(path)

    mesh.save_stl(path, allow_attribute_loss=True)
    assert Mesh3.load_stl(path) is not None


# ================================================================================================
# Primitive constructors
# ================================================================================================


def test_create_cone_uses_the_radius_and_full_height_it_was_given():
    """The two arguments used to reach parry swapped, so a cone came out with them exchanged."""
    cone = Mesh3.create_cone(radius=2.0, height=10.0, steps=32)
    aabb = cone.aabb

    assert aabb.min.x == pytest.approx(-2.0)
    assert aabb.max.x == pytest.approx(2.0)
    assert aabb.min.y == pytest.approx(-5.0)
    assert aabb.max.y == pytest.approx(5.0)


def test_cone_and_cylinder_agree_on_what_height_means():
    cone = Mesh3.create_cone(radius=1.0, height=6.0, steps=16)
    cyl = Mesh3.create_cylinder(radius=1.0, height=6.0, steps=16)

    assert cone.aabb.max.y == pytest.approx(cyl.aabb.max.y)
    assert cone.aabb.min.y == pytest.approx(cyl.aabb.min.y)


def test_primitives_carry_no_attributes():
    for mesh in [
        Mesh3.create_box(1.0, 2.0, 3.0),
        Mesh3.create_sphere(1.0, 8, 8),
        Mesh3.create_cylinder(1.0, 2.0, 8),
        Mesh3.create_cone(1.0, 2.0, 8),
        Mesh3.create_circle(1.0, 8),
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
    mesh = Mesh3.create_box(1.0, 1.0, 1.0)

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


# ---------------------------------------------------------------------------
# Mask-based subsets
# ---------------------------------------------------------------------------

def two_triangles() -> MeshData3:
    """Two disjoint triangles, so a face selection has points it can drop."""
    points = numpy.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [5.0, 1.0, 0.0],
    ])
    faces = numpy.array([[0, 1, 2], [3, 4, 5]], dtype=numpy.uint32)
    data = MeshData3(points, faces)
    data.set_face_labels(numpy.array([11, 22], dtype=numpy.uint32))
    return data


def test_mesh_data_extract_subset_faces_drops_orphaned_points():
    data = two_triangles()
    mask = IndexMask.from_indices([1], 2)

    subset = data.extract_subset_faces(mask)

    assert subset.faces.shape == (1, 3)
    assert subset.points.shape == (3, 3)
    assert numpy.allclose(subset.points[0], [5.0, 0.0, 0.0])
    # Faces are renumbered against the surviving points.
    assert numpy.array_equal(subset.faces, [[0, 1, 2]])
    # Attributes follow the selection.
    assert numpy.array_equal(subset.face_labels, [22])


def test_mesh_data_extract_subset_points_keeps_only_fully_selected_faces():
    data = two_triangles()
    # Select the first triangle's points plus one point of the second.
    mask = IndexMask.from_indices([0, 1, 2, 3], 6)

    subset = data.extract_subset_points(mask)

    assert subset.points.shape == (4, 3)
    assert subset.faces.shape == (1, 3)
    assert numpy.array_equal(subset.face_labels, [11])


def test_mesh_data_unique_masks_are_inverses_over_a_clean_selection():
    data = two_triangles()
    face_mask = IndexMask.from_indices([0], 2)

    point_mask = data.compute_unique_point_mask(face_mask)
    back = data.compute_unique_face_mask(point_mask)

    assert list(point_mask.to_indices()) == [0, 1, 2]
    assert back == face_mask


def test_mesh_data_subset_rejects_a_mask_of_the_wrong_length():
    data = two_triangles()
    with pytest.raises(ValueError):
        data.extract_subset_faces(IndexMask(6))


def test_mesh_data_and_mesh_extract_the_same_subset():
    """The accelerated and plain types wrap the same subset code, so they must agree."""
    data = two_triangles()
    mask = IndexMask.from_indices([1], 2)

    from_data = data.extract_subset_faces(mask)
    from_mesh = data.to_mesh().extract_subset_faces(mask)

    assert numpy.allclose(from_data.points, from_mesh.points)
    assert numpy.array_equal(from_data.faces, from_mesh.faces)


# ---------------------------------------------------------------------------
# Editing
# ---------------------------------------------------------------------------

def test_mesh_data_empty_starts_with_nothing():
    data = MeshData3.empty()

    assert data.is_empty
    assert data.point_count == 0
    assert data.face_count == 0
    assert data.points.shape == (0, 3)


def test_mesh_data_can_be_built_up_from_empty():
    data = MeshData3.empty()

    a = data.push_point(Point3(0.0, 0.0, 0.0))
    b = data.push_point(Point3(1.0, 0.0, 0.0))
    c = data.push_point(Point3(0.0, 1.0, 0.0))
    f = data.push_face((a, b, c))

    assert (a, b, c, f) == (0, 1, 2, 0)
    assert not data.is_empty
    assert data.point_count == 3
    assert data.face_count == 1
    assert numpy.array_equal(data.faces, [[0, 1, 2]])


def test_mesh_data_push_face_rejects_an_index_with_no_point():
    data = MeshData3.empty()
    data.push_point(Point3(0.0, 0.0, 0.0))

    with pytest.raises(ValueError):
        data.push_face((0, 1, 2))


def test_mesh_data_push_point_is_blocked_by_point_attributes():
    """There would be no value to give the attribute for the new point, so the type refuses."""
    data = loaded_mesh_data()

    with pytest.raises(ValueError):
        data.push_point(Point3(9.0, 9.0, 9.0))


def test_mesh_data_set_point_moves_it():
    data = MeshData3(triangle_points(), triangle_faces())

    data.set_point(1, Point3(5.0, 6.0, 7.0))

    assert numpy.allclose(data.points[1], [5.0, 6.0, 7.0])


def test_mesh_data_set_point_invalidates_the_cached_points():
    """The cached numpy array would otherwise keep reporting the old position."""
    data = MeshData3(triangle_points(), triangle_faces())
    before = data.points

    data.set_point(1, Point3(5.0, 6.0, 7.0))

    assert before is not data.points
    assert numpy.allclose(data.points[1], [5.0, 6.0, 7.0])


def test_mesh_data_set_face_replaces_the_indices():
    data = MeshData3(triangle_points(), triangle_faces())

    data.set_face(0, (2, 1, 0))

    assert numpy.array_equal(data.faces, [[2, 1, 0]])


@pytest.mark.parametrize("call", [
    lambda d: d.set_point(9, Point3(0.0, 0.0, 0.0)),
    lambda d: d.set_face(9, (0, 1, 2)),
    lambda d: d.set_face(0, (0, 1, 9)),
    lambda d: d.remove_point(9),
    lambda d: d.remove_face(9),
])
def test_mesh_data_editing_rejects_out_of_range_indices(call):
    data = MeshData3(triangle_points(), triangle_faces())
    with pytest.raises(ValueError):
        call(data)


def test_mesh_data_remove_face_leaves_the_points_alone():
    data = two_triangles()

    data.remove_face(0)

    assert data.face_count == 1
    # The points the removed face referenced become orphans rather than being renumbered.
    assert data.point_count == 6
    assert numpy.array_equal(data.faces, [[3, 4, 5]])
    assert numpy.array_equal(data.face_labels, [22])


def test_mesh_data_remove_point_drops_its_faces_and_shifts_indices():
    data = two_triangles()

    data.remove_point(0)

    assert data.point_count == 5
    assert data.face_count == 1
    # The surviving face's indices all shifted down by one.
    assert numpy.array_equal(data.faces, [[2, 3, 4]])
    assert numpy.array_equal(data.face_labels, [22])


def test_mesh_data_remove_point_invalidates_the_cached_arrays():
    data = two_triangles()
    points_before = data.points
    faces_before = data.faces

    data.remove_point(0)

    assert points_before is not data.points
    assert faces_before is not data.faces
    assert data.points.shape == (5, 3)


# ---------------------------------------------------------------------------
# Computed quantities
# ---------------------------------------------------------------------------

def test_mesh_data_compute_face_quantities():
    data = MeshData3(triangle_points(), triangle_faces())

    normals = data.compute_face_normals()
    areas = data.compute_face_areas()
    centers = data.compute_face_centers()

    assert normals.shape == (1, 3)
    assert numpy.allclose(normals[0], [0.0, 0.0, 1.0])
    assert areas.shape == (1,)
    assert areas[0] == pytest.approx(0.5)
    assert centers.shape == (1, 3)
    assert numpy.allclose(centers[0], numpy.mean(triangle_points(), axis=0))


def test_mesh_data_compute_point_normals():
    data = MeshData3(triangle_points(), triangle_faces())

    normals = data.compute_point_normals()

    assert normals.shape == (3, 3)
    assert numpy.allclose(normals, numpy.tile([0.0, 0.0, 1.0], (3, 1)))


def test_mesh_data_compute_matches_the_accelerated_type():
    """Both types delegate to the same algorithms, so their results have to agree exactly."""
    data = two_triangles()
    mesh = data.to_mesh()

    assert numpy.array_equal(data.compute_face_areas(), mesh.compute_face_areas())
    assert numpy.array_equal(data.compute_face_centers(), mesh.compute_face_centers())
    assert numpy.array_equal(data.compute_point_normals(), mesh.compute_point_normals())


def test_mesh_data_compute_point_normals_reports_a_point_with_no_normal():
    data = MeshData3.empty()
    for p in triangle_points():
        data.push_point(Point3(*p))
    data.push_point(Point3(9.0, 9.0, 9.0))  # orphan, belongs to no face
    data.push_face((0, 1, 2))

    with pytest.raises(ValueError):
        data.compute_point_normals()


# ---------------------------------------------------------------------------
# Primitive parity with Mesh3
# ---------------------------------------------------------------------------

def test_mesh_data_capsule_and_cylinder_between():
    p0, p1 = Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 10.0)

    capsule = MeshData3.create_capsule(p0, p1, 1.0, 16, 8)
    cylinder = MeshData3.create_cylinder_between(p0, p1, 1.0, 16)

    assert capsule.face_count > 0
    assert cylinder.face_count > 0
    # The capsule adds two hemispherical caps, so it is the larger of the two.
    assert capsule.face_count > cylinder.face_count


def test_mesh_data_rect_beam_between():
    beam = MeshData3.create_rect_beam_between(
        Point3(0.0, 0.0, 0.0), Point3(10.0, 0.0, 0.0), 2.0, 4.0
    )

    assert beam.face_count == 12  # a box is 12 triangles
    lengths = beam.points.max(axis=0) - beam.points.min(axis=0)
    assert numpy.allclose(sorted(lengths), [2.0, 4.0, 10.0])


def test_mesh_data_rect_beam_rejects_an_up_along_the_segment():
    with pytest.raises(ValueError):
        MeshData3.create_rect_beam_between(
            Point3(0.0, 0.0, 0.0), Point3(0.0, 0.0, 10.0), 2.0, 4.0, up=Vector3(0.0, 0.0, 1.0)
        )


@pytest.mark.parametrize("name", ["stanford_bunny_res2", "stanford_bunny_res3", "stanford_bunny_res4"])
def test_mesh_data_bunnies_match_the_accelerated_type(name):
    """Both types read the same embedded asset, so the buffers must be identical."""
    data = getattr(MeshData3, name)()
    mesh = getattr(Mesh3, name)()

    assert numpy.allclose(data.points, mesh.points)
    assert numpy.array_equal(data.faces, mesh.faces)
