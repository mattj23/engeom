"""
Tests for the `PointCloud3` bindings.

These focus on what the binding layer is responsible for: the numpy shapes and dtypes crossing the
boundary, the attribute setters, the serialization round trip, and the spatial operations reaching
the Rust side with the right arguments. The geometry itself is tested on the Rust side.

Behaviour `PointCloud3` shares with `MeshData3`, notably the `_in_place` / `_copy` pairs, is tested
against both containers together in `test_mesh_data.py` rather than duplicated here.

You should know that the the class caches a k-d tree between spatial queries and drops it on any
mutation, so a test that mutates and then queries is checking the invalidation rather than the 
query itself.
"""

from __future__ import annotations

import numpy
import pytest

from engeom.common import IndexMask
from engeom.geom3 import Iso3, Point3, PointCloud3


def triangle_points() -> numpy.ndarray:
    return numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


def loaded_cloud_data() -> PointCloud3:
    cloud = PointCloud3(triangle_points())
    cloud.set_point_normals(numpy.tile([0.0, 0.0, 1.0], (3, 1)))
    cloud.set_point_colors(numpy.array([[255, 0, 0], [0, 255, 0], [0, 0, 255]], dtype=numpy.uint8))
    cloud.set_point_stdev(numpy.array([0.001, 0.002, 0.003]))
    return cloud


# ================================================================================================
# Buffers, attributes and serialization
# ================================================================================================


def test_data_buffers_have_the_documented_shapes_and_dtypes():
    cloud = PointCloud3(triangle_points())

    assert cloud.points.shape == (3, 3)
    assert cloud.points.dtype == numpy.float64
    assert len(cloud) == 3


def test_transform_in_place_rotates_the_stored_normals():
    cloud = loaded_cloud_data()

    # A quarter turn about +x maps +z onto -y. The angle is in radians.
    cloud.transform_in_place(Iso3.from_rotation(numpy.pi / 2.0, 1.0, 0.0, 0.0))

    assert cloud.point_normals[0] == pytest.approx([0.0, -1.0, 0.0], abs=1e-12)


@pytest.mark.parametrize("binary", [True, False])
def test_data_ply_round_trip(tmp_path, binary):
    before = loaded_cloud_data()
    path = tmp_path / "cloud.ply"

    before.save_ply(path, binary=binary)
    after = PointCloud3.load_ply(path)

    assert after.points == pytest.approx(before.points)
    assert after.point_normals == pytest.approx(before.point_normals)
    assert numpy.array_equal(after.point_colors, before.point_colors)
    assert after.point_stdev == pytest.approx(before.point_stdev)


# ================================================================================================
# Spatial queries
#
# These are the operations backed by the cached k-d tree.
# ================================================================================================


def test_spatial_queries_do_not_disturb_attributes():
    """Indexing no longer copies the cloud into a second type, so nothing can be lost on the way to
    a spatial query. This replaced a lossy `to_cloud`/`from_cloud` round trip."""
    cloud = loaded_cloud_data()

    kept = cloud.sample_poisson_disk(0.0001)
    assert len(kept) == 3

    sampled = cloud.extract_poisson_sample(0.0001)
    assert sampled.points == pytest.approx(cloud.points)
    assert sampled.point_stdev == pytest.approx(cloud.point_stdev)
    assert numpy.array_equal(sampled.point_colors, cloud.point_colors)


def test_spatial_queries_see_mutations():
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


def test_overlap_accepts_the_same_cloud_twice():
    """The cached tree is built through a shared borrow specifically so that passing one cloud as
    both arguments works rather than raising a borrow error."""
    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    assert cloud.overlap_points_by_reciprocity(cloud, 0.001) == [0, 1]


def test_estimate_normals_recovers_a_plane():
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


def test_estimate_normals_rejects_a_mismatched_reference():
    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    with pytest.raises(ValueError):
        cloud.estimate_normals(numpy.tile([0.0, 0.0, 1.0], (3, 1)), 2.0)


# ================================================================================================
# Voxel reduction
# ================================================================================================


def test_reduce_by_voxel_averages_and_reports_its_work():
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


def test_reduce_by_voxel_reports_coherence_when_normals_are_present():
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


def test_reduce_by_voxel_rejects_a_nonsense_size():
    cloud = PointCloud3(numpy.array([[0.0, 0.0, 0.0]]))
    with pytest.raises(ValueError):
        cloud.reduce_by_voxel(0.0)


# ================================================================================================
# Subsets, size and copying
# ================================================================================================


def test_extract_subset_points_takes_a_mask():
    cloud = loaded_cloud_data()
    mask = IndexMask.from_indices([0, 2], 3)

    sub = cloud.extract_subset_points(mask)

    assert len(sub) == 2
    assert sub.points == pytest.approx(cloud.points[[0, 2]])
    assert sub.point_stdev == pytest.approx([0.001, 0.003])


def test_subset_indices_carries_the_attributes():
    sub = loaded_cloud_data().extract_subset_indices([2, 0])

    assert len(sub) == 2
    assert sub.points[0] == pytest.approx([0.0, 1.0, 0.0])
    assert sub.point_stdev == pytest.approx([0.003, 0.001])
    assert numpy.array_equal(sub.point_colors, [[0, 0, 255], [255, 0, 0]])


def test_reports_its_size_and_bounds():
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


def test_append_in_place_unions_the_attributes():
    cloud = loaded_cloud_data()
    cloud.append_in_place(loaded_cloud_data())

    assert len(cloud) == 6
    assert cloud.point_stdev.shape == (6,)
    assert cloud.point_colors.shape == (6, 3)


def test_append_in_place_rejects_a_mismatch_without_modifying_the_target():
    cloud = loaded_cloud_data()
    other = loaded_cloud_data()
    other.set_point_stdev(None)

    with pytest.raises(ValueError):
        cloud.append_in_place(other)

    assert len(cloud) == 3
    assert cloud.point_stdev == pytest.approx([0.001, 0.002, 0.003])


def test_cloned_is_independent():
    original = loaded_cloud_data()
    copy = original.cloned()

    copy.set_point_stdev(numpy.array([9.0, 9.0, 9.0]))

    assert original.point_stdev == pytest.approx([0.001, 0.002, 0.003])
    assert copy.point_stdev == pytest.approx([9.0, 9.0, 9.0])


def test_cloud_point_flat_round_trips_and_clears():
    cloud = PointCloud3(triangle_points())
    assert cloud.point_flat is None

    cloud.set_point_flat(numpy.array([[0.0, 0.0], [1.5, 0.0], [0.0, 1.5]]))
    assert cloud.point_flat.shape == (3, 2)
    assert cloud.point_flat[2] == pytest.approx([0.0, 1.5])

    with pytest.raises(ValueError):
        cloud.set_point_flat(numpy.zeros((4, 2)))

    cloud.set_point_flat(None)
    assert cloud.point_flat is None
