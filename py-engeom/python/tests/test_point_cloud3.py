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

from pathlib import Path

import numpy
import pytest

from engeom.common import IndexMask
from engeom.geom3 import (
    Iso3,
    Mesh3,
    MeshData3,
    PatchFilter,
    Point3,
    PointCloud3,
    RepairOpts,
)


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

    assert cloud.point_count == 2
    assert len(cloud) == 2
    assert not cloud.is_empty

    aabb = cloud.compute_aabb()
    assert aabb.min.x == pytest.approx(0.0)
    assert aabb.max.z == pytest.approx(3.0)

    blank = PointCloud3.empty()
    assert blank.is_empty
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


# ================================================================================================
# PCD loading
# ================================================================================================

BUNNY_PCD = Path(__file__).parents[3] / "engeom" / "tests" / "data" / "bun_zipper_res4.pcd"


def write_pcd_with_a_nan_normal(path: Path) -> Path:
    """
    Write three points. The middle point has a finite position and the NaN normal that PCL's normal estimation
    writes for a point with too few neighbors.
    """
    path.write_text(
        "FIELDS x y z normal_x normal_y normal_z\n"
        "SIZE 4 4 4 4 4 4\n"
        "TYPE F F F F F F\n"
        "WIDTH 3\n"
        "DATA ascii\n"
        "0 0 0 0 0 1\n"
        "1 0 0 nan nan nan\n"
        "0 1 0 0 1 0\n"
    )
    return path


def test_load_pcd_reads_a_binary_compressed_file():
    cloud = PointCloud3.load_pcd(BUNNY_PCD)

    # The file holds 455 points, two of which have NaN positions and are dropped.
    assert cloud.points.shape == (453, 3)
    assert cloud.points[0] == pytest.approx([-0.0312216, 0.126304, 0.00514924], abs=1e-6)

    # The fixture's colors are derived from each point's index, so any misalignment shows up here.
    colors = cloud.point_colors
    assert colors.dtype == numpy.uint8
    i = numpy.arange(453)
    assert numpy.array_equal(colors, numpy.stack([i % 256, (7 * i) % 256, (13 * i) % 256], axis=1))

    normals = cloud.point_normals
    assert normals.shape == (453, 3)
    assert numpy.linalg.norm(normals, axis=1) == pytest.approx(numpy.ones(453))


def test_load_pcd_refuses_an_invalid_normal_by_default(tmp_path):
    path = write_pcd_with_a_nan_normal(tmp_path / "nan-normal.pcd")

    with pytest.raises(OSError, match="not a valid direction"):
        PointCloud3.load_pcd(path)


def test_load_pcd_can_drop_the_points_with_invalid_normals(tmp_path):
    path = write_pcd_with_a_nan_normal(tmp_path / "nan-normal.pcd")

    cloud = PointCloud3.load_pcd(path, invalid_normals="drop_points")

    assert numpy.array_equal(cloud.points, [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert numpy.array_equal(cloud.point_normals, [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])


def test_load_pcd_can_drop_the_normals(tmp_path):
    path = write_pcd_with_a_nan_normal(tmp_path / "nan-normal.pcd")

    cloud = PointCloud3.load_pcd(path, invalid_normals="drop_normals")

    assert cloud.points.shape == (3, 3)
    assert cloud.point_normals is None


def test_load_pcd_checks_the_invalid_normals_keyword(tmp_path):
    path = write_pcd_with_a_nan_normal(tmp_path / "nan-normal.pcd")

    with pytest.raises(ValueError, match="drop_points"):
        PointCloud3.load_pcd(path, invalid_normals="skip")

    # The choice is keyword-only, so the method refuses a positional string.
    with pytest.raises(TypeError):
        PointCloud3.load_pcd(path, "drop_points")

# ================================================================================================
# Normal orientation
# ================================================================================================


def _sphere_cloud(radius: float = 5.0, spacing: float = 0.3) -> PointCloud3:
    return Mesh3.create_sphere(radius, radius * 0.002).sample_poisson(spacing)


def test_estimate_normals_from_a_viewpoint_faces_it():
    """A point cannot be seen from behind, so every normal must face the viewpoint."""
    cloud = _sphere_cloud()
    viewpoint = numpy.array([0.0, 0.0, 1000.0])

    normals, confidence = cloud.estimate_normals(radius=0.75, viewpoint=viewpoint)

    assert normals.shape == (cloud.point_count, 3)
    assert confidence.shape == (cloud.point_count,)

    toward = viewpoint - cloud.points
    assert numpy.all(numpy.sum(normals * toward, axis=1) >= 0.0)


def test_estimate_normals_from_per_point_viewpoints_can_orient_a_sphere():
    """One external viewpoint per point orients a closed shape outward where one viewpoint cannot."""
    cloud = _sphere_cloud()
    points = cloud.points

    normals, _ = cloud.estimate_normals(radius=0.75, viewpoints=points * 2.0)

    # Outward means agreeing with the direction from the middle of the sphere.
    assert numpy.all(numpy.sum(normals * points, axis=1) > 0.0)


def test_estimate_normals_by_propagation_agrees_with_the_sampled_normals():
    cloud = _sphere_cloud()
    truth = cloud.point_normals

    normals, _ = cloud.estimate_normals(radius=0.75, propagate_k=12)

    agreement = numpy.mean(numpy.sum(normals * truth, axis=1) > 0.0)
    # Propagation seeds the topmost point with an upward normal, which points outward on a sphere.
    # The recovered set should therefore match the sampled normals or reverse all of them.
    assert agreement > 0.99 or agreement < 0.01


def test_estimate_normals_still_takes_must_match_and_radius_positionally():
    """The original call signature predates the orientation arguments and must remain compatible."""
    cloud = _sphere_cloud(5.0, 0.4)
    outward = cloud.points

    normals, confidence = cloud.estimate_normals(outward, 1.0)

    assert normals.shape == (cloud.point_count, 3)
    assert confidence.shape == (cloud.point_count,)
    assert numpy.all(numpy.sum(normals * outward, axis=1) >= 0.0)


def test_estimate_normals_requires_a_radius():
    """`radius` has a default only to make `must_match` optional; omitting it is an error."""
    cloud = _sphere_cloud(2.0, 0.5)

    with pytest.raises(ValueError, match="radius"):
        cloud.estimate_normals(viewpoint=[0.0, 0.0, 10.0])

    with pytest.raises(ValueError, match="radius"):
        cloud.estimate_normals(numpy.tile([0.0, 0.0, 1.0], (cloud.point_count, 1)))


def test_estimate_normals_needs_exactly_one_direction_source():
    cloud = _sphere_cloud(2.0, 0.5)
    up = numpy.tile([0.0, 0.0, 1.0], (cloud.point_count, 1))

    # None given.
    with pytest.raises(ValueError):
        cloud.estimate_normals(radius=0.75)

    # Two given.
    with pytest.raises(ValueError):
        cloud.estimate_normals(radius=0.75, viewpoint=[0.0, 0.0, 10.0], propagate_k=12)

    with pytest.raises(ValueError):
        cloud.estimate_normals(up, 0.75, viewpoint=[0.0, 0.0, 10.0])

    # A neighbor count of zero links a point to nothing.
    with pytest.raises(ValueError):
        cloud.estimate_normals(radius=0.75, propagate_k=0)


def test_estimate_point_spacing_recovers_the_sampling_radius():
    for spacing in (0.2, 0.4):
        cloud = _sphere_cloud(5.0, spacing)
        measured = cloud.estimate_point_spacing()
        assert spacing * 0.98 <= measured < spacing * 1.5

    assert PointCloud3.empty().estimate_point_spacing() == 0.0


# ================================================================================================
# Surface reconstruction
# ================================================================================================


def test_reconstruct_surface_recovers_a_sphere():
    """A smooth, closed sphere has no features for the field to round."""
    radius = 5.0
    spacing = 0.25
    source = Mesh3.create_sphere(radius, radius * 0.002)
    cloud = source.sample_poisson(spacing)

    mesh, report = cloud.reconstruct_surface(spacing)

    assert len(mesh.faces) > 10_000
    assert report.faces == len(mesh.faces)
    assert report.points_used == cloud.point_count
    assert report.points_dropped == 0
    assert report.known_voxels > 0
    assert report.cells_visited > 0
    assert report.repair is not None

    # No orientation pass was requested, so the report contains no orientation results.
    assert report.components is None
    assert report.normals_flipped is None

    # Every vertex must lie within one-tenth of a voxel of the true sphere.
    error = numpy.abs(numpy.linalg.norm(mesh.points, axis=1) - radius)
    assert error.max() < 0.1 * spacing


def test_reconstruct_surface_estimates_and_propagates_normals():
    radius = 5.0
    spacing = 0.25
    cloud = Mesh3.create_sphere(radius, radius * 0.002).sample_poisson(spacing)
    cloud.set_point_normals(None)

    mesh, report = cloud.reconstruct_surface(
        spacing, normal_radius=spacing * 2.5, propagate_k=12
    )

    assert len(mesh.faces) > 10_000
    assert report.components == 1
    assert report.normals_flipped is not None

    error = numpy.abs(numpy.linalg.norm(mesh.points, axis=1) - radius)
    assert error.max() < 0.15 * spacing


def test_reconstruct_surface_estimates_from_a_viewpoint():
    spacing = 0.25
    cloud = Mesh3.create_sphere(5.0, 0.01).sample_poisson(spacing)
    cloud.set_point_normals(None)

    mesh, report = cloud.reconstruct_surface(
        spacing, normal_radius=spacing * 2.5, viewpoint=[0.0, 0.0, 1000.0]
    )

    assert len(mesh.faces) > 1000
    # Viewpoint orientation does not use a graph, so there are no components to report.
    assert report.components is None


def test_reconstruct_surface_needs_normals_from_somewhere():
    spacing = 0.3
    cloud = Mesh3.create_sphere(2.0, 0.01).sample_poisson(spacing)
    cloud.set_point_normals(None)

    with pytest.raises(ValueError, match="no normals"):
        cloud.reconstruct_surface(spacing)


def test_reconstruct_surface_needs_exactly_one_direction_source():
    spacing = 0.3
    cloud = Mesh3.create_sphere(2.0, 0.01).sample_poisson(spacing)

    # Estimation requires an orientation source.
    with pytest.raises(ValueError):
        cloud.reconstruct_surface(spacing, normal_radius=0.75)

    # Estimation accepts only one orientation source.
    with pytest.raises(ValueError):
        cloud.reconstruct_surface(
            spacing, normal_radius=0.75, viewpoint=[0.0, 0.0, 10.0], propagate_k=12
        )

    # An orientation source cannot apply when normal estimation is disabled.
    with pytest.raises(ValueError):
        cloud.reconstruct_surface(spacing, propagate_k=12)


def test_reconstruct_surface_rejects_options_that_cannot_work():
    spacing = 0.3
    cloud = Mesh3.create_sphere(2.0, 0.01).sample_poisson(spacing)

    for voxel_size in (0.0, -1.0):
        with pytest.raises(ValueError):
            cloud.reconstruct_surface(voxel_size)

    # A band thinner than the diagonal of a cell cannot hold a crossing with all eight corners.
    for band_width in (1.0, 1.7):
        with pytest.raises(ValueError):
            cloud.reconstruct_surface(spacing, band_width=band_width)

    # The smallest band that can work.
    mesh, _ = cloud.reconstruct_surface(spacing, band_width=1.75)
    assert len(mesh.faces) > 0

    with pytest.raises(ValueError):
        cloud.reconstruct_surface(spacing, min_normal_confidence=-0.1)

    with pytest.raises(ValueError):
        PointCloud3.empty().reconstruct_surface(spacing)


def test_reconstruct_surface_turns_inside_out_with_reversed_normals():
    """Downstream processing does not correct reversed normals."""
    spacing = 0.25
    cloud = Mesh3.create_sphere(5.0, 0.01).sample_poisson(spacing)
    cloud.set_point_normals(-cloud.point_normals)

    mesh, _ = cloud.reconstruct_surface(spacing)

    # All face normals should point toward the center.
    tri = mesh.points[mesh.faces]
    normals = numpy.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    centroids = tri.mean(axis=1)
    assert numpy.all(numpy.sum(normals * centroids, axis=1) < 0.0)


def test_reconstruct_surface_defaults_to_the_cheaper_repair():
    """Reconstruction omits two unnecessary repair passes. The default and full repair sets must
    produce meshes with the same size, providing evidence for that default."""
    spacing = 0.25
    cloud = Mesh3.create_sphere(5.0, 0.01).sample_poisson(spacing)

    default_mesh, _ = cloud.reconstruct_surface(spacing)
    full_mesh, _ = cloud.reconstruct_surface(spacing, repair=RepairOpts())
    named_mesh, _ = cloud.reconstruct_surface(
        spacing, repair=RepairOpts.assuming_oriented_edges()
    )

    assert len(default_mesh.faces) == len(full_mesh.faces)
    assert len(default_mesh.points) == len(full_mesh.points)
    assert len(named_mesh.faces) == len(default_mesh.faces)


def test_repair_opts_assuming_oriented_edges_drops_only_the_edge_passes():
    cheap = RepairOpts.assuming_oriented_edges()

    assert not cheap.resolve_nonmanifold_edges
    assert not cheap.orient_consistently

    full = RepairOpts()
    assert cheap.drop_degenerate == full.drop_degenerate
    assert cheap.drop_duplicate_faces == full.drop_duplicate_faces
    assert cheap.split_bowtie_vertices == full.split_bowtie_vertices
    assert cheap.drop_isolated_vertices == full.drop_isolated_vertices


def test_reconstruct_surface_accepts_repair_and_patch_options():
    spacing = 0.25
    cloud = Mesh3.create_sphere(5.0, 0.01).sample_poisson(spacing)

    # Python skips repair by disabling every repair pass.
    mesh, report = cloud.reconstruct_surface(spacing, repair=RepairOpts.none())
    assert len(mesh.faces) == report.raw_faces

    filtered, _ = cloud.reconstruct_surface(
        spacing, patch_filter=PatchFilter.keep_largest()
    )
    assert len(filtered.faces) > 0

    smoothed, _ = cloud.reconstruct_surface(spacing, smooth_iterations=2)
    assert len(smoothed.faces) > 0


def test_reconstruct_surface_drops_points_below_the_confidence_threshold():
    spacing = 0.25
    cloud = Mesh3.create_sphere(5.0, 0.01).sample_poisson(spacing)

    # Each point is far from the sphere and the other added points. It is therefore alone within the
    # estimation radius and cannot support a plane fit.
    strays = numpy.array([[40.0 + i * 10.0, 0.0, 0.0] for i in range(5)])
    widened = PointCloud3(numpy.vstack([cloud.points, strays]))

    _, report = widened.reconstruct_surface(
        spacing,
        normal_radius=spacing * 2.5,
        viewpoint=[0.0, 0.0, 1000.0],
        min_normal_confidence=0.01,
    )

    assert report.points_dropped == len(strays)
    assert report.points_used == cloud.point_count

def test_reconstruct_surface_returns_buffers_the_caller_converts():
    """The pipeline returns buffers because it does not query a bounding volume hierarchy."""
    spacing = 0.25
    cloud = Mesh3.create_sphere(5.0, 0.01).sample_poisson(spacing)

    data, _ = cloud.reconstruct_surface(spacing)

    assert isinstance(data, MeshData3)
    assert len(data.faces) > 1000

    # The caller converts the buffers before performing spatial queries.
    mesh = data.to_mesh(is_solid=False)
    assert isinstance(mesh, Mesh3)
    assert len(mesh.faces) == len(data.faces)
    assert mesh.point_closest_to(100.0, 0.0, 0.0) is not None
