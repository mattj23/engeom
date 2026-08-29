import numpy
import pytest
from numpy import linalg
from engeom.geom3 import Mesh3, Iso3, Plane3, PatchFilter
from engeom.common import IndexMask


def test_mesh_offset_points_copy():
    m = Mesh3.create_sphere(1.0, 1.0e-3)
    n = m.offset_points_copy(0.1)

    for v in n.points:
        assert pytest.approx(1.1, abs=1e-5) == linalg.norm(v)


def test_mesh_compute_point_normals_of_a_box_point_along_the_diagonals():
    """
    Angle weighting is invariant to how a flat face is triangulated, so every corner normal of a box
    lands on the diagonal. Area weighting, which this used to do, pulls them off it.
    """
    m = Mesh3.create_box(2.0, 2.0, 2.0)
    normals = m.compute_point_normals()
    diagonal = 1.0 / numpy.sqrt(3.0)

    assert normals.shape == (m.points.shape[0], 3)
    for point, normal in zip(m.points, normals):
        expected = numpy.sign(point) * diagonal
        assert numpy.allclose(normal, expected, atol=1e-12)


def test_mesh_compute_point_normals_is_cached():
    m = Mesh3.create_sphere(1.0, 0.025)
    assert m.compute_point_normals() is m.compute_point_normals()


def test_mesh_counts_and_len():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    assert m.point_count == m.points.shape[0]
    assert m.face_count == m.faces.shape[0]
    assert len(m) == m.point_count


def test_mesh_is_solid_follows_the_constructor():
    assert Mesh3.create_box(1.0, 1.0, 1.0).is_solid
    assert not Mesh3.create_box(1.0, 1.0, 1.0, is_solid=False).is_solid


def test_mesh_new_accepts_is_solid():
    points = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    faces = numpy.array([[0, 1, 2]], dtype=numpy.uint32)

    assert not Mesh3(points, faces).is_solid
    assert Mesh3(points, faces, is_solid=True).is_solid


def _single_triangle() -> Mesh3:
    points = numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    faces = numpy.array([[0, 1, 2]], dtype=numpy.uint32)
    return Mesh3(points, faces)


def test_mesh_attributes_are_none_when_absent():
    m = _single_triangle()
    assert m.point_normals is None
    assert m.point_colors is None
    assert m.point_stdev is None
    assert m.face_colors is None
    assert m.face_labels is None


def test_mesh_point_attributes_round_trip():
    m = _single_triangle()

    m.set_point_normals(numpy.array([[0.0, 0.0, 2.0]] * 3))
    m.set_point_colors(numpy.array([[10, 20, 30]] * 3, dtype=numpy.uint8))
    m.set_point_stdev(numpy.array([0.1, 0.2, 0.3]))

    # Normals are normalized on the way in, so the stored value is the unit direction.
    assert numpy.allclose(m.point_normals, numpy.array([[0.0, 0.0, 1.0]] * 3))
    assert m.point_colors.dtype == numpy.uint8
    assert numpy.array_equal(m.point_colors, numpy.array([[10, 20, 30]] * 3, dtype=numpy.uint8))
    assert numpy.allclose(m.point_stdev, [0.1, 0.2, 0.3])


def test_mesh_face_attributes_round_trip():
    m = _single_triangle()

    m.set_face_colors(numpy.array([[1, 2, 3]], dtype=numpy.uint8))
    m.set_face_labels(numpy.array([7], dtype=numpy.uint32))

    assert numpy.array_equal(m.face_colors, numpy.array([[1, 2, 3]], dtype=numpy.uint8))
    assert m.face_labels.dtype == numpy.uint32
    assert numpy.array_equal(m.face_labels, [7])


def test_mesh_attributes_can_be_cleared():
    m = _single_triangle()
    m.set_face_labels(numpy.array([7], dtype=numpy.uint32))
    assert m.face_labels is not None

    m.set_face_labels(None)
    assert m.face_labels is None


def test_mesh_attribute_setters_reject_a_length_mismatch():
    m = _single_triangle()
    with pytest.raises(ValueError):
        m.set_point_stdev(numpy.array([0.1, 0.2]))


def test_mesh_setting_an_attribute_invalidates_the_cached_arrays():
    """A stale cached array would keep handing back the old attribute after a setter runs."""
    m = _single_triangle()
    m.set_face_labels(numpy.array([7], dtype=numpy.uint32))
    first = m.face_labels

    m.set_face_labels(numpy.array([9], dtype=numpy.uint32))

    assert numpy.array_equal(m.face_labels, [9])
    assert first is not m.face_labels


def test_mesh_flip_normals_in_place_negates_stored_normals():
    m = _single_triangle()
    m.set_point_normals(numpy.array([[0.0, 0.0, 1.0]] * 3))

    m.flip_normals_in_place()

    assert numpy.allclose(m.point_normals, numpy.array([[0.0, 0.0, -1.0]] * 3))


def test_mesh_face_and_point_masks_are_sized_to_the_mesh():
    m = Mesh3.create_box(1.0, 1.0, 1.0)

    assert len(m.face_mask()) == m.face_count
    assert len(m.point_mask()) == m.point_count
    assert not m.face_mask().any()
    assert m.face_mask(True).all()


def test_mesh_extract_subset_faces_by_mask():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    mask = m.face_mask()
    mask[0] = True
    mask[1] = True

    sub = m.extract_subset_faces(mask)

    assert sub.face_count == 2
    # Only the points those two faces reference survive.
    assert sub.point_count == len(m.compute_unique_point_mask(mask).to_indices())


def test_mesh_extract_subset_faces_rejects_a_mask_of_the_wrong_length():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        m.extract_subset_faces(m.point_mask(True))


def test_mesh_extract_subset_faces_matches_the_index_route():
    """The mask and index routes into a subset are the same operation, so they must agree."""
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    mask = m.face_mask()
    for i in (2, 5, 7):
        mask[i] = True

    by_mask = m.extract_subset_faces(mask)
    by_indices = m.extract_subset_faces_from_indices([2, 5, 7])

    assert numpy.array_equal(by_mask.points, by_indices.points)
    assert numpy.array_equal(by_mask.faces, by_indices.faces)


def test_mesh_compute_patch_labels_finds_disconnected_pieces():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    other = Mesh3.create_box(1.0, 1.0, 1.0)
    other.transform_in_place(Iso3.from_translation(10.0, 0.0, 0.0))
    m.append_in_place(other)

    labels = m.compute_patch_labels()

    assert labels.dtype == numpy.uint32
    assert len(labels) == m.face_count
    assert set(numpy.unique(labels)) == {0, 1}

    # Both boxes tessellate the same way, so the faces split evenly between the two patches.
    counts = numpy.bincount(labels)
    assert list(counts) == [m.face_count // 2, m.face_count // 2]

    # Patches are numbered by their lowest-indexed face, so the box appended first is patch 0.
    assert labels[0] == 0
    assert labels[-1] == 1


def test_mesh_compute_patch_labels_marks_masked_faces():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    other = Mesh3.create_box(1.0, 1.0, 1.0)
    other.transform_in_place(Iso3.from_translation(10.0, 0.0, 0.0))
    m.append_in_place(other)

    # Keep only the first box, which is the first half of the faces.
    half = m.face_count // 2
    mask = IndexMask.from_indices(range(half), m.face_count)

    labels = m.compute_patch_labels(mask)

    no_patch = 2**32 - 1
    assert set(numpy.unique(labels)) == {0, no_patch}
    assert numpy.all(labels[: m.face_count // 2] == 0)
    assert numpy.all(labels[m.face_count // 2 :] == no_patch)


def _body_with_flyer() -> Mesh3:
    """A body with a small flyer parked well away from it, as scan data tends to produce."""
    m = Mesh3.create_box(20.0, 20.0, 20.0)
    flyer = Mesh3.create_box(0.5, 0.5, 0.5)
    flyer.transform_in_place(Iso3.from_translation(300.0, 0.0, 0.0))
    m.append_in_place(flyer)
    return m


def test_mesh_remove_small_patches_drops_the_flyer():
    m = _body_with_flyer()
    assert m.compute_patch_labels().max() == 1

    cleaned = m.remove_small_patches(PatchFilter.keep_largest())

    assert cleaned.face_count == 12
    assert cleaned.compute_patch_labels().max() == 0

    # With the flyer gone the bounding box is the body alone.
    lo, hi = cleaned.aabb.min, cleaned.aabb.max
    assert hi.x - lo.x == pytest.approx(20.0, abs=1e-9)


@pytest.mark.parametrize(
    "filter_",
    [
        PatchFilter(min_faces=None, min_area=10.0),
        PatchFilter(min_aabb_diagonal=5.0),
        PatchFilter(min_area_fraction=0.01),
        PatchFilter(keep_largest_n=1),
    ],
)
def test_mesh_remove_small_patches_criteria_each_catch_the_flyer(filter_):
    cleaned = _body_with_flyer().remove_small_patches(filter_)
    assert cleaned.face_count == 12


def test_mesh_patch_mask_reports_what_would_be_kept():
    m = _body_with_flyer()
    mask = m.patch_mask(PatchFilter.keep_largest())

    assert len(mask) == m.face_count
    assert mask.count_true() == 12
    kept = mask.to_indices()
    assert list(kept) == list(range(12))


def test_mesh_remove_small_patches_keeping_everything_is_a_no_op():
    m = _body_with_flyer()
    same = m.remove_small_patches(PatchFilter())

    assert same.face_count == m.face_count
    assert same.point_count == m.point_count


def test_mesh_remove_small_patches_rejects_discarding_everything():
    m = _body_with_flyer()
    with pytest.raises(ValueError):
        m.remove_small_patches(PatchFilter(min_faces=1_000_000))


def test_patch_filter_exposes_its_criteria():
    f = PatchFilter(min_faces=50, min_area_fraction=0.01)

    assert f.min_faces == 50
    assert f.min_area_fraction == 0.01
    assert f.min_area is None
    assert f.min_aabb_diagonal is None
    assert f.keep_largest_n is None
    assert PatchFilter.keep_largest().keep_largest_n == 1


def test_face_select_keep_patches():
    m = _body_with_flyer()
    kept = m.face_select("all").keep_patches(PatchFilter.keep_largest(), "keep").to_mesh()

    assert kept.face_count == 12


def test_mesh_find_points_in_tol_is_indexed_by_the_given_points():
    m = Mesh3.create_box(2.0, 2.0, 2.0)
    points = numpy.array([[0.0, 0.0, 1.0], [0.0, 0.0, 50.0]])

    mask = m.find_points_in_tol(points, 0.1, numpy.pi / 4)

    assert len(mask) == 2
    assert mask[0]
    assert not mask[1]


def test_face_filter_handle_mask_round_trip():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    mask = m.face_mask()
    mask[3] = True
    mask[4] = True

    handle = m.face_select().by_mask(mask, "add")

    assert handle.collect_indices() == [3, 4]
    assert handle.to_mask() == mask
    assert handle.to_mesh().face_count == 2


def test_face_filter_handle_to_mask_does_not_consume_the_handle():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    handle = m.face_select("all")

    first = handle.to_mask()

    assert first.all()
    assert handle.collect_indices() == list(range(m.face_count))


def test_filter_above_plane_splits_a_box():
    m = Mesh3.create_box(2.0, 2.0, 2.0)
    plane = Plane3(0.0, 0.0, 1.0, 0.0)

    any_above = m.face_select().above_plane(plane, False, "add").collect_indices()
    all_above = m.face_select().above_plane(plane, True, "add").collect_indices()

    # Requiring every vertex is strictly stronger than requiring one.
    assert set(all_above) <= set(any_above)
    assert len(all_above) < len(any_above)


def test_filter_vertices_near_point():
    m = Mesh3.create_box(2.0, 2.0, 2.0)
    corner = m.points[0]

    near = m.face_select().vertices_near_point(*corner, 0.1, False, "add").collect_indices()
    far = m.face_select().vertices_near_point(*corner, 100.0, True, "add").collect_indices()

    assert 0 < len(near) < m.face_count
    assert len(far) == m.face_count


def test_filter_expand_dilates_and_erodes():
    m = Mesh3.create_sphere(1.0, 0.025)
    seed = IndexMask.from_indices([0], m.face_count)

    grown = m.face_select().by_mask(seed, "add").expand("add").to_mask()
    assert grown.count_true() > 1
    assert (grown & seed) == seed  # dilation keeps what it started with

    # Erosion needs a border to eat into. A closed sphere has none when everything is selected, so
    # erode the dilated patch instead, which does.
    eroded = m.face_select().by_mask(grown, "add").expand("remove").to_mask()
    assert eroded.count_true() < grown.count_true()


def test_filter_expand_remove_on_a_closed_mesh_selection_is_a_no_op():
    """Erosion works from the unselected side, so a fully selected closed mesh has nothing to erode."""
    m = Mesh3.create_sphere(1.0, 0.025)

    eroded = m.face_select("all").expand("remove").collect_indices()

    assert len(eroded) == m.face_count


def test_filter_expand_n_matches_repeated_expand():
    m = Mesh3.create_sphere(1.0, 0.025)
    seed = IndexMask.from_indices([0], m.face_count)

    once_twice = (m.face_select().by_mask(seed, "add")
                  .expand("add").expand("add").collect_indices())
    n_two = (m.face_select().by_mask(seed, "add")
             .expand_n(2, "add").collect_indices())

    assert once_twice == n_two


def test_filter_expand_respects_the_exclude_mask():
    m = Mesh3.create_sphere(1.0, 0.025)
    seed = IndexMask.from_indices([0], m.face_count)

    free = m.face_select().by_mask(seed, "add").expand("add").to_mask()
    # Exclude everything the unrestricted expansion would have reached, minus the seed itself.
    blocked = m.face_select().by_mask(seed, "add").expand("add", exclude=free - seed).to_mask()

    assert blocked == seed


def test_filter_expand_rejects_a_mask_of_the_wrong_length():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        m.face_select("all").expand("add", exclude=IndexMask(3))


def test_filter_faces_overlap_finds_a_coincident_copy():
    m = Mesh3.create_box(2.0, 2.0, 2.0)
    same = Mesh3.create_box(2.0, 2.0, 2.0)
    apart = Mesh3.create_box(2.0, 2.0, 2.0)
    apart.transform_in_place(Iso3.from_translation(50.0, 0.0, 0.0))

    on_top = m.face_select().faces_overlap(same, 0.1, 0.1, "add").collect_indices()
    nowhere = m.face_select().faces_overlap(apart, 0.1, 0.1, "add").collect_indices()

    assert len(on_top) == m.face_count
    assert nowhere == []


def test_mesh_face_select_seeds():
    m = Mesh3.create_box(1.0, 1.0, 1.0)

    assert m.face_select().collect_indices() == []
    assert m.face_select("none").collect_indices() == []
    assert m.face_select("all").collect_indices() == list(range(m.face_count))


def test_mesh_face_select_accepts_a_mask():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    mask = IndexMask.from_indices([1, 4, 9], m.face_count)

    assert m.face_select(mask).collect_indices() == [1, 4, 9]


def test_mesh_face_select_rejects_a_bad_seed():
    m = Mesh3.create_box(1.0, 1.0, 1.0)

    with pytest.raises(ValueError):
        m.face_select("everything")
    # A mask sized to the points would otherwise be zipped against the faces silently.
    with pytest.raises(ValueError):
        m.face_select(m.point_mask(True))


def test_mesh_measure_deviations_signs_and_shape():
    m = Mesh3.create_box(2.0, 2.0, 2.0, is_solid=False)
    points = numpy.array([[0.0, 0.0, 2.0], [0.0, 0.0, 0.5]])

    values = m.measure_deviations(points, "point")

    assert values.shape == (2,)
    assert values[0] == pytest.approx(1.0, abs=1e-12)
    assert values[1] < 0.0


@pytest.mark.parametrize("mode", ["point", "plane"])
def test_mesh_measure_deviations_matches_the_single_point_measurement(mode):
    """The bulk and single-point paths share a direction helper in the core, so they must agree."""
    m = Mesh3.create_box(2.0, 2.0, 2.0)
    points = numpy.array([
        [0.0, 0.0, 3.0],
        [0.0, 0.0, -3.0],
        [5.0, 5.0, 5.0],
        [0.1, -0.2, 0.05],
    ])

    bulk = m.measure_deviations(points, mode)
    single = [m.measure_point_deviation(*p, mode).value for p in points]

    assert numpy.allclose(bulk, single, atol=1e-15)


def test_mesh_measure_deviations_rejects_an_unknown_mode():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    with pytest.raises(ValueError):
        m.measure_deviations(numpy.array([[0.0, 0.0, 5.0]]), "surface")


def test_mesh_distance_and_face_closest_to():
    m = Mesh3.create_box(2.0, 2.0, 2.0, is_solid=False)

    assert m.distance_closest_to(0.0, 0.0, 2.0) == pytest.approx(1.0, abs=1e-12)
    # The distance is unsigned, unlike measure_deviations.
    assert m.distance_closest_to(0.0, 0.0, 0.5) > 0.0

    face = m.face_closest_to(0.0, 0.0, 5.0)
    assert 0 <= face < m.face_count
    # That face has to be one of the two on the +Z side.
    assert numpy.allclose(m.compute_face_normals()[face], [0.0, 0.0, 1.0])
