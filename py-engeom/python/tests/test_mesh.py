import numpy
import pytest
from numpy import linalg
from engeom.geom3 import Mesh3, Iso3


def test_mesh_offset_points_copy():
    m = Mesh3.create_sphere(1.0, 100, 100)
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
    m = Mesh3.create_sphere(1.0, 20, 20)
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


def test_mesh_compute_patches_finds_disconnected_pieces():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    other = Mesh3.create_box(1.0, 1.0, 1.0)
    other.transform_in_place(Iso3.from_translation(10.0, 0.0, 0.0))
    m.append_in_place(other)

    patches = m.compute_patches()

    assert len(patches) == 2
    assert all(len(p) == m.face_count for p in patches)
    assert sum(p.count_true() for p in patches) == m.face_count


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

    handle = m.face_select_none().by_mask(mask, "add")

    assert handle.collect_indices() == [3, 4]
    assert handle.to_mask() == mask
    assert handle.to_mesh().face_count == 2


def test_face_filter_handle_to_mask_does_not_consume_the_handle():
    m = Mesh3.create_box(1.0, 1.0, 1.0)
    handle = m.face_select_all()

    first = handle.to_mask()

    assert first.all()
    assert handle.collect_indices() == list(range(m.face_count))
