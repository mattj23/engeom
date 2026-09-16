"""
Tests for the row-organized point format (`.tcrpf3`) bindings.

The format preserves the grouping produced by a rasterizing sensor. This grouping permits meshing
by triangulating neighboring rows without surface reconstruction. These tests verify that the
grouping survives. The Rust crate tests coordinate round-tripping.
"""

from __future__ import annotations

import numpy
import pytest

from engeom.geom3 import Mesh3, MeshData3, PointCloud3, RowPointsScan3

ROW_PITCH = 0.1
COL_PITCH = 0.05
TOL = 1e-6


def snapshot_rows(rows: int, cols: int, slope: float = 0.3):
    """
    Create a synthetic scan with the row geometry produced by a snapshot sensor.

    A raster row is a plane through the camera instead of a constant-y line, so its y-coordinate
    varies with depth. Here it varies by several millimeters while adjacent rows are one-tenth of a
    millimeter apart, matching the Gocator data.
    """
    out = []
    for i in range(rows):
        x = numpy.arange(cols, dtype=numpy.float64) * COL_PITCH
        z = 10.0 + 0.5 * (x - 1.0) ** 2
        y = i * ROW_PITCH + slope * z
        out.append(numpy.column_stack([x, y, z]))
    return out


def make_scan(rows: int = 20, cols: int = 50, slope: float = 0.3, **kwargs) -> RowPointsScan3:
    data = snapshot_rows(rows, cols, slope)
    ordinals = numpy.arange(rows, dtype=numpy.uint32)
    return RowPointsScan3.from_rows(data, ordinals, ROW_PITCH, col_pitch=COL_PITCH, **kwargs)


def test_round_trips_through_a_file(tmp_path):
    scan = make_scan(name="pass 1")
    path = tmp_path / "scan.tcrpf3"
    scan.write(path, TOL)

    back = RowPointsScan3.read(path)
    assert back.row_count == scan.row_count
    assert back.point_count == scan.point_count == 20 * 50
    assert back.row_pitch == pytest.approx(ROW_PITCH)
    assert back.col_pitch == pytest.approx(COL_PITCH)
    assert back.along == "x"
    assert back.name == "pass 1"
    numpy.testing.assert_array_equal(back.ordinals, scan.ordinals)

    for i in range(back.row_count):
        numpy.testing.assert_allclose(back.row(i), scan.row(i), atol=TOL)


def test_the_row_arrays_have_the_expected_shape_and_dtype():
    scan = make_scan(rows=4, cols=7)
    row = scan.row(2)
    assert row.shape == (7, 3)
    assert row.dtype == numpy.float64
    assert scan.ordinals.dtype == numpy.uint32


def test_columns_round_trip_when_given(tmp_path):
    rows, cols = 6, 9
    data = snapshot_rows(rows, cols)
    ordinals = numpy.arange(rows, dtype=numpy.uint32)
    columns = [numpy.arange(cols, dtype=numpy.uint32) for _ in range(rows)]

    scan = RowPointsScan3.from_rows(data, ordinals, ROW_PITCH, columns=columns)
    path = tmp_path / "columns.tcrpf3"
    scan.write(path, TOL)

    back = RowPointsScan3.read(path)
    for i in range(rows):
        numpy.testing.assert_array_equal(back.row_columns(i), columns[i])


def test_a_scan_without_columns_reports_none():
    scan = make_scan(rows=3, cols=4)
    assert scan.row_columns(0) is None


def test_metadata_round_trips(tmp_path):
    meta = {
        "engeom.source_frame": "gocator-3x00 left-handed, x negated",
        "exposure_us": 500,
        "captured": 1653436800,
        "calibrated": True,
        "temperature_c": 21.5,
    }
    scan = make_scan(rows=4, cols=6, metadata=meta)
    path = tmp_path / "meta.tcrpf3"
    scan.write(path, TOL)

    back = RowPointsScan3.read(path)
    assert back.metadata == meta
    # A bool must not come back as an int, which it would if `bool` were checked after `int`.
    assert back.metadata["calibrated"] is True


def test_ordinal_gaps_survive(tmp_path):
    data = snapshot_rows(4, 5)
    ordinals = numpy.array([0, 1, 40, 41], dtype=numpy.uint32)
    scan = RowPointsScan3.from_rows(data, ordinals, ROW_PITCH)

    path = tmp_path / "gaps.tcrpf3"
    scan.write(path, TOL)
    numpy.testing.assert_array_equal(RowPointsScan3.read(path).ordinals, ordinals)


def test_mismatched_ordinals_are_rejected():
    data = snapshot_rows(4, 5)
    with pytest.raises(ValueError):
        RowPointsScan3.from_rows(data, numpy.arange(3, dtype=numpy.uint32), ROW_PITCH)


def test_mismatched_column_lengths_are_rejected():
    data = snapshot_rows(2, 5)
    ordinals = numpy.arange(2, dtype=numpy.uint32)
    columns = [numpy.arange(5, dtype=numpy.uint32), numpy.arange(3, dtype=numpy.uint32)]
    with pytest.raises(ValueError):
        RowPointsScan3.from_rows(data, ordinals, ROW_PITCH, columns=columns)


def test_an_unknown_along_axis_is_rejected():
    data = snapshot_rows(2, 5)
    with pytest.raises(ValueError):
        RowPointsScan3.from_rows(
            data, numpy.arange(2, dtype=numpy.uint32), ROW_PITCH, along="z"
        )


def test_non_increasing_ordinals_are_refused_on_write(tmp_path):
    data = snapshot_rows(3, 5)
    ordinals = numpy.array([0, 2, 2], dtype=numpy.uint32)
    scan = RowPointsScan3.from_rows(data, ordinals, ROW_PITCH)

    with pytest.raises(OSError):
        scan.write(tmp_path / "bad.tcrpf3", TOL)


# ================================================================================================
# Loading
# ================================================================================================


def test_loads_as_a_point_cloud(tmp_path):
    scan = make_scan(rows=10, cols=25)
    path = tmp_path / "cloud.tcrpf3"
    scan.write(path, TOL)

    cloud = PointCloud3.load_tc_row_points(path)
    assert cloud.points.shape == (250, 3)


def test_meshes_a_scan_whose_rows_are_not_constant_y(tmp_path):
    """
    A row's strip coordinate comes from its ordinal instead of one of its point coordinates. A row
    whose y-coordinate varies by 3 mm must still mesh with its neighbor 0.1 mm away.
    """
    scan = make_scan(rows=24, cols=60, slope=0.3)
    path = tmp_path / "mesh.tcrpf3"
    scan.write(path, TOL)

    mesh = Mesh3.load_tc_row_points(path)
    assert mesh.faces.shape == (23 * 59 * 2, 3)
    assert mesh.points.shape == (24 * 60, 3)

    # Every face points back toward the sensor, which is what a right-handed frame with +z toward
    # the sensor and ordinals ascending in the sweep direction is supposed to give.
    p = mesh.points[mesh.faces]
    normals = numpy.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    assert (normals[:, 2] > 0).all()


def test_mesh_data_loads_the_same_geometry(tmp_path):
    scan = make_scan(rows=12, cols=20)
    path = tmp_path / "data.tcrpf3"
    scan.write(path, TOL)

    mesh = Mesh3.load_tc_row_points(path)
    data = MeshData3.load_tc_row_points(path)
    numpy.testing.assert_allclose(data.points, mesh.points)
    numpy.testing.assert_array_equal(data.faces, mesh.faces)


def test_a_scan_along_y_meshes_the_same_way_up(tmp_path):
    """A column-major sensor's y-oriented strips flatten to a mirrored plane."""
    rows, cols = 20, 40
    data = []
    for i in range(rows):
        y = numpy.arange(cols, dtype=numpy.float64) * COL_PITCH
        z = 10.0 + 0.5 * (y - 1.0) ** 2
        x = numpy.full(cols, i * ROW_PITCH)
        data.append(numpy.column_stack([x, y, z]))

    scan = RowPointsScan3.from_rows(
        data, numpy.arange(rows, dtype=numpy.uint32), ROW_PITCH, along="y"
    )
    path = tmp_path / "along-y.tcrpf3"
    scan.write(path, TOL)

    mesh = Mesh3.load_tc_row_points(path)
    assert mesh.faces.shape == (19 * 39 * 2, 3)

    p = mesh.points[mesh.faces]
    normals = numpy.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    assert (normals[:, 2] > 0).all()


def test_thinning_drops_points_in_both_directions(tmp_path):
    scan = make_scan(rows=40, cols=100, slope=0.2)
    path = tmp_path / "thin.tcrpf3"
    scan.write(path, TOL)

    full = PointCloud3.load_tc_row_points(path)
    thin = PointCloud3.load_tc_row_points(path, take_every=4)

    assert 0 < len(thin.points) < len(full.points) / 20


def test_smoothing_moves_points_along_z_only(tmp_path):
    scan = make_scan(rows=40, cols=100, slope=0.2)
    path = tmp_path / "smooth.tcrpf3"
    scan.write(path, TOL)

    thin = PointCloud3.load_tc_row_points(path, take_every=4)
    smooth = PointCloud3.load_tc_row_points(
        path, take_every=4, look_scale=1.5, weight_scale=1.0, max_move=1.0
    )

    assert smooth.points.shape == thin.points.shape
    numpy.testing.assert_allclose(smooth.points[:, :2], thin.points[:, :2], atol=1e-12)


def test_partial_smoothing_parameters_are_rejected(tmp_path):
    scan = make_scan(rows=4, cols=6)
    path = tmp_path / "partial.tcrpf3"
    scan.write(path, TOL)

    with pytest.raises(ValueError):
        PointCloud3.load_tc_row_points(path, take_every=4, look_scale=1.5)


def test_a_file_that_is_not_a_scan_is_rejected(tmp_path):
    path = tmp_path / "junk.tcrpf3"
    path.write_bytes(b"NOPE0000\x00\x00")

    with pytest.raises(OSError):
        RowPointsScan3.read(path)


def test_a_capture_that_saw_nothing_is_a_valid_empty_scan(tmp_path):
    """
    A sensor can validly return no measurements, so an empty scan must remain writable and readable.
    It has no faces, and requesting `Mesh3` must return an error without panicking.
    """
    scan = RowPointsScan3.from_rows([], numpy.zeros(0, dtype=numpy.uint32), ROW_PITCH, name="none")
    path = tmp_path / "empty.tcrpf3"
    scan.write(path, TOL)

    back = RowPointsScan3.read(path)
    assert back.row_count == 0
    assert back.point_count == 0
    assert back.name == "none"

    assert PointCloud3.load_tc_row_points(path).points.shape == (0, 3)
    assert MeshData3.load_tc_row_points(path).faces.shape == (0, 3)

    with pytest.raises(ValueError):
        Mesh3.load_tc_row_points(path)
