"""
Tests for the half-edge mesh bindings: repair, both decimation paths, smoothing, and the round trip.

These check the binding layer rather than the algorithms, which are covered on the Rust side. What
matters here is that options reach the core intact, that reports come back readable, and that the
one refusal which protects a tolerance claim survives the boundary.
"""

import pytest

from engeom.geom3 import (
    BestEffortOpts,
    DecimateOpts,
    HalfEdgeMesh3,
    Mesh3,
    RepairOpts,
)

# Loose enough that the coarse bunny actually reduces, and the same on both paths so the two are
# comparable.
TOL = 0.01


@pytest.fixture(scope="module")
def bunny() -> Mesh3:
    return Mesh3.stanford_bunny_res4()


def test_repairing_constructor_reports_what_it_changed(bunny):
    """The bunny is non-manifold, so the repair has something to say and says it."""
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())

    assert he.face_count > 0
    assert he.vertex_count > 0
    assert he.edge_count > 0
    assert not he.repair_report.is_clean()
    assert he.repair_report.faces_rejected_by_ingest == 0


def test_strict_constructor_refuses_a_dirty_mesh(bunny):
    with pytest.raises(ValueError):
        HalfEdgeMesh3(bunny)


def test_strict_constructor_accepts_a_clean_mesh():
    mesh = Mesh3.create_box(2.0, 2.0, 2.0)
    he = HalfEdgeMesh3(mesh)

    assert he.face_count == len(mesh.faces)
    assert he.repair_report.is_clean()


def test_repair_opts_none_disables_every_pass():
    opts = RepairOpts.none()

    assert not opts.drop_degenerate
    assert not opts.drop_duplicate_faces
    assert not opts.resolve_nonmanifold_edges
    assert not opts.split_bowtie_vertices
    assert not opts.orient_consistently
    assert not opts.drop_isolated_vertices


def test_guaranteed_decimation_reduces_and_reports(bunny):
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    before = he.face_count

    report = he.decimate_guaranteed(DecimateOpts(TOL))

    assert report.collapses > 0
    assert report.faces_before == before
    assert report.faces_after < before
    assert report.ratio() == pytest.approx(report.faces_after / report.faces_before)
    assert he.face_count == report.faces_after

    stats = report.stats
    assert stats.evaluations > 0
    assert stats.vetoes() > 0
    assert repr(report)


def test_best_effort_is_the_more_aggressive_path(bunny):
    guaranteed = HalfEdgeMesh3(bunny, repair=RepairOpts())
    guaranteed.decimate_guaranteed(DecimateOpts(TOL))

    best = HalfEdgeMesh3(bunny, repair=RepairOpts())
    best.decimate_best_effort(BestEffortOpts(TOL))

    assert best.face_count < guaranteed.face_count


def test_a_guarantee_cannot_follow_a_best_effort_pass(bunny):
    """
    The refusal that protects the tolerance claim has to survive the binding boundary, since this is
    the whole reason `HalfEdgeMesh3` is the Python surface rather than a thin wrapper over collapses.
    """
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    assert he.error_volume_is_current

    he.decimate_best_effort(BestEffortOpts(TOL))
    assert not he.error_volume_is_current

    with pytest.raises(ValueError, match="best-effort"):
        he.decimate_guaranteed(DecimateOpts(TOL))

    # The other order is fine, as is repeating either path on its own.
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    he.decimate_guaranteed(DecimateOpts(TOL))
    he.decimate_guaranteed(DecimateOpts(TOL))
    he.decimate_best_effort(BestEffortOpts(TOL))


def test_face_count_target(bunny):
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    target = he.face_count // 2

    report = he.decimate_guaranteed(DecimateOpts(1.0, face_count=target))

    assert report.faces_after <= target


def test_ratio_target(bunny):
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    before = he.face_count

    report = he.decimate_guaranteed(DecimateOpts(1.0, ratio=0.5))

    assert report.faces_after <= round(before * 0.5)


def test_the_two_targets_are_mutually_exclusive():
    with pytest.raises(ValueError, match="at most one"):
        DecimateOpts(0.01, face_count=100, ratio=0.5)
    with pytest.raises(ValueError, match="at most one"):
        BestEffortOpts(0.01, face_count=100, ratio=0.5)


def test_options_round_trip_through_the_boundary():
    opts = DecimateOpts(
        0.02,
        ratio=0.25,
        max_normal_deviation=0.5,
        min_aspect=0.05,
        lock_boundary=False,
        boundary_tol=0.001,
        placement="endpoint",
        quadric="probabilistic",
        error_method="affine_face_map",
    )

    assert opts.tol == 0.02
    assert opts.ratio == 0.25
    assert opts.face_count is None
    assert opts.max_normal_deviation == 0.5
    assert opts.min_aspect == 0.05
    assert opts.lock_boundary is False
    assert opts.boundary_tol == 0.001
    assert opts.placement == "endpoint"
    assert opts.quadric == "probabilistic"
    assert opts.error_method == "affine_face_map"


def test_defaults_are_the_core_defaults():
    opts = DecimateOpts(0.02)

    assert opts.face_count is None
    assert opts.ratio is None
    assert opts.lock_boundary is True
    assert opts.boundary_tol is None
    assert opts.placement == "optimal"
    assert opts.quadric == "triangle"
    assert opts.error_method == "projected_overlay"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"placement": "nope"},
        {"quadric": "nope"},
        {"error_method": "nope"},
    ],
)
def test_bad_tokens_are_refused_by_name(kwargs):
    with pytest.raises(ValueError, match="Invalid"):
        DecimateOpts(0.01, **kwargs)


def test_invalid_options_are_refused(bunny):
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())

    with pytest.raises(ValueError):
        he.decimate_guaranteed(DecimateOpts(0.0))
    with pytest.raises(ValueError):
        he.decimate_guaranteed(DecimateOpts(0.01, ratio=1.5))

    # A refused call must not mark the mesh, since nothing was done to it.
    with pytest.raises(ValueError):
        he.decimate_best_effort(BestEffortOpts(-1.0))
    assert he.error_volume_is_current


def test_every_placement_and_quadric_is_accepted(bunny):
    for placement in ("optimal", "midpoint", "endpoint"):
        for quadric in ("triangle", "probabilistic"):
            he = HalfEdgeMesh3(bunny, repair=RepairOpts())
            report = he.decimate_guaranteed(
                DecimateOpts(TOL, placement=placement, quadric=quadric)
            )
            assert report.faces_after <= report.faces_before


def test_every_error_method_is_accepted(bunny):
    for method in ("contraction", "affine_face_map", "projected_overlay"):
        he = HalfEdgeMesh3(bunny, repair=RepairOpts())
        report = he.decimate_guaranteed(DecimateOpts(TOL, error_method=method))
        assert report.faces_after <= report.faces_before


def test_endpoint_placement_keeps_the_points_a_subset(bunny):
    """
    A binding-level check that `placement` is actually reaching the core: with endpoint placement
    every surviving point is one that was in the input.
    """
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    he.decimate_guaranteed(DecimateOpts(0.05, placement="endpoint"))
    out = he.to_mesh()

    original = {tuple(p) for p in bunny.points}
    for point in out.points:
        assert tuple(point) in original


def test_smoothing_moves_points_and_keeps_the_topology(bunny):
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    before = he.to_mesh()

    he.smooth(2)
    after = he.to_mesh()

    assert after.points.shape == before.points.shape
    assert len(after.faces) == len(before.faces)
    assert not (after.points == before.points).all()


def test_zero_smoothing_passes_change_nothing(bunny):
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    before = he.to_mesh()

    he.smooth(0)
    after = he.to_mesh()

    assert (after.points == before.points).all()


def test_round_trip_and_is_solid(bunny):
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())

    surface = he.to_mesh()
    assert len(surface.faces) == he.face_count

    solid = he.to_mesh(is_solid=True)
    assert len(solid.faces) == len(surface.faces)


def test_the_mesh_is_reusable_after_a_conversion(bunny):
    """`to_mesh` builds a copy, so the session stays usable and can be decimated further."""
    he = HalfEdgeMesh3(bunny, repair=RepairOpts())
    first = he.to_mesh()

    he.decimate_guaranteed(DecimateOpts(TOL))
    second = he.to_mesh()

    assert len(second.faces) < len(first.faces)
