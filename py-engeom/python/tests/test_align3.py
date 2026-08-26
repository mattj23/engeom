"""
Tests for the 3-D alignment bindings.

These focus on what the binding layer responsibilities: the numpy arrays crossing the boundary, the options reaching
the Rust side, the shape of the outcome object, and the error contract. The alignment implementation itself is tested
on the Rust side and is not re-checked here.

The error contract is the most complicated part: `points_to_mesh` raises `ValueError` only when there is no
answer at all, meaning the arguments were rejected or the initial solve broke down. A solve that merely exhausts its
evaluation budget returns normally, reporting ``quality == "unconverged"``, because the parameters it leaves behind
are the best it found and are technically valid geometry.
"""

from __future__ import annotations

import numpy
import pytest

from engeom.align3 import AlignParams3, Dof6, points_to_cloud, points_to_mesh
from engeom.geom3 import Iso3, Mesh3, MeshData3, Point3, PointCloud3


def box_mesh() -> Mesh3:
    return Mesh3.create_box(10.0, 5.0, 2.0)


def box_points(mesh: Mesh3, radius: float = 0.5) -> numpy.ndarray:
    """
    Poisson samples on the surface as an ``(n, 3)`` array.

    ``sample_poisson`` returns points and normals interleaved as ``(n, 6)``, so the positions have
    to be sliced out and made contiguous before they cross back over the boundary.
    """
    return numpy.ascontiguousarray(mesh.sample_poisson(radius)[:, :3])


def small_disturbance() -> Iso3:
    return Iso3.from_translation(0.3, -0.2, 0.1)


def large_disturbance() -> Iso3:
    """ Large enough that the solver needs more than a handful of evaluations to undo it. """
    return Iso3.from_translation(3.0, 2.0, 1.0)


def centered_params(points: numpy.ndarray) -> AlignParams3:
    return AlignParams3(center=Point3(*points.mean(axis=0)))


# ==================================================================================================
# The alignment round trip
# ==================================================================================================


def test_recovers_a_known_disturbance():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)

    outcome = points_to_mesh(moved, mesh, centered_params(moved))

    recovered = outcome.alignment.full_transform.transform_points(moved)
    assert numpy.abs(recovered - points).max() < 1e-6


def test_locked_dof_is_not_recovered():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = Iso3.from_translation(0.4, 0.0, 0.0).transform_points(points)

    dof = Dof6(tx=False)
    outcome = points_to_mesh(moved, mesh, AlignParams3(dof=dof))

    # With the local origin and offset both identity, the full transform is the alignment
    # transform, so a locked tx must leave the x translation at exactly zero.
    assert outcome.alignment.full_transform.transform_points(
        numpy.zeros((1, 3))
    )[0, 0] == pytest.approx(0.0, abs=0.0)


def test_residuals_are_one_dimensional_float64():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)

    residuals = points_to_mesh(moved, mesh, centered_params(moved)).alignment.residuals()

    assert residuals.dtype == numpy.float64
    assert residuals.shape == (len(points),)


def test_alignment_component_transforms_are_iso3():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)

    alignment = points_to_mesh(moved, mesh, centered_params(moved)).alignment

    for transform in (
        alignment.full_transform,
        alignment.local_transform,
        alignment.local_origin,
        alignment.offset,
    ):
        assert isinstance(transform, Iso3)


# ==================================================================================================
# The outcome report
# ==================================================================================================


def test_converged_outcome_reports_every_round():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)

    outcome = points_to_mesh(
        moved, mesh, centered_params(moved), refinement_steps=4, sigma_max=0.5
    )

    assert outcome.converged
    assert outcome.quality == "converged"
    assert outcome.refinement_rounds == 4
    # One entry for the initial solve plus one per refinement round.
    assert len(outcome.solves) == 5
    assert outcome.halt is None


def test_solve_terminations_are_snake_case_strings():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)

    solves = points_to_mesh(moved, mesh, centered_params(moved)).solves

    assert all(isinstance(s, str) for s in solves)
    # Not the underlying solver crate's Debug formatting, which would contain braces.
    assert all("{" not in s for s in solves)
    assert all(s.startswith(("converged", "residuals_zero", "orthogonal")) for s in solves)


def test_exhausted_budget_is_reported_not_raised():
    """The central contract: running out of evaluations is a quality report, not an error."""
    mesh = box_mesh()
    points = box_points(mesh)
    moved = large_disturbance().transform_points(points)

    outcome = points_to_mesh(moved, mesh, AlignParams3(), patience=1)

    assert outcome.quality == "unconverged"
    assert not outcome.converged
    assert "lost_patience" in outcome.solves
    # The alignment is still real geometry rather than something degenerate.
    assert numpy.isfinite(
        outcome.alignment.full_transform.transform_points(moved)
    ).all()


def test_refinement_can_be_disabled():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)

    outcome = points_to_mesh(moved, mesh, centered_params(moved), refinement_steps=0)

    assert outcome.refinement_rounds == 0
    assert len(outcome.solves) == 1
    assert outcome.halt is None


def test_repr_mentions_quality_and_rounds():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)

    text = repr(points_to_mesh(moved, mesh, centered_params(moved)))

    assert "AlignOutcome3" in text
    assert "quality" in text


# ==================================================================================================
# Uncertainty
# ==================================================================================================


def test_point_sigma_accepts_ndarray_and_list():
    mesh = box_mesh()
    points = box_points(mesh)
    moved = small_disturbance().transform_points(points)
    params = centered_params(moved)

    as_array = numpy.full(len(moved), 0.01)
    from_array = points_to_mesh(
        moved, mesh, params, point_sigma=as_array, sigma_max=3.0
    )
    from_list = points_to_mesh(
        moved, mesh, params, point_sigma=list(as_array), sigma_max=3.0
    )

    numpy.testing.assert_allclose(
        from_array.alignment.full_transform.transform_points(moved),
        from_list.alignment.full_transform.transform_points(moved),
    )


def test_point_sigma_changes_the_result():
    """A non-uniform sigma must actually reach the solver, not be silently dropped."""
    mesh = box_mesh()
    points = box_points(mesh)
    moved = points.copy()

    # Push one point off the surface so there is something for the weighting to act on.
    bad = len(moved) // 2
    moved[bad] += numpy.array([0.0, 0.0, 1.0])
    params = centered_params(moved)

    uniform = points_to_mesh(moved, mesh, params, refinement_steps=0)

    sigma = numpy.full(len(moved), 0.01)
    sigma[bad] = 100.0
    weighted = points_to_mesh(
        moved, mesh, params, refinement_steps=0, point_sigma=sigma
    )

    clean = [i for i in range(len(points)) if i != bad]
    uniform_dev = numpy.abs(
        uniform.alignment.full_transform.transform_points(moved)[clean] - points[clean]
    ).max()
    weighted_dev = numpy.abs(
        weighted.alignment.full_transform.transform_points(moved)[clean] - points[clean]
    ).max()

    assert weighted_dev < uniform_dev


def test_mesh_point_stdev_reaches_the_solver():
    """A mesh carrying per-vertex uncertainty must have it picked up as target-side sigma."""
    data = MeshData3.create_box(10.0, 5.0, 2.0)
    plain = data.to_mesh()

    # A non-uniform target sigma, since a uniform one is only a global scale on the objective and
    # would leave the minimizer unchanged.
    stdev = numpy.linspace(0.01, 1.0, len(data.points))
    data.set_point_stdev(stdev)
    measured = data.to_mesh()

    points = box_points(plain)
    moved = points.copy()
    moved[len(moved) // 2] += numpy.array([0.0, 0.0, 1.0])
    params = centered_params(moved)

    without = points_to_mesh(moved, plain, params, refinement_steps=0)
    with_stdev = points_to_mesh(moved, measured, params, refinement_steps=0)

    delta = numpy.abs(
        with_stdev.alignment.full_transform.transform_points(moved)
        - without.alignment.full_transform.transform_points(moved)
    ).max()
    assert delta > 1e-9, "the mesh's point_stdev had no effect on the alignment"


# ==================================================================================================
# Validation: the cases that genuinely have no answer
# ==================================================================================================


def test_zero_patience_is_rejected():
    mesh = box_mesh()
    points = box_points(mesh)

    with pytest.raises(ValueError, match="patience"):
        points_to_mesh(points, mesh, AlignParams3(), patience=0)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_sigma_max_is_rejected(bad):
    mesh = box_mesh()
    points = box_points(mesh)

    with pytest.raises(ValueError, match="sigma_max"):
        points_to_mesh(points, mesh, AlignParams3(), sigma_max=bad)


def test_point_sigma_length_is_validated():
    mesh = box_mesh()
    points = box_points(mesh)
    sigma = numpy.full(len(points) - 1, 0.1)

    with pytest.raises(ValueError, match="entries"):
        points_to_mesh(points, mesh, AlignParams3(), point_sigma=sigma)


@pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_point_sigma_values_are_rejected(bad):
    mesh = box_mesh()
    points = box_points(mesh)
    sigma = numpy.full(len(points), 0.1)
    sigma[3] = bad

    with pytest.raises(ValueError):
        points_to_mesh(points, mesh, AlignParams3(), point_sigma=sigma)


def test_wrong_point_shape_is_rejected():
    mesh = box_mesh()

    with pytest.raises(ValueError):
        points_to_mesh(numpy.zeros((10, 2)), mesh, AlignParams3())


# ================================================================================================
# points_to_cloud
# ================================================================================================


def _outward_normals(points: numpy.ndarray) -> numpy.ndarray:
    """
    Outward normals for points on the box surface.

    The box is axis aligned and centered, so a point's face is whichever axis it is closest to the
    extent of, and the normal is that axis signed by which side the point is on.
    """
    half = numpy.array([5.0, 2.5, 1.0])
    axes = numpy.argmax(numpy.abs(points) / half, axis=1)

    normals = numpy.zeros_like(points)
    rows = numpy.arange(len(points))
    normals[rows, axes] = numpy.sign(points[rows, axes])
    return normals


def box_cloud(spacing: float = 0.25) -> PointCloud3:
    """A sample of the box surface carrying the normals a cloud target requires."""
    points = box_points(box_mesh(), spacing)
    cloud = PointCloud3(points)
    cloud.set_point_normals(_outward_normals(points))
    return cloud


def test_points_to_cloud_recovers_a_known_offset():
    cloud = box_cloud()
    points = cloud.points.copy()

    shift = Iso3.from_translation(0.08, -0.05, 0.03)
    moved = numpy.array([list(shift @ Point3(*p)) for p in points])

    outcome = points_to_cloud(moved, cloud, AlignParams3(), max_extrapolation=1.0)

    recovered = outcome.alignment.full_transform
    residuals = outcome.alignment.residuals()
    assert residuals.shape == (len(points),)
    assert numpy.abs(residuals).max() < 1e-3

    # Applying the recovered transform to the moved points should put them back.
    back = numpy.array([list(recovered @ Point3(*p)) for p in moved])
    assert back == pytest.approx(points, abs=1e-3)


def test_points_to_cloud_requires_normals():
    """A cloud without normals cannot supply a tangent plane or a residual sign, so the binding
    must refuse rather than guess."""
    bare = PointCloud3(numpy.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]))
    with pytest.raises(ValueError):
        points_to_cloud(
            numpy.array([[0.0, 0.0, 1.0]]), bare, AlignParams3(), max_extrapolation=1.0
        )


def test_points_to_cloud_rejects_a_nonsense_extrapolation():
    cloud = box_cloud()
    with pytest.raises(ValueError):
        points_to_cloud(cloud.points, cloud, AlignParams3(), max_extrapolation=0.0)
