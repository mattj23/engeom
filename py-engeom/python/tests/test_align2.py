"""Tests for the 2-D alignment bindings."""

from __future__ import annotations

import math

import numpy
import pytest

from engeom.align2 import (
    AlignParams2,
    Dof3,
    multi_curve_adjustment,
    points_to_curve,
    points_to_group,
    points_to_cloud,
)
from engeom.geom2 import Curve2, CurveGroup2, Iso2


def rect_points(x0=-5.0, y0=-3.0, x1=5.0, y1=3.0) -> numpy.ndarray:
    """The four corners of a counter-clockwise rectangle, closed by repeating the first."""
    return numpy.array(
        [[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]], dtype=numpy.float64
    )


def rect_curve() -> Curve2:
    return Curve2(rect_points(), tol=1e-8)


def sampled(curve: Curve2, spacing: float = 0.25) -> numpy.ndarray:
    """Points walked around a curve at a uniform arc length spacing."""
    n = int(curve.length / spacing)
    return numpy.array(
        [[s.point.x, s.point.y] for s in (curve.at_length((k + 0.5) * spacing) for k in range(n))],
        dtype=numpy.float64,
    )


def moved(points: numpy.ndarray, t: Iso2) -> numpy.ndarray:
    return t.transform_points(points)


class TestDof3:
    def test_defaults_are_all_free(self):
        d = Dof3()
        assert (d.tx, d.ty, d.rz) == (True, True, True)
        assert (Dof3.all().tx, Dof3.all().ty, Dof3.all().rz) == (True, True, True)

    def test_locks_are_reported(self):
        d = Dof3(tx=False, ty=True, rz=False)
        assert not d.tx
        assert d.ty
        assert not d.rz


class TestAlignParams2:
    def test_center_and_local_are_mutually_exclusive(self):
        from engeom.geom2 import Point2

        with pytest.raises(ValueError):
            AlignParams2(center=Point2(1.0, 2.0), local=Iso2.identity())

    def test_dof_round_trips(self):
        params = AlignParams2(dof=Dof3(tx=False))
        assert not params.dof.tx
        assert params.dof.ty


class TestPointsToCurve:
    def test_a_disturbed_point_set_is_recovered(self):
        curve = rect_curve()
        pts = sampled(curve)
        disturb = Iso2(0.15, -0.1, 0.02)

        outcome = points_to_curve(moved(pts, disturb), curve, AlignParams2())

        # The recovered transform should undo the disturbance.
        composed = outcome.alignment.full_transform @ disturb
        numpy.testing.assert_allclose(
            composed.as_numpy(), numpy.eye(3), atol=1e-4
        )

    def test_the_outcome_reports_its_solves(self):
        curve = rect_curve()
        pts = sampled(curve)
        outcome = points_to_curve(
            moved(pts, Iso2(0.05, 0.05, 0.0)), curve, AlignParams2()
        )

        assert outcome.quality in {"converged", "unconverged"}
        assert isinstance(outcome.converged, bool)
        assert len(outcome.solves) >= 1
        assert outcome.refinement_rounds == len(outcome.solves) - 1

    def test_residuals_come_back_as_a_float_array(self):
        curve = rect_curve()
        pts = sampled(curve)
        outcome = points_to_curve(
            moved(pts, Iso2(0.05, 0.0, 0.0)), curve, AlignParams2()
        )

        residuals = outcome.alignment.residuals
        assert residuals.dtype == numpy.float64
        assert residuals.shape == (len(pts),)

        mean, std = outcome.alignment.residual_mean_std_dev
        assert math.isfinite(mean)
        assert std >= 0.0

    def test_a_locked_dof_is_not_recovered(self):
        curve = rect_curve()
        pts = sampled(curve)
        disturb = Iso2(0.3, 0.0, 0.0)

        outcome = points_to_curve(
            moved(pts, disturb), curve, AlignParams2(dof=Dof3(tx=False))
        )

        # tx is locked, so the x displacement cannot be undone.
        composed = outcome.alignment.full_transform @ disturb
        assert abs(composed.as_numpy()[0, 2]) > 0.2

    def test_refinement_can_be_switched_off(self):
        curve = rect_curve()
        pts = sampled(curve)
        outcome = points_to_curve(
            moved(pts, Iso2(0.05, 0.0, 0.0)),
            curve,
            AlignParams2(),
            refinement_steps=0,
        )
        assert outcome.refinement_rounds == 0


class TestPointsToGroup:
    def test_a_two_member_body_is_recovered(self):
        left = Curve2(rect_points(-8.0, -2.0, -3.0, 2.0), tol=1e-8)
        right = Curve2(rect_points(3.0, -2.0, 8.0, 2.0), tol=1e-8)
        group = CurveGroup2([left, right])

        pts = numpy.vstack([sampled(left), sampled(right)])
        disturb = Iso2(0.12, -0.08, 0.015)

        outcome = points_to_group(moved(pts, disturb), group, AlignParams2())

        composed = outcome.alignment.full_transform @ disturb
        numpy.testing.assert_allclose(
            composed.as_numpy(), numpy.eye(3), atol=1e-4
        )


class TestPointsToPointSet:
    def test_a_disturbed_point_set_is_recovered(self):
        curve = rect_curve()
        spacing = 0.25
        n = int(curve.length / spacing)
        stations = [curve.at_length((k + 0.5) * spacing) for k in range(n)]

        target_points = numpy.array(
            [[s.point.x, s.point.y] for s in stations], dtype=numpy.float64
        )
        target_normals = numpy.array(
            [[s.normal.x, s.normal.y] for s in stations], dtype=numpy.float64
        )

        disturb = Iso2(0.1, -0.08, 0.015)
        outcome = points_to_cloud(
            moved(target_points, disturb),
            target_points,
            target_normals,
            AlignParams2(),
            1.0,
        )

        composed = outcome.alignment.full_transform @ disturb
        numpy.testing.assert_allclose(
            composed.as_numpy(), numpy.eye(3), atol=1e-3
        )

    def test_a_nonsense_extrapolation_bound_is_rejected(self):
        pts = numpy.array([[0.0, 0.0], [1.0, 0.0]], dtype=numpy.float64)
        normals = numpy.array([[0.0, 1.0], [0.0, 1.0]], dtype=numpy.float64)

        with pytest.raises(ValueError):
            points_to_cloud(pts, pts, normals, AlignParams2(), 0.0)

    def test_mismatched_normals_are_rejected(self):
        pts = numpy.array([[0.0, 0.0], [1.0, 0.0]], dtype=numpy.float64)
        normals = numpy.array([[0.0, 1.0]], dtype=numpy.float64)

        with pytest.raises(ValueError):
            points_to_cloud(pts, pts, normals, AlignParams2(), 1.0)

    def test_a_zero_length_normal_is_rejected(self):
        # A zero row has no direction, and silently normalizing it would feed NaN into the solve.
        pts = numpy.array([[0.0, 0.0], [1.0, 0.0]], dtype=numpy.float64)
        normals = numpy.array([[0.0, 1.0], [0.0, 0.0]], dtype=numpy.float64)

        with pytest.raises(ValueError):
            points_to_cloud(pts, pts, normals, AlignParams2(), 1.0)


class TestMultiCurveAdjustment:
    def test_three_bodies_settle_together(self):
        groups = [CurveGroup2([rect_curve()]) for _ in range(3)]
        initial = [Iso2.identity(), Iso2(0.15, -0.1, 0.01), Iso2(-0.1, 0.12, -0.012)]

        outcome = multi_curve_adjustment(
            groups, 20.0, initial=initial, sample_spacing=0.5
        )

        assert len(outcome) == 3
        assert outcome.quality in {"converged", "unconverged"}

        # Every pair of bodies should end up coincident, which is the same thing as every
        # relative transform being the identity.
        transforms = [a.full_transform for a in outcome.alignments]
        for i in range(3):
            for j in range(3):
                relative = transforms[i].inverse() @ transforms[j]
                numpy.testing.assert_allclose(
                    relative.as_numpy(), numpy.eye(3), atol=1e-4
                )

    def test_fewer_than_two_bodies_is_rejected(self):
        with pytest.raises(ValueError):
            multi_curve_adjustment([CurveGroup2([rect_curve()])], 20.0)

    def test_a_mismatched_initial_list_is_rejected(self):
        groups = [CurveGroup2([rect_curve()]) for _ in range(2)]
        with pytest.raises(ValueError):
            multi_curve_adjustment(groups, 20.0, initial=[Iso2.identity()])

    def test_the_outcome_indexes_its_bodies(self):
        groups = [CurveGroup2([rect_curve()]) for _ in range(2)]
        outcome = multi_curve_adjustment(groups, 20.0, sample_spacing=0.5)

        assert outcome.alignment(0) is not None
        with pytest.raises(IndexError):
            outcome.alignment(5)
