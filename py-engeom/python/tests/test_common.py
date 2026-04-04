"""
Tests for the engeom.common module: AngleDir enum and angle utility functions.
"""
import math

import numpy
import pytest
from engeom.common import (AngleDir, angle_in_direction, shortest_angle_between, angle_signed_pi, angle_to_2pi,
                           signed_compliment_2pi)

TOL = 1e-10
PI = math.pi


# ---------------------------------------------------------------------------
# AngleDir
# ---------------------------------------------------------------------------

def test_angle_dir_variants_exist():
    assert AngleDir.Cw is not None
    assert AngleDir.Ccw is not None


def test_angle_dir_equality():
    assert AngleDir.Cw == AngleDir.Cw
    assert AngleDir.Ccw == AngleDir.Ccw
    assert AngleDir.Cw != AngleDir.Ccw


# ---------------------------------------------------------------------------
# angle_to_2pi
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("radians,expected", [
    (0.0, 0.0),
    (PI, PI),
    (2 * PI, 0.0),  # exactly 2pi wraps to 0
    (-PI, PI),
    (-PI / 2, 3 * PI / 2),
    (3 * PI, PI),  # > 2pi wraps back
    (-2 * PI, 0.0),
])
def test_angle_to_2pi_range(radians, expected):
    assert angle_to_2pi(radians) == pytest.approx(expected, abs=TOL)


def test_angle_to_2pi_preserves_trig():
    for angle in [-7 * PI, -3.5, 0.0, 1.23, 4 * PI, 6.28]:
        r = angle_to_2pi(angle)
        assert 0.0 <= r < 2 * PI + TOL
        assert math.sin(r) == pytest.approx(math.sin(angle), abs=TOL)
        assert math.cos(r) == pytest.approx(math.cos(angle), abs=TOL)


# ---------------------------------------------------------------------------
# angle_signed_pi
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("radians,expected", [
    (0.0, 0.0),
    (PI / 2, PI / 2),
    (2.5 * PI, PI / 2),  # doctest from the Rust source
    (-PI / 2, -PI / 2),
    (3 * PI / 2, -PI / 2),  # > π wraps to negative
])
def test_angle_signed_pi_values(radians, expected):
    assert angle_signed_pi(radians) == pytest.approx(expected, abs=TOL)


def test_angle_signed_pi_range():
    for angle in [-7 * PI, -3.5, 0.0, 1.23, 4 * PI, 6.28]:
        r = angle_signed_pi(angle)
        assert -PI <= r <= PI + TOL
        assert math.sin(r) == pytest.approx(math.sin(angle), abs=TOL)
        assert math.cos(r) == pytest.approx(math.cos(angle), abs=TOL)


# ---------------------------------------------------------------------------
# signed_compliment_2pi
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("degrees,expected", [
    (90.0, -270.0),
    (180.0, -180.0),
    (270.0, -90.0),
    (-91.0, 269.0),
    (-181.0, 179.0),
    (-271.0, 89.0),
])
def test_signed_compliment_2pi_values(degrees, expected):
    assert signed_compliment_2pi(numpy.deg2rad(degrees)) == pytest.approx(numpy.deg2rad(expected), abs=TOL)


# ---------------------------------------------------------------------------
# angle_in_direction
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("start,end,direction,expected", [
    # Quarter turn each way from 0
    (0.0, PI / 2, AngleDir.Ccw, PI / 2),
    (0.0, PI / 2, AngleDir.Cw, 3 * PI / 2),
    # Half-turn — same in both directions
    (0.0, PI, AngleDir.Ccw, PI),
    (0.0, PI, AngleDir.Cw, PI),
    # Full circle (same start and end)
    (0.0, 0.0, AngleDir.Ccw, 0.0),
    (0.0, 0.0, AngleDir.Cw, 0.0),
])
def test_angle_in_direction_values(start, end, direction, expected):
    assert angle_in_direction(start, end, direction) == pytest.approx(expected, abs=TOL)


def test_angle_in_direction_always_positive():
    for start in [-PI, -PI / 2, 0.0, PI / 4, PI]:
        for end in [-PI, -PI / 2, 0.0, PI / 4, PI]:
            for d in (AngleDir.Cw, AngleDir.Ccw):
                assert angle_in_direction(start, end, d) >= -TOL


def test_angle_in_direction_rotation_check():
    """Rotating start by the returned angle in the stated direction must reach end."""
    cases = [
        (0.3, 1.8, AngleDir.Ccw),
        (1.8, 0.3, AngleDir.Ccw),
        (0.3, 1.8, AngleDir.Cw),
        (-2.5, 2.5, AngleDir.Ccw),
    ]
    for start, end, direction in cases:
        arc = angle_in_direction(start, end, direction)
        sign = 1.0 if direction == AngleDir.Ccw else -1.0
        rotated = start + sign * arc
        assert math.sin(rotated) == pytest.approx(math.sin(end), abs=TOL)
        assert math.cos(rotated) == pytest.approx(math.cos(end), abs=TOL)


# ---------------------------------------------------------------------------
# shortest_angle_between
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("start,end,expected", [
    (0.0, PI / 2, PI / 2),  # CCW is shorter → positive
    (0.0, -PI / 4, -PI / 4),  # CW is shorter → negative
    (0.0, PI, PI),  # exactly half — positive by convention (ccw ≤ cw)
    (0.0, 0.0, 0.0),  # same angle
    (PI / 2, -PI / 2, PI),  # half turn
])
def test_shortest_angle_between_values(start, end, expected):
    assert shortest_angle_between(start, end) == pytest.approx(expected, abs=TOL)


def test_shortest_angle_between_magnitude_le_pi():
    for start in [-PI, -PI / 2, 0.0, PI / 3, PI]:
        for end in [-PI, -PI / 2, 0.0, PI / 3, PI]:
            assert abs(shortest_angle_between(start, end)) <= PI


def test_shortest_angle_between_antisymmetric():
    """shortest_angle_between(a, b) == -shortest_angle_between(b, a) when |angle| < π."""
    pairs = [(0.0, 1.0), (0.5, -0.5), (-1.0, 1.0), (0.1, PI - 0.1)]
    for a, b in pairs:
        forward = shortest_angle_between(a, b)
        backward = shortest_angle_between(b, a)
        # Only strictly antisymmetric when not at the ±π boundary
        if abs(abs(forward) - PI) > TOL:
            assert forward + backward == pytest.approx(0.0, abs=TOL)
