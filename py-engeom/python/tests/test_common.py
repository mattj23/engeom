"""
Tests for the engeom.common module: angle direction tokens, angle utility functions, and the
index mask type.
"""
import math

import numpy
import pytest
from engeom.common import (IndexMask, angle_in_direction, shortest_angle_between, angle_signed_pi,
                           angle_to_2pi, signed_compliment_2pi)

TOL = 1e-10
PI = math.pi


# ---------------------------------------------------------------------------
# angle direction tokens
# ---------------------------------------------------------------------------

def test_angle_in_direction_rejects_unknown_token():
    with pytest.raises(ValueError):
        angle_in_direction(0.0, PI / 2, "clockwise")


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
    (0.0, PI / 2, "ccw", PI / 2),
    (0.0, PI / 2, "cw", 3 * PI / 2),
    # Half-turn - same in both directions
    (0.0, PI, "ccw", PI),
    (0.0, PI, "cw", PI),
    # Full circle (same start and end)
    (0.0, 0.0, "ccw", 0.0),
    (0.0, 0.0, "cw", 0.0),
])
def test_angle_in_direction_values(start, end, direction, expected):
    assert angle_in_direction(start, end, direction) == pytest.approx(expected, abs=TOL)


def test_angle_in_direction_always_positive():
    for start in [-PI, -PI / 2, 0.0, PI / 4, PI]:
        for end in [-PI, -PI / 2, 0.0, PI / 4, PI]:
            for d in ("cw", "ccw"):
                assert angle_in_direction(start, end, d) >= -TOL


def test_angle_in_direction_rotation_check():
    """Rotating start by the returned angle in the stated direction must reach end."""
    cases = [
        (0.3, 1.8, "ccw"),
        (1.8, 0.3, "ccw"),
        (0.3, 1.8, "cw"),
        (-2.5, 2.5, "ccw"),
    ]
    for start, end, direction in cases:
        arc = angle_in_direction(start, end, direction)
        sign = 1.0 if direction == "ccw" else -1.0
        rotated = start + sign * arc
        assert math.sin(rotated) == pytest.approx(math.sin(end), abs=TOL)
        assert math.cos(rotated) == pytest.approx(math.cos(end), abs=TOL)


# ---------------------------------------------------------------------------
# shortest_angle_between
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("start,end,expected", [
    (0.0, PI / 2, PI / 2),  # CCW is shorter → positive
    (0.0, -PI / 4, -PI / 4),  # CW is shorter → negative
    (0.0, PI, PI),  # exactly half - positive by convention (ccw ≤ cw)
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


# ---------------------------------------------------------------------------
# IndexMask
# ---------------------------------------------------------------------------

def test_index_mask_starts_uniform_and_reports_length():
    mask = IndexMask(10)
    assert len(mask) == 10
    assert mask.count_true() == 0
    assert not mask.any()
    assert not mask.all()

    filled = IndexMask(10, True)
    assert filled.count_true() == 10
    assert filled.any()
    assert filled.all()


def test_index_mask_zero_length_is_neither_any_nor_nonempty():
    """A zero length mask is a legal degenerate case, not a panic."""
    mask = IndexMask(0)
    assert len(mask) == 0
    assert not mask.any()
    assert mask.all()
    assert mask.to_indices().size == 0


def test_index_mask_get_and_set_by_index():
    mask = IndexMask(5)
    mask[2] = True
    assert mask[2]
    assert not mask[0]
    mask[2] = False
    assert not mask[2]


def test_index_mask_supports_negative_indices():
    mask = IndexMask(5)
    mask[-1] = True
    assert mask[4]
    assert mask[-1]
    assert list(mask.to_indices()) == [4]


@pytest.mark.parametrize("index", [5, 6, -6, -100])
def test_index_mask_out_of_range_index_raises(index):
    mask = IndexMask(5)
    with pytest.raises(IndexError):
        _ = mask[index]
    with pytest.raises(IndexError):
        mask[index] = True


def test_index_mask_from_indices_round_trips():
    mask = IndexMask.from_indices([0, 2, 5], 8)
    assert len(mask) == 8
    assert list(mask.to_indices()) == [0, 2, 5]


def test_index_mask_from_indices_rejects_out_of_bounds():
    with pytest.raises(ValueError):
        IndexMask.from_indices([1, 12], 4)


def test_index_mask_to_indices_is_unsigned():
    """Index arrays cross the boundary as an unsigned dtype, not the numpy default int64."""
    indices = IndexMask.from_indices([1, 3], 8).to_indices()
    assert indices.dtype == numpy.uint64


def test_index_mask_bool_array_round_trips():
    values = numpy.array([True, False, True, True, False])
    mask = IndexMask.from_bool_array(values)
    assert len(mask) == 5
    result = mask.to_bool_array()
    assert result.dtype == numpy.bool_
    assert numpy.array_equal(result, values)


def test_index_mask_union_intersection_and_difference():
    left = IndexMask.from_indices([0, 2, 5], 8)
    right = IndexMask.from_indices([1, 2, 6], 8)

    assert list((left | right).to_indices()) == [0, 1, 2, 5, 6]
    assert list((left & right).to_indices()) == [2]
    assert list((left - right).to_indices()) == [0, 5]
    assert list((~left).to_indices()) == [1, 3, 4, 6, 7]


def test_index_mask_binary_operators_leave_operands_alone():
    left = IndexMask.from_indices([0, 2, 5], 8)
    right = IndexMask.from_indices([1, 2, 6], 8)

    _ = left | right
    _ = left & right
    _ = left - right
    _ = ~left

    assert list(left.to_indices()) == [0, 2, 5]
    assert list(right.to_indices()) == [1, 2, 6]


@pytest.mark.parametrize("op", ["or", "and", "sub"])
def test_index_mask_in_place_operators_modify_the_mask(op):
    mask = IndexMask.from_indices([0, 2, 5], 8)
    other = IndexMask.from_indices([1, 2, 6], 8)

    if op == "or":
        mask |= other
        expected = [0, 1, 2, 5, 6]
    elif op == "and":
        mask &= other
        expected = [2]
    else:
        mask -= other
        expected = [0, 5]

    assert list(mask.to_indices()) == expected


@pytest.mark.parametrize("op", ["or", "and", "sub"])
def test_index_mask_operators_reject_mismatched_lengths(op):
    left = IndexMask(4, True)
    right = IndexMask(5)

    with pytest.raises(ValueError):
        if op == "or":
            _ = left | right
        elif op == "and":
            _ = left & right
        else:
            _ = left - right


def test_index_mask_copy_is_independent():
    original = IndexMask.from_indices([1, 3], 8)
    duplicate = original.copy()
    duplicate[5] = True

    assert list(original.to_indices()) == [1, 3]
    assert list(duplicate.to_indices()) == [1, 3, 5]


def test_index_mask_equality_compares_contents_and_length():
    mask = IndexMask.from_indices([1, 3], 8)
    assert mask == IndexMask.from_indices([1, 3], 8)
    assert mask != IndexMask.from_indices([1, 4], 8)
    assert mask != IndexMask.from_indices([1, 3], 9)
    assert mask != "not a mask"


def test_index_mask_fill():
    mask = IndexMask.from_indices([1, 3], 8)
    mask.fill(True)
    assert mask.all()
    mask.fill(False)
    assert not mask.any()


def test_index_mask_pop_true_drains_in_ascending_order():
    mask = IndexMask.from_indices([1, 4, 6], 8)
    assert [mask.pop_true() for _ in range(3)] == [1, 4, 6]
    assert mask.pop_true() is None
    assert not mask.any()
