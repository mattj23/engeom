from __future__ import annotations
from enum import Enum


class AngleDir(Enum):
    """Enumerates the two possible directions of rotation, clockwise and counter-clockwise."""

    Cw = 0
    """Clockwise rotation (negative angular direction)."""

    Ccw = 1
    """Counter-clockwise rotation (positive angular direction)."""


def angle_in_direction(radians0: float, radians1: float, angle_dir: AngleDir) -> float:
    """
    Returns the positive angle (in radians) needed to rotate `radians0` to `radians1` in the
    given direction. The result is always in the range [0, 2π].

    :param radians0: The starting angle, in radians.
    :param radians1: The ending angle, in radians.
    :param angle_dir: The rotational direction to consider.
    :return: The positive arc length in the given direction, in radians.
    """
    ...


def shortest_angle_between(radians0: float, radians1: float) -> float:
    """
    Returns the signed shortest angular distance from `radians0` to `radians1`.

    A positive result means the shortest path is counter-clockwise; a negative result means
    the shortest path is clockwise. The magnitude is always in the range [0, π].

    :param radians0: The starting angle, in radians.
    :param radians1: The ending angle, in radians.
    :return: The signed shortest angular distance, in radians.
    """
    ...


def angle_signed_pi(radians: float) -> float:
    """
    Re-expresses an angle in the range (-π, π]. Equivalent angles are preserved; -π maps to π.

    :param radians: The angle to re-express, in radians.
    :return: The equivalent angle in (-π, π].
    """
    ...


def angle_to_2pi(radians: float) -> float:
    """
    Re-expresses an angle in the range [0, 2π].

    :param radians: The angle to re-express, in radians.
    :return: The equivalent angle in [0, 2π].
    """
    ...


def signed_compliment_2pi(radians: float) -> float:
    """
    Returns the signed complement of an angle with respect to a full rotation (2π).

    A positive input returns a negative complement, and vice versa — together they sum to ±2π.
    The result is in (-2π, 2π].

    :param radians: The angle to complement, in radians.
    :return: The signed complement angle, in radians.
    """
    ...
