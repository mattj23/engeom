from __future__ import annotations
from typing import Iterable, Literal

import numpy
from numpy.typing import NDArray

AngleDir = Literal["cw", "ccw"]
"""A direction of rotation: ``"cw"`` (clockwise, negative) or ``"ccw"`` (counter-clockwise, positive)."""


class AngleInterval:
    """
    A continuous range of angles, specified by a starting angle and a positive (counter-clockwise)
    included length.

    This is a read-only view returned by other geometric queries (such as
    ``Circle2.intersection_interval`` and ``Arc2.angle_interval``); it is not directly
    constructible from Python. Only the core query surface is exposed here; the full interval
    algebra (overlap, intersection, union, etc.) is not currently bound.
    """

    @property
    def min(self) -> float:
        """The starting angle of the interval, in radians, in the range [0, 2π]."""
        ...

    @property
    def max(self) -> float:
        """
        The ending angle of the interval, in radians. If less than ``min``, the interval wraps
        beyond 2π.
        """
        ...

    @property
    def extent(self) -> float:
        """The angular extent (included angle) of the interval, in radians."""
        ...

    @property
    def center(self) -> float:
        """The angle at the midpoint of the interval, in radians."""
        ...

    def is_empty(self) -> bool:
        """Return whether the interval contains zero angular extent."""
        ...

    def contains_value(self, x: float) -> bool:
        """
        Return whether the given angle (in radians) falls within the interval.

        :param x: the angle to test, in radians.
        :return: True if the angle is within the interval, False otherwise.
        """
        ...


class IndexMask:
    """
    A fixed-length mask of boolean values over a contiguous range of indices, used throughout the
    library to select subsets of an entity's elements (a mesh's faces, a point cloud's points, and
    so on).

    A mask has a length fixed at construction, and every index from ``0`` to ``len - 1`` is either
    ``True`` (selected) or ``False`` (not selected). It is the counterpart to a plain list of
    indices: the mask form is what makes set algebra between two selections cheap, which is why the
    library prefers it.

    The set operations are exposed as Python operators rather than named methods: ``~`` inverts,
    ``|`` unions, ``&`` intersects, and ``-`` takes the difference. The augmented forms (``|=``,
    ``&=``, ``-=``) modify the mask in place instead of allocating a new one. All of the binary
    operations require both masks to be the same length.
    """

    def __init__(self, length: int, value: bool = False):
        """
        Create a mask of the given length, with every index set to the same initial value.

        :param length: the number of indices covered by the mask.
        :param value: the initial value for every index.
        """
        ...

    @staticmethod
    def from_indices(indices: Iterable[int], length: int) -> IndexMask:
        """
        Create a mask of the given length in which only the listed indices are ``True``.

        :param indices: the indices to set to ``True``.
        :param length: the total number of indices covered by the mask.
        :return: the new mask.
        :raises ValueError: if any index is out of bounds for the given length.
        """
        ...

    @staticmethod
    def from_bool_array(values: NDArray[numpy.bool_]) -> IndexMask:
        """
        Create a mask from a 1D numpy array of booleans, one per index.

        :param values: a 1D array with dtype ``bool``. The mask takes its length from the array.
        :return: the new mask.
        """
        ...

    def to_indices(self) -> NDArray[numpy.uint64]:
        """
        Return the indices which are ``True``, in ascending order.

        :return: a 1D array of ``uint64`` indices.
        """
        ...

    def to_bool_array(self) -> NDArray[numpy.bool_]:
        """
        Return the mask as a 1D numpy array of booleans with one entry per index.

        :return: a 1D array of ``bool`` with the same length as the mask.
        """
        ...

    def count_true(self) -> int:
        """
        Return the number of indices which are ``True``.
        """
        ...

    def any(self) -> bool:
        """
        Return whether at least one index is ``True``.

        Note that this is about the mask's contents, not its length: ``len()`` is the size of the
        range the mask covers, and a mask of a thousand indices which are all ``False`` still has a
        length of a thousand.
        """
        ...

    def all(self) -> bool:
        """
        Return whether every index is ``True``. A zero length mask returns ``True``.
        """
        ...

    def copy(self) -> IndexMask:
        """
        Return an independent copy of this mask. Because a mask is mutable, assigning one to
        another name aliases it rather than duplicating it, and this is the way to get a version
        which can be modified without disturbing the original.

        :return: a new mask with the same length and contents.
        """
        ...

    def fill(self, value: bool):
        """
        Set every index in the mask to the same value, in place.

        :param value: the value to set every index to.
        """
        ...

    def pop_true(self) -> int | None:
        """
        Set the lowest ``True`` index to ``False`` and return it, or return ``None`` if no index is
        ``True``. This lets a mask be consumed as a work queue.

        :return: the index which was cleared, or None if the mask had no ``True`` index.
        """
        ...

    def __len__(self) -> int: ...

    def __getitem__(self, index: int) -> bool: ...

    def __setitem__(self, index: int, value: bool): ...

    def __invert__(self) -> IndexMask: ...

    def __or__(self, other: IndexMask) -> IndexMask: ...

    def __and__(self, other: IndexMask) -> IndexMask: ...

    def __sub__(self, other: IndexMask) -> IndexMask: ...

    def __ior__(self, other: IndexMask) -> IndexMask: ...

    def __iand__(self, other: IndexMask) -> IndexMask: ...

    def __isub__(self, other: IndexMask) -> IndexMask: ...


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

    A positive input returns a negative complement, and vice versa - together they sum to ±2π.
    The result is in (-2π, 2π].

    :param radians: The angle to complement, in radians.
    :return: The signed complement angle, in radians.
    """
    ...
