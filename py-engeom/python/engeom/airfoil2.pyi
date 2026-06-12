from __future__ import annotations

from typing import List, Literal, Optional

import numpy as np
from numpy.typing import NDArray

from engeom.geom2 import Circle2, Curve2, Point2, SurfacePoint2, Vector2


class Inscribed:
    """
    An inscribed circle inside an airfoil section.

    The circle center lies on the mean camber line. The two contact points ``p0`` and ``p1``
    are on the section perimeter on opposite sides of the camber line. Their upper/lower
    assignment is ambiguous until the section has been oriented.
    """

    @property
    def c(self) -> Circle2:
        """The inscribed circle (center and radius)."""
        ...

    @property
    def p0(self) -> Point2:
        """First contact point on the section perimeter."""
        ...

    @property
    def p1(self) -> Point2:
        """Second contact point on the section perimeter, opposite side from ``p0``."""
        ...

    @property
    def center(self) -> Point2:
        """The center of the inscribed circle, which lies on the camber line."""
        ...

    @property
    def radius(self) -> float:
        """The radius of the inscribed circle."""
        ...

    def camber_point(self) -> SurfacePoint2:
        """
        Return a ``SurfacePoint2`` at the circle center whose normal direction estimates the
        local camber line tangent. The direction is the unit vector from ``p0`` to ``p1``
        rotated 90 degrees clockwise.
        """
        ...

    def contact_dir(self) -> Vector2:
        """Return the vector from ``p0`` to ``p1``."""
        ...


class AfEdgeGeometry:
    """
    The geometric description of a fitted airfoil edge.

    The ``kind`` property identifies the variant. Depending on ``kind``, specific properties are
    populated and the rest are ``None``:

    - ``"open"``: no geometric data.
    - ``"sharp"``: ``point`` is the apex.
    - ``"square"``: ``corner0`` and ``corner1``.
    - ``"rounded_square"``: ``corner0``, ``corner1``, and ``radius``.
    - ``"full_round"``: ``point`` is the circle center, ``radius`` is the circle radius.
    - ``"blended_round"``: ``point`` is the circle center, ``radius`` is the circle radius.
    """

    @property
    def kind(self) -> Literal["open", "sharp", "square", "rounded_square", "full_round", "blended_round"]:
        """The variant name."""
        ...

    @property
    def point(self) -> Optional[Point2]:
        """
        The apex (``"sharp"``) or circle center (``"full_round"``, ``"blended_round"``).
        ``None`` for all other variants.
        """
        ...

    @property
    def corner0(self) -> Optional[Point2]:
        """First corner point. Populated for ``"square"`` and ``"rounded_square"`` only."""
        ...

    @property
    def corner1(self) -> Optional[Point2]:
        """Second corner point. Populated for ``"square"`` and ``"rounded_square"`` only."""
        ...

    @property
    def radius(self) -> Optional[float]:
        """
        Circle or fillet radius. Populated for ``"rounded_square"``, ``"full_round"``, and
        ``"blended_round"`` only.
        """
        ...


class AfEdge:
    """A fitted airfoil edge: the canonical edge location point and its geometric description."""

    @property
    def point(self) -> Point2:
        """
        The canonical edge location. For example: the apex for ``"sharp"``, the midpoint of the
        flat face for ``"square"`` and ``"rounded_square"``, or the outermost camber-axis
        point for ``"full_round"`` and ``"blended_round"``.
        """
        ...

    @property
    def geometry(self) -> AfEdgeGeometry:
        """The geometric description of the edge shape."""
        ...


class AfEdgeFit:
    """The result of an airfoil edge fitting operation."""

    @property
    def edge(self) -> AfEdge:
        """The fitted edge point and geometry."""
        ...

    @property
    def residuals(self) -> NDArray[np.float64]:
        """Point-to-boundary residuals from the fitting optimization, as a 1-D array."""
        ...

    @property
    def circles(self) -> List[Inscribed]:
        """
        The inscribed circle stack, potentially extended with one additional circle refined
        near the edge (for ``"blended_round"`` fits only).
        """
        ...


def fit_square_edge(section: Curve2, tol: float, circles: List[Inscribed], at_front: bool) -> AfEdgeFit:
    """
    Fit a square (flat) trailing or leading edge to airfoil section data.

    The edge geometry consists of two corner points connected by a straight flat face, with the
    edge point placed at the midpoint between them.

    :param section: the airfoil section curve.
    :param tol: general fitting tolerance.
    :param circles: the inscribed circle stack from the camber line fitting step.
    :param at_front: when ``True``, fit the leading edge; when ``False``, fit the trailing edge.
    :return: an ``AfEdgeFit`` whose edge geometry is of kind ``"square"``.
    :raises ValueError: if fitting fails.
    """
    ...


def fit_rounded_square_edge(section: Curve2, tol: float, circles: List[Inscribed], at_front: bool) -> AfEdgeFit:
    """
    Fit a rounded-square trailing or leading edge to airfoil section data.

    Like ``fit_square_edge``, but the two corners are replaced by circular arc fillets of a
    single optimised radius.

    :param section: the airfoil section curve.
    :param tol: general fitting tolerance.
    :param circles: the inscribed circle stack from the camber line fitting step.
    :param at_front: when ``True``, fit the leading edge; when ``False``, fit the trailing edge.
    :return: an ``AfEdgeFit`` whose edge geometry is of kind ``"rounded_square"``.
    :raises ValueError: if fitting fails.
    """
    ...


def fit_sharp_edge(section: Curve2, tol: float, circles: List[Inscribed], at_front: bool) -> AfEdgeFit:
    """
    Fit a sharp corner trailing or leading edge to airfoil section data.

    The edge geometry is a single apex point where the two airfoil surfaces meet.

    :param section: the airfoil section curve.
    :param tol: general fitting tolerance.
    :param circles: the inscribed circle stack from the camber line fitting step.
    :param at_front: when ``True``, fit the leading edge; when ``False``, fit the trailing edge.
    :return: an ``AfEdgeFit`` whose edge geometry is of kind ``"sharp"``.
    :raises ValueError: if fitting fails.
    """
    ...


def fit_full_round_edge(section: Curve2, tol: float, circles: List[Inscribed], at_front: bool) -> AfEdgeFit:
    """
    Fit a full-round (semicircular) trailing or leading edge to airfoil section data.

    The edge is a single circular arc spanning from the suction-side to the pressure-side
    contact point of the last inscribed circle.

    :param section: the airfoil section curve.
    :param tol: general fitting tolerance.
    :param circles: the inscribed circle stack from the camber line fitting step.
    :param at_front: when ``True``, fit the leading edge; when ``False``, fit the trailing edge.
    :return: an ``AfEdgeFit`` whose edge geometry is of kind ``"full_round"``.
    :raises ValueError: if fitting fails.
    """
    ...


def fit_blended_round_edge(section: Curve2, tol: float, circles: List[Inscribed], at_front: bool) -> AfEdgeFit:
    """
    Fit a blended-round trailing or leading edge to airfoil section data.

    The edge is a three-arc boundary: a blend arc from the suction surface, a central edge
    circle arc, and a blend arc from the pressure surface. Each blend arc is simultaneously
    tangent to its surface tangent line and internally tangent to the central circle.

    The returned inscribed circle stack includes one additional circle refined near the edge.

    :param section: the airfoil section curve.
    :param tol: general fitting tolerance.
    :param circles: the inscribed circle stack from the camber line fitting step.
    :param at_front: when ``True``, fit the leading edge; when ``False``, fit the trailing edge.
    :return: an ``AfEdgeFit`` whose edge geometry is of kind ``"blended_round"``.
    :raises ValueError: if fitting fails.
    """
    ...


def extract_inscribed_circles(section: Curve2, tol: float) -> List[Inscribed]:
    """
    Extract the unambiguous inscribed circles of an airfoil section.

    Circles are returned in a consistent order but may run from leading-to-trailing edge or
    the reverse. The ``p0``/``p1`` contact points are consistently oriented relative to each
    other but their upper/lower surface assignment is also ambiguous until the section is
    oriented.

    The algorithm terminates near the edges when the remaining distance to the farthest edge
    point exceeds 3/8 of the last circle's radius beyond its perimeter. Edge-aware algorithms
    should be used to extend the camber line from there.

    :param section: the airfoil section curve, may be open at one edge but not both.
    :param tol: tolerance controlling circle density; more circles are inserted until the
        interpolation error between neighbors falls below this value. Circle centers are
        resolved to ``tol / 10``.
    :return: a list of ``Inscribed`` circles ordered along the camber line.
    :raises ValueError: if the algorithm fails to find a valid starting crossing line or
        encounters an unrecoverable error during extraction.
    """
    ...
