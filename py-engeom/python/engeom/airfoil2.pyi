from __future__ import annotations

from typing import List, Literal, Optional

import numpy as np
from numpy.typing import NDArray

from engeom.geom2 import Circle2, Curve2, CurveStation2, Point2, SurfacePoint2, Vector2
from engeom.metrology import Distance2


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
    - ``"spline_max_k"``: ``point`` is the osculating circle center, ``radius`` its radius.
    """

    @property
    def kind(self) -> Literal[
        "open", "sharp", "square", "rounded_square", "full_round", "blended_round", "spline_max_k"
    ]:
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


def is_open_at_end(section: Curve2, tol: float, circles: List[Inscribed], at_front: bool) -> bool:
    """
    Determine whether the section is open at the end being processed.

    Whether an edge is open is a property of the section topology rather than something that can
    be measured against a candidate geometry, so this is the test used by ``"auto"`` to
    resolve the open case before any edge fitting is attempted.

    A section is open at this end when it is not a closed curve *and* both lips of the gap lie
    beyond the end of the inscribed circle stack. A section with a gap in the middle of a surface
    is open at neither end.

    :param section: the airfoil section curve.
    :param tol: general fitting tolerance.
    :param circles: the inscribed circle stack from the camber line fitting step.
    :param at_front: when ``True``, test the leading edge; when ``False``, test the trailing edge.
    :return: ``True`` if the section is open at the requested end.
    :raises ValueError: if there are fewer than two inscribed circles.
    """
    ...


def fit_open_edge(section: Curve2, tol: float, circles: List[Inscribed], at_front: bool) -> AfEdgeFit:
    """
    Fit an open trailing or leading edge on the airfoil section data.

    An open edge is one where the airfoil shape is incomplete, most commonly because the
    measurement system did not capture all the way around the airfoil, or because the nominal
    geometry is not a full airfoil (near the root or tip, where it blends into other geometry).

    Because there is no edge geometry to fit, the edge point is found by projecting the end of the
    camber line onto a chord spanning the gap, from the first point of the section curve to the
    last. The camber line is first advanced into the gap by repeatedly stepping halfway to the
    nearer lip and pushing a new refined inscribed circle, until the projected edge point stops
    moving by more than ``tol``.

    :param section: the airfoil section curve, which must be open.
    :param tol: general fitting tolerance.
    :param circles: the inscribed circle stack from the camber line fitting step.
    :param at_front: when ``True``, fit the leading edge; when ``False``, fit the trailing edge.
    :return: an ``AfEdgeFit`` whose edge geometry is of kind ``"open"``, with an empty residual
        array and an inscribed circle stack extended toward the gap.
    :raises ValueError: if the section is closed, or if the gap does not lie ahead of the end of
        the camber line.
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
    single optimized radius.

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


AfSide = Literal["upper", "lower"]
"""
Selects between the upper (suction, convex) and lower (pressure, concave) side of the airfoil.
"""


AfPos = Literal["on_camber", "radius", "edge_offset"]
"""
Selects how a gage point location on an airfoil surface is specified. For all three methods a
positive value is measured from the leading edge and a negative one from the trailing edge.

- ``"on_camber"``: an arc distance along the mean camber line. A line orthogonal to the camber is
  cast at that station and its intersection with the requested side is the point. This is the most
  commonly used method and performs well across the whole length of the airfoil.
- ``"radius"``: the intersection of the side with a circle of the given radius centered on the
  edge point. Best close to an edge, where the angle between the camber line and the surface
  tangent is large.
- ``"edge_offset"``: step from the edge point along the camber tangent *at that edge*, then cast
  orthogonally to that direction. Because the tangent diverges from a curved camber line, this
  only reaches the surface near the edge, which is what it is intended for.
"""


AfEdgeSearch = Literal[
    "auto",
    "open",
    "sharp",
    "square",
    "rounded_square",
    "full_round",
    "blended_round",
    "spline_max_k",
]
"""
Selects which edge geometry to fit at the leading or trailing edge during a geometric analysis.

- ``"auto"``: resolve the edge automatically. If the section is open at this end the edge is
  treated as ``"open"``; otherwise every applicable variant is fit and the lowest-residual result
  is selected.
- ``"open"``: treat the edge as open, meaning the section perimeter does not close at this end.
  No edge geometry is fit; instead the camber line is advanced into the gap and the edge point is
  projected onto a chord spanning the two open lips of the section.
- ``"sharp"``: fit a sharp apex edge.
- ``"square"``: fit a flat-faced edge with two sharp corners.
- ``"rounded_square"``: fit a flat-faced edge with two rounded corners of equal radius.
- ``"full_round"``: fit a full-round edge joined to the surfaces by short tangent segments.
- ``"blended_round"``: fit a full-round edge joined to the surfaces by tangent blending arcs.
- ``"spline_max_k"``: fit a cubic spline and take the point of maximum curvature as the edge.

These tokens match the ``kind`` reported by ``AfEdgeGeometry``, so an edge can be searched for by
the same name it reports back.
"""


class OrientFwdAft:
    """
    Selects the method used to identify which end of the camber line is the leading edge.

    Build one with a variant constructor rather than by calling ``OrientFwdAft()`` directly.
    """

    @staticmethod
    def tmax_fwd() -> OrientFwdAft:
        """
        Identify the leading edge as the end nearer the largest inscribed circle (the maximum
        thickness point). Suitable for typical subsonic airfoils.
        """
        ...

    @staticmethod
    def airflow(x: float, y: float) -> OrientFwdAft:
        """
        Identify the leading edge based on a supplied airflow direction. The end of the camber
        line further upstream of the airflow becomes the leading edge.

        :param x: x-component of the airflow direction.
        :param y: y-component of the airflow direction.
        """
        ...

    @staticmethod
    def fwd(x: float, y: float) -> OrientFwdAft:
        """
        Identify the leading edge using a vector that points from the trailing toward the leading
        edge. The end of the camber line that projects further along this vector becomes the
        leading edge.

        :param x: x-component of the forward direction vector.
        :param y: y-component of the forward direction vector.
        """
        ...


class OrientUpperLower:
    """
    Selects the method used to identify the upper (suction) and lower (pressure) surfaces of
    the airfoil after the forward/aft orientation has been resolved.

    Build one with a variant constructor rather than by calling ``OrientUpperLower()`` directly.
    """

    @staticmethod
    def curvature() -> OrientUpperLower:
        """
        Detect upper and lower from the curvature of the camber line. The more concave side
        becomes the lower (pressure) surface. Will fail or give poor results if the camber line
        is nearly straight.
        """
        ...

    @staticmethod
    def upper(x: float, y: float) -> OrientUpperLower:
        """
        Use a supplied direction vector to identify the upper surface. The side whose contact
        points are further along this vector becomes the upper surface.

        :param x: x-component of the upper-direction vector.
        :param y: y-component of the upper-direction vector.
        """
        ...

    @staticmethod
    def lower(x: float, y: float) -> OrientUpperLower:
        """
        Use a supplied direction vector to identify the lower surface. The side whose contact
        points are further along this vector becomes the lower surface.

        :param x: x-component of the lower-direction vector.
        :param y: y-component of the lower-direction vector.
        """
        ...


class AfGeometry:
    """
    The result of a geometric analysis of an airfoil section. Provides the leading and
    trailing edge fits, the mean camber line, the upper and lower surface curves, and the
    inscribed circle stack used during analysis.

    Construct with :meth:`from_geometric_analysis`.
    """

    @staticmethod
    def from_geometric_analysis(
        section: Curve2,
        general_tol: float,
        fwd_aft: OrientFwdAft,
        upper_lower: OrientUpperLower,
        le_search: AfEdgeSearch,
        te_search: AfEdgeSearch,
    ) -> AfGeometry:
        """
        Run a purely geometric analysis of an airfoil section.

        Internally this method:

        1. Extracts the unambiguous inscribed circles of the section.
        2. Orients the inscribed circles forward/aft using ``fwd_aft``.
        3. Orients them upper/lower using ``upper_lower``.
        4. Fits the leading and trailing edge geometries using ``le_search`` and ``te_search``.
        5. Builds the mean camber line from the edges and circle centers.
        6. Splits the section into upper and lower surface curves.

        :param section: the airfoil section curve.
        :param general_tol: general fitting tolerance used throughout the analysis.
        :param fwd_aft: forward/aft orientation method.
        :param upper_lower: upper/lower orientation method.
        :param le_search: edge-search method for the leading edge.
        :param te_search: edge-search method for the trailing edge.
        :return: a populated :class:`AfGeometry`.
        :raises ValueError: if any step of the analysis fails.
        """
        ...

    @property
    def leading(self) -> AfEdge:
        """The fitted leading edge."""
        ...

    @property
    def trailing(self) -> AfEdge:
        """The fitted trailing edge."""
        ...

    @property
    def camber(self) -> Curve2:
        """
        The mean camber line, oriented so the first point is the leading edge point and the
        last point is the trailing edge point.
        """
        ...

    @property
    def upper(self) -> Curve2:
        """The upper (suction) surface curve."""
        ...

    @property
    def lower(self) -> Curve2:
        """The lower (pressure) surface curve."""
        ...

    @property
    def circles(self) -> List[Inscribed]:
        """
        The inscribed circle stack used in the analysis, oriented leading-to-trailing with
        each circle's ``p0`` on the lower surface and ``p1`` on the upper surface.
        """
        ...

    def max_thickness_circle(self) -> Inscribed:
        """
        Return the inscribed circle with the largest radius. This corresponds to the maximum
        thickness location along the camber line.
        """
        ...

    def af_point(self, side: AfSide, method: AfPos, value: float) -> Optional[CurveStation2]:
        """
        Locate a gage point on one surface of the airfoil.

        :param side: which surface to locate the point on.
        :param method: how ``value`` is interpreted, see ``AfPos``.
        :param value: the position value. Positive is measured from the leading edge, negative
            from the trailing edge.
        :return: the station on the requested surface, or ``None`` when the position does not land
            on it (a camber distance longer than the camber line, a radius whose circle does not
            reach the surface, or an offset whose orthogonal cast misses it).
        :raises ValueError: if ``side`` or ``method`` is not a valid token.
        """
        ...

    def thickness_at(self, method: AfPos, value: float) -> Optional[Distance2]:
        """
        Measure the thickness of the airfoil at a location specified by one of the gage point
        location methods.

        A corresponding pair of points is located on the lower and upper surfaces with the same
        ``method`` and ``value``, and the thickness is the Euclidean distance between them. The
        returned measurement runs from the lower point to the upper point.

        :param method: how ``value`` is interpreted, see ``AfPos``.
        :param value: the position value. Positive is measured from the leading edge, negative
            from the trailing edge.
        :return: the thickness measurement, or ``None`` when the position does not land on both
            surfaces.
        :raises ValueError: if ``method`` is not a valid token.
        """
        ...

    def max_thickness(self) -> Distance2:
        """
        Measure the maximum thickness of the airfoil, taken from the largest inscribed circle.

        The measurement runs between the two contact points of ``max_thickness_circle``, which are
        points measured on the section itself.

        Note that this is **not** the maximum of ``thickness_at``. Sweeping that with
        ``"on_camber"`` measures a chord orthogonal to the camber, which is not required to fit
        inside the section, so on a cambered airfoil its maximum is larger and at a different
        station. This method reports the inscribed circle, which is the conventional definition of
        maximum airfoil thickness.
        """
        ...

    @property
    def circle_array(self) -> NDArray[np.float64]:
        """
        The inscribed circle stack as an ``(N, 3)`` numpy array whose columns are the circle
        center ``x``, the circle center ``y``, and the circle radius.
        """
        ...
