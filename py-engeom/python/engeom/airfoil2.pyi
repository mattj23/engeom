from __future__ import annotations

from typing import List

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
