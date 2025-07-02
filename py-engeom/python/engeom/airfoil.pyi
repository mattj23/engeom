from typing import List

from numpy.typing import NDArray
from enum import Enum
import geom2

type MclOrientEnum = MclOrient.TmaxFwd | MclOrient.DirFwd
type FaceOrientEnum = FaceOrient.Detect | FaceOrient.UpperDir
type EdgeFindEnum = EdgeFind.Open | EdgeFind.OpenIntersect | EdgeFind.Intersect | EdgeFind.RansacRadius
type AfGageEnum = AfGage.OnCamber | AfGage.Radius


class EdgeType(Enum):
    """
    An enumeration of the possible non-geometric types of edges that can be detected on an airfoil cross-section as the
    result of an edge finding operation. When one of these types is returned, it means that the edge finding algorithm
    could not provide more detailed geometric information about the edge.
    """

    Open = 0
    """ Represents an open edge, where the airfoil cross-section is incomplete and/or not closed. """

    Closed = 1
    """ Represents an edge that is closed, where the airfoil cross-section has contiguous vertices through the edge. """


class AfGage:
    """
    An enumeration class used to specify a method for locating a gage points on an airfoil cross-section.
    """

    class OnCamber:
        """
        A gaging method that measures a distance along the mean camber line. A positive distance will be from the
        leading edge towards the trailing edge, and a negative distance will be from the trailing edge towards the
        leading edge.
        """

        def __init__(self, d: float):
            """
            Create a specification for a gage point that is a distance `d` along the mean camber line. A positive
            distance will be from the leading edge towards the trailing edge, and a negative distance will be from the
            trailing edge towards the leading edge.
            :param d: the distance along the mean camber line to find the position
            """
            ...

    class Radius:
        """
        A gaging method that measures by intersection with a circle of a given radius centered on either the
        leading or trailing edge point.  A positive radius indicates that the circle is located on the leading edge
        while a negative radius indicates that the circle is located on the trailing edge.
        """

        def __init__(self, r: float):
            """
            Create a specification for a gage point that is located at the intersection of a circle of radius `r` with
            the airfoil cross-section. A positive radius indicates that the circle is located on the leading edge while
            a negative radius indicates that the circle is located on the trailing edge.
            :param r: the radius of the circle to find the position
            """
            ...


class FaceOrient:
    """
    An enumeration of the possible ways to orient the upper/lower (suction/pressure, convex/concave) faces of an
    airfoil cross-section.
    """

    class Detect:
        """
        In an airfoil with an MCL that exhibits curvature, this will attempt to detect which direction the camber line
        curves and thus identify convex/concave. This will fail if the MCL is straight.
        """

        def __init__(self):
            """
            Create a new face orientation parameter that will attempt to detect the orientation of the faces based on
            the curvature of the mean camber line.
            """
            ...

    class UpperDir:
        """
        This method will orient the faces based on a vector direction provided by the user.
        """

        def __init__(self, x: float, y: float):
            """
            Create a new upper direction parameter. The x and y arguments are components of a direction vector which
            should distinguish the upper (pressure side, convex) face of the airfoil. At the center of the mean camber
            line, an intersection in this direction will be taken with each of the two faces. The intersection that
            is further in the direction of this vector will be considered the upper face of the airfoil, and the other
            will be considered the lower face.

            :param x: the x component of the upper direction vector
            :param y: the y component of the upper direction vector
            """
            ...


class MclOrient:
    """
    An enumeration of the possible ways to orient (to identify which side is the leading edge and which side is the
    trailing edge) the mean camber line of an airfoil.
    """

    class TmaxFwd:
        """
        This method will take advantage of the fact that for most typical subsonic airfoils the maximum thickness point
        is closer to the leading edge than the trailing edge.
        """

        def __init__(self):
            """
            Create a specification for orienting the mean camber line based on which side the maximum thickness point is
            closer to. This method will assume that the maximum thickness point is closer to the leading edge than the
            trailing edge.
            """
            ...

    class DirFwd:
        """
        This method will orient the airfoil based on a vector direction provided by the user.
        """

        def __init__(self, x: float, y: float):
            """
            Create a new forward direction parameter. The x and y arguments are components of a direction vector which
            should distinguish the forward (leading edge) direction of the airfoil. The position of the first and last
            inscribed circle will be projected onto this vector, and the larger result (the one that is more in the
            direction of this vector) will be considered the leading edge of the airfoil.

            For instance, if you know that the airfoil is oriented so that the leading edge will have a smaller x value
            than the trailing edge, `DirFwd(-1, 0)` will correctly orient the airfoil.
            :param x: the x component of the forward direction vector
            :param y: the y component of the forward direction vector
            """
            ...


class EdgeFind:
    """
    An enumeration of the possible techniques to find the leading and/or trailing edge geometry of an airfoil.
    """

    class Open:
        """
        This algorithm will not attempt to find edge geometry, and will simply leave the inscribed circles for the side
        as they are. Use this if you know that the airfoil cross-section is open/incomplete on this side, and you don't
        care to extend the MCL any further.
        """

        def __init__(self):
            """
            Create a specification which assumes that the airfoil cross-section is open/incomplete on this side, and
            makes no attempt to find the edge geometry beyond the unambiguous inscribed circles.
            """
            ...

    class OpenIntersect:
        """
        This algorithm is also for an open edge, but unlike `Open` it will attempt to refine the end of the MCL and
        extend it to intersect the line segment which spans the open gap where the edge should be. This is useful on
        partial cross-sections where you would still like to extend the MCL as much as possible.

        It works by intersecting the end of the inscribed circles camber curve with the open gap in the airfoil
        cross-section, filling and refining more inscribed circles between the last circle and the intersection point,
        and repeating until the location of the end converges to within 1/100th of the general refinement tolerance.
        """

        def __init__(self, max_iter: int):
            """
            This algorithm will attempt to find the edge geometry by intersecting the end of the inscribed circles
            camber curve with the open gap in the airfoil cross-section, then refining the end of the MCL with more
            inscribed circles until the location of the end converges to within 1/100th of the general refinement
            tolerance.

            If the maximum number of iterations is reached before convergence, the method will throw an error instead.

            :param max_iter: the maximum number of iterations to attempt to find the edge geometry
            """
            ...

    class Intersect:
        """
        This algorithm will simply intersect the end of the inscribed circles camber curve with the airfoil
        cross-section. This is the fastest method with the least amount of assumptions, and makes sense for airfoil
        edges where you know the mean camber line has very low curvature in the vicinity of the edge.

        Do not use this method if you know that the airfoil cross-section is open/incomplete on this side, as it will
        throw an error if the MCL does not intersect the cross-section.
        """

        def __init__(self):
            """
            Create a specification which will attempt to find the edge geometry by intersecting the end of the inscribed
            circles camber curve with the airfoil cross-section.  Use on known closed airfoil cross-sections with low
            curvature near the edge.
            """
            ...

    class RansacRadius:
        """
        This technique uses RANSAC (Random Sample Consensus) to find a constant radius leading/trailing edge circle
        that fits the greatest number of points leftover at the edge within the tolerance `in_tol`.

        The method will try `n` different combinations of three points picked at random from the remaining points
        at the edge, construct a circle, and then count the number of points within `in_tol` distance of the circle
        perimeter. The circle with the most points within tolerance will be considered the last inscribed circle.

        The MCL will be extended to this final circle, and then intersected with the airfoil cross-section to find
        the final edge point.
        """

        def __init__(self, in_tol: float, n: int = 500):
            """
            Create a specification which will attempt to find the edge geometry by Random Sample Consensus of a constant
            radius at the edge of the airfoil cross-section. This is useful for airfoils known to have a constant radius
            edge and on section data which is relatively clean and has low noise.

            :param in_tol: the max distance from the circle perimeter for a point to be considered a RANSAC inlier
            :param n: The number of RANSAC iterations to perform
            """
            ...


class InscribedCircle:
    """
    Represents an inscribed circle in an airfoil cross-section. The circle is contained within the airfoil cross-section
    and is tangent to the airfoil section at two points.
    """

    from .geom2 import Circle2, Point2

    @property
    def circle(self) -> Circle2:
        """
        Gets the circle object associated with this inscribed circle.
        :return: The circle entity for the inscribed circle
        """
        ...

    @property
    def contact_a(self) -> Point2:
        """
        Get a contact point of the inscribed circle with one side of the airfoil cross-section. Inscribed circles
        computed together will have a consistent meaning of `a` and `b` sides, but which is the upper or lower surface
        will depend on the ordering of the circles and the coordinate system of the airfoil.

        :return: The first contact point of the inscribed circle with the airfoil cross-section
        """
        ...

    @property
    def contact_b(self) -> Point2:
        """
        Get the other contact point of the inscribed circle with the airfoil cross-section. Inscribed circles computed
        together will have a consistent meaning of `a` and `b` sides, but which is the upper or lower surface will
        depend on the ordering of the circles and the coordinate system of the airfoil.

        :return: The second contact point of the inscribed circle with the airfoil cross-section
        """
        ...


class EdgeResult:
    """
    Represents the results of an airfoil edge detection operation, containing both a point on the airfoil cross-section
    that was detected as the edge, and optional geometric information about the edge depending on the method used.
    """
    from .geom2 import Arc2, Point2

    @property
    def point(self) -> Point2:
        """
        Gets the point on the airfoil cross-section that was detected as the edge.
        :return: The point on the airfoil cross-section that was detected as the edge
        """
        ...

    @property
    def geometry(self) -> EdgeType | Arc2:
        """
        Gets the geometric information about the edge that was detected.

        * This will be an instance of `EdgeType` if the algorithm could not provide more detailed geometric information
        about the edge beyond open/closed.

        * This will be an instance of `Arc2` in the case of constant radius edge detection.

        :return: The geometric information about the edge that was detected
        """
        ...


class AirfoilGeometry:
    """
    This class produces and contains the result of the geometric analysis of an airfoil cross-section. It contains:

    1. The mean camber line

    2. The segregated geometry of the upper (suction, convex) and lower (pressure, concave) faces

    3. The detected position and geometry of the leading and trailing edges

    4. The array of inscribed circles identified during the analysis

    From this information, thicknesses and other measurements can be made on the airfoil.

    See the `from_analyze` method for how to perform the analysis creating an instance of this class.
    """

    from .geom2 import Curve2, Circle2
    from .metrology import Distance2

    @staticmethod
    def from_bytes(data: bytes) -> AirfoilGeometry:
        """
        Create an instance of `AirfoilGeometry` from a byte string containing the airfoil geometry data serialized in
        a msgpack format.
        :param data: The byte string containing the serialized data
        :return: An instance of `AirfoilGeometry`
        """
        ...

    def to_bytes(self) -> bytes:
        """
        Serialize the `AirfoilGeometry` instance to a byte string msgpack representation.
        :return: A byte string containing the serialized data
        """
        ...

    @staticmethod
    def from_analyze(
            section: Curve2,
            refine_tol: float,
            camber_orient: MclOrientEnum,
            leading: EdgeFindEnum,
            trailing: EdgeFindEnum,
            face_orient: FaceOrientEnum,
    ) -> AirfoilGeometry:
        """
        This method attempts to extract the airfoil geometry from the given airfoil cross-section using only the
        geometric information embedded in the section. It is suitable for airfoil cross-sections with clean, low-noise
        or noise-free data, such as those which come from nominal CAD data or from very clean scans/samples of airfoils
        with smooth, continuous surfaces with little to no defects.

        The cross-section data must also be *only* the outer surface of the airfoil, with no internal features or
        points. The vertices should be ordered in a counter-clockwise direction, but the section may be open on one
        side and no particular orientation or position in the XY plane is required.

        Internally, this operation will attempt to:

        1. Extract the unambiguous inscribed circles from the cross-section using an iterative stepping/refinement
        method (see `compute_inscribed_circles` for more detail).

        2. Orient the mean camber line (determine which side is the leading edge and which side is the trailing edge)
        based on the method specified in `camber_orient`.

        3. Extract the exact position (and optionally, the geometry) of the leading and trailing edges based on the
        methods specified in `leading` and `trailing`, respectively.

        4. Identify and extract the upper (suction, convex) and lower (pressure, concave) faces of the airfoil based on
        the method specified in `face_orient`.

        If successful in all of these steps, the method will return an instance of `AirfoilGeometry` containing the
        results of the analysis. If any of the steps fail, the method will raise an error.

        :param section: The curve representing the airfoil cross-section
        :param refine_tol: A general tolerance used in the analysis, typically used to refine results until the error
        or difference between two values is less than this value. It is also used in certain methods as a common
        reference tolerance, where the method will use a fraction of this value as a threshold for convergence or error.
        :param camber_orient: The method to use to orient the mean camber line of the airfoil.
        :param leading: The method to use to detect the leading edge of the airfoil.
        :param trailing: The method to use to detect the trailing edge of the airfoil.
        :param face_orient: The method to identify the upper and lower faces of the airfoil.
        :return: An instance of `AirfoilGeometry` containing the results of the analysis.
        """
        ...

    @property
    def leading(self) -> EdgeResult | None:
        """
        Gets the result of the leading edge detection algorithm.
        """
        ...

    @property
    def trailing(self) -> EdgeResult | None:
        """
        Gets the result of the trailing edge detection algorithm.
        """
        ...

    @property
    def camber(self) -> Curve2:
        """
        Gets the mean camber line of the airfoil cross-section. The curve will be oriented so that the first point is at
        the leading edge of the airfoil and the last point is at the trailing edge.
        :return: The mean camber line of the airfoil cross-section as a `Curve2` object.
        """
        ...

    @property
    def upper(self) -> Curve2 | None:
        """
        The curve representing the upper (suction, convex) side of the airfoil cross-section. The curve will be oriented
        in the same winding direction as the original section, so the first point may be at either the leading or
        trailing edge based on the airfoil geometry and the coordinate system.

        :return: A `Curve2`, or None if there was an issue detecting the leading or trailing edge.
        """
        ...

    @property
    def lower(self) -> Curve2 | None:
        """
        The curve representing the lower (pressure, concave) side of the airfoil cross-section. The curve will be
        oriented in the same winding direction as the original section, so the first point may be at either the leading
        or trailing edge based on the airfoil geometry and the coordinate system.

        :return: A `Curve2`, or None if there was an issue detecting the leading or trailing edge.
        """
        ...

    @property
    def circle_array(self) -> NDArray[float]:
        """
        Returns the list of inscribed circles as a numpy array of shape (N, 3) where N is the number of inscribed
        circles. The first two columns are the x and y coordinates of the circle center, and the third column is the
        radius of the circle.
        :return: A numpy array of shape (N, 3) containing the inscribed circles
        """
        ...

    def get_thickness(self, gage: AfGageEnum) -> Distance2:
        """
        Get the thickness dimension of the airfoil cross-section. The 'a' point of the distance will be on the lower
        surface of the airfoil, and the 'b' point will be on the upper surface.
        :param gage: the gaging method to use
        :return: a `Distance2` object representing the thickness dimension at the specified gage points
        """
        ...

    def get_tmax(self) -> Distance2:
        """
        Get the maximum thickness dimension of the airfoil cross-section. The 'a' point of the distance will be on the
        lower surface of the airfoil, and the 'b' point will be on the upper surface.
        :return: a `Distance2` object representing the maximum thickness dimension
        """
        ...

    def get_tmax_circle(self) -> Circle2:
        """
        Get the circle representing the maximum thickness dimension of the airfoil cross-section.
        :return: a `Circle2` object representing the maximum thickness circle
        """
        ...


def compute_inscribed_circles(section: geom2.Curve2, refine_tol: float) -> List[InscribedCircle]:
    """
    Compute the unambiguous inscribed circles of an airfoil cross-section.

    The cross-section is represented by a curve in the x-y plane. The curve does not need to be closed, but the points
    should be oriented in a counter-clockwise direction and should only contain data from the outer surface of the
    airfoil (internal features/points should not be part of the data).

    The method used to compute these circles is:

    1. We calculate the convex hull of the points in the section and find the longest distance between any two points.
    2. At the center of the longest distance line, we draw a perpendicular line and look for exactly two intersections
       with the section. We assume that one of these is on the upper surface of the airfoil and the other is on the
       lower, though it does not matter which is which.
    3. We fit the maximum inscribed circle whose center is constrained to the line between these two points. The
       location and radius of this circle is refined until it converges to within 1/100th of `refine_tol`.
    4. The inscribed circle has two contact points with the section. The line between these contact points is a good
       approximation of the direction orthogonal to the mean camber line near the circle.  We create a parallel line
       to this one, advancing from the circle center by 1/4 of the circle radius, and looking for exactly two
       intersections with the section.  If we fail, we try again with a slightly less aggressive advancement until we
       either succeed or give up.
    5. We fit the maximum inscribed circle whose center is constrained to the new line, and refine it as in step 3.
    6. We recursively fit inscribed circles between this new circle and the previous one until the error between the
       position and radius of any circle is less than `refine_tol` from the linear interpolation between its next and
       previous neighbors.
    7. We repeat the process from step 4 until the distance between the center of the most recent circle and the
       farthest point in the direction of the next advancement is less than 1/4 of the radius of the most recent
       circle. This terminates the process before we get too close to the leading or trailing edge of the airfoil.
    8. We repeat the process from step 3, but this time in the opposite direction from the first circle. This will
       give us the inscribed circles on the other side of the airfoil.

    When finished, we have a list of inscribed circles from the unambiguous regions (not too close to the leading or
    trailing edges) of the airfoil cross-section. The circles are ordered from one side of the airfoil to the other,
    but the order may be *either* from the leading to the trailing edge *or* vice versa.

    :param section: the curve representing the airfoil cross-section.
    :param refine_tol: a tolerance used when refining the inscribed circles, see description for details.
    :return: a list of inscribed circle objects whose order is contiguous but may be in either direction
    """
    ...
