from __future__ import annotations

from typing import Iterable, Tuple, TypeVar, Iterator, Any

from numpy.typing import NDArray
from engeom.engeom import ResampleEnum

from engeom import geom3

Transformable2 = TypeVar("Transformable2", Vector2, Point2, Iso2, SurfacePoint2)
PointOrVec2 = TypeVar("PointOrVec2", Point2, Vector2)


class Vector2(Iterable[float]):
    """
    A class representing a vector in 2D space. The vector contains an x and y component.  It is iterable and will
    yield the x and y components in order, allowing the Python unpacking operator `*` to be used to compensate for the
    lack of function overloading through some other parts of the library.

    A vector has different semantics than a point when it comes to transformations and some mathematical operations.
    """

    def __iter__(self) -> Iterator[float]:
        pass

    def __init__(self, x: float, y: float):
        """
        Create a 2D vector from the given x and y components.
        :param x: the x component of the vector.
        :param y: the y component of the vector.
        """
        ...

    @property
    def x(self) -> float:
        """
        Access the x component of the vector as a floating point value.
        """
        ...

    @property
    def y(self) -> float:
        """
        Access the y component of the vector as a floating point value.
        """
        ...

    def __rmul__(self, other: float) -> Vector2:
        """
        Multiply the vector by a scalar value. This allows the scalar to be on the left side of the multiplication
        operator.
        :param other: a scalar value to multiply the vector by.
        :return: a new vector that is the result of the multiplication.
        """
        ...

    def __mul__(self, other: float) -> Vector2:
        """
        Multiply the vector by a scalar value.
        :param other:  a scalar value to multiply the vector by.
        :return: a new vector that is the result of the multiplication.
        """
        ...

    def __add__(self, other: PointOrVec2) -> PointOrVec2:
        """
        Add a vector to a point or another vector. Adding a vector to a point will return a new point, while
        adding a vector to a vector will return a new vector.
        :param other: a point or vector to add to the vector.
        :return: a new point or vector that is the result of the addition.
        """
        ...

    def __sub__(self, other: Vector2) -> Vector2:
        """
        Subtract a vector from this vector.
        :param other: the vector to subtract from this vector.
        :return: a new vector that is the result of the subtraction.
        """
        ...

    def __neg__(self) -> Vector2:
        """
        Invert the vector by negating the x and y components.
        :return: a new vector in which the x and y components are negated.
        """
        ...

    def __truediv__(self, other: float) -> Vector2:
        """
        Divide the vector by a scalar value.
        :param other: a scalar value to divide the vector by.
        :return: a new vector that is the result of the division.
        """
        ...

    def as_numpy(self) -> NDArray[float]:
        """
        Create a numpy array of shape (2, ) from the vector.
        """
        ...

    def dot(self, other: Vector2) -> float:
        """
        Compute the dot product of two vectors. The result is a scalar value.
        :param other: the vector to compute the dot product with.
        :return: the scalar dot product of the two vectors.
        """
        ...

    def cross(self, other: Vector2) -> float:
        """
        Compute the cross product of two vectors.
        :param other: the vector to compute the cross product with.
        :return: the scalar cross product of the two vectors.
        """
        ...

    def norm(self) -> float:
        """
        Compute the Euclidian norm (aka magnitude, length) of the vector.
        """
        ...

    def normalized(self) -> Vector2:
        """
        Return a normalized version of the vector. The normalized vector will have the same direction as the original
        vector, but with a magnitude of 1.
        """
        ...

    def angle_to(self, other: Vector2) -> float:
        """
        Compute the smallest angle between two vectors and return it in radians.

        :param other: the vector to compute the angle to.
        :return: the angle between the two vectors in radians.
        """
        ...

    def with_x(self, x: float) -> Vector2:
        """
        Return a new vector with the same y component as this vector, but with the x component set to the given value.
        :param x: the new x component of the vector.
        :return: a new vector with the same y component as this vector, but with the x component set to the given value.
        """
        ...

    def with_y(self, y: float) -> Vector2:
        """
        Return a new vector with the same x component as this vector, but with the y component set to the given value.
        :param y: the new y component of the vector.
        :return: a new vector with the same x component as this vector, but with the y component set to the given value.
        """
        ...


class Point2(Iterable[float]):
    """
    A class representing a point in 2D space. The point contains an x and y component. It is iterable and will yield
    the x and y components in order, allowing the Python unpacking operator `*` to be used to compensate for the lack
    of function overloading through some other parts of the library.

    A point has different semantics than a vector when it comes to transformations and some mathematical operations.
    """

    def __iter__(self) -> Iterator[float]:
        pass

    def __init__(self, x: float, y: float):
        """
        Create a 2D point from the given x and y components.
        :param x: the x component of the point.
        :param y: the y component of the point.
        """
        ...

    @property
    def x(self) -> float:
        """
        Access the x component of the point as a floating point value.
        """
        ...

    @property
    def y(self) -> float:
        """
        Access the y component of the point as a floating point value.
        """
        ...

    @property
    def coords(self) -> Vector2:
        """
        Get the coordinates of the point as a `Vector2` object.
        :return: a `Vector2` object with the same x and y components as the point.
        """
        ...

    def __sub__(self, other: PointOrVec2) -> PointOrVec2:
        """
        Subtract a point or vector from this point. Subtracting a point from a point will return a new vector, while
        subtracting a vector from a point will return a new point.
        :param other: a point or vector to subtract from the point.
        :return: a new point or vector that is the result of the subtraction.
        """
        ...

    def __add__(self, other: Vector2) -> Vector2:
        """
        Add a vector to this point.
        :param other: the vector to add to the point.
        :return: a new point that is the result of the addition.
        """
        ...

    def __mul__(self, other: float) -> Point2:
        """
        Multiply the point's x and y components by a scalar value, returning a new point.
        :param other: the scalar value to multiply the point by.
        :return: a new point that is the result of the multiplication.
        """
        ...

    def __truediv__(self, other) -> Point2:
        """
        Divide the point's x and y components by a scalar value, returning a new point.
        :param other: the scalar value to divide the point by.
        :return: a new point that is the result of the division.
        """
        ...

    def __rmul__(self, other) -> Point2:
        """
        Multiply the point's x and y components by a scalar value, returning a new point. This allows the scalar to be
        on the left side of the multiplication.
        :param other: the scalar value to multiply the point by.
        :return: a new point that is the result of the multiplication.
        """
        ...

    def __neg__(self) -> Point2:
        """
        Invert the point by negating the x and y components.
        :return: a new point in which the x and y components are negated.
        """
        ...

    def as_numpy(self) -> NDArray[float]:
        """
        Create a numpy array of shape (2, ) from the point.
        """
        ...

    @staticmethod
    def mid(a: Point2, b: Point2) -> Point2:
        """
        Return the midpoint between two points. This is the average of the x and y components of the two points.
        """
        ...

    def with_x(self, x: float) -> Point2:
        """
        Return a new point with the same y component as this point, but with the x component set to the given value.
        :param x: the new x component of the point.
        :return: a new point with the same y component as this point, but with the x component set to the given value.
        """
        ...

    def with_y(self, y: float) -> Point2:
        """
        Return a new point with the same x component as this point, but with the y component set to the given value.
        :param y: the new y component of the point.
        :return: a new point with the same x component as this point, but with the y component set to the given value.
        """
        ...


class SurfacePoint2:
    """
    This class is used to represent a surface point in 2D space.

    Surface points are a composite structure that consist of a point in space and a normal direction. Conceptually, they
    come from metrology as a means of representing a point on the surface of an object along with the normal direction
    of the surface at that point. However, they are also isomorphic with the concept of a ray or a parameterized line
    with a direction of unit length, and can be used in that way as well.
    """

    def __init__(self, x: float, y: float, nx: float, ny: float):
        """
        Create a surface point from the given x and y components and the normal vector components. The normal vector
        components will be normalized automatically upon creation.  If the normal vector is the zero vector, an
        exception will be thrown.

        :param x: the x component of the point.
        :param y: the y component of the point.
        :param nx: the x component of the normal vector.
        :param ny: the y component of the normal vector.
        """
        ...

    @property
    def point(self) -> Point2:
        """
        Get the coordinates of the point as a Point2 object.
        :return: a Point2 object
        """
        ...

    @property
    def normal(self) -> Vector2:
        """
        Get the normal of the point as a Vector2 object.
        :return: a Vector2 object
        """
        ...

    def at_distance(self, distance: float) -> Point2:
        """
        Get the point at a distance along the normal from the surface point.
        :param distance: the distance to move along the normal.
        :return: the point at the distance along the normal.
        """
        ...

    def scalar_projection(self, point: Point2) -> float:
        """
        Calculate the scalar projection of a point onto the axis defined by the surface point position and direction.
        Positive values indicate that the point is in the normal direction from the surface point, while negative values
        indicate that the point is in the opposite direction.

        :param point: the point to calculate the projection of.
        :return: the scalar projection of the point onto the normal.
        """
        ...

    def projection(self, point: Point2) -> Point2:
        """
        Calculate the projection of a point onto the axis defined by the surface point position and direction.

        :param point: the point to calculate the projection of.
        :return: the projection of the point onto the plane.
        """
        ...

    def reversed(self) -> SurfacePoint2:
        """
        Return a new surface point with the normal vector inverted, but the position unchanged.
        :return: a new surface point with the inverted normal vector.
        """
        ...

    def planar_distance(self, point: Point2) -> float:
        """
        Calculate the planar (non-normal) distance between the surface point and a point. This is complementary to the
        scalar projection. A point is projected onto the plane defined by the position and normal of the surface point,
        and the distance between the surface point position and the projected point is returned.  The value will always
        be positive.

        :param point: the point to calculate the distance to.
        :return: the planar distance between the surface point and the point.
        """
        ...

    def shift_orthogonal(self, distance: float) -> SurfacePoint2:
        """
        Shift the surface point by a distance orthogonal to the normal vector. The direction of travel is the surface
        point's normal vector rotated 90 degrees clockwise. For instance, if the normal vector is (0, 1), a positive
        distance will move the point to the right and a negative distance will move the point to the left.

        :param distance: the distance to shift the surface point.
        :return: a new surface point shifted by the given distance.
        """
        ...

    def rot_normal(self, angle: float) -> SurfacePoint2:
        """
        Rotate the normal vector of the surface point by a given angle in radians and return a new surface point. The
        position of the surface point is not affected. The angle is positive for counter-clockwise rotation and negative
        for clockwise rotation.

        :param angle: the angle to rotate the normal vector by.
        :return: a new surface point with the rotated normal vector.
        """

    def __mul__(self, other: float) -> SurfacePoint2:
        """
        Multiply the position of the surface point by a scalar value. The normal vector is not affected unless the
        scalar is negative, in which case the normal vector is inverted.
        :param other: the scalar value to multiply the position by.
        :return: a new surface point with the position multiplied by the scalar.
        """
        ...

    def __rmul__(self, other: float) -> SurfacePoint2:
        """
        Multiply the position of the surface point by a scalar value. The normal vector is not affected unless the
        scalar is negative, in which case the normal vector is inverted.
        :param other: the scalar value to multiply the position by.
        :return: a new surface point with the position multiplied by the scalar.
        """
        ...

    def __truediv__(self, other: float) -> SurfacePoint2:
        """
        Divide the position of the surface point by a scalar value. The normal vector is not affected unless the
        scalar is negative, in which case the normal vector is inverted.
        :param other: the scalar value to divide the position by.
        :return: a new surface point with the position divided by the scalar.
        """
        ...

    def __neg__(self) -> SurfacePoint2:
        """
        Invert both the position AND the normal vector of the surface point.
        """
        ...

    def offset(self, offset: Vector2) -> SurfacePoint2:
        """
        Offset the surface point by a given vector. The normal vector is not affected.
        :param offset: the vector to offset the surface point by.
        :return: a new surface point with the position offset by the given vector.
        """
        ...

    def shift(self, distance: float) -> SurfacePoint2:
        """
        Shift the surface point by a given distance along the normal vector. The position of the surface point is
        affected, but the normal vector is not.
        :param distance: the distance to shift the surface point.
        :return: a new surface point with the position shifted by the given distance.
        """
        ...


class Iso2:
    """
    A class representing an isometry in 2D space. An isometry is a transformation that preserves distances and angles,
    also sometimes known as a rigid body transformation. It is composed of a translation and a rotation.

    `Iso2` objects can be used to transform points, vectors, surface points, other isometries, and a number of other
    2D geometric constructs.
    """

    def __init__(self, tx: float, ty: float, r: float):
        """
        Create an isometry from a translation and a rotation. The translation is represented by the x and y components
        of the translation vector. The rotation is represented by the angle in radians, and will be a rotation around
        the origin of the coordinate system.

        In convention with typical transformation matrices, transforming by an isometry constructed this way is the
        equivalent of first rotating by the angle `r` and then translating by the vector `(tx, ty)`.

        :param tx: the x component of the translation vector.
        :param ty: the y component of the translation vector.
        :param r: the angle of rotation in radians around the origin, where a positive value is a counter-clockwise
        rotation.
        """
        ...

    @staticmethod
    def identity() -> Iso2:
        """
        Create the identity isometry.
        """
        ...

    def __matmul__(self, other: Transformable2) -> Transformable2:
        """
        Transform a point, vector, or other transformable object by the isometry using the matrix multiplication
        operator. The transform must be on the right side of the operator, and the object being transformed must be on
        the left side. This is the equivalent of multiplying the object by the isometry matrix.

        When composing multiple isometries together, remember that the order of operations is reversed. For example, if
        you have isometries A, B, and C, and you want to compose them together such that they are the equivalent of
        first applying A, then B, then C, you would write `D = C @ B @ A`.

        :param other: the object to transform.
        :return: an object of the same type as the input, transformed by the isometry.
        """
        ...

    def inverse(self) -> Iso2:
        """
        Get the inverse of the isometry, which is the isometry that undoes the transformation of the original isometry,
        or the isometry that when composed with the original isometry produces the identity isometry.
        """
        ...

    def as_numpy(self) -> NDArray[float]:
        """
        Create a numpy array of shape (3, 3) from the isometry.
        """
        ...

    def transform_points(self, points: NDArray[float]) -> NDArray[float]:
        """
        Transform an array of points using the isometry. The semantics of transforming points are such that the full
        matrix is applied, first rotating the point around the origin and then translating it by the translation vector.

        To transform vectors, use the `transform_vectors` method instead.

        This is an efficient way to transform a large number of points at once, rather than using the `@` operator
        individually on a large number of `Point2` objects.

        :param points: a numpy array of shape (N, 2)
        :return: a numpy array of shape (N, 2) containing the transformed points in the same order as the input.
        """
        ...

    def transform_vectors(self, vectors: NDArray[float]) -> NDArray[float]:
        """
        Transform an array of vectors using the isometry. The semantics of transforming vectors are such that only the
        rotation matrix is applied, and the translation vector is not used. The vectors retain their original
        magnitude, but their direction is rotated by the isometry.

        To transform points, use the `transform_points` method instead.

        This is an efficient way to transform a large number of vectors at once, rather than using the `@` operator
        individually on a large number of `Vector2` objects.

        :param vectors: a numpy array of shape (N, 2)
        :return: a numpy array of shape (N, 2) containing the transformed vectors in the same order as the input.
        """
        ...


class SvdBasis2:
    """
    A class which creates a set of orthonormal basis vectors from a set of points in 2D space. The basis is created
    using a singular value decomposition of the points, and is very similar to the statistical concept of principal
    component analysis.

    The basis can be used to determine the rank of the point set, the variance of the points along the basis vectors,
    and to extract an isometry that will transform points from the world space to the basis space.  It is useful for
    orienting unknown point sets in a consistent way, for finding best-fit lines or planes, and for other similar
    tasks.
    """

    def __init__(
            self,
            points: NDArray[float],
            weights: NDArray[float] | None = None
    ):
        """
        Create a basis from a set of points. The basis will be calculated using a singular value decomposition of the
        points.

        :param points: a numpy array of shape (n, 2) containing the points to calculate the basis from.
        :param weights: a numpy array of shape (n,) containing the weights of the points. If None, all points will be
        weighted equally.
        """
        ...

    def rank(self, tol: float) -> int:
        """
        Retrieve the rank of the decomposition by counting the number of singular values that are
        greater than the provided tolerance.  A rank of 0 indicates that all singular values are
        less than the tolerance, and thus the point set is essentially a single point. A rank of 1
        indicates that the point set is essentially a line. A rank of 2 indicates that the point
        set exists roughly in a plane.

        The singular values do not directly have a clear physical meaning. They are square roots of
        the variance multiplied by the number of points used to compute the basis.  Thus, they can
        be interpreted in relation to each other, and when they are very small.

        This method should be used either when you know roughly what a cutoff tolerance for the
        problem you're working on should be, or when you know the cutoff value should be very
        small.  Otherwise, consider examining the standard deviations of the basis vectors
        instead, as they will be easier to interpret (`basis_stdevs()`).
        :param tol: the tolerance to use when determining the rank.
        :return: the rank of the decomposition.
        """
        ...

    def largest(self) -> Vector2:
        """
        Get the largest singular vector of the basis.
        :return: the largest singular vector.
        """
        ...

    def smallest(self) -> Vector2:
        """
        Get the smallest singular vector of the basis.
        :return: the smallest singular vector.
        """
        ...

    def basis_variances(self) -> NDArray[float]:
        """
        Get the variance of the points along the singular vectors.
        :return: a numpy array of the variance of the points along the singular vectors.
        """
        ...

    def basis_stdevs(self) -> NDArray[float]:
        """
        Get the standard deviation of the points along the singular vectors.
        :return: a numpy array of the standard deviation of the points along the singular vectors.
        """
        ...

    def to_iso2(self) -> Iso2:
        """
        Produce an isometry which will transform from the world space to the basis space.

        For example, if the basis is created from a set of points that lie roughly on an arbitrary line, multiplying
        original points by this isometry will move the points such that all points are aligned with the x-axis.
        :return: the isometry that transforms from the world space to the basis space.
        """
        ...


class CurveStation2:
    """
    A class representing a station along a curve in 2D space. The station is represented by a point on the curve, a
    tangent (direction) vector, and a length along the curve.

    These are created as the result of position finding operations on `Curve2` objects.
    """

    @property
    def point(self) -> Point2:
        """
        Get the point in 2D world space where the station is located.
        :return: the point in 2D world space.
        """
        ...

    @property
    def direction(self) -> Vector2:
        """
        Get the direction vector of the curve at the location of the station. This is the tangent vector of the curve,
        and is typically the direction from the previous vertex to the next vertex.
        :return: the direction vector of the curve at the station.
        """
        ...

    @property
    def normal(self) -> Vector2:
        """
        Get the normal vector of the curve at the location of the station. This is the vector that is orthogonal to the
        direction vector, and is the direction vector at the station rotated by -90 degrees. When the curve represents
        a manifold surface, this vector represents the direction of the surface normal.
        :return: the surface normal vector of the curve at the station.
        """
        ...

    @property
    def direction_point(self) -> SurfacePoint2:
        """
        Get the combined point and direction vector of the curve at the location of the station, returned as a
        `SurfacePoint2` object.
        :return: the combined point and direction vector of the curve at the station.
        """
        ...

    @property
    def surface_point(self) -> SurfacePoint2:
        """
        Get the combined point and normal vector of the curve at the location of the station, returned as a
        `SurfacePoint2` object.
        :return: the combined point and normal vector of the curve at the station.
        """
        ...

    @property
    def index(self) -> int:
        """
        Get the index of the previous vertex on the curve, at or before the station.
        :return: the index of the previous vertex on the curve.
        """
        ...

    @property
    def length_along(self) -> float:
        """
        Get the length along the curve to the station, starting at the first vertex of the curve.
        :return: the length along the curve to the station.
        """
        ...


class Curve2:
    """
    A class representing a curve in 2D space. The curve is defined by a set of vertices and the line segments between
    them (also known as a polyline).

    Because the curve is in 2D space, it also has a concept of a surface normal direction, which is orthogonal to the
    tangent direction of the curve at any point. This normal direction allows a `Curve2` to represent a 2D manifold
    surface boundary, defining the concepts of inside and outside.  It is commonly used to represent the surface of a
    solid body in a 2D cross-section.

    Additionally, the `Curve2` object can be used to represent closed regions by connecting the first and last vertices
    and allowing the curve to be treated as a closed loop. This lets the `Curve2` also represent closed polygons.
    """

    def __init__(
            self,
            vertices: NDArray[float],
            normals: NDArray[float] | None = None,
            tol: float = 1e-6,
            force_closed: bool = False,
            hull_ccw: bool = False,
    ):
        """
        Create a 2d curve from a set of vertices and some additional options.

        It's important to note that in 2d, a curve has a concept of a normal direction, built from the concept of
        inside/outside defined through the winding order of the vertices. This extra information can allow a 2d curve
        to model a manifold surface.

        There are three ways to specify the winding order of the vertices:

        1. Control it manually by passing the vertices array with the rows already organized so that an exterior surface
        is counter-clockwise.

        2. If the vertices represent an exterior shape, pass `hull_ccw=True` to have the constructor automatically
        check the winding order and reverse it if point ordering in the convex hull does not match ordering in the
        original array.

        3. Pass a `normals` array the same size as the `vertices` array, where the normals are non-zero vectors pointed
        in the "outside" direction at each point. The constructor will reverse the winding if the majority of normals
        do not point in the same direction as the winding.

        :param vertices: a numpy array of shape (N, 2) representing the vertices of the curve.
        :param normals: an optional numpy array of shape (N, 2) representing the normals of the curve associated with
        each vertex.
        :param tol: a tolerance value for the curve. If not provided, a default value of 1e-6 is used. This is the
        distance at which two points are considered to be the same.
        :param force_closed: If True, the curve will be closed even if the first and last points are not the same, which
        will be done by adding a new point at the end of the array that is the same as the first point.
        :param hull_ccw: If True, the constructor will check the winding order of the vertices and reverse it if the
        convex hull of the points is not in the same order as the original array. This will do nothing if the `normals`
        parameter is provided.
        """
        ...

    def length(self) -> float:
        """
        Get the total length of the curve as a scalar value.
        :return: the length of the curve.
        """
        ...

    def at_front(self) -> CurveStation2:
        """
        Get the station at the front of the curve.
        :return: the station at the front of the curve.
        """
        ...

    def at_back(self) -> CurveStation2:
        """
        Get the station at the back of the curve.
        :return: the station at the back of the curve.
        """
        ...

    def at_length(self, length: float) -> CurveStation2:
        """
        Get the station at a given length along the curve. Will throw a ValueError if the length is less than zero or
        greater than the length of the curve.
        :param length: the length along the curve.
        :return: the station at the given length.
        """
        ...

    def at_fraction(self, fraction: float) -> CurveStation2:
        """
        Get the station at a given fraction of the length of the curve. Will throw a ValueError if the fraction is less
        than zero or greater than one.
        :param fraction: the fraction of the length of the curve.
        :return: the station at the given fraction.
        """
        ...

    def at_closest_to_point(self, point: Point2) -> CurveStation2:
        """
        Get the station on the curve that is closest to a given point.
        :param point: the point to find the closest station to.
        :return: the station on the curve that is closest to the given point.
        """
        ...

    @property
    def is_closed(self) -> bool:
        """
        Check if the curve is closed.
        :return: True if the curve is closed, False otherwise.
        """
        ...

    def trim_front(self, length: float) -> Curve2:
        """
        Remove the front of the curve by a given length and return a new curve.
        :param length: the length to trim from the front of the curve.
        :return: a new curve with the front trimmed by the given length.
        """
        ...

    def trim_back(self, length: float) -> Curve2:
        """
        Remove the back of the curve by a given length and return a new curve.
        :param length: the length to trim from the back of the curve.
        :return: a new curve with the back trimmed by the given length.
        """
        ...

    def between_lengths(self, l0: float, l1: float) -> Curve2:
        """
        Attempt to get a new curve cut between two lengths along the curve. If the lengths are not valid, a ValueError
        will be thrown.

        If the curve is closed, the lengths will be wrapped around the curve. If the curve is not closed, the value
        of `l0` must be less than `l1`. In either case, the lengths must be within the bounds of the curve.

        :param l0: the start length.
        :param l1: the end length.
        :return: a new curve between the two lengths.
        """
        ...

    def between_lengths_by_control(self, a: float, b: float, control: float) -> Curve2:
        """
        Attempt to get a new curve cut between two lengths along the curve, with a control point that will be used to
        determine which side of the curve to keep. This is primarily helpful on closed curves when you can find a length
        (usually via use of the `at_closest_to_point` method) that is on the side of the curve you want to keep.

        If the lengths are not valid, a ValueError will be thrown.

        :param a: the first length along the curve to cut
        :param b: the second length along the curve to cut
        :param control: a length along the curve that is on a point in the portion of the result that you want to keep
        :return: a new curve between the two lengths
        """

    def reversed(self) -> Curve2:
        """
        Reverse the curve and return a new curve.
        :return: a new curve with the vertices in reverse order.
        """
        ...

    def make_hull(self) -> NDArray[int]:
        """
        Get the vertices of a convex hull of the curve, in counter-clockwise order.
        :return: a numpy array of shape (N, 2) representing the convex hull of the curve.
        """
        ...

    def max_point_in_direction(self, direction: Vector2) -> Tuple[int, Point2]:
        """
        Find the point on the curve that is furthest in a given direction.
        :param direction: the direction to find the furthest point in.
        :return: a tuple of the index of the point and the point itself.
        """
        ...

    def max_distance_in_direction(self, surf_point: SurfacePoint2) -> float:
        """
        Find the maximum scalar projection of all vertices of the curve onto a surface point.
        :param surf_point: the direction to find the furthest point in.
        :return: the maximum scalar projection of all vertices of the curve onto a surface point.
        """
        ...

    @property
    def points(self) -> NDArray[float]:
        """
        Get the points of the curve.
        :return: a numpy array of shape (N, 2) representing the points of the curve.
        """
        ...

    def simplify(self, tol: float) -> Curve2:
        """
        Simplify the curve using the Ramer-Douglas-Peucker algorithm.
        :param tol: the tolerance to use for simplification.
        :return: a new curve with the simplified points.
        """
        ...

    def resample(self, resample: ResampleEnum) -> Curve2:
        """
        Resample the curve using the given resampling method. The resampling method can be one of the following:

        - `Resample.ByCount(count: int)`: resample the curve to have the given number of points.
        - `Resample.BySpacing(distance: float)`: resample the curve to have points spaced by the given distance.
        - `Resample.ByMaxSpacing(distance: float)`: resample the curve to have points spaced by a maximum distance.

        :param resample: the resampling method to use.
        :return: a new curve object with the resampled vertices.
        """
        ...

    def transformed_by(self, transform: Iso2) -> Curve2:
        """
        Transform the curve by the given transform and return a new curve.
        :param transform: the transform to apply to the curve.
        :return: a new curve object with the transformed vertices.
        """
        ...

    def to_3d(self) -> geom3.Curve3:
        """
        Convert the curve to a 3D curve by adding a z-coordinate of 0 to all points.
        :return: a new `Curve3` object representing the curve in 3D space.
        """
        ...

    @property
    def aabb(self) -> Aabb2:
        """
        Get the axis-aligned bounding box of the curve.
        :return: the axis-aligned bounding box of the curve.
        """
        ...

    def offset_vertices(self, offset: float) -> Curve2:
        """
        Create a new curve which is the result of offsetting the vertices of this curve by the
        given offset. The direction of each vertex offset will be the same as the direction of the
        surface normal at the curve station corresponding to that vertex, which is the angle
        bisecting the normals of the two edges that meet at the vertex.  Vertices at the ends of
        the curve (on an open curve) will have the same normal as the edge they are connected to.

        Compared to `offset_segments`, this method will move the vertices of the curve while
        allowing the distance between the bodies of the initial and resulting segments to change.
        Generally speaking, use this method if you primarily care about the vertices and not the
        segments, or if the curvature between adjacent segments is very low.

        :param offset: the distance to offset the vertices by.
        :return: a new curve with the vertices offset by the given distance.
        """
        ...

    def offset_segments(self, offset: float) -> Curve2:
        """
        Create a new curve which is the result of offsetting the segments of this curve by the
        given offset. The direction of the offset is perpendicular to the direction of the segment,
        and a positive offset will move the segment outward from the curve, while a negative offset
        will move it inward.  Outward and inward are defined based on the counter-clockwise winding
        convention.

        Vertices will be moved to the intersection of their adjacent segments.

        Compared to `offset_vertices`, this method will preserve the distance between the segments
        bodies of the initial and resulting curves, while allowing vertices on outside corners to
        get farther from the original as necessary for the segments to be straight lines.

        :param offset: the distance to offset the segments by.
        :return: a new curve with the segments offset by the given distance.
        """
        ...

    def __add__(self, other: Curve2) -> Curve2:
        """
        Concatenate two curves together, returning a new curve that is the result of appending the vertices of the
        second curve to the first curve. Both curves must be open or this will throw an error. The resulting curve
        will be open.

        :param other: the curve to append to this curve.
        :return: a new curve that is the result of concatenating the two curves.
        """
        ...



class Circle2:
    """
    A class representing a circle in 2D space. The circle is defined by a center point and a radius.
    """

    def __init__(self, x: float, y: float, r: float):
        """
        Create a circle from the given center point and radius.
        :param x: the x-coordinate of the center of the circle.
        :param y: the y-coordinate of the center of the circle.
        :param r: the radius of the circle.
        """
        ...

    @property
    def center(self) -> Point2:
        """
        Get the `Point2` at the center of the circle.
        :return: the center of the circle.
        """
        ...

    @property
    def x(self) -> float:
        """
        Get the x-coordinate of the circle.
        :return: the x-coordinate of the circle.
        """
        ...

    @property
    def y(self) -> float:
        """
        Get the y-coordinate of the circle.
        :return: the y-coordinate of the circle.
        """
        ...

    @property
    def r(self) -> float:
        """
        Get the radius of the circle.
        :return: the radius of the circle.
        """
        ...

    @property
    def aabb(self) -> Aabb2:
        """
        Get the axis-aligned bounding box of the circle.
        :return: the axis-aligned bounding box of the circle.
        """
        ...

    def point_at_angle(self, angle: float) -> Point2:
        """
        Get the point on the circle at a given angle.
        :param angle: the angle in radians.
        :return: the point on the circle at the given angle.
        """
        ...

    @staticmethod
    def fitting(points: NDArray[float], guess: Circle2 | None = None, sigma: float | None = None) -> Circle2:
        """
        Fit a circle to a set of points using an unconstrained Levenberg-Marquardt minimization of the sum of
        squared errors between the points and the boundary of the circle.

        The initial guess is used to provide a starting point for the optimization. If no guess is provided, the
        unit circle will be used.

        The sigma parameter is used to weight the points in the optimization. If no sigma is provided, all points
        will be weighted equally, otherwise points beyond `sigma` standard deviations from the mean will be
        assigned a weight of 0.0.
        :param points: the points to fit the circle to.
        :param guess: an optional initial guess for the circle. If None, the unit circle will be used.
        :param sigma: an optional standard deviation to use for weighting the points. If None, all points will be
        weighted equally.
        :return: a new `Circle2` object representing the fitted circle.
        """
        ...

    @staticmethod
    def ransac(points: NDArray[float], tol: float, iterations: int | None = None, min_r: float | None = None,
               max_r: float | None = None) -> Circle2:
        """
        Fit a circle to a set of points using the RANSAC algorithm. The algorithm will randomly sample points from the
        input set and fit a circle to them, then check how many points are within the given tolerance of the fitted
        circle. The best fitting circle will be returned.

        :param points: the points to fit the circle to.
        :param tol: the tolerance for the RANSAC algorithm.
        :param iterations: the number of iterations to run. If None, a default value of 500 will be used.
        :param min_r: the minimum radius of the circle. If None, no minimum will be enforced.
        :param max_r: the maximum radius of the circle. If None, no maximum will be enforced.
        :return: a new `Circle2` object representing the fitted circle.
        """
        ...


class Arc2:
    """
    An arc in 2D space. The arc is defined by a center point, a radius, a start angle, and a sweep angle.

    * The center point and the radius define the circle of which the arc is part.

    * The start angle is the angle in radians from the positive x-axis to the point where the arc begins. A positive
      value is a counter-clockwise rotation, so a start angle of $\\pi / 2$ would start the arc at the top $y=r$ of the
      circle.

    * The sweep angle is the angle in radians that the arc covers, beginning at the starting point. A positive value is
      a counter-clockwise rotation, a negative value is clockwise.
    """

    def __init__(self, x: float, y: float, r: float, start_radians: float, sweep_radians: float):
        """
        Create an arc from the given center point, radius, start angle, and sweep angle.

        :param x: the x-coordinate of the center of the arc.
        :param y: the y-coordinate of the center of the arc.
        :param r: the radius of the arc.
        :param start_radians: the start angle of the arc in radians, which is the angle from the positive x-axis to the
        starting point of the arc. A positive value is a counter-clockwise rotation.
        :param sweep_radians: the sweep angle of the arc in radians, which is the angle that the arc covers, beginning
        at the starting point. A positive value is a counter-clockwise rotation, a negative value is clockwise.
        """

    @property
    def center(self) -> Point2:
        """
        Get the center point of the arc.
        :return: the center of the arc.
        """
        ...

    @property
    def x(self) -> float:
        """
        Get the x-coordinate of the center of the arc.
        :return: the x-coordinate of the arc center.
        """
        ...

    @property
    def y(self) -> float:
        """
        Get the y-coordinate of the center of the arc.
        :return: the y-coordinate of the arc center
        """
        ...

    @property
    def r(self) -> float:
        """
        Get the radius of the arc
        :return: the radius of the arc
        """
        ...

    @property
    def start(self) -> float:
        """
        Get the start angle of the arc, in radians.
        :return: the start angle of the arc in radians.
        """
        ...

    @property
    def sweep(self) -> float:
        """
        Get the sweep angle of the arc, in radians.
        :return: the sweep angle of the arc in radians.
        """
        ...

    @property
    def aabb(self) -> Aabb2:
        """
        Get the axis-aligned bounding box of the arc.
        :return: the axis-aligned bounding box of the arc.
        """
        ...

    @property
    def start_point(self) -> Point2:
        """
        Get the start point of the arc.
        :return: the start point of the arc.
        """
        ...

    @property
    def end_point(self) -> Point2:
        """
        Get the end point of the arc.
        :return: the end point of the arc.
        """
        ...


class Aabb2:
    """
    A class representing an axis-aligned bounding box in 2D space. The bounding box is defined by a minimum point and a
    maximum point, which are the lower-left and upper-right corners of the box, respectively.

    Bounding boxes are typically used for accelerating intersection and distance queries and are used internally inside
    the Rust language `engeom` library for this purpose.  However, they have other useful applications and so are
    exposed here in the Python API.

    Typically, `Aabb2` objects will be retrieved from other `engeom` objects which use them internally, such as curves,
    circles, arcs, etc.  However, they can also be created and manipulated directly.
    """

    def __init__(self, x_min: float, y_min: float, x_max: float, y_max: float):
        """
        Create an axis-aligned bounding box from the minimum and maximum coordinates.

        :param x_min: the minimum x-coordinate of the AABB
        :param y_min: the minimum y-coordinate of the AABB
        :param x_max: the maximum x-coordinate of the AABB
        :param y_max: the maximum y-coordinate of the AABB
        """
        ...

    @staticmethod
    def at_point(x: float, y: float, w: float, h: float | None = None) -> Aabb2:
        """
        Create an AABB centered at a point with a given width and height.
        :param x: the x-coordinate of the center of the AABB.
        :param y: the y-coordinate of the center of the AABB.
        :param w: the width of the AABB.
        :param h: the height of the AABB. If not provided, the AABB will be square.
        :return: a new AABB object.
        """
        ...

    @staticmethod
    def from_points(points: NDArray[float]) -> Aabb2:
        """
        Create an AABB that bounds a set of points. If the point array is empty or the wrong shape, an error will be
        thrown.
        :param points: a numpy array of shape (N, 2) containing the points to bound
        :return: a new AABB object
        """
        ...

    @property
    def min(self) -> Point2:
        """
        Get the minimum point of the AABB.
        :return: the minimum point of the AABB.
        """
        ...

    @property
    def max(self) -> Point2:
        """
        Get the maximum point of the AABB.
        :return: the maximum point of the AABB.
        """
        ...

    @property
    def center(self) -> Point2:
        """
        Get the center point of the AABB.
        :return: the center point of the AABB.
        """
        ...

    @property
    def extent(self) -> Vector2:
        """
        Get the extent of the AABB.
        :return: the extent of the AABB.
        """
        ...

    def expand(self, d: float) -> Aabb2:
        """
        Expand the AABB by a given distance in all directions. The resulting height and
        width will be increased by 2 * d.

        :param d: the distance to expand the AABB by.
        :return: a new AABB object with the expanded bounds.
        """
        ...

    def shrink(self, d: float) -> Aabb2:
        """
        Shrink the AABB by a given distance in all directions. The resulting height and
        width will be decreased by 2 * d.

        :param d: the distance to shrink the AABB by.
        :return: a new AABB object with the shrunk bounds.
        """
        ...

    def merged(self, other: Aabb2) -> Aabb2:
        """
        Merge this AABB with another AABB and return a new AABB.
        :param other: the other AABB to merge with.
        :return: a new AABB object that is the result of merging this AABB with the other AABB.
        """
        ...

    def indices_contained(self, points: NDArray[float]) -> NDArray[int]:
        """
        Get the indices of the points that are contained within the AABB.
        :param points: a numpy array of shape (N, 2) containing the points to check.
        :return: a numpy array of indices of the points that are contained within the AABB.
        """
        ...

    def contains_point(self, point: Point2) -> bool:
        """
        Check if a point is contained within the AABB.
        :param point: the point to check.
        :return: True if the point is contained within the AABB, False otherwise.
        """
        ...
