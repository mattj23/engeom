from __future__ import annotations

from typing import Callable, Iterable, Tuple, TypeVar, Iterator, Any, List

from numpy.typing import NDArray
from engeom.engeom import ResampleEnum, VecDot
from engeom.common import AngleDir, AngleInterval

from engeom import geom3

Transformable2 = TypeVar("Transformable2", Vector2, Point2, Iso2, SurfacePoint2, Segment2, Circle2)
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

    @staticmethod
    def zero() -> Vector2:
        """ Return the zero vector (0, 0). """
        ...

    @staticmethod
    def x_axis() -> Vector2:
        """ Return the unit vector along the X axis (1, 0). """
        ...

    @staticmethod
    def y_axis() -> Vector2:
        """ Return the unit vector along the Y axis (0, 1). """
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

    @staticmethod
    def origin() -> Point2:
        """ Return the origin point (0, 0). """
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

    def shifted_orthogonal(self, distance: float) -> SurfacePoint2:
        """
        Return a new surface point shifted by a distance orthogonal to the normal vector. The direction of travel is the
        surface point's normal vector rotated 90 degrees clockwise. For instance, if the normal vector is (0, 1), a
        positive distance will move the point to the right and a negative distance will move the point to the left.

        :param distance: the distance to shift the surface point.
        :return: a new surface point shifted by the given distance.
        """
        ...

    def normal_rotated(self, angle: float) -> SurfacePoint2:
        """
        Return a new surface point with its normal vector rotated by a given angle in radians. The position of the
        surface point is not affected. The angle is positive for counter-clockwise rotation and negative for clockwise
        rotation.

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

    def transformed_by(self, iso: Iso2) -> SurfacePoint2:
        """
        Transform the surface point by an isometry, moving its position and rotating its normal.
        :param iso: the isometry to apply to the surface point.
        :return: a new surface point transformed by the given isometry.
        """
        ...

    def shifted(self, distance: float) -> SurfacePoint2:
        """
        Shift the surface point by a given distance along the normal vector. The position of the surface point is
        affected, but the normal vector is not.
        :param distance: the distance to shift the surface point.
        :return: a new surface point with the position shifted by the given distance.
        """
        ...

    def to_line(self) -> Line2:
        """
        Convert this surface point to a ``Line2`` whose origin is the surface point's position and
        whose direction is the surface point's normal vector.

        :return: a ``Line2`` through this surface point.
        """
        ...

    def slerp(self, other: SurfacePoint2, t: float) -> SurfacePoint2:
        """
        Spherically interpolate between this surface point and another by parameter `t`. The
        position is linearly interpolated and the normal is spherically interpolated so that it
        remains unit length.

        :param other: the surface point to interpolate towards.
        :param t: the interpolation parameter, where 0.0 returns this point and 1.0 returns `other`.
        :return: a new surface point interpolated between this one and `other`.
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

    @staticmethod
    def from_array(matrix: NDArray[float]) -> Iso2:
        """
        Try to create an isometry from a 3x3 homogeneous transformation matrix. Raises a `ValueError` if the matrix
        does not represent a valid isometry (no shear or scale).
        :param matrix: a numpy array of shape (3, 3) containing the matrix values.
        :return: the isometry represented by the matrix.
        """
        ...

    @staticmethod
    def from_rotation_about(point: Point2, angle: float) -> Iso2:
        """
        Create an isometry representing a rotation around an arbitrary point by a given angle, rather than around
        the origin as with the `Iso2` constructor. Positive angles rotate counter-clockwise.
        :param point: the point in the plane that the isometry will rotate around.
        :param angle: the angle to rotate by, in radians.
        :return: the isometry representing the rotation.
        """
        ...

    @staticmethod
    def from_basis_x(e0: Vector2, origin: Point2 | None = None) -> Iso2:
        """
        Try to create an isometry from a single basis vector and an optional origin. The vector becomes the x-axis
        of the isometry (after being normalized to unit length), and the y-axis is the x-axis rotated 90 degrees
        counter-clockwise. Raises a `ValueError` if the vector cannot be normalized.

        Unlike the 3D `from_basis_xy`/`from_basis_xz`/etc. family, a single vector fully determines a 2D isometry's
        rotation, since there is only one rotational degree of freedom.
        :param e0: the vector that will become the x-axis, will be normalized automatically.
        :param origin: an optional origin point for the isometry, defaults to the world origin if omitted.
        :return: the isometry defined by the basis vector and origin.
        """
        ...

    @staticmethod
    def from_basis_y(e1: Vector2, origin: Point2 | None = None) -> Iso2:
        """
        Try to create an isometry from a single basis vector and an optional origin. The vector becomes the y-axis
        of the isometry (after being normalized to unit length), and the x-axis is the y-axis rotated 90 degrees
        clockwise. Raises a `ValueError` if the vector cannot be normalized.

        Unlike the 3D `from_basis_xy`/`from_basis_xz`/etc. family, a single vector fully determines a 2D isometry's
        rotation, since there is only one rotational degree of freedom.
        :param e1: the vector that will become the y-axis, will be normalized automatically.
        :param origin: an optional origin point for the isometry, defaults to the world origin if omitted.
        :return: the isometry defined by the basis vector and origin.
        """
        ...

    def flipped(self) -> Iso2:
        """
        Return a copy of the isometry rotated 180 degrees in-plane around its own origin.

        This is a rotation, not a mirror/reflection: the origin's location is unchanged, but both the x-axis and
        y-axis directions are reversed together. A 2D isometry cannot represent a reflection (determinant -1), only
        proper rotations, and a 180 degree turn about the implicit out-of-plane axis is the only one that flips both
        in-plane axes at once.
        """
        ...

    @property
    def origin(self) -> Point2:
        """
        The world-space location of the isometry's origin.
        :return: a Point2 at the isometry's translation.
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

    @property
    def x_direction(self) -> Vector2:
        """
        The world-space direction of the isometry's local x-axis.
        :return: a unit Vector2 along the transformed x-axis.
        """
        ...

    @property
    def y_direction(self) -> Vector2:
        """
        The world-space direction of the isometry's local y-axis.
        :return: a unit Vector2 along the transformed y-axis.
        """
        ...

    @property
    def x_axis(self) -> SurfacePoint2:
        """
        A SurfacePoint2 through the isometry's origin along its local x-axis.
        :return: a SurfacePoint2 with origin at the isometry's origin and direction along the x-axis.
        """
        ...

    @property
    def y_axis(self) -> SurfacePoint2:
        """
        A SurfacePoint2 through the isometry's origin along its local y-axis.
        :return: a SurfacePoint2 with origin at the isometry's origin and direction along the y-axis.
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
    def fraction(self) -> float:
        """
        Get the fractional position of the station between the previous vertex and the next vertex
        on the curve, in the range [0, 1].
        :return: the fraction between the previous and next vertex.
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

    def max_dist_in_direction(self, surf_point: SurfacePoint2) -> float:
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

    def new_transformed_by(self, transform: Iso2) -> Curve2:
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

    @staticmethod
    def load_tccurve2(path: str | Path) -> Curve2:
        """
        Load a curve from a tolerance-compressed 2D curve (.tccurve2) file. The tccurve2 format
        stores vertex positions as variable-width integers scaled within the bounding box of the
        point data, using the minimum number of bytes needed to guarantee a round-trip accuracy at
        or below the tolerance that was specified when the file was written. The curve's closed/open
        state and reconstruction tolerance are also stored in the file.

        :param path: the path to the .tccurve2 file to load.
        :return: the curve loaded from the file.
        """
        ...

    def write_tccurve2(self, path: str | Path, tol: float):
        """
        Write the curve to a tolerance-compressed 2D curve (.tccurve2) file. The tolerance controls
        the maximum allowable round-trip position error for any vertex: a smaller tolerance produces
        a more accurate file at the cost of more bytes per vertex, while a larger tolerance allows
        greater compression.

        :param path: the path to the .tccurve2 file to write.
        :param tol: the maximum acceptable round-trip position error, in model units.
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
    def from_fit(points: NDArray[float], weights: NDArray[float] | None = None) -> Circle2:
        """
        Fit a circle to a set of points by ordinary least squares. A closed-form algebraic
        (Kåsa-style) estimate provides the initial guess, which is then refined against the true
        geometric radial residuals with a weighted Levenberg-Marquardt minimization.

        This is not robust to gross outliers; for that, use `from_consensus`.
        :param points: the points to fit the circle to, as an (n, 2) array.
        :param weights: if provided, a length-n array of weights to multiply each point's residual
            by. If None, all points are weighted equally.
        :return: a new `Circle2` object representing the fitted circle.
        :raises ValueError: if there are fewer than three points, or they are collinear.
        """
        ...

    @staticmethod
    def from_consensus(points: NDArray[float], sigma_max: float, min_r: float | None = None,
                       max_r: float | None = None, max_iterations: int | None = None,
                       refinement_steps: int | None = None, confidence: float | None = None,
                       seed: int | None = None) -> Circle2:
        """
        Fit a circle to a set of points robustly using the MAGSAC++ consensus algorithm. Unlike a fixed-threshold
        RANSAC, this takes an upper bound on the inlier noise (`sigma_max`) rather than a hard inlier/outlier
        threshold, and refines each candidate with noise-marginalized iteratively reweighted least squares. It is
        substantially less sensitive to `sigma_max` than RANSAC is to its threshold, as long as `sigma_max` is not
        chosen smaller than the actual noise.

        :param points: the points to fit the circle to.
        :param sigma_max: the upper bound on the expected inlier noise, in the same units as the points.
        :param min_r: the minimum radius of the circle. If None, no minimum will be enforced.
        :param max_r: the maximum radius of the circle. If None, no maximum will be enforced.
        :param max_iterations: the maximum number of minimal-sample iterations. If None, a default of 500 is used.
        :param refinement_steps: the number of iteratively reweighted refinement steps per candidate. If None, a
            default of 4 is used.
        :param confidence: the probability used for adaptive termination. If None, a default of 0.99 is used.
        :param seed: an optional fixed RNG seed for reproducible sampling. If None, a random seed is used.
        :return: a new `Circle2` object representing the fitted circle.
        """
        ...

    def project_point_to_perimeter(self, point: Point2) -> Point2:
        """
        Project a point onto the perimeter of the circle.

        :param point: the point to project.
        :return: the projected point on the perimeter of the circle.
        """
        ...

    def angle_of_point(self, point: Point2) -> float:
        """
        Get the angle of a point on the perimeter of the circle, measured in radians from the positive x-axis, and
        where a positive angle is a counter-clockwise rotation.

        :param point: the point on the perimeter of the circle.
        :return: the angle in radians of the point on the perimeter of the circle.
        """
        ...

    def intersections_with(self, other: Circle2) -> List[Point2]:
        """
        Get the intersection points between this circle and another circle.

        :param other: the other circle to intersect with.
        :return: a list of intersection points. The list will contain 0, 1, or 2 points depending on the number of
        intersections.
        """
        ...

    def distance_to(self, point: Point2) -> float:
        """
        Get the distance from a point to the perimeter of the circle.

        :param point: the point to measure the distance from.
        :return: the distance from the point to the perimeter of the circle. The distance will be negative if the
        point is inside the circle.
        """
        ...

    def contains_point(self, x: float, y: float) -> bool:
        """
        Returns True if the point lies at or inside the boundary of the circle.
        :param x: the x component of the point to test.
        :param y: the y component of the point to test.
        :return: True if the point is inside or on the circle, False otherwise.
        """
        ...

    def tangent_points_to(self, point: Point2) -> List[Point2]:
        """
        Get the tangent points on the circle from a given point outside the circle.

        :param point: the point outside the circle.
        :return: a list of tangent points on the circle. The list will contain 0, 1, or 2 points depending on the
        position of the point relative to the circle.
        """
        ...

    @staticmethod
    def from_point(center: Point2, r: float) -> Circle2:
        """
        Create a circle from a center point and radius.

        :param center: the center of the circle.
        :param r: the radius of the circle.
        :return: a new ``Circle2``.
        """
        ...

    @staticmethod
    def from_3_points(p0: Point2, p1: Point2, p2: Point2) -> Circle2:
        """
        Create a circle passing through three points. Raises ``ValueError`` if the points are collinear.

        :param p0: the first point.
        :param p1: the second point.
        :param p2: the third point.
        :return: a new ``Circle2``.
        """
        ...

    @staticmethod
    def from_tangent_to_corner(corner: Point2, d0: Vector2, d1: Vector2, radius: float) -> Circle2:
        """
        Create a circle tangent to the corner formed by two lines. The corner is defined by a point
        and two direction vectors. Raises ``ValueError`` if the directions are collinear.

        :param corner: the corner point where the two lines meet.
        :param d0: direction vector of the first line.
        :param d1: direction vector of the second line.
        :param radius: radius of the tangent circle.
        :return: a new ``Circle2``.
        """
        ...

    @staticmethod
    def from_tangent_and_point(tangent: Line2, point: Point2) -> Circle2:
        """
        Create a circle tangent to a line and passing through a point.

        :param tangent: the line the circle is tangent to.
        :param point: a point on the circle.
        :return: a new ``Circle2``.
        """
        ...

    def outer_tangents_to(self, other: Circle2) -> Tuple[Segment2, Segment2] | None:
        """
        Compute the two outer tangent segments between this circle and another. Returns ``None`` if
        the circles are concentric.

        The first segment lies in the negative half-space of the line from this center to the other;
        the second lies in the positive half-space.

        :param other: the other circle.
        :return: a tuple of two ``Segment2`` objects, or ``None``.
        """
        ...

    def to_arc(self) -> Arc2:
        """
        Convert the full circle to an ``Arc2`` starting at angle 0 and spanning 2π.

        :return: a new ``Arc2``.
        """
        ...

    def to_partial_arc(self, angle0: float, angle: float) -> Arc2:
        """
        Create a partial arc of the circle.

        :param angle0: the starting angle in radians.
        :param angle: the arc span in radians.
        :return: a new ``Arc2``.
        """
        ...

    def line_direction(self, line: Line2) -> AngleDir:
        """
        Determine whether the circle center is on the clockwise (``"cw"``) or counter-clockwise
        (``"ccw"``) side of a line.

        :param line: the reference line.
        :return: ``"cw"`` or ``"ccw"``.
        """
        ...

    def at_angle(self, angle: float) -> Manifold1Pos2:
        """
        Return the manifold position (point, tangent direction, normal, and arc length) at the
        given angle around the circle.

        :param angle: the angle, in radians.
        :return: the manifold position at that angle.
        """
        ...

    def at_closest_to_point(self, p: Point2) -> Manifold1Pos2:
        """
        Project ``p`` onto the circle's perimeter and return the manifold position there.

        :param p: the point to project.
        :return: the manifold position of the closest point on the perimeter.
        """
        ...

    def intersection_interval(self, other: Circle2) -> AngleInterval | None:
        """
        Return the angular interval, measured on this circle, spanned by its intersection with
        ``other``. Returns ``None`` if the circles do not intersect.

        :param other: the other circle to intersect with.
        :return: the angular interval of the intersection, or ``None``.
        """
        ...

    def transformed_by(self, iso: Iso2) -> Circle2:
        """
        Create a copy of this circle with its center moved by an isometry. The radius is unchanged.

        :param iso: the isometry to transform the circle by.
        :return: a new ``Circle2``.
        """
        ...

    def resized_by(self, delta: float) -> Circle2:
        """
        Create a copy of this circle with ``delta`` added to its radius. The center is unchanged.

        A positive ``delta`` grows the circle and a negative one shrinks it. The result is not
        clamped, so a large enough negative ``delta`` will produce a zero or negative radius.

        :param delta: the amount to add to the circle's radius.
        :return: a new ``Circle2``.
        """
        ...


class Line2:
    """
    A parameterized line in 2D space: ``P(t) = origin + t * direction``.

    The direction is not required to be normalized. Use ``Line2.new_normalize`` or ``normalized()``
    for a unit-speed parameterization where the parameter ``t`` equals arc length from the origin.
    """

    def __init__(self, ox: float, oy: float, dx: float, dy: float):
        """
        Create a line from an origin point ``(ox, oy)`` and a direction vector ``(dx, dy)``.
        The direction is stored as-is and is **not** normalized automatically.

        :param ox: x-coordinate of the origin.
        :param oy: y-coordinate of the origin.
        :param dx: x-component of the direction vector.
        :param dy: y-component of the direction vector.
        """
        ...

    @staticmethod
    def x_axis() -> Line2:
        """Return the X axis: origin at (0, 0), direction (1, 0)."""
        ...

    @staticmethod
    def y_axis() -> Line2:
        """Return the Y axis: origin at (0, 0), direction (0, 1)."""
        ...

    @staticmethod
    def from_points(p1: Point2, p2: Point2) -> Line2:
        """
        Create a line through two points. The direction is ``p2 - p1`` (not normalized).

        :param p1: the origin point of the line.
        :param p2: a second point the line passes through.
        :return: a new ``Line2``.
        """
        ...

    @staticmethod
    def new_normalize(ox: float, oy: float, dx: float, dy: float) -> Line2:
        """
        Create a line from an origin and direction, normalizing the direction so that the
        parameter ``t`` equals arc length from the origin.

        :param ox: x-coordinate of the origin.
        :param oy: y-coordinate of the origin.
        :param dx: x-component of the direction vector (will be normalized).
        :param dy: y-component of the direction vector (will be normalized).
        :return: a new ``Line2`` with a unit-length direction.
        """
        ...

    @staticmethod
    def from_fit(points: NDArray[float], weights: NDArray[float] | None = None) -> Line2:
        """
        Fit a line to a set of points using singular value decomposition, resulting in a
        least-squares fitting. The resulting parameterized line will have its ``t=0`` sitting at
        the center of the SVD result.

        :param points: the points to fit the line to, as an (n, 2) array.
        :param weights: if provided, a length-n array of weights to multiply each point's
            residual by. If None, all points are weighted equally.
        :return: a new ``Line2`` fitted to the points.
        :raises ValueError: if the singular value decomposition fails (e.g. too few points).
        """
        ...

    @staticmethod
    def from_consensus(points: NDArray[float], sigma_max: float, max_iterations: int | None = None,
                       refinement_steps: int | None = None, confidence: float | None = None,
                       seed: int | None = None) -> Line2:
        """
        Fit a line to a set of points robustly using the MAGSAC++ consensus algorithm, rejecting gross outliers.
        Unlike a fixed-threshold RANSAC, this takes an upper bound on the inlier noise (`sigma_max`) rather than a
        hard inlier/outlier threshold, and refines each candidate with noise-marginalized iteratively reweighted
        least squares.

        :param points: the points to fit the line to.
        :param sigma_max: the upper bound on the expected inlier noise, in the same units as the points.
        :param max_iterations: the maximum number of minimal-sample iterations. If None, a default of 500 is used.
        :param refinement_steps: the number of iteratively reweighted refinement steps per candidate. If None, a
            default of 4 is used.
        :param confidence: the probability used for adaptive termination. If None, a default of 0.99 is used.
        :param seed: an optional fixed RNG seed for reproducible sampling. If None, a random seed is used.
        :return: a new ``Line2`` object representing the fitted line.
        """
        ...

    @property
    def origin(self) -> Point2:
        """The origin point of the line."""
        ...

    @property
    def direction(self) -> Vector2:
        """
        The direction vector of the line. May not be unit-length unless the line was created
        with ``new_normalize`` or ``normalized()``.
        """
        ...

    @property
    def normal(self) -> Vector2:
        """
        The unit normal to the line: the direction rotated 90 degrees clockwise. By convention
        this points to the right of the direction of travel, consistent with the outward normal
        on counter-clockwise-wound 2D geometry.
        """
        ...

    def at(self, t: float) -> Point2:
        """
        Evaluate the line at parameter ``t``: returns ``origin + t * direction``.

        :param t: the parameter value.
        :return: the point on the line at ``t``.
        """
        ...

    def scalar_project(self, point: Point2) -> float:
        """
        Return the parameter ``t`` such that ``at(t)`` is the closest point on the line to
        ``point``.  Equivalent to the scalar projection of ``(point - origin)`` onto the
        direction vector.

        :param point: the point to project.
        :return: the parameter value of the closest point on the line.
        """
        ...

    def closest_point(self, point: Point2) -> Point2:
        """
        Return the point on the line closest to ``point``.

        :param point: the query point.
        :return: the nearest point on the line.
        """
        ...

    def distance_to(self, point: Point2) -> float:
        """
        Return the perpendicular (unsigned) distance from ``point`` to the line.

        :param point: the query point.
        :return: the non-negative distance.
        """
        ...

    def signed_distance_to(self, point: Point2) -> float:
        """
        Return the signed perpendicular distance from ``point`` to the line. Positive means the
        point is to the right of the direction of travel (on the normal side); negative means it
        is to the left.

        :param point: the query point.
        :return: the signed distance.
        """
        ...

    def intersect(self, other: Line2) -> Point2 | None:
        """
        Return the intersection point with another line, or ``None`` if the lines are parallel.

        :param other: the other line to intersect with.
        :return: the intersection point, or ``None``.
        """
        ...

    def normalized(self) -> Line2:
        """
        Return a new line with the same origin but a normalized direction, so that the parameter
        ``t`` equals arc length from the origin.
        """
        ...

    def reversed(self) -> Line2:
        """
        Return a new line with the same origin, but with the direction inverted.
        """
        ...

    def offset_by(self, delta_n: float) -> Line2:
        """
        Return a new line parallel to this one, with the origin shifted by ``delta_n`` along the
        normal direction. A positive ``delta_n`` moves the line to the right of the direction of
        travel.

        :param delta_n: the offset distance along the normal.
        :return: a new parallel ``Line2``.
        """
        ...

    def rotated(self, angle: float) -> Line2:
        """
        Return a copy of this line rotated about its own origin by ``angle`` radians. The origin
        is unchanged and the direction is rotated counter-clockwise (a positive ``angle`` rotates
        from the +x axis toward the +y axis); the direction's magnitude is preserved.

        :param angle: the rotation angle in radians, counter-clockwise positive.
        :return: a new rotated ``Line2``.
        """
        ...

    def shifted_origin(self, delta_t: float) -> Line2:
        """
        Return a new line with the origin shifted by ``delta_t`` along the direction vector.

        :param delta_t: the distance to shift the origin along the direction.
        :return: a new ``Line2`` with the shifted origin.
        """
        ...

    def transformed_by(self, iso: Iso2) -> Line2:
        """
        Return a new line with both origin and direction transformed by the given isometry.

        :param iso: the isometry to apply.
        :return: a new transformed ``Line2``.
        """
        ...

    def slerp(self, other: Line2, t: float) -> Line2:
        """
        Return a new line whose origin and direction are spherically interpolated between this
        line and ``other`` by parameter ``t``.

        :param other: the line to interpolate towards.
        :param t: the interpolation parameter, where 0.0 returns this line and 1.0 returns ``other``.
        :return: a new interpolated ``Line2``.
        """
        ...

    def projected_parameter(self, p: Point2) -> float:
        """
        Return the parameter ``t`` at which ``p`` projects onto this line's direction, measured
        from the origin. Unlike ``scalar_project``, this projects onto the raw (possibly
        non-unit-length) direction vector.

        :param p: the point to project.
        :return: the projection parameter.
        """
        ...

    def projected_point(self, p: Point2) -> Point2:
        """
        Return the point on the line at the projected parameter of ``p``. Equivalent to
        ``self.at(self.projected_parameter(p))``.

        :param p: the point to project.
        :return: the projected point on the line.
        """
        ...

    def orthogonal(self) -> Vector2:
        """
        Return the direction vector rotated -90 degrees, typically used as a normal.
        """
        ...

    def signed_projection_dist(self, point: Point2) -> float:
        """
        Return the signed perpendicular distance from ``point`` to this line, measured against
        the ``orthogonal`` direction. Positive values indicate the point is to the right of the
        direction of travel; negative values indicate the left.

        :param point: the query point.
        :return: the signed distance.
        """
        ...

    def intersection_params(self, other: Line2) -> tuple[float, float] | None:
        """
        Return the pair of parameters ``(t_self, t_other)`` at which this line and ``other``
        intersect, or ``None`` if the two directions are parallel.

        :param other: the other line to intersect with.
        :return: the intersection parameters, or ``None``.
        """
        ...

    def winding_direction(self, point: Point2) -> AngleDir:
        """
        Determine the direction that the line winds around ``point``.

        :param point: the reference point.
        :return: ``"cw"`` or ``"ccw"``.
        """
        ...

    def to_iso_from_x(self) -> Iso2:
        """
        Return an isometry whose origin matches the line's origin and whose X direction matches
        the line's direction. Transforming an entity by the *inverse* of this isometry maps its
        relationship with the line back to the origin along the X axis.

        :return: the alignment isometry.
        """
        ...

    def to_iso_from_y(self) -> Iso2:
        """
        Return an isometry whose origin matches the line's origin and whose Y direction matches
        the line's direction. Transforming an entity by the *inverse* of this isometry maps its
        relationship with the line back to the origin along the Y axis.

        :return: the alignment isometry.
        """
        ...

    def to_surface_point(self) -> SurfacePoint2:
        """
        Convert this line to a ``SurfacePoint2`` whose position is the line's origin and whose
        normal vector is the line's (normalized) direction.

        :return: a ``SurfacePoint2`` at this line's origin.
        """
        ...


class CubicSpline2:
    """
    A cubic Bezier curve in 2D space, defined by four control points `p0`, `p1`, `p2`, `p3`.

    The curve is parameterized by a scalar `t`, conventionally in the range `[0, 1]`. At `t = 0`
    the curve passes through `p0` and at `t = 1` it passes through `p3`. The two interior control
    points influence the curve's shape but are not generally interpolated by it.
    """

    def __init__(
            self,
            x0: float,
            y0: float,
            x1: float,
            y1: float,
            x2: float,
            y2: float,
            x3: float,
            y3: float,
    ):
        """
        Create a cubic Bezier curve from the coordinates of its four control points, ordered from
        the start of the curve to the end.

        :param x0: x-coordinate of the first control point (curve start).
        :param y0: y-coordinate of the first control point (curve start).
        :param x1: x-coordinate of the second control point.
        :param y1: y-coordinate of the second control point.
        :param x2: x-coordinate of the third control point.
        :param y2: y-coordinate of the third control point.
        :param x3: x-coordinate of the fourth control point (curve end).
        :param y3: y-coordinate of the fourth control point (curve end).
        """
        ...

    @property
    def p0(self) -> Point2:
        """ The first control point (curve start). """
        ...

    @property
    def p1(self) -> Point2:
        """ The second control point. """
        ...

    @property
    def p2(self) -> Point2:
        """ The third control point. """
        ...

    @property
    def p3(self) -> Point2:
        """ The fourth control point (curve end). """
        ...

    def position(self, t: float) -> Point2:
        """
        Evaluate the curve position at parameter `t`. At `t = 0` this returns `p0`, and at
        `t = 1` it returns `p3`. Values outside `[0, 1]` will extrapolate the underlying
        polynomial.

        :param t: the curve parameter.
        :return: the point on the curve at parameter `t`.
        """
        ...

    def derivative(self, t: float) -> Vector2:
        """
        Evaluate the derivative of the curve at parameter `t`, returning the velocity vector
        (un-normalized). At `t = 0` this is `3 * (p1 - p0)`; at `t = 1` it is `3 * (p3 - p2)`.

        :param t: the curve parameter.
        :return: the derivative vector at parameter `t`.
        """
        ...

    def tangent(self, t: float) -> Vector2:
        """
        Evaluate the unit tangent vector of the curve at parameter `t`. This is the derivative
        normalized to unit length. The result is undefined (will contain NaN entries) at cusps
        where the derivative vanishes.

        :param t: the curve parameter.
        :return: the unit tangent vector at parameter `t`.
        """
        ...

    def normal(self, t: float) -> Vector2:
        """
        Evaluate the unit normal vector of the curve at parameter `t`. The normal is the unit
        tangent rotated 90 degrees clockwise.

        :param t: the curve parameter.
        :return: the unit normal vector at parameter `t`.
        """
        ...

    def curvature(self, t: float) -> float:
        """
        Evaluate the curvature magnitude at parameter `t`. The result is always non-negative;
        its reciprocal is the radius of the osculating circle. Returns NaN at a cusp.

        :param t: the curve parameter.
        :return: the curvature magnitude at parameter `t`.
        """
        ...

    def curvature_circle(self, t: float) -> Circle2 | None:
        """
        Return the circle of curvature (osculating circle) tangent to the spline at parameter `t`:
        the circle matching the curve's position, tangent direction, and curvature there. Its
        radius is the reciprocal of the curvature and its center is the center of curvature, on the
        concave side of the curve.

        Returns `None` where no finite circle exists: a locally straight section (zero curvature,
        infinite radius) or a cusp where the curvature is undefined.

        :param t: the curve parameter.
        :return: the osculating circle at parameter `t`, or `None`.
        """
        ...

    def find_max_curvature(self) -> Tuple[float, float]:
        """
        Find the point of maximum curvature on the curve over the parameter range `[0, 1]`,
        returned as a `(t, curvature)` tuple: the parameter `t` at which the maximum occurs and the
        curvature magnitude there. For a fully degenerate curve whose curvature is undefined
        everywhere, `t` is `0.0` and the curvature is NaN.

        :return: a `(t, curvature)` tuple at the point of maximum curvature.
        """
        ...

    def second_derivative(self, t: float) -> Vector2:
        """
        Evaluate the second derivative of the curve at parameter `t`, returning the acceleration
        vector (un-normalized). At `t = 0` this is `6 * (p0 - 2*p1 + p2)`; at `t = 1` it is
        `6 * (p1 - 2*p2 + p3)`.

        :param t: the curve parameter.
        :return: the second derivative vector at parameter `t`.
        """
        ...

    def polyline(self, tolerance: float) -> NDArray[float]:
        """
        Return an adaptive polyline approximation of the curve as an `Nx2` numpy array of points.
        The linear interpolation between consecutive returned points deviates from the underlying
        spline by no more than `tolerance` (measured as Euclidean distance). Regions where the
        curve is locally close to straight produce widely spaced points; regions of high curvature
        are subdivided more finely. The returned array always starts at `p0` and ends at `p3`.

        :param tolerance: maximum allowed Euclidean deviation between the polyline and the spline.
            Must be positive. Smaller values produce more points.
        :return: an `Nx2` array of the polyline vertices.
        """
        ...

    def query(self) -> CubicSplineQueries2:
        """
        Build a reusable acceleration structure for closest-point queries against this spline.
        Build the structure once and reuse it across many queries rather than rebuilding it per
        query.

        :return: a `CubicSplineQueries2` built from this spline.
        """
        ...

    def project_point(self, point: Point2) -> SplineProjection:
        """
        Find the closest point on the spline to the given point. This builds a temporary query
        structure on each call; for repeated queries against the same spline, build a `query()`
        once and reuse it.

        :param point: the point to project onto the spline.
        :return: a `SplineProjection` holding the closest-point parameter and distance.
        """
        ...

    def line_at(self, t: float) -> Line2:
        """
        Return the position and derivative direction of the curve at parameter `t` in the form of
        a parameterized line.

        :param t: the curve parameter.
        :return: a `Line2` through the curve's position at `t`, in the direction of the derivative.
        """
        ...

    def derivative_roots(self) -> Tuple[List[float], List[float]]:
        """
        Returns the real roots of the derivative of each component of the curve, as an
        `(x_roots, y_roots)` tuple where each is a list of 0, 1, or 2 parameter values, sorted
        ascending. Roots are not filtered to `[0, 1]`.

        :return: the per-axis derivative roots.
        """
        ...

    def find_cusp(self) -> float | None:
        """
        Returns the parameter `t` of a cusp if one exists in `[0, 1]`, otherwise `None`. A cusp is
        a point where the velocity vector vanishes (`B'(t) = 0`).

        :return: the parameter of the cusp, or `None`.
        """
        ...

    def find_curvature_zeros(self) -> List[float]:
        """
        Returns parameter values in `[0, 1]` where the curve's curvature is zero, as a list of 0,
        1, or 2 parameter values.

        :return: the curvature-zero parameters.
        """
        ...

    def find_curvature_maxima(self) -> List[Tuple[float, float]]:
        """
        Returns every local maximum of the curvature over `[0, 1]`, each as a `(t, curvature)`
        tuple, ordered by ascending `t`.

        :return: a list of `(t, curvature)` tuples, one per local maximum.
        """
        ...

    def compute_bounds(self) -> Tuple[Point2, Point2]:
        """
        Returns the corners `(min, max)` of the tight axis-aligned bounding box of the curve over
        the parameter range `[0, 1]`.

        :return: a `(min, max)` tuple of points.
        """
        ...

    def arc_length_between(self, t0: float, t1: float) -> float:
        """
        Returns the arc length of the curve over the parameter range `[t0, t1]`.

        :param t0: the start parameter.
        :param t1: the end parameter.
        :return: the arc length between the two parameters.
        """
        ...

    def arc_length(self) -> float:
        """
        Returns the total arc length of the curve over the parameter range `[0, 1]`.
        """
        ...

    def split(self, t: float) -> Tuple[CubicSpline2, CubicSpline2]:
        """
        Splits the curve at parameter `t` using de Casteljau's algorithm, returning the left and
        right sub-curves. Concatenating the left and right curves reproduces the original.

        :param t: the parameter at which to split, expected to be in `[0, 1]`.
        :return: a `(left, right)` tuple of the two sub-curves.
        """
        ...

    def try_split(self, t: float) -> Tuple[CubicSpline2, CubicSpline2] | None:
        """
        Splits the curve at parameter `t`, returning `None` if `t` is not in `[0, 1]`.

        :param t: the parameter at which to split.
        :return: a `(left, right)` tuple of the two sub-curves, or `None`.
        """
        ...

    @property
    def aabb(self) -> Aabb2:
        """
        Get the axis-aligned bounding box of the curve, computed on demand.
        :return: the axis-aligned bounding box of the curve.
        """
        ...

    def transformed_by(self, iso: Iso2) -> CubicSpline2:
        """
        Return a new curve with all four control points transformed by the given isometry.

        :param iso: the isometry to apply.
        :return: a new transformed `CubicSpline2`.
        """
        ...

    def find_inflections(self) -> List[float]:
        """
        Returns the parameter values of any inflection points of the curve, as a list of 0, 1, or
        2 parameter values. In 2D, an inflection point is where the signed curvature crosses zero.
        Roots are not filtered to `[0, 1]`.

        :return: the inflection parameters.
        """
        ...


class SplineProjection:
    """
    The result of a closest-point query against a cubic spline: the parameter `t` of the closest
    point on the curve and the distance from the queried point to that location. This type is
    shared by both the 2D and 3D spline queries.
    """

    @property
    def t(self) -> float:
        """ The parameter `t` of the closest point on the spline. Recover the point itself with
        `spline.position(t)`. """
        ...

    @property
    def distance(self) -> float:
        """ The distance from the queried point to the closest point on the spline. """
        ...


class CubicSplineQueries2:
    """
    A prebuilt acceleration structure for running repeated closest-point queries against a
    `CubicSpline2`. Build it once and reuse it across many queries.
    """

    def __init__(self, spline: CubicSpline2):
        """
        Build the query acceleration structure for the given spline.

        :param spline: the `CubicSpline2` to build queries against.
        """
        ...

    def project_point(self, point: Point2) -> SplineProjection:
        """
        Find the closest point on the spline to the given point.

        :param point: the point to project onto the spline.
        :return: a `SplineProjection` holding the closest-point parameter and distance.
        """
        ...

    def project_points(self, points: NDArray[float]) -> NDArray[float]:
        """
        Project an `Nx2` array of points onto the spline.

        :param points: an `Nx2` numpy array of points to project.
        :return: an `Nx2` array whose columns are the closest-point parameter `t` and the distance
            to the curve, one row per input point.
        """
        ...


def fit_spline_to_points(
        points: NDArray[float],
        builder: Callable[[NDArray[float]], CubicSpline2],
        initial: NDArray[float],
) -> NDArray[float]:
    """
    Fit a `CubicSpline2` to a set of 2D points using Levenberg-Marquardt optimization and a
    user-supplied builder function.

    Residuals are the unsigned distances from each point to its nearest location on the spline.
    There is no inherent pressure keeping the curve from sliding or growing past the points unless
    the parameterization forbids it, so write the builder so that the desired constraints are
    inherently contained within it.

    :param points: a numpy array of shape ``(N, 2)`` containing the sample points.
    :param builder: a callable that accepts a 1-D ``float64`` numpy array of parameters and returns
        a ``CubicSpline2``. Raise any exception to signal that the parameter vector produces invalid
        geometry; the optimizer treats that step as a failed evaluation.
    :param initial: a 1-D ``float64`` numpy array containing the initial parameter guess. The length
        determines the number of optimization parameters, and the builder must succeed on it.
    :return: a 1-D ``float64`` numpy array of the optimized parameters.
    :raises ValueError: if the initial guess causes the builder to fail, or if the optimizer does
        not converge.
    """
    ...


class Segment2:
    """
    A class representing a line segment in 2D space. The segment is defined by two endpoints.
    """

    def __init__(self, x0: float, y0: float, x1: float, y1: float):
        """
        Create a line segment from two endpoints.
        :param x0: the x-coordinate of the first endpoint.
        :param y0: the y-coordinate of the first endpoint.
        :param x1: the x-coordinate of the second endpoint.
        :param y1: the y-coordinate of the second endpoint.
        """
        ...

    @staticmethod
    def from_fit(points: NDArray[float], weights: NDArray[float] | None = None) -> Segment2:
        """
        Fit a segment to a set of points by ordinary least squares. An infinite line is fit to the points, and the
        segment's endpoints are set to the extreme projections of the points onto that line, so the segment spans
        exactly the range covered by the input.

        This is not robust to gross outliers; for that, use `from_consensus`.
        :param points: the points to fit the segment to, as an (n, 2) array.
        :param weights: if provided, a length-n array of weights to multiply each point's residual by. Weights bias
            the fitted line only; the endpoints are still the extreme projections of every point. If None, all
            points are weighted equally.
        :return: a new ``Segment2`` object representing the fitted segment.
        :raises ValueError: if there are fewer than two distinct points.
        """
        ...

    @staticmethod
    def from_consensus(points: NDArray[float], sigma_max: float, max_iterations: int | None = None,
                       refinement_steps: int | None = None, confidence: float | None = None,
                       seed: int | None = None) -> Segment2:
        """
        Fit a segment to a set of points robustly using the MAGSAC++ consensus algorithm. A robust infinite line is
        estimated (rejecting gross outliers), and the segment's endpoints are set to the extreme projections of the
        *inlier* points onto that line, so outliers influence neither the line nor the segment's extent.

        :param points: the points to fit the segment to.
        :param sigma_max: the upper bound on the expected inlier noise, in the same units as the points.
        :param max_iterations: the maximum number of minimal-sample iterations. If None, a default of 500 is used.
        :param refinement_steps: the number of iteratively reweighted refinement steps per candidate. If None, a
            default of 4 is used.
        :param confidence: the probability used for adaptive termination. If None, a default of 0.99 is used.
        :param seed: an optional fixed RNG seed for reproducible sampling. If None, a random seed is used.
        :return: a new ``Segment2`` object representing the fitted segment.
        """
        ...

    @property
    def a(self) -> Point2:
        """
        Get the first endpoint of the segment.
        :return: the first endpoint of the segment.
        """
        ...

    @property
    def b(self) -> Point2:
        """
        Get the second endpoint of the segment.
        :return: the second endpoint of the segment.
        """
        ...

    @property
    def length(self) -> float:
        """
        Get the length of the segment.
        :return: the length of the segment.
        """
        ...

    @property
    def direction(self) -> Vector2:
        """
        Get the direction vector of the segment, which is a unit vector pointing from the first endpoint to the second
        endpoint.
        :return: the direction vector of the segment.
        """
        ...

    @property
    def aabb(self) -> Aabb2:
        """
        Get the axis-aligned bounding box of the segment.
        :return: the axis-aligned bounding box of the segment.
        """
        ...

    def to_line(self) -> Line2:
        """
        Convert the segment to an infinite line passing through its two endpoints.
        :return: a Line2 passing through the segment's endpoints.
        """
        ...

    def transformed_by(self, iso: Iso2) -> Segment2:
        """
        Return a new segment with both endpoints transformed by the given isometry.
        :param iso: the isometry to apply.
        :return: a new transformed Segment2.
        """
        ...

    def at(self, t: float) -> Point2:
        """
        Return the point at parameter t, where t=0 is the first endpoint and t=1 is the second.
        :param t: the parameter value.
        :return: the point on the segment at t.
        """
        ...

    def reversed(self) -> Segment2:
        """
        Return a new segment with the endpoints reversed.
        :return: a new Segment2 with a and b swapped.
        """
        ...

    def scalar_projection(self, other: Point2) -> float:
        """
        Calculate the scalar projection of a point onto the segment, where 0.0 represents the
        first endpoint and 1.0 represents the second endpoint. The result can be any finite value,
        including negative ones or ones greater than one.
        :param other: the point to project onto the segment.
        :return: the scalar projection parameter.
        """
        ...

    def closest_point(self, other: Point2) -> Point2:
        """
        Return the closest point on the segment to the given point, clamped to the segment's
        endpoints.
        :param other: the point to find the closest location to.
        :return: the closest point on the segment.
        """
        ...

    def offset_by(self, d: float) -> Segment2:
        """
        Create a new segment shifted by distance d in the direction of the segment's normal
        vector. The normal vector is the direction vector rotated 90 degrees clockwise.
        :param d: the distance to shift the segment along its normal vector.
        :return: a new shifted Segment2.
        """
        ...

    def normal(self) -> Vector2:
        """
        Return the unit normal of the segment: the direction vector rotated 90 degrees clockwise.
        :return: the unit normal vector.
        """
        ...

    def at_t(self, t: float) -> Manifold1Pos2:
        """
        Return the manifold position (point, tangent direction, normal, and arc length) at
        parameter t, where t=0 is the first endpoint and t=1 is the second.
        :param t: the parameter value.
        :return: the manifold position at t.
        """
        ...

    def intersects_other(self, other: Segment2) -> bool:
        """
        Determine whether this segment intersects another segment, considering only the bounded
        extent of both segments (not their infinite line extensions).
        :param other: the other segment to test against.
        :return: True if the segments intersect, False otherwise.
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
    def r(self) -> float:
        """
        Get the radius of the arc
        :return: the radius of the arc
        """
        ...

    @property
    def angle0(self) -> float:
        """
        Get the start angle of the arc, in radians.
        :return: the start angle of the arc in radians.
        """
        ...

    @property
    def angle(self) -> float:
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
    def a(self) -> Point2:
        """
        Get the start point of the arc.
        :return: the start point of the arc.
        """
        ...

    @property
    def b(self) -> Point2:
        """
        Get the end point of the arc.
        :return: the end point of the arc.
        """
        ...

    @property
    def circle(self) -> Circle2:
        """
        Get the circle that has the segment as its diameter.
        :return: the circle that has the segment as its diameter.
        """
        ...

    def to_points(self, tol: float) -> NDArray:
        """
        Sample the arc into a sequence of points such that the maximum deviation between any chord
        and the true arc is at most ``tol``.

        :param tol: the maximum allowable chordal deviation.
        :return: an ``(n, 2)`` array of points along the arc, including both endpoints.
        """
        ...

    @staticmethod
    def from_circle(circle: Circle2, angle0: float, angle: float) -> Arc2:
        """
        Create an arc from a circle, a start angle, and a sweep angle.

        :param circle: the circle of which the arc is a part.
        :param angle0: the angle in radians (with respect to the x-axis) at which the arc starts.
        :param angle: the angle in radians which the arc sweeps through, beginning at ``angle0``.
            A positive value indicates a counter-clockwise sweep, while a negative value indicates
            a clockwise sweep.
        :return: a new ``Arc2``.
        """
        ...

    @staticmethod
    def from_ends(start: Point2, end: Point2, center: Point2, clockwise: bool) -> Arc2:
        """
        Create an arc from a start point, an end point, a center point, and a winding direction.

        :param start: the starting point of the arc, on the perimeter of the circle.
        :param end: the ending point of the arc, on the perimeter of the circle.
        :param center: the center point of the arc.
        :param clockwise: whether the arc sweeps clockwise (True) or counter-clockwise (False)
            from ``start`` to ``end``.
        :return: a new ``Arc2``.
        :raises ValueError: if ``start`` and ``center`` are not consistent with ``end`` and
            ``center`` being equidistant.
        """
        ...

    @staticmethod
    def from_point_angle(center: Point2, radius: float, point: Point2, angle: float) -> Arc2:
        """
        Create an arc from a center point, a radius, a point on the perimeter, and an included
        angle starting at the point.

        :param center: the arc center point.
        :param radius: the arc radius.
        :param point: a point on the perimeter of the arc at which the arc starting point should
            be located.
        :param angle: the angle in radians which the arc sweeps through, beginning at the point.
            A positive value indicates a counter-clockwise sweep, while a negative value indicates
            a clockwise sweep.
        :return: a new ``Arc2``.
        """
        ...

    @staticmethod
    def from_3_points(p0: Point2, p1: Point2, p2: Point2) -> Arc2:
        """
        Create an arc from three points. The arc will begin at the first point, pass through the
        second point, and end at the third point.

        :param p0: the starting point of the arc.
        :param p1: a point on the arc, between the start and end points.
        :param p2: the ending point of the arc.
        :return: a new ``Arc2``.
        """
        ...

    def length(self) -> float:
        """
        Return the arc length: the radius times the absolute value of the sweep angle.
        """
        ...

    def point_at_angle(self, angle: float) -> Point2:
        """
        Return the point on the arc's circle at ``angle0 + angle``.

        :param angle: the angle, in radians, measured from the arc's start angle.
        :return: the point at that angle.
        """
        ...

    def point_at_fraction(self, fraction: float) -> Point2:
        """
        Return the point at the given fraction of the arc's sweep, where 0.0 is the start and 1.0
        is the end.

        :param fraction: the fraction of the sweep angle.
        :return: the point at that fraction.
        """
        ...

    def point_at_length(self, length: float) -> Point2:
        """
        Return the point at the given arc length from the start of the arc.

        :param length: the arc length from the start.
        :return: the point at that arc length.
        """
        ...

    def is_ccw(self) -> bool:
        """Return whether the arc sweeps counter-clockwise (a positive sweep angle)."""
        ...

    def angle_interval(self) -> AngleInterval:
        """
        Return the angular interval spanned by the arc, starting at ``angle0`` and extending for
        ``angle`` radians.
        """
        ...

    def is_theta_on_arc(self, theta: float) -> bool:
        """
        Return whether the given absolute angle, in radians, falls within the arc's angular span.

        :param theta: the angle to test, in radians.
        :return: True if the angle is on the arc, False otherwise.
        """
        ...

    def theta_to_fraction(self, theta: float) -> float:
        """
        Convert an absolute angle, in radians, to the fraction of the arc's sweep at which it
        occurs, where 0.0 is the start and 1.0 is the end.

        :param theta: the angle to convert, in radians.
        :return: the fraction of the sweep at that angle.
        """
        ...

    def at_fraction(self, fraction: float) -> Point2:
        """
        Return the point at the given fraction of the arc's sweep. Equivalent to
        ``point_at_fraction``.

        :param fraction: the fraction of the sweep angle.
        :return: the point at that fraction.
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
        Expand the AABB by a given distance in all directions. The resulting height and width will be increased by
        2 * d.

        :param d: the distance to expand the AABB by.
        :return: a new AABB object with the expanded bounds.
        """
        ...

    def shrink(self, d: float) -> Aabb2:
        """
        Shrink the AABB by a given distance in all directions. The resulting height and width will be decreased by
        2 * d.

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


class Manifold1Pos2:
    """
    A position along a 1-manifold (boundary) embedded in 2D space.

    Instances are returned by spatial query methods on ``Boundary2``, such as ``at_length``, ``at_start``, ``at_end``,
    and ``at_closest_to_point``.  They carry the full geometric state at that location: where the point is, which way
    the manifold is traveling, and what the outward surface normal is.
    """

    @property
    def l(self) -> float:
        """
        The arc-length coordinate of this position measured from the start of the boundary.
        :return: the length along the boundary to this position.
        """
        ...

    @property
    def point(self) -> Point2:
        """
        The 2D world-space position of this manifold location.
        :return: the point in 2D space.
        """
        ...

    @property
    def direction(self) -> Vector2:
        """
        The unit tangent vector of the manifold at this position, pointing in the direction of increasing arc length.
        :return: the tangent direction vector.
        """
        ...

    @property
    def normal(self) -> Vector2:
        """
        The unit surface normal at this position.  By the counter-clockwise winding convention this is the tangent
        direction rotated 90 degrees clockwise (pointing to the right of travel).
        :return: the surface normal vector.
        """
        ...

    @property
    def surface_point(self) -> SurfacePoint2:
        """
        The position and surface normal combined as a ``SurfacePoint2``.
        :return: a ``SurfacePoint2`` at this manifold position.
        """
        ...

    @property
    def direction_line(self) -> Line2:
        """
        A ``Line2`` through this position whose direction matches the manifold tangent.
        :return: a ``Line2`` aligned with the tangent direction at this position.
        """
        ...


class BoundaryData2:
    """
    A mutable builder for 2D boundary geometry.

    ``BoundaryData2`` accumulates line segments and arcs in order, maintaining the continuity constraint that each
    new element must begin where the previous one ended.  Once all elements have been added, call ``to_boundary()``
    to produce a queryable ``Boundary2``.

    Boundaries may be *open* (a path with distinct start and end points) or *closed* (a loop that wraps back onto
    itself).

    **Open boundary**: provide a starting point::

        data = BoundaryData2(x=0.0, y=0.0)
        data.add_seg_xy(1.0, 0.0)

    **Closed boundary**: no starting point needed::

        data = BoundaryData2(closed=True)
        data.add_seg_xy(1.0, 0.0)
        data.add_seg_xy(2.0, 1.0)
        data.add_seg_xy(1.5, 2.0)
    """

    def __init__(self, x: float | None = None, y: float | None = None, closed: bool = False):
        """
        Create a new boundary data builder.

        For an open boundary, supply ``x`` and ``y`` as the starting point. For a closed boundary, set
        ``closed=True`` and omit ``x``/``y``. Alternately, use the ``new_open`` or ``new_closed`` static methods for
        more explicit construction.

        :param x: x-coordinate of the open boundary's start point.
        :param y: y-coordinate of the open boundary's start point.
        :param closed: if ``True``, create a closed boundary (``x`` and ``y`` are ignored).
        :raises ValueError: if ``closed=False`` and ``x`` or ``y`` is not provided.
        """
        ...

    @staticmethod
    def new_open(x: float, y: float) -> BoundaryData2:
        """
        Create an open boundary starting at ``(x, y)``.

        :param x: x-coordinate of the start point.
        :param y: y-coordinate of the start point.
        :return: a new open ``BoundaryData2``.
        """
        ...

    @staticmethod
    def new_closed() -> BoundaryData2:
        """
        Create an empty closed boundary.

        :return: a new closed ``BoundaryData2``.
        """
        ...

    def add_seg_xy(self, x: float, y: float) -> int:
        """
        Append a line segment whose end point is ``(x, y)``.  The segment begins at the end point of the previous
        element (or the boundary start point if no elements have been added yet).

        :param x: x-coordinate of the segment end point.
        :param y: y-coordinate of the segment end point.
        :return: the integer element id assigned to this segment.
        """
        ...

    def add_arc_xy(self, cx: float, cy: float, ex: float, ey: float, clockwise: bool) -> int:
        """
        Append an arc defined by a center point and an end point.  The arc begins at the end point of the previous
        element and sweeps through the given angle direction to ``(ex, ey)``.

        Note that the distance from the previous end point to the center must be the same as the distance from the
        center to the end point, or an error will be thrown when the geometry is constructed.

        :param cx: x-coordinate of the arc center.
        :param cy: y-coordinate of the arc center.
        :param ex: x-coordinate of the arc end point.
        :param ey: y-coordinate of the arc end point.
        :param clockwise: if ``True``, the arc sweeps clockwise; otherwise counter-clockwise.
        :return: the integer element id assigned to this arc.
        """
        ...

    def add_corner_fillets(
            self,
            points: list[tuple[float, float]],
            radius: float,
    ) -> list[int]:
        """
        Append a sequence of line segments joined by circular arc fillets.

        For ``n`` corner points, ``n`` segments and ``n-1`` arc fillets are created.  At least two points must be
        provided.  The method raises ``ValueError`` if the fillet radius is too large for any corner, or if the
        boundary is empty when called.

        :param points: a list of ``(x, y)`` tuples defining the segment end-points / corners.
        :param radius: the fillet radius at each interior corner.
        :return: the list of element ids created (segments and arcs interleaved).
        :raises ValueError: if fewer than two points are given or any fillet is geometrically
            invalid.
        """
        ...

    def transform_by(self, iso: Iso2) -> None:
        """
        Apply an isometry to all geometry stored in this builder in place.

        :param iso: the isometry to apply.
        """
        ...

    def is_closed(self) -> bool:
        """
        Return ``True`` if this is a closed boundary.
        :return: ``True`` for closed, ``False`` for open.
        """
        ...

    def __len__(self) -> int:
        """Return the number of elements currently in the builder."""
        ...

    def to_boundary(self) -> Boundary2:
        """
        Build a queryable ``Boundary2`` from the accumulated elements.

        :return: a new ``Boundary2``.
        :raises ValueError: if the geometry is inconsistent (e.g. duplicate points).
        """
        ...


class Boundary2:
    """
    A queryable 2D boundary composed of line segments and arcs.

    ``Boundary2`` instances are produced by ``BoundaryData2.to_boundary()``. Once created the geometry is immutable.
    The boundary may be open or closed.

    Spatial queries return ``Manifold1Pos2`` objects that carry the full geometric state (position, tangent direction,
    and surface normal) at the queried location.
    """

    def is_closed(self) -> bool:
        """
        Return ``True`` if this boundary forms a closed loop.
        :return: ``True`` for closed, ``False`` for open.
        """
        ...

    def length(self) -> float:
        """
        Return the total arc length of the boundary.
        :return: the total length.
        """
        ...

    def at_length(self, length: float) -> Manifold1Pos2:
        """
        Return the manifold position at a given arc-length along the boundary.

        :param length: the arc length, which must be in ``[0, self.length()]``.
        :return: the manifold position at ``length``.
        :raises ValueError: if ``length`` is outside the valid range.
        """
        ...

    def at_start(self) -> Manifold1Pos2:
        """
        Return the manifold position at the start of the boundary (``l = 0``).
        :return: the manifold position at the start.
        """
        ...

    def at_end(self) -> Manifold1Pos2:
        """
        Return the manifold position at the end of the boundary (``l = length()``).
        :return: the manifold position at the end.
        """
        ...

    def at_closest_to_point(self, point: Point2) -> Tuple[int, Manifold1Pos2]:
        """
        Find the point on the boundary closest to a query point.

        :param point: the query point.
        :return: a tuple of ``(element_id, manifold_position)`` where ``element_id`` is the integer id of the element
            that contains the closest point.
        """
        ...

    def at_lengths(self, lengths: NDArray[float]) -> NDArray[float]:
        """
        Evaluate the boundary at multiple arc-length positions in a single call.

        For each length value, the method returns the position, tangent direction, and surface normal at that location
        on the manifold.  The result is an ``(N, 6)`` array whose columns are arranged as:

        ======  ====================================
        Col 0   x-coordinate of the point
        Col 1   y-coordinate of the point
        Col 2   x-component of the unit tangent
        Col 3   y-component of the unit tangent
        Col 4   x-component of the surface normal
        Col 5   y-component of the surface normal
        ======  ====================================

        The tangent points in the direction of increasing arc length.  The normal is the tangent rotated 90 degrees
        clockwise (consistent with the counter-clockwise winding convention, so it points to the *right* of travel /
        outward from a CCW-wound closed boundary).

        :param lengths: a 1-D ``float64`` numpy array of arc-length positions.  Every value must lie in
            ``[0, self.length()]``; any out-of-range value raises ``ValueError``.
        :return: a numpy array of shape ``(N, 6)`` where ``N = len(lengths)``.
        :raises ValueError: if any length value is outside the valid range, or if ``lengths`` is not a contiguous
            1-D array.
        """
        ...

    def to_points(self, tol: float) -> NDArray[float]:
        """
        Sample the boundary as an array of 2D points.

        The points are placed densely enough that the chord error between adjacent points and the true boundary
        geometry does not exceed ``tol``.  On straight segments this has no effect; on arcs it controls the angular
        step size.

        :param tol: the maximum chord deviation from the true geometry.
        :return: a numpy array of shape ``(N, 2)`` containing the sampled points.
        :raises ValueError: if sampling fails.
        """
        ...

    @property
    def aabb(self) -> Aabb2:
        """
        The axis-aligned bounding box of the boundary.
        :return: the bounding box.
        """
        ...


def fit_boundary_to_points(
        points: NDArray[float],
        builder: Callable[[NDArray[float]], BoundaryData2],
        initial: NDArray[float],
        ignore_ends: bool = False,
) -> NDArray[float]:
    """
    Fit a boundary to a set of 2D points using Levenberg-Marquardt optimization.

    Residuals are the unsigned distances from each point to its nearest location on the boundary.

    :param points: a numpy array of shape ``(N, 2)`` containing the sample points.
    :param builder: a callable that accepts a 1-D ``float64`` numpy array of parameters and returns a ``BoundaryData2``.
        Raise any exception to signal that the parameter vector produces invalid geometry; the optimizer treats that
        step as a failed evaluation.
    :param initial: a 1-D ``float64`` numpy array containing the initial parameter guess. The length determines the
        number of optimization parameters.
    :param ignore_ends: if ``True``, samples that project onto the two endpoints of an open boundary contribute zero
        residual (useful when boundary extent is unknown).
    :return: a 1-D ``float64`` numpy array of the optimized parameters.
    :raises ValueError: if the initial guess causes the builder to fail, or if the optimizer does not converge.
    """
    ...


def fit_boundary_to_surface_points(
        points: NDArray[float],
        builder: Callable[[NDArray[float]], BoundaryData2],
        initial: NDArray[float],
        weight_mode: VecDot,
        ignore_ends: bool = False,
) -> NDArray[float]:
    """
    Fit a boundary to a set of 2D surface points using Levenberg-Marquardt optimization.

    ``points`` is an ``(N, 4)`` array whose columns are ``[x, y, nx, ny]``: the first two columns are the sample
    positions and the last two are the outward normal vectors (normalized internally).  Residuals are the unsigned
    distances from each sample to its nearest location on the boundary, weighted by the dot product of the sample normal
    with the boundary normal. The ``weight_mode`` parameter controls how the dot product is applied:

    * ``"as_is"``     : raw dot product (can downweight or negate antiparallel samples).
    * ``"abs"``       : absolute value (de-weights orthogonal normals, ignores direction).
    * ``"clamp_pos"`` : clamped to ``[0, 1]`` (ignores samples facing away from boundary).

    I suggest that you use ``"clamp_pos"`` if you know that the boundary normal and the surface point normals are
    definitely facing the correct directions and want to take advantage of the additional filtering.  Otherwise, I
    recommend using ``"abs"`` if you aren't sure that the boundary will be facing the correct way, since it won't
    penalize the optimizer for "flipping" the boundary normal during optimization.

    :param points: a numpy array of shape ``(N, 4)`` with columns ``[x, y, nx, ny]``.
    :param builder: a callable that accepts a 1-D ``float64`` numpy array of parameters and returns a ``BoundaryData2``.
    :param initial: a 1-D ``float64`` numpy array containing the initial parameter guess.
    :param weight_mode: the ``VecDot`` mode to use for normal-based weighting.
    :param ignore_ends: if ``True``, samples projecting onto the endpoints of an open boundary contribute zero residual.
    :return: a 1-D ``float64`` numpy array of the optimized parameters.
    :raises ValueError: if ``points`` does not have exactly 4 columns, if the initial guess causes the builder to fail,
        or if the optimizer does not converge.
    """
    ...


def rot90(dir: AngleDir) -> Iso2:
    """
    Returns an isometry representing a 90-degree rotation in the given direction.
    A CCW direction produces a positive (counter-clockwise) rotation; CW produces a negative one.

    :param dir: The direction of rotation.
    :return: An isometry representing the rotation.
    """
    ...


def rot270(dir: AngleDir) -> Iso2:
    """
    Returns an isometry representing a 270-degree rotation in the given direction.
    A CCW direction produces a negative (clockwise) rotation; CW produces a positive one.

    :param dir: The direction of rotation.
    :return: An isometry representing the rotation.
    """
    ...


def signed_angle(v1: Vector2, v2: Vector2) -> float:
    """
    Returns the signed angle from `v1` to `v2`, in radians, in the range (-π, π].
    Positive means the shortest rotation from `v1` to `v2` is counter-clockwise; negative means clockwise.

    :param v1: The reference vector.
    :param v2: The vector to measure the angle to.
    :return: The signed angle in radians.
    """
    ...


def directed_angle(v1: Vector2, v2: Vector2, direction: AngleDir) -> float:
    """
    Returns the angle from `v1` to `v2` measured in the given rotational direction, in radians, in the range [0, 2π].

    :param v1: The reference vector.
    :param v2: The vector to measure the angle to.
    :param direction: The direction in which to measure the angle.
    :return: The angle in radians.
    """
    ...


def convex_hull_2d(points: NDArray[float]) -> NDArray[float]:
    """
    Compute the convex hull of a set of 2D points, returning the hull vertex coordinates.

    The returned array contains the coordinates of the hull vertices in counter-clockwise order.
    To get the indices of those vertices in the original array instead, use ``convex_hull_idx``.

    :param points: a numpy array of shape (N, 2) containing the 2D points.
    :return: a numpy array of shape (M, 2) containing the hull vertex coordinates in
        counter-clockwise order.
    """
    ...


def convex_hull_idx(points: NDArray[float]) -> NDArray[int]:
    """
    Compute the convex hull of a set of 2D points, returning the indices of the hull vertices.

    The returned indices index into the input ``points`` array and are ordered counter-clockwise.
    To get the coordinates of the hull vertices directly instead, use ``convex_hull_2d``.

    :param points: a numpy array of shape (N, 2) containing the 2D points.
    :return: a 1-D numpy array of indices into the input array which constitute the convex hull,
        in counter-clockwise order.
    """
    ...
