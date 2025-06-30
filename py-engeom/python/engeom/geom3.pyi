from __future__ import annotations

from pathlib import Path
from typing import Tuple, Iterable, List, TypeVar, Iterator, Any

import numpy
from numpy.typing import NDArray
import engeom
import metrology

Transformable3 = TypeVar("Transformable3", Vector3, Point3, Plane3, Iso3, SurfacePoint3)
PointOrVector3 = TypeVar("PointOrVector3", Vector3, Point3)


class Vector3(Iterable[float]):
    """
    A class representing a vector in 3D space. The vector is represented by its x, y, and z components. It is
    iterable and will yield the x, y, and z components in that order, allowing the Python unpacking operator `*` to be
    used to compensate for the lack of function overloading in other parts of the library.

    A vector supports a number of mathematical operations, including addition, subtraction, scalar multiplication,
    dot and cross products, and normalization. It also supports transformation by isometry.

    Vectors have different semantics than points when it comes to transformations and some mathematical operations. Be
    sure to use the type which matches the conceptual use of the object in your code.
    """

    def __init__(self, x: float, y: float, z: float):
        """
        Create a vector in 3D space by specifying the x, y, and z components.
        :param x: the x component of the vector
        :param y: the y component of the vector
        :param z: the z component of the vector
        """
        ...

    @property
    def x(self) -> float:
        """
        Get the x component of the vector as a floating point value.
        """
        ...

    @property
    def y(self) -> float:
        """
        Get the y component of the vector as a floating point value.
        """
        ...

    @property
    def z(self) -> float:
        """
        Get the z component of the vector as a floating point value.
        """
        ...

    def __iter__(self) -> Iterator[float]:
        ...

    def __rmul__(self, other: float) -> Vector3:
        """
        Multiply the vector by a scalar value. This allows the scalar to be on the left side of the multiplication
        operator.
        :param other: the scalar value to multiply the vector by.
        :return: a new vector that is the result of the multiplication.
        """
        ...

    def __add__(self, other: PointOrVector3) -> PointOrVector3:
        """
        Add a vector to another vector or a point. Adding a vector to a point will return a new point, and adding a
        vector to a vector will return a new vector.
        :param other: the other vector or point to add to this vector.
        :return: a new vector or point that is the result of the addition.
        """
        ...

    def __sub__(self, other: Vector3) -> Vector3:
        """
        Subtract another vector from this vector.
        :param other: the other vector to subtract from this vector.
        :return: a new vector that is the result of the subtraction.
        """
        ...

    def __neg__(self) -> Vector3:
        """
        Invert the vector by negating all of its components.
        :return: a new vector in which the x, y, and z components are negated.
        """
        ...

    def __mul__(self, other: float) -> Vector3:
        """
        Multiply the vector by a scalar value.
        :param other: the scalar value to multiply the vector by.
        :return: a new vector that is the result of the multiplication.
        """
        ...

    def __truediv__(self, other: float) -> Vector3:
        """
        Divide the vector by a scalar value.
        :param other: the scalar value to divide the vector by.
        :return: a new vector that is the result of the division.
        """
        ...

    def as_numpy(self) -> NDArray[float]:
        """
        Create a numpy array of shape (3, ) from the vector.
        """
        ...

    def dot(self, other: Vector3) -> float:
        """
        Calculate the dot product of this vector with another vector.
        :param other: the other vector to calculate the dot product with.
        :return: the dot product of the two vectors.
        """
        ...

    def cross(self, other: Vector3) -> Vector3:
        """
        Calculate the cross product of this vector with another vector.
        :param other: the other vector to calculate the cross product with.
        :return: the cross product of the two vectors.
        """
        ...

    def norm(self) -> float:
        """
        Calculate the Euclidian norm (aka magnitude, length) of the vector.
        :return: the length of the vector as a floating point value.
        """

    def normalized(self) -> Vector3:
        """
        Return a normalized version of the vector. The normalized vector will have the same direction as the original
        vector, but will have a length of 1.
        :return: a new vector that has unit length
        """

    def angle_to(self, other: Vector3) -> float:
        """
        Calculate the smallest angle between this vector and another vector and return it in radians.
        :param other: the other vector to calculate the angle to.
        :return: the angle between the two vectors in radians.
        """
        ...

    def with_x(self, x: float) -> Vector3:
        """
        Return a new vector with the same y and z components as this vector, but with the x component set to the
        specified value.
        :param x: the new x component of the vector.
        :return: a new vector with the specified x component.
        """
        ...

    def with_y(self, y: float) -> Vector3:
        """
        Return a new vector with the same x and z components as this vector, but with the y component set to the
        specified value.
        :param y: the new y component of the vector.
        :return: a new vector with the specified y component.
        """
        ...

    def with_z(self, z: float) -> Vector3:
        """
        Return a new vector with the same x and y components as this vector, but with the z component set to the
        specified value.
        :param z: the new z component of the vector.
        :return: a new vector with the specified z component.
        """
        ...


class Point3(Iterable[float]):
    """
    A class representing a point in 3D space. The point is represented by its x, y, and z coordinates. It is iterable
    and will yield the x, y, and z coordinates in that order, allowing the Python unpacking operator `*` to be used to
    compensate for the lack of function overloading in other parts of the library.

    A point supports a number of mathematical operations, including addition and subtraction with vectors, subtraction
    with other points, and scaling by a scalar value. It also supports transformation by isometry.

    Points have different semantics than vectors when it comes to transformations and some mathematical operations. Be
    sure to use the type which matches the conceptual use of the object in your code.
    """

    def __init__(self, x: float, y: float, z: float):
        """
        Create a point in 3D space by specifying the x, y, and z coordinates.

        :param x: the x coordinate of the point
        :param y: the y coordinate of the point
        :param z: the z coordinate of the point
        """
        ...

    @property
    def x(self) -> float:
        """
        Get the x coordinate of the point as a floating point value.
        """
        ...

    @property
    def y(self) -> float:
        """
        Get the y coordinate of the point as a floating point value.
        """
        ...

    @property
    def z(self) -> float:
        """
        Get the z coordinate of the point as a floating point value.
        """
        ...

    def __iter__(self) -> Iterator[float]:
        ...

    @property
    def coords(self) -> Vector3:
        """
        Get the coordinates of the point as a Vector3 object.
        :return: a Vector3 object
        """
        ...

    def __sub__(self, other: PointOrVector3) -> PointOrVector3:
        """
        Subtract a vector from a point to get a new point, or subtract a point from a point to get a new vector.
        :param other: the other point or vector to subtract from this point.
        :return: a new point or vector that is the result of the subtraction.
        """
        ...

    def __add__(self, other: Vector3) -> Vector3:
        """
        Add a vector to a point to get a new point.
        :param other: the vector to add to this point.
        :return: a new point that is the result of the addition.
        """
        ...

    def __neg__(self) -> Point3:
        """
        Invert the point by negating all of its components.
        :return: a new point in which the x, y, and z components are negated
        """
        ...

    def __mul__(self, other: float) -> Point3:
        """
        Multiply the coordinates of the point by a scalar value.
        :param other: the scalar value to multiply the point by.
        :return: a new point that is the result of the multiplication.
        """
        ...

    def __rmul__(self, other: float) -> Point3:
        """
        Multiply the coordinates of the point by a scalar value. This allows the scalar to be on the left side of the
        multiplication operator.
        :param other: the scalar value to multiply the point by.
        :return: a new point that is the result of the multiplication.
        """
        ...

    def __truediv__(self, other: float) -> Point3:
        """
        Divide the coordinates of the point by a scalar value.
        :param other: the scalar value to divide the point by.
        :return: a new point that is the result of the division.
        """
        ...

    def as_numpy(self) -> NDArray[float]:
        """
        Create a numpy array of shape (2, ) from the point.
        """
        ...

    @staticmethod
    def mid(a: Point3, b: Point3) -> Point3:
        """
        Return the midpoint between two points. This is the average of the x, y, and z coordinates of the two points.
        """
        ...

    def with_x(self, x: float) -> Point3:
        """
        Return a new point with the same y and z coordinates as this point, but with the x coordinate set to the
        specified value.
        :param x: the new x coordinate of the point.
        :return: a new point with the specified x coordinate.
        """
        ...

    def with_y(self, y: float) -> Point3:
        """
        Return a new point with the same x and z coordinates as this point, but with the y coordinate set to the
        specified value.
        :param y: the new y coordinate of the point.
        :return: a new point with the specified y coordinate.
        """
        ...

    def with_z(self, z: float) -> Point3:
        """
        Return a new point with the same x and y coordinates as this point, but with the z coordinate set to the
        specified value.
        :param z: the new z coordinate of the point.
        :return: a new point with the specified z coordinate.
        """
        ...


class SurfacePoint3:
    """
    This class is used to represent a surface point in 3D space.

    Surface points are a composite structure that consist of a point in space and a normal direction. Conceptually, they
    come from metrology as a means of representing a point on the surface of an object along with the normal direction
    of the surface at that point. However, they are also isomorphic with the concept of a ray or a parameterized line
    with a direction of unit length, and can be used in that way as well.
    """

    def __init__(self, x: float, y: float, z: float, nx: float, ny: float, nz: float):
        """
        Create a surface point in 3D space by specifying the x, y, and z coordinates of the point, as well as the x, y,
        and z components of the normal vector.  The normal components will be normalized before being stored, so they
        do not need to be scaled to unit length before being passed to this constructor.

        :param x: the x coordinate of the point
        :param y: the y coordinate of the point
        :param z: the z coordinate of the point
        :param nx: the x component of the normal vector (will be normalized after construction)
        :param ny: the y component of the normal vector (will be normalized after construction)
        :param nz: the z component of the normal vector (will be normalized after construction)
        """
        ...

    @property
    def point(self) -> Point3:
        """
        Get the coordinates of the point as a Point3 object.
        :return: a Point3 object
        """
        ...

    @property
    def normal(self) -> Vector3:
        """
        Get the normal of the point as a Vector3 object.
        :return: a Vector3 object
        """
        ...

    def at_distance(self, distance: float) -> Point3:
        """
        Get the point at a distance along the normal from the surface point.
        :param distance: the distance to move along the normal.
        :return: the point at the distance along the normal.
        """
        ...

    def scalar_projection(self, point: Point3) -> float:
        """
        Calculate the scalar projection of a point onto the axis defined by the surface point position and direction.
        Positive values indicate that the point is in the normal direction from the surface point, while negative values
        indicate that the point is in the opposite direction.

        :param point: the point to calculate the projection of.
        :return: the scalar projection of the point onto the normal.
        """
        ...

    def projection(self, point: Point3) -> Point3:
        """
        Calculate the projection of a point onto the axis defined by the surface point position and direction.

        :param point: the point to calculate the projection of.
        :return: the projection of the point onto the plane.
        """
        ...

    def reversed(self) -> SurfacePoint3:
        """
        Return a new surface point with the normal vector inverted, but the position unchanged.
        :return: a new surface point with the inverted normal vector.
        """
        ...

    def planar_distance(self, point: Point3) -> float:
        """
        Calculate the planar (non-normal) distance between the surface point and a point. This is complementary to the
        scalar projection. A point is projected onto the plane defined by the position and normal of the surface point,
        and the distance between the surface point position and the projected point is returned.  The value will always
        be positive.

        :param point: the point to calculate the distance to.
        :return: the planar distance between the surface point and the point.
        """
        ...

    def get_plane(self) -> Plane3:
        """
        Get the plane defined by the surface point.
        :return: the plane defined by the surface point.
        """
        ...

    def __mul__(self, other: float) -> SurfacePoint3:
        """
        Multiply the position of the surface point by a scalar value. The normal vector is not affected unless the
        scalar is negative, in which case the normal vector is inverted.
        :param other:
        :return:
        """
        ...

    def __rmul__(self, other: float) -> SurfacePoint3:
        """
        Multiply the position of the surface point by a scalar value. The normal vector is not affected unless the
        scalar is negative, in which case the normal vector is inverted.
        :param other:
        :return:
        """
        ...

    def __truediv__(self, other: float) -> SurfacePoint3:
        """
        Divide the position of the surface point by a scalar value. The normal vector is not affected unless the
        scalar is negative, in which case the normal vector is inverted.
        :param other:
        :return:
        """
        ...

    def __neg__(self) -> SurfacePoint3:
        """
        Invert both the position AND the normal vector of the surface point.
        """
        ...

    def shift(self, offset: float) -> SurfacePoint3:
        """
        Return a new surface point with the position shifted by the specified distance along the normal vector. The
        normal vector is unchanged.
        :param offset: the distance to shift the position by along the normal vector.
        :return:
        """
        ...


class Iso3:
    """
    A class representing an isometry in 3D space. An isometry is a transformation that preserves distances and angles,
    and is also sometimes known as a rigid body transformation.  It is composed of a translation and a rotation, with
    the rotation part being internally represented by a unit quaternion.

    `Iso3` objects can be used to transform 3D points, vectors, surface points, other isometries, and a few other types
    of objects. They can also be inverted and decomposed.
    """

    def __init__(self, matrix: NDArray[float]):
        """
        Attempt to create an isometry from a 4x4 matrix in the form of a numpy array. If the matrix is not a valid
        isometry (it is not orthogonal, it has scale or shear, etc.), an exception will be raised.

        Use this method if you explicitly have a known matrix to convert to an isometry, otherwise consider using a
        composition of the `from_translation` and `from_rotation` methods.

        :param matrix: a numpy array of shape (4, 4) containing the matrix representation of the isometry.
        """
        ...

    @staticmethod
    def identity() -> Iso3:
        """ Return the identity isometry. """
        ...

    @staticmethod
    def from_translation(x: float, y: float, z: float) -> Iso3:
        """
        Create an isometry representing a translation by the specified x, y, and z components.
        :param x: the x component of the translation.
        :param y: the y component of the translation.
        :param z: the z component of the translation.
        :return: an isometry containing only a translation component
        """
        ...

    @staticmethod
    def from_rotation(angle: float, ax: float, ay: float, az: float) -> Iso3:
        """
        Create an isometry representing a rotation around an axis defined by a vector direction and the origin. The
        components of the direction will be automatically normalized before the rotation applied.

        When looking down the axis of rotation (the axis is pointing towards the observer), the rotation will be
        counter-clockwise.

        :param angle: the angle to rotate by in radians.
        :param ax: the x component of the rotation axis.
        :param ay: the y component of the rotation axis.
        :param az: the z component of the rotation axis.
        :return: the isometry representing the rotation.
        """
        ...

    @staticmethod
    def from_xyzwpr(x: float, y: float, z: float, w: float, p: float, r: float) -> Iso3:
        """
        Create an isometry from the specified translation and rotation angles in yaw, pitch, and roll format, following
        the convention typically used in robotics. The angles are specified in degrees.
        :param x:
        :param y:
        :param z:
        :param w:
        :param p:
        :param r:
        :return:
        """
        ...

    def to_xyzwpr(self) -> List[float]:
        """
        Convert the isometry to a list of translation and rotation angles in yaw, pitch, and roll format, following the
        convention typically used in robotics. The angles are returned in degrees.
        :return: a list of 6 floats representing the translation and rotation angles.
        """
        ...

    def __matmul__(self, other: Transformable3) -> Transformable3:
        """
        Multiply another object by the isometry, transforming it and returning a new object of the same type.
        :param other: an object of one of the transformable types
        :return: a new object of the same type as the input object, transformed by the isometry.
        """
        ...

    def inverse(self) -> Iso3:
        """
        Get the inverse of the isometry. The inverse is the isometry that will undo the transformation of the original
        isometry, or the isometry that when applied to the original isometry will return the identity isometry.
        """
        ...

    def transform_points(self, points: NDArray[float]) -> NDArray[float]:
        """
        Transform an array of points using the isometry. The semantics of transforming points are such that the full
        matrix is applied, first rotating the point around the origin and then translating it by the translation vector.

        To transform vectors, use the `transform_vectors` method instead.

        This is an efficient way to transform a large number of points at once, rather than using the `@` operator
        individually on a large number of `Point3` objects.

        :param points: a numpy array of shape (N, 3)
        :return: a numpy array of shape (N, 3) containing the transformed points in the same order as the input.
        """
        ...

    def transform_vectors(self, vectors: NDArray[float]) -> NDArray[float]:
        """
        Transform an array of vectors using the isometry. The semantics of transforming vectors are such that only the
        rotation matrix is applied, and the translation vector is not used. The vectors retain their original
        magnitude, but their direction is rotated by the isometry.

        To transform points, use the `transform_points` method instead.

        This is an efficient way to transform a large number of vectors at once, rather than using the `@` operator
        individually on a large number of `Vector3` objects.

        :param vectors: a numpy array of shape (N, 3)
        :return: a numpy array of shape (N, 3) containing the transformed vectors in the same order as the input.
        """
        ...

    def as_numpy(self) -> NDArray[float]:
        """
        Return a copy of the 4x4 matrix representation of the isometry.
        """
        ...

    def flip_around_x(self) -> Iso3:
        """
        Return a new isometry that flips the isometry 180° around the x-axis. The origin of the isometry will be
        preserved, but the y and z axes will point in the opposite directions.
        :return: a new isometry that is the result of the flip.
        """
        ...

    def flip_around_y(self) -> Iso3:
        """
        Return a new isometry that flips the isometry 180° around the y-axis. The origin of the isometry will be
        preserved, but the x and z axes will point in the opposite directions.
        :return: a new isometry that is the result of the flip.
        """
        ...

    def flip_around_z(self) -> Iso3:
        """
        Return a new isometry that flips the isometry 180° around the z-axis. The origin of the isometry will be
        preserved, but the x and y axes will point in the opposite directions.
        :return: a new isometry that is the result of the flip.
        """
        ...

    @property
    def origin(self) -> Point3:
        """
        Get the origin of the isometry as a Point3 object.
        :return: a Point3 object representing the origin of the isometry.
        """
        ...

    def translation(self) -> Iso3:
        """
        Return the translation component of the isometry as a separate isometry.
        """
        ...

    def rotation(self) -> Iso3:
        """
        Return the rotation component of the isometry as a separate isometry.
        """
        ...

    @staticmethod
    def from_rx(angle: float) -> Iso3:
        """
        Create an isometry representing a rotation around the x-axis by the specified angle in radians.
        :param angle: the angle to rotate by in radians.
        :return: an isometry containing only a rotation component
        """
        ...

    @staticmethod
    def from_ry(angle: float) -> Iso3:
        """
        Create an isometry representing a rotation around the y-axis by the specified angle in radians.
        :param angle: the angle to rotate by in radians.
        :return: an isometry containing only a rotation component
        """
        ...

    @staticmethod
    def from_rz(angle: float) -> Iso3:
        """
        Create an isometry representing a rotation around the z-axis by the specified angle in radians.
        :param angle: the angle to rotate by in radians.
        :return: an isometry containing only a rotation component
        """
        ...

    @staticmethod
    def from_basis_xy(e0: Vector3, e1: Vector3, origin: Point3 | None = None) -> Iso3:
        """
        Create an isometry from two vectors defining the X and Y axes.

        This method creates an isometry where:
        - The X axis aligns with e0 (normalized)
        - The Y axis aligns with e1 (normalized and adjusted to be orthogonal to e0)
        - The Z axis is calculated as the cross product of X and Y

        :param e0: Vector defining the X axis direction
        :param e1: Vector defining the approximate Y axis direction (will be orthogonalized)
        :param origin: Optional point to use as the origin of the isometry (defaults to world origin)
        :return: A new Iso3 representing the coordinate system
        :raises: ValueError if the vectors are parallel or too small to create a valid basis
        """
        ...

    @staticmethod
    def from_basis_xz(e0: Vector3, e2: Vector3, origin: Point3 | None = None) -> Iso3:
        """
        Create an isometry from two vectors defining the X and Z axes.

        This method creates an isometry where:
        - The X axis aligns with e0 (normalized)
        - The Z axis aligns with e2 (normalized and adjusted to be orthogonal to e0)
        - The Y axis is calculated as the cross product of Z and X

        :param e0: Vector defining the X axis direction
        :param e2: Vector defining the approximate Z axis direction (will be orthogonalized)
        :param origin: Optional point to use as the origin of the isometry (defaults to world origin)
        :return: A new Iso3 representing the coordinate system
        :raises: ValueError if the vectors are parallel or too small to create a valid basis
        """
        ...

    @staticmethod
    def from_basis_yz(e1: Vector3, e2: Vector3, origin: Point3 | None = None) -> Iso3:
        """
        Create an isometry from two vectors defining the Y and Z axes.

        This method creates an isometry where:
        - The Y axis aligns with e1 (normalized)
        - The Z axis aligns with e2 (normalized and adjusted to be orthogonal to e1)
        - The X axis is calculated as the cross product of Y and Z

        :param e1: Vector defining the Y axis direction
        :param e2: Vector defining the approximate Z axis direction (will be orthogonalized)
        :param origin: Optional point to use as the origin of the isometry (defaults to world origin)
        :return: A new Iso3 representing the coordinate system
        :raises: ValueError if the vectors are parallel or too small to create a valid basis
        """
        ...

    @staticmethod
    def from_basis_yx(e1: Vector3, e0: Vector3, origin: Point3 | None = None) -> Iso3:
        """
        Create an isometry from two vectors defining the Y and X axes.

        This method creates an isometry where:
        - The Y axis aligns with e1 (normalized)
        - The X axis aligns with e0 (normalized and adjusted to be orthogonal to e1)
        - The Z axis is calculated as the cross product of X and Y

        :param e1: Vector defining the Y axis direction
        :param e0: Vector defining the approximate X axis direction (will be orthogonalized)
        :param origin: Optional point to use as the origin of the isometry (defaults to world origin)
        :return: A new Iso3 representing the coordinate system
        :raises: ValueError if the vectors are parallel or too small to create a valid basis
        """
        ...

    @staticmethod
    def from_basis_zx(e2: Vector3, e0: Vector3, origin: Point3 | None = None) -> Iso3:
        """
        Create an isometry from two vectors defining the Z and X axes.

        This method creates an isometry where:
        - The Z axis aligns with e2 (normalized)
        - The X axis aligns with e0 (normalized and adjusted to be orthogonal to e2)
        - The Y axis is calculated as the cross product of Z and X

        :param e2: Vector defining the Z axis direction
        :param e0: Vector defining the approximate X axis direction (will be orthogonalized)
        :param origin: Optional point to use as the origin of the isometry (defaults to world origin)
        :return: A new Iso3 representing the coordinate system
        :raises: ValueError if the vectors are parallel or too small to create a valid basis
        """
        ...

    @staticmethod
    def from_basis_zy(e2: Vector3, e1: Vector3, origin: Point3 | None = None) -> Iso3:
        """
        Create an isometry from two vectors defining the Z and Y axes.

        This method creates an isometry where:
        - The Z axis aligns with e2 (normalized)
        - The Y axis aligns with e1 (normalized and adjusted to be orthogonal to e2)
        - The X axis is calculated as the cross product of Y and Z

        :param e2: Vector defining the Z axis direction
        :param e1: Vector defining the approximate Y axis direction (will be orthogonalized)
        :param origin: Optional point to use as the origin of the isometry (defaults to world origin)
        :return: A new Iso3 representing the coordinate system
        :raises: ValueError if the vectors are parallel or too small to create a valid basis
        """
        ...


class SvdBasis3:
    """
    A class which creates a set of orthonormal basis vectors from a set of points in 3D space. The basis is created
    using a singular value decomposition of the points, and is very similar to the statistical concept of principal
    component analysis.

    The basis can be used to determine the rank of the point set, the variance of the points along the basis vectors,
    and to extract an isometry that will transform points from the world space to the basis space.  It is useful for
    orienting unknown point sets in a consistent way, for finding best-fit lines or planes, and for other similar
    tasks.
    """

    def __init__(self, points: NDArray[float], weights: NDArray[float] | None = None):
        """
        Create a basis from a set of points. The basis will be calculated using a singular value decomposition of the
        points.

        :param points: a numpy array of shape (n, 3) containing the points to calculate the basis from.
        :param weights: a numpy array of shape (n, ) containing the weights of the points. If None, all points will be
        weighted equally.
        """
        ...

    def to_iso3(self) -> Iso3:
        """
        Produce an isometry which will transform from the world space to the basis space.

        For example, if the basis is created from a set of points that lie in an arbitrary plane, transforming the
        original points by this isometry will move the points such that all points lie on the XY plane.
        :return: the isometry that transforms from the world space to the basis space.
        """
        ...

    def largest(self) -> Vector3:
        """
        Return the largest normalized basis vector.
        :return: a Vector3 object containing the largest basis vector.
        """
        ...

    def smallest(self) -> Vector3:
        """
        Return the smallest normalized basis vector.
        :return: a Vector3 object containing the smallest basis vector.
        """
        ...

    def basis_variances(self) -> NDArray[float]:
        """
        Return the variances of the basis vectors.
        :return: a numpy array of shape (3, ) containing the variances of the basis vectors.
        """
        ...

    def basis_stdevs(self) -> NDArray[float]:
        """
        Return the standard deviations of the basis vectors.
        :return: a numpy array of shape (3, ) containing the standard deviations of the basis vectors.
        """
        ...

    def rank(self, tol: float) -> int:
        """
        Retrieve the rank of the decomposition by counting the number of singular values that are
        greater than the provided tolerance.  A rank of 0 indicates that all singular values are
        less than the tolerance, and thus the point set is essentially a single point. A rank of 1
        indicates that the point set is essentially a line. A rank of 2 indicates that the point
        set exists roughly in a plane.  The maximum rank is 3, which indicates that the point set
        cannot be reduced to a lower dimension.

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


class Plane3:
    """
    A class representing a plane in 3D space. The plane is represented by a unit normal vector and a distance from the
    origin along the normal vector.
    """

    def __init__(self, a: float, b: float, c: float, d: float):
        """
        Create a plane from the equation ax + by + cz + d = 0.
        :param a: the x value of the unit normal vector.
        :param b: the y value of the unit normal vector.
        :param c: the z value of the unit normal vector.
        :param d: the distance from the origin along the normal vector.
        """
        ...

    def inverted_normal(self) -> Plane3:
        """
        Return a new plane with the normal vector inverted.
        :return: a new plane with the inverted normal vector.
        """
        ...

    def signed_distance_to_point(self, point: Point3) -> float:
        """
        Calculate the signed distance from the plane to a point. The distance will be positive if the point is on the
        same side of the plane as the normal vector, and negative if the point is on the opposite side.
        :param point: the point to calculate the distance to.
        :return: the signed distance from the plane to the point.
        """
        ...

    def project_point(self, point: Point3) -> Point3:
        """
        Project a point onto the plane. The projected point will be the closest point on the plane to the input point.
        :param point: the point to project.
        :return: the projected point.
        """
        ...

    def intersection_distance(self, sp: SurfacePoint3) -> float | None:
        """

        :param sp:
        :return:
        """
        ...


class Mesh:
    """
    A class holding an unstructured, 3-dimensional mesh of triangles.
    """

    def __init__(
            self,
            vertices: NDArray[float],
            faces: NDArray[numpy.uint32],
            merge_duplicates: bool = False,
            delete_degenerate: bool = False
    ):
        """
        Create an engeom mesh from vertices and triangles.  The vertices should be a numpy array of shape (n, 3), while
        the triangles should be a numpy array of shape (m, 3) containing the indices of the vertices that make up each
        triangle. The triangles should be specified in counter-clockwise order when looking at the triangle from the
        front/outside.

        !!! tip
            If you get an error `TypeError: argument 'faces': 'ndarray' object cannot be converted to 'PyArray<T, D>'`,
            make sure to convert the faces array to an unsigned integer type, e.g. `numpy.uint32`.

        :param vertices: a numpy array of shape (n, 3) containing the vertices of the mesh.
        :param faces: a numpy array of shape (m, 3) containing the triangles of the mesh, should be uint.
        :param merge_duplicates: merge duplicate vertices and triangles
        :param delete_degenerate: delete degenerate triangles
        """
        ...

    @property
    def aabb(self) -> Aabb3:
        """ Get the axis-aligned bounding box of the mesh. """
        ...

    @staticmethod
    def load_stl(path: str | Path, merge_duplicates: bool = False, delete_degenerate: bool = False) -> Mesh:
        """
        Load a mesh from an STL file. This will return a new mesh object containing the vertices and triangles from the
        file.  Optional parameters can be used to control the behavior of the loader when handling duplicate vertices/
        triangles and degenerate triangles.

        !!! note
            The STL loader will automatically merge duplicate vertices due to the extreme redundancy of the STL format,
            but the `merge_duplicates` flag will control whether to merge vertices when new meshes are appended to the
            existing mesh.

        :param path: the path to the STL file to load.
        :param merge_duplicates: merge duplicate vertices and triangles. If None, the default behavior is to do nothing
        :param delete_degenerate: delete degenerate triangles. If None, the default behavior is to do nothing
        :return: the mesh object containing the data from the file.
        """
        ...

    def write_stl(self, path: str | Path):
        """
        Write the mesh to an STL file. This will write the vertices and triangles of the mesh to the file in binary
        format.

        :param path: the path to the STL file to write.
        """
        ...

    def cloned(self) -> Mesh:
        """
        Will return a copy of the mesh. This is a copy of the data, so modifying the returned mesh will not modify the
        original mesh.

        :return: an independent copy of the mesh.
        """

    def transform_by(self, iso: Iso3):
        """
        Transforms the vertices of the mesh by an isometry. This will modify the mesh in place.  Any copies made of
        the vertices will no longer match the mesh after this operation.
        :param iso: the isometry to transform the mesh by.
        """
        ...

    def append(self, other: Mesh):
        """
        Append another mesh to this mesh. This will add the vertices and triangles from the other mesh to this mesh,
        changing this one and leaving the other one unmodified.

        The `merge_duplicates` and `delete_degenerate` flags will control whether to merge duplicate vertices/triangles
        and delete degenerate triangles when appending the mesh.

        :param other: the mesh to append to this mesh, will not be modified in this operation
        """
        ...

    @property
    def vertices(self) -> NDArray[float]:
        """
        Will return an immutable view of the vertices of the mesh as a numpy array of shape (n, 3).
        :return: a numpy array of shape (n, 3) containing the vertices of the mesh.
        """
        ...

    @property
    def face_normals(self) -> NDArray[float]:
        """
        Will return an immutable view of the face normals of the mesh as a numpy array of shape (m, 3), where m is the
        number of triangles in the mesh.
        :return: a numpy array of shape (m, 3) containing the normals of the triangles of the mesh.
        """
        ...

    @property
    def vertex_normals(self) -> NDArray[float]:
        """
        Will return an immutable view of the vertex normals of the mesh as a numpy array of shape (n, 3), where n is the
        number of vertices in the mesh.  If a vertex has no faces, the normal will be (0, 0, 0), otherwise the normal
        will have been averaged from the normals of the faces that share the vertex.
        :return: a numpy array of shape (n, 3) containing the normals of the vertices of the mesh.
        """
        ...

    @property
    def faces(self) -> NDArray[numpy.uint32]:
        """
        Will return an immutable view of the triangles of the mesh as a numpy array of shape (m, 3).
        :return: a numpy array of shape (m, 3) containing the triangles of the mesh.
        """
        ...

    def split(self, plane: Plane3) -> Tuple[Mesh | None, Mesh | None]:
        """
        Split the mesh by a plane. The plane will divide the mesh into two possible parts and return them as two new
        objects.  If the part lies entirely on one side of the plane, the other part will be `None`.

        :param plane: the plane to split the mesh by.

        :return: a tuple of two optional meshes, the first being that on the negative side of the plane, the second
        being that on the positive side of the plane.
        """
        ...

    def deviation(self, points: NDArray[float], mode: engeom.DeviationMode) -> NDArray[float]:
        """
        Calculate the deviation between a set of points and their respective closest points on the mesh surface.
        There are two possible modes of computing the distance, specified using the `DeviationMode` enum.  The two
        modes are essentially the same except for how they treat points which are beyond the edge of the closest face.

        - `DeviationMode.Point`: The deviation is calculated as the direct distance from the test point to the closest
        point on the face.

        - `DeviationMode.Plane`: The deviation is calculated as the distance from the test point to the plane of the
        face on which the closest point lies. This allows for points that are slightly beyond the edge of the closest
        face to have a deviation which would be the same as if the edge of the face extended to beyond the test point.

        In both cases, the deviation will be positive if the point is outside the surface and negative if the point is
        inside the surface.

        This is a means of efficiently calculating the deviation of a large number of points from a mesh surface. For
        a single point, the `measure_point_deviation` method will provide a `Distance3` result, which is more suitable
        for visualization and reporting.

        :param points: a numpy array of shape (n, 3) containing the points to calculate the deviation for.
        :param mode: the mode to calculate the deviation in.
        :return: a numpy array of shape (n, ) containing the deviation for each point.
        """
        ...

    def sample_poisson(self, radius: float) -> NDArray[float]:
        """
        Sample the surface of the mesh using a Poisson disk sampling algorithm. This will return a numpy array of points
        and their normals that are approximately evenly distributed across the surface of the mesh. The radius parameter
        controls the minimum distance between points.

        Internally, this algorithm will first re-sample each triangle of the mesh with a dense array of points at a
        maximum distance of radius/2, before applying a random poisson disk sampling algorithm to thin the resampled
        points. This means that the output points are not based on the mesh vertices, so large triangles will not be
        under-represented and small triangles will not be over-represented.

        :param radius: the minimum distance between points.
        :return: a numpy array of shape (n, 6) containing the sampled points.
        """
        ...

    def section(self, plane: Plane3, tol: float | None = None) -> List[Curve3]:
        """
        Calculate and return the intersection curves between the mesh and a plane.

        :param plane: The plane to intersect the mesh with.
        :param tol: The curve tolerance to use when constructing the intersection curves. See the `Curve3` class
        initializer for more information on the tolerance parameter.
        :return: a list of `Curve3` objects representing the intersection curves.
        """
        ...

    def face_select_none(self) -> FaceFilterHandle:
        """
        Start a filter operation on the faces of the mesh beginning with no faces selected. This will return a filter
        object that can be used to further add or remove faces from the selection.

        :return: a filter object for the triangles of the mesh.
        """
        ...

    def face_select_all(self) -> FaceFilterHandle:
        """
        Start a filter operation on the faces of the mesh beginning with all faces selected. This will return a filter
        object that can be used to further add or remove faces from the selection.

        :return: a filter object for the triangles of the mesh.
        """
        ...

    def separate_patches(self) -> List[Mesh]:
        """
        Separate the mesh into connected patches. This will return a list of new mesh objects, each containing one
        connected patch of the original mesh. These objects will be clones of the original mesh, so modifying them will
        have no effect on the original mesh.
        :return: a list of new mesh objects containing the connected patches.
        """

    def create_from_indices(self, indices: List[int]) -> Mesh:
        """
        Create a new mesh from a list of triangle indices. This will build a new mesh object containing only the
        triangles (and their respective vertices) identified by the given list of indices.  Do not allow duplicate
        indices in the list.
        :param indices: the triangle indices to include in the new mesh
        :return: a new mesh object containing only the specified triangles
        """
        ...

    def measure_point_deviation(self, x: float, y: float, z: float,
                                dist_mode: engeom.DeviationMode) -> metrology.Distance3:
        """
        Compute the deviation of a point from this mesh's surface and return it as a measurement object.

        The deviation is the distance from the point to its closest projection onto the mesh using the specified
        distance mode.  The direction of the measurement is the direction between the point and the projection,
        flipped so that it points in the same direction as the mesh surface normal.

        If the distance is less than a very small floating point epsilon, the direction will be taken directly from the
        mesh surface normal.

        The first point `.a` of the measurement is the point on the mesh, and the second point `.b` is the test point
        that was given as an argument.

        There are two possible modes of computing the distance, specified using the `DeviationMode` enum.  The two
        modes are essentially the same except for how they treat points which are beyond the edge of the closest face.

        - `DeviationMode.Point`: The deviation is calculated as the direct distance from the test point to the closest
        point on the face.

        - `DeviationMode.Plane`: The deviation is calculated as the distance from the test point to the plane of the
        face on which the closest point lies. This allows for points that are slightly beyond the edge of the closest
        face to have a deviation which would be the same as if the edge of the face extended to beyond the test point.

        In both cases, the deviation will be positive if the point is outside the surface and negative if the point is
        inside the surface.

        This method is appropriate for measuring the deviation at a few points of interest, returning a rich object
        that contains features to aid in visualization or analysis.  For bulk measurement of large numbers of points,
        use the `deviation` method instead.

        :param x: the x component of the point to measure
        :param y: the y component of the point to measure
        :param z: the z component of the point to measure
        :param dist_mode: the deviation mode to use
        :return: a `Distance3` object containing the deviation measurement
        """

    def boundary_first_flatten(self) -> NDArray[float]:
        """
        This method will perform a conformal mapping of the mesh to the XY plane using the boundary-first flattening
        algorithm developed by Crane et al.  This mapping attempts to preserve angles from the original mesh to the
        flattened mesh, and is useful for applications such as texture mapping or transformation to an image/raster
        space for analysis.

        There are a number of limitations to this method based on the implementation:

        * There can be no non-manifold edges in the mesh.  Non-manifold edges are edges that have more than two faces
          connected to them. If there are non-manifold edges, the method will raise an exception.

        * There must be a single patch (sets of faces connected by common edges) in the mesh. If there are multiple
          patches, the method will raise an exception.

        * There can be only one boundary loop in the mesh, meaning that there can be no holes. If there are holes, the
          method will raise an exception.

        The method will return a numpy array of shape (n, 2) containing the flattened vertices of the mesh in the XY
        plane. There is no specific orientation or position guarantee to the output vertices, so they may need to be
        transformed, scaled, and/or rotated to fit a specific application.

        The 2D vertices in the output will be in the exact same order as those in the mesh and will have a 1:1
        correspondence by index, meaning that the `faces` array from the mesh also describes the triangles in the
        flattened output.

        :return: a numpy array of shape (n, 2) containing the flattened vertices of the mesh.
        """

    def surface_closest_to(self, x: float, y: float, z: float) -> SurfacePoint3:
        """
        Find the closest point on the surface of the mesh to a given point in space, returning the point and normal
        in the form of a `SurfacePoint3` object.
        :param x: the x coordinate of the point to find the closest point to
        :param y: the y coordinate of the point to find the closest point to
        :param z: the z coordinate of the point to find the closest point to
        :return: a `SurfacePoint3` object containing the closest point and normal
        """
        ...

    def visual_outline(
            self,
            facing: Vector3,
            max_edge_length: float,
            corner_angle: float | None = None
    ) -> Tuple[NDArray[float], NDArray[numpy.uint8]]:
        """
        Capture the edges of a visual outline of the mesh, used to draw a line diagram of the mesh in a 2D space. The
        returned results will consist of two numpy arrays.  The first will be a floating point array of shape (N, 6)
        where N is the number of edges in the outline. The first three columns will be the start point of the edge, and
        the last three columns will be the end point of the edge. The second array will be a numpy array of shape (N,)
        containing a 0 or 1 for each edge at the associated index. A 0 indicates that the edge is unobstructed in the
        view direction, while a 1 indicates that the edge is obstructed by the mesh.

        :param facing: A vector with the direction to look at the mesh
        :param max_edge_length: The maximum length of an edge to be included in the outline. Edges longer than this
        will be broken up into smaller edges.
        :param corner_angle: The minimum angle between two adjacent faces for the common edge to be considered a corner
        and included in the outline. If None, the default value is 45 degrees.
        :return: a tuple of two numpy arrays. The first array is the outline edges, and the second array is a mask
        indicating whether the edge is obstructed or not.
        """
        ...

    def convex_hull(self) -> Mesh:
        """
        Calculate the convex hull of the mesh. The convex hull is the smallest convex shape that contains all the
        vertices of the mesh. This will return a new mesh object containing the vertices and triangles of the convex
        hull. This method will not modify the original mesh.
        :return: a new mesh object containing the convex hull of the original mesh.
        """
        ...

    @staticmethod
    def create_box(length: float, width: float, height: float) -> Mesh:
        """
        Creates a box with the center at the origin and the specified length, width, and height

        :param length: the size of the box along the X-axis
        :param width: the size of the box along the Y-axis
        :param height: the size of the box along the Z-axis
        :return: a new `Mesh` object representing the box
        """
        ...

    @staticmethod
    def create_cylinder(radius: float, height: float, steps: int) -> Mesh:
        """
        Creates a cylinder with a radius and height. The cylinder will be centered at the origin and oriented along the
        Y-axis.

        :param radius: the radius of the cylinder
        :param height: the size of the cylinder along the Y-axis
        :param steps: the number of subdivisions to create vertices around the cylinder. The more steps the smoother the
        cylinder will be.
        :return: a new `Mesh` object representing the cylinder
        """
        ...

    @staticmethod
    def create_sphere(radius: float, n_theta: int, n_phi: int) -> Mesh:
        """
        Creates a sphere with a radius. The sphere will be centered at the origin. The step counts `n_theta` and `n_phi`
        will determine the smoothness of the sphere in the radial (n_theta) and polar (n_phi) directions. The poles
        will be located at Y=+radius and Y=-radius, and the equator will lie in the XZ plane.

        :param radius: the radius of the sphere
        :param n_theta: the number of subdivisions to create vertices around the sphere in the theta direction
        :param n_phi: the number of subdivisions to create vertices around the sphere in the phi direction
        :return: a new `Mesh` object representing the sphere
        """
        ...

    @staticmethod
    def create_cone(radius: float, height: float, steps: int) -> Mesh:
        """
        Creates a cone with a radius and height. The cone will be centered at the origin and oriented so that the
        point of the cone is located at Y=height/2 and the base is located at Y=-height/2.

        :param radius: the radius of the base of the cone
        :param height: the size of the cone along the Y-axis
        :param steps: the number of subdivisions to create vertices around the cone. The more steps the smoother the
        cone will be.
        :return: a new `Mesh` object representing the cone
        """
        ...

    @staticmethod
    def create_capsule(p0: Point3, p1: Point3, radius: float, n_theta: int, n_phi: int) -> Mesh:
        """
        Creates a capsule shape between two points with a specified radius. The capsule will be centered between the two
        points and oriented along the line connecting them. The step counts `n_theta` and `n_phi` will determine the
        smoothness of the sphere in the radial (n_theta) and polar (n_phi) directions.

        :param p0: the first point of the capsule
        :param p1: the second point of the capsule
        :param radius: the radius of the capsule
        :param n_theta: the number of subdivisions to create vertices around the sphere in the theta direction
        :param n_phi: the number of subdivisions to create vertices around the sphere in the phi direction
        :return: a new `Mesh` object representing the capsule
        """
        ...

    @staticmethod
    def create_cylinder_between(p0: Point3, p1: Point3, radius: float, steps: int) -> Mesh:
        """
        Creates a cylinder between two points with a specified radius. The cylinder will be centered between the two
        points and oriented along the line connecting them.

        :param p0: the first point of the cylinder
        :param p1: the second point of the cylinder
        :param radius: the radius of the cylinder
        :param steps: the number of subdivisions to create vertices around the cylinder. The more steps the smoother the
        cylinder will be.
        :return: a new `Mesh` object representing the cylinder
        """
        ...

    @staticmethod
    def create_rect_beam_between(p0: Point3, p1: Point3, width: float, height: float,
                                 up: Vector3 | None = None) -> Mesh:
        """
        Create a rectangular cross-sectioned prism between two points with a specified width and height. The prism will
        be centered between the two points and oriented along the line connecting them. The up vector's projection onto
        the line connecting the two end points will determine the direction of the height of the prism. If None, the
        height will be aligned with the projection of the Z-axis.

        If the up vector is parallel to the line connecting the two points, an error will be thrown.

        :param p0: the first point of the prism
        :param p1: the second point of the prism
        :param width: the width of the prism
        :param height: the height of the prism
        :param up: the up vector to use for the height direction. If None, the Z-axis will be used
        :return: a new `Mesh` object representing the prism
        """
        ...


class FaceFilterHandle:
    """
    A class that acts as a handle to a filtering (selection/deselection) operation of faces on a mesh.

    A filtering operation is started using the `face_select_all` or `face_select_none` methods on a `Mesh` object, and
    then further filtering operations can be done on the handle to select or deselect faces based on various criteria.

    Once finished, the handle can be finalized into a list of the final indices of the triangles that passed the filter,
    or used to directly create a new mesh containing only the filtered triangles.
    """

    def collect(self) -> List[int]:
        """
        Finalize the handle by collecting the final indices of the triangles that passed the filter.
        :return: a list of the final indices of the triangles that passed the filter.
        """
        ...

    def create_mesh(self) -> Mesh:
        """
        Create a new mesh from the filtered triangles. This will build a new mesh object containing only the triangles
        (and their respective vertices) that are still retained in the filter.
        :return: a new mesh object containing only the filtered triangles.
        """
        ...

    def facing(self, x: float, y: float, z: float, angle: float, mode: engeom.SelectOp) -> FaceFilterHandle:
        """
        Add, remove, or keep only the faces whose normals are facing a given direction within a certain angle tolerance.

        This method will alter the filter handle object in place and return `self` to allow for the use of a fluent-like
        interface if desired.

        :param x: the x component of the direction to check against
        :param y: the y component of the direction to check against
        :param z: the z component of the direction to check against
        :param angle: the maximum angle in radians between the face normal and the filter direction
        :param mode: the operation to perform on the faces, one of `SelectOp.Add`, `SelectOp.Remove`, or `SelectOp.Keep`
        :return: the altered filter handle object
        """
        ...

    def near_mesh(
            self,
            other: Mesh,
            all_points: bool,
            distance_tol: float,
            mode: engeom.SelectOp,
            planar_tol: float | None = None,
            angle_tol: float | None = None,
    ) -> FaceFilterHandle:
        """
        Add, remove, or keep only the faces that are within a certain distance of their closest projection onto another
        mesh. The distance can require that all three vertices of the triangle are within the tolerance, or just one.

        There are two additional optional tolerances that can be applied.

        1. A planar tolerance, which checks the distance of the vertex projected onto the plane of the reference mesh
           triangle and looks at how far it is from the projection point. This is useful to filter out triangles
           that go past the edge of the reference mesh.

        2. An angle tolerance, which checks the angle between the normal of the current triangle and the normal of the
           reference triangle. This is useful to filter out triangles that are not facing the same direction as the
           reference mesh.

        :param other: the mesh to use as a reference
        :param all_points: if True, all points of the triangle must be within the tolerance, if False, only one point
        :param distance_tol: the maximum distance between the triangle and its projection onto the reference mesh
        :param mode: the operation to perform on the faces, one of `SelectOp.Add`, `SelectOp.Remove`, or `SelectOp.Keep`
        :param planar_tol: the maximum in-plane distance between the triangle and its projection onto the reference mesh
        :param angle_tol: the maximum angle between the normals of the triangle and the reference mesh
        """
        ...


class MeshCollisionSet:
    """
    A class holding a set of meshes for collision detection. This class is used to perform collision detection between
    a set of moving and stationary meshes in 3D space.
    """

    def __init__(self):
        """
        Creates an empty collision set.
        """
        ...

    def add_stationary(self, mesh: Mesh) -> int:
        """
        Add a stationary mesh to the collision set. This mesh will be used as a reference for collision detection.
        :param mesh: the mesh to add to the collision set.
        :return: the ID of the mesh in the collision set.
        """
        ...

    def add_moving(self, mesh: Mesh) -> int:
        """
        Add a moving mesh to the collision set. This mesh will be used to check for collisions against the stationary
        meshes in the set.
        :param mesh: the mesh to add to the collision set.
        :return: the ID of the mesh in the collision set.
        """
        ...

    def add_exception(self, id1: int, id2: int):
        """
        Add an exception to the collision set. This will prevent the two meshes from being checked for collisions.
        :param id1: the ID of the first mesh.
        :param id2: the ID of the second mesh.
        """
        ...

    def check_all(self, transforms: List[Tuple[int, Iso3]], stop_at_first: bool) -> List[Tuple[int, int]]:
        """
        Check all moving meshes against all stationary meshes for collisions. This will return a list of tuples
        containing the IDs of the two meshes that are colliding.
        :return: a list of tuples containing the IDs of the colliding meshes.
        """
        ...


class CurveStation3:
    """
    A class representing a station along a curve in 3D space. The station is represented by a point on the curve, a
    tangent (direction) vector, and a length along the curve.

    These are created as the result of position finding operations on `Curve3` objects.
    """

    @property
    def point(self) -> Point3:
        """
        Get the point in 3D world space where the station is located.
        :return: the point in 3D world space
        """
        ...

    @property
    def direction(self) -> Vector3:
        """
        Get the direction vector of the curve at the location of the station. This is the tangent vector of the curve,
        and is typically the direction from the previous vertex to the next vertex.
        :return: the direction vector of the curve at the station.
        """
        ...

    @property
    def direction_point(self) -> SurfacePoint3:
        """
        A `SurfacePoint3` object representing the point on the curve and the curve's tangent/direction vector.
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


class Curve3:
    """
    A class representing a polyline in 3D space. The curve is represented by a set of vertices and the lien segments
    between them (also known as a polyline).

    !!! note
        Because this curve is a simplicial 1-complex in 3D space it cannot divide space the way a `Curve2` can in 2D
        space. As a result, it lacks the concept of left/right or inside/outside, and so does not have all the same
        features that exist in the 2D curve class.
    """

    def __init__(self, vertices: NDArray[float], tol: float = 1.0e-6):
        """
        Create a curve from a set of vertices. The vertices should be a numpy array of shape (n, 3).

        :param vertices: a numpy array of shape (n, 3) containing the vertices of the curve.
        :param tol: the inherent tolerance of the curve; points closer than this distance will be considered the same.
        """
        ...

    def clone(self) -> Curve3:
        """
        Will return a copy of the curve. This is a copy of the data, so modifying the returned curve will not modify
        the original curve.

        :return: a copy of the curve.
        """
        ...

    def length(self) -> float:
        """
        Return the total length of the curve in the units of the vertices.

        :return: the length of the curve.
        """
        ...

    @property
    def points(self) -> NDArray[float]:
        """
        Will return an immutable view of the vertices of the mesh as a numpy array of shape (n, 3).
        :return: a numpy array of shape (n, 3) containing the vertices of the mesh.
        """
        ...

    def at_length(self, length: float) -> CurveStation3:
        """
        Return a station along the curve at the given length. The length is measured from the start of the curve to the
        station. If the length is greater than the length of the curve or less than 0, an error will be raised.

        :param length: the length along the curve to return the station at.
        :return: a `CurveStation3` object representing the station along the curve.
        """
        ...

    def at_fraction(self, fraction: float) -> CurveStation3:
        """
        Return a station along the curve at the given fraction of the length of the curve. If the fraction is greater
        than 1 or less than 0, an error will be raised.

        :param fraction: the fraction of the length of the curve to return the station at.
        :return: a `CurveStation3` object representing the station along the curve.
        """
        ...

    def at_closest_to_point(self, point: Point3) -> CurveStation3:
        """
        Return a station along the curve at the closest point to the given point. The station will be the point on the
        curve that is closest to the given point.

        :param point: the point to find the closest station to.
        :return: a `CurveStation3` object representing the station along the curve.
        """
        ...

    def at_front(self) -> CurveStation3:
        """
        Return a station at the front of the curve. This is equivalent to calling `at_length(0)`.

        :return: a `CurveStation3` object representing the station at the front of the curve.
        """
        ...

    def at_back(self) -> CurveStation3:
        """
        Return a station at the back of the curve. This is equivalent to calling `at_length(length)`.

        :return: a `CurveStation3` object representing the station at the back of the curve.
        """
        ...

    def resample(self, resample: engeom.ResampleEnum) -> Curve3:
        """
        Resample the curve using the given resampling method. The resampling method can be one of the following:

        - `Resample.ByCount(count: int)`: resample the curve to have the given number of points.
        - `Resample.BySpacing(distance: float)`: resample the curve to have points spaced by the given distance.
        - `Resample.ByMaxSpacing(distance: float)`: resample the curve to have points spaced by a maximum distance.

        :param resample: the resampling method to use.
        :return: a new curve object with the resampled vertices.
        """
        ...

    def simplify(self, tolerance: float) -> Curve3:
        """
        Simplify the curve using the Ramer-Douglas-Peucker algorithm. This will remove vertices from the curve that are
        within the given tolerance of the line between the previous and next vertices.

        :param tolerance: the tolerance to use when simplifying the curve.
        :return: a new curve object with the simplified vertices.
        """
        ...

    def transformed_by(self, iso: Iso3) -> Curve3:
        """
        Transform the curve by an isometry. This will return a new curve object with the transformed vertices.

        :param iso: the isometry to transform the curve by.
        :return: a new curve object with the transformed vertices.
        """
        ...


class Aabb3:
    """
    A class representing an axis-aligned bounding box in 3D space. The bounding box is defined by a minimum point and a
    maximum point, which are the lower-left and upper-right corners of the box, respectively.

    Bounding boxes are typically used for accelerating intersection and distance queries and are used internally inside
    the Rust language `engeom` library for this purpose.  However, they have other useful applications and so are
    exposed here in the Python API.

    Typically, `Aabb3` objects will be retrieved from other `engeom` objects which use them internally, such as
    `Curve3` and `Mesh` entities.  However, they can also be created and manipulated directly.
    """

    def __init__(self, x_min: float, y_min: float, z_min: float, x_max: float, y_max: float, z_max: float):
        """
        Create an axis-aligned bounding box from the minimum and maximum coordinates.
        :param x_min: the minimum x coordinate of the box.
        :param y_min: the minimum y coordinate of the box.
        :param z_min: the minimum z coordinate of the box.
        :param x_max: the maximum x coordinate of the box.
        :param y_max: the maximum y coordinate of the box.
        :param z_max: the maximum z coordinate of the box.
        """
        ...

    @property
    def min(self) -> Point3:
        """
        Get the minimum point of the AABB.
        :return: the minimum point of the AABB.
        """
        ...

    @property
    def max(self) -> Point3:
        """
        Get the maximum point of the AABB.
        :return: the maximum point of the AABB.
        """
        ...

    @property
    def center(self) -> Point3:
        """
        Get the center point of the AABB.
        :return: the center point of the AABB.
        """
        ...

    @property
    def extent(self) -> Vector3:
        """
        The extent of the box (equivalent to `self.max - self.min`).
        :return: A vector representing the extent of the box.
        """
        ...

    @staticmethod
    def at_point(x: float, y: float, z: float, w: float, h: float | None = None, l: float | None = None) -> Aabb3:
        """
        Create an AABB centered at a point with a given width and height.
        :param x: The x-coordinate of the center of the AABB.
        :param y: The y-coordinate of the center of the AABB.
        :param z: The z-coordinate of the center of the AABB.
        :param w: The width (x extent) of the AABB.
        :param h: The height (y extent) of the AABB. If not provided, it will have the same value as the width.
        :param l: The length (z extent) of the AABB. If not provided, it will have the same value as the width.
        :return: A new axis aligned bounding box object.
        """
        ...

    @staticmethod
    def from_points(points: NDArray[float]) -> Aabb3:
        """
        Create an AABB that bounds a set of points. If the point array is empty or the wrong shape, an error will be
        thrown.
        :param points: a numpy array of shape (N, 2) containing the points to bound
        :return: a new AABB object
        """
        ...

    def expand(self, d: float) -> Aabb3:
        """
        Expand the AABB by a given distance in all directions. The resulting height and
        width will be increased by 2 * d.

        :param d: the distance to expand the AABB by.
        :return: a new AABB object with the expanded bounds.
        """
        ...

    def shrink(self, d: float) -> Aabb3:
        """
        Shrink the AABB by a given distance in all directions. The resulting height and
        width will be decreased by 2 * d.

        :param d: the distance to shrink the AABB by.
        :return: a new AABB object with the shrunk bounds.
        """
        ...

    def merged(self, other: Aabb3) -> Aabb3:
        """
        Merge this AABB with another AABB. The resulting AABB will be the smallest AABB that contains both AABBs.

        :param other: the other AABB to merge with.
        :return: a new AABB object representing the merged bounds.
        """
        ...

    def indices_contained(self, points: NDArray[float]) -> NDArray[int]:
        """
        Get the indices of the points that are contained within this AABB. The points should be a numpy array of shape
        (n, 3) where n is the number of points.

        :param points: a numpy array of shape (n, 3) containing the points to check.
        :return: a numpy array of shape (m,) containing the indices of the points that are contained within this AABB,
        where m is the number of points contained within the AABB.
        """
        ...

    def contains_point(self, point: Point3) -> bool:
        """
        Check if a point is contained within this AABB.

        :param point: the point to check.
        :return: True if the point is contained within the AABB, False otherwise.
        """
        ...


class RayBundle3:
    """
    A class representing a bundle of rays in 3D space. The rays are represented by a set of origins and directions,
    which are used to define the rays.

    This class is used for ray tracing operations, such as intersection tests with meshes or other geometric objects.
    """

    def __init__(self, array: NDArray[float]):
        """
        Create a ray bundle from a numpy array of shape (n, 6) containing the origins and directions of the rays.
        :param array: a numpy array of shape (n, 6) containing the origins and directions of the rays.
        """
        ...

    def intersect_mesh(self, mesh: Mesh, iso: Iso3 | None = None, angle: float | None = None) -> NDArray[float]:
        """

        :param angle:
        :param mesh:
        :param iso:
        :return:
        """
        ...


class PointCloud:
    """

    """

    def __init__(self, points: NDArray[float], normals: NDArray[float] | None = None,
                 colors: NDArray[numpy.uint8] | None = None):
        ...

    @property
    def points(self) -> NDArray[float]:
        """
        Get the points of the point cloud as a numpy array of shape (n, 3).
        :return: a numpy array of shape (n, 3) containing the points of the point cloud.
        """
        ...

    @property
    def normals(self) -> NDArray[float] | None:
        """
        Get the normals of the point cloud as a numpy array of shape (n, 3). If no normals were provided, this will
        return None.
        :return: a numpy array of shape (n, 3) containing the normals of the point cloud, or None if no normals were
        provided.
        """
        ...

    @property
    def colors(self) -> NDArray[numpy.uint8] | None:
        """
        Get the colors of the point cloud as a numpy array of shape (n, 3). If no colors were provided, this will
        return None.
        :return: a numpy array of shape (n, 3) containing the colors of the point cloud, or None if no colors were
        provided.
        """
        ...

    def cloned(self) -> PointCloud:
        """
        Create a copy of the point cloud. This will return a new `PointCloud` object with the same points, normals, and
        colors as the original.

        :return: a new `PointCloud` object with the same points, normals, and colors as the original.
        """
        ...

    @staticmethod
    def load_lptf3(path: str | Path, take_every : int | None = None) -> PointCloud:
        """
        Load a point cloud from a LPTF3 file. The LPTF3 format is a binary format used to store measurements from a
        triangulation-based laser profile sensor.

        :param path: the path to the LPTF3 file to load.
        :param take_every: if provided, this will take every nth row from the file, and will be used to estimate the
        spacing of the points along the profile to also sample at a rough distance
        :return: a `PointCloud` object containing the points, normals, and colors from the file.
        """
        ...

    def append(self, other: PointCloud) -> PointCloud:
        """
        Append another point cloud to this one. The points, normals, and colors from the other point cloud will be
        added to this point cloud.

        Will throw an error if the other point cloud has a different combination of normals and colors than this one.

        :param other: the other point cloud to append.
        :return: a new `PointCloud` object containing the combined points, normals, and colors.
        """
        ...

    def sample_poisson_disk(self, radius: float) -> list[int]:
        """
        Sample a subset of points from the point cloud using a Poisson disk sampling algorithm. This will return a list
        of indices of the points that were preserved. The points will be selected such that no two points are closer
        than the given radius.

        :param radius: the minimum distance between sampled points.
        :return: a list of indices of the points that were selected.
        """
        ...

    def create_from_indices(self, indices: list[int]) -> PointCloud:
        """
        Create a new point cloud from a subset of the points in this point cloud, specified by the given indices.
        The normals and colors will also be subsetted to match the points.

        :param indices: a list of indices to select from the point cloud.
        :return: a new `PointCloud` object containing the selected points, normals, and colors.
        """
        ...

    def create_from_poisson_sample(self, radius: float) -> PointCloud:
        """
        Create a new point cloud from a Poisson disk sampling of the points in this point cloud. The points will be
        selected such that no two points are closer than the given radius.

        :param radius: the minimum distance between sampled points.
        :return: a new `PointCloud` object containing the sampled points, normals, and colors.
        """
        ...

    def transform_by(self, iso: Iso3) -> PointCloud:
        """
        Transform the point cloud by an isometry. This will return a new `PointCloud` object with the transformed
        points, normals, and colors.

        :param iso: the isometry to transform the point cloud by.
        :return: a new `PointCloud` object with the transformed points, normals, and colors.
        """
        ...

    def overlap_points_by_reciprocity(self, other: PointCloud, max_distance: float) -> list[int]:
        """
        Find the indices of points in this point cloud that "overlap" with points in another point
        cloud by looking for reciprocity in the closest point in each direction.

        For each point in this point cloud "p_this", we will find the closest point in the other
        point cloud "p_other".  Then we take "p_other" and find the closest point to it in this
        point cloud, "p_recip".

        In an ideally overlapping point cloud, "p_recip" should be the same as "p_this".  We will
        use a maximum distance tolerance to determine if "p_recip" is close enough to "p_this" that
        "p_this" is considered to be overlapping with the other point cloud.

        :param other: the other point cloud to check for overlap.
        :param max_distance: the maximum distance to consider a point as overlapping.
        :return: a list of indices of points in this point cloud that overlap with points in the other point cloud.
        """
        ...
    
    def overlap_mesh_by_reciprocity(self, mesh: Mesh, max_distance: float) -> list[int]:
        """
        Find the indices of points in this point cloud that "overlap" with triangles in a mesh by looking for
        reciprocity in the closest point in each direction.

        For each point in this point cloud "p_this", we will find the closest point on the surface of the mesh 
        "p_other".  Then we will take "p_other" and find the closest point in the point cloud, "p_recip".

        In an ideally overlapping point cloud, "p_this" should be the same as "p_recip".  We will use a maximum 
        distance tolerance instead to determine if "p_recip" is close enough to "p_this" that "p_this" is 
        considered to be overlapping with the mesh.

        :param mesh: the mesh to check for overlap.
        :param max_distance: the maximum distance to consider a point as overlapping.
        :return: a list of indices of points in this point cloud that overlap with triangles in the mesh.
        """
        ...
