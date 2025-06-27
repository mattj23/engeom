class Distance2:
    """
    Represents a distance between two points in 2D space.
    """

    from .geom2 import Point2, Vector2, SurfacePoint2
    from .geom3 import Iso3

    def __init__(self, a: Point2, b: Point2, direction: Vector2 | None = None):
        """
        Initialize a new 2D distance object.
        :param a: The first point in the distance calculation.
        :param b: The second point in the distance calculation.
        :param direction: The direction of the distance calculation. If not provided, the direction will be calculated
        automatically as the vector `b - a`.
        """
        ...

    @property
    def a(self) -> Point2:
        """
        Get the first point in the distance calculation.
        :return: The first point
        """
        ...

    @property
    def b(self) -> Point2:
        """
        Get the second point in the distance calculation.
        :return: The second point
        """
        ...

    @property
    def direction(self) -> Vector2:
        """
        Get the direction of the distance calculation.
        :return: The direction vector
        """
        ...

    @property
    def value(self) -> float:
        """
        Get the signed distance scalar value.
        :return: the signed distance
        """
        ...

    @property
    def center(self) -> SurfacePoint2:
        """
        Get a center surface point, located halfway between the `a` and `b` points and with a normal facing the
        `direction` vector.
        :return: the center surface point
        """
        ...

    def to_3d(self, iso: Iso3) -> Distance3:
        """
        Convert this 2D distance to a 3D distance by adding a zero Z component to the points and direction and then
        transforming them using the provided isometry.
        :param iso: The isometry to transform the entity by after adding a zero Z component.
        :return: The 3D distance object
        """
        ...


class Distance3:
    """
    Represents a distance between two points in 3D space.
    """

    from .geom3 import Point3, Vector3, SurfacePoint3, Iso3

    def __init__(self, a: Point3, b: Point3, direction: Vector3 | None = None):
        """
        Initialize a new 3D distance object.
        :param a: The first point in the distance calculation.
        :param b: The second point in the distance calculation.
        :param direction: The direction of the distance calculation. If not provided, the direction will be calculated
        automatically as the vector `b - a`.
        """
        ...

    @property
    def a(self) -> Point3:
        """
        Get the first point in the distance calculation.
        :return: The first point
        """
        ...

    @property
    def b(self) -> Point3:
        """
        Get the second point in the distance calculation.
        :return: The second point
        """
        ...

    @property
    def direction(self) -> Vector3:
        """
        Get the direction of the distance calculation.
        :return: The direction vector
        """
        ...

    @property
    def value(self) -> float:
        """
        Get the signed distance scalar value.
        :return: the signed distance
        """
        ...

    @property
    def center(self) -> SurfacePoint3:
        """
        Get a center surface point, located halfway between the `a` and `b` points and with a normal facing the
        `direction` vector.
        :return: the center surface point
        """
        ...

    def to_2d(self, iso: Iso3) -> Distance2:
        """
        Convert this 3D distance to a 2D distance by transforming the points and direction using the provided isometry
        and then removing the Z component.

        :param iso: The isometry to transform the entity by before removing the Z component.
        :return: The 2D distance object
        """
        ...

