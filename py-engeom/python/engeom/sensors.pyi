from __future__ import annotations

from .geom3 import Point3, Mesh, Iso3, Vector3, PointCloud


class LaserProfile:
    def __init__(
            self,
            emitter_z: float,
            detector_y: float,
            detector_z: float,
            volume_width: float,
            volume_z_min: float,
            volume_z_max: float,
            resolution: int,
            angle_limit: float | None = None,
    ):
        """
        Create the base geometry of a laser profile line sensor, which emits a laser line into a
        scene and detects the reflection of that line to triangulate the distance to points on a
        surface.
       
        The general coordinate system is specified in X and Z. The center of the detection volume
        is at the origin, with the laser line ranging from the -X direction to the +X direction.
        The +Z direction points directly up towards the emitter.  The +Y direction is orthogonal to
        laser line and is typically the direction which the sensor will be panned.
       
        The geometry is specified with the following assumptions:
          - The laser line is emitted from a point directly on the +Z axis, with no offset in
            the X or Y direction.
          - The detector is not offset in the X direction, and can be specified with a Y and
            Z offset from the center of the detection volume.
          - The detection volume is trapezoidal, and its flat top and bottom are specified by a
            maximum and minimum Z value.
          - The detection volume's with is specified at Z=0, and is symmetrical around X=0.
       
        # Arguments
       
        :param emitter_z: The Z coordinate of the laser emitter. This is the height from the volume
          center where the laser fans into a triangle.
        :param detector_y: The Y coordinate of the detector's optical center. This is the out-of-plane
          offset from the plane of the laser line.
        :param detector_z: The Z coordinate of the detector's optical center. This is the height from
          the volume center where the detector's optical center is located.
        :param volume_width: The width of the detection volume at Z=0. The volume is assumed to be
          symmetrical around the X axis, ranging from -volume_width/2 to +volume_width/2.
        :param volume_z_min: The minimum Z value of the detection volume. This is the bottom of the
          trapezoidal volume, the farthest distance from the emitter where the sensor will still
          return points.
        :param volume_z_max: The maximum Z value of the detection volume. This is the top of the
          trapezoidal volume, the closest distance to the emitter where the sensor will still
          return points.
        :param resolution: The number of rays to cast across the laser line. This is the number of
          points that will be returned in the point cloud.
        :param angle_limit: An optional angle limit in radians. If specified, the sensor will only
          return a point if the angle between the surface normal at the point and the detector is
          less than this limit.
        """
        ...
    
        

    def get_points(self, target: Mesh, obstruction: Mesh | None, iso: Iso3) -> PointCloud:
        """

        :param target:
        :param obstruction:
        :param iso:
        :return:
        """
        ...


class PanningLaserProfile:
    def __init__(self, laser_line: LaserProfile, y_step: float, steps: int):
        """
        :param laser_line:
        :param y_step:
        :param steps:
        """
        ...

    def get_points(self, target: Mesh, obstruction: Mesh | None, iso: Iso3) -> PointCloud:
        """
        :param target:
        :param obstruction:
        :param iso:
        :return:
        """
        ...
