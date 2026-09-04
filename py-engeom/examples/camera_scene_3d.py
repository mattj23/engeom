"""
    This example shows the simulated shot from `camera_scene_preview.py` from an external
    viewpoint, revealing spatial relationships that the 2D image cannot show.

    The image shows the recorded data, while the 3D view explains that data through three elements:

    * The camera's coordinate frame, which makes the convention explicit. The
      blue Z axis runs into the scene along the optical axis, and the green Y axis points down in
      the image, which is the computer vision convention this module uses.
    * The view volume, as two nested frustums. The outer one is everything the sensor can see
      between two arbitrary depths; the inner one is the slab that is actually in focus. Seeing
      its thickness against the part depth demonstrates an important distinction: a working
      distance specifies a plane, while the part must fit within the in-focus region.
    * The simulated point cloud, which is where the camera's rays actually landed. It is denser
      where the surface faces the camera and thins out where the surface turns away, because the
      rays are spaced evenly on the sensor rather than on the part. The resulting bias means that
      the cloud is not a uniform sample of the surface.

    The plotter's own camera is set to look from the simulated camera's position, which requires
    converting between the two conventions this codebase uses: the camera pose is camera-to-world
    with +Y down and +Z into the scene, while a plotting view is world-to-view with +Y up and +Z
    toward the viewer. They differ by a half turn about X and an inversion.
"""

import numpy
from pyvista import Plotter

from _common import DATA_DIR
from engeom.geom3 import Iso3, Mesh3, Point3, Vector3
from engeom.plot.pyvista.convert import LineBuilder
from engeom.sensors import Camera, Sensor, ViewBuffer, look_at

SENSOR = Sensor(2448, 2048, 0.00345)
F_NUMBER = 4.0
WORKING_DISTANCE = 400.0
FIELD_WIDTH = 160.0
PREVIEW_SCALE = 0.06

MAX_BLUR_PX = 2.0


def frustum_lines(corners) -> LineBuilder:
    """
    Build the twelve edges of a view volume from the eight corners.

    `frustum_corners` returns the near face first, then the far face, each running from the top
    left of the image clockwise, so the two faces are loops over the first and second groups of
    four and the connecting edges join matching indices.
    """
    builder = LineBuilder()
    for face in (0, 4):
        for i in range(5):
            builder.add(corners[face + i % 4])
        builder.skip()
    for i in range(4):
        builder.add(corners[i])
        builder.add(corners[i + 4])
        builder.skip()
    return builder


def apex_lines(corners, eye) -> LineBuilder:
    """The four sight lines running from the camera out through the corners of the view volume."""
    builder = LineBuilder()
    for i in range(4):
        builder.add(eye)
        builder.add(corners[i + 4])
        builder.skip()
    return builder


def main():
    mesh = Mesh3.load_tcmesh(DATA_DIR / "engine-blade.tcmesh")
    camera = Camera.from_field_width(SENSOR, F_NUMBER, WORKING_DISTANCE, FIELD_WIDTH)

    center = mesh.aabb.center
    eye = Point3(center.x, center.y - WORKING_DISTANCE, center.z)
    iso = look_at(eye, center, Vector3(0.0, 0.0, 1.0))

    near_focus, far_focus = camera.depth_of_field_px(MAX_BLUR_PX)
    print(f"lens          {camera.lens.focal_length:.1f} mm f/{F_NUMBER:.0f} at "
          f"{WORKING_DISTANCE:.0f} mm")
    print(f"  in focus      {near_focus:.1f} .. {far_focus:.1f} mm "
          f"({far_focus - near_focus:.1f} mm deep)")
    print(f"  part depth    {mesh.aabb.extent.y:.1f} mm along the optical axis")
    if mesh.aabb.extent.y > far_focus - near_focus:
        print("  the part is deeper than the focus slab, so one shot cannot hold all of it sharp")

    # Use a coarse render because this view illustrates the geometry rather than image detail
    view = ViewBuffer(camera.rescaled(PREVIEW_SCALE), iso, mesh)
    cloud = view.to_point_cloud()
    print(f"\ncloud         {cloud.point_count} points from a "
          f"{view.width} x {view.height} render")

    plotter = Plotter(shape=(1, 2), window_size=(1500, 750))

    depth = mesh.aabb.extent.y
    outer = camera.frustum_corners(iso, WORKING_DISTANCE - depth, WORKING_DISTANCE + depth)
    focus = camera.frustum_corners(iso, near_focus, far_focus)

    # ------------------------------------------------------------------------------------------
    # Left: show the setup from an external viewpoint
    # ------------------------------------------------------------------------------------------
    plotter.subplot(0, 0)
    plotter.engeom.draw_mesh(mesh, color="lightgray", opacity=0.9)

    plotter.add_mesh(frustum_lines(outer).build(), color="steelblue", line_width=2)
    plotter.add_mesh(frustum_lines(focus).build(), color="orange", line_width=5)
    plotter.add_mesh(apex_lines(outer, eye).build(), color="steelblue", line_width=1,
                     opacity=0.5)
    plotter.engeom.draw_point_cloud(cloud, color="firebrick", point_size=2.0)
    plotter.engeom.draw_coordinate_system(iso, length=WORKING_DISTANCE * 0.12, line_width=4)
    plotter.engeom.draw_point(eye, color="black", point_size=14.0)
    plotter.add_text(
        "blue: sight lines and framed volume\n"
        "orange: the slab that is actually in focus\n"
        "red: where the rays landed",
        font_size=9,
    )

    # Off to one side and above, so the frustum reads as a volume and the thinness of the focus
    # slab against the depth of the part is visible. Looking down the optical axis would collapse
    # both into a rectangle.
    plotter.engeom.view_from((-1.0, -1.0, 0.45), up=(0.0, 0.0, 1.0))

    # Framed to include the camera itself, so the standoff is visible against the part. The part
    # ends up small, accurately showing that the working distance is many times its depth.
    # The eye is added to the corner array rather than passed on its own: a single point has no
    # extent, and `fit_to` refuses an extent with no size.
    with_eye = numpy.vstack([outer, [[eye.x, eye.y, eye.z]]])
    plotter.engeom.fit_to(mesh, with_eye, pad=0.05)
    plotter.add_axes()

    # ------------------------------------------------------------------------------------------
    # Right: the same scene from the camera's own pose, which checks the convention conversion
    # ------------------------------------------------------------------------------------------
    plotter.subplot(0, 1)
    plotter.engeom.draw_mesh(mesh, color="lightgray")
    plotter.engeom.draw_point_cloud(cloud, color="firebrick", point_size=3.0)
    plotter.add_text("from the camera's own pose", font_size=9)

    # A plotting view is world-to-view with +Y up and +Z toward the viewer, while a camera pose is
    # camera-to-world with +Y down and +Z into the scene, so the two differ by a half turn about X
    # and an inversion. If this is right, the part is oriented here the way the rendered image
    # showed it; if the half turn were dropped, it would appear upside down.
    plotter.engeom.view_pose(Iso3.from_rx(numpy.pi) @ iso.inverse())
    plotter.engeom.fit_to(mesh, pad=0.02)

    plotter.show()


if __name__ == "__main__":
    main()
