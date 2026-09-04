"""
    This example simulates what a camera on a robot records from a part and shows the
    three things that decide whether a shot is worth taking.

    A rendered view casts one ray per pixel through the scene. Each pixel that strikes the part
    records three quantities:

    * its distance, which affects the lens,
    * its defocus blur at that distance, in pixels of this camera, and
    * the surface incidence angle, which helps determine whether even perfectly sharp data is
      trustworthy.

    The three panels show these quantities. The last panel includes the part outline to make the
    shape legible. The outline comes from the mesh silhouette, creases, and boundaries rather than
    from an image edge filter. It is tested for hidden lines against the geometry and then
    projected, so faint lines represent actual edges hidden from this pose.

    Note the depth range against the depth of field printed in the console. A blade is deep enough
    that no single aperture holds all of it in focus. The blur map reveals this limitation, while
    the working distance alone does not.

    The view is deliberately rendered at reduced resolution. `Camera.rescaled` keeps the field of
    view and coarsens the sampling, so a preview costs a fraction of a full frame and frames the
    part identically.
"""

import numpy
from matplotlib.pyplot import Axes, Figure, figure, show

from _common import DATA_DIR
from engeom.geom3 import Mesh3, Point3, Vector3
from engeom.plot.matplotlib import AxesHelper, TraceBuilder
from engeom.sensors import WAVELENGTH_550NM_MM, Camera, Sensor, ViewBuffer, look_at

SENSOR = Sensor(2448, 2048, 0.00345)
F_NUMBER = 4.0
WORKING_DISTANCE = 400.0

# Wide enough that the whole blade fits: 160 mm across the sensor's long axis is 134 mm across its
# short one, and the part is 127 mm tall. Framing the whole part is what makes the outline
# readable; a tighter shot would cut the blade off at the top of the frame.
FIELD_WIDTH = 160.0

# Render at an eighth of the sensor's resolution. Blur is measured in pixels of whichever camera
# produced the view, so the preview reports smaller numbers for the same physical defocus; the
# console prints both to make the difference explicit.
PREVIEW_SCALE = 0.125

# What counts as usable data: in focus to within two pixels, and not seen edge on
MAX_BLUR_PX = 2.0
MAX_INCIDENCE = numpy.radians(60.0)


def draw_image(helper: AxesHelper, values, title: str, label: str, cmap: str):
    """Show a per-pixel map as an image, with the color scale beside it."""
    height, width = values.shape

    # The extent makes the image span the same continuous pixel coordinates the camera projects
    # into, where pixel `i` covers `[i, i + 1)`. Without it, imshow centers pixel `i` on `i` and
    # anything drawn over the image sits half a pixel off.
    image = helper.ax.imshow(values, cmap=cmap, extent=(0, width, height, 0))
    helper.draw_colorbar(mappable=image, label=label, fraction=0.046, pad=0.04)
    helper.ax.set_title(title, fontsize=10)
    helper.ax.set_xticks([])
    helper.ax.set_yticks([])


def main():
    mesh = Mesh3.load_tcmesh(DATA_DIR / "engine-blade.tcmesh")
    extent = mesh.aabb.extent
    print(f"part          {extent.x:.1f} x {extent.y:.1f} x {extent.z:.1f} mm, "
          f"{len(mesh.faces)} faces")

    # Calculate the lens from the framing, as in the lens-selection example
    camera = Camera.from_field_width(SENSOR, F_NUMBER, WORKING_DISTANCE, FIELD_WIDTH)
    near, far = camera.depth_of_field_px(MAX_BLUR_PX)
    print(f"lens          {camera.lens.focal_length:.1f} mm f/{F_NUMBER:.0f} focused at "
          f"{WORKING_DISTANCE:.0f} mm")
    print(f"  object pixel  {camera.pixel_size * 1000:.1f} um")
    print(f"  in focus      {near:.1f} .. {far:.1f} mm "
          f"(a {far - near:.1f} mm slab at {MAX_BLUR_PX:.0f} px)")
    print(f"  Airy disk     {camera.airy_disk_px(WAVELENGTH_550NM_MM):.2f} px")

    # Look at the middle of the part from the working distance, along -y, with world +z up in
    # the image. `look_at` builds the camera-to-world pose in the camera convention.
    center = mesh.aabb.center
    eye = Point3(center.x, center.y - WORKING_DISTANCE, center.z)
    iso = look_at(eye, center, Vector3(0.0, 0.0, 1.0))

    preview = camera.rescaled(PREVIEW_SCALE)
    view = ViewBuffer(preview, iso, mesh)

    depth = view.depth
    if view.hit_count == 0:
        raise SystemExit("the part is not in frame from this pose")

    d_min, d_max = float(numpy.nanmin(depth)), float(numpy.nanmax(depth))
    print(f"\nview          {view.width} x {view.height} px "
          f"(1/{1 / PREVIEW_SCALE:.0f} resolution), {view.hit_fraction * 100:.1f}% of the frame")
    print(f"  depth range   {d_min:.1f} .. {d_max:.1f} mm, a {d_max - d_min:.1f} mm spread")

    # The same physical blur reads differently on the two cameras, because a preview pixel is
    # eight times wider. Both numbers are correct about their own sensor.
    worst = d_min if abs(d_min - WORKING_DISTANCE) > abs(d_max - WORKING_DISTANCE) else d_max
    print(f"  worst blur    {preview.coc_px_at(worst):.2f} preview px "
          f"= {camera.coc_px_at(worst):.2f} px at full resolution")

    # Evaluate two quantities. Relative to every mesh face, the mask measures how much coverage
    # one shot contributes. Relative to only the visible faces, it measures shot quality: how much
    # of the in-frame data is trustworthy.
    seen = view.face_mask.count_true
    usable = view.face_mask_usable(max_coc_px=MAX_BLUR_PX, max_incidence=MAX_INCIDENCE)
    print(f"  faces seen    {seen} of {view.target_face_count} "
          f"({seen / view.target_face_count * 100:.0f}% of the part)")
    print(f"  usable        {usable.count_true} of those seen "
          f"({usable.count_true / seen * 100:.0f}% of this shot is trustworthy)")

    # A tall part on a landscape sensor leaves the sides of the frame empty. Rotating the camera
    # about its optical axis, by passing a different `up` to `look_at`, would recover most of it.
    fill_limit = min(extent.x, extent.z) / max(extent.x, extent.z)
    print(f"  note          the part is {extent.z / extent.x:.1f}x taller than wide, so a "
          f"landscape frame cannot exceed about {fill_limit * 100:.0f}% fill in this orientation")

    # The outline of the part as this camera sees it, already in pixel coordinates
    segments, codes = preview.project_outline(mesh, iso)
    print(f"  outline       {len(segments)} segments, "
          f"{int((codes == 0).sum())} visible and {int((codes == 1).sum())} hidden")

    # ------------------------------------------------------------------------------------------
    # Three panels sharing one figure
    # ------------------------------------------------------------------------------------------
    fig: Figure = figure(figsize=(15, 5.5))
    axes = fig.subplots(1, 3)

    # skip_aspect because imshow sets its own square aspect, and the helper's datalim adjustment
    # would conflict with it
    helpers = [AxesHelper(ax, skip_aspect=True) for ax in axes]

    draw_image(helpers[0], depth, "depth", "mm", "viridis")
    # Converted into full resolution pixels, because that is the number the real shot is judged
    # by. A preview pixel is 1/PREVIEW_SCALE times wider, so the blur measured in them is smaller
    # by that factor; dividing reverses the scaling.
    draw_image(helpers[1], view.coc_px / PREVIEW_SCALE, "defocus blur (full-res px)", "px",
               "magma")
    draw_image(helpers[2], numpy.degrees(view.incidence), "incidence", "deg", "cividis")

    # The outline over the incidence panel. Hidden edges are drawn faintly so the shape reads as
    # a solid object rather than a wireframe.
    ax: Axes = axes[2]
    for kind, style in [(1, dict(color="white", linewidth=0.5, alpha=0.35)),
                        (0, dict(color="white", linewidth=1.0))]:
        trace = TraceBuilder()
        for x0, y0, x1, y1 in segments[codes == kind]:
            trace.add_segment((x0, y0), (x1, y1))
        if trace.xs:
            ax.plot(*trace.xy, **style)

    # imshow already set the limits from the extent; restore them after the outline, which can
    # reach outside the frame where the part is clipped by the sensor
    ax.set_xlim(0, view.width)
    ax.set_ylim(view.height, 0)

    fig.suptitle(
        f"{camera.lens.focal_length:.0f} mm f/{F_NUMBER:.0f} at {WORKING_DISTANCE:.0f} mm  |  "
        f"{view.hit_fraction * 100:.0f}% frame fill  |  "
        f"{usable.count_true / seen * 100:.0f}% of what was seen is usable",
        fontsize=11,
    )
    fig.tight_layout()
    show()


if __name__ == "__main__":
    main()
