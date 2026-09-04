"""
    This example answers the first question in a machine-vision project, before a scene or robot
    is available: which lens should you use?

    The initial inputs are the required width of the part in one frame and the approximate camera
    distance. The remaining values follow from the optics. This example evaluates them in decision
    order:

    1. The focal length is not a free choice once the field width and working distance are fixed.
       `Camera.from_field_width` solves for it.
    2. The f-number is a free choice that trades two limits against each other. Stopping down
       increases the depth of field, which helps a part with relief remain sharp. It also widens
       the Airy disk. Once that disk covers more than one pixel, the sensor cannot resolve finer
       detail regardless of focus. The best aperture is where these two curves cross, and the
       crossing depends on the sensor.
    3. You can then determine whether the pixels are small enough for the feature of interest.
       This calculation uses object-space pixel size rather than megapixel count.

    The plot shows defocus blur against object distance for several apertures, with the pixel and
    diffraction limits drawn as horizontal lines to show the crossing directly.
"""

import numpy
from matplotlib.pyplot import Axes, Figure, figure, show

from engeom.sensors import WAVELENGTH_550NM_MM, Camera, Sensor

# A common 5MP machine vision sensor: 2448 x 2048 pixels on a 3.45 micron pitch. Everything in
# this example is in millimeters.
SENSOR = Sensor(2448, 2048, 0.00345)

# What we have to see, and from how far away
FIELD_WIDTH = 90.0
WORKING_DISTANCE = 250.0

# The smallest feature that has to be measurable, and how many pixels across it needs to be
FEATURE_SIZE = 0.5
PIXELS_PER_FEATURE = 4.0


def main():
    # ------------------------------------------------------------------------------------------
# Calculate the focal length from the framing requirements
    # ------------------------------------------------------------------------------------------
    reference = Camera.from_field_width(SENSOR, 4.0, WORKING_DISTANCE, FIELD_WIDTH)
    width, height = reference.footprint

    print(f"sensor            {SENSOR.width_px} x {SENSOR.height_px} px, "
          f"{SENSOR.pitch * 1000:.2f} um pitch, {SENSOR.width:.2f} x {SENSOR.height:.2f} mm")
    print(f"to see {FIELD_WIDTH:.0f} mm at {WORKING_DISTANCE:.0f} mm you need a "
          f"{reference.lens.focal_length:.1f} mm lens")
    print(f"  field of view   {width:.1f} x {height:.1f} mm "
          f"({numpy.degrees(reference.horizontal_fov):.1f} deg horizontal)")
    print(f"  object pixel    {reference.pixel_size * 1000:.1f} um")

    # Is that fine enough for the feature we care about?
    needed = FEATURE_SIZE / PIXELS_PER_FEATURE
    verdict = "enough" if reference.pixel_size <= needed else "NOT enough"
    print(f"  a {FEATURE_SIZE} mm feature at {PIXELS_PER_FEATURE:.0f} px across needs "
          f"{needed * 1000:.1f} um pixels: {verdict}")

    # ------------------------------------------------------------------------------------------
    # Evaluate the aperture tradeoff
    # ------------------------------------------------------------------------------------------
    print(f"\n{'f/':>6} {'DOF at 1px':>12} {'DOF at 2px':>12} {'Airy':>8}  limited by")
    print("  " + "-" * 52)

    f_numbers = [1.4, 2.0, 2.8, 4.0, 5.6, 8.0, 11.0, 16.0, 22.0]
    for f_number in f_numbers:
        cam = Camera.from_field_width(SENSOR, f_number, WORKING_DISTANCE, FIELD_WIDTH)
        airy = cam.airy_disk_px(WAVELENGTH_550NM_MM)

        # Below one pixel the sensor is the limit; above it, diffraction is
        limit = "pixels" if airy < 1.0 else "diffraction"
        print(f"{f_number:>6.1f} {cam.total_depth_of_field_px(1.0):>9.2f} mm "
              f"{cam.total_depth_of_field_px(2.0):>9.2f} mm {airy:>6.2f} px  {limit}")

    # The crossing is where the Airy disk first covers a whole pixel
    crossing = next(
        (f for f in numpy.arange(1.0, 32.0, 0.01)
         if Camera.from_field_width(SENSOR, float(f), WORKING_DISTANCE, FIELD_WIDTH)
         .airy_disk_px(WAVELENGTH_550NM_MM) >= 1.0),
        None,
    )
    print(f"\ndiffraction overtakes the pixels at about f/{crossing:.1f}; stopping down past that "
          f"buys depth of field at the cost of resolution")

    # ------------------------------------------------------------------------------------------
    # Plot blur against distance to show the tradeoff
    # ------------------------------------------------------------------------------------------
    fig: Figure = figure(figsize=(10, 6))
    ax: Axes = fig.subplots()

    distances = numpy.linspace(WORKING_DISTANCE - 25.0, WORKING_DISTANCE + 25.0, 400)
    for f_number in [2.0, 4.0, 8.0, 16.0]:
        cam = Camera.from_field_width(SENSOR, f_number, WORKING_DISTANCE, FIELD_WIDTH)

        # coc_px_at returns None behind the lens plane, which cannot happen over this range
        blur = numpy.array([cam.coc_px_at(float(z)) for z in distances])
        airy = cam.airy_disk_px(WAVELENGTH_550NM_MM)

        line, = ax.plot(distances, blur, label=f"f/{f_number:.1f}  (Airy {airy:.2f} px)")

        # The diffraction floor for this aperture: no amount of focus gets below it
        ax.axhline(airy, color=line.get_color(), linestyle=":", linewidth=1.0, alpha=0.7)

    ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0,
               label="one pixel")
    ax.set_xlabel("object distance (mm)")
    ax.set_ylabel("defocus blur (pixels)")
    ax.set_ylim(0.0, 8.0)
    ax.set_title(f"{reference.lens.focal_length:.0f} mm lens on a "
                 f"{SENSOR.pitch * 1000:.2f} um pitch sensor, focused at "
                 f"{WORKING_DISTANCE:.0f} mm\n"
                 "dotted lines are each aperture's diffraction floor")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    show()


if __name__ == "__main__":
    main()
