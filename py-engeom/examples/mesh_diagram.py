"""
    This example builds a technical illustration of a 3D mesh: an outline drawing with dimensions
    and callouts, of the sort you would put in a document, rather than a shaded render of the sort
    you would rotate around on screen.

    That is what `ViewPort3` is for. It draws 3D entities onto an ordinary matplotlib `Axes` in
    parallel projection, so the result is a vector figure that composes with everything else
    matplotlib can do, including saving straight to PDF or SVG. There is no perspective, no
    lighting, and no depth sorting beyond what each method does for itself, which is good enough
    for a line drawing and absolutely the wrong way to try to visualize a scene. Use the PyVista
    helper for the latter.
"""

import numpy
from matplotlib.pyplot import Axes, Figure, figure, show

from _common import DATA_DIR
from engeom.geom3 import Iso3, Mesh3, Plane3, Point3, Vector3
from engeom.metrology import Distance3
from engeom.plot.matplotlib import AxesHelper, TraceBuilder


def build_view(mesh: Mesh3) -> Iso3:
    """
    Build the isometry that transforms world space into the image plane, where +X is to the right,
    +Y is up, and +Z points into the page.

    The mesh sits wherever the scanner left it, so the view first brings its center to the origin
    and then rotates. Composing in that order means the rotations pivot about the part rather than
    about some distant world origin, which is what you want when aiming a view by hand.
    """
    center = mesh.aabb.center
    recenter = Iso3.from_translation(-center.x, -center.y, -center.z)

    # Roll the part so its long axis lies across the page, then tip it back to get a three-quarter
    # view with a little of the top showing.
    tip = Iso3.from_rotation(numpy.deg2rad(-70.0), 1.0, 0.0, 0.0)
    spin = Iso3.from_rotation(numpy.deg2rad(-25.0), 0.0, 1.0, 0.0)
    return spin @ tip @ recenter


def ground_grid(mesh: Mesh3, view_port, spacing: float = 10.0) -> TraceBuilder:
    """
    Build a grid on the plane under the part, to give the eye something to read the orientation
    against.

    A grid is a pile of disjoint line segments, and drawing each as its own matplotlib artist would
    be both slow and awkward to restyle. `TraceBuilder` exists for exactly this: it accumulates
    runs of points with `None` breaks between them, so the whole grid ends up as a single `Line2D`
    that can be styled, hidden, or re-ordered as one thing.
    """
    box = mesh.aabb.expand(spacing)
    floor = mesh.aabb.min.z
    xs = numpy.arange(box.min.x, box.max.x + spacing, spacing)
    ys = numpy.arange(box.min.y, box.max.y + spacing, spacing)

    trace = TraceBuilder()
    for x in xs:
        trace.add_segment(view_port.view @ Point3(x, ys[0], floor),
                          view_port.view @ Point3(x, ys[-1], floor))
    for y in ys:
        trace.add_segment(view_port.view @ Point3(xs[0], y, floor),
                          view_port.view @ Point3(xs[-1], y, floor))
    return trace


def main():
    mesh = Mesh3.load_g3d(DATA_DIR / "stud-bolt.g3d")
    print(f"Loaded {len(mesh.points)} vertices, {len(mesh.faces)} faces")

    # Diagrams want no axis decorations, and `hide_axes` turns off the ticks, labels, and spines.
    fig: Figure = figure(figsize=(12, 5))
    ax: Axes = fig.subplots()
    helper = AxesHelper(ax, hide_axes=True)

    # Ask the helper for a viewport. It draws onto the same axes the helper wraps, so 2D and 3D
    # entities can be mixed freely in one figure.
    view = helper.viewport(build_view(mesh))

    # The grid goes down first so everything else sits over it. This is the one place we reach past
    # the helpers to matplotlib directly, because `TraceBuilder` holds raw coordinates rather than
    # any `engeom` entity, and `helper.ax` is always available for that.
    grid = ground_grid(mesh, view)
    helper.ax.plot(*grid.xy, color="0.85", linewidth=0.6, zorder=0)

    # The outline drawing itself. Silhouette edges and, with `corner_angle` given, the sharp
    # interior corners are found for this view direction and classified as visible or hidden. The
    # two classes get separate styles, which is what makes a line drawing readable: hidden detail
    # stays faint instead of competing with the profile.
    view.draw_mesh_outline(
        mesh,
        visible_kwargs={"color": "black", "linewidth": 1.2, "label": "outline"},
        hidden_kwargs={"color": "black", "linewidth": 0.5, "alpha": 0.15},
        corner_angle=numpy.deg2rad(35.0),
        max_edge_len=1.0,
    )

    # A section curve through the middle of the shank, cut with a plane and drawn as a `Curve3`.
    # `section_with_plane` returns a list because a plane can cut a mesh in several places at once.
    center = mesh.aabb.center
    section_plane = Plane3.from_point_normal(center.x, center.y, center.z, 1.0, 0.0, 0.0)
    view.draw_curve(*mesh.section_with_plane(section_plane, tol=1e-3),
                    color="tab:blue", linewidth=1.5, label="section")

    # The overall length, as a dimension. Note that the label reports the true 3D distance: the
    # arrow is foreshortened by the projection, as it must be, but a dimension whose number
    # disagreed with the geometry it measures would be worse than useless.
    low = mesh.aabb.min
    high = mesh.aabb.max
    standoff = Vector3(0.0, 0.0, -6.0)
    view.draw_distance(
        Distance3(Point3(low.x, center.y, low.z) + standoff,
                  Point3(high.x, center.y, low.z) + standoff,
                  Vector3(1.0, 0.0, 0.0)),
        template="{value:.2f} mm",
        label_place="inside",
        fontsize=11,
    )

    # Two callouts. `find_mesh_edge_point` walks out from the part center in a 2D *view* direction
    # and returns the point on the mesh it reaches, so a leader lands on the silhouette rather than
    # somewhere inside the part. The direction is in view coordinates, not world ones, which is what
    # makes it useful for placing annotations: "the top of the part as drawn" is what you mean.
    #
    # `offset_2d` shifts the label in data units after projection, so these are millimeters on the
    # page rather than typographic points.
    view.draw_labeled_point(view.find_mesh_edge_point(0.3, 1.0, mesh), "rounded end",
                            offset_2d=(4.0, 9.0), arrow=True, box=True, fontsize=10)
    view.draw_labeled_point(view.find_mesh_edge_point(-1.0, 0.2, mesh), "cut end",
                            offset_2d=(-11.0, 6.0), arrow=True, box=True, fontsize=10)

    # A coordinate frame for orientation. Any axis pointing nearly straight into or out of the page
    # would project to a meaningless stub, so those are skipped automatically.
    view.draw_coordinate_system(Iso3.from_translation(low.x - 6.0, center.y, low.z - 4.0), 12.0)

    ax.legend(loc="upper right")
    ax.set_title("ViewPort3: parallel-projection diagram of a scanned threaded stud")
    fig.tight_layout()

    # Called last, once everything is drawn and the layout is settled. `set_aspect("equal")` alone
    # keeps a true aspect by shrinking the plot inside its box, leaving a band of dead space down
    # one side; this widens whichever axis has room to spare instead, so the figure stays full and
    # the geometry keeps its shape. Nothing already in view is ever cropped.
    helper.fill_available_space()

    show()


if __name__ == "__main__":
    main()
