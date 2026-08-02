"""
    This example builds a technical illustration of a 3D mesh: an outline drawing with sections,
    dimensions and callouts, of the sort you would put in a document, rather than a shaded render of
    the sort you would rotate around on screen.

    That is what `ViewPort3` is for. It draws 3D entities onto an ordinary matplotlib `Axes` in
    parallel projection, so the result is a vector figure that composes with everything else
    matplotlib can do, including saving straight to PDF or SVG. There is no perspective, no
    lighting, and no depth sorting beyond what each method does for itself, which is good enough
    for a line drawing and absolutely the wrong way to try to visualize a scene. Use the PyVista
    helper for the latter.

    The subject of this example is half of the mesh of a threaded stud, scanned on a Zeiss ATOS 5
    in an arbitrary coordinate system.  The first part of this example uses the `SvdBasis3` entity
    and some quick mean radius tests to orient it to a coordinate system. Then we take some
    cross-sections and some measurements and display them all in the viewport.
"""

import numpy
from matplotlib.pyplot import Axes, Figure, figure, show

from _common import DATA_DIR
from engeom.geom3 import Circle3, Iso3, Mesh3, Plane3, Point3, SvdBasis3, Vector3
from engeom.metrology import Distance3
from engeom.plot.matplotlib import AxesHelper, TraceBuilder

# Where to cut the round cross-section, as a distance in from the flat cut end. Far enough in to
# clear the edge of the scan, close enough to stay on the smooth unthreaded shank.
CROSS_SECTION_INSET = 4.0


def align_to_axis(mesh: Mesh3) -> Iso3:
    """
    Build the isometry that moves the bolt from wherever the scanner left it onto the origin, with
    its axis along X.

    `SvdBasis3` is essentially a principal component analysis of the vertices: the direction of
    greatest variance becomes the first basis vector, and for a part far longer than it is wide that
    direction is the axis. `to_iso3` hands back the transform into that basis, which drops the part
    on the origin with its longest direction along X for free.

    What the decomposition cannot give you is which *way* along the axis to point, because the sign
    of an SVD basis vector is arbitrary. That has to come from something you know about the part.
    Here the flat cut end is fatter than the tapered end, so comparing the mean radius of the two
    ends settles it.
    """
    to_basis = SvdBasis3(mesh.points).to_iso3()
    points = to_basis.transform_points(mesh.points)

    # Compare the mean radius about the axis within the outer fifth at each end.
    x = points[:, 0]
    radius = numpy.hypot(points[:, 1], points[:, 2])
    low, high = x.min(), x.max()
    span = high - low
    fat_end_is_positive = (radius[x > high - span * 0.2].mean()
                           > radius[x < low + span * 0.2].mean())

    if fat_end_is_positive:
        # Spin 180 degrees about Z to put the cut end on the left, where the drawing wants it.
        return Iso3.from_rotation(numpy.pi, 0.0, 0.0, 1.0) @ to_basis
    return to_basis


def build_view() -> Iso3:
    """
    Build the isometry that transforms world space into the image plane, where +X is to the right,
    +Y is up, and +Z points into the page.

    Because the part has already been aligned onto the origin, this is nothing but a pair of
    rotations. Tip the part back to show a little of the top, then spin it for a three-quarter view.
    """
    tip = Iso3.from_rotation(numpy.deg2rad(-70.0), 1.0, 0.0, 0.0)
    spin = Iso3.from_rotation(numpy.deg2rad(-25.0), 0.0, 1.0, 0.0)
    return spin @ tip


def ground_grid(mesh: Mesh3, view: Iso3, spacing: float = 10.0) -> TraceBuilder:
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
        trace.add_segment(view @ Point3(x, ys[0], floor), view @ Point3(x, ys[-1], floor))
    for y in ys:
        trace.add_segment(view @ Point3(xs[0], y, floor), view @ Point3(xs[-1], y, floor))
    return trace


def main():
    scanned = Mesh3.load_g3d(DATA_DIR / "stud-bolt.g3d")
    print(f"Loaded {len(scanned.points)} vertices, {len(scanned.faces)} faces")

    # Put the part on the origin with its axis along X. Everything after this is expressed in that
    # frame, which is the whole point of doing it first.
    mesh = scanned.transform_copy(align_to_axis(scanned))
    low, high = mesh.aabb.min, mesh.aabb.max
    print(f"Aligned extents: x {low.x:.3f} to {high.x:.3f}, "
          f"y {low.y:.3f} to {high.y:.3f}, z {low.z:.3f} to {high.z:.3f}")

    # ---------------------------------------------------------------------------------------------
    # Sections and the fitted circle
    # ---------------------------------------------------------------------------------------------

    # A round cross-section, taken just inside the flat cut end where the shank is unthreaded and the
    # scan is clean. `section_with_plane` returns a list because one plane can cut a mesh in several
    # places at once; here it is a single closed loop.
    cut_x = low.x + CROSS_SECTION_INSET
    cross_section = mesh.section_with_plane(Plane3.from_point_normal(cut_x, 0, 0, 1, 0, 0),
                                            tol=1e-3)[0]

    # Fit a circle to that loop with the MAGSAC++ consensus algorithm rather than a plain
    # least-squares fit. `sigma_max` is an upper bound on the expected noise rather than a hard
    # inlier threshold, so it is far less touchy than the threshold a classic RANSAC needs. The
    # `_planar` variant projects onto the best-fit plane first, which is the correct choice here
    # because the section is planar by definition. A fixed `seed` keeps the result reproducible,
    # which isn't necessary but is a feature I'm showing off.
    circle: Circle3 = Circle3.from_consensus_planar(cross_section.points, sigma_max=0.15, seed=7)

    # How far the real surface departs from the fitted circle, which is the roundness of the shank.
    # We're going to use numpy here to show how, outside of the Rust ecosystem, we do have to be
    # thoughtful about using vectorization to get reasonable performance. Python will never be
    # Rust, but we don't have to completely trade away our ability to do large computations for
    # Python's convenience.
    radii = numpy.linalg.norm(cross_section.points - circle.center.as_numpy(), axis=1)
    form_error = numpy.abs(radii - circle.r)
    print(f"Cross-section at x = {cut_x:.3f}: fitted diameter {2 * circle.r:.4f}, form error "
          f"max {form_error.max():.4f}, rms {numpy.sqrt((form_error ** 2).mean()):.4f}")

    # A longitudinal section straight down the axis, which cuts the thread profile. This one comes
    # back as two curves, one for each side of the part. You wouldn't inherently know that ahead
    # of time, because how many curves come back depends on the connectivity of the triangles that
    # the section plane passes through.
    axial_sections = mesh.section_with_plane(Plane3.from_point_normal(0, 0, 0, 0, 1, 0), tol=1e-3)
    print(f"Axial section: {len(axial_sections)} curves")

    # ---------------------------------------------------------------------------------------------
    # The drawing
    # ---------------------------------------------------------------------------------------------

    # Diagrams usually want no axis decorations, and `hide_axes` turns off the ticks, labels, and spines.
    fig: Figure = figure(figsize=(12, 5))
    ax: Axes = fig.subplots()
    helper = AxesHelper(ax, hide_axes=True)

    # Ask the helper for a viewport. It draws onto the same axes the helper wraps, so 2D and 3D
    # entities can be mixed freely in one figure.
    view = helper.viewport(build_view())

    # The grid goes down first so everything else sits over it. This is the one place we reach past
    # the helpers to matplotlib directly, because `TraceBuilder` holds raw coordinates rather than
    # any `engeom` entity, and `helper.ax` is always available for that.
    helper.ax.plot(*ground_grid(mesh, view.view).xy, color="0.85", linewidth=0.6, zorder=0)

    # The outline drawing itself. Silhouette edges and, with `corner_angle` given, the sharp interior
    # corners are found for this view direction and classified as visible or hidden. The two classes
    # get separate styles, which is what makes a line drawing readable: hidden detail stays faint
    # instead of competing with the profile.
    view.draw_mesh_outline(
        mesh,
        visible_kwargs={"color": "black", "linewidth": 1.2, "label": "outline"},
        hidden_kwargs={"color": "black", "linewidth": 0.5, "alpha": 0.15},
        corner_angle=numpy.deg2rad(35.0),
        max_edge_len=1.0,
    )

    # The thread profile, as projected polylines. Only the first curve carries the legend label, so
    # that the pair produces one entry rather than two identical ones.  If you didn't want any legend
    # labels you would just make one call to `draw_curve(*axial_sections, ...)`, this just a trick to
    # manage the legend clutter.
    view.draw_curve(axial_sections[0], color="magenta", linewidth=1.0, label="axial section")
    view.draw_curve(*axial_sections[1:], color="magenta", linewidth=1.0)

    # The measured cross-section, and the circle fitted to it. A circle whose plane is not facing the
    # viewer projects to an ellipse, which `draw_circle` computes exactly rather than approximating
    # with a polyline, so the fit lands right on top of the section it came from.
    #
    # The fit is given an explicit `zorder` because it is a patch and the section is a line, and
    # matplotlib draws all patches beneath all lines regardless of the order they were added in.
    # Without it the fit would be hidden underneath the very curve it is meant to be compared to.
    view.draw_curve(cross_section, color="tab:blue", linewidth=3.0, alpha=0.45,
                    label="cross-section")
    view.draw_circle(circle, color="tab:red", linewidth=1.4, linestyle="--", zorder=5,
                     label=f"consensus fit, ⌀{2 * circle.r:.3f}")

    # ---------------------------------------------------------------------------------------------
    # Dimensions and callouts
    # ---------------------------------------------------------------------------------------------

    # The overall length. Note that the label reports the true 3D distance: the arrow is
    # foreshortened by the projection, as it must be, but a dimension whose number disagreed with the
    # geometry it measures would be worse than useless.
    standoff = Vector3(0.0, 0.0, -6.0)
    view.draw_distance(
        Distance3(Point3(low.x, 0.0, low.z) + standoff, Point3(high.x, 0.0, low.z) + standoff,
                  Vector3(1.0, 0.0, 0.0)),
        template="{value:.2f} mm overall",
        label_place="inside",
        fontsize=11,
    )

    # The fitted diameter, dimensioned across the circle in its own section plane. This one is
    # heavily foreshortened in this view, which reiterates the point above: the arrow is short,
    # the number is not.
    across = Vector3(0.0, 0.0, 1.0) * circle.r
    view.draw_distance(
        Distance3(circle.center - across, circle.center + across, Vector3(0.0, 0.0, 1.0)),
        template="⌀{value:.3f}",
        label_place="outside_rev",
        fontsize=10,
    )

    # Callouts. `find_mesh_edge_point` walks out from the part center in a 2D *view* direction and
    # returns the point on the mesh it reaches, so a leader lands on the silhouette rather than
    # somewhere inside the part. The direction is in view coordinates, not world ones, which is what
    # makes it useful for placing annotations: "the top of the part as drawn" is what you mean.
    #
    # `offset_2d` shifts the label in data units after projection, so these are millimeters on the
    # page rather than typographic points.
    view.draw_labeled_point(view.find_mesh_edge_point(-1.0, 0.2, mesh), "flat cut end",
                            offset_2d=(-12.0, 6.0), arrow=True, box=True, fontsize=10)
    view.draw_labeled_point(view.find_mesh_edge_point(0.3, 1.0, mesh), "rounded end",
                            offset_2d=(4.0, 9.0), arrow=True, box=True, fontsize=10)

    # A coordinate frame for orientation, sitting on the origin the alignment put the part on. Any
    # axis pointing nearly straight into or out of the page would project to a meaningless stub, so
    # the draw function will automatically skip them. Try playing with the view direction if you
    # want to see it in action.
    view.draw_coordinate_system(Iso3.identity(), 12.0)

    ax.legend(loc="upper right", fontsize=9)
    ax.set_title("ViewPort3: sections and a consensus circle fit on a scanned threaded stud")
    fig.tight_layout()

    # Called last, once everything is drawn and the layout is settled. `set_aspect("equal")` alone
    # keeps a true aspect by shrinking the plot inside its box, leaving a band of dead space down
    # one side; this widens whichever axis has room to spare instead, so the figure stays full and
    # the geometry keeps its shape. Nothing already in view is ever cropped.
    helper.fill_available_space()

    show()


if __name__ == "__main__":
    main()
