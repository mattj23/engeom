"""
    This example shows how coloring a mesh can answer the questions that drive an investigation:
    not merely "what does the part look like?" but "where is it out, and which faces are those?"

    Four panels, sharing one camera:

    1. Deviations of a measured point set from the nominal surface, colored through the same GOM
       color map the Matplotlib helper uses, with the out-of-range colors carried across so points
       beyond the scale are called out rather than clamped.
    2. A face selection highlighted on the mesh itself, in a single actor, so the highlight cannot
       z-fight with the surface it sits on.
    3. The feature edges, which read a part's shape far better than drawing every triangle does.
    4. The same mesh drawn twice at different poses without copying its points, which is what makes
       a before-and-after of a transform cheap.
"""

import math

import numpy
from pyvista import Plotter

from _common import DATA_DIR
from engeom.align3 import AlignParams3, points_to_mesh
from engeom.geom3 import Iso3, Mesh3
from engeom.plot.matplotlib.colors import GOM_CMAP


def main():
    mesh = Mesh3.load_tcmesh(DATA_DIR / "engine-blade.tcmesh")

    # Stand in for a measurement. Sampling the surface and moving the points rigidly would not do:
    # the alignment inverts a rigid motion exactly, and every deviation would come back zero. So the
    # points are bowed, growing with height like a blade that has crept in service, and then put
    # back on the nominal surface by an alignment. What survives that is real form error, which is
    # what an inspection is actually looking at.
    sampled = mesh.sample_poisson(2.0).points
    height = (sampled[:, 2] - mesh.aabb.min.z) / mesh.aabb.extent.z
    bowed = sampled.copy()
    bowed[:, 1] += 0.6 * height ** 2

    result = points_to_mesh(bowed, mesh, AlignParams3())
    aligned = result.alignment.full_transform.transform_points(bowed)
    deviations = mesh.measure_deviations(aligned, "plane")

    print(f"{len(aligned)} points, deviation from {deviations.min():.4f} to {deviations.max():.4f}, "
          f"rms {numpy.sqrt((deviations ** 2).mean()):.4f}")

    # The faces pointing towards +Y, as an `IndexMask`. The same mask type comes back from
    # `pick_faces` when a selection is made in the render window, so this is interchangeable with
    # one picked by hand.
    towards_y = mesh.face_select().facing(0.0, 1.0, 0.0, math.pi / 5, "add").to_mask()
    print(f"{towards_y.count_true} of {mesh.face_count} faces point towards +Y")

    plotter = Plotter(shape=(1, 4), window_size=(1600, 500))

    plotter.subplot(0, 0)
    plotter.engeom.draw_mesh(mesh, color="lightgray", opacity=0.4)
    # `clim` fixes the scale so the colors mean the same thing between runs, and the color map's
    # own over and under colors mark anything outside it.
    plotter.engeom.draw_point(aligned, scalars=deviations, cmap=GOM_CMAP, clim=(-0.2, 0.2),
                              point_size=8)
    plotter.add_text("deviation of measured points", font_size=10)

    plotter.subplot(0, 1)
    plotter.engeom.draw_mesh(mesh, color="lightgray", highlight=towards_y,
                             highlight_color="crimson")
    plotter.add_text("a face selection, one actor", font_size=10)

    plotter.subplot(0, 2)
    plotter.engeom.draw_mesh(mesh, color="white", opacity=0.25)
    plotter.engeom.draw_feature_edges(mesh, angle=40.0, color="black", line_width=1.0)
    plotter.add_text("feature edges", font_size=10)

    plotter.subplot(0, 3)
    plotter.engeom.draw_mesh(mesh, color="lightgray", opacity=0.5)
    # `pose` moves the actor rather than the mesh, so this second copy costs no geometry at all and
    # the mesh's own points are untouched.
    plotter.engeom.draw_mesh(mesh, color="steelblue", opacity=0.5,
                             pose=Iso3.from_translation(60.0, 0.0, 0.0))
    plotter.add_text("one mesh, two poses", font_size=10)

    plotter.link_views()
    plotter.view_isometric()
    plotter.show()


if __name__ == '__main__':
    main()
