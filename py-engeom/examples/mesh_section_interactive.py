"""
    `mesh_section.py` cuts a mesh at one plane worked out ahead of time. This example turns the
    same operation into an instrument: the plane is a widget you drag, and the section is recut and
    redrawn as it moves.

    Two things make that work, and neither is specific to sectioning:

    * The widget reports a `Plane3`, so the callback is ordinary `engeom` code rather than anything
      to do with the render window.
    * Redrawing under a fixed actor `name` replaces the previous actor instead of adding another,
      which is what keeps a callback that fires on every mouse move from filling the scene with
      dead geometry.

    The same file also switches on face picking. Drag a rubber band over the part with the left
    mouse button and the faces you enclose come back as an `IndexMask` over the mesh, which is the
    same currency `Mesh3.face_select` produces and `extract_subset_faces` consumes. That is the
    round trip worth understanding: something chosen by hand in the window becomes a value the Rust
    side can act on.

    This one needs a real window, so unlike the other examples there is nothing to see if it is run
    off-screen.
"""

from pyvista import Plotter

from _common import DATA_DIR
from engeom.geom3 import Mesh3


def main():
    mesh = Mesh3.load_tcmesh(DATA_DIR / "engine-blade.tcmesh")

    plotter = Plotter()

    # Several translucent surfaces overlap here, the part and the widget's own plane among them, and
    # they are drawn in the order they were added rather than by depth. `plotter.enable_depth_peeling()`
    # sorts that out properly and is worth adding if the transparency looks wrong, but it needs real
    # graphics hardware: on a software renderer it can leave the translucent geometry drawing as
    # nothing at all, so it is left out here rather than making the example fragile.
    plotter.engeom.draw_mesh(mesh, color="white", opacity=0.35, name="part")

    def recut(plane):
        """Recut the mesh wherever the widget has been dragged to."""
        try:
            curves = mesh.section_with_plane(plane)
        except Exception:
            # The plane can be dragged clear of the part, where there is nothing to cut. Leaving
            # the last section on screen would be a lie, so it is removed instead.
            plotter.remove_actor("section")
            return
        plotter.engeom.draw_curve(*curves, color="crimson", line_width=4, name="section")

    plotter.engeom.plane_widget(recut, normal=(0.0, 0.0, 1.0))

    def report(mask):
        """Report whatever was just selected with the rubber band."""
        if mask.count_true == 0:
            print("nothing selected")
            return
        picked = mesh.extract_subset_faces(mask)
        area = picked.compute_face_areas().sum()
        print(f"selected {mask.count_true} faces, {area:.2f} mm^2, "
              f"in {len(picked.separate_patches())} connected patch(es)")

    plotter.engeom.pick_faces(mesh, report)

    plotter.add_axes()
    plotter.add_text("drag the plane to recut; rubber-band select faces with the left button",
                     font_size=9)
    plotter.show()


if __name__ == '__main__':
    main()
