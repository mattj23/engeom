from pyvista import Plotter

from _common import DATA_DIR
from engeom.geom3 import Mesh3, Iso3, Plane3


def main():
    # Load a mesh of a small turbine blade to demonstrate. The blade mesh is dimensioned in millimeters and is roughly
    # aligned with the +Z direction pointing in the stacking axis and +X pointing in the engine axis direction towards
    # the front.
    mesh = Mesh3.load_tcmesh(DATA_DIR / "engine-blade.tcmesh")

    # We'll create a plane that is parallel to the XY plane and passes through the Z coordinate of the mesh's AABB
    # center. We'll then use the `section_with_plane` method to extract the curves that intersect the plane. They
    # come back as a `CurveGroup3`, which unpacks like a sequence.
    plane = Plane3.xy().offset_by(mesh.aabb.center.z)
    curves = mesh.section_with_plane(plane)

    # Finally, we'll plot the mesh, the plane we cut it with, and the curves that came back. A
    # `Plane3` has no origin or size of its own, so `draw_plane` needs to be told how much of it to
    # show. Left as it is here, the extent comes from whatever is already in the scene, which is
    # why the mesh is drawn first; the plane is then clipped to a box around it and drawn slightly
    # oversized so that it reads as passing through the part rather than stopping at its surface.
    plotter = Plotter()
    plotter.engeom.draw_mesh(mesh, color="white", opacity=0.5)
    plotter.engeom.draw_plane(plane, color="steelblue", opacity=0.25)
    plotter.engeom.draw_curve(*curves, color="black")
    plotter.add_axes()
    plotter.show()


if __name__ == '__main__':
    main()
