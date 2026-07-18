from pyvista import Plotter

from _common import DATA_DIR
from engeom.geom3 import Mesh, Iso3, Plane3
from engeom.plot import PyvistaPlotterHelper


def main():
    # Load a mesh of a small turbine blade to demonstrate. The blade mesh is dimensioned in millimeters and is roughly
    # aligned with the +Z direction pointing in the stacking axis and +X pointing in the engine axis direction towards
    # the front.
    mesh = Mesh.load_umesh(DATA_DIR / "engine-blade.umesh.gz")

    # We'll create a plane that is parallel to the XY plane and passes through the Z coordinate of the mesh's AABB
    # center. We'll then use the `section` method to extract the curves that intersect the plane.
    plane = Plane3.xy().offset_by(mesh.aabb.center.z)
    curves = mesh.section(plane)

    # Finally, we'll plot the original points, the aligned points, and the original mesh.
    plotter = Plotter()
    helper = PyvistaPlotterHelper(plotter)
    helper.mesh(mesh, color="white", opacity=0.5)
    helper.curves(*curves, color="black")
    plotter.add_axes()
    plotter.show()


if __name__ == '__main__':
    main()
