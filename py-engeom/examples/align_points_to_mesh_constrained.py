import math

from engeom.engeom import DeviationMode
from pyvista import Plotter

from engeom.align import points_to_mesh, Dof6
from engeom.geom3 import Mesh, Iso3
from engeom.plot import PyvistaPlotterHelper


def main():
    # Load a mesh of a small turbine blade to demonstrate. The blade mesh is dimensioned in millimeters and is roughly
    # aligned with the +Z direction pointing in the stacking axis and +X pointing in the engine axis direction towards
    # the front.
    mesh = Mesh.create_box(10, 5, 2)
    sample_points = mesh.sample_poisson(0.5)

    disturb = Iso3.from_translation(0.5, 1, 1) @ Iso3.from_rotation(-math.pi / 12, 1, 1, 1)
    to_align = disturb.transform_points(sample_points[:, :3])

    # Now we perform the alignment. If the result is successful, we'll get an `Iso3` back, otherwise the call to
    dof = Dof6(tx=False, ty=True, tz=True, rx=True, ry=True, rz=True)
    result = points_to_mesh(to_align, mesh, mode=DeviationMode.Point, dof=dof)
    aligned = result.transform_points(to_align)

    # Finally, we'll plot the original points, the aligned points, and the original mesh.
    plotter = Plotter()
    helper = PyvistaPlotterHelper(plotter)
    helper.mesh(mesh, color="white")
    plotter.add_points(to_align, point_size=5, color="red")
    plotter.add_points(aligned, point_size=5, color="green")
    plotter.add_axes()
    plotter.add_text("Original points are in red, aligned points are in green", font_size=10, font="courier")
    plotter.show()


if __name__ == '__main__':
    main()
