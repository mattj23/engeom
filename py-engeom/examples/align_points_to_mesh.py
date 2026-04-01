import math

from engeom.engeom import SelectOp, DeviationMode
from pyvista import Plotter
from engeom.geom3 import Mesh, Iso3
from engeom.align import points_to_mesh
from engeom.plot import PyvistaPlotterHelper

from _common import DATA_DIR


def main():
    # Load a mesh of a small turbine blade to demonstrate. The blade mesh is dimensioned in millimeters and is roughly
    # aligned with the +Z direction pointing in the stacking axis and +X pointing in the engine axis direction towards
    # the front.
    mesh = Mesh.load_umesh(DATA_DIR / "engine-blade.umesh.gz")

    # We're going to grab a sub-mesh consisting of the faces that are rougly pointed towards +Y and then generate
    # points from it using a poisson disk sample with radius of 2mm. Note that the sample points include mesh normals,
    # so the 3d coordinates are in `sample_points[:, :3]`.
    sub_mesh = mesh.face_select_none().facing(0, 1, 0, math.pi/4, SelectOp.Add).create_mesh()
    sample_points = sub_mesh.sample_poisson(2)

    # We'll clone the points and move them away from the original mesh to have the alignment do something. In
    # reality, the points would have come from some other source, like a 3d scanner.
    disturb = Iso3.from_translation(-100, 150, 0) @ Iso3.from_rotation(math.pi / 6, 1, 1, 1)
    to_align = disturb.transform_points(sample_points[:, :3])

    # Now we perform the alignment. If the result is succesful we'll get an `Iso3` back, otherwise the call to
    # `points_to_mesh` will throw an error. We'll create a new set of points by transforming `to_align` to represent
    # where the points are with the alignment applied. The alignment itself does not mutate any objects.
    result = points_to_mesh(to_align, mesh, Iso3.identity(), DeviationMode.Point)
    aligned = result.transform_points(to_align)

    # Finally, we'll plot the original points, the aligned points, and the original mesh.
    plotter = Plotter()
    helper = PyvistaPlotterHelper(plotter)
    helper.add_mesh(mesh, color="white")
    plotter.add_points(to_align, point_size=5, color="red")
    plotter.add_points(aligned, point_size=5, color="green")
    plotter.add_axes()
    plotter.add_text("Original points are in red, aligned points are in green", font_size=10, font="courier")
    plotter.show()


if __name__ == '__main__':
    main()
