import math

from pyvista import Plotter

from engeom.align3 import AlignParams3, Dof6, points_to_mesh
from engeom.geom3 import Mesh3, Iso3
from engeom.plot.pyvista import PlotterHelper


def main():
    # A simple box stands in for a real part here, so that the effect of locking a degree of freedom is easy to see
    # against the box's flat faces. The points to align are a poisson disk sample of the box's own surface.
    mesh = Mesh3.create_box(10, 5, 2)
    sample_points = mesh.sample_poisson(0.5).points

    disturb = Iso3.from_translation(0.5, 1, 1) @ Iso3.from_rotation(-math.pi / 12, 1, 1, 1)
    to_align = disturb.transform_points(sample_points)

    # Alignment parameters carry the degrees of freedom the solver is allowed to use. Locking `tx` leaves translation
    # along the world X axis untouched, so the 0.5mm of the disturbance along X survives the alignment while
    # everything else is corrected.
    params = AlignParams3(dof=Dof6(tx=False))
    result = points_to_mesh(to_align, mesh, params)
    aligned = result.alignment.full_transform.transform_points(to_align)

    # Finally, we'll plot the original points, the aligned points, and the original mesh. The other
    # examples reach the drawing methods through `plotter.engeom`, which is the same helper attached
    # to the plotter. Building one directly, as here, is the way to draw on a PyVista older than
    # 0.48, which has no mechanism for attaching it.
    plotter = Plotter()
    helper = PlotterHelper(plotter)
    helper.draw_mesh(mesh, color="white")
    plotter.add_points(to_align, point_size=5, color="red")
    plotter.add_points(aligned, point_size=5, color="green")
    plotter.add_axes()
    plotter.add_text("Original points are in red, aligned points are in green", font_size=10, font="courier")
    plotter.show()


if __name__ == '__main__':
    main()
