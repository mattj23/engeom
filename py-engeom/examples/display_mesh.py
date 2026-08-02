from pyvista import Plotter
from engeom.geom3 import Mesh3
from engeom.plot.pyvista import PlotterHelper


def main():
    mesh = Mesh3.stanford_bunny_res4()
    plotter = Plotter()
    helper = PlotterHelper(plotter)
    helper.draw_mesh(mesh)
    plotter.show()


if __name__ == '__main__':
    main()
