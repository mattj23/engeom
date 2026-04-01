from pyvista import Plotter
from engeom.geom3 import Mesh
from engeom.plot import PyvistaPlotterHelper


def main():
    mesh = Mesh.stanford_bunny_res4()
    plotter = Plotter()
    helper = PyvistaPlotterHelper(plotter)
    helper.add_mesh(mesh)
    plotter.show()


if __name__ == '__main__':
    main()
