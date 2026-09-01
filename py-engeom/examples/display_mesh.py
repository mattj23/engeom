from pyvista import Plotter
from engeom.geom3 import Mesh3


def main():
    mesh = Mesh3.stanford_bunny_res4()

    # Note that nothing from `engeom.plot` is imported here. `engeom` declares its PyVista helper as
    # a plugin, so PyVista attaches it to every plotter as `plotter.engeom` and imports it the first
    # time that attribute is read.
    plotter = Plotter()
    plotter.engeom.draw_mesh(mesh)
    plotter.show()


if __name__ == '__main__':
    main()
