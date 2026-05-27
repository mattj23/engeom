import numpy

from _common import DATA_DIR
from engeom.geom2 import Curve2, Circle2
from matplotlib.pyplot import figure, Figure, Axes, show
from engeom.plot import MatplotlibAxesHelper


def main():
    # Load the 2d curve of a small turbine blade airfoil cross-section.
    curve = Curve2.load_tccurve2(DATA_DIR / "airfoil-0.tccurve2")

    fig: Figure = figure(figsize=(10, 6))
    ax: Axes = fig.subplots()
    helper = MatplotlibAxesHelper(ax)

    helper.plot_curve(curve, color="black")
    show()


if __name__ == '__main__':
    main()
