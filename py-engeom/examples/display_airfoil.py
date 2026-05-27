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

    circles = []
    p0s = []
    p1s = []
    with open("/tmp/circles.txt", "r") as handle:
        for line in handle.readlines():
            x, y, r, p0x, p0y, p1x, p1y = [float(q) for q in line.split()]
            circles.append(Circle2(x, y, r))
            p0s.append([p0x, p0y])
            p1s.append([p1x, p1y])
    p0s = numpy.array(p0s)
    p1s = numpy.array(p1s)

    for c in circles[1:]:
        helper.plot_circle(c, color="red")

    ax.scatter(p0s[:, 0], p0s[:, 1])
    ax.scatter(p1s[:, 0], p1s[:, 1])


    helper.plot_curve(curve, color="black")
    show()


if __name__ == '__main__':
    main()
