"""
    This example shows the creation and plotting of a simple open 2D boundary. We'll use the direct `BoundaryData2`
    interface to create the boundary, and then plot it using Matplotlib. The boundary consists of a line segment,
    followed by an arc, and then another line segment.
"""
import matplotlib.pyplot
from matplotlib.pyplot import figure, Figure, Axes

from engeom.geom2 import BoundaryData2
from engeom.plot.matplotlib import AxesHelper


def main():
    # We'll create an open boundary, which we do by providing the x, y coordinates of the first point. If we wanted
    # a closed boundary we would not provide the points, but instead would set the `closed=True` argument in the
    # initalizer. There are also the more explicit `BoundaryData2.new_open(x, y)` and `BoundaryData2.new_closed()`
    # static constructor methods.
    data = BoundaryData2(0, 0)

    # Now we'll add the segment, the arc, and the final segment
    data.add_seg_xy(1, 0)
    data.add_arc_xy(1, 0.5, 1, 1, False)
    data.add_seg_xy(0, 1)

    # We use the `to_boundary()` method to generate actual queriable geometry from the boundary definition
    boundary = data.to_boundary()

    # Now we'll create a matplotlib figure and an engeom helper to plot the geometry
    fig: Figure = figure(figsize=(8, 5))
    ax: Axes = fig.subplots()
    helper = AxesHelper(ax)

    # We'll use the helper's boundary drawing function to plot the geometry. The `tol` argument controls how densely
    # points are generated over arc geometry. The rest of the arguments are passthroughs to the matplotlib `plot`
    # function.
    helper.boundary(boundary, color="black", linewidth=2.0, linestyle="--", tol=0.001)

    # We'll additionally plot the boundary's normals. This is useful for debugging when trying to make sure that the
    # winding direction correctly matches your intention.
    helper.boundary_normals(boundary, 10, 0.1, color="red")

    # Finally, because the annotation arrows used to draw the normals don't automatically adjust the plot bounds,
    # we'll grab the `Aabb2` of the boundary, expand it by 0.5 in all directions, and then use it to set the plot
    # boundary through the helper's `.set_bounds` method.
    helper.set_bounds(boundary.aabb.expand(0.5))

    matplotlib.pyplot.show()


if __name__ == '__main__':
    main()
