"""
    This example demonstrates fitting a boundary to a set of noisy 2D points using `fit_boundary_to_points`. We fit a
    slot shape (a rectangle with semicircular ends) to points sampled from a known slot with added Gaussian noise.

    This is representative of a 2D cross-section metrology workflow where you want to recover the intent of the nominal
    geometry from measured surface points, sometimes for the purpose of direct measurement, other times just to create
    a reference for segmenting data in preparation for more intense measurement methods.

    The slot is parameterized by four values:
      params[0]  cx        x-coordinate of the slot centre
      params[1]  cy        y-coordinate of the slot centre
      params[2]  half_len  half the distance between the two arc centres
      params[3]  radius    radius of the semicircular ends (= half the slot width)

    The builder function creates a closed CCW-wound BoundaryData2 from those four values. Then `fit_boundary_to_points`
    runs a Levenberg-Marquardt minimization that drives the boundary towards the sample points.
"""

import numpy
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, Axes

from engeom.geom2 import BoundaryData2, fit_boundary_to_points
from engeom.plot.matplotlib import AxesHelper


def build_slot(params: NDArray) -> BoundaryData2:
    """
    Build a closed slot (rectangle + two semicircles) from params = [cx, cy, half_len, radius]. The boundary is wound
    counter-clockwise so that surface normals point outward, but typically a slot would have the normals facing
    inward.

    This is the function that will get injected into `fit_boundary_to_points`, which will then try to optimize the
    parameters based on the geometry this creates.
    """
    cx, cy, half_len, radius = params[0], params[1], params[2], params[3]

    # A closed boundary has no explicit starting point; the first element's start is
    # the implicit wrap-around from the last element.
    data = BoundaryData2(closed=True)

    # Bottom segment: (cx-half_len, cy-radius) → (cx+half_len, cy-radius)
    data.add_seg_xy(cx + half_len, cy - radius)

    # Right semicircle: CCW arc around (cx+half_len, cy) up to (cx+half_len, cy+radius)
    data.add_arc_xy(cx + half_len, cy, cx + half_len, cy + radius, False)

    # Top segment: (cx+half_len, cy+radius) → (cx-half_len, cy+radius)
    data.add_seg_xy(cx - half_len, cy + radius)

    # Left semicircle: CCW arc around (cx-half_len, cy) down to (cx-half_len, cy-radius)
    data.add_arc_xy(cx - half_len, cy, cx - half_len, cy - radius, False)

    return data


def main():
    rng = numpy.random.default_rng(24601)

    # We'll start by making the true slot geometry so that we can generate points off of it and add random noise to
    # them. Normally you'd be working with data taken off of some sort of measurement system.
    # ----------------------------------------------------------------------------------------------------------------
    true_params = numpy.array([2.0, 1.0, 3.0, 1.5])
    true_boundary = build_slot(true_params).to_boundary()

    lengths = numpy.linspace(0, true_boundary.length, 200, endpoint=False)
    clean_points = true_boundary.at_lengths(lengths)[:, :2]

    # Add Gaussian noise with σ = 0.05 (5 % of the slot radius).
    noisy_points = clean_points + rng.normal(0.0, 0.05, clean_points.shape)

    # Now we're going to perform the fitting. We'll start with an intial guess that's intentially displaced and
    # not sized correctly. The fitting function will minimize the deviation.
    # ----------------------------------------------------------------------------------------------------------------
    initial_guess = numpy.array([0.0, 0.0, 2.0, 1.0])

    result = fit_boundary_to_points(noisy_points, build_slot, initial_guess)
    fitted_boundary = build_slot(result).to_boundary()

    print("Parameter       True     Initial   Fitted")
    print(f"  cx          {true_params[0]:7.3f}  {initial_guess[0]:7.3f}  {result[0]:7.3f}")
    print(f"  cy          {true_params[1]:7.3f}  {initial_guess[1]:7.3f}  {result[1]:7.3f}")
    print(f"  half_len    {true_params[2]:7.3f}  {initial_guess[2]:7.3f}  {result[2]:7.3f}")
    print(f"  radius      {true_params[3]:7.3f}  {initial_guess[3]:7.3f}  {result[3]:7.3f}")

    # Finally, for fun, we're going to plot the results
    # ----------------------------------------------------------------------------------------------------------------
    fig: Figure = plt.figure(figsize=(9, 6))
    ax: Axes = fig.subplots()
    helper = AxesHelper(ax)

    # `draw_point` takes the whole (n, 2) array directly, so the columns don't need splitting out.
    helper.draw_point(noisy_points, markersize=4.0, color="steelblue", zorder=3,
                      label=f"Sample points (σ = 0.05, n = {len(noisy_points)})")

    helper.draw_boundary(true_boundary, color="green", linewidth=2.0, linestyle="--", tol=0.005, label="True boundary")
    helper.draw_boundary(fitted_boundary, color="firebrick", linewidth=2.0, linestyle="-", tol=0.005, label="Fitted boundary")

    ax.legend()
    ax.set_title("fit_boundary_to_points: slot fitting")
    helper.set_bounds(true_boundary.aabb.expand(0.75))

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
