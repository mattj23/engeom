"""
    This example demonstrates fitting a 2D cubic Bezier spline to a set of noisy points using
    `fit_spline_to_points`, and then using the spline's closest-point query to measure how far a
    handful of probe points lie from the fitted curve.

    This is representative of a 2D profile metrology workflow: you have measured surface points that
    follow some smooth curve, and you want to recover a clean analytical curve through them. Once you
    have the curve, the query structure lets you project arbitrary points onto it (for deviation
    measurement, segmenting, etc.).

    A single cubic Bezier segment is defined by four control points `p0, p1, p2, p3`. Here we treat
    the profile endpoints as known and hold them fixed, letting the optimizer move only the two
    interior control points. The parameter vector is therefore four values:

      params[0]  p1x   x-coordinate of the second control point
      params[1]  p1y   y-coordinate of the second control point
      params[2]  p2x   x-coordinate of the third control point
      params[3]  p2y   y-coordinate of the third control point

    The builder function turns those four values (plus the fixed endpoints) into a `CubicSpline2`,
    and `fit_spline_to_points` runs a Levenberg-Marquardt minimization driving the curve towards the
    sample points. The residual for each point is its distance to the closest point on the curve.

    A note on what gets "recovered": the control points of a cubic Bezier are not uniquely
    determined by points sampled along it, because the curve can be re-parameterized without
    changing its shape. So the fitted interior control points will generally differ from the ones
    that generated the data even though the fitted *curve* tracks the true curve closely. That is
    why this example reports the fit residual and overlays the curves rather than comparing control
    points directly.
"""

from typing import Callable

import numpy
from numpy.typing import NDArray
from matplotlib.pyplot import Figure, Axes, figure, show as show_figure

from engeom.geom2 import Aabb2, CubicSpline2, Point2, Segment2, fit_spline_to_points
from engeom.plot.matplotlib import AxesHelper


def make_builder(p0: Point2, p3: Point2) -> Callable[[NDArray], CubicSpline2]:
    """
    Create the builder function that `fit_spline_to_points` will call. The endpoints are captured
    here and held fixed, so the optimizer only sees the four interior-control-point parameters.
    """

    def build_spline(params: NDArray) -> CubicSpline2:
        return CubicSpline2(
            p0.x, p0.y,
            params[0], params[1],
            params[2], params[3],
            p3.x, p3.y,
        )

    return build_spline


def main():
    rng = numpy.random.default_rng(24601)

    # Start by making the true spline so that we can generate points off of it and add random noise.
    # Normally you'd be working with data taken off of some sort of measurement system.
    # ----------------------------------------------------------------------------------------------
    true_spline = CubicSpline2(0.0, 0.0, 1.0, 2.0, 2.0, 2.0, 3.0, 0.0)

    ts = numpy.linspace(0.0, 1.0, 60)
    clean_points = numpy.array([[true_spline.position(t).x, true_spline.position(t).y] for t in ts])

    # Add Gaussian noise with σ = 0.02.
    noisy_points = clean_points + rng.normal(0.0, 0.02, clean_points.shape)

    # Perform the fitting. The profile endpoints are treated as known and pinned, and the interior
    # control points start from a straight line between them (an intentionally poor guess that the
    # optimizer has to bend into shape). The initial guess must be a float array.
    # ----------------------------------------------------------------------------------------------
    p0 = Point2(true_spline.p0.x, true_spline.p0.y)
    p3 = Point2(true_spline.p3.x, true_spline.p3.y)
    builder = make_builder(p0, p3)

    initial_guess = numpy.array([
        p0.x + (p3.x - p0.x) / 3.0, p0.y + (p3.y - p0.y) / 3.0,
        p0.x + 2.0 * (p3.x - p0.x) / 3.0, p0.y + 2.0 * (p3.y - p0.y) / 3.0,
    ])

    result = fit_spline_to_points(noisy_points, builder, initial_guess)
    fitted_spline = builder(result)

    # Report how well the fitted curve matches: the residual of the sample points, and how far the
    # fitted curve sits from the original true curve. The query structure makes both easy.
    # ----------------------------------------------------------------------------------------------
    queries = fitted_spline.query()

    sample_dist = queries.project_points(noisy_points)[:, 1]
    rms = float(numpy.sqrt(numpy.mean(sample_dist ** 2)))

    dense_true = numpy.array([[true_spline.position(t).x, true_spline.position(t).y]
                              for t in numpy.linspace(0.0, 1.0, 400)])
    curve_dev = queries.project_points(dense_true)[:, 1]

    print(f"Sample residual : rms = {rms:.4f}, max = {sample_dist.max():.4f}")
    print(f"Curve deviation : max = {curve_dev.max():.4f}  (fitted curve vs. true curve)")

    # Now use the closest-point query to measure a few probe points against the fitted curve. Build
    # the query structure once and reuse it for every probe.
    # ----------------------------------------------------------------------------------------------
    probes = numpy.array([
        [1.5, 2.0],
        [2.5, 1.2],
        [0.5, 0.4],
    ])

    print("\nProbe point        t       distance")
    projections = []
    for probe in probes:
        proj = queries.project_point(Point2(probe[0], probe[1]))
        foot = fitted_spline.position(proj.t)
        projections.append((probe, foot))
        print(f"  ({probe[0]:5.2f}, {probe[1]:5.2f})    {proj.t:6.3f}   {proj.distance:8.4f}")

    # Finally, plot the results.
    # ----------------------------------------------------------------------------------------------
    fig: Figure = figure(figsize=(9, 6))
    ax: Axes = fig.subplots()
    helper = AxesHelper(ax)

    # `draw_point` takes the whole (n, 2) array directly, so the columns don't need splitting out.
    helper.draw_point(noisy_points, markersize=4.0, color="steelblue", zorder=3,
                      label=f"Sample points (σ = 0.02, n = {len(noisy_points)})")

    # `draw_spline` flattens each spline into a polyline itself, to the tolerance given, so there's
    # no need to call `polyline` and unpack the columns by hand.
    helper.draw_spline(true_spline, tol=0.001, color="green", linewidth=2.0, linestyle="--",
                       label="True spline")
    helper.draw_spline(fitted_spline, tol=0.001, color="firebrick", linewidth=2.0, linestyle="-",
                       label="Fitted spline")

    # Draw the closest-point projections as leader lines from each probe to its foot on the curve.
    # A `Segment2` is the natural way to say "a line between these two points".
    helper.draw_segment(*(Segment2(probe[0], probe[1], foot.x, foot.y) for probe, foot in projections),
                        color="black", linewidth=1.0, zorder=4)

    # All the probe markers go into one artist, which is also what gives them a single legend entry
    # instead of needing the label suppressed on all but the first.
    helper.draw_point(probes, marker="o", markersize=6.0, color="darkorange", zorder=5,
                      label="Probe points")

    ax.legend()
    ax.set_title("fit_spline_to_points + closest-point query")

    # Frame everything that was drawn, probes included, or the leader to the probe above the curve
    # runs off the top of the view.
    bounds = Aabb2.from_points(noisy_points).merged(Aabb2.from_points(probes)).expand(0.5)
    helper.set_bounds(bounds)

    fig.tight_layout()
    show_figure()


if __name__ == "__main__":
    main()
