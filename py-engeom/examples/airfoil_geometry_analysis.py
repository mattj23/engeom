from _common import DATA_DIR
from engeom.airfoil2 import AfGeometry, OrientFwdAft, OrientUpperLower
from engeom.geom2 import Curve2
from engeom.plot.matplotlib import AxesHelper
from matplotlib.pyplot import Axes, Figure, figure, show


def main():
    # This is a sample airfoil cross-section from a small metal hot section blade found on ebay some years ago.  The
    # original data is from a scan of the blade on a Zeiss (formerly GOM) ATOS 5 3D scanner, and the cross-section was
    # taken after some mesh smoothing and error handling. Units are in millimeters.
    section = Curve2.load_tccurve2(DATA_DIR / "airfoil-0.tccurve2")

    # We're going to establish the baseline airfoil geometry (mean camber line, upper/lower surfaces, and
    # leading/trailing points and geometries) purely from the geometry, since I do not have a CAD model for this part.
    # Because this data is from an extremely high quality scan, this will work. Typically, you should use the
    # pure geometric analysis only on nominal data, and then use the hybrid analysis on cross-sections from
    # actual parts.
    #
    # For parameters, we'll allow the algorithm to detect the leading edge orientation by seeing which side the
    # maximum thickness (T-max) is biased towards.  The blade also has plenty of curvature, so we can use the
    # curvature based method to automatically identify the upper and lower skins. We'll allow the auto search to
    # try to identify the leading and trailing edge geometries.
    geom = AfGeometry.from_geometric_analysis(
        section=section,
        general_tol=1e-3,
        fwd_aft=OrientFwdAft.tmax_fwd(),
        upper_lower=OrientUpperLower.curvature(),
        le_search="auto",
        te_search="auto",
    )

    # For plotting, we'll use matplotlib, and `engeom`'s helper tools.
    fig: Figure = figure(figsize=(12, 6))
    ax: Axes = fig.subplots()
    helper = AxesHelper(ax)

    # We'll draw the inscribed circles first, making them thin and slighty transparent
    for c in geom.circles:
        helper.plot_circle(c.c, color="gray", alpha=0.4, linewidth=0.5)

    # Let's plot the tmax circle with a different color
    helper.plot_circle(geom.max_thickness_circle().c, color="red")

    # Now we'll draw the upper, lower, and camber curves. Upper will be red and lower blue.
    helper.plot_curve(geom.upper, color="lightsalmon", linewidth=3.0, label="upper")
    helper.plot_curve(geom.lower, color="cornflowerblue", linewidth=3.0, label="lower")
    helper.plot_curve(geom.camber, color="lightseagreen", linewidth=1.0, linestyle="--", label="camber")

    # We'll draw the original section over top of the curves so that we can verify that the upper and lower surfaces
    # match the original section data.
    helper.plot_curve(section, color="black", linewidth=1.0, linestyle="--")

    # Finally, we'll plot the leading and trailing edge points, and annotate them.
    le = geom.leading.point
    te = geom.trailing.point
    ax.scatter([le.x, te.x], [le.y, te.y], color="black", s=20, zorder=5)
    ax.annotate("LE", (le.x, le.y), textcoords="offset points", xytext=(6, 6))
    ax.annotate("TE", (te.x, te.y), textcoords="offset points", xytext=(6, 6))

    print(f"Leading edge: kind={geom.leading.geometry.kind}, point=({le.x:.4f}, {le.y:.4f})")
    print(f"Trailing edge: kind={geom.trailing.geometry.kind}, point=({te.x:.4f}, {te.y:.4f})")

    ax.legend(loc="best")
    show()


if __name__ == "__main__":
    main()
