from _common import DATA_DIR
from engeom.airfoil2 import AfGeometry, OrientFwdAft, OrientUpperLower
from engeom.geom2 import Curve2
from engeom.plot.matplotlib import AxesHelper
from matplotlib.pyplot import Axes, Figure, figure, show


def main():
    # -------------------------------------------------------------------------------------------------------------
    # Note, this example does not use the specialized drawing tool for airfoils, but deliberately draws all the
    # geometry manually as a way of demonstrating the different parts of the airfoil analysis results.
    #
    # Everything between the inscribed circles and the LE/TE annotations below is what `draw_airfoil` does in a
    # single call, with the same default styling (see the docstring for more information on control):
    #     helper.draw_airfoil(geom)
    # -------------------------------------------------------------------------------------------------------------

    # This is a sample airfoil cross-section from a small metal hot section blade I bought on ebay some years ago. The
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
        helper.draw_circle(c.c, color="gray", alpha=0.4, linewidth=0.5)

    # Let's plot the tmax circle with a different color
    helper.draw_circle(geom.max_thickness_circle().c, color="red")

    # Now we'll draw the upper, lower, and camber curves. Upper will be red and lower blue.
    helper.draw_curve(geom.upper, color="lightsalmon", linewidth=3.0, label="upper")
    helper.draw_curve(geom.lower, color="cornflowerblue", linewidth=3.0, label="lower")
    helper.draw_curve(geom.camber, color="lightseagreen", linewidth=1.0, linestyle="--", label="camber")

    # We'll draw the original section over top of the curves so that we can verify that the upper and lower surfaces
    # match the original section data.
    helper.draw_curve(section, color="black", linewidth=1.0, linestyle="--")

    # Finally, we'll plot the leading and trailing edge points, and annotate them.
    le = geom.leading.point
    te = geom.trailing.point
    ax.scatter([le.x, te.x], [le.y, te.y], color="black", s=20, zorder=5)
    ax.annotate("LE", (le.x, le.y), textcoords="offset points", xytext=(6, 6))
    ax.annotate("TE", (te.x, te.y), textcoords="offset points", xytext=(6, 6))

    # Now we'll add a pair of thickness gage dimensions, one near each edge. These use the "on_camber" method, which
    # locates the gage point by walking a given arc distance along the mean camber line and then casting orthogonally
    # to the camber to find the two surfaces. It's the most commonly used of the three location methods and behaves
    # well anywhere along the airfoil. A positive value is measured from the leading edge and a negative one from the
    # trailing edge, so the two stations below sit 3 mm aft of the LE and 5 mm forward of the TE on a camber line
    # about 42.8 mm long. Both return None if the station doesn't land on both surfaces, which is worth checking
    # rather than assuming, since a value longer than the camber line will simply run off the end.
    le_gage = geom.thickness_at("on_camber", 3.0)
    te_gage = geom.thickness_at("on_camber", -5.0)

    # The measurements are ordinary `Distance2` objects, so the standard dimension drawing handles them. We push the
    # labels out to opposite sides so they don't collide with the section, using "outside_rev" at the leading edge to
    # place its label forward instead of aft.
    if le_gage is not None:
        helper.draw_distance(le_gage, label_place="outside_rev", template="LE+3.0: {value:.3f}")
    if te_gage is not None:
        helper.draw_distance(te_gage, label_place="outside", template="TE-5.0: {value:.3f}")

    print(f"Leading edge: kind={geom.leading.geometry.kind}, point=({le.x:.4f}, {le.y:.4f})")
    print(f"Trailing edge: kind={geom.trailing.geometry.kind}, point=({te.x:.4f}, {te.y:.4f})")
    print(f"Camber length: {geom.camber.length():.4f}, max thickness: {geom.max_thickness().value:.4f}")
    if le_gage is not None:
        print(f"Thickness 3.0 aft of the LE: {le_gage.value:.4f}")
    if te_gage is not None:
        print(f"Thickness 5.0 forward of the TE: {te_gage.value:.4f}")

    ax.legend(loc="best")
    show()


if __name__ == "__main__":
    main()
