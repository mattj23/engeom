# Camber Line

Most airfoil cross-section analysis begins with the detection of the mean camber line. The mean camber line is a
theoretical line which runs from an intersection at the center of the leading edge to an intersection at the center of
the trailing edge, all while staying equidistant from the upper and lower surfaces of the cross-section.

## Importance of the Mean Camber Line for Airfoil Analysis

The mean camber line (MCL) is very closely related to the geometric concept of
the [medial axis](https://en.wikipedia.org/wiki/Medial_axis), which is part of an alternate way to represent the shape
of a 2D object by defining a "skeleton" like graph structure and a set of widths at each point along the skeleton. In
the case of an airfoil, the MCL can be thought of as a "skeleton" of an airfoil cross-section, but extended forward and
backwards to intersect with the leading and trailing edges.

For smooth airfoil geometry (having no sharp corners or discontinuities) the MCL and the medial axis are the same
outside the leading and trailing edge regions, and the MCL is redundant with the medial axis within those regions.
Geometry like this is common in the idealized airfoil shapes created during design, and is typical of CAD and other
design geometry.

On idealized geometry, the mean camber line is usually very easy to compute and will provide a large amount of
information about the airfoil shape, including its orientation, it's direction of curvature, the location of the leading
and trailing edges, and the location and value of maximum thickness. At the very least this calculation, which can be
done with very little more than the set of points defining the cross-section (and a few extra hints on non-typical
airfoils), is an important tool for establishing a frame of reference to locate other measurements and features, even if
the MCL is not used for any other purpose.

Actual, as-built airfoils tend to have non-negligible bumps, pits, and other high curvature regions which occur during
normal manufacturing processes and will be captured during measurement. Depending on the quality of the manufacturing
and the measurement system, these features may prevent the reliable calculation of MCL from just the naive cross-section
points. However, by using data from the nominal cross-section and nominal MCL to orient and provide a frame of
reference, actual MCLs can be calculated even on very warped/distorted airfoils with noisy measurement data.

Finally, because of the shape-representing properties of the MCL which it inherits from the medial axis concept, the
combination of MCL and thickness distribution of an as-built airfoil can be thought of as an alternate, reduced
representation of that airfoil which preserves more aerodynamically relevant features (like thicknesses, curvatures,
orientation, and size) while filtering out some of the measurement and manufacturing surface artifacts that might
otherwise frustrate measurement techniques.

To summarize, the MCL is often the first step in airfoil analysis because it:

1. Calculates a number of aerodynamically relevant features of the airfoil
2. Provides an unambiguous frame of reference by which the location of other features and measurements can be
   established
3. Provides a reduced representation of an as-built airfoil which is more aerodynamically relevant and less sensitive to
   measurement artifacts

## Engeom 