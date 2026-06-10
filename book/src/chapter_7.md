# Airfoils (2D)

The `engeom` library has an entire module dedicated to performing dimensional analysis on 2D airfoil cross-sections.

In the aerospace industry, common measurements taken on airfoils include:

- Thickness measurements through the airfoil, especially around the leading/trailing edges and
  the position of maximum thickness, usually located by position along the mean camber line or
  in local directions at the edges.
- Chord lengths, often specified by a variety of different methods.
- Form and profile measurements, often with multiple tolerance regions, sometimes requiring
  special rules for partially constrained floating zones
- Section position and angle, measured from nominal references
- Leading and trailing edge position and shape, such as leading edge radius or trailing edge
  trim position, etc.

> [!IMPORTANT]
> It's important to start with the understanding that there's a huge difference between
> running airfoil shape analysis tools on nominal section data, such as that exported from a CAD
> system or a mathematical representation of design geometry, and running the same tools on
> measured section data that came from a system like a 3D scanner or CMM.  Actual data will have
> both noise from the measurement system and actual defects/roughness from the manufacturing
> process. This brings a whole collection of problems that range from the practical to the
> philosophical, all which need to be addressed at some level.

## Typical Analysis Sequence

Before any airfoil-specific measurements can be taken on a section, the core airfoil geometry needs to be established. 
The mean camber line needs to be identified, the leading and/or trailing edge points need to be located, and the upper
and lower sides of the airfoil need to be partitioned.

In the `engeom` library there are two main ways of going about this, depending on whether or not you have nominal
reference geometry to serve as a starting point.

- If you _don't_ have nominal reference geometry it is either because you don't have a CAD model, or because the section
  you're working with _is_ the nominal reference geometry. In this case, the airfoil geometry is established using
  purely geometric algorithms.
- If you _do_ have nominal reference geometry, the `engeom` library has tools to take a hybrid approach which uses the
  nominal geometry to establish were features should be, and then performs more robust local searching to capture the
  actual features.

## Geometry Extraction

Pure geometric approach:

1. Find the unambiguous inscribed circles
2. Order the inscribed circles from leading edge to trailing edge
3. Orient the inscribed circles so that their first contact point is on the lower surface and the second contact point 
   is on the upper surface
4. Identify the leading edge geometry
5. Identify the trailing edge geometry
6. Construct the full mean camber line through the inscribed circle centers and the leading/trailing edge geometry
7. Partition the upper and lower surfaces using the mean camber line and the leading/trailing

Hybrid approach:
1. Best fit the actual section to the nominal reference geometry
2. Transfer a central inscribed circle from the nominal geometry to the actual section, then align them
3. 


## Mean Camber Line

All measurements, airfoil or otherwise, require some method of specifying _where_ the measurement should be taken. On 
older/legacy airfoil components this was sometimes done without the need for any special understanding of airfoil 
geometry.  Over the past few decades, however, drawings have shifted to measurements that require isolating the upper 
and lower surfaces and identifying the leading/trailing edges in order to specify measurement locations and tolerance
regions.

The main landmark for identifying places on an airfoil is the mean camber line (MCL), which is the curve that runs from 
the leading edge to the trailing edge and is equidistant from the upper and lower surfaces. The mean camber line is a 
portion of the [medial axis](https://en.wikipedia.org/wiki/Medial_axis), a geometric concept with a strong mathematical
foundation.

In the `engeom` library, any airfoil-specific measurement tools first require the establishment of the MCL. Once 
established, the rest of the airfoil can be partitioned into upper and lower surfaces, the leading and trailing edges 
identified, and any point along the airfoil can be specified using the MCL and the side of interest.






