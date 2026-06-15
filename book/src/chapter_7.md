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

## Airfoil Geometry Primer

An airfoil is an aerodynamic body with a shape designed to create lift in a flowing fluid. Aircraft wings are the most
obvious example of airfoils, but airfoils are found all throughout industry in applications large and small, terrestrial
and otherwise, anywhere that the flow of liquids and gasses must be controlled or managed. With their ubiquity comes
the need to measure, analyze, and control their geometry.

Typically, the main purpose of an airfoil is to create some specific aerodynamic condition, and its performance is
critical to the overall performance of the system it belongs to. Because aerodynamic analysis is expensive, the 
aerodynamic designers will use control over geometric features as a proxy for control over aerodynamic performance, 
applying limits that are based on expected performance losses from various geometric defects.

Over the years, these measurements have become more sophisticated to act as better proxies for the underlying 
aerodynamic characteristics, and so the tools required to make these measurements have become more sophisticated as 
a result. 

The current (and historical) standard for airfoil analysis is to take measurements on two-dimensional cross-sections of 
the airfoil in planes that roughly have no out-of-plane component of fluid flow, and then analyze the geometry of the
2D section perimeter.  

> [!NOTE] 
> Some organizations are starting to develop methodologies and tools for 3D airfoil analysis, but the 2D section-based 
> approach is still the most commonly used and broadly understood.

### Edges, Sides, and Curvature

Airfoil shapes vary wildly across applications, but in general they have some common geometric properties that can be
used to talk about them in a general way.

- With the notable exception of a spinning cylinder, airfoils are always longer in one direction than the other. They 
  range from elliptical or teardrop shaped to thin and flat.
- Airfoils can be thought of as having two distinct "edges", which form the ends of the longer direction, and two
  "sides" which connect the edges.
- Airfoils generally are thinner near the edges and thicker near the middle. It is very rare for an airfoil to have 
  more than one local thickness maximum, and in most cases the thickness distribution grows, reaches a single 
  maximum that is both local and global, and then tapers back down until it reaches the other edge.
- Airfoils may be curved along their long direction.  This curvature is called camber, and it usually (but not always)
  is in a single direction. When an airfoil is curved, one side of the airfoil becomes more convex than the other.

### The Forward/Aft and Upper/Lower Directions

All airfoils have two conceptual directions which are the basis for a formal definition of its geometry.  

The first is the forward vs aft direction, which describes the general orientation of fluid flow across the airfoil. 
Fluid is coming from the forward direction of the airfoil and flowing towards the aft direction. The edge of the airfoil 
in the forward direction is usually called the leading edge, and the edge in the aft direction called the trailing edge.

The second conceptual direction describes which side of the airfoil is in the direction of the lift force. There are 
a number of different conventions for this naming, varying by organization and application:

- Upper, Suction, or Convex: the side of the airfoil that is in the direction of the lift force. On cambered airfoils 
  this is the side that is geometrically _more_ convex, hence that name. It also is the side with the lower pressure in 
  the fluid flow, leading it to be called the suction side by some organizations.  On traditional aircraft wings which 
  use lift to counteract gravity, this is side of the wing facing away from the earth, hence the name upper.

- Lower, Pressure, or Concave: the side of the airfoil that is opposite the lift force. It is the side which generates 
  the higher pressure in the fluid flow, pushing on the airfoil to generate the lift force, and leading to the name 
  "pressure side". In traditional aircraft wings which use lift to counteract gravity, this is the side of the wing 
  facing the earth, hence the name lower.  Finally, on cambered airfoils it is common to call this the concave side,
  but although this side will be less convex than the other side, it may not _actually_ be concave, depending on the 
  ratio of thickness to camber.

> [!NOTE]
> There is no universal convention for the upper/lower, suction/pressure, or convex/concave naming. Convex/concave 
> naming doesn't make sense for symmetric airfoils or those with low curvature.  The upper/lower distinction only makes
> physical sense for aircraft wings generating lift, and can be confusing for airfoils designed to produce downforce, or
> for airfoils in rotating or vertical applications. The suction/pressure naming convention is the most universally
> applicable, since all airfoils designed to produce lift have one side that sustains higher pressure than the other, 
> but it seems to be the least commonly used. 

These directions are easy to define in unambiguous parts of the airfoil, such as near the middle, but the exact 
transition between upper and lower sides and the position of the _most_ forward and _most_ aft points requires the
leading and trailing edges to be formally defined.

### The Mean Camber Line (MCL)

The mean camber line (MCL) is a curve that runs through the middle of the airfoil, equidistant from the upper and lower 
surfaces.  It is mostly a subset of the mathematical [medial axis](https://en.wikipedia.org/wiki/Medial_axis) of the 
airfoil, and it provides a means of formally defining the leading and trailing edges, and thus the partitioning of the 
upper and lower sides.

Every point on the mean camber line is equally distant from the closest point on the upper and lower surfaces, although
some additional information is usually required at the leading and trailing edges.

- Mathematically, the mean camber line should start at the point of maximum curvature near the leading edge, and end at 
  the point of maximum curvature near the trailing edge.  For airfoils that have a maximum curvature point at both 
  edges, such as a sharp edge or a minimum radius edge, the mean camber line definition is unambiguous.

- For airfoils that have a leading and/or trailing edge that has a finite arc length of a fixed radius (such as a full 
  round), the mean camber line is only unambiguous up to the center point of that arc, after which it could turn
  towards any point on the arc and still satisfy the equidistant condition.  In this case, the convention is that the 
  direction of the mean camber line at the arc center is preserved and the mean camber line is extended in a straight 
  line out to the arc perimeter.

- Finally, for manufacturing reasons, it's common for some airfoils have square trailing edges. Sometimes these include 
  small corner fillets.  In this case the mean camber line is mathematically ambiguous near the edge, but upper and 
  lower sides are clearly defined by the square edge. The usual convention is that the mean camber line ends at the
  midpoint of the square edge.

There are a variety of different methods to find the mean camber line, but historically the most common is to find the
inscribed circles of the airfoil, which are the circles that are tangent to both the upper and lower surfaces.  The 
inscribed circles are unambiguous until they get close to the leading and trailing edges, and the centers of these 
circles will lie on the mean camber line.  Once the unambiguous inscribed circles are found, the leading and trailing
edges can be identified using methods specific to the geometry of the edges.

### Chord 

The chord of an airfoil is usually defined as a straight line from the leading edge point to the trailing edge point.
The chord is used often used to measure the angle of the airfoil, to measure the length of the airfoil, and under
certain conventions (especially older ones) to measure the location of features along the airfoil.

There is a legacy method for measuring the chord length that you may still encounter, sometimes called the "caliper 
method".  It was used on highly cambered airfoils in the era of physical layout inspections.  In this method, the 
concave side of the airfoil is rested against a straightedge, making contact somewhere near the leading and trailing
edges. The tips of a pair of caliper jaws was placed in contact with the straightedge...forcing the jaws to be 
perpendicular to the straightedge...and then closed until contacting the airfoil.

### Airfoil Locations

Most modern measurements taken on airfoil require some mechanism for specifying and locating points on the surface of 
the airfoil cross-section.  These points are then used for things such as measuring thickness, defining GD&T zones, 
defining alignment points, checking positions and angles, etc.

There are a number of different methods used to specify locations on the airfoil, which have varied by organization,
industry sector, and decade.  Some of the commonly used methods include:

- _Position along the mean camber line._ In this method, a point is specified by advancing along the mean camber line
  by a given distance, starting at either the leading or trailing edge point. Then a line of intersection is made 
  orthogonal to the MCL at that position, and the intersection of that line with the side of interest (upper or lower)
  defines the point.  In modern contexts this is the most commonly used method, and performs equally well across the
  entire length of the airfoil, without compromises on highly curved airfoil shapes.

- _Position along the chord line._ This is an older convention that may still be encountered. It is nearly identical to
  the mean camber line position method, but uses the chord instead.  Because the chord is a straight line, this method
  was historically easier to perform.  The point is specified as a distance along the chord line, starting at either
  the leading or trailing edge, and then a line is created orthogonal to the chord line, defining the point of interest
  as the intersection between the orthogonal line and the side of interest (upper or lower).  

- _Offset from leading or trailing edge._ This method has often been used for points near the leading or trailing edge. 
  It involves starting at the leading or trailing edge point, offsetting a specified distance in a specified direction, 
  and then finding the intersection between a line perpendicular to the offset direction and the side of interest 
  (upper or lower).  The direction of the offset may be specified in a number of different ways, most commonly as a 
  nominal vector from design geometry or as the tangent direction of the mean camber line at the edge.  In legacy
  contexts, this method could circumvent the need to find the mean camber line.

- _Intersection with radius_. When used, this method is usually reserved for locating points very close to the leading
  or trailing edges. In this method, a circle of a specified radius is created, centered at the leading or trailing edge
  point, and the point of interest is defined as the intersection between that circle and the side of interest (upper 
  or lower).  This method has an advantage over the previous methods in that it is better at resolving points near the
  edge of the airfoil where the angle between the mean camber line and the surface tangent is large.  It has a
  disadvantage farther from the edges, especially on highly curved airfoils, where the upper and lower points specified
  by the same radius circle can be very different distances along the airfoil.

### Airfoil Thickness

The thickness of an airfoil is the distance between the upper and lower surfaces.  To measure thickness at any discrete
location, a pair of corresponding points on the upper and lower surfaces needs to be identified.  These points are 
identified using one of the methods described in the previous section, and the thickness is measured as the Euclidean
distance between them.

Commonly, airfoil thickness is controlled directly at some combination of the leading edge, the trailing edge, and 
the position of maximum thickness.  Typically, leading and trailing edge thicknesses are each measured at a single pair 
of gage points defined by one of the methods discussed above, while the maximum thickness is found by checking all 
possible pairs of points that would be associated through either the mean camber line or the chord line methods.

The maximum thickness mean camber line method is equivalent to finding the largest inscribed circle of the airfoil,
which in turn can be done without needing to find the entire mean camber line for most airfoil shapes. As a result, 
the camber line/inscribed circle has been used to define maximum thickness even in legacy contexts, and it is unusual
to encounter maximum thickness defined perpendicular to the chord line except on very straight airfoils.


## Pre-Measurement Geometry Establishment

Before any airfoil-specific measurements can be taken on a section, the core airfoil geometry needs to be established. 
The mean camber line needs to be identified, the leading and/or trailing edge points need to be located, and the upper
and lower sides of the airfoil need to be partitioned.

In the `engeom` library there are two main ways of going about this, depending on whether you have nominal reference 
geometry to serve as a starting point.

- If you _don't_ have nominal reference geometry it is either because you don't have a CAD model, or because the section
  you're working with _is_ the nominal reference geometry. In this case, the airfoil geometry is established using
  purely geometric algorithms.
- If you _do_ have nominal reference geometry, the `engeom` library has tools to take a hybrid approach which uses the
  nominal geometry to establish were features should be, and then performs more robust local searching to capture the
  actual features.

### Purely Geometric Analysis

1. Find the unambiguous inscribed circles
2. Order the inscribed circles from leading edge to trailing edge
3. Orient the inscribed circles so that their first contact point is on the lower surface and the second contact point 
   is on the upper surface
4. Identify the leading edge geometry
5. Identify the trailing edge geometry
6. Construct the full mean camber line through the inscribed circle centers and the leading/trailing edge geometry
7. Partition the upper and lower surfaces using the mean camber line and the leading/trailing

### Hybrid Analysis with Nominal Reference Geometry


