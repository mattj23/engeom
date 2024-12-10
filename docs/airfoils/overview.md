# Airfoil Overview

One of the first feature targets of the `engeom` library was a set of tools to perform analysis of airfoil cross-sections.

Measurement of airfoils is an important task in the aerospace propulsion industry, and the ecosystem is dominated by closed-source/proprietary software which is either very expensive, closely guarded, or both.  These tools also tend to have categories of geometries where they perform poorly or unexpectedly, with limited ability to adjust the methodology or tune parameters to correct the results.

The `engeom` library contains a set of airfoil cross-section analysis tools which I developed over roughly fifteen years to perform analysis on highly curved turbine airfoil data produced by industrial structured light scanning systems.  These tools were originally developed in Python and C++, and have been ported to Rust for `engeom`.  A major goal of the re-implementation in Rust is to provide a set of tools which are open-source, transparent, and have an open architecture for modification and extension by users.  

The airfoil tools in `engeom` consist of a general architecture for extracting the mean camber line, thickness distribution, and edge geometry/conditions of nominal airfoil cross-sections, a set of primitive types and traits that form a skeleton for implementing replacement/alternative algorithms, and a set of helper tools, utilities, and partial algorithms to prevent users from having to re-implement common tasks while building alternate or custom airfoil analysis tools.


