# Discrete Domain

Conceptually, a domain is the set of values which a variable may take the value of.  For example, in a simple function $y=f(x)$ in the $XY$ plane, the domain $x$ is something like the real number line $\mathbb{R}$.  A *discrete* domain is like a general domain, except that it is a finite set.

There is a `DiscreteDomain` struct in `engeom::common` which is used to represent a discrete domain of real numbers as a set of finite `f64` values kept in a sorted order.  This is basically a wrapper around a `Vec<f64>` which is guaranteed to contain only finite values and only in ascending order, enabling the use of algorithms which rely on these properties to be used.

A `DiscreteDomain` will dereference automatically as a `[f64]`, allowing it to be used in place of a `Vec` or array of `f64` values.


