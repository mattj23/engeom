# Project Conventions Overview

I'm writing this document in July 2026 (current version is 0.3.5) after a few years of expanding the library and using it across a large number of applications. I feel that the library is getting to the point where it would benefit from some consistency, and that I've seen enough use cases to be in a position to understand at least some of the consequences of the choices made to achieve it.

This document is a place to start keeping track of conventions used in the project.  If, in the long term, it needs to be something more formal I'll change it then. For now, I'll start accumulating these design decisions here and see how it goes.

# Definitions

| Term                | Definition                                                                                                                                                                              |
|---------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Numeric Primitive   | A fundamental numeric entity used to compose more complicated concepts and entities (ex: points, vectors, scalars, etc.)                                                                |
| Geometric Primitive | A small, foundational entity representing a simple geometric construct. Geometric primitives are composed of numeric primitives and can be composed to solve more complicated problems. |

# Naming Conventions

One of the things I struggled the most with while writing and using the library was coming up with consistency in how entities and functions are named. This section is far from exhaustive, but it captures the handful of decisions I've made so far.

## User-Facing Entity Struct Names

For now, these are the guidelines for user-facing structs that represent some conceptual entity (geometric primitives, for example):

- Shorter names are better, but they should accurately represent the concept. An example of where I _failed_ at this was the `Curve2`/`Curve3` entities, which should probably have names that reference that they are polylines, not curves.
- If it is possible that a similar concept exists in a different number of dimensions, use the suffix `2` or `3` in the entity name, _even if `engeom` does not implement a version for the other dimension_.  Only entities that are conceptually identical across dimensions _or_ are a generic implementation from which the 2D/3D versions are derived should be named without a dimension number suffix.

> [!NOTE]
> **Known exception: the `Dof` family.** `Dof3` (in `geom2::align2`) and `Dof6` (in `geom3::align3`) are numbered by their degree-of-freedom count, not by dimension. This collides visually with the suffix rule above: `Dof3` lives in the 2D module but its `3` refers to `tx`/`ty`/`rz`. The count is the more useful thing to convey for these types, so the exception stands deliberately. Note to future self: don't get clever and "fix" `Dof3` to `Dof2`.

## General Conventions for Function Names

These conventions apply to the naming of all functions, including trait and struct implementations.

### Checked/Unchecked Sibling Methods: `try_` and `_unchecked`

Because `Result<Self, E>` and `Option<Self>` are `#[must_use]` and already visible in the signature, there doesn't feel like much benefit in universally prefixing methods that return these types with the word `try_`. For any operation where it's obvious that a failure is a normal possible result (such as an intersection check, or a creation method that's impossible with degenerate input) the function should just return a `Result` or `Option` and leave `try_` out of the function name.

However, if there is sibling method with the same name that doesn't return a fallible type, add the prefix `try_` to the name to distinguish the one returning `Result` or `Option`.

Importantly, if a method can either panic or continue in a state where invalid/unexpected behavior will result, _especially if that behavior is non-obvious and/or deferred_, the method should be suffixed with `_unchecked` to indicate that the caller is responsible for the state of the input.  For example, if a struct depends on being initialized with a sorted `Vec<f64>`, you may want to have a construction method that skips the expensive sort in case the caller already knows that the list is sorted.  This should be named something like `::new_unchecked(a: Vec<f64>)` to indicate that there's something different about this function.

### Verb prefixes in function names

1. When there is an obviously correct, specific verb for the operation the function is doing, use that as the prefix for the function.  Examples include: `find_`, `fit_`, `solve_`, `measure_`, `estimate_`, `extract_`, and so on.

2. _If_ there is no obviously correct verb like #1, but the function is performing a non-trivial computation that is building an output, prefer the prefix `compute_`. Examples I have been guilty of using that should be replaced with `compute_` are: `build_`, `make_`, `calc_`, and `get_`.

3. As per the Rust API guidelines **do not** use the prefix `get_` unless there's one obvious thing to get, and any function named `get_` should only be a pass-through to a reference or value, it shouldn't be hiding computation inside it.

## Struct Function (Method) Names

These are conventions for methods that create, derive, or transform struct instances.  The goal is to have a naming scheme that's clear, consistent, reasonably familiar to users of Rust, and not too verbose.  There were a number of cases where the tradeoffs weren't easy to choose between, so I've written down the reasoning behind some of the decisions so that they can be revisited in the future.

### Simple constructors: `new(...)`

Use `new` when the arguments map roughly 1:1 to the struct's internal fields, e.g. `Circle2::new(x, y, r)`, `Plane3::new(normal, d)`. This is normal Rust practice.

### Derived/complex constructors: `from_<description>(...)`

Use a `from_` name when a constructor builds `Self` from other non-trivial input, but the methodology itself is too complicated or ambiguous to justify a `From<T>` implementation.

Follow `from_` with what the input is (ex. `Circle2::from_3_points`, `Circle3::from_point_normal`).  

If it needs further disambiguation from another method with the same type of input, put something distinguishing about the method at the end of the name.  *However*, if multiple methods can be used with the same input, consider making the method an argument instead of having a separate method.

It makes sense to reverse the input vs method naming convention in certain cases when the input and method are closely related and the method is more recognizable than the input.  For example, `from_least_squares` or `from_fit` is preferable to `from_points_least_squares`.

### Static shape factories: `create_<shape>(...)`

The mesh and point cloud types have a family of static constructors which tessellate a named primitive shape from its dimensions: `Mesh3::create_box`, `Mesh3::create_sphere`, `MeshData3::create_cylinder`, and so on. These keep the `create_` prefix rather than becoming `from_box`/`from_sphere`, because the family is large, mirrored across several types, and reads better with a verb than with `from_` followed by a shape that isn't really the "input" in the sense the rule above means.

`create_` is **reserved** for this. It always means: an associated function with no receiver, which builds a whole new entity out of parameters. A method taking `&self` must never use it — that was the old `Mesh3::create_from_mask`, which read as a constructor but actually returned a subset of an existing mesh, and is now `extract_subset_faces`.

### Methods returning a modified copy: bare past-participle, no prefix

A method that takes `&self` and returns a new, modified `Self` (as opposed to mutating in place) gets a bare past-participle name: `rotated`, `reversed`, `normalized`, `transformed_by`.

Where a verb doesn't have a natural participle, try to find a different name that reads correctly both grammatically and semantically.

This rule is aimed at *operations* on an entity, such as geometric transformations. It does **not** apply to builder-style field overrides on configuration/parameter types, where the standard Rust `with_<field>` idiom is preferred: `AlignParams2::with_tx(1.0)`, `AlignSurfMatch2::with_sigma(0.05)`. Forcing those into participle form produces names that read poorly and diverge from what Rust users expect.

Tentatively, on heavy objects which will primarily use modify-in-place operations, I'm considering using a `_copy` suffix to make expensive operations distinguishable with a quick glance.

### Methods that modify in place

This is still a work in progress, but as of now (July 2026), my tentative plan is to try the following:

- Methods that modify a heavy object in place will be suffixed with `_in_place` to make clear what they're doing.  In Rust this is mostly for readability, because the compiler makes it hard to accidentally modify something you didn't expect to.  However, one of my goals is to keep the names as close as possible in the Python bindings, and the suffix will be distinguishing in that language.

- For now, we'll avoid in-place mutation on lightweight objects like geometric primitives.  If in the future there's a compelling reason to re-introduce them, we can revisit it then.

> [!NOTE]
> While it isn't decided yet, it may make sense to enforce that a struct _either_ has no modify-in-place functions, _or_ it has `_in_place` and `_copy` to distinguish between them.

### Conversions to a different type: `to_` / `into_` / `as_` 

We're going to follow the standard Rust API guidelines:

- `to_`: expensive conversions of borrowed to borrowed, borrowed to owned (via clone), or owned to owned.
- `into_`: a consuming conversion to a different type
- `as_`: a cheap reference-to-reference reinterpretation 

### Conforming pass-through on foreign types

Several types that are part of engeom's public surface aren't ours: `Iso2`, `Iso3`, `Point2`/`3`, `Vector2`/`3` and friends are aliases over `nalgebra`/`parry` types, so their inherent methods are named by someone else and can't be renamed. Where those names violate the conventions above, we do have the option of adding a conforming name to the corresponding extension trait (`IsoExtensions2`, `IsoExtensions3`, ...) that just forwards to the foreign method.

This is worth doing when *all* the following hold:

- The foreign name actually violates a rule here, rather than merely being unfamiliar.
- The operation is part of the everyday surface that users will reach for, not an obscure corner.
- The extension trait already owns that part of the namespace, so the addition completes an existing family instead of starting a redundant parallel one.

The worked example is `IsoExtensions2/3::from_translation` and `from_rotation`, added July 2026. nalgebra spells these `Iso3::translation(x, y, z)` and `Iso3::rotation(axisangle)`, which read as accessors rather than the constructors they are. `IsoExtensions3` already held twelve `from_*` constructors (`from_array`, six `from_basis_*`, `from_rot_axis`, `from_rx`/`ry`/`rz`, `from_z_arbitrary_xy`), so the two most basic constructors were a hole in that family: autocompleting `Iso3::from_` suggested no way to build a plain translation.  Additionally, nalgebra uses the name `translation` for both the constructor *and* the isometry's public field, so `iso.translation.x` and `Iso3::translation(..)` sit side by side meaning different things; `from_translation` can't be confused with the field. 

# Implementation Conventions

## Geometric Primitives

The idea of a geometric primitive is to be a lightweight entity that maps to a fundamental geometric concept.  They should be composed of numeric primitives and, in some cases, other geometric primitives.

- Internally, they should be implemented with the minimum viable representation; for instance, a plane should be represented by a normal and an origin offset (`f64`x4) instead of a point and a normal (`f64`x6) even though the latter is valid.  

- Fields on geometric primitive structs should be publicly accessible.  Accessor methods for convenience are OK for ergonomic reasons.

- Geometric primitive structs must not contain/cache any values that are derived from its own fields.  If the instance owner mutates the value of one of the public fields, it must not be allowed to leave the struct in an incoherent state.

- If derived values are required, such as an `Aabb2`/`3` for some primitives, the struct should provide a method to compute the derived value, and ideally the name and/or documentation should indicate that the value is computed, especially if the computation is non-trivial and the user may want to cache the result.

# Python Binding Conventions

The Python bindings in `py-engeom` have historically tried to mirror the Rust API directly, representing Rust enums and option structs with custom Python classes and enums. After a few years of use I've found this cumbersome, both for the ergonomics of writing Python and for discoverability, and because Python classes don't map cleanly onto Rust enums (especially the data-carrying ones). Rather than mirror the Rust type one-to-one, the convention now is to choose the Python representation based on the _shape_ of the argument.

## Choosing how an argument crosses into Python

| When a Rust argument enum/struct is...                                                         | Python technique                        | Example                                                   |
|------------------------------------------------------------------------------------------------|-----------------------------------------|-----------------------------------------------------------|
| an enum whose variants carry **no payload**, or **one common payload** of the same meaning     | `Literal[...]` string argument          | `AngleDir`, `SelectOp`, `VecDot`, `DistMode`, `Smoothing` |
| an enum where each variant is **uniquely identifiable by which parameter(s) are supplied**     | implicit-by-presence keyword arguments  | `Resample`, `AlignOrigin`, `Lptf3Load`                    |
| an enum with **differing or colliding per-variant payloads** that presence cannot disambiguate | classmethod variant constructors        | `OrientFwdAft`, `OrientUpperLower`                        |
| a **product-type** config/params bundle of independent fields, not a sum type                  | plain keyword constructor with defaults | `Dof6`, `GAPParams`, `AlignParams3`                       |

### No-payload or common-payload enums: `Literal` discriminator plus value

Enums whose variants carry no data are just label selectors. Expose them as lowercase string tokens typed with `Literal[...]` in the stub, e.g. `AngleDir` (`Cw`/`Ccw`) becomes `dir: Literal["cw", "ccw"]`. On the Rust side accept a `&str` and parse it, returning a `PyValueError` whose message lists the valid tokens. Prefer `Literal[...]` over a bare `str` so that IDE autocomplete and static type-checking still work.

When every variant carries the same single payload with the same meaning, the variant name is really just a method label and the payload is one shared argument. Expose the label as a `Literal` string and the payload as an ordinary argument, e.g. `Smoothing` (`Gaussian(f64)`/`Quadratic(f64)`/`Cubic(f64)`) becomes `wrapped_rust_function(method="gaussian", width=...)`.

### Presence-distinguishable enums: implicit-by-presence keyword arguments

When each variant can be identified by _which_ parameters are supplied (distinct argument names, distinct types, or variants that nest as supersets of each other), flatten the variants into optional keyword arguments and let their presence select the variant. Validate at runtime and raise a `ValueError` that names the mutually-exclusive set ("supply exactly one of ...") or the all-or-nothing group ("if any of ... is given, all must be given"). For example `Resample` becomes `count=` / `spacing=` / `max_spacing=` with exactly one required, and `Lptf3Load` collapses into a `take_every=` argument plus an all-or-nothing smoothing parameter group. Where it materially helps a caller, express the valid combinations with `@overload` in the stub.

### Colliding-payload enums: classmethod variant constructors

When variants carry differing payloads, or payloads of the same shape but different meaning, that presence cannot disambiguate, keep a single wrapper class with one named constructor per variant. In PyO3 this is a `#[pyclass]` struct holding `inner: <RustEnum>`, with one `#[staticmethod]` per variant and a single `From<Wrapper> for <RustEnum>` that just unwraps `inner`. For example `OrientFwdAft.airflow(x, y)` and `OrientFwdAft.tmax_fwd()`. Avoid the PyO3 complex-enum representation with nested variant classes, and avoid the `Variant {}` empty-struct trick used to make unit variants callable from Python.

```python
wrapped_rust_function(MyVariant.method1(x=1.0, y=2.0, z=3.0))
wrapped_rust_function(MyVariant.method2())
wrapped_rust_function(MyVariant.method3(true, n=1000))
```

### Product-type config bundles: plain keyword constructor

Structs that bundle independent fields (a product type, not a sum type) get an ordinary keyword `__init__` with sensible defaults, rather than being hidden behind factory staticmethods. `Dof6` is the reference example: `Dof6(tx=True, ty=True, ..., rz=True)`.

> [!IMPORTANT]
> The `.pyi` need to be kept up-to-date with the binding signatures. I have a simple automated name tester in `py-engeom/python/tests/test_stub_drift.py`, which we'll see if it proves to be useful in the long term.

# Plotting Helper Conventions

The plotting helpers in `py-engeom/python/engeom/plot/` are pure Python wrappers that draw `engeom` entities onto some third-party plotting object. They aren't bindings, so nothing above about crossing the Rust boundary applies to them, but they *are* public API and had drifted badly.

These conventions came out of the matplotlib overhaul (August 2026). They're written down because the PyVista helper still has to be done and most of the decisions should apply.

## API Surface shape

- `draw_` prefixes anything that adds an artist.** A method without the prefix configures or queries instead. This splits the API surface so any editor with autocomplete gives you a makeshift table of contents for drawing functions, which seem to be the ones that you spend the most time jumping around between during ordinary use.  Note: when doing the `pyvista` helper, I may consider using `add_` to match the existing convention instead of `draw_`...at this point I'm still undecided.

- Names are singular and take varargs.** `draw_circle(*circles)`, not `draw_circles`. One rule with no per-method judgement, chosen because `draw_circle(c)` is 95% of calls and has to read naturally. Drawing nothing returns an empty list rather than raising, so a computed and possibly empty collection needs no special case at the call site.

> [!NOTE]
> **Known exception: `draw_normals`.** It's plural because the varargs are *sources* (curves or boundaries) and the entities drawn are the normals sampled along them; one call on one source draws `count` arrows. The name describes what lands on the plot, which is the thing the user is looking for in autocomplete. `draw_point` is a different sort of exception: it's singular and varargs per the rule, but additionally accepts a single `(n, 2)` array, since that's the form the rest of the library hands point sets back in. The two are told apart by `numpy.ndim`, which reports 0 for an `engeom` primitive, 1 for a loose coordinate, and 2 for an array.

- Every `draw_*` returns its artists, annotated.** A varargs method returns one artist per entity in the order given; a composite returns every artist it added, in draw order. This is what lets a caller drop down to the host library for anything the helper doesn't expose, and it's the reason `helper.ax` alone isn't sufficient, since finding an artist again after the fact is worse than being handed it. Two deliberate departures: `draw_point` returns a single artist because all its points go into one, and `draw_airfoil` returns a `dict` keyed by element name because with six named parts, indices into a flat list would shift whenever an element was suppressed.

- A more specific verb wins where the host has a separate entry point for it.** `fill_curve` stays `fill_curve` rather than becoming `draw_curve(fill=True)`, because filling goes through `Axes.fill` and stroking through `Axes.plot`, and the two accept disjoint keyword arguments, so a flag switching between them would be a discoverability trap. Where the host takes a real flag on one artist, it stays a flag: `draw_circle(..., fill=False)` maps onto `Circle(fill=)`.  This isn't ideal, but if I find a better way I may change it.

- No `get_`,** same as the Rust rule. `get_3d_viewport` became `viewport`.

## API Structure

- Ideally, wrap instead of subclass, and keep the host object public.** `helper.ax` and `helper.pv` are the documented escape hatch. 

- The helper is named `<HostType>Helper`,** after the object it wraps rather than the backend: `AxesHelper`, `PlotterHelper`. The backend prefix was redundant with the host-type noun once the module path already named the backend. The `Helper` suffix stays because `Axes` and `Plotter` are names users have in scope in the same file, and shadowing them would be hostile.

- One public module per backend** under `engeom.plot`, so the import statement *is* the dependency declaration. Each backend imports its dependency at the top of its `__init__` behind a single guard that re-raises with a message naming the package and the install command. `engeom.plot` itself stays backend-neutral and must not pull either one in, and there's a test asserting that via `sys.modules` in a subprocess. Naming a submodule `matplotlib` is safe: Python 3 absolute imports mean `import matplotlib` inside it resolves to the real package.

## Arguments and styling

- Name the common styling arguments; keep `**kwargs` open.** Each `draw_*` spells out the handful of arguments that get used constantly (`color`, `linewidth`, `linestyle`, `alpha`, `label`, and the patch/text equivalents) as real keyword parameters so editors complete them, and still forwards an open untyped `**kwargs` for everything else the host accepts. I considered PEP 692 `Unpack[TypedDict]` for this and rejected it: it *closes* the set to type checkers, so every valid-but-unlisted matplotlib argument would be flagged as an error.

> [!IMPORTANT]
> Naming a styling argument means an unsupplied one arrives as `None`, and forwarding that is **not** the same as omitting it. `Axes.plot(color=None)` overrides matplotlib's color cycle, so two `draw_curve` calls would come out the same color. Route named styling arguments through `_style.merge_style`, which drops the `None`s. There is no collision risk with `**kwargs`: a named argument binds to its parameter and never reaches the dict, so Python rejects a duplicate before the merge runs.  The downside is that you can't deliberately insert a `None`, which may have some consequences for certain methods.  I'm not sure what the best thing to do is here.

- Where possible, composite draw methods get one argument per element, doing double duty as toggle and style.** `False` suppresses the element, `True` or `None` accepts the defaults, and a dict of keyword arguments is merged **over** the defaults rather than replacing them, because restyling one property must not silently discard the rest of an element's designed appearance. `_style.element_style` resolves this; note that it tests `value is False` rather than falsiness, so that an empty dict means "defaults" instead of "suppress". `draw_airfoil` is the worked example. The alternative, a separate `camber=True` plus `camber_kwargs={...}` per element, doubles the signature for no gain.

## Documentation and tests

Every public member carries a complete Sphinx docstring: `:param:` for every argument including `kwargs`, `:return:`, and `:raises:` where it can raise. These helpers have no `.pyi`, so the docstring and the annotations *are* the reference; `mkdocs.yml` sets `docstring_style: sphinx`.

The examples in `py-engeom/examples/` are the real acceptance test for a plotting API, and are worth running and *looking at* after any change here.
