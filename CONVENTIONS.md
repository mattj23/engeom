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

## Coherence with Foreign Library Entity Names

This example came up while working on mesh decimation. The `alum` library uses the spelling "decimater" due to its OpenMesh origins, while the spelling "decimator" is more comfortable to Americans and throws off less spell-checkers.

In cases where a library has an undesirable or unfamiliar spelling of an entity, whether or not it may be renamed in interacting parts of `engeom` depends on how public facing the entity is.  An inconsistent spelling is worse than an unfamiliar one, so if `engeom` requires that a user is interacting with that entity directly as a normal part of the associated workflow, it's better to keep the spelling consistent.

In the case of "decimater" vs "decimator", `alum`'s `alum::Decimater` is not directly needed for normal uses of the decimation machinery, and `alum`'s half edge `PolyMeshT` is wrapped by our `HalfEdgeMesh3`.  Because `alum` is mostly contained, and in congruence with wrapping `alum` in the event that it would ever need to be replaced, I chose to use the spelling `Decimator` within `engeom`.

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

## Properties vs. methods

A zero-argument accessor is a `#[getter]` when it names **inherent state or a cheap derived value** of the object: a noun like `length`, `normal`, `point_count`, `a`/`b`, or an `is_`/`has_` predicate. It stays a plain method when it takes arguments, when it can fail, or when it has to build something substantial on each call. If it builds something substantial that *is* inherent state, cache it and make it a property anyway; `Mesh3.points` is the reference example, holding the built numpy array in an `Option<Py<PyArray2<f64>>>` field so that repeated reads cost nothing.

Arity is what separates a property from a method of the same name, and it does so cleanly: `Plane3.normal` and `CubicSpline2.normal(t)` are different things, so they are spelled differently, and no rule is needed to reconcile them.

The same concept must be reached the same way on every type that has it. This drifted badly before it was noticed: `mesh.point_count` was a property while `cloud.point_count()` was a method, `segment.length` was a property while `curve.length()` was a method, and so on across eight names and seventeen accessors, all of it unintentional. `py-engeom/python/tests/test_accessor_consistency.py` now fails on any accessor that is a property on one type and a zero-argument method on another.

Note that the stub drift test cannot catch this. The `.pyi` files were accurate the whole time, faithfully recording a `def` wherever the Rust said method; what disagreed was the runtime surface with itself.

That first check only compares names appearing on two or more types, which is how `Cone3.base_center` stayed a method: nothing else has a `base_center` to disagree with. A second check catches it from the other direction, by requiring every zero-argument method whose name reads as a noun to be either a property or a recorded decision in `DELIBERATE_METHODS`, with the reason written beside it. The decisions made so far:

| Kept a method | Why |
|---|---|
| `Curve2/3.at_front`, `at_back`, `Boundary2.at_start`, `at_end` | the rest of the `at_` family takes arguments, and splitting it would read worse than either extreme |
| `CubicSpline2/3.arc_length` | Gauss-Legendre quadrature, not a stored value |
| `Iso2.flipped`, `Iso3.flipped_around_*`, `Circle3/Plane3.normal_reversed` | bare past-participles returning a modified copy |
| `Iso2/3.inverse` | a derived transform rather than state, and the two dimensions agree |
| `Vector2/3.norm` | `norm()` is the spelling in numpy and in nalgebra, which `Vector3` wraps |
| `IndexMask.any`, `all` | mirror `ndarray.any()`/`.all()`; `count_true` is a property, naming a quantity rather than mirroring numpy |
| `AfGeometry.max_thickness`, `max_thickness_circle` | search the inscribed circle stack |

> [!NOTE]
> One accessor is deliberately left copying on each read: `Align2.residuals` and `Align3.residuals` build a fresh numpy array every time, because caching it would cost those types the `Clone` they need to be extractable from a Python object. It is the smaller cost in any case, since `AlignOutcome.alignment` clones the whole alignment, residuals included, every time it is read.

> [!IMPORTANT]
> The `.pyi` need to be kept up-to-date with the binding signatures. I have a simple automated name tester in `py-engeom/python/tests/test_stub_drift.py`, which we'll see if it proves to be useful in the long term.

# Plotting Helper Conventions

The plotting helpers in `py-engeom/python/engeom/plot/` are pure Python wrappers that draw `engeom` entities onto some third-party plotting object. They aren't bindings, so nothing above about crossing the Rust boundary applies to them, but they *are* public API and had drifted badly.

These conventions came out of the matplotlib overhaul (August 2026), and the PyVista helper was brought onto them afterwards. Where the two backends differ it is because the host libraries differ: argument names follow the host (`linewidth` in matplotlib, `line_width` in PyVista), while anything that isn't the host's word for something is the same in both.

## API surface shape

- **`draw_` prefixes anything that adds an artist.** A method without the prefix configures or queries instead. This splits the API surface so an editor's autocomplete provides a makeshift table of contents for the drawing functions used most often. The PyVista helper also uses `draw_`, rather than the `add_` that would match PyVista's own `add_mesh` and `add_points`, so the convention is consistent across backends and `draw_circle` means the same thing whichever helper is in hand.

- **Names are singular and take varargs.** Use `draw_circle(*circles)`, not `draw_circles`. This avoids per-method judgment while allowing the common `draw_circle(c)` call to read naturally. Drawing nothing returns an empty list rather than raising, so a computed and possibly empty collection needs no special case at the call site.

> [!NOTE]
> **Known exception: `draw_normals`.** It's plural because the varargs are *sources* (curves or boundaries) and the entities drawn are the normals sampled along them; one call on one source draws `count` arrows. The name describes what lands on the plot, which is the thing the user is looking for in autocomplete. `draw_point` is a different sort of exception: it's singular and varargs per the rule, but additionally accepts a single `(n, 2)` array, since that's the form the rest of the library hands point sets back in. The two are told apart by `numpy.ndim`, which reports 0 for an `engeom` primitive, 1 for a loose coordinate, and 2 for an array.

> [!NOTE]
> **Known exceptions on the PyVista side: `draw_mesh` and `draw_point_cloud`.** Both take one entity rather than varargs, because their `scalars`, `highlight` and `use_colors` arguments describe that one entity's per-point or per-face data and could not be shared across several. `draw_point` is singular-with-varargs but returns a single actor, since every point given goes into one, and it also accepts an `(n, 3)` array or a `PointCloud3` for the same reason its matplotlib counterpart accepts an array.

- **Two of PyVista's arguments cannot simply be repeated across a varargs call.** An actor `name` replaces any existing actor of that name, so reusing one would leave only the last entity drawn; it is suffixed with the entity's index when more than one is drawn. A legend `label` would produce one identical entry per entity, so only the first carries it and the group gets a single legend entry. Both are handled in one place, `_per_entity`.

- **Every `draw_*` method returns its artists and annotates the return type.** A varargs method returns one artist per entity in the order given; a composite returns every artist it added, in draw order. This lets callers use the host library directly for anything the helper does not expose. Merely exposing `helper.ax` is insufficient because recovering an artist afterward is less useful than returning it directly. There are two deliberate departures: `draw_point` returns a single artist because all its points go into one, and `draw_airfoil` returns a `dict` keyed by element name because, with six named parts, indices into a flat list would shift whenever an element was suppressed.

- **A more specific verb wins where the host has a separate entry point for it.** `fill_curve` stays `fill_curve` rather than becoming `draw_curve(fill=True)`, because filling goes through `Axes.fill` and stroking through `Axes.plot`. The two accept disjoint keyword arguments, so a flag switching between them would make the accepted options harder to discover. Where the host takes a real flag on one artist, it stays a flag: `draw_circle(..., fill=False)` maps onto `Circle(fill=)`. This is not ideal, but it can change if a better convention emerges.

- **No `get_`,** following the Rust rule. `get_3d_viewport` became `viewport`.

## API structure

- **Prefer wrapping to subclassing, and keep the host object public.** `helper.ax` and `helper.pv` are the documented escape hatches.

- **The PyVista helper is also attached to the plotter as `plotter.engeom`.** Wrapping is the right internal structure, but it leaves two very similar objects in scope at every call site. Mixing `helper.draw_mesh(...)` with `plotter.add_points(...)` makes them appear unrelated, while the attachment reflects that the helper belongs to that plotter. It uses PyVista's component mechanism (0.48+), which constructs the helper lazily on first access and caches it on the plotter. `PlotterHelper(plotter)` still works and provides access on older PyVista versions.

  Registering against `pyvista.BasePlotter` rather than `Plotter` is deliberate, so that the Qt and background plotters from `pyvistaqt` get it too. It is declared as a `pyvista.plotter_components` entry point in `py-engeom/pyproject.toml`, so `plotter.engeom` resolves in a session that never imported `engeom.plot.pyvista`. Attaching to someone else's class has to be undoable, so `register()` and `unregister()` are public and tested.

- **The one subclass, `EngeomPlotter`, exists only so editors can see the accessor.** An attribute attached at runtime is invisible to a type checker reading `pyvista.Plotter`, undermining the discoverability provided by the `draw_` prefix. The subclass declares `engeom: PlotterHelper` and nothing else, using an annotation rather than an assignment so it does not shadow the descriptor that does the work. This is not a license to subclass for behavior.

- **The helper is named `<HostType>Helper`,** after the object it wraps rather than the backend: `AxesHelper`, `PlotterHelper`. A backend prefix would be redundant once the module path names the backend. The `Helper` suffix remains because users commonly have `Axes` and `Plotter` in scope in the same file, and shadowing those names would be hostile.

- **Use one public module per backend** under `engeom.plot`, so the import statement *is* the dependency declaration. Each backend imports its dependency at the top of its `__init__` behind a single guard that re-raises with a message naming the package and installation command. `engeom.plot` itself remains backend-neutral and must not import either dependency; a subprocess test asserts this through `sys.modules`. Naming a submodule `matplotlib` is safe because Python 3 absolute imports ensure that `import matplotlib` inside it resolves to the real package.

## Arguments and styling

- **Name the common styling arguments, but keep `**kwargs` open.** Each `draw_*` method spells out the frequently used arguments (`color`, `linewidth`, `linestyle`, `alpha`, `label`, and the patch/text equivalents) as real keyword parameters so editors can complete them. It also forwards an open, untyped `**kwargs` for everything else the host accepts. PEP 692's `Unpack[TypedDict]` was considered and rejected because it *closes* the set to type checkers, causing every valid but unlisted Matplotlib argument to be flagged as an error.

> [!IMPORTANT]
> Naming a styling argument means an unsupplied one arrives as `None`, and forwarding that is **not** the same as omitting it. `Axes.plot(color=None)` overrides matplotlib's color cycle, so two `draw_curve` calls would come out the same color. Route named styling arguments through `_style.merge_style`, which drops the `None`s. There is no collision risk with `**kwargs`: a named argument binds to its parameter and never reaches the dict, so Python rejects a duplicate before the merge runs.  The downside is that you can't deliberately insert a `None`, which may have some consequences for certain methods.  I'm not sure what the best thing to do is here.

- **Where possible, composite draw methods have one argument per element, serving as both toggle and style.** `False` suppresses the element, `True` or `None` accepts the defaults, and a dictionary of keyword arguments is merged **over** the defaults rather than replacing them. This ensures that restyling one property does not silently discard the rest of an element's designed appearance. `_style.element_style` resolves this; it tests `value is False` rather than falsiness, so an empty dictionary means "defaults" instead of "suppress." `draw_airfoil` is the worked example. The alternative—a separate `camber=True` plus `camber_kwargs={...}` for each element—doubles the signature for no gain.

## Unbounded entities

A `Plane3` has no origin, orientation or size, and a `Line3` has no ends, so neither can be drawn without deciding how much of it to show. The PyVista helper takes an `extent` argument for this: an `Aabb3`, anything carrying a bounding box, an `(n, 3)` array of points, or `None` to mean everything already drawn into the active renderer. The entity is then clipped to that box, so the drawn polygon follows the shape of the region rather than being an arbitrary square floating in it, and a `pad` fraction grows the box so a plane cutting through a part protrudes past it instead of stopping flush with its surface.

`None` being the default is the point of the whole thing: draw the part, then draw the plane through it, with nothing to work out at the call site. An empty scene is refused rather than defaulted, because PyVista reports a two-unit cube for a renderer with nothing in it and sizing an entity against that placeholder would draw something plausible and wrong.

An entity that misses the extent raises, naming its position in the argument list, rather than being skipped. Skipping would break the rule that a varargs method returns one artist per entity in the order given, by shifting every later index.

The Matplotlib `ViewPort3` predates this and requires an explicit `center` and `size` on every call. Its `draw_plane` and `draw_line` should move to the same `extent` argument, minus the scene-bounds default, which has no equivalent on an axes.

## The view convention

Both backends describe a view as an `Iso3` transforming world space into the image plane, where **+X is to the right, +Y is up, and +Z points out of the image towards the viewer**. The third of those is not a free choice: with +X right and +Y up, a rigid transform has to have +Z coming towards the viewer, because the other arrangement is left-handed and therefore not a rotation at all. It is also what VTK and OpenGL use for a camera.

`PlotterHelper.camera_pose` returns a view in this convention and `view_pose` takes one, so a view found by dragging the render window can be replayed later or handed to `AxesHelper.viewport` to draw a line diagram from the same angle. VTK does not keep its camera's up vector square to the view direction, so `camera_pose` orthogonalizes before building the frame; without that, `Iso3` rejects the matrix outright.

`ViewPort3` had this documented backward for a while, saying +Z pointed into the image. Only the sentence was wrong: every drawing it produced, including the hidden-line classification in `draw_mesh_outline`, was already using the convention above. This distinction is worth knowing if a saved view looks mirrored, because flipping the isometry is not the fix.

## Coloring, and the order of specificity

Where an entity can be colored several ways, the arguments are ranked by specificity, and the most specific wins: an explicit `color`, then `scalars`, then a `highlight`, and finally the entity's stored colors through `use_colors`. `scalars` and `highlight` cannot be combined because each sets every element's color; asking for both raises rather than silently choosing one. Giving an explicit `color` precedence over stored colors required a deliberate choice: PyVista will not accept a color and scalars together, and the argument written by the caller is more specific than the colors already stored on the entity.

A `highlight` is drawn as per-element colors on the one actor rather than as a second actor over the selection, so it cannot z-fight with the surface beneath it and brings no scalar bar with it. It takes an `IndexMask`, a boolean array, or indices; a mask knows how many elements it covers, which is what tells a selection of faces from one of points without the caller saying which.

## Interaction

A widget or picker takes a callback and passes it `engeom` entities rather than the host library's raw values. The plane widget reports a `Plane3` instead of a normal-and-origin pair, surface picking reports a `SurfacePoint3` instead of a bare coordinate, and face picking reports an `IndexMask` over the mesh's own faces instead of the renumbered cells of an extracted dataset. The callback can therefore be ordinary `engeom` code with no knowledge of the render window, and its values are the same types already accepted by the Rust side.

Face picking works because `to_polydata` records each face's original index on every dataset it builds, in the `FACE_ID` array. PyVista's filters carry cell data through, whereas VTK's own original-cell-ids are discarded by the extraction a rubber-band selection performs, so that stamp is the only thing relating a selection back to the mesh it came from.

Anything the helper switches on it also switches off: the component's `__plotter_close__` releases the pickers and widgets when the plotter closes. PyVista calls that on every close and a plotter can be closed more than once, so it has to be safe to run again with nothing left to do.

## Documentation and tests

Every public member carries a complete Sphinx docstring: `:param:` for every argument including `kwargs`, `:return:`, and `:raises:` where it can raise. These helpers have no `.pyi`, so the docstring and the annotations *are* the reference; `mkdocs.yml` sets `docstring_style: sphinx`.

The examples in `py-engeom/examples/` are the real acceptance test for a plotting API, and are worth running and *looking at* after any change here.
