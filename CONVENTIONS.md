# Project Conventions

I'm writing this document in July 2026 (current version is 0.3.5) after a few years of expanding the library and using it across a large number of applications. I feel that the library is getting to the point where it would benefit from some consistency, and that I've seen enough use cases to be in a position to understand at least some of the consequences of the choices made to achieve it.

This document is a place to start keeping track of conventions used in the project.  If, in the long term, it needs to be something more formal I'll change it then. For now, I'll start accumulating these design decisions here and see how it goes.

## Definitions

| Term                | Definition                                                                                                                                                                              |
|---------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Numeric Primitive   | A fundamental numeric entity used to compose more complicated concepts and entities (ex: points, vectors, scalars, etc.)                                                                |
| Geometric Primitive | A small, foundational entity representing a simple geometric construct. Geometric primitives are composed of numeric primitives and can be composed to solve more complicated problems. |

## Naming Conventions

One of the things I struggled the most with while writing and using the library was coming up with consistency in how entities and functions are named. This section is far from exhaustive, but it captures the handful of decisions I've made so far.

### User-Facing Entity Struct Names

For now, these are the guidelines for user-facing structs that represent some conceptual entity (geometric primitives, for example):

- Shorter names are better, but they should accurately represent the concept. An example of where I _failed_ at this was the `Curve2`/`Curve3` entities, which should probably have names that reference that they are polylines, not curves.
- If it is possible that a similar concept exists in a different number of dimensions, use the suffix `2` or `3` in the entity name, _even if `engeom` does not implement a version for the other dimension_.  Only entities that are conceptually identical across dimensions _or_ are a generic implementation from which the 2D/3D versions are derived should be named without a dimension number suffix.

### General Conventions for Function Names

These conventions apply to the naming of all functions, including trait and struct implementations.

#### Checked/Unchecked Sibling Methods: `try_` and `_unchecked`

Because `Result<Self, E>` and `Option<Self>` are `#[must_use]` and already visible in the signature, there doesn't feel like much benefit in universally prefixing methods that return these types with the word `try_`. For any operation where it's obvious that a failure is a normal possible result (such as an intersection check, or a creation method that's impossible with degenerate input) the function should just return a `Result` or `Option` and leave `try_` out of the function name.

However, if there is sibling method with the same name that doesn't return a fallible type, add the prefix `try_` to the name to distinguish the one returning `Result` or `Option`.

Importantly, if a method can either panic or continue in a state where invalid/unexpected behavior will result, _especially if that behavior is non-obvious and/or deferred_, the method should be suffixed with `_unchecked` to indicate that the caller is responsible for the state of the input.  For example, if a struct depends on being initialized with a sorted `Vec<f64>`, you may want to have a construction method that skips the expensive sort in case the caller already knows that the list is sorted.  This should be named something like `::new_unchecked(a: Vec<f64>)` to indicate that there's something different about this function.

#### Verb prefixes in function names

1. When there is an obviously correct, specific verb for the operation the function is doing, use that as the prefix for the function.  Examples include: `find_`, `fit_`, `solve_`, `measure_`, `estimate_`, `extract_`, and so on.

2. _If_ there is no obviously correct verb like #1, but the function is performing a non-trivial computation that is building an output, prefer the prefix `compute_`. Examples I have been guilty of using that should be replaced with `compute_` are: `build_`, `make_`, `calc_`, and `get_`.

3. As per the Rust API guidelines **do not** use the prefix `get_` unless there's one obvious thing to get, and any function named `get_` should only be a pass-through to a reference or value, it shouldn't be hiding computation inside it.

### Struct Function (Method) Names

These are conventions for methods that create, derive, or transform struct instances.  The goal is to have a naming scheme that's clear, consistent, reasonably familiar to users of Rust, and not too verbose.  There were a number of cases where the tradeoffs weren't easy to choose between, so I've written down the reasoning behind some of the decisions so that they can be revisited in the future.

#### Simple constructors: `new(...)`

Use `new` when the arguments map roughly 1:1 to the struct's internal fields, e.g. `Circle2::new(x, y, r)`, `Plane3::new(normal, d)`. This is normal Rust practice.

#### Derived/complex constructors: `from_<description>(...)`

Use a `from_` name when a constructor builds `Self` from other non-trivial input, but the methodology itself is too complicated or ambiguous to justify a `From<T>` implementation.

Follow `from_` with what the input is (ex. `Circle2::from_3_points`, `Circle3::from_point_normal`).  

If it needs further disambiguation from another method with the same type of input, put something distinguishing about the method at the end of the name.  *However*, if multiple methods can be used with the same input, consider making the method an argument instead of having a separate method.

It makes sense to reverse the input vs method naming convention in certain cases when the input and method are closely related and the method is more recognizable than the input.  For example, `from_least_squares` or `from_fit` is preferable to `from_points_least_squares`.

#### Methods returning a modified copy: bare past-participle, no prefix

A method that takes `&self` and returns a new, modified `Self` (as opposed to mutating in place) gets a bare past-participle name: `rotated`, `reversed`, `normalized`, `transformed_by`.

Where a verb doesn't have a natural participle, try to find a different name that reads correctly both grammatically and semantically.

Tentatively, on heavy objects which will primarily use modify-in-place operations, I'm considering using a `_copy` suffix to make expensive operations distinguishable with a quick glance.

#### Methods that modify in place

This is still a work in progress, but as of now (July 2026), my tentative plan is to try the following:

- Methods that modify a heavy object in place will be suffixed with `_in_place` to make clear what they're doing.  In Rust this is mostly for readability, because the compiler makes it hard to accidentally modify something you didn't expect to.  However, one of my goals is to keep the names as close as possible in the Python bindings, and the suffix will be distinguishing in that language.

- For now, we'll avoid in-place mutation on lightweight objects like geometric primitives.  If in the future there's a compelling reason to re-introduce them, we can revisit it then.

> [!NOTE]
> While it isn't decided yet, it may make sense to enforce that a struct _either_ has no modify-in-place functions, _or_ it has `_in_place` and `_copy` to distinguish between them.

#### Conversions to a different type: `to_` / `into_` / `as_` 

We're going to follow the standard Rust API guidelines:

- `to_`: expensive conversions of borrowed to borrowed, borrowed to owned (via clone), or owned to owned.
- `into_`: a consuming conversion to a different type
- `as_`: a cheap reference-to-reference reinterpretation 

## Implementation Conventions

### Geometric Primitives

The idea of a geometric primitive is to be a lightweight entity that maps to a fundamental geometric concept.  They should be composed of numeric primitives and, in some cases, other geometric primitives.

- Internally, they should be implemented with the minimum viable representation; for instance, a plane should be represented by a normal and an origin offset (`f64`x4) instead of a point and a normal (`f64`x6) even though the latter is valid.  

- Fields on geometric primitive structs should be publicly accessible.  Accessor methods for convenience are OK for ergonomic reasons.

- Geometric primitive structs must not contain/cache any values that are derived from its own fields.  If the instance owner mutates the value of one of the public fields, it must not be allowed to leave the struct in an incoherent state.

- If derived values are required, such as an `Aabb2`/`3` for some primitives, the struct should provide a method to compute the derived value, and ideally the name and/or documentation should indicate that the value is computed, especially if the computation is non-trivial and the user may want to cache the result.
