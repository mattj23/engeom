"""
Guards that an accessor spelled the same way on two types is reached the same way on both.

The bindings drifted into a state where `mesh.point_count` was a property while
`cloud.point_count()` was a method, `segment.length` was a property while `curve.length()` was a
method, and so on for eight names across seventeen accessors. Nothing caught it: the `.pyi` stubs
were accurate, faithfully recording a `def` wherever the Rust said method, so `test_stub_drift.py`
compares the two and finds them in agreement. What was inconsistent was the runtime surface with
itself.

The rule this enforces is that a zero-argument accessor naming inherent state, a cheap derived
value, or an `is_`/`has_` predicate is a `#[getter]`. An accessor that takes arguments is a method
and is left alone here, which is what allows `CubicSpline2.normal(t)` to coexist with
`Plane3.normal` without a special case: they differ in arity, so they are different things.
"""

from __future__ import annotations

import collections
import importlib
import inspect

import pytest

MODULES = [
    "engeom.airfoil2",
    "engeom.align2",
    "engeom.align3",
    "engeom.common",
    "engeom.geom2",
    "engeom.geom3",
    "engeom.metrology",
    "engeom.raster2",
    "engeom.sensors",
]


def _takes_arguments(value) -> bool:
    """
    Whether a bound method needs anything beyond `self`.

    A signature that cannot be read is treated as taking arguments, so that an accessor this cannot
    inspect is left out of the comparison rather than reported as a false inconsistency.
    """
    try:
        parameters = list(inspect.signature(value).parameters.values())
    except (TypeError, ValueError):
        return True
    return any(p.name != "self" for p in parameters)


def _surface() -> dict[str, list[tuple[str, str]]]:
    """
    Map every public accessor name to the `(class, kind)` pairs exposing it.

    `kind` is `"property"` for a getter and `"method"` for a zero-argument method. Anything taking
    arguments is left out, since a method and a property cannot be spelled the same way if one of
    them needs an argument.
    """
    found: dict[str, list[tuple[str, str]]] = collections.defaultdict(list)
    for module_name in MODULES:
        module = importlib.import_module(module_name)
        for class_name in dir(module):
            cls = getattr(module, class_name)
            if not inspect.isclass(cls):
                continue
            for name in dir(cls):
                if name.startswith("_"):
                    continue
                value = inspect.getattr_static(cls, name, None)
                kind = type(value).__name__
                if kind == "getset_descriptor":
                    found[name].append((class_name, "property"))
                elif kind in ("method_descriptor", "function") and not _takes_arguments(value):
                    found[name].append((class_name, "method"))
    return found


def test_the_modules_under_test_all_import():
    """ Guard against the whole check passing vacuously because a module name went stale. """
    for module_name in MODULES:
        assert importlib.import_module(module_name) is not None


def test_the_surface_is_big_enough_to_be_meaningful():
    """ Guard against the introspection silently finding nothing and the comparison passing. """
    assert len(_surface()) > 100


def test_no_accessor_is_a_property_on_one_type_and_a_method_on_another():
    surface = _surface()
    split = {
        name: entries
        for name, entries in surface.items()
        if len({kind for _, kind in entries}) > 1
    }
    if split:
        report = []
        for name, entries in sorted(split.items()):
            properties = sorted(c for c, k in entries if k == "property")
            methods = sorted(c for c, k in entries if k == "method")
            report.append(f"  {name}: property on {properties}, method on {methods}")
        pytest.fail(
            "These accessors are reached differently depending on the type, so a caller has to "
            "remember which is which:\n" + "\n".join(report) +
            "\n\nAdd `#[getter]` to the odd ones out in py-engeom/src, and `@property` to the "
            "matching entries in the .pyi stubs."
        )


# Prefixes that mark a zero-argument method as doing something rather than naming something. A
# method whose name starts with one of these is a verb, a conversion, or a modified copy, so it is
# a method by the naming conventions and is not a candidate for being a property.
_VERB_PREFIXES = {
    "as", "boundary", "build", "clone", "cloned", "collect", "compute", "derivative", "empty",
    "estimate", "extract", "find", "flip", "garbage", "get", "into", "load", "make", "merged",
    "normalized", "patch", "pop", "query", "reversed", "sample", "save", "separate", "split",
    "to", "try",
}

# Zero-argument methods whose names read as nouns or predicates, and which are deliberately not
# properties anyway. Each is a decision that was made once and should not have to be made again by
# whoever next reads the surface, so the reason is recorded beside it.
DELIBERATE_METHODS = {
    # The `at_` family is a set of queries, and the rest of it takes arguments: `at_length(l)`,
    # `at_fraction(f)`, `at_closest_to_point(p)`. Splitting the family so that two members are
    # properties and three are calls would read worse than either extreme.
    "Boundary2.at_end",
    "Boundary2.at_start",
    "Curve2.at_back",
    "Curve2.at_front",
    "Curve3.at_back",
    "Curve3.at_front",
    # Gauss-Legendre quadrature over the parameter range, so a real computation rather than a
    # value the spline is holding.
    "CubicSpline2.arc_length",
    "CubicSpline3.arc_length",
    # Bare past-participles returning a modified copy, which the naming conventions make methods.
    "Circle3.normal_reversed",
    "Iso2.flipped",
    "Iso3.flipped_around_x",
    "Iso3.flipped_around_y",
    "Iso3.flipped_around_z",
    "Plane3.normal_reversed",
    # A derived transform rather than state the isometry holds, and the two dimensions agree.
    "Iso2.inverse",
    "Iso3.inverse",
    # Stored rather than derived here, but `inverse` has to be reached the same way as on the
    # isometries.
    "PlanarMap.inverse",
    # `norm()` is the spelling in numpy, nalgebra, and essentially every linear algebra library,
    # and `Vector3` wraps nalgebra's. Following the surrounding ecosystem wins here.
    "Vector2.norm",
    "Vector3.norm",
    # Likewise `ndarray.any()` and `ndarray.all()`, which is what a mask is most often used
    # alongside. `count_true` is a property, since it names a quantity rather than mirroring numpy.
    "IndexMask.all",
    "IndexMask.any",
    # Both search the inscribed circle stack rather than reading a stored value.
    "AfGeometry.max_thickness",
    "AfGeometry.max_thickness_circle",
}


def _noun_shaped_methods() -> set[str]:
    """Zero-argument methods whose names do not begin with a verb prefix."""
    found = set()
    for module_name in MODULES:
        module = importlib.import_module(module_name)
        for class_name in dir(module):
            cls = getattr(module, class_name)
            if not inspect.isclass(cls):
                continue
            for name in dir(cls):
                if name.startswith("_") or name.split("_")[0] in _VERB_PREFIXES:
                    continue
                value = inspect.getattr_static(cls, name, None)
                if type(value).__name__ not in ("method_descriptor", "function"):
                    continue
                if not _takes_arguments(value):
                    found.add(f"{class_name}.{name}")
    return found


def test_every_noun_shaped_method_is_a_recorded_decision():
    """
    The cross-type check above cannot see an accessor that is a method on every type that has it,
    which is how `Cone3.base_center` and a few dozen others stayed methods long after the same
    concept had become a property elsewhere. This catches that case from the other direction: any
    zero-argument method that reads as a noun has to be either a property or a recorded decision.
    """
    found = _noun_shaped_methods()

    undecided = found - DELIBERATE_METHODS
    assert not undecided, (
        "These zero-argument methods read as nouns, so they look like properties that were never "
        f"made into one: {sorted(undecided)}. Either add `#[getter]` in py-engeom/src and "
        "`@property` in the matching stub, or add the name to DELIBERATE_METHODS with the reason."
    )

    stale = DELIBERATE_METHODS - found
    assert not stale, (
        f"These are recorded as deliberate methods but no longer exist as such: {sorted(stale)}. "
        "Remove them from DELIBERATE_METHODS."
    )
