"""
    This is an experimental test suite to help try to keep the stubs matched to the actual bindings.

    It iterates through the runtime-discovered modules (skipping airfoil2 for now) and checks what's actually
    in them against the contents of of the stubs. Right now all it does is check the names, so if something
    exists but shouldn't, or doesn't exist but should, it will throw an error.

    The double-under methods are ignored.
"""

from __future__ import annotations

import ast
import inspect
import importlib
from pathlib import Path

import pytest

STUB_DIR = Path(__file__).resolve().parent.parent / "engeom"

# Stub stems to skip for now. The airfoil modules are not yet stabilized and their stubs are
# known to lag; drop these from the exclusion set once they are cleaned up.
EXCLUDE = {"airfoil", "airfoil2"}

# A few modules we always expect to discover, so an empty/broken glob can't make the whole
# parametrized check pass vacuously.
EXPECTED = {"engeom.common", "engeom.geom2", "engeom.geom3", "engeom.engeom"}


def discover_modules():
    """
    Derive the (module, stub filename) pairs from the ``.pyi`` files present next to the package.

    Each ``<name>.pyi`` maps to the importable module ``engeom.<name>``, except ``engeom.pyi``
    which stubs the top-level native module ``engeom.engeom``.
    """
    modules = []
    for stub in sorted(STUB_DIR.glob("*.pyi")):
        stem = stub.stem
        if stem in EXCLUDE:
            continue

        # Manual override for the main module
        module_name = "engeom.engeom" if stem == "engeom" else f"engeom.{stem}"

        modules.append((module_name, stub.name))
    return modules


MODULES = discover_modules()


def test_module_discovery_is_sane():
    """ Guard against the stub glob silently returning nothing or missing a core module. """
    discovered = {m[0] for m in MODULES}
    assert MODULES, f"No .pyi stubs discovered under {STUB_DIR}"
    missing = EXPECTED - discovered
    assert not missing, f"Expected modules not discovered (check EXCLUDE / stub dir): {sorted(missing)}"


@pytest.mark.parametrize("module_name,stub_name", MODULES, ids=[m[0] for m in MODULES])
def test_stub_matches_runtime(module_name, stub_name):
    module = importlib.import_module(module_name)
    tree = ast.parse((STUB_DIR / stub_name).read_text())

    rt_classes, rt_functions = runtime_surface(module)
    stub_classes, stub_functions, stub_bodies, stub_imported = stub_surface(tree)
    declared = stub_classes | stub_functions | stub_imported

    problems = []

    for name in sorted(rt_classes - declared):
        problems.append(f"class `{name}` exists at runtime but is missing from {stub_name}")
    for name in sorted(stub_classes - rt_classes):
        problems.append(f"class `{name}` declared in {stub_name} but does not exist at runtime")
    for name in sorted(rt_functions - declared):
        problems.append(f"function `{name}` exists at runtime but is missing from {stub_name}")
    for name in sorted(stub_functions - rt_functions):
        problems.append(f"function `{name}` declared in {stub_name} but does not exist at runtime")

    for cls_name in sorted(rt_classes & stub_classes):
        cls = getattr(module, cls_name)
        rt_members = runtime_class_members(cls)
        stub_members = stub_class_members(stub_bodies[cls_name])
        for m in sorted(rt_members - stub_members):
            problems.append(f"`{cls_name}.{m}` exists at runtime but is missing from {stub_name}")
        for m in sorted(stub_members - rt_members):
            problems.append(f"`{cls_name}.{m}` declared in {stub_name} but does not exist at runtime")

    assert not problems, "Stub drift detected:\n" + "\n".join(f"  - {p}" for p in problems)


def runtime_surface(module):
    """Return (classes, functions) name-sets for the public runtime surface of a module."""
    classes, functions = set(), set()
    for name in dir(module):
        if not _public(name):
            continue
        obj = getattr(module, name)
        if inspect.isclass(obj):
            classes.add(name)
        elif inspect.isroutine(obj):
            functions.add(name)
    return classes, functions


def runtime_class_members(cls) -> set[str]:
    """Non-dunder names defined directly on a runtime class (methods, properties, nested types)."""
    return {n for n in vars(cls) if _public(n) and not _is_dunder(n)}


def stub_surface(tree: ast.Module):
    """Return (classes, functions, class_bodies, imported) declared at a stub's top level.

    ``imported`` holds names brought in via ``import`` / ``from ... import`` - a stub may
    legitimately re-export a class defined in another module (e.g. ``SplineProjection`` in
    ``geom3``), so those names count as "declared" for the runtime->stub direction.
    """
    classes, functions, bodies, imported = set(), set(), {}, set()
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and _public(node.name):
            classes.add(node.name)
            bodies[node.name] = node
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and _public(node.name):
            functions.add(node.name)
        elif isinstance(node, ast.ImportFrom):
            imported.update(a.asname or a.name for a in node.names)
        elif isinstance(node, ast.Import):
            imported.update((a.asname or a.name).split(".")[0] for a in node.names)
    return classes, functions, bodies, imported


def stub_class_members(node: ast.ClassDef) -> set[str]:
    """Non-dunder member names declared in a stub class body (deduped, so overloads collapse)."""
    members = set()
    for item in node.body:
        if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
            name = item.name
        elif isinstance(item, ast.ClassDef):
            name = item.name
        elif isinstance(item, ast.AnnAssign) and isinstance(item.target, ast.Name):
            name = item.target.id
        elif isinstance(item, ast.Assign) and len(item.targets) == 1 and isinstance(item.targets[0], ast.Name):
            name = item.targets[0].id
        else:
            continue
        if _public(name) and not _is_dunder(name):
            members.add(name)
    return members


def _is_dunder(name: str) -> bool:
    return name.startswith("__") and name.endswith("__")


def _public(name: str) -> bool:
    return not name.startswith("_")

