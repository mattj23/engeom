"""
Helpers for drawing `engeom` entities into a PyVista `Plotter` object.

Importing this module attaches the helper to every PyVista plotter under the name `engeom`, so that
a plotter and the helper that draws onto it are one object rather than two sitting side by side in
the same scope:

```python
import pyvista
plotter = pyvista.Plotter()
plotter.engeom.draw_mesh(mesh, color="white")
plotter.add_axes()
plotter.show()
```

Every method that adds an actor is prefixed `draw_`, so `plotter.engeom.draw_<tab>` lists
everything that can be drawn. `PlotterHelper(plotter)` still builds a helper directly, for code
that would rather hold one, and `helper.pv` is still the plotter it wraps.

The accessor uses PyVista's own plugin mechanism, which arrived in PyVista 0.48. On an older
PyVista nothing is attached and `PlotterHelper(plotter)` is the way in. Because the attachment is
declared as a `pyvista.plotter_components` entry point, `plotter.engeom` also resolves in a session
that never imported `engeom.plot.pyvista`; PyVista imports it the first time the attribute is
looked up.

To detach it, call `pyvista.unregister_plotter_component("engeom")` or the `unregister` below, and
`register` puts it back. Nothing else about the helper changes either way.

## Reaching past the helper

The helper covers drawing `engeom` entities and little else. Everything a `Plotter` can do is
still on the plotter, and a few of those are worth knowing about because they matter for the kinds
of scene this draws:

* `enable_depth_peeling` sorts overlapping transparent surfaces properly. Worth switching on
  whenever more than one translucent thing is in the scene, which planes and shelled parts both
  produce; without it they are drawn in the order they were added.
* `enable_eye_dome_lighting` shades a point cloud by depth, which makes its shape readable where
  flat points are just a haze.
* `add_camera_orientation_widget` puts an interactive orientation cube in the corner, and
  `show_grid` and `add_axes` give the eye something to measure against.
* `Plotter(shape=(1, 2))` with `subplot` and `link_views` puts two states of the same part side by
  side under one camera, which is the quickest way to see what an operation changed.
* `screenshot` and `export_html` get a scene out of the window and into a document.

Anything PyVista accepts can also be passed through a draw method: each one forwards an open
`**kwargs` to the PyVista call named in its docstring. `to_polydata` and `to_mesh3` are there for
the rest of PyVista, so an `engeom` mesh can go through a filter or a widget PyVista offers and
come back, for example `plotter.add_mesh_clip_plane(to_polydata(mesh))`.

Requires `pyvista`, which is an optional dependency: `pip install engeom[pyvista]`.
"""

try:
    # Probe for the dependency before importing the submodules, so that a missing PyVista is
    # reported as itself rather than as a confusing failure deeper in the package.
    import pyvista
except ImportError as _e:
    raise ImportError(
        "engeom.plot.pyvista requires the `pyvista` package, which is not installed. "
        "Install it with `pip install engeom[pyvista]`."
    ) from _e

from .convert import FACE_ID, POINT_ID, LineBuilder, to_mesh3, to_polydata
from .extent import clip_line_to_aabb, clip_plane_to_aabb, resolve_extent
from .helper import PlotterHelper

__all__ = [
    "FACE_ID",
    "POINT_ID",
    "EngeomPlotter",
    "LineBuilder",
    "PlotterHelper",
    "clip_line_to_aabb",
    "clip_plane_to_aabb",
    "register",
    "resolve_extent",
    "to_mesh3",
    "to_polydata",
    "unregister",
]

# The attribute the helper is reached through, on the plotter and on the typed subclass below.
_COMPONENT = "engeom"


def _supported() -> bool:
    """
    Whether this PyVista has the component mechanism the accessor is built on, added in 0.48.
    """
    return hasattr(pyvista, "register_plotter_component")


def _attached() -> bool:
    """
    Whether the accessor is currently attached.

    This reads the class dictionary rather than asking
    `pyvista.registered_plotter_components`, because that function forces PyVista to resolve every
    pending entry point, which would import this module while it is still being imported.
    """
    return _COMPONENT in vars(pyvista.BasePlotter)


def register() -> bool:
    """
    Attach `PlotterHelper` to every PyVista plotter as the `engeom` namespace.

    This runs when the module is imported, so it rarely needs calling directly. The exception is
    putting the accessor back after `unregister`: re-importing the module will not do it, since
    Python does not run a module body twice.

    The component is attached to `pyvista.BasePlotter`, so it reaches every kind of plotter,
    including the Qt plotters from `pyvistaqt`. It is constructed the first time it is read on a
    given plotter and cached there, so `plotter.engeom` is the same helper every time and costs
    nothing after the first access.

    :return: True if the accessor was attached, or False on a PyVista older than 0.48, which has no
        mechanism to attach it to.
    """
    if not _supported():
        return False

    # Re-registering warns unless it is declared as deliberate, which it is whenever this is called
    # to restore the accessor after it was detached.
    pyvista.register_plotter_component(_COMPONENT, override=_attached())(PlotterHelper)
    return True


def unregister() -> bool:
    """
    Detach the `engeom` accessor from PyVista's plotters.

    Equivalent to `pyvista.unregister_plotter_component("engeom")`, and offered here so that the
    way out is documented next to the way in. Helpers already constructed on existing plotters are
    unaffected, and `PlotterHelper(plotter)` keeps working.

    :return: True if the accessor was attached and has been removed, False if there was nothing to
        remove.
    """
    if not (_supported() and _attached()):
        return False

    pyvista.unregister_plotter_component(_COMPONENT)
    return True


class EngeomPlotter(pyvista.Plotter):
    """
    A `pyvista.Plotter` that declares the `engeom` accessor, for editors and type checkers.

    The accessor is attached to every plotter at runtime, but it is attached dynamically, so a type
    checker reading `pyvista.Plotter` has no way to know that `plotter.engeom` exists and an editor
    will not complete it. This subclass adds nothing except the declaration, so that code which
    wants completion can build one of these instead:

    ```python
    from engeom.plot.pyvista import EngeomPlotter
    plotter = EngeomPlotter()
    plotter.engeom.draw_mesh(mesh)
    ```

    The declaration is an annotation and not an assignment, so it does not shadow the descriptor
    that does the work; a plain `pyvista.Plotter` behaves identically at runtime.
    """

    engeom: PlotterHelper


register()
