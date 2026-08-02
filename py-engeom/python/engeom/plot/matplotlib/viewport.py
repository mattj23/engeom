"""
Parallel projection of 3D entities onto a 2D Matplotlib `Axes`, for generating diagrams and
illustrations of 3D geometry.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy
from matplotlib.artist import Artist
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Polygon

from engeom.geom2 import Point2, Vector2
from engeom.geom3 import Aabb3, Circle3, Curve3, Iso3, Line3, Mesh3, Plane3, Point3, PointCloud, Vector3
from engeom.metrology import Distance2, Distance3

from .._common import LabelPlace
from .._coerce import PointLike, to_point2, to_point3
from ._style import merge_style
from .trace import TraceBuilder

if TYPE_CHECKING:
    from .axes import AxesHelper

# Below this ratio between a projected direction and its 3D original, the entity is considered to be
# viewed edge-on to the point where the projection carries no usable information.
_DEGENERATE = 1.0e-9


def _plane_basis(normal: Vector3) -> tuple[Vector3, Vector3]:
    """
    Build an orthonormal pair of vectors spanning the plane with the given normal.

    The pair is deterministic but otherwise arbitrary, since a plane normal alone does not fix a
    rotation about itself. The world axis least aligned with the normal is used as the seed, which
    keeps the construction numerically well conditioned for every possible normal.

    :param normal: the plane normal. Need not be unit length.
    :return: two orthonormal vectors, both perpendicular to `normal` and to each other.
    """
    n = normal.normalized()
    seed = min((Vector3.x_axis(), Vector3.y_axis(), Vector3.z_axis()), key=lambda a: abs(n.dot(a)))
    u = n.cross(seed).normalized()
    return u, n.cross(u)


class ViewPort3:
    """
    Draws 3D entities onto a 2D Matplotlib ``Axes`` in parallel projection, for generating diagrams
    and illustrations of 3D geometry.

    This is not a renderer. There is no perspective, no lighting, and no depth sorting beyond what
    the individual methods do for themselves, and entities are drawn in call order. It is intended
    for line diagrams, not for visualizing a scene; use the PyVista helper for that.

    Obtain one from `AxesHelper.viewport` rather than constructing it directly, and note that it
    draws onto that helper's axes, so 2D and 3D entities can be mixed in one figure.

    Every method that adds an artist is prefixed ``draw_``, matching `AxesHelper`.
    """

    def __init__(self, view: Iso3, helper: AxesHelper):
        """
        :param view: the isometry transforming 3D space into the 2D image plane, where +X is to the
            right, +Y is up, and +Z points into the image plane.
        :param helper: the helper whose axes to draw onto.
        """
        self.view = view
        self.helper = helper

    def find_mesh_edge_point(self, view_x: float, view_y: float, mesh: Mesh3) -> Point3:
        """
        Find a point on the visual perimeter of a mesh in a given 2D view direction, so that a
        label can be attached to the silhouette without its leader crossing the mesh itself.

        This is a rough approximation: it projects far outward from the mesh centroid along the
        requested direction and takes the closest point on the mesh to that far sample. It is good
        enough for placing an annotation but is not an exact silhouette query.

        :param view_x: the x component of the direction to search in, in view coordinates.
        :param view_y: the y component of the direction to search in, in view coordinates.
        :param mesh: the mesh to find a perimeter point on.
        :return: a point on the mesh surface, in the mesh's own 3D coordinates rather than view
            coordinates.
        """
        sample_point = mesh.aabb.center + self.view.inverse() @ (
                Vector3(view_x, view_y, 0) * 100 * mesh.aabb.extent.norm())
        return mesh.point_closest_to(*sample_point)

    def draw_line(self, line: Line3, t: float = 1.0, t0: float | None = None, color: str = "black",
                  linewidth: float = 1.0, linestyle: str = "-", alpha: float = 1.0,
                  label: str | None = None, **kwargs) -> Line2D:
        """
        Draw a segment of an infinite 3D line, projected into the view.

        A `Line3` has no endpoints, so a parameter range has to be chosen. By default the segment is
        drawn symmetrically about the line origin, from ``-t`` to ``t``; give `t0` to draw from
        there to `t` instead.

        :param line: the line to draw.
        :param t: the parameter value to stop drawing at.
        :param t0: the parameter value to start drawing at. If None, ``-t`` is used, giving a
            segment centered on the line origin. Note that ``0.0`` is honored, and produces a ray
            from the origin rather than a centered segment.
        :param color: the color of the line.
        :param linewidth: the width of the line.
        :param linestyle: the line style.
        :param alpha: the opacity of the line, from 0.0 to 1.0.
        :param label: the legend label. If None, the line is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: the ``Line2D`` that was drawn.
        """
        kwargs = merge_style(kwargs, label=label)
        l = self.view @ line
        p0 = l.at(t0 if t0 is not None else -t)
        p1 = l.at(t)
        return self.helper.ax.plot([p0.x, p1.x], [p0.y, p1.y], color=color, linewidth=linewidth,
                                   linestyle=linestyle, alpha=alpha, **kwargs)[0]

    def draw_curve(self, *curves: Curve3, color: str | None = None, linewidth: float | None = None,
                   linestyle: str | None = None, alpha: float | None = None,
                   label: str | None = None, **kwargs) -> list[Line2D]:
        """
        Draw one or more 3D curves as projected polylines, each as its own separate line.

        The curve's own vertices are projected directly, so the drawn polyline is exactly the
        projection of the polyline the curve already is, with no resampling.

        :param curves: the curves to draw.
        :param color: the line color. If None, the axes' color cycle supplies one.
        :param linewidth: the line width in points. If None, the Matplotlib default applies.
        :param linestyle: the line style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the line is fully opaque.
        :param label: the legend label. If None, the line is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per curve, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)
        drawn = []
        for curve in curves:
            points = self.view.transform_points(curve.points)
            drawn.append(self.helper.ax.plot(points[:, 0], points[:, 1], **kwargs)[0])
        return drawn

    def draw_point_cloud(self, *clouds: PointCloud, marker: str = ".", markersize: float = 1.0,
                         color: str | None = None, alpha: float | None = None,
                         label: str | None = None, **kwargs) -> list[Line2D]:
        """
        Draw one or more point clouds as projected markers.

        Points are drawn in storage order with no depth sorting, so a cloud with a far side will
        show it through the near side. For a cloud dense enough that this matters, a small marker
        with some transparency reads better than trying to hide the back.

        Per-point colors carried by the cloud are ignored, since a single Matplotlib line artist has
        one color. Reach for the PyVista helper if the colors are the point of the plot.

        :param clouds: the point clouds to draw.
        :param marker: the Matplotlib marker style for each point.
        :param markersize: the marker size in points, defaulting small because clouds are dense.
        :param color: the marker color. If None, the axes' color cycle supplies one.
        :param alpha: the opacity from 0.0 to 1.0. If None, the markers are fully opaque.
        :param label: the legend label. If None, the cloud is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per cloud, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, alpha=alpha, label=label)
        drawn = []
        for cloud in clouds:
            points = self.view.transform_points(cloud.points)
            drawn.append(self.helper.ax.plot(points[:, 0], points[:, 1], linestyle="none",
                                             marker=marker, markersize=markersize, **kwargs)[0])
        return drawn

    def draw_circle(self, *circles: Circle3, fill: bool = False, color: str | None = None,
                    edgecolor: str | None = None, facecolor: str | None = None,
                    linewidth: float | None = None, linestyle: str | None = None,
                    alpha: float | None = None, label: str | None = None,
                    **kwargs) -> list[Ellipse]:
        """
        Draw one or more 3D circles as the ellipses they project to.

        A circle not facing the viewer projects to an ellipse under parallel projection, and that
        ellipse is computed exactly rather than approximated by a polyline: the projection maps the
        circle's plane basis to a 2x2 matrix whose singular values are the ellipse's semi-axes. A
        circle seen exactly edge-on therefore collapses to a line, which is correct.

        :param circles: the circles to draw.
        :param fill: whether to fill the ellipses.
        :param color: sets both the edge and face color at once. If None, the Matplotlib default
            applies. `edgecolor` and `facecolor` override this for the side they name.
        :param edgecolor: the outline color. If None, `color` or the Matplotlib default applies.
        :param facecolor: the fill color, which has no effect unless `fill` is True. If None,
            `color` or the Matplotlib default applies.
        :param linewidth: the outline width in points. If None, the Matplotlib default applies.
        :param linestyle: the outline style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the ellipse is fully opaque.
        :param label: the legend label. If None, the ellipse is not labeled.
        :param kwargs: any other keyword argument accepted by the Matplotlib ``Ellipse`` patch.
        :return: one ``Ellipse`` patch per circle, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, edgecolor=edgecolor, facecolor=facecolor,
                             linewidth=linewidth, linestyle=linestyle, alpha=alpha, label=label)
        patches = []
        for circle in circles:
            in_view = self.view @ circle
            u, v = _plane_basis(in_view.normal)
            # The circle is the image of the unit circle under [r*u | r*v]; dropping z leaves a 2x2
            # whose singular values are the projected semi-axes and whose left singular vectors give
            # their directions.
            matrix = numpy.array([[u.x, v.x], [u.y, v.y]]) * in_view.r
            left, semi_axes, _ = numpy.linalg.svd(matrix)
            angle = numpy.rad2deg(numpy.arctan2(left[1, 0], left[0, 0]))
            center = in_view.center
            patch = Ellipse((center.x, center.y), 2 * semi_axes[0], 2 * semi_axes[1], angle=angle,
                            fill=fill, **kwargs)
            patches.append(self.helper.ax.add_patch(patch))
        return patches

    def draw_plane(self, *planes: Plane3, center: PointLike, size: float, fill: bool = False,
                   color: str | None = None, edgecolor: str | None = None,
                   facecolor: str | None = None, linewidth: float | None = None,
                   linestyle: str | None = None, alpha: float | None = None,
                   label: str | None = None, **kwargs) -> list[Polygon]:
        """
        Draw one or more infinite planes as bounded square quads, projected into the view.

        A `Plane3` has no extent, so both where to draw it and how big to draw it have to be given.
        `center` is projected onto each plane to place the quad, meaning the same anchor point can
        be reused for several planes without having to sit on any of them. The quad's rotation about
        the plane normal is deterministic but arbitrary, since a plane does not carry one.

        A quad's projection is exact, unlike a circle's, because parallel projection maps straight
        edges to straight edges.

        :param planes: the planes to draw.
        :param center: the point to center each quad on, projected onto the plane first. Any
            coordinate the helpers accept.
        :param size: the edge length of the square quad, in data units.
        :param fill: whether to fill the quads.
        :param color: sets both the edge and face color at once. If None, the Matplotlib default
            applies. `edgecolor` and `facecolor` override this for the side they name.
        :param edgecolor: the outline color. If None, `color` or the Matplotlib default applies.
        :param facecolor: the fill color, which has no effect unless `fill` is True. If None,
            `color` or the Matplotlib default applies.
        :param linewidth: the outline width in points. If None, the Matplotlib default applies.
        :param linestyle: the outline style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. A filled plane is usually worth making partly
            transparent so that geometry behind it stays visible.
        :param label: the legend label. If None, the quad is not labeled.
        :param kwargs: any other keyword argument accepted by the Matplotlib ``Polygon`` patch.
        :return: one ``Polygon`` patch per plane, in the order given.
        :raises ValueError: if `size` is not positive.
        """
        if size <= 0.0:
            raise ValueError(f"size must be positive, got {size}")
        kwargs = merge_style(kwargs, color=color, edgecolor=edgecolor, facecolor=facecolor,
                             linewidth=linewidth, linestyle=linestyle, alpha=alpha, label=label)
        anchor = to_point3(center)
        half = size / 2.0
        patches = []
        for plane in planes:
            origin = plane.project_point(anchor)
            u, v = _plane_basis(plane.normal)
            corners = [origin + u * half + v * half, origin - u * half + v * half,
                       origin - u * half - v * half, origin + u * half - v * half]
            projected = [self.view @ corner for corner in corners]
            patch = Polygon([(p.x, p.y) for p in projected], closed=True, fill=fill, **kwargs)
            patches.append(self.helper.ax.add_patch(patch))
        return patches

    def draw_aabb(self, *boxes: Aabb3, color: str | None = None, linewidth: float | None = None,
                  linestyle: str | None = None, alpha: float | None = None,
                  label: str | None = None, **kwargs) -> list[Line2D]:
        """
        Draw one or more 3D bounding boxes as projected wireframes.

        All twelve edges of a box go into a single artist, using `TraceBuilder` to break between
        them, so the box restyles as one thing. Every edge is drawn, including the three at the
        back, since hiding them would need a depth test this viewport does not do.

        :param boxes: the bounding boxes to draw.
        :param color: the line color. If None, the axes' color cycle supplies one.
        :param linewidth: the line width in points. If None, the Matplotlib default applies.
        :param linestyle: the line style, such as ``"--"``. If None, the Matplotlib default applies.
        :param alpha: the opacity from 0.0 to 1.0. If None, the wireframe is fully opaque.
        :param label: the legend label. If None, the box is not labeled.
        :param kwargs: any other keyword argument accepted by ``Axes.plot``.
        :return: one ``Line2D`` per box, in the order given.
        """
        kwargs = merge_style(kwargs, color=color, linewidth=linewidth, linestyle=linestyle,
                             alpha=alpha, label=label)
        drawn = []
        for box in boxes:
            axes = ((box.min.x, box.max.x), (box.min.y, box.max.y), (box.min.z, box.max.z))
            corners = {(i, j, k): Point3(axes[0][i], axes[1][j], axes[2][k])
                       for i in (0, 1) for j in (0, 1) for k in (0, 1)}
            trace = TraceBuilder()
            for a in sorted(corners):
                for b in sorted(corners):
                    # One differing index means the two corners share an edge; `a < b` keeps each
                    # of the twelve edges from being walked in both directions.
                    if a < b and sum(1 for x, y in zip(a, b) if x != y) == 1:
                        trace.add_segment(self.view @ corners[a], self.view @ corners[b])
            drawn.append(self.helper.ax.plot(*trace.xy, **kwargs)[0])
        return drawn

    def draw_distance(
            self,
            distance: Distance3,
            side_shift: float = 0.0,
            template: str = "{value:.3f}",
            fontsize: int = 10,
            label_place: LabelPlace = "outside",
            label_offset: float | None = None,
            fontname: str | None = None,
            scale_value: float = 1.0,
    ) -> list[Artist]:
        """
        Draw a `Distance3` as a dimension, by projecting it into the view and handing it to
        `AxesHelper.draw_distance`, so that leaders, arrows, and label placement all behave exactly
        as they do for a 2D dimension.

        **The label reports the true 3D distance, not the projected one.** Projection foreshortens
        any measurement not parallel to the image plane, and a dimension whose label disagreed with
        the geometry it measures would be worse than useless. The arrow is drawn foreshortened, as
        it must be, but the number is the real one.

        :param distance: the 3D measurement to draw.
        :param side_shift: offset the leader line ends by this many data units, orthogonally to the
            projected measurement direction, with positive values shifting to the right.
        :param template: the format string for the label. The measured value is substituted for the
            ``value`` key.
        :param fontsize: the font size of the label.
        :param label_place: where to put the label relative to the anchor points. See `LabelPlace`.
        :param label_offset: how far the label sits from its anchor, in data units. If None, a
            default derived from the font size is used.
        :param fontname: the name of the font to use for the label. If None, the Matplotlib default
            is used.
        :param scale_value: a factor applied to the true 3D value before it is formatted into the
            label, for displaying the measurement in different units.
        :return: every artist drawn, in the same order as `AxesHelper.draw_distance` returns them.
        :raises ValueError: if the measurement direction is so close to the view direction that the
            dimension projects to nothing and cannot be drawn.
        :raises ValueError: if `label_place` is not one of the valid placement tokens.
        """
        a = self.view @ distance.a
        b = self.view @ distance.b
        direction = self.view @ distance.direction
        flat_direction = Vector2(direction.x, direction.y)

        if flat_direction.norm() <= _DEGENERATE * direction.norm():
            raise ValueError(
                "the measurement direction is parallel to the view direction, so the dimension "
                "projects to a point and cannot be drawn; view the geometry from another angle"
            )

        flat = Distance2(Point2(a.x, a.y), Point2(b.x, b.y), flat_direction.normalized())
        if distance.value != 0.0 and abs(flat.value) <= _DEGENERATE * abs(distance.value):
            raise ValueError(
                "the measurement projects to a zero-length dimension and cannot be drawn; view the "
                "geometry from another angle"
            )

        # `draw_distance` formats whatever the 2D value is, which foreshortening has shrunk. Fold
        # the recovery factor into the scale so the label shows the true 3D distance. A measurement
        # that is genuinely zero in 3D needs no recovery, and labels as zero either way.
        correction = 1.0 if flat.value == 0.0 else distance.value / flat.value
        return self.helper.draw_distance(
            flat,
            side_shift=side_shift,
            template=template,
            fontsize=fontsize,
            label_place=label_place,
            label_offset=label_offset,
            fontname=fontname,
            scale_value=scale_value * correction,
        )

    def draw_dimension_arrow(
            self,
            point0: PointLike,
            point1: PointLike,
            label: str,
            leader_shift: PointLike | None = None,
            label_shift: PointLike | None = None,
            linewidth: float = 1.5,
            fontsize: int = 14,
            color: str = "black",
            label_position: float = 0.5,
    ) -> list[Artist]:
        """
        Draw a dimension between two 3D points as a double-headed arrow with a boxed label,
        optionally standing the dimension off from the geometry it measures.

        When `leader_shift` is given, the dimension arrow is drawn at the shifted position and
        dotted leader lines are added connecting it back to the two original points, which is the
        usual way to keep a dimension clear of the part it refers to.

        :param point0: the first anchor point, as any coordinate the helpers accept.
        :param point1: the second anchor point.
        :param label: the text to show in the label box.
        :param leader_shift: a 3D offset applied to the dimension arrow, moving it away from the
            anchor points. If given and non-zero, dotted leaders back to the anchors are drawn too.
        :param label_shift: a further 3D offset applied to the label alone, for nudging the text
            clear of the arrow.
        :param linewidth: the width of the dimension arrow. Leaders are drawn at three quarters of
            this width.
        :param fontsize: the font size of the label.
        :param color: the color of the arrows and the label text.
        :param label_position: where the label sits along the dimension, as a fraction from
            `point0` to `point1`.
        :return: every artist drawn, in draw order: the dimension arrow, the label, then the two
            dotted leaders if a `leader_shift` was given.
        """
        ldr_shift = to_point3(leader_shift).coords if leader_shift is not None else Vector3.zero()
        lbl_shift = to_point3(label_shift).coords if label_shift is not None else Vector3.zero()
        p0 = to_point3(point0)
        p1 = to_point3(point1)
        p0_end = p0 + ldr_shift
        p1_end = p1 + ldr_shift
        label_pos = p0 + (p1 - p0) * label_position + ldr_shift + lbl_shift

        artists: list[Artist] = [
            self.helper.draw_arrow(self.view @ p0_end, self.view @ p1_end, color=color, linewidth=linewidth,
                                   arrow="<->"),
            self.helper.draw_text(label, self.view @ label_pos, fontsize=fontsize, ha="center", va="center",
                                  color=color,
                                  bbox=dict(boxstyle="round", ec="black", fc="w", lw=linewidth)),
        ]

        if ldr_shift.norm() > 0.0:
            artists.append(self.helper.draw_arrow(self.view @ p0, self.view @ p0_end, color=color,
                                                  linewidth=linewidth * 0.75, arrow="-", linestyle="dotted"))
            artists.append(self.helper.draw_arrow(self.view @ p1, self.view @ p1_end, color=color,
                                                  linewidth=linewidth * 0.75, arrow="-", linestyle="dotted"))
        return artists

    def draw_labeled_point(
            self,
            point: PointLike,
            label: str,
            offset_3d: PointLike | None = None,
            offset_2d: PointLike | None = None,
            fontsize: int = 14,
            color: str = "black",
            marker: str = "o",
            marker_size: float = 5.0,
            weight: str = "normal",
            arrow: bool = False,
            arrow_style: str = "->",
            linewidth: float = 1.5,
            linestyle: str = "-",
            arrow_props: dict | None = None,
            box: bool = False,
    ) -> list[Artist]:
        """
        Draw a marker at a 3D point with a text label beside it, optionally with a leader arrow
        from the label back to the point.

        The label can be offset in either space, and the two compose. `offset_3d` moves the label
        in 3D before projection, so it stays anchored to the geometry as the view rotates;
        `offset_2d` moves it on the page after projection, so it always displaces in the same
        screen direction. Use the former to place a label relative to the part, the latter to keep
        text from colliding in a particular view.

        :param point: the point to label, as any coordinate the helpers accept.
        :param label: the label text.
        :param offset_3d: an offset applied to the label position in 3D, before projection.
        :param offset_2d: an offset applied to the label position in 2D, after projection.
        :param fontsize: the font size of the label.
        :param color: the color of the label, marker, and leader arrow.
        :param marker: the Matplotlib marker style for the point.
        :param marker_size: the marker size in points. Set to 0.0 to draw the label with no marker.
        :param weight: the font weight of the label.
        :param arrow: whether to draw a leader arrow from the label back to the point. Only useful
            in combination with an offset, since otherwise the two coincide.
        :param arrow_style: the Matplotlib arrow style, used when `arrow` is True.
        :param linewidth: the width of the leader arrow.
        :param linestyle: the line style of the leader arrow.
        :param arrow_props: additional keyword arguments to pass to the arrow.
        :param box: whether to draw a rounded box around the label text.
        :return: every artist drawn, in draw order: the leader arrow if requested, the label, then
            the marker if `marker_size` is positive.
        """
        p = to_point3(point)
        o3 = Vector3.zero() if offset_3d is None else to_point3(offset_3d).coords
        o2 = Vector3.zero() if offset_2d is None else Vector3(*to_point2(offset_2d), 0)
        text_position = self.view @ (p + o3) + o2

        artists: list[Artist] = []
        if arrow:
            artists.append(self.helper.draw_arrow(text_position, self.view @ p, arrow=arrow_style, color=color,
                                                  linewidth=linewidth, linestyle=linestyle,
                                                  **(arrow_props or {})))
        text_params = {
            "fontsize": fontsize,
            "color": color,
            "va": "center",
            "ha": "center",
            "weight": weight,
        }
        if box:
            text_params["bbox"] = dict(boxstyle="round", ec=color, fc="w", lw=linewidth)

        artists.append(self.helper.draw_text(label, text_position, **text_params))

        if marker_size > 0.0:
            artists.append(self.helper.draw_point(self.view @ p, marker=marker, markersize=marker_size, color=color))
        return artists

    def draw_coordinate_system(
            self,
            cs: Iso3,
            length: float,
            linewidth: float = 2.0,
            label_offset: float = 0.1,
            fontsize: int = 14,
    ) -> list[Artist]:
        """
        Draw a coordinate frame as three labeled arrows: X in red, Y in green, and Z in blue.

        An axis pointing nearly straight into or out of the view would project to almost nothing
        and render as a meaningless stub, so any axis whose projected length falls below a tenth of
        `length` is skipped entirely. Looking down an axis therefore yields two arrows, not three.

        :param cs: the coordinate frame to draw, as an isometry.
        :param length: the length of each axis arrow, in data units.
        :param linewidth: the width of the axis arrows.
        :param label_offset: how far beyond the arrow tip to place its label, as a fraction of
            `length`.
        :param fontsize: the font size of the axis labels.
        :return: every artist drawn: an arrow and a label for each visible axis, so between zero
            and six artists depending on the view direction.
        """

        def _visible(_v: Point3, _len: float) -> bool:
            return _v.with_z(0).coords.norm() > _len * 0.1

        o_view = self.view @ cs.origin
        to_draw = [("x", "red", Vector3.x_axis()),
                   ("y", "green", Vector3.y_axis()),
                   ("z", "blue", Vector3.z_axis())]

        artists: list[Artist] = []
        for text, color, vector in to_draw:
            p_view: Point3 = self.view @ cs @ (Point3.origin() + vector * length)
            if _visible(p_view, length):
                artists.append(self.helper.draw_arrow(o_view, p_view, color=color, linewidth=linewidth))
                text_pos = self.view @ cs @ (Point3.origin() + vector * length * (1 + label_offset))
                artists.append(self.helper.draw_text(f"${text}$", text_pos, color=color, fontsize=fontsize,
                                                     ha="center", va="center"))
        return artists

    def draw_mesh_outline(
            self,
            mesh: Mesh3,
            visible_kwargs: dict | None = None,
            hidden_kwargs: dict | None = None,
            no_hidden=False,
            max_edge_len=10.0,
            corner_angle=None
    ) -> list[Line2D]:
        """
        Draw a mesh as a line drawing of its silhouette and corner edges, in the manner of a
        technical illustration rather than a shaded render.

        Edges are classified as visible or hidden with respect to the view direction and drawn with
        separate styles, so that hidden detail can be shown faintly or suppressed. All the edges of
        one class go into a single artist, using `TraceBuilder` to break between disjoint runs.

        :param mesh: the mesh to outline.
        :param visible_kwargs: keyword arguments for the visible edges, passed to ``Axes.plot``. If
            None, ``{'color': 'black', 'linewidth': 1.0}`` is used. Pass an empty dict to fall back
            to Matplotlib's own defaults instead.
        :param hidden_kwargs: keyword arguments for the hidden edges. If None,
            ``{'color': 'black', 'linewidth': 0.5, 'alpha': 0.125}`` is used.
        :param no_hidden: if True, hidden edges are not drawn at all and only one artist is
            returned.
        :param max_edge_len: the maximum length of an outline edge, in data units. Longer edges are
            subdivided, which matters because a long straight edge on a curved surface would
            otherwise cut across the silhouette.
        :param corner_angle: the minimum angle between adjacent faces, in radians, for their shared
            edge to count as a corner and be included in the outline. If None, only silhouette
            edges are drawn and interior corners are omitted.
        :return: the visible-edge ``Line2D``, followed by the hidden-edge ``Line2D`` unless
            `no_hidden` was set.
        """
        visible = TraceBuilder()
        hidden = TraceBuilder()

        if visible_kwargs is None:
            visible_kwargs = {'color': 'black', 'linewidth': 1.0}
        if hidden_kwargs is None:
            hidden_kwargs = {'color': 'black', 'linewidth': 0.5, 'alpha': 0.125}

        points, edge_types = mesh.compute_visual_outline(self.view.inverse() @ Vector3.z_axis(), max_edge_len, corner_angle)
        p0s = self.view.transform_points(points[:, :3])
        p1s = self.view.transform_points(points[:, 3:])
        for edge_type, (p0, p1) in zip(edge_types, zip(p0s, p1s)):
            if edge_type == 0:
                visible.add_segment(p0, p1)
            elif not no_hidden:
                hidden.add_segment(p0, p1)

        lines = [self.helper.ax.plot(*visible.xy, **visible_kwargs)[0]]
        if not no_hidden:
            lines.append(self.helper.ax.plot(*hidden.xy, **hidden_kwargs)[0])
        return lines
