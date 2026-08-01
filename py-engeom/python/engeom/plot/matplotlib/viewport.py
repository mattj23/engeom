"""
Parallel projection of 3D entities onto a 2D Matplotlib `Axes`, for generating diagrams and
illustrations of 3D geometry.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from engeom.geom3 import Iso3, Line3, Mesh3, Point3, Vector3

from .._coerce import PointLike, to_point2, to_point3
from .trace import TraceBuilder

if TYPE_CHECKING:
    from .axes import AxesHelper


class ViewPort3:
    """
    A helper class that helps draw diagrams of 3d objects being rendered in parallel projection onto a 2d view,
    such as for generating diagrams and illustrations.
    """

    def __init__(self, view: Iso3, helper: AxesHelper):
        self.view = view
        self.helper = helper

    def mesh_edge_point_in_dir(self, view_x: float, view_y: float, mesh: Mesh3) -> Point3:
        """ A quick and rough method for finding a point on the visual perimeter of a mesh in a 2d viewport
        direction to label something without the arrow overlapping the actual mesh. There's probably a better way
        to do this."""
        sample_point = mesh.aabb.center + self.view.inverse() @ (
                Vector3(view_x, view_y, 0) * 100 * mesh.aabb.extent.norm())
        return mesh.point_closest_to(*sample_point)

    def line(self, line: Line3, t: float = 1.0, t0: float | None = None, color: str = "black",
             linewidth: float = 1.0, linestyle: str = "-", alpha: float = 1.0, **kwargs):
        """
        Draws a line in the 2d view using the provided view. If `t0` is `None`, the line will be drawn from `-t` to
        `t`. Otherwise, the line will be drawn from `t0` to `t`.
        :param line: The line to draw.
        :param t: The ending value to draw the line at. If `t0` is `None`, the line will be drawn from `-t` to `t`.
        :param t0: The starting value to draw the line at. If `t0` is `None`, the line will be drawn from `-t` to `t`.
        :param color: The color of the line.
        :param linewidth: The width of the line.
        :param linestyle: The style of the line.
        :param alpha: The transparency of the line.
        :param kwargs: Additional keyword arguments to pass to the line plot function.
        """
        l = self.view @ line
        p0 = l.at(t0 or -t)
        p1 = l.at(t)
        self.helper.ax.plot([p0.x, p1.x], [p0.y, p1.y], color=color, linewidth=linewidth, linestyle=linestyle,
                            alpha=alpha, **kwargs)

    def dimension_arrow(
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
    ):
        ldr_shift = to_point3(leader_shift).coords if leader_shift is not None else Vector3.zero()
        lbl_shift = to_point3(label_shift).coords if label_shift is not None else Vector3.zero()
        p0 = to_point3(point0)
        p1 = to_point3(point1)
        p0_end = p0 + ldr_shift
        p1_end = p1 + ldr_shift
        label_pos = p0 + (p1 - p0) * label_position + ldr_shift + lbl_shift

        self.helper.arrow(self.view @ p0_end, self.view @ p1_end, color=color, linewidth=linewidth, arrow="<->")
        self.helper.text(label, self.view @ label_pos, fontsize=fontsize, ha="center", va="center", color=color,
                         bbox=dict(boxstyle="round", ec="black", fc="w", lw=linewidth))

        if ldr_shift.norm() > 0.0:
            self.helper.arrow(self.view @ p0, self.view @ p0_end, color=color, linewidth=linewidth * 0.75,
                              arrow="-", linestyle="dotted")
            self.helper.arrow(self.view @ p1, self.view @ p1_end, color=color, linewidth=linewidth * 0.75,
                              arrow="-", linestyle="dotted")

    def labeled_point(
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
    ):
        """
        Generates a labeled point representation on a 3D or 2D space. The labeled point is defined
        by its position, label text, offsets for 3D and 2D rendering adjustments, and additional
        visual attributes.

        :param point: The coordinates of the point to be labeled, specified as a type compatible with PointLike.
        :param label: The label text to display near the point.
        :param offset_3d: The 3D offset value to position the label relative to the point in the 3D rendered space.
        :param offset_2d: The 2D offset value to position the label relative to the point in the 2D rendered space.
        :param fontsize: Font size of the label. Defaults to 14.
        :param color: The color of the label text. Defaults to "black".
        :param marker: The marker style used to visually indicate the point. Defaults to "o".
        :param marker_size: Size of the marker representing the point. Defaults to 5.0.
        :param weight: The font weight of the label text. Defaults to "normal".
        :param arrow: Whether to display an arrow from the point to the label. Defaults to False.
        :param arrow_style: The style of the arrow, if arrow is True. Defaults to "->".
        :param linewidth: The width of the arrow line, if arrow is True. Defaults to 1.5.
        :param linestyle: The style of the arrow line, if arrow is True. Defaults to "-".
        :param arrow_props: Extra keyword arguments to pass to the arrow function.
        :param box: Whether to draw a box around the label. Defaults to False.
        """
        p = to_point3(point)
        o3 = Vector3.zero() if offset_3d is None else to_point3(offset_3d).coords
        o2 = Vector3.zero() if offset_2d is None else Vector3(*to_point2(offset_2d), 0)
        text_position = self.view @ (p + o3) + o2

        if arrow:
            self.helper.arrow(text_position, self.view @ p, arrow=arrow_style, color=color, linewidth=linewidth,
                              linestyle=linestyle, **(arrow_props or {}))
        text_params = {
            "fontsize": fontsize,
            "color": color,
            "va": "center",
            "ha": "center",
            "weight": weight,
        }
        if box:
            text_params["bbox"] = dict(boxstyle="round", ec=color, fc="w", lw=linewidth)

        self.helper.text(label, text_position, **text_params)

        if marker_size > 0.0:
            self.helper.points(self.view @ p, marker=marker, markersize=marker_size, color=color)

    def coordinate_system(
            self,
            cs: Iso3,
            length: float,
            linewidth: float = 2.0,
            label_offset: float = 0.1,
            fontsize: int = 14,
    ):
        """
        Draws a coordinate system in the 2d view using the provided view. The coordinate system is drawn as a set
        of red, green, and blue arrows representing the x, y, and z axes respectively. If one of the arrows
        aligns too closely to the view direction, it will be hidden.

        :param cs: The coordinate system to draw
        :param length: The length of the arrows representing the axes
        :param linewidth: The width of the arrows representing the axes
        :param label_offset: The additional scale factor for length to determine where the labels are placed.
        :param fontsize: The font size of the labels
        """

        def _visible(_v: Point3, _len: float) -> bool:
            return _v.with_z(0).coords.norm() > _len * 0.1

        o_view = self.view @ cs.origin
        to_draw = [("x", "red", Vector3.x_axis()),
                   ("y", "green", Vector3.y_axis()),
                   ("z", "blue", Vector3.z_axis())]

        for text, color, vector in to_draw:
            p_view: Point3 = self.view @ cs @ (Point3.origin() + vector * length)
            if _visible(p_view, length):
                self.helper.arrow(o_view, p_view, color=color, linewidth=linewidth)
                text_pos = self.view @ cs @ (Point3.origin() + vector * length * (1 + label_offset))
                self.helper.text(f"${text}$", text_pos, color=color, fontsize=fontsize, ha="center", va="center")

    def mesh_outline(
            self, mesh:
            Mesh3,
            visible_kwargs: dict | None = None,
            hidden_kwargs: dict | None = None,
            no_hidden=False,
            max_edge_len=10.0,
            corner_angle=None
    ):
        """
        Draws the outline of a mesh in the 2d view using the provided view and helper. The default parameters for
        the visible edges are {'color': 'black', 'linewidth': 1.0} and for the hidden edges are
        {'color': 'black', 'linewidth': 0.5, alpha=0.125}.

        :param mesh: The mesh to draw the outline of.
        :param visible_kwargs: Optional keyword arguments for the visible edges of the mesh.
        :param hidden_kwargs: Optional keyword arguments for the hidden edges of the mesh.
        :param no_hidden: If True, only the visible edges will be drawn.
        :param max_edge_len: Max length of edges for the outline, edges longer than this will be broken up
        :param corner_angle: Minimum angle between two adjacent faces for the common edge to be considered a corner and included in the outline
        """
        visible = TraceBuilder()
        hidden = TraceBuilder()

        visible_kwargs = visible_kwargs or {'color': 'black', 'linewidth': 1.0}
        hidden_kwargs = hidden_kwargs or {'color': 'black', 'linewidth': 0.5, 'alpha': 0.125}
        assert visible_kwargs is not None
        assert hidden_kwargs is not None

        points, edge_types = mesh.compute_visual_outline(self.view.inverse() @ Vector3.z_axis(), max_edge_len, corner_angle)
        p0s = self.view.transform_points(points[:, :3])
        p1s = self.view.transform_points(points[:, 3:])
        for edge_type, (p0, p1) in zip(edge_types, zip(p0s, p1s)):
            if edge_type == 0:
                visible.add_segment(p0, p1)
            elif not no_hidden:
                hidden.add_segment(p0, p1)

        self.helper.ax.plot(*visible.xy, **visible_kwargs)
        if not no_hidden:
            self.helper.ax.plot(*hidden.xy, **hidden_kwargs)
