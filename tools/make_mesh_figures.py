"""Regenerate the mesh-generation figures for docs/fem/mesh.md.

Four figures, each answering one question the page asks the reader to decide:

  mesh_element_nodes.png   What the five element types are, and where their
                           nodes sit — the local numbering the mesh arrays use.
  mesh_tri_vs_quad.png     What the element-type choice looks like on a real
                           section: the same layered slope as tri3 and quad4.
  mesh_zone_size.png       What a per-zone Size does: a thin weak layer that the
                           global size cannot fit elements across, resolved by a
                           Size on that zone alone.
  mesh_refine_features.png What refine-near-features does: elements shrink along
                           the reinforcement lines and grow back away from them.

Every panel is built by meshing a committed sample file at a stated element
size, so the figures are reproducible and stay honest as the mesher changes.

Run from the repo root:

    PYTHONPATH=. python3 tools/make_mesh_figures.py
"""
import contextlib
import io
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Patch

from xslope.fileio import load_slope_data
from xslope.mesh import (build_mesh_from_polygons, extract_constraint_line_geometry,
                         get_material_polygons)
from xslope.plot import get_material_color

IMAGES = "docs/fem/images"

DPI = 200
PAD_FRAC = 0.035        # uniform cushion around the domain (house frame spec)
FIG_W_IN = 9.5
TITLE_H_IN = 0.32
CAPTION_H_IN = 0.34
LEGEND_H_IN = 0.50


# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------

def _quiet(fn, *args, **kwargs):
    """gmsh and the loaders are chatty; the figure is the output, not the log."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _corners(mesh, i):
    """Corner nodes of element ``i`` — 4 for a quad, 3 for a triangle. Midside
    nodes change neither the outline nor the shape statistics."""
    n = 4 if mesh["element_types"][i] in (4, 8, 9) else 3
    return mesh["nodes"][mesh["elements"][i][:n].astype(int)]


def _load(model):
    """Sample file -> (polygons, constraint lines, material names, ground width)."""
    root = _root()
    data = _quiet(load_slope_data, os.path.join(root, model))
    lines, _, _ = _quiet(extract_constraint_line_geometry, data)
    polygons = _quiet(get_material_polygons, data, reinf_lines=lines or None)
    gx = [x for x, _ in data["ground_surface"].coords]
    names = [m.get("name") or f"Material {i + 1}"
             for i, m in enumerate(data.get("materials") or [])]
    return polygons, lines, names, max(gx) - min(gx)


def _root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _extent(mesh):
    nodes = mesh["nodes"]
    return (nodes[:, 0].min(), nodes[:, 0].max(),
            nodes[:, 1].min(), nodes[:, 1].max())


def _draw(ax, mesh, names, window=None, lw=0.25):
    """Material-colored element patches with hairline edges, equal aspect.

    ``window`` limits the view to (x0, x1, y0, y1) for a detail panel; otherwise
    the whole mesh is shown with a uniform cushion. Returns legend handles."""
    mats = np.asarray(mesh["element_materials"])
    handles = []
    for mid in sorted(set(mats.tolist())):
        polys = [_corners(mesh, i) for i in np.flatnonzero(mats == mid)]
        # Mesh material IDs are 1-based (gmsh); the palette is 0-based, so the
        # zone colors match the Inputs view and every other xslope plot.
        color = get_material_color(int(mid) - 1)
        ax.add_collection(PolyCollection(polys, facecolors=[color], alpha=0.55,
                                         edgecolors="0.25", linewidths=lw))
        label = names[mid - 1] if 0 <= mid - 1 < len(names) else f"Material {mid}"
        handles.append(Patch(facecolor=color, alpha=0.55, edgecolor="0.25",
                             label=label))
    if window is None:
        x0, x1, y0, y1 = _extent(mesh)
        px, py = PAD_FRAC * (x1 - x0), PAD_FRAC * (y1 - y0)
        ax.set_xlim(x0 - px, x1 + px)
        ax.set_ylim(y0 - py, y1 + py)
    else:
        ax.set_xlim(window[0], window[1])
        ax.set_ylim(window[2], window[3])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("0.6")
    return handles


def _stacked_figure(panels, names, out, extras=None, windows=None):
    """A column of equal-aspect mesh panels, each with a bold title and a caption.

    ``panels`` is [(title, caption, mesh), ...]; every panel is drawn at the same
    scale, so the comparison is between the meshes and not between two zooms.
    ``extras`` is an optional per-panel callable(ax) for overlays. ``windows``
    optionally restricts every panel to one (x0, x1, y0, y1) view."""
    left, right = 0.03, 0.99
    axes_w_in = FIG_W_IN * (right - left)
    if windows is None:
        x0, x1, y0, y1 = _extent(panels[0][2])
        span_x = (x1 - x0) * (1 + 2 * PAD_FRAC)
        span_y = (y1 - y0) * (1 + 2 * PAD_FRAC)
    else:
        span_x = windows[1] - windows[0]
        span_y = windows[3] - windows[2]
    panel_h_in = axes_w_in * span_y / span_x
    row_h_in = panel_h_in + CAPTION_H_IN + TITLE_H_IN
    fig_h_in = len(panels) * row_h_in + LEGEND_H_IN + 0.15
    fig = plt.figure(figsize=(FIG_W_IN, fig_h_in), dpi=DPI)

    handles = []
    for i, (title, caption, mesh) in enumerate(panels):
        bottom = (LEGEND_H_IN + (len(panels) - 1 - i) * row_h_in
                  + CAPTION_H_IN) / fig_h_in
        ax = fig.add_axes((left, bottom, right - left, panel_h_in / fig_h_in))
        handles = _draw(ax, mesh, names, window=windows)
        if extras:
            extras[i](ax)
        ax.set_title(title, fontsize=11, fontweight="bold", pad=6)
        fig.text(left + (right - left) / 2,
                 bottom - CAPTION_H_IN / fig_h_in * 0.62, caption,
                 ha="center", va="center", fontsize=8.5, color="0.25")

    fig.legend(handles=handles, loc="lower center", ncol=max(1, len(handles)),
               frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.10 / fig_h_in))
    _save(fig, out)


def _save(fig, out):
    path = os.path.join(_root(), IMAGES, out)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"wrote {IMAGES}/{out}")


def _stats(mesh):
    return f"{len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements"


# --------------------------------------------------------------------------
# 1. element types and local node numbering
# --------------------------------------------------------------------------

TRI = np.array([(0.0, 0.0), (1.0, 0.0), (0.5, 0.88)])
QUAD = np.array([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])


def _element_schematic(kind):
    """Corner outline plus the node coordinates in xslope's local order.

    Corners come first, counterclockwise; then the midside nodes in edge order
    (0-1, 1-2, ...); then, for quad9, the center. This is the order the mesh
    ``elements`` array stores, so the labels are the array columns."""
    base = TRI if kind.startswith("tri") else QUAD
    nodes = [tuple(p) for p in base]
    if kind in ("tri6", "quad8", "quad9"):
        n = len(base)
        nodes += [tuple((base[i] + base[(i + 1) % n]) / 2) for i in range(n)]
    if kind == "quad9":
        nodes.append(tuple(base.mean(axis=0)))
    return base, np.array(nodes)


def figure_element_nodes():
    kinds = [("tri3", "3-node triangle"), ("tri6", "6-node triangle"),
             ("quad4", "4-node quadrilateral"), ("quad8", "8-node quadrilateral"),
             ("quad9", "9-node quadrilateral")]
    fig_w, panel_h = 10.0, 1.95
    fig = plt.figure(figsize=(fig_w, panel_h + 0.55), dpi=DPI)
    for i, (kind, title) in enumerate(kinds):
        ax = fig.add_axes((i / len(kinds) + 0.012, 0.10,
                           1 / len(kinds) - 0.024, panel_h / (panel_h + 0.55) * 0.86))
        outline, nodes = _element_schematic(kind)
        ring = np.vstack([outline, outline[:1]])
        ax.fill(ring[:, 0], ring[:, 1], color="#cfe3ee", alpha=0.85, zorder=0)
        ax.plot(ring[:, 0], ring[:, 1], color="0.2", lw=1.4, zorder=1)
        ax.scatter(nodes[:, 0], nodes[:, 1], s=42, color="#c0392b", zorder=3)
        center = nodes[:len(outline)].mean(axis=0)
        for n, (x, y) in enumerate(nodes):
            # Push each label away from the element center so it never sits on
            # an edge; the center node of quad9 has nowhere to go but up.
            d = np.array([x, y]) - center
            if np.hypot(*d) < 1e-9:
                off = np.array([0.0, 0.16])
            else:
                off = 0.20 * d / np.hypot(*d)
            ax.annotate(str(n), (x + off[0], y + off[1]), ha="center", va="center",
                        fontsize=9.5, fontweight="bold", color="#7b241c")
        ax.set_title(f"{kind}\n{title}", fontsize=9.5, fontweight="bold", pad=4,
                     linespacing=1.35)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-0.34, 1.34)
        ax.set_ylim(-0.30, 1.38)
        ax.axis("off")
    _save(fig, "mesh_element_nodes.png")


# --------------------------------------------------------------------------
# 2. triangles vs quadrilaterals on the same section
# --------------------------------------------------------------------------

TRIQUAD_MODEL = "docs/fem/files/xslope_simple_mult_layers_fem.xlsx"


def figure_tri_vs_quad():
    polygons, _, names, width = _load(TRIQUAD_MODEL)
    target = width / 60
    panels = []
    for kind, title in (("tri3", "tri3 — 3-node triangles"),
                        ("quad4", "quad4 — 4-node quadrilaterals")):
        mesh = _quiet(build_mesh_from_polygons, polygons, target_size=target,
                      element_type=kind)
        panels.append((title, _stats(mesh), mesh))
    _stacked_figure(panels, names, "mesh_tri_vs_quad.png")
    print(f"   target element size {target:.3f}")


# --------------------------------------------------------------------------
# 3. a per-zone Size on a thin layer
# --------------------------------------------------------------------------

SIZE_MODEL = "docs/seep/files/xslope_earth_dam1.xlsx"
SIZE_ZONE_MAT = 1          # the central clay core (0-based mat_id)
SIZE_TARGET = 1.1          # global target element size
SIZE_ZONE = 0.55           # Size declared on the core alone


def _median_edge(mesh, mat_id):
    """Median element edge length inside one material zone.

    Edge length rather than area, so the number in the caption is directly
    comparable with the element size that was requested."""
    mats = np.asarray(mesh["element_materials"])
    lengths = []
    for i in np.flatnonzero(mats == mat_id + 1):     # mesh material IDs are 1-based
        v = _corners(mesh, i)
        lengths.append(np.linalg.norm(np.roll(v, -1, axis=0) - v, axis=1).mean())
    return float(np.median(lengths))


def figure_zone_size():
    polygons, _, names, _ = _load(SIZE_MODEL)
    refined = [dict(p, size=SIZE_ZONE) if p["mat_id"] == SIZE_ZONE_MAT else p
               for p in polygons]
    other = 1 - SIZE_ZONE_MAT

    panels = []
    for title, polys in ((f"Global target size {SIZE_TARGET:g}, no zone Size",
                          polygons),
                         (f"Size = {SIZE_ZONE:g} on the core", refined)):
        mesh = _quiet(build_mesh_from_polygons, polys, target_size=SIZE_TARGET,
                      element_type="tri3")
        caption = (f"{_stats(mesh)}   |   median element edge "
                   f"{_median_edge(mesh, SIZE_ZONE_MAT):.2f} in the core, "
                   f"{_median_edge(mesh, other):.2f} in the shell")
        panels.append((title, caption, mesh))
    _stacked_figure(panels, names, "mesh_zone_size.png")
    for t, c, _ in panels:
        print(f"   {t:38s} {c}")


# --------------------------------------------------------------------------
# 4. refine near features
# --------------------------------------------------------------------------

REFINE_MODEL = "docs/fem/files/xslope_reinforce_fem.xlsx"
REFINE_FACTOR = 3.0


def figure_refine_features():
    polygons, lines, names, width = _load(REFINE_MODEL)
    target = width / 45

    def overlay(ax):
        segs = [np.asarray(ln) for ln in lines]
        ax.add_collection(LineCollection(segs, colors="#111111", linewidths=1.6,
                                         zorder=6))

    panels, extras = [], []
    for factor, title in ((None, "Off (default) — one element size everywhere"),
                          (REFINE_FACTOR,
                           f"refine_factor = {REFINE_FACTOR:g} — elements shrink "
                           "near the reinforcement")):
        mesh = _quiet(build_mesh_from_polygons, polygons, target_size=target,
                      element_type="tri3", lines=lines or None,
                      refine_factor=factor)
        panels.append((title, _stats(mesh), mesh))
        extras.append(overlay)
    _stacked_figure(panels, names, "mesh_refine_features.png", extras=extras)
    print(f"   target element size {target:.3f}")


def main():
    figure_element_nodes()
    figure_tri_vs_quad()
    figure_zone_size()
    figure_refine_features()


if __name__ == "__main__":
    main()
