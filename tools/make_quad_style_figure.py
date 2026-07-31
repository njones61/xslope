"""Regenerate the quadrilateral-meshing-style comparison figure for docs/fem/mesh.md.

The page's "Quadrilateral Meshing Styles" section asks the reader to choose
between the two styles the mesh builder offers, so it shows them: the same
section, at the same requested element size, meshed each way, with the numbers
that separate them printed under each panel.

The section is the levee with a grouted foundation
(``docs/seep/files/xslope_levee_poly.xlsx``) because it carries both outcomes at
once — its three foundation blocks are swept as exact grids under
``quad_style='structured'``, while its 5:1 trapezoidal embankment is not
sweepable and stays free-meshed. That is the figure's argument: "structured"
means structured where the shape allows.

Mesh configuration (the Studio defaults, so the panels are what a reader
reproduces by opening the file and pressing Build):

    element_type = 'quad4'
    target_size  = (ground-surface x-extent) / 100   = 0.580
    quad_style   = 'free' | 'structured'

Metrics under each panel use the same definitions as ``test/quad_mesh_check.py``
(element aspect ratio = longest corner edge / shortest corner edge; element size
= sqrt(area)), so the figure and the standing mesh checks measure one thing.

Figure sizing follows the house frame spec: equal aspect at one common scale for
both panels, a uniform cushion around the domain, and the material legend in
reserved space rather than over the mesh.

Run from the repo root:

    PYTHONPATH=. python3 tools/make_quad_style_figure.py
"""
import contextlib
import io
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch

from xslope.fileio import load_slope_data
from xslope.mesh import build_mesh_from_polygons, get_material_polygons
from xslope.plot import get_material_color

MODEL = "docs/seep/files/xslope_levee_poly.xlsx"
OUT = "docs/fem/images/quad_styles_levee.png"

STYLES = [("free", "Free (recommended)"),
          ("structured", "Structured where possible")]

DIVISIONS = 100        # target size = ground-surface width / DIVISIONS
FIG_W_IN = 10.0
PAD_FRAC = 0.035       # uniform cushion around the domain (frame spec)
LEGEND_H_IN = 0.55     # reserved strip for the material legend
CAPTION_H_IN = 0.42    # reserved strip under each panel for its metrics line
TITLE_H_IN = 0.30
DPI = 200


def _quiet(fn, *args, **kwargs):
    """gmsh and the loaders are chatty; the figure is the output, not the log."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        return fn(*args, **kwargs)


def _corners(mesh, i):
    """The corner nodes of element ``i`` — 4 for a quad, 3 for a triangle. Midside
    nodes are ignored: they change neither the outline nor the shape statistics."""
    n = 4 if mesh["element_types"][i] in (4, 8, 9) else 3
    return mesh["nodes"][mesh["elements"][i][:n].astype(int)]


def _metrics(mesh, target_size):
    """Element count, quadrilateral share, AR95 and delivered element size.

    Same definitions as test/quad_mesh_check.py::_metrics."""
    ars, sizes = [], []
    for i in range(len(mesh["elements"])):
        v = _corners(mesh, i)
        lengths = np.linalg.norm(np.roll(v, -1, axis=0) - v, axis=1)
        if lengths.min() <= 0:
            continue
        ars.append(lengths.max() / lengths.min())
        x, y = v[:, 0], v[:, 1]
        area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
        sizes.append(np.sqrt(area))
    ars = np.asarray(ars)
    types = np.asarray(mesh["element_types"])
    return {
        "n": len(mesh["elements"]),
        "quad_pct": float(100.0 * np.mean(np.isin(types, (4, 8, 9)))),
        "ar95": float(np.percentile(ars, 95)),
        "med_over_req": float(np.median(sizes) / target_size),
    }


def _caption(m):
    return (f"{m['n']} elements   |   {m['quad_pct']:.1f}% quadrilateral   |   "
            f"AR95 {m['ar95']:.2f}   |   median element "
            f"{m['med_over_req']:.2f} x requested")


def _draw(ax, mesh, names):
    """Material-colored element patches with hairline edges, equal aspect."""
    nodes = mesh["nodes"]
    mats = np.asarray(mesh["element_materials"])
    handles = []
    for mid in sorted(set(mats.tolist())):
        polys = [_corners(mesh, i) for i in np.flatnonzero(mats == mid)]
        # Mesh material IDs are 1-based (gmsh); the palette is 0-based, so the
        # zone colors here match the Inputs view and every other xslope plot.
        color = get_material_color(int(mid) - 1)
        ax.add_collection(PolyCollection(polys, facecolors=[color], alpha=0.55,
                                         edgecolors="0.25", linewidths=0.25))
        label = names[mid - 1] if 0 <= mid - 1 < len(names) else f"Material {mid}"
        handles.append(Patch(facecolor=color, alpha=0.55, edgecolor="0.25",
                             label=label))
    x0, x1 = nodes[:, 0].min(), nodes[:, 0].max()
    y0, y1 = nodes[:, 1].min(), nodes[:, 1].max()
    px, py = PAD_FRAC * (x1 - x0), PAD_FRAC * (y1 - y0)
    ax.set_xlim(x0 - px, x1 + px)
    ax.set_ylim(y0 - py, y1 + py)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color("0.6")
    return handles, (x1 - x0 + 2 * px), (y1 - y0 + 2 * py)


def main(repo_root=None):
    root = repo_root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data = _quiet(load_slope_data, os.path.join(root, MODEL))
    polygons = _quiet(get_material_polygons, data)
    gx = [x for x, _ in data["ground_surface"].coords]
    target = (max(gx) - min(gx)) / DIVISIONS
    names = [m.get("name") or f"Material {i + 1}"
             for i, m in enumerate(data.get("materials") or [])]

    meshes = [(label, _quiet(build_mesh_from_polygons, polygons,
                             target_size=target, element_type="quad4",
                             quad_style=key))
              for key, label in STYLES]

    # Equal-aspect panels at one common scale: the axes width is the figure width
    # less the margins, and the panel height follows the domain's own aspect.
    left, right = 0.03, 0.99
    axes_w_in = FIG_W_IN * (right - left)
    # Both styles mesh the same domain, so one panel's extent sets the height of
    # every panel and the two are drawn at the same scale.
    nodes = meshes[0][1]["nodes"]
    span_x = (nodes[:, 0].max() - nodes[:, 0].min()) * (1 + 2 * PAD_FRAC)
    span_y = (nodes[:, 1].max() - nodes[:, 1].min()) * (1 + 2 * PAD_FRAC)
    panel_h_in = axes_w_in * span_y / span_x
    row_h_in = panel_h_in + CAPTION_H_IN + TITLE_H_IN
    fig_h_in = len(meshes) * row_h_in + LEGEND_H_IN + 0.15
    fig = plt.figure(figsize=(FIG_W_IN, fig_h_in), dpi=DPI)

    handles = []
    for i, (label, mesh) in enumerate(meshes):
        # Rows run top to bottom; the legend strip sits at the very bottom.
        bottom = (LEGEND_H_IN + (len(meshes) - 1 - i) * row_h_in
                  + CAPTION_H_IN) / fig_h_in
        ax = fig.add_axes((left, bottom, right - left, panel_h_in / fig_h_in))
        handles, _, _ = _draw(ax, mesh, names)
        ax.set_title(label, fontsize=11, fontweight="bold", pad=6)
        fig.text(left + (right - left) / 2,
                 bottom - CAPTION_H_IN / fig_h_in * 0.62,
                 _caption(_metrics(mesh, target)),
                 ha="center", va="center", fontsize=8.5, color="0.25")

    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.10 / fig_h_in))

    out = os.path.join(root, OUT)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=DPI, facecolor="white")
    plt.close(fig)
    print(f"wrote {OUT}  (target element size {target:.3f})")
    for label, mesh in meshes:
        m = _metrics(mesh, target)
        print(f"  {label:30s} {m['n']:5d} elements  {m['quad_pct']:5.1f}% quad  "
              f"AR95 {m['ar95']:.2f}  size {m['med_over_req']:.2f}x")


if __name__ == "__main__":
    main()
