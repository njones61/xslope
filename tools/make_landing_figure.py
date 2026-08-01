"""Render the landing-page figure: three analysis modes from one model.

One deterministic run of the committed Johnson Reservoir sample
(``docs/seep/files/xslope_johnson_res.xlsx``, with its companion mesh and
seepage pore pressures) produces four panels, all drawn by the package's own
plotting functions, and composes them into a single image:

  docs/images/landing_three_modes.png

    A   The problem as it is defined — geometry, profile lines, reservoir and
        the water load it exerts, the trial circle, over the mesh built from
        that same geometry (``plot_inputs``).
    B   Finite element seepage: total-head contours and the phreatic surface
        (``plot_seep_solution``), with the computed discharge in the title.
    C   Limit equilibrium: the critical circle found by the automated search,
        with its slices (``plot_solution``), titled with Spencer's factor of
        safety.
    D   Finite element strength reduction: the viscoplastic shear strain of the
        captured at-failure field (``plot_fem_results``), titled with the SSRM
        factor of safety.

Three arrows run from panel A to B, C and D: one problem definition, three
analyses, no re-meshing and no manual transfer of the pore pressure field.

Every number printed on the figure is computed by the run, never hard-coded, so
the figure cannot drift from what the package does. The run is deterministic:
the committed mesh and seepage companion files fix the discretization, the
circular search is a deterministic grid refinement, and the SSRM bisection is
seeded from the file's own bracket. It takes about a minute.

Run from the repo root:

    PYTHONPATH=. python3 tools/make_landing_figure.py
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from xslope.fileio import load_slope_data                     # noqa: E402
from xslope.seep import build_seep_data, run_seepage_analysis  # noqa: E402
from xslope.slice import generate_slices                       # noqa: E402
from xslope.solve import spencer                               # noqa: E402
from xslope.search import circular_search                      # noqa: E402
from xslope.fem import build_fem_data, solve_ssrm              # noqa: E402
from xslope.plot import plot_inputs, plot_solution             # noqa: E402
from xslope.plot_seep import plot_seep_solution                # noqa: E402
from xslope.plot_fem import plot_fem_results                   # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = os.path.join(REPO, "docs", "seep", "files", "xslope_johnson_res.xlsx")
OUT = os.path.join(REPO, "docs", "images", "landing_three_modes.png")

METHOD = "spencer"
PANEL_FIGSIZE = (8.6, 3.15)     # every panel rendered at one size and dpi …
PANEL_DPI = 150                 # … so the four images compose on a common scale

C_INK = "#1f2933"
C_MUTED = "#5a6673"
C_ARROW = "#1f6fb2"


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _panel(draw):
    """Run one plotting call on a fresh figure and return it as a PIL image."""
    fig = plt.figure(figsize=PANEL_FIGSIZE)
    _quiet(draw, fig)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=PANEL_DPI, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


# ---------------------------------------------------------------- the analyses

def run_all():
    """One load of the sample file, three analyses, four panels."""
    t0 = time.time()
    slope_data = load_slope_data(XLSX)
    mesh = slope_data["mesh"]
    print(f"  model: {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements")

    # A — the problem definition, over the mesh built from the same geometry
    panel_a = _panel(lambda fig: plot_inputs(
        slope_data, fig=fig, mode="lem", mat_table=False,
        show_title=False, show_legend=True, legend_ncol=6))

    # B — finite element seepage
    seep_data = build_seep_data(mesh, slope_data)
    seep_solution = _quiet(run_seepage_analysis, seep_data)
    flowrate = float(seep_solution["flowrate"])
    print(f"  seepage: q = {flowrate:.3f} (converged={seep_solution['converged']})"
          f"  [{time.time() - t0:.0f}s]")
    panel_b = _panel(lambda fig: plot_seep_solution(
        seep_data, seep_solution, fig=fig, levels=14, phreatic=True,
        flowlines=True, base_mat=3, mesh=False, pad_frac=0.04,
        show_title=False, show_legend=False))  # base_mat 3 = foundation; the
        # default (1, the k=1.0 shell) sizes the flow net to ~zero channels

    # C — limit equilibrium on the critical circle from the automated search
    fs_cache, converged, _path, _cache = _quiet(circular_search, slope_data, METHOD)
    best = fs_cache[0]
    circle = {k: float(best[k]) for k in ("Xo", "Yo", "Depth")}
    circle["R"] = circle["Yo"] - circle["Depth"]   # Depth is an elevation
    ok, res = generate_slices(slope_data, circle=circle, num_slices=40)
    assert ok, res
    slice_df, failure_surface = res
    ok, lem = spencer(slice_df)
    assert ok, lem
    fs_lem = float(lem["FS"])
    print(f"  LEM: {METHOD} FS = {fs_lem:.3f} on the critical circle "
          f"(search converged={converged})  [{time.time() - t0:.0f}s]")
    def _draw_lem(fig):
        plot_solution(slope_data, slice_df, failure_surface, lem, fig=fig,
                      seep_contours=False, show_title=False, show_legend=False)
        # The base-stress arrows push the auto y-range well below the model;
        # crop back to the section so this panel is drawn at the same scale as
        # the other three.
        fig.axes[0].set_ylim(-18, 205)

    panel_c = _panel(_draw_lem)

    # D — finite element strength reduction on the same mesh and pore pressures
    # These solver settings are the ones the sample's own locked test tag uses
    # (docs/seep/seep_slope.md, type=fem_ssrm), so the factor of safety on the
    # figure is the same number the regression suite holds this file to.
    fem_data = _quiet(build_fem_data, slope_data, mesh)
    ssrm = _quiet(solve_ssrm, fem_data, F_min=1.0, F_max=1.6, tolerance=0.01,
                  max_iterations=16000, debug_level=0)
    fs_ssrm = float(ssrm["FS"])
    print(f"  SSRM: FS = {fs_ssrm:.3f} ({ssrm['failure_criterion']})"
          f"  [{time.time() - t0:.0f}s]")
    panel_d = _panel(lambda fig: plot_fem_results(
        fem_data, ssrm["last_solution"], plot_type="shear_strain",
        fs=fs_ssrm, failure_solution=ssrm.get("failure_solution"),
        fig=fig, show_title=False))

    unit = "ft" if slope_data.get("unit_system", "imperial") == "imperial" else "m"
    titles = [
        "One problem definition — geometry, materials, water",
        f"1 · Finite element seepage — q = {flowrate:.2f} {unit}³/day per {unit}",
        f"2 · Limit equilibrium — Spencer, FS = {fs_lem:.2f}",
        f"3 · FEM strength reduction — FS = {fs_ssrm:.2f}",
    ]
    subtitles = [
        "profile lines, reservoir and its water load, trial circle, over the mesh",
        "flow net — head contours, flow lines, and the phreatic surface",
        "critical circle from the automated search, with its slices",
        "shear strain at failure — the mechanism emerges, unassumed",
    ]
    return [panel_a, panel_b, panel_c, panel_d], titles, subtitles

# ------------------------------------------------------------- the composition

def compose(panels, titles, subtitles, out_path):
    """Lay the four panels out 2x2, title each, and draw the three arrows."""
    fig_w = 16.0                      # inches; the page renders it at 1200 px
    margin = 0.014                    # left/right breathing room
    gap_x = 0.072                     # the gutter the arrows live in
    col_w = (1.0 - 2 * margin - gap_x) / 2
    title_h_in = 0.62                 # space reserved above each panel
    row_gap_in = 0.30
    top_in = 0.12
    bottom_in = 0.78                  # two-line footer

    # Every panel is rendered at one figsize/dpi, but tight bounding boxes and
    # colorbars leave the widths slightly unequal. Scale each to the column
    # width and let its own aspect set its height, then give each row the
    # height of its taller panel.
    heights_in = [col_w * fig_w * p.height / p.width for p in panels]
    row_h = [max(heights_in[0], heights_in[1]), max(heights_in[2], heights_in[3])]
    fig_h = (top_in + sum(row_h) + 2 * title_h_in + row_gap_in + bottom_in)

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    boxes = []
    y_cursor_in = fig_h - top_in
    for row in (0, 1):
        y_cursor_in -= title_h_in
        title_y = y_cursor_in / fig_h
        band_top_in = y_cursor_in
        for col in (0, 1):
            i = 2 * row + col
            h_in = heights_in[i]
            x0 = margin if col == 0 else margin + col_w + gap_x
            y0 = (band_top_in - h_in) / fig_h
            ax = fig.add_axes([x0, y0, col_w, h_in / fig_h])
            ax.imshow(panels[i], aspect="auto", interpolation="lanczos")
            ax.set_xticks([]); ax.set_yticks([])
            for side in ax.spines.values():
                side.set_color("#c9d2da")
                side.set_linewidth(0.8)
            fig.text(x0, title_y + 0.30 / fig_h, titles[i], ha="left", va="bottom",
                     fontsize=15, color=C_INK, fontweight="bold")
            fig.text(x0, title_y + 0.07 / fig_h, subtitles[i], ha="left", va="bottom",
                     fontsize=11.5, color=C_MUTED)
            boxes.append((x0, y0, col_w, h_in / fig_h))
        y_cursor_in = band_top_in - row_h[row] - row_gap_in

    ax0, ay0, aw, ah = boxes[0]

    def arrow(p0, p1, rad=0.0):
        fig.patches.append(FancyArrowPatch(
            p0, p1, transform=fig.transFigure, figure=fig,
            arrowstyle="-|>", mutation_scale=24, linewidth=2.2,
            color=C_ARROW, connectionstyle=f"arc3,rad={rad}",
            shrinkA=2, shrinkB=2, zorder=5))

    # One arrow from the problem definition into each of the three analyses.
    arrow((ax0 + aw + 0.005, ay0 + ah / 2), (boxes[1][0] - 0.005, ay0 + ah / 2))
    x_down = ax0 + aw * 0.82          # clear of the titles it runs past
    arrow((x_down, ay0 - 0.004), (x_down, boxes[2][1] + boxes[2][3] + 0.004))
    arrow((ax0 + aw + 0.005, ay0 + ah * 0.18),
          (boxes[3][0] - 0.005, boxes[3][1] + boxes[3][3] * 0.62), rad=-0.22)

    fig.text(margin, 0.014,
             "Johnson Reservoir dam: seepage, limit equilibrium and finite element strength "
             "reduction from a single XSLOPE input file, on one mesh,\nwith the computed pore "
             "pressures passed straight to both stability analyses.",
             ha="left", va="bottom", fontsize=12, color=C_MUTED, linespacing=1.35)

    fig.savefig(out_path, dpi=150, facecolor="white")
    plt.close(fig)
    print("wrote", os.path.relpath(out_path, REPO),
          Image.open(out_path).size)


def main():
    print("rendering the landing figure from", os.path.relpath(XLSX, REPO))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    panels, titles, subtitles = run_all()
    compose(panels, titles, subtitles, OUT)


if __name__ == "__main__":
    main()
