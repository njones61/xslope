"""Render the figures for the rapid-drawdown-from-a-transient-solution worked
example (docs/lem/rapid.md, "Rapid Drawdown from a Transient Solution").

One deterministic script solves the Johnson Reservoir rapid-drawdown workbook
(``docs/lem/files/xslope_johnson_res_rapid.xlsx``) end to end and renders two
figures into ``docs/lem/rapid_images/``:

  johnson_res_rapid_stages.png   The two staged pore-pressure fields the three-stage
                                 method consumes, stacked: **Stage 1** (full pool,
                                 the transient frame at ``stage_1`` = t = 0 day) and
                                 **Stage 2** (end of drawdown, the frame at
                                 ``stage_2`` = t = 50 day). Each panel is the total-head
                                 field with material zones and the phreatic surface, and
                                 the critical upstream slip circle drawn over it. Between
                                 the panels the reservoir has been drawn all the way down,
                                 but the low-permeability core still holds an elevated
                                 head pocket — the retained pore pressure that governs
                                 Stage 2.
  johnson_res_rapid_solution.png The limit-equilibrium solution: the slope, the material
                                 zones, the critical upstream circle and its slices, and
                                 the full-pool head contours, titled with the rapid-
                                 drawdown factor of safety.

The whole pipeline runs in memory — build the mesh, march the transient seepage
solution, stage the two frames with ``stage_transient_for_drawdown`` (no seep.csv /
seep2.csv companions), then run the three-stage ``rapid_drawdown`` on the critical
circle. Deterministic: same mesh, same solve, same frames.

Run from the repo root:
  PYTHONPATH=. python3 tools/make_johnson_rapid_figures.py
"""

from __future__ import annotations

import contextlib
import io
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from xslope.fileio import load_slope_data                              # noqa: E402
from xslope.mesh import build_mesh_from_polygons, get_material_polygons  # noqa: E402
from xslope.seep import (build_seep_data, build_tseep_data,            # noqa: E402
                         run_transient_seepage, stage_transient_for_drawdown,
                         transient_frame_index, _transient_frame_solution)
from xslope.slice import generate_slices                              # noqa: E402
from xslope.advanced import rapid_drawdown                            # noqa: E402
from xslope.plot import plot_solution                                 # noqa: E402
from xslope.plot_seep import plot_seep_solution                      # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
XLSX = os.path.join(REPO_ROOT, "docs", "lem", "files", "xslope_johnson_res_rapid.xlsx")
IMG = os.path.join(REPO_ROOT, "docs", "lem", "rapid_images")

METHOD = "spencer"
TARGET = 10.0          # gmsh characteristic length -> the 1024-node tri3 mesh
PANEL_SIZE = (8.6, 3.05)


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def solve():
    sd = load_slope_data(XLSX)
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(sd), TARGET, "tri3")
    seep = build_seep_data(mesh, sd)
    sol = _quiet(run_transient_seepage, seep, build_tseep_data(sd), verbose=False)
    sd = stage_transient_for_drawdown(sd, sol)
    sd["mesh"] = mesh
    ok, res = generate_slices(sd, circle=sd["circles"][0], num_slices=40, debug=False)
    assert ok, res
    slice_df, failure_surface = res
    ok, results = rapid_drawdown(slice_df.copy(), METHOD, debug_level=0)
    assert ok, results
    return sd, seep, sol, slice_df, failure_surface, results


def _stage_panel(seep, frsol, arc, title):
    """One staged pore-pressure field with the slip circle overlaid -> PNG bytes."""
    fig = plt.figure(figsize=PANEL_SIZE)
    _quiet(plot_seep_solution, seep, frsol, fig=fig, levels=14,
           phreatic=True, flowlines=False, vectors=False, mesh=False,
           pad_frac=0.04, cmap="Spectral_r", show_title=False, show_legend=False)
    ax = fig.axes[0]
    ax.plot(arc[:, 0], arc[:, 1], color="#c1121f", lw=2.2, zorder=6,
            label="critical upstream circle")
    ax.set_title(title, fontsize=12, pad=5)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=145)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def fig_stages(sd, seep, sol, failure_surface):
    uncon = bool(sol.get("unconfined", False))
    arc = np.asarray(failure_surface.coords)
    ts = sd["tseep"]
    panels = []
    for t, tag in [(ts["stage_1"], "Stage 1 — full pool"),
                   (ts["stage_2"], "Stage 2 — end of drawdown")]:
        fr = sol["frames"][transient_frame_index(sol, float(t))]
        frsol = _transient_frame_solution(seep, fr["head"], fr["u"], fr.get("phi"),
                                          fr.get("inflow"), fr.get("outflow"),
                                          uncon, time=fr["time"])
        panels.append(_stage_panel(seep, frsol, arc,
                                   f"{tag}   (t = {fr['time']:g} day)"))
    w = min(p.width for p in panels)
    panels = [p if p.width == w else p.resize((w, round(p.height * w / p.width)))
              for p in panels]
    out = Image.new("RGB", (w, sum(p.height for p in panels)), "white")
    y = 0
    for p in panels:
        out.paste(p, (0, y)); y += p.height
    path = os.path.join(IMG, "johnson_res_rapid_stages.png")
    out.save(path)
    print("wrote", os.path.relpath(path, REPO_ROOT))


def fig_solution(sd, slice_df, failure_surface, results):
    fig = plt.figure(figsize=(11, 5.2))
    _quiet(plot_solution, sd, slice_df, failure_surface, results,
           fig=fig, seep_contours=True, show_title=True, show_legend=True)
    path = os.path.join(IMG, "johnson_res_rapid_solution.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(path, REPO_ROOT))


def main():
    print("solving Johnson Reservoir rapid-drawdown worked example …")
    sd, seep, sol, slice_df, failure_surface, results = solve()
    print(f"  {len(seep['nodes'])} nodes, {len(sol['frames'])} frames, "
          f"{METHOD}: rapid FS = {results['FS']:.3f}  "
          f"(Stage 1 = {results['stage1_FS']:.3f}, Stage 2 = {results['stage2_FS']:.3f}, "
          f"Stage 3 = {results['stage3_FS']:.3f})")
    fig_stages(sd, seep, sol, failure_surface)
    fig_solution(sd, slice_df, failure_surface, results)
    print("done")


if __name__ == "__main__":
    main()
