"""Render the figures for the transient reservoir-drawdown seepage samples
(docs/seep/samples.md).

One deterministic script serves both transient samples — each is solved once
(serial) and rendered into the figures its page carries, in docs/seep/images/:

  <sample>_flownet.png   The time-stamped transient series: a handful of panels
                         stacked vertically (full pool, early / mid / late
                         drawdown, and recovery frames toward quasi-equilibrium),
                         each the full field on identical framing (equal aspect,
                         uniform cushion) with a 't = X day' title. Every panel
                         shows the material zone fills (shell / core / foundation)
                         under filled head contours + the phreatic surface, so the
                         core's retained head pocket is visible against the drained
                         shells, PLUS the instantaneous reservoir/tailwater water
                         level for that frame -- the pool drops through the series
                         alongside the lagging phreatic (the pedagogical win). Flow
                         lines are OFF on every panel: a flow net requires divergence-
                         free through-flow, which a transient storage-release state is
                         not, so no stream function exists and equal-drop flow channels
                         have no meaning -- the caveat on the transient page. Velocity
                         vectors read a frame's instantaneous flow direction in Studio.
  <sample>_history.png   NOT built here. The time history -- reservoir level,
                         phreatic elevation and exit point above, boundary
                         inflow and outflow below -- is drawn by
                         xslope.plot_seep.plot_transient_history, and the
                         shipped figure is rebuilt by its own generator beside
                         the image, docs/seep/images/tseep_history.py. One
                         producer per figure.
  <sample>_inputs.png    Drawn only for a sample whose page carries one. The
                         plot_inputs() view of the model: geometry, material
                         zones, and the seep boundary conditions with the v18
                         reservoir symbology -- the submerged-only reservoir face is
                         drawn distinctly from a fixed-head boundary, and its
                         schedule shows as two waterlines (the full-pool level at
                         t = 0 and the drawn-down level at t = end), each with the
                         standard apex-down water symbol.

The flow-net panels are rendered one-per-figure through the unchanged
``plot_seep_solution`` (so they inherit the repo frame conventions) and stacked
losslessly with Pillow. Deterministic: same mesh, same solve, same frames.

Run from the repo root:
  PYTHONPATH=. python3 tools/make_earth_dam_tseep_figures.py            # both
  PYTHONPATH=. python3 tools/make_earth_dam_tseep_figures.py johnson    # one
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

from xslope.fileio import load_slope_data                       # noqa: E402
from xslope.mesh import get_material_polygons, build_mesh_from_polygons  # noqa: E402
from xslope.seep import (build_seep_data, build_tseep_data,     # noqa: E402
                         run_transient_seepage, transient_frame_index,
                         _transient_frame_solution)
from xslope.plot_seep import plot_seep_solution                 # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO_ROOT, "docs", "seep", "files")
IMG = os.path.join(REPO_ROOT, "docs", "seep", "images")


# --------------------------------------------------------------------------- #
#  Per-sample rendering configuration
# --------------------------------------------------------------------------- #
# target       gmsh characteristic length -> a coarse-but-core-resolving tri3 mesh.
# panels       (time, label) pairs to stack in the flow-net series figure.
# panel_size   per-panel figsize (in), sized to the domain aspect.
SAMPLES = {
    "earth_dam": dict(
        xlsx="xslope_earth_dam_tseep.xlsx",
        stem="earth_dam_tseep",
        target=110.0 / 64.0,
        panel_size=(8.6, 2.55),
        panels=[(0, "full pool"), (15, "early drawdown"), (30, "mid drawdown"),
                (47, "end of drawdown (max lag)"), (120, "recovery"),
                (360, "quasi-equilibrium")],
        # SEEP-3 builds this model input by input, so the sample entry carries
        # the solved series alone and no inputs panel of its own.
        inputs=False,
    ),
    "johnson": dict(
        xlsx="xslope_johnson_res_tseep.xlsx",
        stem="johnson_res_tseep",
        target=10.0,
        panel_size=(8.6, 3.05),
        panels=[(0, "full pool"), (35, "mid drawdown"),
                (50, "end of drawdown (max lag)"), (150, "early recovery"),
                (400, "recovery"), (1000, "quasi-equilibrium")],
    ),
}


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _solve(cfg):
    sd = load_slope_data(os.path.join(FILES, cfg["xlsx"]))
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(sd),
                  cfg["target"], "tri3")
    seep = build_seep_data(mesh, sd)
    sol = _quiet(run_transient_seepage, seep, build_tseep_data(sd), verbose=False)
    return sd, seep, sol


def _panel_png(cfg, seep, frsol, t, label, vmin, vmax):
    """Render one transient frame to a fixed-size PNG (bytes) via plot_seep_solution.

    Flow lines are OFF on every panel: a flow net requires divergence-free through-
    flow, which a transient storage-release state is not (its flow is sourced from
    released storage), so a stream function — and the equal-drop flow channels it
    would draw — does not exist for these frames. Even the full-pool t = 0 frame is
    treated like the rest for series-wide consistency. Each panel instead shows the
    material zones under filled head contours, the phreatic surface, and the
    instantaneous reservoir/tailwater WATER LEVEL for that frame (show_bc_levels) —
    so the pool visibly drops through the series alongside the lagging phreatic."""
    fig = plt.figure(figsize=cfg["panel_size"])
    _quiet(plot_seep_solution, seep, frsol, fig=fig, levels=12,
           phreatic=True, flowlines=False, vectors=False,
           mesh=False, pad_frac=0.04, cmap="Spectral_r",
           vmin=vmin, vmax=vmax,
           show_title=False, show_legend=False, show_bc_levels=True)
    title = f"t = {t:g} day" + (f"   — {label}" if label else "")
    fig.axes[0].set_title(title, fontsize=12, pad=5)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=145)          # no tight bbox -> identical width
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def fig_flownet(cfg, sd, seep, sol):
    uncon = bool(sol.get("unconfined", False))
    # Shared color scale across the whole series: the series-wide head range (which the
    # full-pool t = 0 frame dominates). Passing it to every panel pins all colorbars to
    # one range, so the drawdown reads as one continuous story — late frames fade toward
    # uniform instead of re-normalizing into a per-frame bullseye.
    all_head = np.concatenate([np.asarray(fr["head"], float) for fr in sol["frames"]])
    vmin, vmax = float(all_head.min()), float(all_head.max())
    panels = []
    for t, label in cfg["panels"]:
        fr = sol["frames"][transient_frame_index(sol, float(t))]
        frsol = _transient_frame_solution(seep, fr["head"], fr["u"], fr.get("phi"),
                                          fr.get("inflow"), fr.get("outflow"),
                                          uncon, time=fr["time"])
        panels.append(_panel_png(cfg, seep, frsol, fr["time"], label, vmin, vmax))
    w = min(p.width for p in panels)
    panels = [p if p.width == w else p.resize((w, round(p.height * w / p.width)))
              for p in panels]
    out = Image.new("RGB", (w, sum(p.height for p in panels)), "white")
    y = 0
    for p in panels:
        out.paste(p, (0, y))
        y += p.height
    path = os.path.join(IMG, f"{cfg['stem']}_flownet.png")
    out.save(path)
    print("wrote", os.path.relpath(path, REPO_ROOT))


def fig_inputs(cfg, sd):
    """Render the plot_inputs() seep view to <stem>_inputs.png.

    Showcases the v18 BC symbology: the submerged-only ``reservoir`` face is drawn
    distinctly from a fixed-head boundary, and because its level is a tseep series
    the reservoir surface is shown as two waterlines -- the full-pool level at
    t = 0 and the drawn-down level at t = end -- each with the standard apex-down
    water symbol. Geometry only (no mesh), so the boundary conditions read clearly.
    """
    from xslope.plot import plot_inputs
    xs = [x for x, _ in sd["ground_surface"].coords]
    ys = [y for _, y in sd["ground_surface"].coords]
    aspect = (max(ys) - min(ys)) / (max(xs) - min(xs))
    fig = plt.figure(figsize=(10.0, max(3.2, 10.0 * aspect + 2.0)))
    _quiet(plot_inputs, sd, fig=fig, mode="seep", show_title=False,
           show_legend=True, frame="content", pad_frac=0.04)
    path = os.path.join(IMG, f"{cfg['stem']}_inputs.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(path, REPO_ROOT))


def main(argv):
    names = argv[1:] or list(SAMPLES)
    for name in names:
        cfg = SAMPLES[name]
        print(f"solving transient sample '{name}' …")
        sd, seep, sol = _solve(cfg)
        print(f"  {len(seep['nodes'])} nodes, {len(sol['frames'])} saved frames, "
              f"converged={sol['converged']}")
        if cfg.get("inputs", True):
            fig_inputs(cfg, sd)
        fig_flownet(cfg, sd, seep, sol)
    print("done")


if __name__ == "__main__":
    main(sys.argv)
