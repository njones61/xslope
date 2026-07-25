"""Render the figures for the transient earth-dam drawdown sample
(docs/seep/samples.md, "Earth Dam — Reservoir Drawdown (Transient)").

Solves ``docs/seep/files/xslope_earth_dam_tseep.xlsx`` once (serial, ~10 s on a
coarse tri3 mesh) and renders two figures into docs/seep/images/:

  earth_dam_tseep_flownet.png   The time-stamped flow-net series: six panels
                                stacked vertically (full pool, early / mid / late
                                drawdown, and two recovery frames toward quasi-
                                equilibrium), each the full field on identical
                                framing (equal aspect, uniform cushion) with a
                                't = X day' title. Every panel shows filled head
                                contours + the phreatic surface; flow lines are
                                drawn only where a stream function exists (the
                                steady full-pool frame). Once the reservoir face
                                turns to an exit face the flow is storage-release,
                                not divergence-free through-flow, so no stream
                                function is defined and the panel is contours +
                                phreatic only -- the instantaneous-flow-net caveat
                                on the transient page.
  earth_dam_tseep_history.png   Two stacked history plots: (top) reservoir level,
                                the phreatic elevation at an upstream station, and
                                the upstream-face exit point vs time -- the
                                phreatic lag and exit-point migration; (bottom)
                                boundary inflow vs outflow -- inflow falls to zero
                                as the face drains while outflow spikes on the
                                released storage and decays toward the drained
                                steady state.

The flow-net panels are rendered one-per-figure through the unchanged
``plot_seep_solution`` (so they inherit the repo frame conventions) and stacked
losslessly with Pillow. Deterministic: same mesh, same solve, same frames.

Run from the repo root:  PYTHONPATH=. python3 tools/make_earth_dam_tseep_figures.py
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
XLSX = os.path.join(REPO_ROOT, "docs", "seep", "files", "xslope_earth_dam_tseep.xlsx")
IMG = os.path.join(REPO_ROOT, "docs", "seep", "images")
TARGET_SIZE = 110.0 / 64.0            # coarse tri3, ~490 nodes
PANEL_TIMES = [0, 15, 30, 47, 120, 360]
PANEL_LABEL = {0: "full pool", 15: "early drawdown", 30: "mid drawdown",
               47: "end of drawdown (max lag)", 120: "recovery",
               360: "quasi-equilibrium"}
FACE_SLOPE = 42.0 / 18.0              # upstream face: x = (42/18) * y, toe (0,0)->(42,18)
STATION_X = 30.0                      # upstream monitoring station


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _solve():
    sd = load_slope_data(XLSX)
    mesh = _quiet(build_mesh_from_polygons, get_material_polygons(sd),
                  TARGET_SIZE, "tri3")
    seep = build_seep_data(mesh, sd)
    sol = _quiet(run_transient_seepage, seep, build_tseep_data(sd), verbose=False)
    return sd, seep, sol


def _pool(sd, t):
    ts = sd["tseep"]
    return float(np.interp(t, ts["times"], ts["series"]["pool"]))


def _panel_png(seep, frsol, t):
    """Render one flow-net frame to a fixed-size PNG (bytes) via plot_seep_solution."""
    phi = np.asarray(frsol.get("phi"))
    has_stream = phi.size and np.ptp(phi) > 1e-9      # flow lines only where phi exists
    fig = plt.figure(figsize=(8.6, 2.55))
    _quiet(plot_seep_solution, seep, frsol, fig=fig, levels=12,
           phreatic=True, flowlines=bool(has_stream), vectors=False,
           mesh=False, pad_frac=0.04, cmap="Spectral_r",
           show_title=False, show_legend=False)
    lab = PANEL_LABEL.get(int(round(t)), "")
    title = f"t = {t:g} day" + (f"   — {lab}" if lab else "")
    fig.axes[0].set_title(title, fontsize=12, pad=5)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=145)          # no tight bbox -> identical width
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def fig_flownet(sd, seep, sol):
    uncon = bool(sol.get("unconfined", False))
    panels = []
    for t in PANEL_TIMES:
        fr = sol["frames"][transient_frame_index(sol, float(t))]
        frsol = _transient_frame_solution(seep, fr["head"], fr["u"], fr.get("phi"),
                                          fr.get("inflow"), fr.get("outflow"),
                                          uncon, time=fr["time"])
        panels.append(_panel_png(seep, frsol, fr["time"]))
    w = min(p.width for p in panels)
    panels = [p if p.width == w else p.resize((w, round(p.height * w / p.width)))
              for p in panels]
    out = Image.new("RGB", (w, sum(p.height for p in panels)), "white")
    y = 0
    for p in panels:
        out.paste(p, (0, y))
        y += p.height
    path = os.path.join(IMG, "earth_dam_tseep_flownet.png")
    out.save(path)
    print("wrote", os.path.relpath(path, REPO_ROOT))


def fig_history(sd, seep, sol):
    nodes = seep["nodes"]
    frames = sol["frames"]
    t = np.array([f["time"] for f in frames])
    pool = np.array([_pool(sd, tt) for tt in t])
    inflow = np.array([f.get("inflow", np.nan) for f in frames])
    outflow = np.array([f.get("outflow", np.nan) for f in frames])

    def _saturated_top(mask, head):
        idx = np.where(mask)[0]
        sat = idx[head[idx] - nodes[idx, 1] >= -0.1]
        return float(np.max(nodes[sat, 1])) if sat.size else np.nan

    col = np.abs(nodes[:, 0] - STATION_X) < 2.5
    face = (np.abs(nodes[:, 0] - FACE_SLOPE * nodes[:, 1]) < 0.6) & (nodes[:, 1] <= 18.1)
    phr = np.array([_saturated_top(col, f["head"]) for f in frames])
    exitpt = np.array([_saturated_top(face, f["head"]) for f in frames])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.4, 6.4), sharex=True)

    ax1.plot(t, pool, "-", color="#2b7bb0", lw=2.0, marker="o", ms=4,
             label="reservoir level $h(t)$")
    ax1.plot(t, phr, "--", color="#1f4e79", lw=1.8, marker="s", ms=4,
             label=f"phreatic elev. at x = {STATION_X:g} ft")
    ax1.plot(t, exitpt, ":", color="#e08214", lw=1.8, marker="^", ms=5,
             label="exit point on upstream face")
    ax1.axvspan(2, 47, color="#f2e6cf", alpha=0.6, zorder=0)
    ax1.annotate("45-day\ndrawdown", xy=(24.5, 16.5), ha="center", va="top",
                 fontsize=8.5, color="#9a5c14")
    ax1.set_ylabel("elevation  (ft)")
    ax1.set_ylim(0, 20)
    ax1.set_title("Phreatic surface lags the reservoir; the exit point migrates "
                  "down the face", fontsize=10.5)
    ax1.grid(alpha=0.25)
    ax1.legend(loc="upper right", fontsize=8.5)

    ax2.plot(t, inflow, "-", color="#2b7bb0", lw=1.8, marker="o", ms=4,
             label="boundary inflow")
    ax2.plot(t, outflow, "-", color="#0f8a86", lw=2.0, marker="D", ms=4,
             label="boundary outflow")
    ax2.axvspan(2, 47, color="#f2e6cf", alpha=0.6, zorder=0)
    ax2.set_xlabel("time  (day)")
    ax2.set_ylabel("flow rate  (ft$^3$/day per ft)")
    ax2.set_title("Inflow → 0 as the face drains; outflow spikes on released "
                  "storage, then decays", fontsize=10.5)
    ax2.grid(alpha=0.25)
    ax2.legend(loc="upper right", fontsize=8.5)

    fig.tight_layout()
    path = os.path.join(IMG, "earth_dam_tseep_history.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", os.path.relpath(path, REPO_ROOT))


def main():
    print("solving transient earth-dam drawdown …")
    sd, seep, sol = _solve()
    print(f"  {len(seep['nodes'])} nodes, {len(sol['frames'])} saved frames, "
          f"converged={sol['converged']}")
    fig_flownet(sd, seep, sol)
    fig_history(sd, seep, sol)
    print("done")


if __name__ == "__main__":
    main()
