"""Render the figures for the transient reservoir-drawdown seepage samples
(docs/seep/samples.md).

One deterministic script serves both transient samples — each is solved once
(serial) and rendered into two figures in docs/seep/images/:

  <sample>_flownet.png   The time-stamped flow-net series: a handful of panels
                         stacked vertically (full pool, early / mid / late
                         drawdown, and recovery frames toward quasi-equilibrium),
                         each the full field on identical framing (equal aspect,
                         uniform cushion) with a 't = X day' title. Every panel
                         shows the material zone fills (shell / core / foundation)
                         under filled head contours + the phreatic surface, so the
                         core's retained head pocket is visible against the drained
                         shells. Flow lines are drawn only where a stream function
                         exists (the steady full-pool frame). Once the reservoir
                         face turns to an exit face the flow is storage-release, not
                         divergence-free through-flow, so no stream function is
                         defined and the panel is zones + contours + phreatic only
                         -- the instantaneous-flow-net caveat on the transient page.
  <sample>_history.png   Two stacked history plots: (top) reservoir level, the
                         phreatic elevation at an upstream station, and the
                         upstream-face exit point vs time -- the phreatic lag and
                         exit-point migration; (bottom) boundary inflow vs outflow
                         -- inflow falls to zero as the face drains while outflow
                         spikes on the released storage and decays toward the
                         drained steady state.

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
# station_x    upstream monitoring column for the phreatic-lag history trace.
# face         a callable y -> x mapping the submerged upstream face (exit-point
#              track), plus the (y_lo, y_hi) band over which it applies.
# elev_lim     y-limits for the history elevation axis.
SAMPLES = {
    "earth_dam": dict(
        xlsx="xslope_earth_dam_tseep.xlsx",
        stem="earth_dam_tseep",
        target=110.0 / 64.0,
        panel_size=(8.6, 2.55),
        panels=[(0, "full pool"), (15, "early drawdown"), (30, "mid drawdown"),
                (47, "end of drawdown (max lag)"), (120, "recovery"),
                (360, "quasi-equilibrium")],
        station_x=30.0,
        station_halfwidth=2.5,
        face=lambda y: (42.0 / 18.0) * y,     # toe (0,0) -> crest heel (42,18)
        face_band=(0.0, 18.1),
        face_tol=0.6,
        elev_lim=(0, 20),
        drawdown=(2, 47),
    ),
    "johnson": dict(
        xlsx="xslope_johnson_res_tseep.xlsx",
        stem="johnson_res_tseep",
        target=10.0,
        panel_size=(8.6, 3.05),
        panels=[(0, "full pool"), (35, "mid drawdown"),
                (50, "end of drawdown (max lag)"), (150, "early recovery"),
                (400, "recovery"), (1000, "quasi-equilibrium")],
        station_x=280.0,
        station_halfwidth=15.0,
        # submerged upstream face: heel (200,100) -> pool line (320,160), slope 0.5
        face=lambda y: 200.0 + 2.0 * (y - 100.0),
        face_band=(100.0, 160.1),
        face_tol=3.0,
        elev_lim=(95, 170),
        drawdown=(5, 50),
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


def _pool(sd, t):
    ts = sd["tseep"]
    return float(np.interp(t, ts["times"], ts["series"]["pool"]))


def _panel_png(cfg, seep, frsol, t, label):
    """Render one flow-net frame to a fixed-size PNG (bytes) via plot_seep_solution."""
    phi = np.asarray(frsol.get("phi"))
    has_stream = phi.size and np.ptp(phi) > 1e-9      # flow lines only where phi exists
    fig = plt.figure(figsize=cfg["panel_size"])
    _quiet(plot_seep_solution, seep, frsol, fig=fig, levels=12,
           phreatic=True, flowlines=bool(has_stream), vectors=False,
           mesh=False, pad_frac=0.04, cmap="Spectral_r",
           show_title=False, show_legend=False)
    title = f"t = {t:g} day" + (f"   — {label}" if label else "")
    fig.axes[0].set_title(title, fontsize=12, pad=5)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=145)          # no tight bbox -> identical width
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).convert("RGB")


def fig_flownet(cfg, sd, seep, sol):
    uncon = bool(sol.get("unconfined", False))
    panels = []
    for t, label in cfg["panels"]:
        fr = sol["frames"][transient_frame_index(sol, float(t))]
        frsol = _transient_frame_solution(seep, fr["head"], fr["u"], fr.get("phi"),
                                          fr.get("inflow"), fr.get("outflow"),
                                          uncon, time=fr["time"])
        panels.append(_panel_png(cfg, seep, frsol, fr["time"], label))
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


def fig_history(cfg, sd, seep, sol):
    nodes = seep["nodes"]
    frames = sol["frames"]
    t = np.array([f["time"] for f in frames])
    pool = np.array([_pool(sd, tt) for tt in t])
    inflow = np.array([f.get("inflow", np.nan) for f in frames])
    outflow = np.array([f.get("outflow", np.nan) for f in frames])
    station_x = cfg["station_x"]
    lo, hi = cfg["face_band"]

    def _saturated_top(mask, head):
        idx = np.where(mask)[0]
        sat = idx[head[idx] - nodes[idx, 1] >= -0.1]
        return float(np.max(nodes[sat, 1])) if sat.size else np.nan

    col = np.abs(nodes[:, 0] - station_x) < cfg["station_halfwidth"]
    face = (np.abs(nodes[:, 0] - cfg["face"](nodes[:, 1])) < cfg["face_tol"]) & \
           (nodes[:, 1] >= lo) & (nodes[:, 1] <= hi)
    phr = np.array([_saturated_top(col, f["head"]) for f in frames])
    exitpt = np.array([_saturated_top(face, f["head"]) for f in frames])

    d0, d1 = cfg["drawdown"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.4, 6.4), sharex=True)

    ax1.plot(t, pool, "-", color="#2b7bb0", lw=2.0, marker="o", ms=4,
             label="reservoir level $h(t)$")
    ax1.plot(t, phr, "--", color="#1f4e79", lw=1.8, marker="s", ms=4,
             label=f"phreatic elev. at x = {station_x:g} ft")
    ax1.plot(t, exitpt, ":", color="#e08214", lw=1.8, marker="^", ms=5,
             label="exit point on upstream face")
    ax1.axvspan(d0, d1, color="#f2e6cf", alpha=0.6, zorder=0)
    ax1.annotate("45-day\ndrawdown", xy=(0.5 * (d0 + d1),
                 cfg["elev_lim"][0] + 0.92 * (cfg["elev_lim"][1] - cfg["elev_lim"][0])),
                 ha="center", va="top", fontsize=8.5, color="#9a5c14")
    ax1.set_ylabel("elevation  (ft)")
    ax1.set_ylim(*cfg["elev_lim"])
    ax1.set_title("Phreatic surface lags the reservoir; the exit point migrates "
                  "down the face", fontsize=10.5)
    ax1.grid(alpha=0.25)
    ax1.legend(loc="upper right", fontsize=8.5)

    ax2.plot(t, inflow, "-", color="#2b7bb0", lw=1.8, marker="o", ms=4,
             label="boundary inflow")
    ax2.plot(t, outflow, "-", color="#0f8a86", lw=2.0, marker="D", ms=4,
             label="boundary outflow")
    ax2.axvspan(d0, d1, color="#f2e6cf", alpha=0.6, zorder=0)
    ax2.set_xlabel("time  (day)")
    ax2.set_ylabel("flow rate  (ft$^3$/day per ft)")
    ax2.set_title("Inflow → 0 as the face drains; outflow spikes on released "
                  "storage, then decays", fontsize=10.5)
    ax2.grid(alpha=0.25)
    ax2.legend(loc="upper right", fontsize=8.5)

    fig.tight_layout()
    path = os.path.join(IMG, f"{cfg['stem']}_history.png")
    fig.savefig(path, dpi=150)
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
        fig_flownet(cfg, sd, seep, sol)
        fig_history(cfg, sd, seep, sol)
    print("done")


if __name__ == "__main__":
    main(sys.argv)
