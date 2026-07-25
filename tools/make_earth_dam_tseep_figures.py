"""Render the figures for the transient reservoir-drawdown seepage samples
(docs/seep/samples.md).

One deterministic script serves both transient samples — each is solved once
(serial) and rendered into two figures in docs/seep/images/:

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
  <sample>_history.png   Two stacked history plots: (top) reservoir level, the
                         phreatic elevation at an interior station, and the
                         upstream-face exit point vs time -- the phreatic lag and
                         exit-point migration; (bottom) boundary inflow vs outflow
                         -- inflow falls to zero as the face drains while outflow
                         spikes on the released storage and decays toward the
                         drained steady state. Inflow is blue, outflow a
                         contrasting dark red so the two are unmistakable.

                         The "exit point on upstream face" trace reports the top of
                         the upstream SEEPAGE face -- the highest still-draining
                         (saturated) point ABOVE the current pool. Once the pool
                         stabilises and the interior drains, no exposed face node
                         seeps, so the trace CLAMPS to the pool waterline: below the
                         line the face is submerged (held at reservoir head) and
                         nothing exits, so no exit point exists there. The interior
                         phreatic station can legitimately settle a little BELOW the
                         drawn-down pool -- at quasi-equilibrium the water table is a
                         through-flow surface sloping from the low pool toward the
                         downstream exit -- and the panel annotates that where it
                         happens so the trace does not read as an error.
  <sample>_inputs.png    The plot_inputs() view of the model: geometry, material
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

    def _exit_point(mask, head, pool_level):
        """Top of the upstream SEEPAGE face: the highest exposed (above-pool)
        face node that is still saturated (pressure head >= 0), i.e. still
        draining. Once the reservoir stabilises and the interior drains, no
        exposed face node seeps, so the diagnostic CLAMPS to the pool waterline
        -- below the line the face is submerged (held at reservoir head) and
        nothing exits, so reporting a point there would be an artifact."""
        idx = np.where(mask)[0]
        exposed = idx[nodes[idx, 1] > pool_level + 1e-6]
        seeping = exposed[head[exposed] - nodes[exposed, 1] >= -0.1]
        if seeping.size:
            return max(pool_level, float(np.max(nodes[seeping, 1])))
        return pool_level

    col = np.abs(nodes[:, 0] - station_x) < cfg["station_halfwidth"]
    face = (np.abs(nodes[:, 0] - cfg["face"](nodes[:, 1])) < cfg["face_tol"]) & \
           (nodes[:, 1] >= lo) & (nodes[:, 1] <= hi)
    phr = np.array([_saturated_top(col, f["head"]) for f in frames])
    exitpt = np.array([_exit_point(face, f["head"], pool[i])
                       for i, f in enumerate(frames)])

    d0, d1 = cfg["drawdown"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.4, 6.4), sharex=True)

    ax1.plot(t, pool, "-", color="#2b7bb0", lw=2.0, marker="o", ms=4,
             label="reservoir level $h(t)$")
    ax1.plot(t, phr, "--", color="#1f4e79", lw=1.8, marker="s", ms=4,
             label=f"phreatic elev. at x = {station_x:g} ft")
    ax1.plot(t, exitpt, ":", color="#e08214", lw=1.8, marker="^", ms=5,
             label="exit point (top of seepage face)")
    # The interior phreatic station can settle a little BELOW the drawn-down pool:
    # at quasi-equilibrium the water table is a through-flow surface sloping from
    # the (low) pool toward the downstream exit, so an interior column drains below
    # the upstream level. Annotate it where it happens so the trace is not misread.
    if np.isfinite(phr[-1]) and phr[-1] < pool[-1] - 0.15:
        ax1.annotate("interior surface drains\ntoward the downstream exit",
                     xy=(t[-1], phr[-1]), xytext=(0.60, 0.26),
                     textcoords="axes fraction", fontsize=8.0, color="#1f4e79",
                     ha="left", va="center",
                     arrowprops=dict(arrowstyle="-|>", color="#1f4e79", lw=0.9))
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
    ax2.plot(t, outflow, "-", color="#c0392b", lw=2.0, marker="D", ms=4,
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
        fig_inputs(cfg, sd)
        fig_flownet(cfg, sd, seep, sol)
        fig_history(cfg, sd, seep, sol)
    print("done")


if __name__ == "__main__":
    main(sys.argv)
