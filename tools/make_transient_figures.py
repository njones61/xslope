"""Render the schematic / analytic illustration figures for docs/seep/transient.md.

None of these figures runs the seepage solver — they are pure matplotlib schematics
and closed-form curves, so they are deterministic and cheap. Five figures land in
docs/seep/images/:

  transient_reservoir_bc.png    Two-panel drawdown schematic of the submerged-only
                                Dirichlet ("reservoir") boundary: as the reservoir
                                level falls, face nodes cross from held reservoir
                                head (submerged) to exit face (unsubmerged) and the
                                still-saturated soil drains back through the exposed
                                upper face.
  transient_isochrones.png      The dissipation picture: a sudden boundary head
                                change diffusing into a 1-D column, head-vs-depth
                                isochrones at log-spaced times (erfc solution, the
                                same closed form the phase-1 analytical lock uses).
  transient_storage_models.png  The storage coefficient S(psi) for the two families
                                as implemented in storage_capacity_vec: the
                                linear-front / Gardner boxcar (Ss floor + Sy/|h0|
                                band) and the smooth van Genuchten capacity peak.
  transient_series_semantics.png  A head time series showing the breakpoint rules:
                                defined points, linear interpolation, hold-constant
                                beyond the ends, and a repeated-time vertical step.
  transient_frame_schedule.png  A timing-diagram of the saved-frame schedule as the
                                union of the save_interval grid, explicit save_times,
                                the stage times, and every series breakpoint.

Run from the repo root:  python tools/make_transient_figures.py
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch
from matplotlib.lines import Line2D
from scipy.special import erfc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from xslope.seep import storage_capacity_vec, KR_VG  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "docs", "seep", "images")

# ---- shared palette -------------------------------------------------------
C_WATER_FILL = "#cfe7f5"
C_WATER_LINE = "#2b7bb0"
C_DAM_FILL = "#e9dfc6"
C_DAM_EDGE = "#7a6a4c"
C_GROUND = "#d9cdb0"
C_SUBMERGED = "#1f6fb2"          # held-head nodes (blue, filled)
C_EXIT = "#e08214"              # exit-face nodes (orange, open)
C_PHREATIC = "#1f4e79"
C_SEEP = "#0f8a86"              # seepage-out arrows (teal)


# ===========================================================================
# Figure 1 — submerged-only reservoir boundary (two-panel drawdown schematic)
# ===========================================================================

# Embankment geometry, shared by both panels (upstream face is the key polyline).
US_TOE = (2.0, 0.0)
US_CREST = (9.0, 7.0)            # upstream face is 45 deg: x = 2 + y
DS_CREST = (13.0, 7.0)
DS_TOE = (20.0, 0.0)
FOUND_L, FOUND_R = -5.5, 21.5
DAM_XY = [US_TOE, US_CREST, DS_CREST, DS_TOE]

# nodes marched up the upstream face (elevations); x = 2 + y on the 45-deg face
FACE_Y = np.array([0.4, 1.2, 2.0, 2.8, 3.6, 4.4, 5.2, 6.0, 6.8])


def _face_x(y):
    return US_TOE[0] + (np.asarray(y) - US_TOE[1]) * \
        (US_CREST[0] - US_TOE[0]) / (US_CREST[1] - US_TOE[1])


def _water_symbol(ax, x, y, size=0.45, z=7):
    """Standard apex-down water-surface triangle: the BOTTOM TIP touches the
    waterline (the tip marks the level), body above it, ticks below."""
    tri = np.array([[x - size, y + 1.25 * size], [x + size, y + 1.25 * size], [x, y]])
    ax.add_patch(Polygon(tri, closed=True, facecolor="white",
                         edgecolor=C_WATER_LINE, lw=1.2, zorder=z))
    ax.plot([x - 0.55 * size, x + 0.55 * size], [y - 0.40 * size] * 2,
            color=C_WATER_LINE, lw=1.0, zorder=z)
    ax.plot([x - 0.28 * size, x + 0.28 * size], [y - 0.70 * size] * 2,
            color=C_WATER_LINE, lw=1.0, zorder=z)


def _phreatic_curve(ctrl_x, ctrl_y, n=120):
    """Smooth (cubic-fit) monotone-ish phreatic line through control points."""
    c = np.polyfit(ctrl_x, ctrl_y, 3)
    xs = np.linspace(min(ctrl_x), max(ctrl_x), n)
    return xs, np.polyval(c, xs)


def _draw_panel(ax, h, phreatic_ctrl, seepage_band, title, wsym_x, hlabel):
    ax.set_aspect("equal")
    ax.set_xlim(FOUND_L, FOUND_R)
    ax.set_ylim(-1.4, 10.2)
    ax.axis("off")
    ax.set_title(title, fontsize=11.5, pad=6)

    # ground / foundation
    ax.add_patch(Polygon([[FOUND_L, 0], [FOUND_R, 0], [FOUND_R, -1.4],
                          [FOUND_L, -1.4]], closed=True,
                         facecolor=C_GROUND, edgecolor="none", zorder=0))
    ax.plot([FOUND_L, FOUND_R], [0, 0], color="#8a7c5c", lw=1.1, zorder=1)

    # reservoir water (fill up to level h against the upstream face)
    hx = float(_face_x(h))
    water = [[FOUND_L, h], [hx, h], [US_TOE[0], 0.0], [FOUND_L, 0.0]]
    ax.add_patch(Polygon(water, closed=True, facecolor=C_WATER_FILL,
                         edgecolor="none", zorder=1.5))
    ax.plot([FOUND_L, hx], [h, h], color=C_WATER_LINE, lw=1.4, zorder=4)
    _water_symbol(ax, wsym_x, h)
    ax.annotate(hlabel, xy=(FOUND_L + 0.4, h + 0.35), fontsize=10.5,
                color=C_WATER_LINE, ha="left", va="bottom")

    # embankment body
    ax.add_patch(Polygon(DAM_XY, closed=True, facecolor=C_DAM_FILL,
                         edgecolor=C_DAM_EDGE, lw=1.6, zorder=2))
    # emphasise the upstream face polyline (the face drawn in full)
    ax.plot([US_TOE[0], US_CREST[0]], [US_TOE[1], US_CREST[1]],
            color=C_DAM_EDGE, lw=3.0, solid_capstyle="round", zorder=3)

    # interior phreatic surface (dashed)
    px, py = _phreatic_curve(*phreatic_ctrl)
    ax.plot(px, py, ls=(0, (5, 3)), color=C_PHREATIC, lw=1.8, zorder=5)

    # seepage-out arrows through the exposed upper face (t2 only)
    if seepage_band is not None:
        nrm = np.array([-1.0, 1.0]) / np.sqrt(2.0)   # outward face normal (up-left)
        for yb in seepage_band:
            base = np.array([float(_face_x(yb)), yb])
            tip = base + 1.15 * nrm
            ax.add_patch(FancyArrowPatch(base, tip, arrowstyle="-|>",
                                         mutation_scale=13, lw=1.7,
                                         color=C_SEEP, zorder=6))
        ytop = max(seepage_band)
        ax.annotate("seepage out\n(soil still saturated)",
                    xy=tuple(np.array([float(_face_x(ytop)), ytop]) + 0.8 *
                             np.array([-1, 1]) / np.sqrt(2)),
                    xytext=(2.6, 6.6), fontsize=9.0, color=C_SEEP,
                    ha="center", va="bottom",
                    arrowprops=dict(arrowstyle="-", color=C_SEEP, lw=0.9))

    # nodes on the face: submerged (held head) vs unsubmerged (exit face)
    fy = FACE_Y
    fx = _face_x(fy)
    sub = fy <= h + 1e-9
    ax.plot(fx[sub], fy[sub], "o", ms=8, mfc=C_SUBMERGED, mec="white",
            mew=1.0, zorder=8)
    ax.plot(fx[~sub], fy[~sub], "o", ms=8, mfc="white", mec=C_EXIT,
            mew=1.8, zorder=8)


def fig_reservoir_bc():
    # stacked vertically so each panel spans the full docs column — side-by-side
    # rendered the text unreadably small at page width (Norm)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.6, 9.6))

    # t1 — high reservoir; phreatic enters at the waterline, no exposed seepage
    _draw_panel(
        ax1, h=5.6,
        phreatic_ctrl=([7.6, 9.0, 13.0, 18.0], [5.6, 5.5, 3.6, 1.3]),
        seepage_band=None,
        title=r"$t = t_1$   (early drawdown)",
        wsym_x=4.0, hlabel=r"$h(t_1)$")

    # t2 — lower reservoir; interior phreatic LAGS above it and drains out the face
    _draw_panel(
        ax2, h=2.4,
        phreatic_ctrl=([6.6, 9.0, 13.0, 18.0], [4.6, 4.9, 3.5, 1.1]),
        seepage_band=[2.9, 3.6, 4.3],
        title=r"$t = t_2$   (later — reservoir drawn down)",
        wsym_x=-0.5, hlabel=r"$h(t_2)$")

    # inset: falling head series h(t), tying the panels to the tseep concept
    # top-right sky: keeps the waterline + h(t1) label clear (Norm noted the
    # overlap in the side-by-side layout); inset x-margin leaves its own
    # y-tick side unclipped
    axi = ax1.inset_axes([0.70, 0.66, 0.28, 0.30])
    tt = np.array([0.0, 1.0, 6.0, 10.0])
    hh = np.array([6.0, 6.0, 2.0, 2.0])
    axi.plot(tt, hh, color=C_WATER_LINE, lw=1.8)
    axi.plot([1.5, 5.5], [5.6, 2.4], "o", ms=5, color="#123",
             mfc=C_WATER_LINE, mec="white", mew=0.8, zorder=5)
    for tx, hy, lab in [(1.5, 5.6, r"$t_1$"), (5.5, 2.4, r"$t_2$")]:
        axi.annotate(lab, xy=(tx, hy), xytext=(tx + 0.25, hy + 0.6),
                     fontsize=9, color="#123")
    axi.set_title(r"reservoir level $h(t)$", fontsize=9, pad=3)
    axi.annotate(r"$t\rightarrow$", xy=(0.98, 0.02), xycoords="axes fraction",
                 ha="right", va="bottom", fontsize=8.5, color="#555")
    axi.set_xticks([]); axi.set_yticks([])
    axi.set_ylim(1.0, 7.0)
    for s in axi.spines.values():
        s.set_color("#999")

    # shared legend below both panels
    handles = [
        Line2D([], [], ls="none", marker="o", ms=9, mfc=C_SUBMERGED,
               mec="white", mew=1.0, label="held at reservoir head $h(t)$  (submerged)"),
        Line2D([], [], ls="none", marker="o", ms=9, mfc="white", mec=C_EXIT,
               mew=1.8, label="exit face — free to drain  (unsubmerged)"),
        Line2D([], [], ls=(0, (5, 3)), color=C_PHREATIC, lw=1.8,
               label="interior phreatic surface"),
        FancyArrowPatch((0, 0), (1, 0), arrowstyle="-|>", mutation_scale=12,
                        color=C_SEEP, lw=1.7, label="seepage out of exposed face"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=10.5,
               frameon=False, bbox_to_anchor=(0.5, 0.005),
               handletextpad=0.5, columnspacing=1.8)
    fig.suptitle("Submerged-only reservoir boundary: the exit point migrates "
                 "down the face as the level falls", fontsize=13, y=0.985)
    fig.tight_layout(rect=(0, 0.065, 1, 0.955))
    _save(fig, "transient_reservoir_bc.png")


# ===========================================================================
# Figure 2 — dissipation isochrones (erfc column)
# ===========================================================================

def fig_isochrones():
    # h(z,t) = h0 + dH*erfc(z/(2 sqrt(D t))) — the phase-1 analytical lock's form.
    h0, dH, D = 10.0, 1.0, 1.0
    z = np.linspace(0.0, 3.5, 400)
    times = [0.02, 0.06, 0.18, 0.5, 1.5]
    cols = plt.cm.plasma(np.linspace(0.12, 0.82, len(times)))

    fig, ax = plt.subplots(figsize=(6.6, 5.3))
    # initial condition (t = 0): step at the boundary
    ax.plot([h0, h0], [z[-1], 0.0], color="0.6", lw=1.3, ls=(0, (4, 3)))
    ax.plot([h0, h0 + dH], [0.0, 0.0], color="0.6", lw=1.3, ls=(0, (4, 3)))
    ax.annotate("t = 0\n(initial step)", xy=(h0 + 0.02, 0.08), fontsize=8.5,
                color="0.4", ha="left", va="top")

    for t, c in zip(times, cols):
        h = h0 + dH * erfc(z / (2.0 * np.sqrt(D * t)))
        ax.plot(h, z, color=c, lw=2.0, label=f"t = {t:g}")

    # arrow crosses the family steeply (Norm's marked direction); label sits in
    # the open space past the tip so it never overlaps the curves
    ax.annotate("", xy=(10.42, 2.02), xytext=(10.20, 0.62),
                arrowprops=dict(arrowstyle="-|>", color="0.35", lw=1.4))
    ax.text(10.46, 2.12, "increasing t", fontsize=10, color="0.25",
            ha="left", va="top")
    ax.set_xlim(9.97, 11.03)
    ax.set_ylim(3.5, 0.0)                     # depth increases downward
    ax.set_xlabel("hydraulic head  $h$")
    ax.set_ylabel("depth into column  $z$   (stepped boundary at $z=0$)")
    ax.set_title("A transient solution is a family of head fields in time\n"
                 "(boundary step diffusing into a column — erfc)", fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=9, title="saved times")
    fig.tight_layout()
    _save(fig, "transient_isochrones.png")


# ===========================================================================
# Figure 3 — storage coefficient S(psi) for the two families
# ===========================================================================

def fig_storage_models():
    Ss, Sy, h0 = 1.0e-4, 0.20, -2.0        # representative (1/m, -, m)
    a_vg, n_vg = 1.6, 2.2                   # van Genuchten (1/m, -)

    psi = np.linspace(-3.0, 0.5, 700)
    S_vg = storage_capacity_vec(psi, Ss, Sy, h0, vg_a=a_vg, vg_n=n_vg,
                                model=np.full(psi.shape, KR_VG))

    # linear-front / Gardner boxcar, drawn with explicit vertical jumps at h0 and 0.
    # Every level comes from storage_capacity_vec itself, so the figure cannot drift
    # from the solver: dry floor, draining band, saturated Ss.
    def _S_lf(p):
        return float(storage_capacity_vec(np.array([p]), Ss, Sy, h0)[0])

    dry, band, sat = _S_lf(-2.5), _S_lf(-1.0), _S_lf(0.25)
    lf_x = [-3.0, h0, h0, 0.0, 0.0, 0.5]
    lf_y = [dry, dry, band, band, sat, sat]

    fig, ax = plt.subplots(figsize=(7.4, 4.9))
    ax.axvspan(h0, 0.0, color="#f2e6cf", zorder=0)          # draining band shade
    ax.plot(lf_x, lf_y, color=C_EXIT, lw=2.2,
            label="linear-front / Gardner")
    ax.plot(psi, S_vg, color=C_WATER_LINE, lw=2.2,
            label="van Genuchten")
    ax.axhline(Ss, color="0.45", lw=1.0, ls=(0, (4, 3)))

    ax.set_yscale("log")
    ax.set_xlim(-3.0, 0.5)
    ipk = int(np.nanargmax(S_vg))
    peak = float(S_vg[ipk])
    # headroom: the vG-peak annotation is two text lines anchored (va="bottom")
    # at 2x the peak, so the top must clear ~a further half-decade above that
    top = max(band, peak) * 7.0
    ax.set_ylim(min(dry, float(np.nanmin(S_vg))) * 0.4, top)

    # Ss, which applies in the saturated zone only (dashed line, labelled at the left)
    ax.annotate(r"$S_s$ (saturated zone only)", xy=(-2.92, Ss * 1.45),
                fontsize=9.3, color="0.35", ha="left", va="bottom")
    # draining band, seated low inside the shaded band (shading shows the extent)
    ax.annotate(r"$S_y/|h_0|$ draining band" + "\n" + r"$h_0<\psi<0$",
                xy=(-1.0, 4.0e-4), fontsize=9.3, color="#9a5c14",
                ha="center", va="center")
    # van Genuchten capacity peak, labelled from above
    ax.annotate("van Genuchten peak\n(near air entry)",
                xy=(psi[ipk], peak), xytext=(psi[ipk] - 0.15, peak * 2.0),
                fontsize=9.3, color=C_WATER_LINE, ha="center", va="bottom",
                arrowprops=dict(arrowstyle="-|>", color=C_WATER_LINE, lw=1.1))
    # saturated zone (bottom-right)
    ax.annotate("saturated\n$\\psi\\geq0$", xy=(0.28, Ss * 3.4), fontsize=8.6,
                color="0.4", ha="center", va="center")
    # residual floor carried by the drained end of both curves
    ax.annotate("residual floor", xy=(-2.55, dry * 1.5), fontsize=8.6,
                color="0.4", ha="center", va="bottom")

    ax.set_xlabel(r"pressure head  $\psi$")
    ax.set_ylabel(r"storage coefficient  $S(\psi)$   [1/length]")
    ax.set_title("Storage capacity of the two unsaturated families", fontsize=11.5)
    ax.legend(loc="upper left", fontsize=9.3)
    ax.grid(alpha=0.22, which="both")
    fig.tight_layout()
    _save(fig, "transient_storage_models.png")


# ===========================================================================
# Figure 4 — time-series breakpoint semantics
# ===========================================================================

def fig_series_semantics():
    # defined points (a repeated time at t=4 makes the step)
    tp = [1.0, 4.0, 4.0, 7.0]
    vp = [10.0, 7.0, 3.0, 3.0]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    # held-constant continuations (dashed) beyond first / last entry
    ax.plot([-0.4, 1.0], [10.0, 10.0], color=C_WATER_LINE, lw=1.6,
            ls=(0, (5, 3)))
    ax.plot([7.0, 8.4], [3.0, 3.0], color=C_WATER_LINE, lw=1.6, ls=(0, (5, 3)))
    # the series itself (linear segments + the vertical step)
    ax.plot(tp, vp, color=C_WATER_LINE, lw=2.2)
    # defined points as markers
    ax.plot(tp, vp, "o", ms=8, mfc="white", mec=C_WATER_LINE, mew=1.8,
            zorder=5)

    ax.annotate("defined points", xy=(1.0, 10.0), xytext=(1.4, 11.0),
                fontsize=9.5, color="0.2",
                arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.0))
    ax.annotate("linear interpolation", xy=(2.5, 8.5), xytext=(3.2, 9.6),
                fontsize=9.5, color="0.2",
                arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.0))
    ax.annotate("repeated time = step", xy=(4.0, 5.0), xytext=(4.5, 6.4),
                fontsize=9.5, color="#b5460f",
                arrowprops=dict(arrowstyle="-|>", color="#b5460f", lw=1.2))
    ax.annotate("held constant\nbefore first entry", xy=(0.2, 10.0),
                xytext=(-0.3, 7.9), fontsize=8.8, color="0.35", ha="left")
    ax.annotate("held constant\nafter last entry", xy=(7.8, 3.0),
                xytext=(6.0, 1.4), fontsize=8.8, color="0.35", ha="left")

    ax.set_xlim(-0.5, 8.5)
    ax.set_ylim(0.5, 11.8)
    ax.set_xlabel("time  $t$")
    ax.set_ylabel("series value  $h$")
    ax.set_title("Reading a tseep series: linear between points, "
                 "held beyond the ends,\nvertical at a repeated time", fontsize=10.5)
    ax.grid(alpha=0.22)
    fig.tight_layout()
    _save(fig, "transient_series_semantics.png")


# ===========================================================================
# Figure 5 — saved-frame schedule (union timing diagram)
# ===========================================================================

def fig_frame_schedule():
    T = 20.0
    grid = [4, 8, 12, 16, 20]              # save_interval grid
    save_times = [3, 11]
    stages = {"stage_1": 2, "stage_2": 14}
    breaks = [1, 7, 14]                    # 14 coincides with stage_2 (dedup)

    union = sorted(set([0] + grid + save_times +
                       list(stages.values()) + breaks))

    lanes = [
        (4, grid, "save_interval grid", "#6a6a6a", "|", 12, False),
        (3, save_times, "save_times", "#2b7bb0", "D", 7, True),
        (2, list(stages.values()), "stage_1 / stage_2", "#b5460f", "^", 9, True),
        (1, breaks, "series breakpoints", "#5a3f9c", "x", 8, True),
    ]

    fig, ax = plt.subplots(figsize=(8.6, 4.1))
    # faint vertical guides at every union time (shows alignment + dedup)
    for u in union:
        ax.plot([u, u], [-0.15, 4.4], color="0.85", lw=0.8, zorder=0)

    for y, ts, lab, col, mk, ms, filled in lanes:
        ax.plot([0, T], [y, y], color="0.8", lw=0.8, zorder=1)
        if mk == "|":
            ax.plot(ts, [y] * len(ts), marker="|", ls="none", color=col,
                    ms=15, mew=2.0, zorder=3)
        else:
            ax.plot(ts, [y] * len(ts), marker=mk, ls="none", color=col,
                    ms=ms, mfc=(col if filled else "white"), mew=1.6, zorder=3)
        ax.annotate(lab, xy=(-0.4, y), ha="right", va="center", fontsize=9.5,
                    color=col)

    # stage flag labels
    for lab, tx in stages.items():
        ax.annotate(lab.replace("_", r"$_") + "$", xy=(tx, 2), xytext=(tx, 2.42),
                    ha="center", fontsize=8.2, color="#b5460f")

    # merged lane
    ax.plot([0, T], [0, 0], color="0.25", lw=1.2, zorder=1)
    ax.plot(union, [0] * len(union), marker="|", ls="none", color="0.1",
            ms=18, mew=2.4, zorder=3)
    ax.annotate("saved frames\n(union)", xy=(-0.4, 0), ha="right", va="center",
                fontsize=9.5, fontweight="bold", color="0.1")
    ax.annotate(r"$t=0$ always saved", xy=(0, 0), xytext=(0.3, -0.62),
                fontsize=8.5, color="0.3", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.0))
    ax.annotate("coincident times\nde-duplicated", xy=(14, 0.0),
                xytext=(15.4, -0.7), fontsize=8.3, color="0.3", ha="left",
                arrowprops=dict(arrowstyle="-|>", color="0.4", lw=1.0))

    # time axis
    ax.annotate("", xy=(T + 0.6, -1.05), xytext=(0, -1.05),
                arrowprops=dict(arrowstyle="-|>", color="0.3", lw=1.2))
    ax.annotate("time  $t$", xy=(T + 0.6, -1.05), xytext=(T + 0.7, -1.05),
                fontsize=9.5, color="0.3", va="center")

    ax.set_xlim(-4.6, T + 2.2)
    ax.set_ylim(-1.5, 4.9)
    ax.axis("off")
    ax.set_title("The saved-frame schedule is the union of every source, "
                 "de-duplicated and sorted", fontsize=11)
    fig.tight_layout()
    _save(fig, "transient_frame_schedule.png")


def _save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", os.path.relpath(path))


if __name__ == "__main__":
    for fn in (fig_reservoir_bc, fig_isochrones, fig_storage_models,
               fig_series_semantics, fig_frame_schedule):
        fn()
