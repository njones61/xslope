"""Render the reinforcement illustration figures for the two reinforcement pages.

Two figures land in docs/lem/images/, two in docs/fem/images/:

  reinf_envelope.png    The LEM capacity envelope T(s) along a reinforcement
                        line, for the four end conditions the page describes:
                        friction-only free ends, a plated end, a fully anchored
                        end (Lp = 0), and a line shorter than Lp1 + Lp2. Every
                        curve is evaluated with the shipping
                        xslope.fileio.reinforce_available_tension, so the figure
                        and the solver cannot disagree.
  reinf_direction.png   How the crossing-point force enters slice equilibrium
                        for Dir = Tangent and Dir = Axial: the angle psi from
                        horizontal, and the resolution into components normal
                        and tangential to the slice base.
  reinf_bar_law.png     The truss-bar constitutive law: tension only, the
                        elastic slope EA/L, and the three post-peak models
                        (Tres blank, Tres entered, Tres = 0).
  reinf_bond_slip.png   The optional bond-slip envelope against the fixed
                        pullout ramp, on a line whose overburden grows from one
                        end to the other.

All four are closed-form schematics; nothing is solved.

Run from the repo root:

    PYTHONPATH=. python3 tools/make_reinforcement_figures.py
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, Polygon

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from xslope.fileio import reinforce_available_tension

HERE = os.path.dirname(os.path.abspath(__file__))
LEM_OUT = os.path.join(HERE, "..", "docs", "lem", "images")
FEM_OUT = os.path.join(HERE, "..", "docs", "fem", "images")

# ---- shared palette (matches tools/make_fem_overview_figures.py) ----------
C_INK = "#1f2933"
C_MUTED = "#6b7785"
C_SOIL = "#f6e3d0"
C_SOIL_EDGE = "#1f2933"
C_ACCENT = "#1f6fb2"        # the primary curve / the fixed ramp
C_FORCE = "#c0392b"         # the reinforcement force P
C_GREEN = "#2f7d5b"         # anchorage / bond envelope
C_ORANGE = "#c05621"        # the contrasting curve
C_GRID = "#d7dde3"


def _finish(fig, out_dir, name):
    path = os.path.abspath(os.path.join(out_dir, name))
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", path)


def _tidy(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(C_MUTED)
    ax.tick_params(colors=C_MUTED, labelsize=9)
    ax.grid(True, color=C_GRID, linewidth=0.6)
    ax.set_axisbelow(True)


# ===========================================================================
# Figure 1 — the LEM capacity envelope
# ===========================================================================

def fig_capacity_envelope():
    Tmax = 100.0

    cases = [
        dict(title="Free ends — friction only\n$T_{end}=0$, both ends taper",
             L=10.0, lp1=2.5, lp2=2.5, t1=0.0, t2=0.0),
        dict(title="Plate at end 1\n$T_{end1}>0$, friction adds to it",
             L=10.0, lp1=2.5, lp2=2.5, t1=35.0, t2=0.0),
        dict(title="End 1 fully anchored\n$L_{p1}=0$",
             L=10.0, lp1=0.0, lp2=2.5, t1=0.0, t2=0.0),
        dict(title="Line shorter than $L_{p1}+L_{p2}$\nno anchorage",
             L=3.5, lp1=2.5, lp2=2.5, t1=0.0, t2=0.0),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 6.6))

    for ax, cs in zip(axes.ravel(), cases):
        L = cs["L"]
        s = np.linspace(0.0, L, 801)
        T = np.array([reinforce_available_tension(si, L - si, Tmax,
                                                  cs["lp1"], cs["lp2"],
                                                  cs["t1"], cs["t2"])
                      for si in s])

        y_top = Tmax * 1.52
        ax.axhline(Tmax, color=C_MUTED, linewidth=1.0, linestyle=(0, (5, 4)))
        ax.text(L * 0.985, Tmax + 3, "$T_{max}$", fontsize=9.5, color=C_MUTED,
                ha="right", va="bottom")

        ax.fill_between(s, 0, T, color=C_ACCENT, alpha=0.13, zorder=1)
        ax.plot(s, T, color=C_ACCENT, linewidth=2.2, zorder=3)

        # pullout lengths, drawn as spans in the clear band above the envelope
        for lp, at_end1, y_span, label in ((cs["lp1"], True, 0.94, "$L_{p1}$"),
                                           (cs["lp2"], False, 0.83, "$L_{p2}$")):
            if lp <= 0:
                continue
            x0, x1 = ((0.0, min(lp, L)) if at_end1
                      else (L - min(lp, L), L))
            y = y_span * y_top
            ax.annotate("", xy=(x0, y), xytext=(x1, y),
                        arrowprops=dict(arrowstyle="<->", color=C_ORANGE,
                                        linewidth=1.1))
            ax.text((x0 + x1) / 2, y - 4, label, fontsize=9.5, color=C_ORANGE,
                    ha="center", va="top")

        if cs["t1"] > 0:
            ax.plot([0], [cs["t1"]], marker="o", markersize=6, color=C_GREEN,
                    zorder=5)
            ax.annotate("$T_{end1}$", xy=(0, cs["t1"]),
                        xytext=(L * 0.10, cs["t1"] + 26), fontsize=9.5,
                        color=C_GREEN,
                        arrowprops=dict(arrowstyle="->", color=C_GREEN,
                                        linewidth=1.0))
        if cs["lp1"] == 0:
            ax.plot([0], [Tmax], marker="o", markersize=6, color=C_GREEN,
                    zorder=5)
            ax.annotate("full capacity\nat the end", xy=(0, Tmax),
                        xytext=(L * 0.09, Tmax * 0.60), fontsize=9,
                        color=C_GREEN,
                        arrowprops=dict(arrowstyle="->", color=C_GREEN,
                                        linewidth=1.0))
        if T.max() < Tmax - 1e-9:
            i = int(np.argmax(T))
            ax.plot([s[i]], [T[i]], marker="o", markersize=6, color=C_FORCE,
                    zorder=5)
            ax.annotate("ramps cross below $T_{max}$", xy=(s[i], T[i]),
                        xytext=(s[i], T[i] + 34), fontsize=9, color=C_FORCE,
                        ha="center",
                        arrowprops=dict(arrowstyle="->", color=C_FORCE,
                                        linewidth=1.0))

        ax.set_xlim(0, L)
        ax.set_ylim(0, y_top)
        ax.set_title(cs["title"], fontsize=10, color=C_INK)
        ax.set_xlabel("distance along the line   $s$", fontsize=9.5)
        ax.set_ylabel("available tension   $T$", fontsize=9.5)
        _tidy(ax)

    fig.tight_layout(h_pad=2.6, w_pad=2.2)
    _finish(fig, LEM_OUT, "reinf_envelope.png")


# ===========================================================================
# Figure 2 — Dir = Tangent against Dir = Axial
# ===========================================================================

def _slice_panel(ax, alpha_deg, psi_deg, line_deg, title, note):
    """One slice with a reinforcement line crossing its base at r."""
    a = np.radians(alpha_deg)
    ta = np.tan(a)

    def base_y(x):
        return x * ta

    # --- slice body: vertical sides, base at alpha, flat-ish top
    xl, xr = -1.75, 1.15
    ytop = 2.35
    ax.add_patch(Polygon([[xl, base_y(xl)], [xr, base_y(xr)],
                          [xr, ytop + 0.16 * xr], [xl, ytop + 0.16 * xl]],
                         closed=True, facecolor=C_SOIL, edgecolor=C_SOIL_EDGE,
                         linewidth=1.6, zorder=2))

    # --- the slip surface, running past the slice
    xs = np.array([-3.5, 3.0])
    ax.plot(xs, base_y(xs), color=C_INK, linewidth=2.0, zorder=3)
    ax.text(-2.05, base_y(-2.05) - 0.16, "slip surface", fontsize=9,
            color=C_INK, ha="left", va="top", rotation=alpha_deg,
            rotation_mode="anchor")

    # --- alpha at the left, from a horizontal dashed reference
    ax.plot([-3.3, -2.15], [base_y(-3.3), base_y(-3.3)], color=C_MUTED,
            linewidth=1.0, linestyle=(0, (3, 3)), zorder=3)
    ax.add_patch(Arc((-3.3, base_y(-3.3)), 1.5, 1.5, theta1=0.0,
                     theta2=alpha_deg, color=C_MUTED, linewidth=1.1, zorder=4))
    ax.text(-2.42, base_y(-3.3) + 0.16, r"$\alpha$", fontsize=12,
            color=C_MUTED, ha="center", va="bottom")

    # --- the reinforcement line through the crossing point r = (0, 0)
    lr = np.radians(line_deg)
    ux, uy = np.cos(lr), np.sin(lr)
    ax.plot([-3.3 * ux, 2.55 * ux], [-3.3 * uy, 2.55 * uy], color=C_GREEN,
            linewidth=2.4, zorder=4)
    ax.text(-3.25 * ux, -3.25 * uy + 0.12, "reinforcement line", fontsize=8.5,
            color=C_GREEN, ha="left", va="bottom", rotation=line_deg,
            rotation_mode="anchor")

    # --- P at psi from horizontal, applied at r
    p = np.radians(psi_deg)
    Plen = 1.75
    px, py = Plen * np.cos(p), Plen * np.sin(p)
    ax.plot([0, 1.15], [0, 0], color=C_MUTED, linewidth=1.0,
            linestyle=(0, (3, 3)), zorder=5)
    ax.add_patch(FancyArrowPatch((0, 0), (px, py), arrowstyle="-|>",
                                 mutation_scale=16, linewidth=2.2,
                                 color=C_FORCE, zorder=7))
    ax.text(px + 0.10, py + 0.10, "$P$", fontsize=13, color=C_FORCE,
            ha="left", va="bottom")
    ax.add_patch(Arc((0, 0), 1.55, 1.55, theta1=min(0.0, psi_deg),
                     theta2=max(0.0, psi_deg), color=C_FORCE, linewidth=1.2,
                     zorder=6))
    ax.text(0.93, 0.30 * np.sign(psi_deg if psi_deg else 1) - 0.02,
            r"$\psi$", fontsize=12, color=C_FORCE, ha="left", va="center")
    ax.plot([0], [0], marker="o", markersize=8, color=C_FORCE, zorder=8)
    ax.text(-0.16, 0.06, "$r$", fontsize=12, color=C_INK, ha="right",
            va="bottom")

    # --- components normal and tangential to the base
    tn = np.array([np.cos(a), np.sin(a)])          # base tangent, up-slope
    nn = np.array([-np.sin(a), np.cos(a)])         # base normal, outward
    Pt = Plen * np.cos(p - a)
    Pn = Plen * np.sin(a - p)
    if abs(Pn) > 1e-9:
        ax.add_patch(FancyArrowPatch((0, 0), tuple(Pt * tn), arrowstyle="-|>",
                                     mutation_scale=12, linewidth=1.4,
                                     color=C_FORCE, linestyle=(0, (4, 2.5)),
                                     zorder=6))
        ax.add_patch(FancyArrowPatch((0, 0), tuple(-Pn * nn), arrowstyle="-|>",
                                     mutation_scale=12, linewidth=1.4,
                                     color=C_FORCE, linestyle=(0, (4, 2.5)),
                                     zorder=6))
        ax.plot([Pt * tn[0], px, -Pn * nn[0]], [Pt * tn[1], py, -Pn * nn[1]],
                color=C_FORCE, linewidth=0.9, linestyle=(0, (2, 3)), zorder=5)
        e = -Pn * nn
        ax.text(e[0] + 0.06, e[1] - 0.16, r"$P\sin(\alpha-\psi)$", fontsize=9.5,
                color=C_FORCE, ha="center", va="top")
        e = Pt * tn
        ax.text(e[0] + 0.26, e[1] + 0.14, r"$P\cos(\alpha-\psi)$", fontsize=9.5,
                color=C_FORCE, ha="left", va="bottom")
    else:
        ax.text(px + 0.16, py - 0.14,
                r"$P\cos(\alpha-\psi)=P$" + "\n" + r"$P\sin(\alpha-\psi)=0$",
                fontsize=9.5, color=C_FORCE, ha="left", va="top")

    ax.set_title(title, fontsize=11, color=C_INK)
    ax.text(0.5, -0.045, note, transform=ax.transAxes, fontsize=9.5,
            color=C_INK, ha="center", va="top")
    ax.set_xlim(-3.6, 4.35)
    ax.set_ylim(-2.1, 3.0)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def fig_force_direction():
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.0))

    _slice_panel(axes[0], alpha_deg=25.0, psi_deg=25.0, line_deg=0.0,
                 title="Dir = Tangent   (flexible reinforcement)",
                 note=r"$\psi=\alpha$: the force turns with the sliding mass,"
                      "\nnot with the line — no component normal to the base")
    _slice_panel(axes[1], alpha_deg=25.0, psi_deg=-18.0, line_deg=-18.0,
                 title="Dir = Axial   (rigid supports)",
                 note=r"$\psi$ = the inclination of the line: the normal"
                      "\ncomponent adds frictional resistance on the base")

    fig.tight_layout(w_pad=1.0)
    _finish(fig, LEM_OUT, "reinf_direction.png")


# ===========================================================================
# Figure 3 — the truss-bar constitutive law
# ===========================================================================

def fig_bar_law():
    T_allow = 1.0
    T_res = 0.45
    d_y = 1.0                      # elongation at first yield
    d_end = 3.4

    fig, ax = plt.subplots(figsize=(8.0, 4.6))

    # compression side: the bar carries nothing
    ax.plot([-1.35, 0], [0, 0], color=C_INK, linewidth=2.6, zorder=4)
    ax.text(-1.30, 0.10, "tension only —\na compressed bar\ncarries nothing",
            fontsize=9, color=C_INK, ha="left", va="bottom")

    # elastic branch, slope EA/L
    d = np.array([0.0, d_y])
    ax.plot(d, d * T_allow / d_y, color=C_INK, linewidth=2.6, zorder=4)
    ax.annotate("slope $EA/L$", xy=(0.66 * d_y, 0.66 * T_allow),
                xytext=(0.02 * d_y, 0.72 * T_allow), fontsize=9.5, color=C_INK,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color=C_INK, linewidth=1.0))

    # the three post-peak models
    ax.plot([d_y, d_end], [T_allow, T_allow], color=C_ACCENT, linewidth=2.6,
            zorder=5, label="$T_{res}$ blank — elastic-perfectly-plastic (default)")
    ax.plot([d_y, d_y, d_end], [T_allow, T_res, T_res], color=C_ORANGE,
            linewidth=2.4, linestyle=(0, (6, 3)), zorder=5,
            label="$T_{res}$ entered — peak-residual")
    ax.plot([d_y, d_y, d_end], [T_allow, 0.0, 0.0], color=C_FORCE,
            linewidth=2.2, linestyle=(0, (1.6, 2.4)), zorder=5,
            label="$T_{res}=0$ — brittle rupture")

    ax.axhline(T_allow, color=C_MUTED, linewidth=1.0, linestyle=(0, (5, 4)),
               zorder=2)
    ax.axhline(T_res, color=C_MUTED, linewidth=1.0, linestyle=(0, (5, 4)),
               zorder=2)

    ax.annotate("the drop is applied to a converged\nstate, and the model re-solved",
                xy=(d_y, 0.5 * (T_allow + T_res)),
                xytext=(d_y + 0.30, 0.72 * T_allow), fontsize=9,
                color=C_ORANGE, ha="left", va="top",
                arrowprops=dict(arrowstyle="->", color=C_ORANGE, linewidth=1.0))

    ax.set_xlim(-1.45, d_end + 0.05)
    ax.set_ylim(-0.10, T_allow * 1.42)
    ax.set_xlabel(r"elongation   $\delta$", fontsize=9.5)
    ax.set_ylabel("axial force   $T$", fontsize=9.5)
    ax.set_xticks([0])
    ax.set_xticklabels(["0"])
    ax.set_yticks([0, T_res, T_allow])
    ax.set_yticklabels(["0", "$T_{res}$", "$T_{allow}$"])
    _tidy(ax)
    ax.legend(loc="upper left", fontsize=9, frameon=False)
    fig.tight_layout()
    _finish(fig, FEM_OUT, "reinf_bar_law.png")


# ===========================================================================
# Figure 4 — the bond-slip envelope against the fixed pullout ramp
# ===========================================================================

def fig_bond_slip():
    L = 8.0                        # line length, m
    Tmax = 100.0                   # kN/m
    Lp = 2.0                       # fixed pullout length at both ends, m
    gamma = 20.0                   # kN/m3
    perim = 2.0                    # geotextile: friction on both faces
    c_bond = 0.0                   # kPa
    phi_bond = np.radians(28.35)

    s = np.linspace(0.0, L, 1601)
    depth = 0.5 + (6.0 - 0.5) * s / L          # overburden grows into the fill
    sig_n = gamma * depth
    q = perim * (c_bond + sig_n * np.tan(phi_bond))

    from_1 = np.concatenate([[0.0], np.cumsum(np.diff(s) * 0.5 * (q[1:] + q[:-1]))])
    from_2 = from_1[-1] - from_1
    T_bond = np.minimum(Tmax, np.minimum(from_1, from_2))

    T_fixed = np.array([reinforce_available_tension(si, L - si, Tmax, Lp, Lp)
                        for si in s])

    fig, (ax_q, ax_t) = plt.subplots(2, 1, figsize=(8.6, 6.4), sharex=True,
                                     gridspec_kw=dict(height_ratios=[1.0, 1.55]))

    ax_q.fill_between(s, 0, q, color=C_GREEN, alpha=0.13)
    ax_q.plot(s, q, color=C_GREEN, linewidth=2.2)
    ax_q.text(0.25, q.max() * 1.34,
              r"$q(s)=P\,(c_{bond}+\sigma_n(s)\tan\phi_{bond})$",
              fontsize=10, color=C_GREEN, ha="left", va="top")
    ax_q.set_ylabel("bond rate   $q$", fontsize=9.5)
    ax_q.set_ylim(0, q.max() * 1.42)
    _tidy(ax_q)

    ax_t.plot(s, T_fixed, color=C_ACCENT, linewidth=2.2,
              linestyle=(0, (6, 3)), label="fixed ramp, $L_p$ at each end")
    ax_t.plot(s, T_bond, color=C_GREEN, linewidth=2.4,
              label="bond envelope, $\\int q\\,ds$ from each end")
    ax_t.axhline(Tmax, color=C_MUTED, linewidth=1.0, linestyle=(0, (5, 4)))
    ax_t.text(L, Tmax + 3, "$T_{max}$", fontsize=9.5, color=C_MUTED,
              ha="right", va="bottom")

    ax_t.annotate("shallow end: force develops\nmore slowly than the ramp allows",
                  xy=(1.15, np.interp(1.15, s, T_bond)),
                  xytext=(3.35, Tmax * 0.72), fontsize=9, color=C_INK,
                  ha="left", va="center",
                  arrowprops=dict(arrowstyle="->", color=C_INK, linewidth=1.0))
    ax_t.annotate("deep end: the thick overburden\ndevelops it in a fraction of $L_p$",
                  xy=(7.35, np.interp(7.35, s, T_bond)),
                  xytext=(6.55, Tmax * 0.40), fontsize=9, color=C_INK,
                  ha="right", va="center",
                  arrowprops=dict(arrowstyle="->", color=C_INK, linewidth=1.0))

    ax_t.set_xlim(0, L)
    ax_t.set_ylim(0, Tmax * 1.30)
    ax_t.set_xlabel("distance along the line   $s$   (end 1 at the face,"
                    " end 2 under the fill)", fontsize=9.5)
    ax_t.set_ylabel("available tension   $T_{allow}$", fontsize=9.5)
    _tidy(ax_t)
    ax_t.legend(loc="lower center", fontsize=9, frameon=False, ncol=2)

    fig.tight_layout(h_pad=1.2)
    _finish(fig, FEM_OUT, "reinf_bond_slip.png")


if __name__ == "__main__":
    fig_capacity_envelope()
    fig_force_direction()
    fig_bar_law()
    fig_bond_slip()
