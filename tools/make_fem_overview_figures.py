"""Render the illustration figures for docs/fem/overview.md.

Four figures land in docs/fem/images/ with the fem_ov_ prefix:

  fem_ov_viscoplastic_loop.png  Flow diagram of the viscoplastic algorithm: one
                                factorization, the per-Gauss-point yield check, the
                                body-load correction, and the two convergence tests
                                with the hybrid verdict that follows them.
  fem_ov_tension_cutoff.png     The Mohr-Coulomb envelope in the tensile quadrant —
                                the implicit apex at c/tan(phi), its invariance under
                                strength reduction, and the Rankine cap t_cut that
                                removes it.
  fem_ov_k0_initial.png         Initial lateral effective stress with depth: the
                                gravity turn-on's nu/(1-nu) coefficient against
                                stated K0 values.
  fem_ov_ssrm_sweep.png         Strength-reduction sweep on the Griffiths & Lane
                                Example 1 sample file — viscoplastic displacement
                                against F, marked by whether the trial reached
                                equilibrium, with the bisection's factor of safety.

The first three are closed-form schematics. The fourth runs the solver on the
committed sample docs/fem/files/xslope_griffiths1.xlsx (about a minute).

Run from the repo root:  python tools/make_fem_overview_figures.py
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "docs", "fem", "images")
SAMPLE = os.path.join(HERE, "..", "docs", "fem", "files", "xslope_griffiths1.xlsx")

# ---- shared palette -------------------------------------------------------
C_INK = "#1f2933"
C_MUTED = "#6b7785"
C_BOX = "#eef3f7"
C_BOX_EDGE = "#7a97ad"
C_ACCENT = "#1f6fb2"        # the primary/current curve
C_ACCENT2 = "#c05621"       # the reduced / contrasting curve
C_GREEN = "#2f7d5b"         # stable / equilibrium
C_RED = "#b03a2e"           # failed
C_GRID = "#d7dde3"


def _finish(fig, name):
    path = os.path.abspath(os.path.join(OUT, name))
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
# Figure 1 — the viscoplastic loop
# ===========================================================================

def fig_viscoplastic_loop():
    fig, ax = plt.subplots(figsize=(8.2, 8.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 11.6)
    ax.axis("off")

    def box(x, y, w, h, text, fill=C_BOX, edge=C_BOX_EDGE, fontsize=9.5, style="round,pad=0.02,rounding_size=0.12"):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle=style,
                                    facecolor=fill, edgecolor=edge, linewidth=1.2))
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize,
                color=C_INK, linespacing=1.45)

    def arrow(x1, y1, x2, y2, color=C_MUTED, text=None, tx=0.0, ha="left"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=13, linewidth=1.2,
                                     color=color, shrinkA=1, shrinkB=1))
        if text:
            ax.text((x1 + x2) / 2 + tx, (y1 + y2) / 2, text, fontsize=8.5,
                    color=color, ha=ha, va="center")

    cx = 4.6

    box(cx, 11.0, 7.4, 0.85,
        "Assemble the elastic stiffness $[K]$ from $[D_e]$ and factorize it — once")
    arrow(cx, 10.57, cx, 10.15)
    box(cx, 9.7, 7.4, 0.85,
        "Elastic solution   $[K]\\{U\\} = \\{F\\}_{\\mathrm{applied}}$,   $\\{\\varepsilon^{vp}\\} = 0$")
    arrow(cx, 9.27, cx, 8.85)

    # --- iteration band
    ax.add_patch(FancyBboxPatch((0.35, 3.05), 8.5, 5.72,
                                boxstyle="round,pad=0.02,rounding_size=0.15",
                                facecolor="#f7fafc", edgecolor=C_MUTED,
                                linewidth=1.0, linestyle=(0, (5, 4))))
    ax.text(0.55, 8.56, "viscoplastic iteration", fontsize=9, color=C_MUTED,
            style="italic", ha="left", va="center")

    box(cx, 7.85, 7.2, 0.9,
        "At every Gauss point:  $\\{\\sigma\\} = [D_e](\\{\\varepsilon\\} - \\{\\varepsilon^{vp}\\})$\n"
        "on the elastic strain, effective stress")
    arrow(cx, 7.4, cx, 6.98)
    box(cx, 6.5, 7.2, 0.95,
        "Yield function $f$ (Mohr-Coulomb, invariant form)\n"
        "$f > 0$:  $\\Delta\\varepsilon^{vp} = f\\,\\frac{\\partial Q}{\\partial\\sigma}\\,\\Delta t$   accumulated")
    arrow(cx, 6.02, cx, 5.6)
    box(cx, 5.1, 7.2, 0.9,
        "Body-load correction\n"
        "$\\{F\\} = \\{F\\}_{\\mathrm{applied}} + \\sum_e \\int [B]^T[D_e]\\{\\varepsilon^{vp}\\}\\,dA$")
    arrow(cx, 4.65, cx, 4.23)
    box(cx, 3.75, 7.2, 0.85,
        "Re-solve $[K]\\{U\\} = \\{F\\}$ — same factorization, back-substitution only")

    # feedback route up the left side, outside the boxes
    ax.plot([1.0, 0.72, 0.72], [3.75, 3.75, 7.85], color=C_MUTED, linewidth=1.2,
            solid_capstyle="round", zorder=1)
    ax.add_patch(FancyArrowPatch((0.72, 7.85), (1.0, 7.85), arrowstyle="-|>",
                                 mutation_scale=13, linewidth=1.2, color=C_MUTED))
    ax.text(0.52, 5.7, "not yet in equilibrium", fontsize=8.5, color=C_MUTED,
            rotation=90, ha="center", va="center")

    arrow(cx, 3.3, cx, 2.86)
    box(cx, 2.35, 7.4, 1.05,
        "Equilibrium requires BOTH\n"
        "$\\max|\\Delta U| / \\max|U| < \\mathrm{tol}$   and   "
        "$\\max_i |\\mathbf{r}_i| / |\\mathbf{f}^{\\,grav}_i| < \\mathrm{force\\_tol}$",
        fill="#eaf2ea", edge=C_GREEN)

    arrow(2.3, 1.83, 2.3, 1.42, color=C_GREEN)
    ax.text(2.15, 1.63, "met", fontsize=8.5, color=C_GREEN, ha="right", va="center")
    box(2.3, 0.95, 3.6, 0.85, "CONVERGED\nthe slope stands at this $F$",
        fill="#eaf2ea", edge=C_GREEN, fontsize=9)

    arrow(6.9, 1.83, 6.9, 1.42, color=C_RED)
    ax.text(7.05, 1.63, "budget spent", fontsize=8.5, color=C_RED, ha="left", va="center")
    box(6.9, 0.95, 3.9, 0.85,
        "Displacement history classified\nFAILED / STABLE_STUCK / AMBIGUOUS",
        fill="#f9ecea", edge=C_RED, fontsize=9)

    _finish(fig, "fem_ov_viscoplastic_loop.png")


# ===========================================================================
# Figure 2 — the tensile quadrant: apex, reduction, Rankine cap
# ===========================================================================

def fig_tension_cutoff():
    c, phi = 20.0, 35.0                 # kPa, degrees
    F = 2.0
    tanphi = np.tan(np.radians(phi))
    apex = -c / tanphi                  # sigma' at tau = 0  (about -28.6 kPa)
    t_cut = 5.0                         # a stated tensile strength, kPa

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3), sharey=True)

    for ax, cap in zip(axes, (None, t_cut)):
        s = np.linspace(-40, 60, 400)
        env = c + s * tanphi
        env_r = c / F + s * (tanphi / F)

        ax.plot(s[env >= 0], env[env >= 0], color=C_ACCENT, linewidth=2.0)
        ax.plot(s[env_r >= 0], env_r[env_r >= 0], color=C_ACCENT2, linewidth=2.0,
                linestyle=(0, (6, 3)))
        ax.axhline(0, color=C_MUTED, linewidth=0.9)
        ax.axvline(0, color=C_MUTED, linewidth=0.9)

        ax.plot([apex], [0], marker="o", markersize=7, color=C_INK, zorder=5)
        ax.annotate("apex  $\\sigma'_t = -c/\\tan\\phi$\nunmoved by $F$",
                    xy=(apex, 0.4), xytext=(apex + 3, 13),
                    fontsize=9, color=C_INK, ha="left",
                    arrowprops=dict(arrowstyle="->", color=C_INK, linewidth=1.0))

        ax.text(26, 34, "$\\tau = c + \\sigma'\\tan\\phi$",
                color=C_ACCENT, fontsize=9.5, ha="left", va="top", zorder=6)
        ax.text(44, 19.5, "reduced by $F$",
                color=C_ACCENT2, fontsize=9.5, ha="left", va="top", zorder=6)

        if cap is None:
            ax.fill_between([apex, 0], 0, 46, color="#f3d9d4", alpha=0.75, zorder=0)
            ax.text((apex) / 2, 40, "tension carried\nby the criterion alone",
                    fontsize=9, color=C_RED, ha="center", va="center")
            ax.set_title("No cutoff — the Griffiths & Lane convention",
                         fontsize=10.5, color=C_INK)
        else:
            # tension is negative on this axis, so a tensile strength T sits at -T
            ax.axvline(-cap, color=C_GREEN, linewidth=1.8)
            ax.axvline(-cap / F, color=C_GREEN, linewidth=1.4, linestyle=(0, (4, 3)))
            ax.fill_between([apex, -cap], 0, 46, color="#e6e9ec", alpha=0.9, zorder=0)
            ax.text(-cap + 1.5, 40.5, "$-T$", fontsize=9, color=C_GREEN,
                    ha="left", va="center")
            ax.text(-cap / F + 1.5, 35.5, "$-T/F$", fontsize=9, color=C_GREEN,
                    ha="left", va="center")
            ax.text((apex - cap) / 2 - 1, 40, "removed by the\nRankine cap",
                    fontsize=9, color=C_INK, ha="center", va="center")
            ax.set_title("With a stated tensile strength (mat!t_cut)",
                         fontsize=10.5, color=C_INK)

        ax.set_xlim(-40, 60)
        ax.set_ylim(0, 46)
        ax.set_xlabel("effective normal stress $\\sigma'$   (tension negative)", fontsize=9.5)
        _tidy(ax)

    axes[0].set_ylabel("shear stress $\\tau$", fontsize=9.5)
    fig.text(0.5, -0.03,
             "c = 20 kPa, $\\phi$ = 35°, strength reduced by F = 2",
             ha="center", fontsize=9, color=C_MUTED)
    fig.tight_layout()
    _finish(fig, "fem_ov_tension_cutoff.png")


# ===========================================================================
# Figure 3 — initial lateral stress: gravity turn-on vs stated K0
# ===========================================================================

def _label_along(ax, k, gamma, depth, text, color, fontsize=9.0):
    """Write text along the line sigma_h = k*gamma*z at the given depth."""
    p0 = ax.transData.transform((0.0, 0.0))
    p1 = ax.transData.transform((k * gamma * 12.0, 12.0))
    ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
    ax.text(k * gamma * depth, depth, text, color=color, fontsize=fontsize,
            rotation=ang, rotation_mode="anchor", ha="center", va="bottom",
            zorder=5)


def fig_k0_initial():
    gamma = 19.0          # kN/m3
    z = np.linspace(0, 12, 200)
    sv = gamma * z

    fig, ax = plt.subplots(figsize=(7.6, 5.2))

    # The gravity turn-on's coefficient is nu/(1-nu): a band, not a value, because
    # nu is chosen for reasons of its own.
    k_lo, k_hi = 0.2 / 0.8, 0.4 / 0.6
    ax.fill_betweenx(z, k_lo * sv, k_hi * sv, color="#f7e6dc", zorder=0)
    ax.plot(0.3 / 0.7 * sv, z, color=C_ACCENT2, linewidth=1.8, linestyle=(0, (6, 3)))
    for k0 in (1.0, 1.5):
        ax.plot(k0 * sv, z, color=C_ACCENT, linewidth=2.0)

    ax.set_xlim(0, 400)
    ax.set_ylim(12, 0)
    ax.set_xlabel("initial lateral effective stress $\\sigma'_h$   (kPa)", fontsize=9.5)
    ax.set_ylabel("depth below level ground   (m)", fontsize=9.5)
    _tidy(ax)

    _label_along(ax, 0.3 / 0.7, gamma, 6.6,
                 "gravity turn-on,  $\\nu$ = 0.3", C_ACCENT2)
    _label_along(ax, k_hi, gamma, 3.4, "$\\nu$ = 0.4", C_ACCENT2, fontsize=8.5)
    _label_along(ax, k_lo, gamma, 10.6, "$\\nu$ = 0.2", C_ACCENT2, fontsize=8.5)
    _label_along(ax, 1.0, gamma, 7.6,
                 "$K_0$ = 1.0   compacted fill, RS2", C_ACCENT)
    _label_along(ax, 1.5, gamma, 6.4,
                 "$K_0$ = 1.5   overconsolidated clay", C_ACCENT)

    ax.text(0.985, 0.965,
            "solid: at-rest initialization — $\\sigma'_h = K_0\\sigma'_v$, stated with the soil\n"
            "dashed band: gravity turn-on — $\\sigma'_h$ follows from Poisson's ratio",
            transform=ax.transAxes, fontsize=9.0, color=C_INK, va="top", ha="right",
            bbox=dict(boxstyle="round,pad=0.45", facecolor="white",
                      edgecolor=C_GRID))

    fig.tight_layout()
    _finish(fig, "fem_ov_k0_initial.png")


# ===========================================================================
# Figure 4 — strength-reduction sweep on the Griffiths & Lane Example 1 sample
# ===========================================================================

def fig_ssrm_sweep():
    import contextlib
    import io

    from xslope.fileio import load_slope_data
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    from xslope.fem import build_fem_data, solve_fem, solve_ssrm

    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(os.path.abspath(SAMPLE))
        mesh = build_mesh_from_polygons(get_material_polygons(slope_data),
                                        target_size=6, element_type='tri6')
        fem_data = build_fem_data(slope_data, mesh)

        F_values = np.round(np.arange(1.00, 1.81, 0.05), 3)
        disp, stable = [], []
        for F in F_values:
            sol = solve_fem(fem_data, F=float(F), max_iterations=4000,
                            max_disp_factor=None)
            u = sol["displacements"] - sol["displacements_elastic"]
            disp.append(float(np.max(np.abs(u))))
            stable.append(bool(sol["stable"]))

        # A fixed global grid makes the reported factor independent of the starting
        # bracket, so the figure is reproducible from any bracket.
        ssrm = solve_ssrm(fem_data, F_min=1.0, F_max=1.8, tolerance=0.02,
                          grid=0.025, max_iterations=4000,
                          capture_failure_state=False)
        FS = float(ssrm["FS"])

    disp = np.array(disp)
    stable = np.array(stable)

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.plot(F_values, disp, color=C_MUTED, linewidth=1.2, zorder=1)
    ax.scatter(F_values[stable], disp[stable], s=52, facecolor=C_GREEN,
               edgecolor=C_GREEN, zorder=3)
    ax.scatter(F_values[~stable], disp[~stable], s=52, facecolor="white",
               edgecolor=C_RED, linewidth=1.6, zorder=3)

    ax.axvline(FS, color=C_INK, linewidth=1.3, linestyle=(0, (5, 4)), zorder=2)
    ax.annotate(f"bisection: FS = {FS:.2f}", xy=(FS, disp.max() * 0.92),
                xytext=(FS - 0.02, disp.max() * 0.92), fontsize=10,
                color=C_INK, ha="right", va="center")

    ax.text(F_values[0] + 0.01, disp.max() * 0.30,
            "filled: equilibrium reached\n(both convergence tests met)",
            fontsize=9, color=C_GREEN, ha="left", va="center")
    ax.text(F_values[-1] - 0.01, disp.max() * 0.42,
            "open: no equilibrium —\ndisplacement runs away",
            fontsize=9, color=C_RED, ha="right", va="center")

    ax.set_xlabel("strength reduction factor $F$", fontsize=9.5)
    ax.set_ylabel("maximum viscoplastic displacement   (m)", fontsize=9.5)
    _tidy(ax)
    fig.tight_layout()
    _finish(fig, "fem_ov_ssrm_sweep.png")


if __name__ == "__main__":
    fig_viscoplastic_loop()
    fig_tension_cutoff()
    fig_k0_initial()
    fig_ssrm_sweep()
