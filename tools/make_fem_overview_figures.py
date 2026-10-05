"""Render the illustration figures for docs/fem/overview.md.

Three figures land in docs/fem/images/ with the fem_ov_ prefix:

  fem_ov_viscoplastic_loop.png  Flow diagram of the viscoplastic algorithm: one
                                factorization, the per-Gauss-point yield check, the
                                body-load correction, and the two convergence tests
                                with the hybrid verdict that follows them.
  fem_ov_k0_initial.png         Initial lateral effective stress with depth: the
                                gravity turn-on's nu/(1-nu) coefficient against
                                stated K0 values.
  fem_ov_ssrm_sweep.png         Strength-reduction sweep on the Griffiths & Lane
                                Example 1 sample file — viscoplastic displacement
                                against F, marked by whether the trial reached
                                equilibrium, with the bisection's factor of safety.

The first two are closed-form schematics. The third runs the solver on the
committed sample docs/fem/files/xslope_griffiths1.xlsx (about a minute).

The tension-cutoff schematic is exported from the private repo's
drawings/docs/fem/tension_cutoff/ source; this script never overwrites it.

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
FEM01_CURVE_TOLERANCE = 0.02  # Overview figure only; tutorial tolerance is unchanged.

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
    gamma = 19.0          # kN/m3: same uniform dry column as the previous figure
    z = np.linspace(0, 12, 200)
    sv = gamma * z

    fig = plt.figure(figsize=(4.9, 6.1))
    gs = fig.add_gridspec(2, 2, height_ratios=(1.1, 2.0),
                          left=0.12, right=0.98, bottom=0.19, top=0.91,
                          wspace=0.36, hspace=0.24)

    def element(ax, title, relation, color):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_title(title, fontsize=8.3, color=C_INK, pad=9)
        ax.add_patch(plt.Rectangle((0.1, 0.02), 0.8, 0.88,
                                   facecolor="#e5d2ad", alpha=0.3, edgecolor="none"))
        ax.plot([0.1, 0.9], [0.9, 0.9], color=C_INK, linewidth=0.9)
        ax.add_patch(plt.Rectangle((0.39, 0.34), 0.26, 0.24,
                                   facecolor="white", edgecolor=C_INK, linewidth=0.9))
        for tail, tip in (((0.52, 0.77), (0.52, 0.59)),
                          ((0.52, 0.14), (0.52, 0.33))):
            ax.annotate("", xy=tip, xytext=tail,
                        arrowprops=dict(arrowstyle="-|>", color=C_INK, linewidth=0.9))
        for tail, tip in (((0.17, 0.46), (0.38, 0.46)),
                          ((0.87, 0.46), (0.66, 0.46))):
            ax.annotate("", xy=tip, xytext=tail,
                        arrowprops=dict(arrowstyle="-|>", color=color, linewidth=0.9))
        ax.text(0.52, 0.81, r"$\sigma'_v$ from overburden",
                fontsize=7.5, ha="center", color=C_INK)
        ax.text(0.84, 0.53, r"$\sigma'_h$", fontsize=8, ha="center", color=color)
        ax.annotate("", xy=(0.1, 0.46), xytext=(0.1, 0.9),
                    arrowprops=dict(arrowstyle="<->", color=C_INK, linewidth=0.7))
        ax.text(0.045, 0.66, "$z$", fontsize=8, ha="right", color=C_INK)
        ax.text(0.52, -0.03, relation, fontsize=8.3, ha="center", color=color)

    element(fig.add_subplot(gs[0, 0]),
            "Gravity turn-on\nElastic, zero lateral strain",
            r"$\sigma'_h=\frac{\nu}{1-\nu}\,\sigma'_v$", C_ACCENT2)
    element(fig.add_subplot(gs[0, 1]),
            "At-rest initialization\nSpecified $K_0$",
            r"$\sigma'_h=K_0\,\sigma'_v$", C_ACCENT)

    for column in (0, 1):
        ax = fig.add_subplot(gs[1, column])
        ax.plot(sv, z, color=C_INK, linewidth=1.5, label=r"$|\sigma'_v|$")
        if column == 0:
            k_lo, k_hi = 0.2 / 0.8, 0.4 / 0.6
            ax.fill_betweenx(z, k_lo * sv, k_hi * sv, color="#f7e6dc", zorder=0,
                             label=r"$\nu=0.2$–$0.4$")
            ax.plot(0.3 / 0.7 * sv, z, color=C_ACCENT2, linewidth=1.4,
                    linestyle=(0, (6, 3)), label=r"$|\sigma'_h|,\ \nu=0.3$")
        else:
            for k0, style in ((0.5, (0, (5, 3))), (1.5, "solid")):
                ax.plot(k0 * sv, z, color=C_ACCENT, linewidth=1.1, linestyle=style,
                        label=rf"$|\sigma'_h|,\ K_0={k0:.1f}$")
        ax.set_xlim(0, 400)
        ax.set_ylim(12, 0)
        ax.set_xticks((0, 200, 400))
        ax.set_yticks((0, 4, 8, 12))
        ax.set_xlabel("Stress magnitude (kPa)", fontsize=7.5)
        if column == 0:
            ax.set_ylabel("Depth below level ground (m)", fontsize=7.5)
        _tidy(ax)
        ax.tick_params(labelsize=7.2)
        ax.legend(loc="upper right", fontsize=6.8, framealpha=1, edgecolor=C_GRID,
                  borderpad=0.4, handlelength=1.5, labelspacing=0.45)

    fig.text(0.31, 0.105, r"Gravity band: $\nu=0.2$–$0.4$",
             ha="center", fontsize=7.4, color=C_INK)
    fig.text(0.31, 0.075, r"Shown line: $\nu=0.3$",
             ha="center", fontsize=7.4, color=C_INK)
    fig.text(0.55, 0.025, r"$\gamma=19$ kN/m$^3$; both vertical profiles are the same.",
             ha="center", fontsize=7.5, color=C_INK)
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
            # early_failure off: this curve is the displacement each trial
            # REACHES, so a failing trial has to run its budget to reach it.
            sol = solve_fem(fem_data, F=float(F), max_iterations=4000,
                            max_disp_factor=None, early_failure=False)
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


def fig_capture_comparison():
    """Reload Tutorial W-3's saved standing and failure fields; never solve."""
    from io import BytesIO
    import json
    from PIL import Image
    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data, import_fem_solution
    from xslope.plot_fem import plot_fem_results

    stem = os.path.join(HERE, "..", "docs", "tutorials", "files",
                        "xslope_johnson_res_solved")
    data = load_slope_data(stem + ".xlsx")
    fem_data = build_fem_data(data)
    standing = import_fem_solution(fem_data, stem)
    failure = standing.get("failure_solution")
    if not standing.get("converged") or not failure or failure.get("converged"):
        raise ValueError("The committed tutorial must hold both distinct field states")
    with open(stem + "_fem_meta.json") as stream:
        meta = json.load(stream)
    panels = []
    for state in ("converged", "failure"):
        fig, ax = plot_fem_results(
            fem_data, standing, plot_type="deformation", field_state=state,
            failure_solution=failure, fs=meta["FS"], ssrm_record=meta)
        print(f"{state}: {ax.get_title()}", flush=True)
        buffer = BytesIO()
        fig.savefig(buffer, dpi=200, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        buffer.seek(0)
        panels.append(Image.open(buffer).convert("RGB"))

    # Preserve each standard render pixel for pixel. Only add white padding
    # and a separator; never replace titles, legends or automatic scales.
    gap = 16
    combined = Image.new("RGB", (max(p.width for p in panels),
                                  sum(p.height for p in panels) + gap), "white")
    y = 0
    for panel in panels:
        x = (combined.width - panel.width) // 2
        combined.paste(panel, (x, y))
        assert combined.crop((x, y, x + panel.width, y + panel.height)).tobytes() == panel.tobytes()
        y += panel.height + gap
    path = os.path.abspath(os.path.join(OUT, "fem_capture_comparison.png"))
    combined.save(path)
    print("wrote", path, "from saved W-3 fields; no solve", flush=True)


def run_fem1_curve():
    """Figure-only FEM-1 search at the wider tolerance; launch through the gate."""
    import json
    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data, solve_ssrm
    from tools.make_tutorial_figures import (
        FEM01_DONE, FEM01_CRITERION, FEM01_F_MIN, FEM01_F_MAX,
        FEM01_MAX_ITERATIONS, _fem01_mesh)

    tolerance = FEM01_CURVE_TOLERANCE
    model = load_slope_data(FEM01_DONE)
    fem_data = build_fem_data(model, _fem01_mesh(model))
    result = solve_ssrm(
        fem_data, F_min=FEM01_F_MIN, F_max=FEM01_F_MAX, tolerance=tolerance,
        max_iterations=FEM01_MAX_ITERATIONS, failure_criterion=FEM01_CRITERION,
        capture_failure_state=False, debug_level=1)
    record = {key: result.get(key) for key in (
        "FS", "final_interval", "interval_width", "trials", "fs_is_lower_bound", "summary")}
    record.update(benchmark="FEM-1-overview-tol-0.02", analysis="ssrm",
                  file="docs/tutorials/files/xslope_ssrm_embankment.xlsx",
                  tolerance=tolerance, F_min=FEM01_F_MIN, F_max=FEM01_F_MAX,
                  max_iter=FEM01_MAX_ITERATIONS, failure_criterion=FEM01_CRITERION,
                  target_size=model["target_size"], element_type=model["element_type"],
                  unit_system=model.get("unit_system"), figure_only=True)
    path = os.path.join(OUT, "fem01_ssrm_curve_fem_meta.json")
    with open(path, "w") as stream:
        json.dump(record, stream, indent=2,
                  default=lambda value: value.tolist() if isinstance(value, np.ndarray)
                  else value.item() if isinstance(value, np.generic) else str(value))
        stream.write("\n")
    print("wrote figure-only record", path, flush=True)
    fig_displacement_curves()


def fig_displacement_curves():
    """Replay only the figure-only FEM-1 record, without a new solve."""
    import json
    from xslope.fileio import load_slope_data
    from xslope.plot_fem import plot_ssrm_curve

    model = load_slope_data(os.path.join(HERE, "..", "docs", "tutorials", "files",
                                        "xslope_ssrm_embankment.xlsx"))
    with open(os.path.join(OUT, "fem01_ssrm_curve_fem_meta.json")) as stream:
        record = json.load(stream)
    fig, ax = plt.subplots(figsize=(9, 5))
    plot_ssrm_curve(ax, record, fem_data=model)
    fig.tight_layout()
    print("FEM-1", ax.get_title(), ax.get_ylabel(),
          [t.get_text() for t in ax.get_legend().get_texts()], flush=True)
    path = os.path.abspath(os.path.join(OUT, "fem01_ssrm_curve.png"))
    fig.savefig(path, dpi=150, facecolor="white")
    plt.close(fig)
    print("wrote", path, "from saved record; no solve", flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["run-fem1-curve"]:
        run_fem1_curve()
    elif sys.argv[1:] == ["displacement-curves"]:
        fig_displacement_curves()
    elif sys.argv[1:] == ["capture-comparison"]:
        fig_capture_comparison()
    elif sys.argv[1:]:
        raise SystemExit("usage: make_fem_overview_figures.py "
                         "[run-fem1-curve|displacement-curves|capture-comparison]")
    else:
        fig_viscoplastic_loop()
        fig_k0_initial()
        fig_ssrm_sweep()
