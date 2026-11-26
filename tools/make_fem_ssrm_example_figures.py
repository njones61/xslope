"""Build the two-layer clay FEM sample and render the strength-reduction example figures for
docs/fem/overview.md and docs/fem/solver.md.

The sample, docs/fem/files/xslope_simple_mult_layers_clay_fem.xlsx, is the limit-equilibrium sample
docs/lem/files/xslope_simple_mult_layers.xlsx (a 20 ft clay embankment, c = 400 psf, on a 10 ft clay
foundation, c = 800 psf, both phi = 0, over a firm base) with the three finite element inputs added on
the mat sheet: E = 150,000 psf (embankment) and 300,000 psf (foundation), nu = 0.40, and t_cut = 0
(no tension) on both.

Runs, one process each, cached in .figure_cache/fem_ssrm_examples/ (about two minutes each):

  base     SSRM as Studio's Run FEM makes it: tri6 at a 2 ft target size, Non-convergence, F from 1 to 2,
           tolerance 0.01, 12,000 iterations per trial. Keeps every trial's element stresses, viscoplastic
           shear strain, and largest-displacement history.
  tension  the same with t_cut blank, so the clays may carry tension; keeps every trial's stresses.
  e01      every E times 0.1.
  e10      every E times 10.
  spencer  Spencer's circular search from the file's starting circles.

Figures, in docs/fem/images/ at 200 dpi, each shown at half its pixel width:

  fem_ssrm_lem_fem.png   Spencer's critical circle above the strength-reduction shear strain.
  fem_ssrm_yielding.png  Shear stress over strength, and shear strain, at three standing trials.
  fem_ssrm_trials.png    The search's trials in order: largest displacement against F, standing or not.
  fem_ssrm_history.png   Largest displacement against iteration, one trial that stands and one that fails.
  fem_ssrm_tension.png   The soil in tension when tension is allowed, and one element's Mohr circle.

The factors of safety and displacements the docs quote are checked against EXPECTED every time the
figures are drawn; a solver change that moves one stops the script, so the pages get updated with it.

Run from the repo root:

    python tools/make_fem_ssrm_example_figures.py --sample          # write the sample workbook
    python tools/make_fem_ssrm_example_figures.py --run base        # and tension, e01, e10, spencer
    python tools/make_fem_ssrm_example_figures.py                   # figures from the cached runs
"""
import contextlib
import io
import json
import os
import pickle
import shutil
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

LEM_SAMPLE = os.path.join(ROOT, "docs", "lem", "files", "xslope_simple_mult_layers.xlsx")
SAMPLE = os.path.join(ROOT, "docs", "fem", "files", "xslope_simple_mult_layers_clay_fem.xlsx")
OUT = os.path.join(ROOT, "docs", "fem", "images")
CACHE = os.path.join(ROOT, ".figure_cache", "fem_ssrm_examples")

# mat sheet rows 11 and 12: (E psf, nu, t_cut)
FEM_INPUTS = {11: (150000.0, 0.40, 0.0), 12: (300000.0, 0.40, 0.0)}
TARGET = 2.0
SSRM = {"F_min": 1.0, "F_max": 2.0, "tolerance": 0.01, "max_iterations": 12000,
        "max_iterations_ceiling": 50000, "failure_criterion": "non_convergence",
        "capture_failure_state": True}
F_SHOWN = (1.0, 1.125)          # the two standing trials drawn beside the last one in fem_ssrm_yielding.png
F_HISTORY_STAND = 1.125         # the standing trial drawn in fem_ssrm_history.png
MOHR_AT = (47.0, 17.0)          # the element whose Mohr circle fem_ssrm_tension.png draws
FRAME = (-12.0, 72.0, -10.0, 20.0)   # the section's frame in the field figures (x0, x1, y0, y1), ft

# What the docs quote; drawing the figures checks every one (to the digits quoted)
EXPECTED = {"base_FS": 1.223, "tension_FS": 1.293, "e01_FS": 1.207, "e10_FS": 1.223, "spencer_FS": 1.245,
            "base_disp_F1": 0.153, "e01_disp_F1": 1.53, "e10_disp_F1": 0.0153,
            "at_limit": {"1.000": 13, "1.125": 22, "1.219": 40},
            "elements": 1581, "nodes": 3306}

# docs palette (tools/make_fem_overview_figures.py)
C_INK = "#1f2933"
C_MUTED = "#6b7785"
C_ACCENT = "#1f6fb2"
C_ACCENT2 = "#c05621"
C_GREEN = "#2f7d5b"
C_RED = "#b03a2e"
C_GRID = "#d7dde3"
C_PALE = "#dbe7f1"


def quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return fn(*a, **k)


# ----------------------------------------------------------------------------- the sample and its variants
def write_workbook(dst, scale_E=1.0, tcut_blank=False):
    """The LEM sample with E, nu, and t_cut added (E scaled, t_cut left blank on request)."""
    import openpyxl
    shutil.copy2(LEM_SAMPLE, dst)
    wb = openpyxl.load_workbook(dst)
    ws = wb["mat"]
    assert [ws[f"{c}10"].value for c in "LMN"] == ["t_cut", "E", "nu"], "mat sheet columns moved"
    for row, (E, nu, tcut) in FEM_INPUTS.items():
        for col in "LMN":
            v = ws[f"{col}{row}"].value
            assert not (isinstance(v, str) and v.startswith("=")), f"mat!{col}{row} holds a formula"
        ws[f"L{row}"] = None if tcut_blank else tcut
        ws[f"M{row}"] = E * scale_E
        ws[f"N{row}"] = nu
    wb.save(dst)
    return dst


def load(path):
    from xslope.fileio import load_slope_data
    return quiet(load_slope_data, path)


def build_mesh(sd):
    """The mesh Studio's Build Mesh makes at a 2 ft target size (studio/runners.py): Refine thin zones adds a
    thin-zone refinement only where it finds a thin zone, and this section has none."""
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    from studio.runners import thin_zone_refinement
    polys = get_material_polygons(sd)
    plan = thin_zone_refinement(polys, TARGET, materials=sd.get("materials"))
    kw = {"refine_factor": 3.0, "refine_features": ["thin_zones"]} if plan else {}
    return quiet(build_mesh_from_polygons, polys, target_size=TARGET, element_type="tri6", **kw)


# ----------------------------------------------------------------------------- runs
def run_ssrm(label, scale_E=1.0, tcut_blank=False):
    """One SSRM, keeping each trial's fields and the largest-displacement record solve_fem samples."""
    import xslope.fem as fem
    work = tempfile.mkdtemp(prefix=f"fem_ssrm_{label}_")
    sd = load(write_workbook(os.path.join(work, os.path.basename(SAMPLE)), scale_E, tcut_blank))
    mesh = build_mesh(sd)
    fd = quiet(fem.build_fem_data, sd, mesh)
    hist, cur, fields = [], {"h": None}, []
    ef0, sf0 = fem._early_failure, fem.solve_fem

    def ef(disp_hist, *a, **k):
        cur["h"] = disp_hist
        return ef0(disp_hist, *a, **k)

    def sf(*a, **k):
        prev, cur["h"] = cur["h"], None
        sol = sf0(*a, **k)
        F = float(k.get("F", a[1] if len(a) > 1 else np.nan))
        arr = lambda key: None if sol.get(key) is None else np.asarray(sol.get(key), dtype=float)
        fields.append({"F": F, "verdict": sol.get("verdict"), "iterations": int(sol.get("iterations", 0)),
                       "max_displacement": float(sol.get("max_displacement", np.nan)),
                       "stresses": arr("stresses"), "vp_shear_strain": arr("vp_shear_strain")})
        if cur["h"] is not None:
            hist.append({"F": F, "iterations": int(sol.get("iterations", 0)),
                         "max_displacement": float(sol.get("max_displacement", np.nan)),
                         "sample_every": int(fem._HYBRID_SAMPLE_EVERY), "disp": [float(v) for v in cur["h"]]})
        cur["h"] = prev
        return sol

    fem._early_failure, fem.solve_fem = ef, sf
    t0 = time.time()
    try:
        res = quiet(fem.solve_ssrm, fd, **SSRM)
    finally:
        fem._early_failure, fem.solve_fem = ef0, sf0
    trials = [{"F": float(t["F"]), "stable": bool(t["stable"]), "iterations": int(t["iterations"]),
               "max_displacement": float(t["max_displacement"])} for t in res["trials"]]
    out = {"label": label, "FS": float(res["FS"]), "final_interval": [float(v) for v in res["final_interval"]],
           "trials": trials, "fields": fields, "hist": hist, "mesh": mesh, "seconds": time.time() - t0,
           "scale_E": scale_E, "tcut_blank": tcut_blank}
    shutil.rmtree(work, ignore_errors=True)
    print(f"{label}: FS {out['FS']:.4f} {out['final_interval']} in {out['seconds']:.0f} s; "
          f"{len(mesh['elements'])} elements, {len(mesh['nodes'])} nodes; "
          f"trials {[(round(t['F'], 4), t['stable']) for t in trials]}")
    return out


def run_spencer():
    from xslope.search import circular_search
    sd = load(SAMPLE)
    fs_cache, conv, path, cc = quiet(circular_search, sd, "spencer", diagnostic=False)
    best = fs_cache[0]
    surf = best.get("failure_surface")
    out = {"FS": float(best["FS"]), "converged": bool(conv),
           "surface": [list(map(float, p)) for p in surf.coords] if surf is not None else None}
    print(f"spencer: FS {out['FS']:.4f}, converged {conv}")
    return out


RUNS = {"base": lambda: run_ssrm("base"), "tension": lambda: run_ssrm("tension", tcut_blank=True),
        "e01": lambda: run_ssrm("e01", scale_E=0.1), "e10": lambda: run_ssrm("e10", scale_E=10.0),
        "spencer": run_spencer}


def cached(label):
    with open(os.path.join(CACHE, f"{label}.pkl"), "rb") as f:
        return pickle.load(f)


# ----------------------------------------------------------------------------- drawing helpers
def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import xslope.plot_fem  # noqa: F401  (sets its rcParams on import: import before ours)
    plt.rcParams.update({"font.size": 12, "axes.titlesize": 12, "axes.labelsize": 12, "xtick.labelsize": 11,
                         "ytick.labelsize": 11, "legend.fontsize": 11, "axes.edgecolor": C_INK,
                         "axes.labelcolor": C_INK, "xtick.color": C_INK, "ytick.color": C_INK, "text.color": C_INK,
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    return plt


def finish(fig, name):
    import matplotlib.pyplot as plt
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=200, bbox_inches="tight", pad_inches=0.06, facecolor="white")
    plt.close(fig)
    print("wrote", os.path.relpath(path, ROOT))


def fem_data(run, sd):
    from xslope.fem import build_fem_data
    return quiet(build_fem_data, sd, run["mesh"])


def corners(mesh):
    nodes = np.asarray(mesh["nodes"], float)
    return np.asarray([e[:3] for e in mesh["elements"]], int), nodes


def outline(ax, sd, lw=1.1):
    from xslope.mesh import get_material_polygons
    for poly in get_material_polygons(sd):
        P = np.asarray(poly["coords"], dtype=float)
        P = np.vstack([P, P[:1]])
        ax.plot(P[:, 0], P[:, 1], color=C_INK, lw=lw, zorder=5)


def frame(ax, box=FRAME):
    x0, x1, y0, y1 = box
    pad = 0.035 * (x1 - x0)
    ax.set_xlim(x0 - pad, x1 + pad)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_aspect("equal")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def fields_by_F(run):
    return {round(f["F"], 7): f for f in run["fields"]}


def principal(st):
    """In-plane principal stresses and radius, compression positive."""
    sx, sy, txy = st[:, 0], st[:, 1], st[:, 2]
    R = np.sqrt(((sx - sy) / 2) ** 2 + txy ** 2)
    m = (sx + sy) / 2
    return m + R, m - R, R, m


def stress_level(st, mats, materials, F):
    """Shear stress over the reduced strength on the largest of the three Mohr circles, as the solver's yield check
    reads it: plane strain carries sigma_z, recovered from the stored deviatoric invariant (column 3) as
    |sigma_z - in-plane center| = sqrt(svm^2 - 3 R^2). Its sign is not stored; where phi > 0 needs the circle's center,
    sigma_z is taken on the side elastic plane strain puts it, toward zero from the in-plane center."""
    s1, s3, R, m = principal(st)
    dz = np.sqrt(np.maximum(st[:, 3] ** 2 - 3 * R ** 2, 0))
    out = dz > R
    sz = m - np.sign(m) * dz
    far = np.where(m >= 0, s1, s3)
    t = np.where(out, np.abs(far - sz) / 2, R)
    centre = np.where(out, (far + sz) / 2, m)
    m0 = int(mats.min())
    c = np.array([float(materials[k - m0].get("c") or 0) for k in mats]) / F
    phi = np.radians([float(materials[k - m0].get("phi") or 0) for k in mats])
    phr = np.arctan(np.tan(phi) / F)
    return np.clip(t / np.maximum(centre * np.sin(phr) + c * np.cos(phr), 1e-9), 0, 1)


def nodal(tris, n_nodes, v):
    z, n = np.zeros(n_nodes), np.zeros(n_nodes)
    for k in range(3):
        np.add.at(z, tris[:, k], v)
        np.add.at(n, tris[:, k], 1)
    return np.where(n > 0, z / np.maximum(n, 1), 0)


# ----------------------------------------------------------------------------- figures
def fig_lem_fem(plt, sd, fd, base, spencer):
    """Spencer's critical circle on the section, above the strength-reduction shear strain at the last trial that
    stood (no surface given)."""
    from matplotlib.collections import PolyCollection
    from xslope.mesh import get_material_polygons
    from xslope.plot import get_material_color
    from xslope.plot_fem import plot_shear_strain_contours
    lo = base["final_interval"][0]
    g = fields_by_F(base)[round(lo, 7)]["vp_shear_strain"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.0, 6.4))
    polys = get_material_polygons(sd)
    a1.add_collection(PolyCollection([np.asarray(p["coords"], float) for p in polys],
                                     facecolors=[get_material_color(i) for i in range(len(polys))],
                                     edgecolors=C_INK, linewidths=1.0, alpha=0.55))
    P = np.asarray(spencer["surface"])
    a1.plot(P[:, 0], P[:, 1], color=C_ACCENT2, lw=2.6, zorder=6)
    a1.text(62, 12, "embankment\nc = 400 psf, φ = 0", ha="center", va="center", fontsize=11)
    a1.text(30, -5, "foundation, c = 800 psf, φ = 0", ha="center", va="center", fontsize=11)
    a1.set_title(f"Limit equilibrium: Spencer's critical circle, F = {spencer['FS']:.3f}", loc="left")
    with contextlib.redirect_stdout(io.StringIO()):
        m, _ = plot_shear_strain_contours(a2, fd, {"vp_shear_strain": g}, False, False, single_panel=True)
    a2.set_title("")
    a2.set_xlabel("")
    a2.set_ylabel("")
    a2.set_title(f"Strength reduction, no surface given: F = {base['FS']:.3f}", loc="left")
    for ax in (a1, a2):
        frame(ax)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    for ax, mm in ((a1, None), (a2, m)):
        cax = make_axes_locatable(ax).append_axes("right", size="2.5%", pad=0.1)
        if mm is None:
            cax.axis("off")
        else:
            from matplotlib.ticker import MaxNLocator
            cb = fig.colorbar(mm, cax=cax)
            cb.locator = MaxNLocator(5)
            cb.update_ticks()
            cb.set_label("shear strain")
    fig.tight_layout(h_pad=1.2)
    finish(fig, "fem_ssrm_lem_fem.png")


def fig_yielding(plt, sd, fd, base, strain_top=0.02):
    """Rows: three standing trials. Left: shear stress over strength (orange at the limit). Right: XSLOPE's
    viscoplastic shear strain on one scale for all three, topping out at strain_top."""
    import matplotlib.tri as mtri
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from xslope.plot_fem import plot_shear_strain_contours
    tris, nodes = corners(fd)
    tri = mtri.Triangulation(nodes[:, 0], nodes[:, 1], tris)
    mats = np.asarray(fd["element_materials"], int)
    levels = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.97, 1.0001]
    cmap = ListedColormap(list(plt.get_cmap("Blues")(np.linspace(0.08, 0.75, 10))) + [C_ACCENT2])
    norm = BoundaryNorm(levels, cmap.N)
    F_all = list(F_SHOWN) + [base["final_interval"][0]]
    fbF = fields_by_F(base)
    fig, axes = plt.subplots(3, 2, figsize=(10.0, 6.6))
    at_limit, largest = {}, {}
    for r, F in enumerate(F_all):
        f = fbF[round(F, 7)]
        sl = stress_level(f["stresses"], mats, sd["materials"], F)
        key = f"{F:.3f}"
        at_limit[key] = int(round(100 * float((sl >= 0.97).mean())))
        cs = axes[r, 0].tricontourf(tri, nodal(tris, len(nodes), sl), levels=levels, cmap=cmap, norm=norm)
        g = f["vp_shear_strain"]
        largest[key] = float(g.max())
        with contextlib.redirect_stdout(io.StringIO()):
            ss, _ = plot_shear_strain_contours(axes[r, 1], fd, {"vp_shear_strain": g}, False, False,
                                               single_panel=True, vmin=0.0, vmax=strain_top)
        gm = g.max()
        axes[r, 0].set_title(f"F = {F:.3f}: {at_limit[key]}% of elements at the limit", loc="left")
        axes[r, 1].set_title("")
        axes[r, 1].set_title(f"F = {F:.3f}: largest shear strain {gm:.3f}" if gm < 0.1
                             else f"F = {F:.3f}: largest shear strain {gm:.2f}", loc="left")
        for ax in axes[r]:
            ax.set_xlabel("")
            ax.set_ylabel("")
            outline(ax, sd)
            frame(ax)
            if r < 2:
                ax.tick_params(labelbottom=False)
    fig.tight_layout(h_pad=1.6, w_pad=1.5)
    fig.subplots_adjust(bottom=0.17)
    for col, (mappable, label, ticks, labels) in enumerate((
            (cs, "shear stress / strength", [0, 0.5, 0.9, 1.0], ["0", "0.5", "0.9", "1"]),
            (ss, "shear strain", [0, strain_top / 2, strain_top], ["0", f"{strain_top / 2:g}", f"{strain_top:g}+"]))):
        pos = axes[2, col].get_position()
        cax = fig.add_axes([pos.x0 + 0.1 * pos.width, 0.06, 0.8 * pos.width, 0.018])
        cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", ticks=ticks)
        cb.ax.set_xticklabels(labels)
        cb.set_label(label)
    finish(fig, "fem_ssrm_yielding.png")
    return at_limit, largest


def fig_trials(plt, base):
    """The search's trials in order: largest displacement against F. A standing trial's point is where it came to
    rest; a failing trial's is where it was when stopped, still moving (the arrow)."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    lo, hi = base["final_interval"]
    FS = base["FS"]
    T = [(t["F"], t["stable"], t["max_displacement"]) for t in base["trials"]]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.0), gridspec_kw={"width_ratios": [1.25, 1]})
    zoom = (1.175, 1.26)
    for ax, (x0, x1) in ((a1, (0.95, 2.08)), (a2, zoom)):
        ax.axvspan(lo, hi, color=C_GRID, lw=0, zorder=0)
        for i, (F, ok, d) in enumerate(T, 1):
            if not (x0 <= F <= x1):
                continue
            if ok:
                ax.plot(F, d, "o", color=C_GREEN, ms=8, zorder=3)
            else:
                ax.plot(F, d, "o", mfc="white", mec=C_RED, mew=1.8, ms=8, zorder=3)
                ax.annotate("", xy=(F, d + 0.30), xytext=(F, d + 0.07),
                            arrowprops=dict(arrowstyle="-|>", color=C_RED, lw=1.4))
            if ax is a2 or not (zoom[0] <= F <= zoom[1]):
                dx = (x1 - x0) * 0.018
                left = ax is a2 and ok and F > lo - 1e-9
                if ax is a1:
                    ax.text(F, d + (0.36 if not ok else 0.10), str(i), fontsize=10, color=C_MUTED,
                            ha="center", va="bottom")
                else:
                    ax.text(F - dx if left else F + dx, d - 0.05, str(i), fontsize=10, color=C_MUTED,
                            ha="right" if left else "left", va="top")
        ax.set_xlim(x0, x1)
        ax.set_ylim(0, 2.75)
        ax.grid(color=C_GRID, lw=0.6)
        ax.set_xlabel("trial F")
    a1.set_ylabel("largest displacement (ft)")
    a1.add_patch(Rectangle((zoom[0], 0.02), zoom[1] - zoom[0], 2.0, fill=False, ec=C_MUTED, lw=0.9, ls=":", zorder=4))
    a1.text(zoom[1] + 0.02, 1.95, "close up →", fontsize=10, color=C_MUTED, ha="left", va="center")
    a2.axvline(FS, color=C_INK, ls="--", lw=1.2)
    a2.text(FS + 0.002, 2.6, f"FS = {FS:.3f}", color=C_INK, ha="left", va="top")
    a2.set_xticks([1.18, 1.20, 1.22, 1.24, 1.26])
    a1.set_title(f"All {len(T)} trials, numbered in order", loc="left")
    a2.set_title(f"The final bracket, {lo:.4f} to {hi:.4f}", loc="left")
    a1.legend([Line2D([], [], ls="", marker="o", color=C_GREEN, ms=8),
               Line2D([], [], ls="", marker="o", mfc="white", mec=C_RED, mew=1.8, ms=8)],
              ["stands: where it came to rest", "fails: where it was stopped,\nstill moving"],
              loc="lower right", bbox_to_anchor=(1.0, 0.06), handletextpad=0.3)
    fig.tight_layout(w_pad=2.0)
    finish(fig, "fem_ssrm_trials.png")


def fig_history(plt, base):
    """Largest displacement against iteration: a trial that stands levels off; the first trial above the factor of
    safety is still moving when its 12,000 iterations run out."""
    lo, hi = base["final_interval"]
    H = {round(h["F"], 7): h for h in base["hist"]}
    hs, hf = H[round(F_HISTORY_STAND, 7)], H[round(hi, 7)]
    assert abs(hs["disp"][-1] - hs["max_displacement"]) < 1e-3, "the standing trial's record ends off its final value"
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    xf = (np.arange(len(hf["disp"])) + 1) * hf["sample_every"]
    ax.plot(xf, hf["disp"], color=C_RED, lw=2.0)
    ax.plot([hf["iterations"]], [hf["max_displacement"]], "o", mfc="white", mec=C_RED, mew=1.8, ms=8)
    xs = (np.arange(len(hs["disp"])) + 1) * hs["sample_every"]
    ax.plot(xs, hs["disp"], color=C_GREEN, lw=2.0)
    ax.plot([hs["iterations"]], [hs["max_displacement"]], "o", color=C_GREEN, ms=8)
    ax.set_xscale("log")
    ax.set_xlim(8, 16000)
    ax.set_ylim(0, 0.8)
    ax.set_xticks([10, 100, 1000, 10000])
    ax.set_xticklabels(["10", "100", "1,000", "10,000"])
    ax.minorticks_off()
    ax.set_xlabel("iteration (log scale)")
    ax.set_ylabel("largest displacement (ft)")
    ax.grid(color=C_GRID, lw=0.6)
    ax.text(hs["iterations"] * 1.25, hs["max_displacement"] - 0.03,
            f"F = {F_HISTORY_STAND:.3f}: levels off,\nequilibrium at {hs['iterations']} iterations",
            color=C_GREEN, va="top", ha="left")
    ax.text(hf["iterations"] * 0.85, hf["max_displacement"] + 0.05,
            f"F = {hi:.3f}: still moving\nafter {hf['iterations']:,} iterations", color=C_RED, va="bottom", ha="right")
    fig.tight_layout()
    finish(fig, "fem_ssrm_history.png")
    return {"stand_iterations": hs["iterations"], "fail_iterations": hf["iterations"]}


def fig_tension(plt, sd, fd, tension):
    """Left: tension allowed (t_cut blank), the last trial that stood, the elements whose minor principal stress is
    tensile in red. Right: the Mohr circle of the element nearest MOHR_AT under the flat phi = 0 line, with the
    cutoff a t_cut of zero would put at zero normal stress."""
    from matplotlib.collections import PolyCollection
    from matplotlib.patches import Rectangle
    from xslope.plot import get_material_color
    lo = tension["final_interval"][0]
    st = fields_by_F(tension)[round(lo, 7)]["stresses"]
    s1, s3, R, m = principal(st)
    ten = s3 < -1.0
    tris, nodes = corners(fd)
    cen = nodes[tris].mean(axis=1)
    k = int(np.argmin(np.hypot(cen[:, 0] - MOHR_AT[0], cen[:, 1] - MOHR_AT[1])))
    mats = np.asarray(fd["element_materials"], int)
    m0 = mats.min()
    c_emb = float(sd["materials"][int(mats[k] - m0)]["c"])
    cF = c_emb / lo
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 3.6), gridspec_kw={"width_ratios": [1.55, 1]})
    base_cols = [get_material_color(int(t - m0)) for t in mats]
    cols = [C_RED if t else (*bc[:3], 0.35) for t, bc in zip(ten, base_cols)]
    a1.add_collection(PolyCollection([nodes[t] for t in tris], facecolors=cols, edgecolors="none"))
    outline(a1, sd)
    a1.plot(*cen[k], marker="o", mfc="white", mec=C_INK, mew=1.4, ms=6, zorder=7)
    a1.annotate("the circle's element", xy=cen[k], xytext=(cen[k][0] + 6, cen[k][1] - 8), fontsize=10,
                arrowprops=dict(arrowstyle="-", color=C_INK, lw=0.8, shrinkA=2, shrinkB=5))
    frame(a1)
    a1.set_title(f"Tension allowed, F = {lo:.3f}: the soil in tension (red)", loc="left")
    # Mohr circle: one scale on both axes so the circle is round
    lo_s, hi_s = -450.0, 600.0
    a2.add_patch(Rectangle((lo_s, 0), -lo_s, 1.45 * cF, color=C_PALE, lw=0, zorder=0))
    a2.text(lo_s + 15, 1.38 * cF, "tension, σ′ < 0", fontsize=10, color=C_MUTED, va="top")
    a2.axhline(cF, color=C_INK, lw=2.0, zorder=2)
    a2.text(hi_s - 10, cF + 12, f"τ$_f$ = c/F = {cF:.0f} psf (φ = 0)", fontsize=10, ha="right", va="bottom")
    th = np.linspace(0, np.pi, 121)
    ctr, rad = (s1[k] + s3[k]) / 2, (s1[k] - s3[k]) / 2
    a2.plot(ctr + rad * np.cos(th), rad * np.sin(th), color=C_ACCENT2, lw=2.2, zorder=3)
    a2.plot([s3[k], s1[k]], [0, 0], "o", color=C_INK, ms=5, zorder=4)
    a2.text(s3[k], -22, f"σ′$_3$ = {s3[k]:.0f}".replace("-", "−"), ha="center", va="top", fontsize=10)
    a2.text(s1[k], -22, f"σ′$_1$ = {s1[k]:.0f}".replace("-", "−"), ha="center", va="top", fontsize=10)
    a2.plot([0, 0], [0, 1.45 * cF], color=C_RED, lw=2.6, zorder=3)
    a2.text(12, 1.38 * cF, "t$_{cut}$ = 0", color=C_RED, fontsize=10, va="top", fontweight="bold")
    a2.axhline(0, color=C_INK, lw=0.9)
    a2.set_xlim(lo_s, hi_s)
    a2.set_ylim(-75, 1.5 * cF)
    a2.set_aspect("equal")
    a2.set_xlabel("σ′ (psf, compression positive)")
    a2.set_ylabel("τ (psf)")
    a2.set_title("Its Mohr circle", loc="left")
    fig.tight_layout(w_pad=2.0)
    finish(fig, "fem_ssrm_tension.png")
    return {"F_shown": lo, "cF": cF, "s1": float(s1[k]), "s3": float(s3[k]), "x": float(cen[k][0]),
            "y": float(cen[k][1]), "tension_max": float(-s3.min()), "tension_fraction": float(ten.mean())}


# ----------------------------------------------------------------------------- checks
def check(name, got, want):
    digits = len(str(want).split(".")[1]) if "." in str(want) else 0
    if round(got, digits) != round(want, digits):
        raise SystemExit(f"{name}: {got:.6g} no longer rounds to the docs' {want}: update the pages and EXPECTED")


def disp_at(run, F):
    return next(t["max_displacement"] for t in run["trials"] if abs(t["F"] - F) < 1e-9)


def figures():
    plt = style()
    base, tension, e01, e10 = (cached(k) for k in ("base", "tension", "e01", "e10"))
    spencer = cached("spencer")
    sd = load(SAMPLE)
    fd = fem_data(base, sd)
    check("elements", len(base["mesh"]["elements"]), EXPECTED["elements"])
    check("nodes", len(base["mesh"]["nodes"]), EXPECTED["nodes"])
    for name, run in (("base", base), ("tension", tension), ("e01", e01), ("e10", e10)):
        check(f"{name}_FS", run["FS"], EXPECTED[f"{name}_FS"])
    check("spencer_FS", spencer["FS"], EXPECTED["spencer_FS"])
    for name, run in (("base", base), ("e01", e01), ("e10", e10)):
        check(f"{name}_disp_F1", disp_at(run, 1.0), EXPECTED[f"{name}_disp_F1"])
    fig_lem_fem(plt, sd, fd, base, spencer)
    at_limit, largest = fig_yielding(plt, sd, fd, base)
    for k, v in EXPECTED["at_limit"].items():
        check(f"at_limit {k}", at_limit[k], v)
    fig_trials(plt, base)
    hist = fig_history(plt, base)
    tsd = load(write_workbook(os.path.join(tempfile.mkdtemp(), "t.xlsx"), tcut_blank=True))
    ten = fig_tension(plt, tsd, fem_data(tension, tsd), tension)
    summary = {"base": {"FS": base["FS"], "interval": base["final_interval"], "trials": base["trials"]},
               "tension": {"FS": tension["FS"], **ten}, "e01": {"FS": e01["FS"], "disp_F1": disp_at(e01, 1.0)},
               "e10": {"FS": e10["FS"], "disp_F1": disp_at(e10, 1.0)}, "base_disp_F1": disp_at(base, 1.0),
               "spencer": spencer["FS"], "at_limit": at_limit, "largest_strain": largest, "history": hist}
    print(json.dumps(summary, indent=1, default=float))


def main():
    if "--sample" in sys.argv:
        print("wrote", os.path.relpath(write_workbook(SAMPLE), ROOT))
        return
    if "--run" in sys.argv:
        label = sys.argv[sys.argv.index("--run") + 1]
        out = RUNS[label]()
        os.makedirs(CACHE, exist_ok=True)
        with open(os.path.join(CACHE, f"{label}.pkl"), "wb") as f:
            pickle.dump(out, f, protocol=4)
        return
    figures()


if __name__ == "__main__":
    main()
