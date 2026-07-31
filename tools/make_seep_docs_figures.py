"""Render the figures for the seepage documentation pages.

Three pages are served, each figure prefixed by the page it belongs to. Every
figure is built from a committed sample input in ``docs/seep/files/`` and from
the same solver and plotting code the package ships, so a figure cannot drift
from the behaviour it illustrates.

docs/seep/overview.md
  seep_ov_bcs.png        The three boundary-condition types on the earth-dam
                         section, drawn with ``plot_seep_data``'s own symbology:
                         specified head (blue squares), exit face (red circles)
                         and specified flux (green triangles). The sample file
                         carries the head and exit-face boundaries; a rainfall
                         flux is added on the crest here so all three types are
                         visible in one figure.
  seep_ov_kr_models.png  The three relative-conductivity models as the solver
                         evaluates them (``kr_relative_vec``): the linear front
                         at several kr0/h0, and van Genuchten / Gardner curves
                         for representative soils.
  seep_ov_flownet.png    The flow net of that same earth-dam section: total-head
                         contours, stream-function flow lines, and the phreatic
                         surface, with the flow-net elements labelled.
  seep_ov_variables.png  The four contourable nodal fields of one solution
                         (total head, pore pressure, velocity magnitude,
                         hydraulic gradient magnitude).

docs/seep/transient.md   (one transient march of the earth-dam drawdown sample,
                          solved once and shared by all three figures)
  seep_tr_phreatic_lag.png   Phreatic surfaces at successive saved times on one
                             section, against the falling reservoir levels: the
                             lag that makes the answer depend on when the slope
                             is examined.
  seep_tr_storage_zones.png  Where each branch of the storage coefficient acts
                             at one drawdown frame — the saturated Ss zone, the
                             draining band, and the drained soil below it.
  seep_tr_stage_frames.png   The two frames a rapid-drawdown analysis consumes:
                             the pore-pressure field at stage_1 and at stage_2.

docs/seep/seep_slope.md
  seep_sl_u_mapping.png     Nodal pore pressures interpolated onto the slice
                            bases of one circle, and the resulting u along the
                            failure surface.
  seep_sl_piezo_vs_seep.png The same circle's base pore pressures from a
                            piezometric line and from the seepage solution.

Deterministic: fixed meshes, fixed solves, fixed frames. Run from the repo root:

    PYTHONPATH=. python3 tools/make_seep_docs_figures.py          # all
    PYTHONPATH=. python3 tools/make_seep_docs_figures.py overview # one page
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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from xslope.fileio import load_slope_data                              # noqa: E402
from xslope.mesh import (build_polygons, build_mesh_from_polygons,      # noqa: E402
                         get_material_polygons, interpolate_at_point)
from xslope.seep import (build_seep_data, build_tseep_data,             # noqa: E402
                         run_seepage_analysis, run_transient_seepage,
                         transient_frame_index, kr_relative_vec,
                         KR_LF, KR_VG, KR_GARD)
from xslope.plot_seep import plot_seep_data, plot_seep_solution         # noqa: E402
from xslope.plot import get_material_color                              # noqa: E402
from xslope.slice import generate_slices                                # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO, "docs", "seep", "files")
IMG = os.path.join(REPO, "docs", "seep", "images")

DPI = 150


def _quiet(fn, *args, **kwargs):
    """Run a chatty solver/plot call with its stdout swallowed."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return fn(*args, **kwargs)


def _save(fig, name, pad=0.12):
    """Save a figure. ``pad`` is generous where a colorbar label sits at the
    figure edge — matplotlib's tight box cuts a rotated label that overhangs the
    canvas, and the extra margin keeps it whole."""
    path = os.path.join(IMG, name)
    fig.savefig(path, dpi=DPI, bbox_inches="tight", pad_inches=pad,
                facecolor="white")
    plt.close(fig)
    print(f"  wrote {name}")


def _stack(names, out, pad=8):
    """Stack rendered panels vertically into one image (equal width)."""
    imgs = [Image.open(os.path.join(IMG, n)) for n in names]
    w = max(im.width for im in imgs)
    imgs = [im if im.width == w else
            im.resize((w, round(im.height * w / im.width)), Image.LANCZOS)
            for im in imgs]
    h = sum(im.height for im in imgs) + pad * (len(imgs) - 1)
    canvas = Image.new("RGB", (w, h), "white")
    y = 0
    for im in imgs:
        canvas.paste(im, ((w - im.width) // 2, y))
        y += im.height + pad
    canvas.save(os.path.join(IMG, out))
    for n in names:
        os.remove(os.path.join(IMG, n))
    print(f"  wrote {out}")


# =========================================================================== #
#  docs/seep/overview.md
# =========================================================================== #

DAM = "xslope_earth_dam1.xlsx"
DAM_TARGET = 1.6          # gmsh characteristic length -> ~1400-node tri3 mesh
# A rainfall flux on the dam crest, added to the sample's own head / exit-face
# boundaries so the figure can show all three boundary-condition types.
CREST_FLUX = {"coords": [(51.0, 22.0), (59.0, 22.0)], "flux": 0.02}


def _dam_model(with_flux=False):
    data = _quiet(load_slope_data, os.path.join(FILES, DAM))
    if with_flux:
        data["seepage_bc"]["specified_fluxes"] = [dict(CREST_FLUX)]
    polys = _quiet(build_polygons, data)
    mesh = _quiet(build_mesh_from_polygons, polys, DAM_TARGET, "tri3")
    seep = _quiet(build_seep_data, mesh, data)
    return data, mesh, seep


def fig_bcs():
    _, _, seep = _dam_model(with_flux=True)
    fig = plt.figure(figsize=(9.5, 3.4))
    _quiet(plot_seep_data, seep, show_bc=True, alpha=0.45, show_title=False, fig=fig)
    ax = fig.axes[0]
    ax.annotate("specified head\n(upstream pool, h = 18)", xy=(21, 9), xytext=(6, 26),
                fontsize=8, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.25"))
    ax.annotate("specified flux\n(rainfall on crest)", xy=(55, 22), xytext=(46, 33),
                fontsize=8, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.25"))
    ax.annotate("exit face\n(seepage face)", xy=(84, 11), xytext=(86, 29),
                fontsize=8, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.25"))
    ax.annotate("specified head\n(tailwater, h = 2)", xy=(107.5, 1), xytext=(106, 18),
                fontsize=8, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.25"))
    _save(fig, "seep_ov_bcs.png")


def fig_kr_models():
    psi = np.linspace(-10.0, 0.5, 600)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.6))

    ax = axes[0]
    for kr0, h0, c in [(0.01, -1.0, "#1f6fb2"), (0.01, -5.0, "#0f8a86"),
                       (0.001, -2.0, "#b2562b")]:
        kr = kr_relative_vec(psi, np.full_like(psi, kr0), np.full_like(psi, h0))
        ax.semilogy(psi, kr, color=c, lw=1.8,
                    label=f"kr0 = {kr0:g}, h0 = {h0:g}")
    ax.set_title("Linear front (unsat = lf)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left")

    ax = axes[1]
    vg = [("Sand (a = 4.42, n = 2.68)", 4.42, 2.68, "#1f6fb2"),
          ("Loam (a = 1.10, n = 1.56)", 1.10, 1.56, "#0f8a86"),
          ("Clay (a = 0.24, n = 1.09)", 0.24, 1.09, "#7a4fb2")]
    for label, a, n, c in vg:
        kr = kr_relative_vec(psi, np.zeros_like(psi), np.full_like(psi, -1.0),
                             np.full_like(psi, a), np.full_like(psi, n),
                             np.full(psi.shape, KR_VG))
        ax.semilogy(psi, kr, color=c, lw=1.8, label=label)
    kr = kr_relative_vec(psi, np.zeros_like(psi), np.full_like(psi, -1.0),
                         np.full_like(psi, 1.0), np.full_like(psi, 2.0),
                         np.full(psi.shape, KR_GARD))
    ax.semilogy(psi, kr, color="#b2562b", lw=1.8, ls="--",
                label="Gardner (a = 1, n = 2)")
    ax.set_title("van Genuchten / Gardner (unsat = vg, gard)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left")

    for ax in axes:
        ax.axvline(0.0, color="0.4", lw=0.8, ls=":")
        ax.set_xlabel("pressure head $\\psi$  (suction to the left)")
        ax.set_ylabel("relative conductivity $k_r$")
        ax.set_ylim(5e-5, 2.0)
        ax.set_xlim(psi[0], psi[-1])
        ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    _save(fig, "seep_ov_kr_models.png")


def fig_flownet_and_variables():
    _, _, seep = _dam_model()
    sol = _quiet(run_seepage_analysis, seep, tol=1e-4)

    # --- annotated flow net -------------------------------------------------
    fig = plt.figure(figsize=(9.5, 3.6))
    _quiet(plot_seep_solution, seep, sol, levels=15, base_mat=2,
           fill_contours=False, phreatic=True, mesh=False, flowlines=True, fig=fig)
    ax = fig.axes[0]
    note = dict(fontsize=8, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85),
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.25"))
    ax.annotate("phreatic surface", xy=(45.5, 17.8), xytext=(22, 20.5), **note)
    ax.annotate("total-head contour", xy=(53, 8), xytext=(24, 11), **note)
    ax.annotate("flow line", xy=(70, 5.5), xytext=(88, 17), **note)
    _save(fig, "seep_ov_flownet.png")

    # --- the four contourable fields ---------------------------------------
    panels = [("head", 'variable="head"'), ("u", 'variable="u"'),
              ("v_mag", 'variable="v_mag"'), ("i_mag", 'variable="i_mag"')]
    names = []
    for variable, label in panels:
        fig = plt.figure(figsize=(9.0, 2.6))
        _quiet(plot_seep_solution, seep, sol, variable=variable, levels=15,
               fill_contours=True, alpha=0.9, phreatic=True, mesh=False,
               flowlines=False, show_legend=False, show_title=False, fig=fig)
        ax = fig.axes[0]
        ax.text(0.005, 0.95, label, transform=ax.transAxes, fontsize=9,
                family="monospace", va="top", ha="left")
        name = f"seep_ov_var_{variable}.png"
        _save(fig, name, pad=0.3)
        names.append(name)
    _stack(names, "seep_ov_variables.png")


# =========================================================================== #
#  docs/seep/transient.md
# =========================================================================== #

TDAM = "xslope_earth_dam_tseep.xlsx"
TDAM_TARGET = 110.0 / 64.0
LAG_TIMES = [0.0, 15.0, 30.0, 47.0, 120.0, 360.0]
ZONE_TIME = 47.0


def _transient_run():
    data = _quiet(load_slope_data, os.path.join(FILES, TDAM))
    polys = _quiet(build_polygons, data)
    mesh = _quiet(build_mesh_from_polygons, polys, TDAM_TARGET, "tri3")
    seep = _quiet(build_seep_data, mesh, data)
    tseep = _quiet(build_tseep_data, data)
    sol = _quiet(run_transient_seepage, seep, tseep, verbose=False)
    return data, mesh, seep, tseep, sol


def _phreatic_xy(nodes, elements, element_types, head, xs, n_probe=240):
    """Elevation of the psi = 0 surface at each x.

    A vertical scan rather than a bisection: the domain's top is not the column's
    top on a sloping section, so the scan keeps only the probes that fall inside
    the mesh and interpolates between the deepest positive and the next negative
    pressure head.
    """
    y_lo, y_hi = float(nodes[:, 1].min()), float(nodes[:, 1].max())
    ys = np.linspace(y_hi - 1e-6, y_lo + 1e-6, n_probe)     # top down
    out = np.full(len(xs), np.nan)
    for k, x in enumerate(xs):
        prev_y = prev_p = None
        for yy in ys:
            val, found = interpolate_at_point(nodes, elements, element_types,
                                              head, (x, yy), return_found=True,
                                              signed=True)
            if not found:
                continue
            p = val - yy
            if p >= 0.0:
                if prev_p is None:                # saturated to the surface
                    out[k] = yy
                else:                             # crossing between the probes
                    f = (0.0 - prev_p) / (p - prev_p)
                    out[k] = prev_y + f * (yy - prev_y)
                break
            prev_y, prev_p = yy, p
    return out


def _dam_outline(ax, data, mesh, fill=False):
    """Section outline, optionally with the material zones washed in behind it."""
    if fill:
        for poly in get_material_polygons(data):
            xy = np.asarray(poly["coords"], dtype=float)
            ax.fill(xy[:, 0], xy[:, 1],
                    color=get_material_color(int(poly["mat_id"])), alpha=0.18,
                    lw=0, zorder=0)
    nodes = mesh["nodes"]
    gs = np.array(data["ground_surface"].coords)
    ax.plot(gs[:, 0], gs[:, 1], color="0.25", lw=1.2)
    ax.plot([nodes[:, 0].min(), nodes[:, 0].max()],
            [nodes[:, 1].min(), nodes[:, 1].min()], color="0.25", lw=1.2)


def fig_transient(run):
    data, mesh, seep, tseep, sol = run
    nodes, elements, et = mesh["nodes"], mesh["elements"], mesh["element_types"]
    times = np.asarray(sol["times"])
    series_t, series_v = tseep["series"]["pool"]

    # --- phreatic lag -------------------------------------------------------
    xs = np.linspace(nodes[:, 0].min() + 0.5, nodes[:, 0].max() - 0.5, 120)
    fig, ax = plt.subplots(figsize=(9.5, 3.4))
    _dam_outline(ax, data, mesh, fill=True)
    cmap = plt.get_cmap("viridis")
    for j, t in enumerate(LAG_TIMES):
        i = transient_frame_index(sol, t)
        head = np.asarray(sol["frames"][i]["head"])
        ph = _phreatic_xy(nodes, elements, et, head, xs)
        c = cmap(j / max(len(LAG_TIMES) - 1, 1) * 0.9)
        ax.plot(xs, ph, color=c, lw=1.8, label=f"t = {t:g} day")
        pool = float(np.interp(t, series_t, series_v))
        ax.plot([nodes[:, 0].min(), 42.0 * (pool / 18.0)], [pool, pool],
                color=c, lw=1.0, ls=":")
    ax.set_aspect("equal")
    ax.set_xlabel("x (ft)")
    ax.set_ylabel("elevation (ft)")
    ax.legend(fontsize=8, frameon=False, ncol=3, loc="upper right")
    ax.set_title("Phreatic surface through a reservoir drawdown "
                 "(dotted: pool level at that time)", fontsize=10)
    fig.tight_layout()
    _save(fig, "seep_tr_phreatic_lag.png")

    # --- storage zones ------------------------------------------------------
    i = transient_frame_index(sol, ZONE_TIME)
    head = np.asarray(sol["frames"][i]["head"])
    psi = head - nodes[:, 1]
    h0 = seep["h0_by_mat"][np.asarray(mesh["element_materials"]) - 1]
    tri = np.asarray(elements)[:, :3]
    psi_e = psi[tri].mean(axis=1)
    zone = np.where(psi_e >= 0.0, 0, np.where(psi_e > h0, 1, 2))
    colors = ["#1f6fb2", "#7fc4a8", "#e8d9b5"]
    labels = [r"saturated: $S = S_s$",
              r"draining band: $S = S_y/|h_0|$",
              r"drained: residual floor"]
    fig, ax = plt.subplots(figsize=(9.5, 3.2))
    ax.tripcolor(nodes[:, 0], nodes[:, 1], tri, facecolors=zone.astype(float),
                 cmap=matplotlib.colors.ListedColormap(colors), vmin=-0.5, vmax=2.5,
                 edgecolors="none")
    ph = _phreatic_xy(nodes, elements, et, head, xs)
    ax.plot(xs, ph, color="black", lw=2.0, label="phreatic surface")
    _dam_outline(ax, data, mesh)
    ax.set_aspect("equal")
    ax.set_xlabel("x (ft)")
    ax.set_ylabel("elevation (ft)")
    handles = [Patch(facecolor=c, label=l) for c, l in zip(colors, labels)]
    handles.append(Line2D([0], [0], color="black", lw=2.0, label="phreatic surface"))
    ax.legend(handles=handles, fontsize=8, frameon=False, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.28))
    ax.set_title(f"Storage zones at t = {ZONE_TIME:g} day", fontsize=10)
    fig.tight_layout()
    _save(fig, "seep_tr_storage_zones.png")

    # --- the two rapid-drawdown stage frames --------------------------------
    stages = [(tseep["stage_1"], "stage_1"), (tseep["stage_2"], "stage_2")]
    u_all = np.concatenate([np.asarray(sol["frames"][transient_frame_index(sol, t)]["u"])
                            for t, _ in stages])
    vmin, vmax = float(u_all.min()), float(u_all.max())
    names = []
    for t, label in stages:
        i = transient_frame_index(sol, t)
        frame = sol["frames"][i]
        fig = plt.figure(figsize=(8.6, 2.9))
        _quiet(plot_seep_solution, seep, frame, variable="u", levels=15,
               fill_contours=True, alpha=0.9, phreatic=True, mesh=False,
               flowlines=False, show_legend=False, vmin=vmin, vmax=vmax, fig=fig)
        ax = fig.axes[0]
        ax.text(0.01, 0.94, label, transform=ax.transAxes,
                fontsize=9, va="top", ha="left")
        name = f"seep_tr_stage_{label}.png"
        _save(fig, name, pad=0.3)
        names.append(name)
    _stack(names, "seep_tr_stage_frames.png")


# =========================================================================== #
#  docs/seep/seep_slope.md
# =========================================================================== #

SEEP_KEY = "xslope_rface_SEEP_KEY.xlsx"
PIEZO_KEY = "xslope_rface_PIEZO_KEY.xlsx"
NUM_SLICES = 40


def _rface_slices(xlsx):
    data = _quiet(load_slope_data, os.path.join(FILES, xlsx))
    ok, res = _quiet(generate_slices, data, circle=data["circles"][0],
                     num_slices=NUM_SLICES)
    if not ok:
        raise RuntimeError(f"generate_slices failed on {xlsx}: {res}")
    slice_df, surface = res
    return data, slice_df, surface


def fig_seep_slope():
    data, slice_df, surface = _rface_slices(SEEP_KEY)
    mesh = data["mesh"]
    u_nodal = np.asarray(data["seep_u"], dtype=float)
    xb = slice_df["x_c"].to_numpy()
    yb = slice_df["y_cb"].to_numpy()
    ub = slice_df["u"].to_numpy()

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 5.6),
                             gridspec_kw=dict(height_ratios=[1.55, 1.0]))

    ax = axes[0]
    nodes = np.asarray(mesh["nodes"])
    tri = np.asarray(mesh["elements"])[:, :3]
    top = 500.0 * np.ceil(float(np.max(u_nodal)) / 500.0)
    lev = np.linspace(0.0, top, 11)
    cf = ax.tricontourf(nodes[:, 0], nodes[:, 1], tri, np.maximum(u_nodal, 0.0),
                        levels=lev, cmap="Spectral_r", alpha=0.85)
    ax.tricontour(nodes[:, 0], nodes[:, 1], tri, u_nodal, levels=[0.0],
                  colors="black", linewidths=2.0)
    sx, sy = surface.xy
    ax.plot(sx, sy, color="0.15", lw=1.0, ls="--", zorder=5)
    for x, y_top, y_bot in zip(xb, slice_df["y_ct"].to_numpy(), yb):
        ax.plot([x, x], [y_bot, y_top], color="0.35", lw=0.4, zorder=4)
    ax.scatter(xb, yb, c=ub, cmap="Spectral_r", vmin=lev[0], vmax=lev[-1],
               s=26, edgecolor="black", linewidth=0.5, zorder=6)
    gs = np.array(data["ground_surface"].coords)
    ax.plot(gs[:, 0], gs[:, 1], color="0.2", lw=1.2)
    ax.set_aspect("equal")
    ax.set_ylabel("elevation (ft)")
    ax.set_title("Nodal pore pressures interpolated to the slice bases", fontsize=10)
    ax.legend(handles=[
        Line2D([0], [0], color="black", lw=2.0, label="phreatic surface (u = 0)"),
        Line2D([0], [0], color="0.15", lw=1.0, ls="--", label="failure surface"),
        Line2D([0], [0], color="black", lw=0, marker="o", markersize=6,
               markerfacecolor="#e0e0e0", label="slice base centre"),
    ], fontsize=8, frameon=False, loc="lower left", ncol=3,
        bbox_to_anchor=(0.0, -0.02))

    ax = axes[1]
    ax.plot(xb, ub, color="#1f6fb2", lw=1.4, marker="o", markersize=3.5)
    ax.set_xlabel("x (ft)")
    ax.set_ylabel("u at slice base (psf)")
    ax.grid(alpha=0.25, lw=0.5)

    fig.tight_layout()
    # Colorbar last, sized to the (equal-aspect) drawn height of the top panel.
    fig.canvas.draw()
    axes[0].apply_aspect()
    pos = axes[0].get_position()
    axes[1].set_xlim(axes[0].get_xlim())
    cax = fig.add_axes([pos.x1 + 0.012, pos.y0, 0.014, pos.height])
    cb = fig.colorbar(cf, cax=cax)
    cb.set_label("pore pressure u (psf)")
    _save(fig, "seep_sl_u_mapping.png", pad=0.3)

    # --- piezo vs seepage on the same circle --------------------------------
    data_p, slice_p, _ = _rface_slices(PIEZO_KEY)
    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    ax.plot(slice_p["x_c"], slice_p["u"], color="#b2562b", lw=1.6, marker="o",
            markersize=3.5, label="piezometric line")
    ax.plot(xb, ub, color="#1f6fb2", lw=1.6, marker="s", markersize=3.5,
            label="seepage solution")
    ax.set_xlabel("x (ft)")
    ax.set_ylabel("u at slice base (psf)")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=9, frameon=False)
    ax.set_title("Base pore pressures on the same circle, from the two "
                 "pore-pressure models", fontsize=10)
    fig.tight_layout()
    _save(fig, "seep_sl_piezo_vs_seep.png")


# =========================================================================== #

def main(which):
    if which in ("all", "overview"):
        print("overview.md figures")
        fig_bcs()
        fig_kr_models()
        fig_flownet_and_variables()
    if which in ("all", "transient"):
        print("transient.md figures (one transient march)")
        fig_transient(_transient_run())
    if which in ("all", "seep_slope"):
        print("seep_slope.md figures")
        fig_seep_slope()


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all")
