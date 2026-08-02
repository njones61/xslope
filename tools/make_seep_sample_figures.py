"""Regenerate the steady seepage flow-net figures for docs/seep/samples.md.

These six ``docs/seep/images/<name>_solution.png`` panels had no committed
generator (they were one-off manual saves), so they drifted from the page's
locked flowrates as the solver evolved. This script is their deterministic,
committed source: each figure is produced by the *same* mesh-and-solve code path
that ``run_tests.py::run_seep_test`` uses to compute the page's ``type=seep``
lock, so the flowrate printed on the figure's subtitle is — by construction —
the locked value.

Solve configuration (identical to run_seep_test's defaults, which is what every
``type=seep`` lock tag for these files was recorded from):

    element_type = 'tri3'
    target_size  = (ground-surface x-extent) / 120
    build_seep_data(mesh, slope_data)          # default seep_bc=1
    run_seepage_analysis(seep, tol=1e-4, max_iter=400)

Figure sizing (one coherent rule for all six, per the house frame spec —
equal aspect, uniform cushion, full-height colorbar, legend in reserved space,
all supplied unchanged by plot_seep_solution): a *fixed figure width* and *dpi*
so every panel comes out the same pixel width, with the height following each
domain's own aspect ratio at equal scale. ``bbox_inches='tight'`` then crops the
height to the domain-driven content, so a wide-thin blanket renders short and a
compact dam renders taller — but all at one consistent width and scale, instead
of the old artifacts' 1411x338…3634x1462 grab-bag.

Run from the repo root:

    PYTHONPATH=. python3 tools/make_seep_sample_figures.py
"""
import io
import os
import contextlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import build_seep_data, run_seepage_analysis
from xslope.plot_seep import plot_seep_solution

OUT = "docs/seep/images"

# Coherent sizing rule (see module docstring): a fixed figure width in inches and
# a fixed dpi give every panel the same pixel width; the height tracks the domain
# aspect so the flow net stays at one consistent scale across all six figures.
FIG_W_IN = 11.0     # figure width (inches) -> ~FIG_W_IN*DPI px wide after tight crop
AXES_W_FRAC = 0.80  # fraction of the width the equal-aspect axes span (rest: y-axis
                    # margin on the left, colorbar + label on the right)
EXTRA_H_IN = 2.05   # headroom the two-line title + below-axes legend + margins need
MIN_H_IN = 2.6      # floor so a very flat domain (clay blanket) still lays out cleanly
DPI = 200

# (file stem, phreatic, base_mat):
#   phreatic=True draws the free surface for the unconfined (partially saturated)
#   dams; phreatic=False for the fully saturated / confined problems, where the
#   zero-pressure contour has no physical meaning.
#   base_mat is the reference material for the flow-net's flow-channel count
#   (plot_seep_solution scales the number of flow lines to q/(k_base·Δh) so the
#   cells read as curvilinear squares in that zone). It does NOT affect the
#   flowrate — only the flow-line density. Choose the zone that controls the
#   through-flow so the net is readable: the low-k core for the zoned dams, the
#   foundation for Johnson (whose near-cutoff core, k≈0.001, would otherwise
#   demand hundreds of flow lines), and the single/through-flow layer otherwise.
#   Matches the flow-net densities of the pre-existing (Norm-approved) figures.
CASES = [
    ("clay_blanket",  False, 1),   # saturated: sheetpile + clay blanket (1 material)
    ("sea_trench",    False, 1),   # saturated: trench, silt is the through-flow layer
    ("earth_dam1",    True,  2),   # partially saturated: shell + core; base = core
    ("earth_dam1_vg", True,  2),   # same dam, van Genuchten model; base = core
    ("johnson_res",   True,  3),   # shell/core/foundation; base = foundation
    ("earth_dam2",    True,  3),   # shell/filter/core; base = core
]


def _figsize(nodes):
    """One coherent size for every panel: fixed width, height from the domain's
    equal-aspect ratio (computed from the meshed node extents, which are what
    plot_seep_solution frames)."""
    dom_w = float(nodes[:, 0].max() - nodes[:, 0].min())
    dom_h = float(nodes[:, 1].max() - nodes[:, 1].min())
    axes_w_in = FIG_W_IN * AXES_W_FRAC
    axes_h_in = axes_w_in * (dom_h / dom_w)
    return (FIG_W_IN, max(MIN_H_IN, axes_h_in + EXTRA_H_IN))


def seep_sample(stem, phreatic, base_mat):
    xlsx = f"docs/seep/files/xslope_{stem}.xlsx"
    slope_data = load_slope_data(xlsx)

    # Mesh + solve EXACTLY as run_tests.py::run_seep_test does, so the figure's
    # flowrate is the page's locked value.
    x_coords = [x for x, _ in slope_data["ground_surface"].coords]
    target_size = (max(x_coords) - min(x_coords)) / 120
    with contextlib.redirect_stdout(io.StringIO()):
        polygons = get_material_polygons(slope_data)
        mesh = build_mesh_from_polygons(polygons, target_size, "tri3")
        seep_data = build_seep_data(mesh, slope_data)
        # 1000, not the solver default of 400: earth_dam2 needs 427 sweeps of the
        # exit face to settle on this mesh, and its test tag carries the same
        # ceiling. Every other sample converges in far fewer, so the higher cap
        # changes no figure but the one that would otherwise not be drawn.
        solution = run_seepage_analysis(seep_data, tol=1e-4, max_iter=1000)

    if not solution.get("converged", True):
        raise RuntimeError(f"{stem}: seepage solve did not converge")

    figsize = _figsize(mesh["nodes"])
    path = os.path.join(OUT, f"{stem}_solution.png")

    # mesh=False: these are flow-net figures — element edges muddy the head fill
    # and chop the contour/flow lines into a dashed look (same convention as
    # benchmarks/make_benchmark_figures.py::seep_case and make_gw_figures.py).
    saved = []

    def _show(*a, **k):
        plt.gcf().savefig(path, dpi=DPI, bbox_inches="tight")
        saved.append(path)
        plt.close("all")

    orig = plt.show
    plt.show = _show
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            plot_seep_solution(seep_data, solution, figsize=figsize, levels=20,
                               base_mat=base_mat, fill_contours=True,
                               phreatic=phreatic, mesh=False)
        if not saved:
            plt.gcf().savefig(path, dpi=DPI, bbox_inches="tight")
            plt.close("all")
    finally:
        plt.show = orig

    q = solution["flowrate"]
    print(f"  {stem}: q = {q:.4f}  ->  {path}")
    return q


if __name__ == "__main__":
    for stem, phreatic, base_mat in CASES:
        seep_sample(stem, phreatic, base_mat)
    print("\nseep sample figures generated")
