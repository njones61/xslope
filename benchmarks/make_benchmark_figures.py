"""Generate the docs images for the verification benchmark sample problems.

Renders input-geometry and solution figures for the benchmark models into
docs/{lem,seep,fem}/images/, for the corresponding samples.md entries.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/make_benchmark_figures.py
"""
import io
import contextlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from xslope.fileio import load_slope_data
from xslope.mesh import get_material_polygons, build_mesh_from_polygons
from xslope.seep import build_seep_data, run_seepage_analysis
from xslope.search import circular_search, noncircular_search
from xslope.plot import plot_inputs, plot_solution
from xslope.plot_seep import plot_seep_solution


def capture(path, fn, *args, **kwargs):
    """Run a plotting function, saving whatever it shows to `path`."""
    saved = []

    def _show(*a, **k):
        fig = plt.gcf()
        fig.savefig(path, dpi=200, bbox_inches="tight")
        saved.append(path)
        plt.close("all")

    orig = plt.show
    plt.show = _show
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn(*args, **kwargs)
        if not saved:  # function didn't call show(); save current figure
            plt.gcf().savefig(path, dpi=200, bbox_inches="tight")
            plt.close("all")
    finally:
        plt.show = orig
    print("wrote", path)


def lem_case(xlsx, img_prefix, surface):
    sd = load_slope_data(xlsx)
    capture(f"docs/lem/images/{img_prefix}_inputs.png", plot_inputs, sd, mode="lem")
    with contextlib.redirect_stdout(io.StringIO()):
        if surface == "circular":
            fs_cache, _, _, _ = circular_search(sd, "spencer", num_slices=50, diagnostic=False)
        else:
            fs_cache, _, _ = noncircular_search(sd, "spencer", num_slices=50, diagnostic=False)
    crit = fs_cache[0]
    capture(f"docs/lem/images/{img_prefix}_solution.png", plot_solution,
            sd, crit["slices"], crit["failure_surface"], crit["solver_result"])
    print(f"  {img_prefix}: Spencer FS = {crit['solver_result']['FS']:.3f}")


def seep_case(xlsx, img_prefix, target_size, base_mat=1, phreatic=True):
    # phreatic=False for confined problems: the zero-pressure contour drawn as a
    # "phreatic surface" has no physical meaning in fully saturated confined flow.
    sd = load_slope_data(xlsx)
    with contextlib.redirect_stdout(io.StringIO()):
        polys = get_material_polygons(sd)
        mesh = build_mesh_from_polygons(polys, target_size=target_size, element_type="tri6")
        seep = build_seep_data(mesh, sd, seep_bc=1)
        sol = run_seepage_analysis(seep, tol=1e-8)
    capture(f"docs/seep/images/{img_prefix}_solution.png", plot_seep_solution,
            seep, sol, base_mat=base_mat, levels=20, fill_contours=True, phreatic=phreatic)
    print(f"  {img_prefix}: q = {sol['flowrate']:.4f}")


def fem_case(xlsx, img_prefix):
    sd = load_slope_data(xlsx)
    capture(f"docs/fem/images/{img_prefix}_inputs.png", plot_inputs, sd, mode="fem")


if __name__ == "__main__":
    lem_case("docs/lem/files/xslope_acads_simple.xlsx", "acads_simple", "circular")
    lem_case("docs/lem/files/xslope_arai_tagyo.xlsx", "arai_tagyo", "circular")
    lem_case("docs/lem/files/xslope_acads_weak_layer.xlsx", "acads_weak_layer", "non_circular")
    seep_case("docs/seep/files/xslope_confined_radial.xlsx", "confined_radial", 1.0, phreatic=False)
    seep_case("docs/seep/files/xslope_sheetpile_s50.xlsx", "sheetpile_s50", 1.0, phreatic=False)
    fem_case("docs/fem/files/xslope_griffiths6_full.xlsx", "griffiths6_full")
    fem_case("docs/fem/files/xslope_griffiths6_dry.xlsx", "griffiths6_dry")
    print("\nbenchmark figures generated")
