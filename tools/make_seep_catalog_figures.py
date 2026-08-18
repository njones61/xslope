"""Regenerate the solved figures on docs/seep/samples.md that have producers.

Each figure re-runs its sample the same way the page's regression tag does
(width/120 tri3 mesh, the file's own declared units), so the drawn title and
the tag lock the same number.

    python tools/make_seep_catalog_figures.py

Covered here: sample #1's solution (clay_blanket_solution.png). The remaining
catalog images predate this script and do not yet have producers.
"""

import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(_HERE)
OUT = os.path.join(REPO, "docs", "seep", "images")


def clay_blanket_solution():
    from xslope.fileio import load_slope_data
    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)
    from xslope.plot_seep import plot_seep_solution
    from xslope.seep import build_seep_data, run_seepage_analysis

    sd = load_slope_data(os.path.join(REPO, "docs", "seep", "files",
                                      "xslope_clay_blanket.xlsx"))
    polygons = get_material_polygons(sd)
    xs = [x for poly in polygons for x, _ in poly["coords"]]
    target = (max(xs) - min(xs)) / 120           # the seep tag's own default
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            polygons, target, "tri3", size_regions=extract_size_regions(sd))
        seep_data = build_seep_data(mesh, sd)
        solution = run_seepage_analysis(seep_data, tol=1e-4)
    fig = plot_seep_solution(seep_data, solution, mesh=False)
    path = os.path.join(OUT, "clay_blanket_solution.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(path, REPO),
          "| q = %.3f" % solution["flowrate"])


if __name__ == "__main__":
    clay_blanket_solution()
