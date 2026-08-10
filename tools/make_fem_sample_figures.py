# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Redraw the figures on the FEM sample pages from the runs committed beside the
models.

The Griffiths & Lane figures on the verification pages are redrawn by
``benchmarks/make_griffiths_figures.py`` and the sample-page inputs panels for
two benchmark models by ``benchmarks/make_benchmark_figures.py``. The pile
sample's figures on ``docs/fem/samples.md`` were covered by neither: they were
grabbed by hand from a Studio session and nothing could redraw them, so when
``tools/make_fem_docs_sidecars.py`` rebuilt that sample's companions from the
page's own ``fem_ssrm`` tag the pictures stayed behind. The mesh panel still
announced the 1,350 elements of the session's mesh where the tag builds 1,521,
and the results panels were titled at a trial factor of 1.38 that belongs to no
run the page reports.

Nothing is solved here. Each figure is drawn from the sidecars
``make_fem_docs_sidecars.py`` writes, so this script and the page describe one
run, and a rebuild of the companions is followed by a rebuild of the pictures.

The results panels are drawn the way the verification figures and the report
draw them: the last converged field is passed with the at-failure snapshot
beside it, and the panels trace the developed mechanism.

The pile summary both fields produce is printed at the end. The page quotes that
block, and quoting it from a field is the whole reason the numbers on the page
have to say which field they came from.

The reinforcement and non-circular samples' figures are still hand grabs with no
producer; only the pile sample is covered here.

    PYTHONPATH=. python3 tools/make_fem_sample_figures.py
"""
import contextlib
import io
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib                                                    # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                      # noqa: E402

from xslope.fem import print_pile_summary                            # noqa: E402
from xslope.fileio import load_slope_data                            # noqa: E402
from xslope.plot import plot_inputs                                  # noqa: E402
from xslope.plot_fem import plot_fem_data, plot_fem_results          # noqa: E402
from xslope.report import solutions_from_sidecars                    # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(ROOT, "docs", "fem", "files")
OUT = os.path.join(ROOT, "docs", "fem", "images")
DPI = 200

#: ``stem -> figure prefix`` for the samples this script draws.
CASES = {"xslope_piles_fem": "piles_fem"}

#: Figures on the pile sample's section that this script does NOT draw, each with
#: the reason. ``piles_fem_results_no_pile`` is the same slope with the piles
#: taken out — a different model, with no input file and no companions committed
#: for it, so there is nothing here to draw it from.
NOT_DRAWN = {
    "piles_fem_results_no_pile": "a model with the piles removed, which is not "
                                 "committed and ships no solved run",
}


def _save(prefix, kind, note=""):
    path = os.path.join(OUT, f"{prefix}_{kind}.png")
    plt.gcf().savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close("all")
    print(f"  wrote {os.path.relpath(path, ROOT)}{note}")


def draw(stem, prefix):
    path = os.path.join(FILES, f"{stem}.xlsx")
    with contextlib.redirect_stdout(io.StringIO()):
        slope_data = load_slope_data(path)
        bundle = solutions_from_sidecars(path, slope_data, None).get("fem")
    if not bundle:
        raise RuntimeError(f"{stem}: no solved run beside the model — run "
                           f"tools/make_fem_docs_sidecars.py first")
    fem_data, solution = bundle["fem_data"], bundle["solution"]
    mesh = slope_data.get("mesh") or {}

    print(f"{stem}: FS = {bundle.get('FS'):.4f}, "
          f"{len(mesh.get('nodes', []))} nodes, "
          f"{len(mesh.get('elements', []))} elements")

    # frame="content" — equal aspect with the axes BOX shrunk to the section,
    # which is what the report's own model figure uses. The default "fill"
    # framing keeps the box and expands the data limits to match it, and on a
    # section 110 ft wide and 30 ft tall that puts 30 ft of empty ground above
    # and below the slope. The mesh underlay is left on: it is what the page's
    # figure has always shown, and the mesh figure of its own is much further
    # down the page.
    with contextlib.redirect_stdout(io.StringIO()):
        plot_inputs(slope_data, mode="fem", frame="content")
    _save(prefix, "inputs")

    with contextlib.redirect_stdout(io.StringIO()):
        plot_fem_data(fem_data)
    _save(prefix, "mesh", f"  ({len(mesh.get('elements', []))} elements)")

    with contextlib.redirect_stdout(io.StringIO()):
        plot_fem_results(fem_data, solution, fs=bundle.get("FS"),
                         failure_solution=bundle.get("failure_solution"))
    _save(prefix, "results")

    for name, why in NOT_DRAWN.items():
        if name.startswith(prefix):
            print(f"  {name}.png not drawn: {why}")

    # What the page's summary block reads, from each field it could be read at.
    for label, field in (("last converged", solution),
                         ("at failure", bundle.get("failure_solution"))):
        if field is None:
            continue
        print(f"\n--- print_pile_summary, {label} ---")
        print_pile_summary(fem_data, field)


def main():
    os.makedirs(OUT, exist_ok=True)
    for stem, prefix in CASES.items():
        draw(stem, prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())
