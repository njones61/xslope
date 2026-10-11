"""No standard plot draws Matplotlib's "Ignoring fixed y limits" warning.

The plots that frame a section with ``set_aspect('equal', adjustable='datalim')``
— the inputs plot, the solution plot, the mesh plot and Studio's editor previews —
also set the window that fit starts from. Set with ``set_ylim``/``set_xlim``, that
window is FIXED, and the data-lim aspect overrides it on every draw with
"Ignoring fixed y limits to fulfill fixed data aspect with adjustable data limits."
The line reached the assistant's tool results and the dock, and filled the console.
Set as bounds (``set_ybound``/``set_xbound``) it is the same starting window, left
adjustable, so every figure is drawn exactly as before and nothing is logged.

Renders each standard plot on shipped models at several figure shapes (Studio's
canvases are sized to their panes, so the shape varies) and fails on any record of
that message.

Run: python test/plot_aspect_limits_check.py
"""
import contextlib
import io
import logging
import os
import sys
import warnings

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

NEEDLE = "to fulfill fixed data aspect"
SHAPES = ((11, 6), (6, 6), (5, 7), (14, 5))
_hits = []


class _Catch(logging.Handler):
    def emit(self, record):
        if NEEDLE in record.getMessage():
            _hits.append(record.getMessage())


def _load(rel):
    from xslope.fileio import load_slope_data
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return load_slope_data(os.path.join(_REPO, rel))


def _drawn(label, make, out):
    """Build one figure with ``make(fig)``, draw it as the kernel saves it, and
    record a failure if the warning was logged at any point."""
    del _hits[:]
    fig = plt.figure(figsize=make.size)
    try:
        make(fig)
        fig.savefig(io.BytesIO(), dpi=80, bbox_inches="tight")
    finally:
        plt.close(fig)
    if _hits:
        out.append(f"{label} at {make.size}: {_hits[0]!r} ({len(_hits)} time(s))")


def _sized(fn, size):
    fn.size = size
    return fn


def check_lem_plots():
    from xslope import plot as P
    from xslope.search import run_lem_analysis
    out = []
    sd = _load("docs/lem/files/xslope_simple_mult_layers.xlsx")
    with contextlib.redirect_stdout(io.StringIO()):
        b = run_lem_analysis(sd, "spencer", analysis="auto_search",
                             surface="circular", announce=False)
    s = b["search"]
    for size in SHAPES:
        _drawn("plot_inputs", _sized(lambda f: P.plot_inputs(sd, fig=f), size), out)
        _drawn("plot_solution", _sized(lambda f: P.plot_solution(
            sd, b["slice_df"], b["failure_surface"], b["results"], fig=f), size), out)
        _drawn("plot_solution with slice numbers", _sized(lambda f: P.plot_solution(
            sd, b["slice_df"], b["failure_surface"], b["results"], fig=f,
            slice_numbers=True), size), out)
        _drawn("plot_circular_search_results", _sized(
            lambda f: P.plot_circular_search_results(
                sd, s["fs_cache"], s["search_path"], s.get("circle_cache"), fig=f),
            size), out)
    return out


def check_seep_and_fem_plots():
    from xslope import plot as P
    from xslope.mesh import build_mesh_from_polygons, get_material_polygons
    out = []
    seep = _load("docs/seep/files/xslope_earth_dam1.xlsx")
    fem = _load("docs/fem/files/xslope_griffiths1.xlsx")
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        xs = [x for x, _ in fem["ground_surface"].coords]
        mesh = build_mesh_from_polygons(get_material_polygons(fem),
                                        target_size=(max(xs) - min(xs)) / 30,
                                        element_type="tri6")
    for size in SHAPES:
        _drawn("plot_inputs (seepage)", _sized(
            lambda f: P.plot_inputs(seep, mode="seep", fig=f), size), out)
        _drawn("plot_inputs (finite element)", _sized(
            lambda f: P.plot_inputs(fem, mode="fem", fig=f), size), out)
        _drawn("plot_mesh", _sized(
            lambda f: P.plot_mesh(mesh, materials=fem["materials"], fig=f), size), out)
    return out


def check_editor_preview():
    """Studio's editor previews finish their axes the way plot_inputs does."""
    from studio.editors import _finish_preview_axes
    out = []
    sd = _load("docs/lem/files/xslope_noncircular.xlsx")

    def preview(fig):
        ax = fig.add_subplot(111)
        for pl in sd["profile_lines"]:
            xs, ys = zip(*pl["coords"])
            ax.plot(xs, ys)
        _finish_preview_axes(ax)

    for size in SHAPES + ((6, 3.5),):
        _drawn("an editor preview", _sized(preview, size), out)
    return out


CHECKS = [
    ("LEM inputs, solution and search plots", check_lem_plots),
    ("seepage and finite element plots", check_seep_and_fem_plots),
    ("Studio editor previews", check_editor_preview),
]


def _listen():
    log = logging.getLogger("matplotlib.axes._base")
    if not any(isinstance(h, _Catch) for h in log.handlers):
        log.addHandler(_Catch())


def run(verbose=False):
    """Every leg's failures as one list (empty when no plot logs the warning)."""
    _listen()
    failures = []
    for name, fn in CHECKS:
        fs = fn()
        if verbose:
            print(f"  {name:42s} {'ok' if not fs else 'FAIL (%d)' % len(fs)}")
        failures += [f"{name}: {f}" for f in fs]
    return failures


def main():
    print("plot aspect/limits checks:")
    failures = run(verbose=True)
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print("  - " + f)
        return 1
    print("\nNo standard plot draws the fixed-limits warning.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
