"""Solve the two SIGMA/W wall models and write the companions the docs read.

The pair on docs/verification/geostudio.md ("Slope stabilization with a wall")
is the corpus's only structural-beam-in-a-slope benchmark, and its docs entry
reports the wall's internal forces, not just a factor of safety. Those come from
the solved field, so the field is committed beside the model:

    {stem}_fem_nodes.csv / _fem_elements.csv     the converged field
    {stem}_fem_meta.json                         the run record
    {stem}_fem_failure_*.csv / .json             the at-failure mechanism
    {stem}_fem_piles.csv (+ _fem_failure_piles)  the per-beam-element forces

so the figures redraw and the page's numbers can be re-read without a solve, and
``xslope.report.solutions_from_sidecars`` finds a complete finite-element bundle
beside each model.

The mesh is NOT rebuilt here. Both models carry ``u = 'seep'``, so their stability
run has to sit on the very mesh the seepage field was written on -- the mesh
companion the corpus builder emitted (benchmarks/geostudio/build_gs2_wall.py) --
and that is the mesh the suite's own tag path picks up too.

Everything else comes from the page's ``fem_ssrm`` tag: the bracket, the
tolerance, the iteration cap and the factor of safety the run must land on. A run
that misses its lock is reported and nothing is written, so a companion file can
never disagree with the number printed beside it.

Solves run on the pure-NumPy reference kernel, the path the locks are defined on.
Deterministic: same mesh, same bracket, same trials.

    PYTHONPATH=. python3 tools/make_gs2_wall_sidecars.py                 # both
    PYTHONPATH=. python3 tools/make_gs2_wall_sidecars.py gs2_wall        # one
    PYTHONPATH=. python3 tools/make_gs2_wall_sidecars.py --bracket 1.15 1.65
"""
import argparse
import contextlib
import io
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import run_tests as RT                                                  # noqa: E402
import xslope.fem as _fem                                               # noqa: E402
from xslope.fem import (build_fem_data, export_fem_solution,            # noqa: E402
                        solve_ssrm, ssrm_run_record)
from xslope.fileio import load_slope_data                               # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FILES = os.path.join(REPO_ROOT, "docs", "verification", "files", "geostudio")
PAGE = os.path.join(REPO_ROOT, "docs", "verification", "geostudio.md")

MODELS = ["gs2_wall_none", "gs2_wall"]


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _tag(stem):
    """The one ``fem_ssrm`` tag on the page that names this model."""
    want = f"{stem}.xlsx"
    hits = [t for t in RT.parse_test_tags(PAGE)
            if t.get("type") == "fem_ssrm"
            and os.path.basename(str(t.get("file", ""))) == want]
    if len(hits) != 1:
        raise RuntimeError(f"geostudio.md: {len(hits)} fem_ssrm tags name "
                           f"{want}, expected exactly one")
    return hits[0]


def build(stem, bracket=None, tolerance=None, max_iter=None):
    path = os.path.join(FILES, f"{stem}.xlsx")
    tag = None
    if bracket is None or tolerance is None:
        tag = _tag(stem)
    f_min, f_max = bracket if bracket else (float(tag["f_min"]),
                                            float(tag["f_max"]))
    tol = tolerance if tolerance is not None else float(tag["tolerance"])
    iters = int(max_iter if max_iter is not None
                else (tag["max_iter"] if tag else 16000))

    slope_data = load_slope_data(path)
    mesh = slope_data.get("mesh")
    if mesh is None:
        raise RuntimeError(f"{stem}: no mesh companion — run "
                           f"benchmarks/geostudio/build_gs2_wall.py first")
    fem_data = _quiet(build_fem_data, slope_data, mesh)

    t0 = time.time()
    options = {"F_min": f_min, "F_max": f_max, "tolerance": tol}
    with RT._force_fast_kernel(_fem, False):
        result = _quiet(solve_ssrm, fem_data, debug_level=0,
                        capture_failure_state=True, max_iterations=iters,
                        **options)
    if not result.get("converged"):
        raise RuntimeError(f"{stem}: SSRM did not converge")
    fs = result["FS"]

    # The tag is the truth. A run that lands outside its lock writes nothing:
    # a companion that disagrees with the number printed beside it is exactly
    # the defect this script exists to prevent.
    if tag is not None:
        expected = float(tag["expected_fs"])
        if abs(fs - expected) > tol:
            raise RuntimeError(
                f"{stem}: solved FS = {fs:.6f}, tag expects {expected:.3f} "
                f"+/- {tol:g} — nothing written")

    field = result["last_solution"]
    meta = ssrm_run_record(result, fem_data, options)
    meta.update({"analysis_type": "ssrm", "FS": fs, "element_type": "tri6",
                 "max_iter": iters, "F": field.get("F")})
    if tag is not None:
        meta["expected_fs"] = float(tag["expected_fs"])
    _quiet(export_fem_solution, fem_data, field, os.path.join(FILES, stem),
           meta=meta, failure_solution=result.get("failure_solution"))

    print(f"{stem}: FS = {fs:.4f}"
          + (f" (tag {float(tag['expected_fs']):.3f})" if tag else "")
          + f", {len(mesh['nodes'])} nodes, {len(mesh['elements'])} elements, "
          f"{fem_data.get('n_pile_elements', 0)} beam elements, "
          f"{len(result['trials'])} trials, {time.time() - t0:.0f}s")
    return fs


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("models", nargs="*", default=None,
                    help="model stems to build (default: both)")
    ap.add_argument("--bracket", nargs=2, type=float, metavar=("F_MIN", "F_MAX"),
                    help="override the tag's bracket (use before the tag exists)")
    ap.add_argument("--tolerance", type=float,
                    help="override the tag's bisection tolerance")
    ap.add_argument("--max-iter", type=int, help="override the iteration cap")
    args = ap.parse_args(argv)

    for stem in (args.models or MODELS):
        build(stem, bracket=args.bracket, tolerance=args.tolerance,
              max_iter=args.max_iter)


if __name__ == "__main__":
    main(sys.argv[1:])
