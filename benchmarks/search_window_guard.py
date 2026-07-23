"""Guard for circular_search's optional search-window constraints.

The adaptive circular search accepts four optional limits — ``center_box``,
``entry_range``, ``exit_range`` and ``tangent_depth`` (Slide2's "search limits"
semantics; see xslope.search.circular_search). This guard asserts the two
properties those limits must satisfy:

  1. INVARIANCE. With every limit left at its default (None) the search is
     byte-for-byte the historical unconstrained search — a change here that
     perturbs the default path would silently move every locked FS in the corpus.

  2. LOCAL-MINIMUM SELECTION. On the RS2-61 benched geometry (Cheng, Lansivaara &
     Wei 2007), which has competing minima by construction, the UNCONSTRAINED
     Spencer search settles on the GLOBAL minimum (~1.338, the mid-lower-face
     circle), while the SAME search CONFINED to the upper-face window (crest entry,
     first-bench exit, bench-elevation tangent — Fig. 5 of RS2 Verification Part 3
     Problem 61) settles on the distinct upper-face LOCAL minimum (~1.437). The two
     must differ by a clear margin, proving the window redirected the search rather
     than being ignored or clamping onto the same surface.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/search_window_guard.py
Exits non-zero on any failure.
"""
import io
import contextlib
import sys

from xslope.fileio import load_slope_data
from xslope.search import circular_search

RS2_61 = "docs/verification/files/rocscience/rs2_61a.xlsx"
# Three diverse existing circular_search benchmarks used for the invariance check.
INVARIANCE = [
    ("docs/verification/files/rocscience/rs2_60a.xlsx", "spencer", 40),
    ("docs/verification/files/rocscience/rs2_61a.xlsx", "spencer", 40),
    ("docs/verification/files/rocscience/vp002.xlsx", "bishop", 40),
]
# Upper-face local-minimum window (RS2-61 Case 3), read from Fig. 5:
#   entry on the crest bench (x 42.91-54.32), exit at the first bench (x 26-31.94),
#   tangent bottoming out near the bench elevation (~19.86).
CASE3_WINDOW = dict(entry_range=(42, 54), exit_range=(23, 32), tangent_depth=(16, 22))

GLOBAL_FS = 1.338
LOCAL_FS = 1.437
FS_TOL = 0.01


def _search(path, method, ns, **kw):
    d = load_slope_data(path)
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _conv, _sp, _cc = circular_search(
            d, method, num_slices=ns, seed="circles", **kw)
    if not fs_cache or fs_cache[0]["FS"] >= 9999:
        return None
    return fs_cache[0]["FS"]


def main():
    failures = []

    # (1) Invariance: passing the constraints explicitly as None must reproduce the
    # no-argument call to full float precision.
    for path, method, ns in INVARIANCE:
        base = _search(path, method, ns)
        with_none = _search(path, method, ns, center_box=None, entry_range=None,
                            exit_range=None, tangent_depth=None)
        name = path.split("/")[-1]
        if base is None or with_none is None or base != with_none:
            failures.append(f"invariance {name}: {base!r} != {with_none!r}")
        else:
            print(f"[invariance] {name:16s} {method:8s} FS={base:.10f}  (None == default)")

    # (2) Local-minimum selection on RS2-61.
    fs_global = _search(RS2_61, "spencer", 40)
    fs_local = _search(RS2_61, "spencer", 40, **CASE3_WINDOW)
    print(f"[select] RS2-61 unconstrained -> global FS={fs_global:.4f} (published Slide2 1.336)")
    print(f"[select] RS2-61 upper-face window -> local FS={fs_local:.4f} (published Slide2 1.443)")

    if fs_global is None or abs(fs_global - GLOBAL_FS) > FS_TOL:
        failures.append(f"global min: {fs_global} not ~{GLOBAL_FS}")
    if fs_local is None or abs(fs_local - LOCAL_FS) > FS_TOL:
        failures.append(f"local min: {fs_local} not ~{LOCAL_FS}")
    if fs_global is not None and fs_local is not None and fs_local - fs_global < 0.05:
        failures.append(f"window did not redirect the search: local {fs_local} "
                        f"vs global {fs_global} (< 0.05 apart)")

    if failures:
        print("\nFAILED:")
        for f in failures:
            print("  -", f)
        return 1
    print("\nOK: search-window constraints are invariant by default and select the "
          "RS2-61 upper-face local minimum when applied.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
