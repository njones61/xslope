"""The non-circular search may not bend the surface back down into the mass.

With every point free, the search used to notch the surface until Spencer returned
a root far below the critical circle. Started from eight points on VP17's critical
circle (F = 1.341), it settled on a surface bending back 84 degrees at one vertex,
where Spencer read 0.456 and Janbu, on the same slices, 2.04. Spencer's root there
passes every test Spencer applies to a root. The search now refuses any trial
surface that bends back more than ``MAX_KINK_BACK`` (10 degrees) at a vertex,
before it is sliced or solved, and the same start settles at 1.337 with Janbu
within 2.5%.

Four legs:

* the guarded search from that start lands at its circle, with Janbu agreeing and
  no bend past the limit;
* the same search with the guard off (``max_kink_back=None``) still takes the false
  minimum, which is what makes the first leg a test of the guard rather than of
  the start. If this leg ever fails, the search or Spencer changed and the
  guard's measurements (search.MAX_KINK_BACK) need taking again;
* a surface entered with a 36.5-degree bend at one vertex is searched rather than
  refused, and the search never takes that vertex past its entered bend, or any
  other vertex past 10 degrees;
* ``kinks_back`` reads the same bend whichever way round a surface is entered,
  reads 0 for a surface bending only the other way, and reads 0 beside a segment
  of no length.

The guarded legs take a few seconds; the unguarded one takes about a minute.
"""
import contextlib
import io
import os
import sys
import warnings

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)

from xslope import solve  # noqa: E402
from xslope.fileio import load_slope_data  # noqa: E402
from xslope.search import MAX_KINK_BACK, kinks_back, noncircular_search  # noqa: E402

VP17 = os.path.join(ROOT, "docs", "verification", "files", "rocscience", "vp017.xlsx")
NUM_SLICES = 50
CIRCLE_FS = 1.3408           # Spencer's critical circle on VP17, 50 slices

#: Eight points on that circle, evenly spaced in x, the ends on the ground.
START = [(4.889, 5.0), (6.699, 4.454), (8.51, 4.263), (10.32, 4.41),
         (12.131, 4.909), (13.941, 5.825), (15.752, 7.323), (17.562, 10.0)]


def _search(points, max_kink_back=MAX_KINK_BACK):
    sd = load_slope_data(VP17)
    sd["non_circ"] = [{"X": x, "Y": y, "Movement": "Free"} for x, y in points]
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        cache, _conv, _path = noncircular_search(sd, "spencer", num_slices=NUM_SLICES,
                                                 max_kink_back=max_kink_back)
    return cache[0] if cache else None


def _janbu(best):
    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = solve.janbu(best["slices"].copy())
    return res["FS"] if ok else None


def check_guarded_search_lands_at_its_circle():
    best = _search(START)
    if best is None or not np.isfinite(best["FS"]):
        return ["the guarded search from VP17's circle found no surface"]
    fails = []
    fs = best["FS"]
    if not 0.98 * CIRCLE_FS <= fs <= 1.005 * CIRCLE_FS:
        fails.append(f"guarded search F = {fs:.4f}, expected at or just under the "
                     f"circle's {CIRCLE_FS}")
    fj = _janbu(best)
    if fj is None or abs(fj - fs) > 0.05 * fs:
        fails.append(f"guarded search: Spencer {fs:.4f} against Janbu {fj} on the "
                     f"same slices")
    worst = float(kinks_back(best["points"]).max())
    if worst > MAX_KINK_BACK + 1e-9:
        fails.append(f"guarded search bends back {worst:.1f} degrees")
    return fails


def check_unguarded_search_still_takes_the_false_minimum():
    best = _search(START, max_kink_back=None)
    if best is None or not best["FS"] < 0.7 * CIRCLE_FS:
        got = None if best is None else round(best["FS"], 4)
        return [f"with the guard off the search returned F = {got}, not the false "
                f"minimum (0.456): the guard is no longer what this check tests; "
                f"re-measure search.MAX_KINK_BACK"]
    return []


def check_an_entered_bend_is_kept_as_its_own_limit():
    pts = list(START)
    pts[3] = (10.32, 5.2)                       # bends back 36.5 degrees at vertex 3
    entered = kinks_back(pts)
    if not entered[2] > 30 or entered[[0, 1, 3, 4, 5]].max() > 1e-9:
        return [f"the test surface's bends are {np.round(entered, 1).tolist()}"]
    best = _search(pts)
    if best is None or not np.isfinite(best["FS"]):
        return ["a surface entered with a 36.5-degree bend was refused"]
    got = kinks_back(best["points"])
    limit = np.maximum(MAX_KINK_BACK, entered)
    if np.any(got > limit + 1e-9):
        return [f"the search bent the surface back {np.round(got, 1).tolist()} "
                f"against limits {np.round(limit, 1).tolist()}"]
    return []


def check_the_bend_measure():
    fails = []
    notch = [(0, 0), (1, 1), (2, 1), (3, 3)]
    if not np.allclose(kinks_back(notch), [45.0, 0.0]):
        fails.append(f"a 45-degree bend reads {kinks_back(notch)}")
    if not np.allclose(kinks_back(notch[::-1]), kinks_back(notch)):
        fails.append("the bend reads differently with the points entered right to left")
    if np.any(kinks_back([(0, 0), (1, -1), (2, -1), (3, 1)]) != 0):
        fails.append("a surface bending only upward reads a bend back")
    if np.any(kinks_back([(0, 0), (1, 1), (1, 1), (2, 0)]) != 0):
        fails.append("a vertex beside a segment of no length reads a bend")
    return fails


def run():
    """Failures as a list, for run_tests.py."""
    fails = []
    for check in (check_the_bend_measure, check_guarded_search_lands_at_its_circle,
                  check_an_entered_bend_is_kept_as_its_own_limit,
                  check_unguarded_search_still_takes_the_false_minimum):
        try:
            fails += check()
        except Exception as exc:          # a crash is a failure, not a pass
            fails.append(f"{check.__name__} raised {exc!r}")
    return fails


if __name__ == "__main__":
    failures = run()
    if failures:
        print("FAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("Non-circular bend-back guard checks passed.")
