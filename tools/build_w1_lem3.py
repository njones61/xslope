"""Build W-1's broken start file: LEM-3's layered slope with three planted faults.

The W-1 diagnosis session hands the assistant a model that is wrong and says only
that the factor of safety looks wrong. For that session to mean anything the
breakage has to be *exactly* known, so the file is built here rather than edited
by hand: LEM-3's own completed model
(``docs/lem/files/xslope_simple_mult_layers.xlsx``) is loaded, three faults are
written into it, and the result is saved through
:func:`xslope.fileio.save_slope_data_to_xlsx` — the writer Studio's Save As uses,
so the destination is rebuilt from a fresh blank template and no cell is edited in
place.

The three faults, each a mistake a person actually makes:

1. **The materials are in the wrong order.** The foundation is entered first and
   the embankment second, while the profile lines still reference material 1 for
   the ground surface and material 2 for the layer beneath it. The row order is
   what fixes the Mat IDs, so the section comes out upside down: the strong
   foundation clay on top, the fill underneath.
2. **The embankment unit weight is 13 pcf, not 130** — a dropped zero.
3. **The maximum depth is -100, not -10** — the rigid rock put ten times too
   deep, so the search can reach mechanisms the site does not have.

Nothing else differs from LEM-3. The builder proves that by reloading what it
wrote, undoing the three faults in memory, and re-solving: the result has to be
LEM-3's published Spencer answer, 1.244.

Run:  PYTHONPATH=. python3 tools/build_w1_lem3.py
"""

from __future__ import annotations

import contextlib
import copy
import io
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

import matplotlib                                                   # noqa: E402
matplotlib.use("Agg")

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx  # noqa: E402
from xslope.search import run_lem_analysis                          # noqa: E402

#: LEM-3's completed model — the one the tutorial page hands the reader, and the
#: project every other W-1 session opens.
SOURCE = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_mult_layers.xlsx")
DEST = os.path.join(REPO_ROOT, "docs/tutorials/files/w1_diagnose_start.xlsx")

#: LEM-3's published Spencer answer, on the circle tangent to the contact.
PUBLISHED_FS = 1.244
FS_TOL = 0.002

#: Fault 2: the embankment's unit weight, as typed and as it should be.
BROKEN_GAMMA = 13.0
GOOD_GAMMA = 130.0

#: Fault 3: the maximum depth, as typed and as it should be.
BROKEN_MAX_DEPTH = -100.0
GOOD_MAX_DEPTH = -10.0


def _solve(slope_data):
    """Spencer with a search, muted — the FS and the circle it was found on."""
    sd = copy.deepcopy(slope_data)
    with contextlib.redirect_stdout(io.StringIO()):
        bundle = run_lem_analysis(sd, "spencer", analysis="auto_search",
                                  announce=False)
    results = bundle.get("results")
    if not isinstance(results, dict):
        raise SystemExit("Spencer found no solvable surface: %s"
                         % (bundle.get("failure") or "no solution"))
    search = bundle.get("search") or {}
    # The search ranks its trials in ``fs_cache``; the first entry is the
    # critical circle it reports (the same one Studio's result header names).
    circle = (search.get("fs_cache") or [None])[0]
    return results["FS"], circle


def _load(path):
    with contextlib.redirect_stdout(io.StringIO()):
        return load_slope_data(path)


def _break(slope_data):
    """LEM-3's model with the three faults written into it."""
    sd = copy.deepcopy(slope_data)

    names = [m["name"] for m in sd["materials"]]
    if names != ["embankment", "foundation"]:
        raise SystemExit("expected LEM-3's two materials in fill-then-foundation "
                         "order, found %r" % (names,))

    # 1. The rows swap; the profile lines' Mat IDs do not, so line 1 (the ground
    #    surface) now names the foundation and line 2 (the contact) the fill.
    sd["materials"] = [sd["materials"][1], sd["materials"][0]]

    # 2. The embankment — now the second row — loses a zero.
    embankment = next(m for m in sd["materials"] if m["name"] == "embankment")
    embankment["gamma"] = BROKEN_GAMMA

    # 3. The rigid base goes ten times too deep.
    sd["max_depth"] = BROKEN_MAX_DEPTH
    return sd


def _repair(slope_data):
    """The three faults undone — nothing else touched."""
    sd = copy.deepcopy(slope_data)
    sd["materials"] = [sd["materials"][1], sd["materials"][0]]
    embankment = next(m for m in sd["materials"] if m["name"] == "embankment")
    embankment["gamma"] = GOOD_GAMMA
    sd["max_depth"] = GOOD_MAX_DEPTH
    return sd


def build():
    source = _load(SOURCE)
    broken = _break(source)
    save_slope_data_to_xlsx(broken, DEST)

    # Everything below is measured on what was WRITTEN, not on what was held in
    # memory: the file is what the session opens.
    back = _load(DEST)

    names = [m["name"] for m in back["materials"]]
    if names != ["foundation", "embankment"]:
        raise SystemExit("the material order did not survive the round trip: %r"
                         % (names,))
    gamma = back["materials"][1]["gamma"]
    if gamma != BROKEN_GAMMA:
        raise SystemExit("the embankment unit weight did not survive the round "
                         "trip: %r" % (gamma,))
    if back["max_depth"] != BROKEN_MAX_DEPTH:
        raise SystemExit("the maximum depth did not survive the round trip: %r"
                         % (back["max_depth"],))
    if [line["mat_id"] for line in back["profile_lines"]] != \
            [line["mat_id"] for line in source["profile_lines"]]:
        raise SystemExit("the profile lines' material references changed; the "
                         "fault is the material ORDER, with the lines left alone")

    fs_broken, circle_broken = _solve(back)
    fs_fixed, circle_fixed = _solve(_repair(back))

    print("wrote %s" % os.path.relpath(DEST, REPO_ROOT))
    print("  broken   Spencer FS = %.4f   %s" % (fs_broken, _circle(circle_broken)))
    print("  repaired Spencer FS = %.4f   %s" % (fs_fixed, _circle(circle_fixed)))

    if abs(fs_fixed - PUBLISHED_FS) > FS_TOL:
        raise SystemExit("undoing the three faults gives %.4f, not LEM-3's "
                         "published %.3f — the file carries a fourth difference"
                         % (fs_fixed, PUBLISHED_FS))
    if abs(fs_broken - fs_fixed) < 0.05:
        raise SystemExit("the broken model answers %.4f, within 0.05 of the "
                         "repaired %.4f — the faults have to be visible in the "
                         "factor of safety for the session to start from one"
                         % (fs_broken, fs_fixed))
    return fs_broken, fs_fixed


def _circle(circle):
    """One trial circle as ``Xo / Yo / R / Depth`` — the radius from the tangent
    elevation, which is what the cache carries."""
    if not circle:
        return ""
    try:
        xo, yo, depth = (float(circle["Xo"]), float(circle["Yo"]),
                         float(circle["Depth"]))
    except (KeyError, TypeError, ValueError):
        return ""
    return "(Xo %.2f, Yo %.2f, R %.2f, Depth %.2f)" % (xo, yo, yo - depth, depth)


if __name__ == "__main__":
    build()
