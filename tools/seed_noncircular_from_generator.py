"""Seed Sample 7's failure surface from the weak-zone generator.

``xslope_noncircular.xlsx`` carried a hand-drawn four-point surface whose track
ran mid-seam at y = -5, entering 10 ft beyond the toe. It is a legal surface, but
it is a worse *seed*: a non-circular search starts from the surface the file
carries, and from that one Spencer settles at 1.7388, while the same search
seeded from the surface :func:`xslope.generators.generate_noncircular_surface`
builds on the same model reaches 1.6563 — the lower minimum, and the one the page
should report. The difference is where the track sits (10% of the seam's
thickness above its base, y = -5.8, rather than in its middle) and where the ends
land (at the toe and on the crest, rather than out on the flat ground).

So the file's surface becomes the generator's surface: this writes the
generator's exact output — four points, each with its explicit Y and Movement —
into the ``non-circ`` table. Sample 7 then ships the seed the search wants, and
the surface Tutorial LEM-5 teaches is the same one the generator builds, so the
page's Excel, assistant and Studio paths all arrive at one model.

Written through :func:`xslope.fileio.save_slope_data_to_xlsx`, the same writer
Studio's **Save As** uses: the destination is rebuilt from a fresh copy of the
current blank template and every input category is written from ``slope_data``,
so no cell is edited in place and no formula is overwritten. The round trip is
asserted rather than assumed, in both directions that matter — the table reads
back point for point, and the generator run on the *written* file reproduces it.

Run:  PYTHONPATH=. python3 tools/seed_noncircular_from_generator.py
"""

from __future__ import annotations

import contextlib
import io
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx   # noqa: E402
from xslope.generators import generate_noncircular_surface           # noqa: E402
from xslope.slice import generate_slices                             # noqa: E402
from xslope.solve import solve_selected                              # noqa: E402

TARGET = os.path.join(REPO_ROOT, "docs/lem/files/xslope_noncircular.xlsx")

#: The zone the generator has to seed on, and the point count it has to build.
#: A surface seeded on another zone, or one carrying stations the flat seam does
#: not need, is not the surface Sample 7 and Tutorial LEM-5 are written about.
ZONE = "Soft Clay"
POINTS = 4

#: Spencer's factor of safety on the written surface, solved as it stands (no
#: search). The page quotes it as what a generated shape is worth held still.
SPENCER_FS = 2.0216


def _spencer(sd, non_circ):
    ok, res = generate_slices(sd, non_circ=non_circ, num_slices=40)
    if not ok:
        raise SystemExit("generate_slices failed: %s" % (res,))
    slice_df, _ = res
    result = solve_selected("spencer", slice_df)
    if isinstance(result, str):
        raise SystemExit("spencer failed: %s" % result)
    return result["FS"]


def _points(surface):
    return [(round(p["X"], 10), round(p["Y"], 10), p["Movement"]) for p in surface]


def seed():
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(TARGET)
        out = generate_noncircular_surface(sd, report=True)

    if not out["surface"]:
        raise SystemExit("the weak-zone generator built nothing: %s" % out["reason"])
    if not out["confident"]:
        raise SystemExit("the generator found no clearly weakest zone; it must not "
                         "be guessing on this model")
    if out["zone"].name != ZONE:
        raise SystemExit("the generator seeded on %r, not %r" % (out["zone"].name, ZONE))
    if len(out["surface"]) != POINTS:
        raise SystemExit("the generator built %d points, expected %d"
                         % (len(out["surface"]), POINTS))
    for p in out["surface"]:
        if p.get("Y") is None or p.get("Movement") in (None, ""):
            raise SystemExit("a generated point is missing its Y or Movement: %r" % (p,))

    sd["non_circ"] = [dict(p) for p in out["surface"]]

    with contextlib.redirect_stdout(io.StringIO()):
        fs = _spencer(sd, sd["non_circ"])
    if abs(fs - SPENCER_FS) > 5e-4:
        raise SystemExit("the generated surface solves at Spencer %.4f, not the "
                         "%.4f this file is locked at" % (fs, SPENCER_FS))

    save_slope_data_to_xlsx(sd, TARGET)

    # Round trip: the file has to read back as the model that was written, and the
    # generator run on the written file has to reproduce the table it now holds.
    with contextlib.redirect_stdout(io.StringIO()):
        back = load_slope_data(TARGET)
        again = generate_noncircular_surface(back, report=True)
    if _points(back["non_circ"]) != _points(sd["non_circ"]):
        raise SystemExit("round trip changed the surface:\n  wrote %r\n  read  %r"
                         % (sd["non_circ"], back["non_circ"]))
    if _points(again["surface"]) != _points(back["non_circ"]):
        raise SystemExit("the generator no longer reproduces the file's surface:\n"
                         "  file      %r\n  generated %r"
                         % (back["non_circ"], again["surface"]))
    for key in ("profile_lines", "max_depth", "unit_system", "k_seismic",
                "tcrack_depth"):
        if repr(back[key]) != repr(sd[key]):
            raise SystemExit("round trip changed %s:\n  wrote %r\n  read  %r"
                             % (key, sd[key], back[key]))
    if [(m["name"], m["gamma"], m["c"], m["phi"], m["u"]) for m in back["materials"]] != \
       [(m["name"], m["gamma"], m["c"], m["phi"], m["u"]) for m in sd["materials"]]:
        raise SystemExit("round trip changed the material table")

    print("-> %s" % os.path.relpath(TARGET, REPO_ROOT))
    print("   %s" % out["summary"])
    for p in back["non_circ"]:
        print("   X = %-20.10g Y = %-6g %s" % (p["X"], p["Y"], p["Movement"]))
    print("   Spencer on it, held still: %.4f" % fs)
    return TARGET


if __name__ == "__main__":
    seed()
