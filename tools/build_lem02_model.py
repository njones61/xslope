"""Build Tutorial LEM-2's completed model: the LEM-1 embankment plus a crest surcharge.

LEM-2 does not build a slope. It starts from the model
[LEM-1](../docs/tutorials/lem01_simple_embankment.md) leaves behind —
``docs/lem/files/xslope_simple_embankment.xlsx``, the sample-page model — and adds
the one thing the tutorial is about: a 750 psf surcharge across a 10 ft strip of
the crest, from x = 25 to x = 35. Everything else is LEM-1's, untouched.

The load is the one LEM Sample Problem 1 already catalogues as its variant (a);
that sample's own file (``xslope_simple_embankment_mods.xlsx``) carries three
further changes on top of it — a 3 ft water-filled tension crack and a 10 ft
submergence — so it is not the state a tutorial about surcharges can start the
reader from. Hence a file of its own, holding exactly one change from LEM-1.

Written through :func:`xslope.fileio.save_slope_data_to_xlsx`, the same writer
Studio's **Save As** uses: the destination is rebuilt from a fresh copy of the
current blank template and every input category is written from ``slope_data``, so
no cell is edited in place and no formula is overwritten. The round trip is
asserted here rather than assumed — the file is read back and compared against the
model that was written.

The tutorial's other states (the line load, the load direction on the face, the
seismic coefficient, the design sweep) are edits the page has the reader make to
this model. They are not files: LEM-1 does the same with its tension crack.

Run:  PYTHONPATH=. python3 tools/build_lem02_model.py
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

SOURCE = os.path.join(REPO_ROOT, "docs/lem/files/xslope_simple_embankment.xlsx")
DEST = os.path.join(REPO_ROOT, "docs/lem/files/xslope_crest_surcharge.xlsx")

#: The surcharge: a uniform 750 psf over the crest between x = 25 and x = 35.
#: Two points is the whole load — the intensity is constant along it, and the
#: strip's ends are where the two points sit.
SURCHARGE = [{"X": 25.0, "Y": 20.0, "Normal": 750.0},
             {"X": 35.0, "Y": 20.0, "Normal": 750.0}]

#: ``normal`` is the loader's default for a blank Direction cell, and on the level
#: crest it is identical to ``vertical`` — the loaded line is horizontal, so the
#: perpendicular IS straight down. Written explicitly all the same, because the
#: tutorial's face-load section turns on the two words meaning different things.
DIRECTION = "normal"


def build():
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(SOURCE)

    if sd.get("dloads"):
        raise SystemExit("%s already carries a distributed load; LEM-2's model is "
                         "LEM-1's plus exactly one." % os.path.basename(SOURCE))

    sd["dloads"] = [list(SURCHARGE)]
    sd["dload_dirs"] = [DIRECTION]

    save_slope_data_to_xlsx(sd, DEST)

    # Round trip: the file has to read back as the model that was written.
    with contextlib.redirect_stdout(io.StringIO()):
        back = load_slope_data(DEST)
    if back["dloads"] != [list(SURCHARGE)]:
        raise SystemExit("round trip lost the surcharge: %r" % (back["dloads"],))
    if back["dload_dirs"] != [DIRECTION]:
        raise SystemExit("round trip lost the direction: %r" % (back["dload_dirs"],))
    for key in ("materials", "profile_lines", "circles", "max_depth",
                "unit_system", "k_seismic", "tcrack_depth"):
        if repr(back[key]) != repr(sd[key]):
            raise SystemExit("round trip changed %s:\n  wrote %r\n  read  %r"
                             % (key, sd[key], back[key]))

    print("-> %s" % os.path.relpath(DEST, REPO_ROOT))
    print("   %s + %g psf over x %g..%g (Direction %s)"
          % (os.path.basename(SOURCE), SURCHARGE[0]["Normal"],
             SURCHARGE[0]["X"], SURCHARGE[-1]["X"], DIRECTION))
    return DEST


if __name__ == "__main__":
    build()
