"""Correct the pore-pressure options in Sample 7's model, ``xslope_noncircular.xlsx``.

The shipped file carried ``u = none`` on Dense Sand (material 4). That is a
transcription defect, not a modeling choice: Dense Sand is an effective-stress
soil (c' = 0, phi' = 37) lying entirely below the piezometric line at elevation
-2, so its strength must be read against the pore pressure the line produces —
``u = piezo``. Sand Fill (material 1) is the same kind of soil above the line,
where the piezometric reading is zero; it is set to ``piezo`` as well so that
every drained sand reads the same source and the file states the physics rather
than an accident of geometry. Soft Clay stays ``none``: its strength is
undrained (Su = 200 psf, phi = 0), a total-stress strength that no pore
pressure enters.

Neither change moves a number. Every base of the taught surface, of every
searched surface, and of every trial the sample's tags lock lies in soils whose
``u`` reading is unchanged — verified by running the full battery (single
surface all methods, non-circular search all methods, movement variants, track
elevations, the circular search on generated circles, and the weak-zone
generator) on the file before and after: identical to four decimals throughout.
This script asserts the taught surface's Spencer FS as its own guard.

Written through :func:`xslope.fileio.save_slope_data_to_xlsx`, the same writer
Studio's **Save As** uses: the destination is rebuilt from a fresh copy of the
current blank template and every input category is written from ``slope_data``,
so no cell is edited in place and no formula is overwritten. The round trip is
asserted rather than assumed.

Run:  PYTHONPATH=. python3 tools/fix_noncircular_u.py
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
from xslope.slice import generate_slices                             # noqa: E402
from xslope.solve import solve_selected                              # noqa: E402

TARGET = os.path.join(REPO_ROOT, "docs/lem/files/xslope_noncircular.xlsx")

#: The corrected u option per material, by name. The drained sands read the
#: piezometric line; the undrained clay reads nothing.
U_OPTIONS = {"Sand Fill": "piezo", "Sand": "piezo",
             "Soft Clay": "none", "Dense Sand": "piezo"}

#: The taught surface's Spencer factor of safety, which the u correction must
#: not move (no base of that surface lies in a soil whose reading changed).
SPENCER_FS = 1.7886


def _spencer(sd):
    ok, res = generate_slices(sd, non_circ=sd["non_circ"], num_slices=40)
    if not ok:
        raise SystemExit("generate_slices failed: %s" % (res,))
    slice_df, _ = res
    result = solve_selected("spencer", slice_df)
    if isinstance(result, str):
        raise SystemExit("spencer failed: %s" % result)
    return result["FS"]


def fix():
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(TARGET)

    for mat in sd["materials"]:
        want = U_OPTIONS.get(mat["name"])
        if want is None:
            raise SystemExit("unexpected material %r in %s"
                             % (mat["name"], os.path.basename(TARGET)))
        mat["u"] = want

    with contextlib.redirect_stdout(io.StringIO()):
        fs = _spencer(sd)
    if abs(fs - SPENCER_FS) > 5e-4:
        raise SystemExit("the u correction moved the taught surface: Spencer "
                         "%.4f, expected %.4f" % (fs, SPENCER_FS))

    save_slope_data_to_xlsx(sd, TARGET)

    # Round trip: the file has to read back as the model that was written.
    with contextlib.redirect_stdout(io.StringIO()):
        back = load_slope_data(TARGET)
    for mat in back["materials"]:
        if mat["u"] != U_OPTIONS[mat["name"]]:
            raise SystemExit("round trip lost a u option: %s reads %r"
                             % (mat["name"], mat["u"]))
    for key in ("profile_lines", "non_circ", "max_depth", "unit_system",
                "k_seismic", "tcrack_depth"):
        if repr(back[key]) != repr(sd[key]):
            raise SystemExit("round trip changed %s:\n  wrote %r\n  read  %r"
                             % (key, sd[key], back[key]))
    if [(m["name"], m["gamma"], m["c"], m["phi"], m["u"])
            for m in back["materials"]] != \
       [(m["name"], m["gamma"], m["c"], m["phi"], m["u"])
            for m in sd["materials"]]:
        raise SystemExit("round trip changed the material table")

    print("-> %s" % os.path.relpath(TARGET, REPO_ROOT))
    for mat in back["materials"]:
        print("   %-10s u = %s" % (mat["name"], mat["u"]))
    print("   Spencer on the taught surface: %.4f (unchanged)" % fs)
    return TARGET


if __name__ == "__main__":
    fix()
