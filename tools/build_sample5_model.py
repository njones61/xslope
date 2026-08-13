"""Build LEM Sample Problem 5's model: three layers over a piezometric line.

``docs/lem/files/xslope_method_slices_problem.xlsx`` — the section of
[samples §5](../docs/lem/samples.md) and the model
[Tutorial LEM-4](../docs/tutorials/lem04_water_in_the_slope.md) builds.

The section is a 44 ft slope in three soils on a rigid base at elevation 0, with a
piezometric line falling from elevation 80 behind the crest to the toe and running
along the ground from there out. That geometry, that water line and the specified
circle at (195, 150) are the published exercise's, and this script pins all three:
they are declared here as constants and checked against the workbook before
anything is written, so the file cannot drift from what the two pages describe.

What this script sets is the **strength and unit-weight table**. Both pages turn on
water, and a water page needs a section whose critical mechanism is one that water
acts on. The materials here put a firm embankment (soil 1) and a firm middle layer
(soil 2) over a soft, nearly frictionless foundation clay (soil 3, c = 900 psf,
phi = 12 deg): the section's weakest mechanism is a deep foundation circle cutting
into soil 3, the piezometric line runs above most of its base, and the water is
therefore worth something on the surface the search actually finds. Every material
carries **both** unit weights — a gamma above the water table and a gamma_sat below
it — so the weight split the pages describe is live in the file as it ships rather
than an edit a reader has to make to see it.

Written through :func:`xslope.fileio.save_slope_data_to_xlsx`, the same writer
Studio's **Save As** uses: the destination is rebuilt from a fresh copy of the
current blank template and every input category is written from ``slope_data``, so
no cell is edited in place and no formula is overwritten. The round trip is
asserted here rather than assumed — the file is read back and compared against the
model that was written.

The companion file ``xslope_method_slices_problem2.xlsx`` is the same section under
a three-seed search; it is built by this script too, from the same materials, so
the two cannot disagree about the soil.

Run:  PYTHONPATH=. python3 tools/build_sample5_model.py
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

FILES = os.path.join(REPO_ROOT, "docs", "lem", "files")
SINGLE = os.path.join(FILES, "xslope_method_slices_problem.xlsx")
MULTI = os.path.join(FILES, "xslope_method_slices_problem2.xlsx")

#: The soil table, top down. Every strength is effective-stress (a c' and a phi')
#: read against the piezometric line, and every material carries both unit
#: weights: ``gamma`` above the line, ``gamma_sat`` below it. Both are TOTAL unit
#: weights — a saturated soil is a few pcf heavier than the same soil moist, and
#: neither column is ever a buoyant weight.
MATERIALS = [
    dict(name="soil 1", gamma=125.0, gamma_sat=130.0, c=400.0, phi=30.0),
    dict(name="soil 2", gamma=122.0, gamma_sat=127.0, c=600.0, phi=28.0),
    dict(name="soil 3", gamma=115.0, gamma_sat=118.0, c=900.0, phi=12.0),
]

#: What this script must not change, and checks it has not: the section, the water
#: line, the rigid base and the declared unit system.
MAX_DEPTH = 0.0
UNIT_SYSTEM = "imperial"
PROFILE_LINES = [
    (0, [(0.0, 84.0), (150.0, 84.0), (174.7, 64.0)]),
    (1, [(0.0, 64.0), (174.7, 64.0), (204.3, 40.0)]),
    (2, [(0.0, 40.0), (320.0, 40.0)]),
]
PIEZO_LINE = [(0.0, 80.0), (75.0, 79.0), (112.0, 76.0), (140.0, 70.0),
              (170.0, 61.0), (189.5, 52.0), (204.3, 40.0), (320.0, 40.0)]

#: The starting circles each file carries. The single-seed file states the
#: exercise's own circle; the companion states three seeds spread across the
#: section, which is what a global search is given.
CIRCLES = {
    SINGLE: [dict(Xo=195.0, Yo=150.0, Depth=18.1, R=131.9)],
    MULTI: [dict(Xo=175.0, Yo=100.0, Depth=0.0, R=100.0),
            dict(Xo=175.0, Yo=120.0, Depth=40.0, R=80.0),
            dict(Xo=175.0, Yo=120.0, Depth=64.0, R=56.0)],
}


def _check(path, sd):
    """Refuse to write a file whose geometry is not the one declared above."""
    def bad(what, got, want):
        raise SystemExit("%s: %s is %r, not the %r this script pins"
                         % (os.path.basename(path), what, got, want))

    if sd["max_depth"] != MAX_DEPTH:
        bad("the maximum depth", sd["max_depth"], MAX_DEPTH)
    if sd["unit_system"] != UNIT_SYSTEM:
        bad("the unit system", sd["unit_system"], UNIT_SYSTEM)
    got = [(ln["mat_id"], [tuple(map(float, p)) for p in ln["coords"]])
           for ln in sd["profile_lines"]]
    if got != PROFILE_LINES:
        bad("the section", got, PROFILE_LINES)
    if [tuple(map(float, p)) for p in sd["piezo_line"]] != PIEZO_LINE:
        bad("the piezometric line", sd["piezo_line"], PIEZO_LINE)
    want = CIRCLES[path]
    got = [{k: c[k] for k in ("Xo", "Yo", "Depth", "R")} for c in sd["circles"]]
    if got != want:
        bad("the starting circles", got, want)
    if len(sd["materials"]) != len(MATERIALS):
        bad("the material count", len(sd["materials"]), len(MATERIALS))


def build(path):
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(path)
    _check(path, sd)

    for mat, spec in zip(sd["materials"], MATERIALS):
        mat.update(spec, option="mc", u="piezo", ru=0.0)

    save_slope_data_to_xlsx(sd, path)

    # Round trip: the file has to read back as the model that was written.
    with contextlib.redirect_stdout(io.StringIO()):
        back = load_slope_data(path)
    _check(path, back)
    for got, spec in zip(back["materials"], MATERIALS):
        for key, want in list(spec.items()) + [("option", "mc"), ("u", "piezo")]:
            if got[key] != want:
                raise SystemExit("%s: round trip changed %s of %s: wrote %r, read %r"
                                 % (os.path.basename(path), key, spec["name"],
                                    want, got[key]))
    for key in ("profile_lines", "circles", "piezo_line", "max_depth",
                "unit_system", "water_loads", "dloads", "non_circ"):
        if repr(back[key]) != repr(sd[key]):
            raise SystemExit("%s: round trip changed %s:\n  wrote %r\n  read  %r"
                             % (os.path.basename(path), key, sd[key], back[key]))

    print("-> %s" % os.path.relpath(path, REPO_ROOT))
    for spec in MATERIALS:
        print("   %-7s g=%g gsat=%g c=%g phi=%g u=piezo"
              % (spec["name"], spec["gamma"], spec["gamma_sat"],
                 spec["c"], spec["phi"]))
    return path


if __name__ == "__main__":
    for target in (SINGLE, MULTI):
        build(target)
