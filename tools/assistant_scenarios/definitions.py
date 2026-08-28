"""The scenarios — thirty conversations, and what each one has to get right.

Coverage is by FAMILY rather than by file: building (from a description and from a
drawing), editing (geometry on both geometry sources, loads, materials,
reinforcement, piles), running each engine (LEM, seepage, FEM), the parametric and
reliability modes, answering without running, diagnosing a model broken on
purpose, and the two behaviors that are neither — asking when a request is
ambiguous, and catching a number handed over in the wrong unit.

The eight W-1 sessions are here verbatim, prompts unchanged, so a scored run is
directly comparable with the recorded ones the tutorial was built from.
"""

from __future__ import annotations

import os
import re

from .core import Scenario, load_model, repo
from . import faults as F
from . import scorers as S

# --------------------------------------------------------------------------- #
# Models
# --------------------------------------------------------------------------- #
REINF = repo("docs/tutorials/files/xslope_reinforced_slope.xlsx")
SURCHARGE = repo("docs/lem/files/xslope_crest_surcharge.xlsx")
LAYERS = repo("docs/lem/files/xslope_eight_layers.xlsx")
GEOGRID = repo("docs/inputs/slope/xslope_reinf.xlsx")
PILES = repo("docs/tutorials/files/xslope_pile_wall.xlsx")
POLY = repo("docs/lem/files/xslope_sloping_bottom.xlsx")
DAM = repo("docs/inputs/slope/xslope_dam.xlsx")
SEEP_DAM = repo("docs/seep/files/xslope_earth_dam1.xlsx")
SSRM = repo("docs/tutorials/files/xslope_ssrm_embankment.xlsx")
DESIGN = repo("docs/lem/files/xslope_design.xlsx")
HASSIOTIS = repo("docs/lem/files/xslope_hassiotis.xlsx")
RELIABILITY = repo("docs/inputs/slope/xslope_reliability.xlsx")
RAPID = repo("docs/inputs/slope/xslope_rapid.xlsx")
ROCK = repo("docs/tutorials/files/xslope_rock_slope.xlsx")

LEM08_SKETCH = repo("docs/tutorials/images/lem08_problem_sketch.png")
LEM02_SKETCH = repo("docs/tutorials/images/lem02_problem_sketch.png")


# --------------------------------------------------------------------------- #
# Reading a model back
# --------------------------------------------------------------------------- #
def ground(sd):
    """The ground surface as a list of (x, y), whichever source built it."""
    gs = sd.get("ground_surface")
    if gs is None:
        return []
    try:
        return [(float(x), float(y)) for x, y in gs.coords]
    except Exception:
        return [(float(x), float(y)) for x, y in gs]


def face(sd):
    """``(height, run_per_rise, toe_x, crest_x)`` of the steepest descending run.

    Read off the ground surface rather than off a profile line, so it answers the
    same question for a polygon-native model and for one built from profile lines.
    """
    pts = ground(sd)
    if len(pts) < 2:
        return (0.0, 0.0, 0.0, 0.0)
    best = max(range(len(pts) - 1),
               key=lambda i: abs(pts[i + 1][1] - pts[i][1]))
    (x1, y1), (x2, y2) = pts[best], pts[best + 1]
    rise = abs(y2 - y1)
    run = abs(x2 - x1)
    lo, hi = (x1, x2) if y1 < y2 else (x2, x1)
    return (rise, (run / rise) if rise else 0.0, lo, hi)


def material(sd, index=0, name=None):
    mats = sd.get("materials") or []
    if name is not None:
        for mat in mats:
            if str(mat.get("name", "")).strip().lower() == name.lower():
                return mat
        return {}
    return mats[index] if index < len(mats) else {}


def near(value, want, rel=0.02, absolute=None):
    if value is None:
        return False
    tol = absolute if absolute is not None else rel * max(1.0, abs(want))
    return abs(float(value) - float(want)) <= tol


def reinforcement(sd):
    return sd.get("reinforcement_lines") or sd.get("reinforce_lines") or []


# --------------------------------------------------------------------------- #
# Predicates used by more than one scenario
# --------------------------------------------------------------------------- #
def units_are(system):
    def predicate(sd):
        got = str(sd.get("unit_system") or "").lower()
        return got == system, "unit_system = %r" % got
    return predicate


def materials_are(*specs):
    """Each spec is ``(gamma, c, phi)``; order-independent, tolerance 2%."""
    def predicate(sd):
        mats = sd.get("materials") or []
        if len(mats) != len(specs):
            return False, "%d material(s), wanted %d" % (len(mats), len(specs))
        unmatched = list(specs)
        for mat in mats:
            hit = next((s for s in unmatched
                        if near(mat.get("gamma"), s[0])
                        and near(mat.get("c"), s[1], absolute=max(1.0, 0.02 * s[1]))
                        and near(mat.get("phi"), s[2], absolute=0.6)), None)
            if hit is None:
                return False, "no spec matches %s (γ %s, c %s, φ %s)" % (
                    mat.get("name"), mat.get("gamma"), mat.get("c"), mat.get("phi"))
            unmatched.remove(hit)
        return True, "%d material(s) match" % len(specs)
    return predicate


def slope_face_is(height, run_per_rise, rel=0.05):
    def predicate(sd):
        h, ratio, _lo, _hi = face(sd)
        ok = near(h, height, rel=rel) and near(ratio, run_per_rise, rel=rel)
        return ok, "face %.1f high at %.2f:1 (wanted %g at %g:1)" % (
            h, ratio, height, run_per_rise)
    return predicate


def extents_beyond(multiple=2.0):
    """Flat ground at least ``multiple`` × the slope height past toe and crest."""
    def predicate(sd):
        pts = ground(sd)
        if not pts:
            return False, "no ground surface"
        height, _r, lo, hi = face(sd)
        left = lo - min(p[0] for p in pts)
        right = max(p[0] for p in pts) - hi
        want = multiple * height
        ok = left >= want - 1e-6 and right >= want - 1e-6
        return ok, "%.0f left / %.0f right of the face, wanted %.0f each" % (
            left, right, want)
    return predicate


def has_starting_circles(at_least=2):
    def predicate(sd):
        circles = sd.get("circles") or []
        return len(circles) >= at_least, "%d starting circle(s)" % len(circles)
    return predicate


def circle_tangent_to(depth, tol=0.5):
    def predicate(sd):
        depths = [c.get("Depth") for c in sd.get("circles") or []]
        hit = [d for d in depths if d is not None and abs(float(d) - depth) <= tol]
        return bool(hit), "circle depths %s (wanted one at %g)" % (
            [None if d is None else round(float(d), 2) for d in depths], depth)
    return predicate


def base_below_toe(depth, tol=0.3):
    """``max_depth`` sits ``depth`` below the toe, whatever datum was chosen.

    Where the elevation zero goes is the modeller's choice — toe at 0 with the
    rock at -6, or rock at 0 with the toe at +6 — and both are correct. The
    invariant is the thickness under the toe.
    """
    def predicate(sd):
        pts = ground(sd)
        toe = min(y for _x, y in pts) if pts else None
        got = sd.get("max_depth")
        if toe is None or got is None:
            return False, "no ground surface or no max_depth"
        return abs((toe - float(got)) - depth) <= tol, \
            "max_depth %g is %g below the toe at %g (wanted %g)" % (
                got, toe - float(got), toe, depth)
    return predicate


def circle_tangent_to_base(tol=0.5):
    """One starting circle bottoms on ``max_depth`` — the base the model declares."""
    def predicate(sd):
        base = sd.get("max_depth")
        depths = [c.get("Depth") for c in sd.get("circles") or []]
        hit = [d for d in depths
               if d is not None and base is not None
               and abs(float(d) - float(base)) <= tol]
        return bool(hit), "circle depths %s against a base at %s" % (
            [None if d is None else round(float(d), 2) for d in depths], base)
    return predicate


def has_water():
    def predicate(sd):
        piezo = len(sd.get("piezo_line") or [])
        using = [m.get("name") for m in sd.get("materials") or []
                 if str(m.get("u", "")).lower() in ("piezo", "seep")]
        ok = piezo >= 2 and bool(using)
        return ok, "%d piezometric point(s), read by %s" % (piezo, using or "nobody")
    return predicate


def ponded_load(gamma_water, surface_elev, rel=0.10):
    """A hydrostatic distributed load standing against the slope."""
    def predicate(sd):
        rows = sd.get("dloads") or []
        pressures = [p.get("Normal") for row in rows for p in row
                     if p.get("Normal") is not None]
        if not pressures:
            return False, "no distributed load at all — the water was left out"
        deepest = max(float(v) for v in pressures)
        want = gamma_water * surface_elev
        ok = near(deepest, want, rel=rel)
        return ok, "deepest water pressure %.0f, hydrostatic value %.0f" % (
            deepest, want)
    return predicate


# --------------------------------------------------------------------------- #
# Criterion bundles
# --------------------------------------------------------------------------- #
def discipline(completions=6, tokens=400_000):
    return [S.finished(), S.no_latex(), S.no_exploration(),
            S.completions_at_most(completions), S.tokens_at_most(tokens)]


# --------------------------------------------------------------------------- #
# The scenarios
# --------------------------------------------------------------------------- #
SCENARIOS = []


def add(scenario):
    SCENARIOS.append(scenario)
    return scenario


# --- build from a description ---------------------------------------------- #
add(Scenario(
    "build_desc_simple", "build", None,
    ["Build a model of a 30 ft high slope at 2:1 (horizontal:vertical) in one "
     "soil: total unit weight 125 pcf, c = 400 psf, phi = 20 degrees, no water. "
     "US customary units (ft, psf, pcf). Add starting circles. Don't run "
     "anything yet."],
    discipline(completions=5) + [
        S.model_matches("US customary units", units_are("imperial")),
        S.model_matches("one soil at 125 / 400 / 20",
                        materials_are((125.0, 400.0, 20.0))),
        S.model_matches("30 ft face at 2:1", slope_face_is(30.0, 2.0)),
        S.model_matches("ground extends 2H past toe and crest", extents_beyond(2.0)),
        S.model_matches("starting circles seeded", has_starting_circles(2)),
        S.solves(),
        S.helper_not_used("run_lem"),
    ],
    what="a slope described in one sentence, built and left unsolved"))

add(Scenario(
    "build_desc_layered", "build", None,
    ["Build this in SI units (m, kN/m3, kPa). A 10 m high embankment on a 2:1 "
     "face sits on a 6 m thick soft clay foundation over rock; put the bottom of "
     "the model at the top of the rock. Embankment: gamma = 20, c' = 0, "
     "phi' = 34. Clay: gamma = 17, c = 25, phi = 0. The water table is at the "
     "original ground surface — enter it as a piezometric line and have both "
     "materials read it. Add starting circles, including one tangent to the base "
     "of the clay. Don't run anything yet."],
    discipline(completions=6) + [
        S.model_matches("SI units", units_are("si")),
        S.model_matches("two materials at 20/0/34 and 17/25/0",
                        materials_are((20.0, 0.0, 34.0), (17.0, 25.0, 0.0))),
        S.model_matches("10 m face at 2:1", slope_face_is(10.0, 2.0)),
        S.model_matches("piezometric line the materials read", has_water()),
        S.model_matches("bottom of the model at the top of rock",
                        base_below_toe(6.0)),
        S.model_matches("a circle tangent to the base of the clay",
                        circle_tangent_to_base(tol=0.5)),
        S.solves(),
    ],
    what="two layers, a water table and a declared base, in SI"))

add(Scenario(
    "build_desc_ponded", "build", None,
    ["Build a 20 ft high slope on a 3:1 face in one clay — total unit weight "
     "120 pcf, c = 800 psf, phi = 0 — with a reservoir standing against the face "
     "at 15 ft above the toe. US customary units. Add starting circles. Don't "
     "run anything yet."],
    discipline(completions=6) + [
        S.model_matches("US customary units", units_are("imperial")),
        S.model_matches("one clay at 120 / 800 / 0",
                        materials_are((120.0, 800.0, 0.0))),
        S.model_matches("20 ft face at 3:1", slope_face_is(20.0, 3.0)),
        S.model_matches("the reservoir is a hydrostatic load on the face",
                        ponded_load(62.4, 15.0)),
        S.model_matches("starting circles seeded", has_starting_circles(2)),
        S.solves(),
    ],
    what="standing water — the rule that it is always an explicit load"))


# --- build from a drawing --------------------------------------------------- #
add(Scenario(
    "w1_build_from_image", "build", None,
    [("Build this model. Use the dimensions and properties on the drawing. "
      "Unit system: US customary (ft, psf, pcf). Add a starting circle and run "
      "Spencer with a search.", LEM08_SKETCH)],
    discipline(completions=6) + [
        S.model_matches("US customary units", units_are("imperial")),
        S.model_matches("the two soils on the drawing",
                        materials_are((130.0, 300.0, 37.0), (130.0, 0.0, 37.0))),
        S.model_matches("six geogrid layers, 800 lb/ft, 4 ft above each other",
                        lambda sd: _geogrids(sd, 6, 800.0, 4.0)),
        S.model_matches("240 psf on the crest",
                        lambda sd: _crest_load(sd, 240.0)),
        S.solves("spencer"),
        S.fs_matches(method="spencer"),
        S.claims_grounded(),
        S.fs_close_to(1.587, rel=0.05,
                      label="within 5% of the published Spencer answer"),
    ],
    timeout_s=1200,
    what="W-1's build-from-drawing prompt, unchanged"))

add(Scenario(
    "build_image_surcharge", "build", None,
    [("Build this model from the drawing. Unit system: US customary (ft, psf, "
      "pcf). Add starting circles and stop there — don't run anything.",
      LEM02_SKETCH)],
    discipline(completions=6) + [
        S.model_matches("US customary units", units_are("imperial")),
        S.model_matches("one soil at 125 / 500 / 0",
                        materials_are((125.0, 500.0, 0.0))),
        S.model_matches("a 750 psf strip on the crest",
                        lambda sd: _crest_load(sd, 750.0, width=10.0)),
        S.model_matches("20 ft face at 1:1", slope_face_is(20.0, 1.0)),
        S.solves(),
        S.fs_close_to(0.918, rel=0.10,
                      label="within 10% of the published answer"),
    ],
    timeout_s=1200,
    what="a second drawing — a crest surcharge read off dimension arrows"))


def _geogrids(sd, count, tmax, spacing):
    lines = reinforcement(sd)
    if len(lines) != count:
        return False, "%d reinforcement line(s), wanted %d" % (len(lines), count)
    tmaxes = [line.get("t_max", line.get("Tmax")) for line in lines]
    if not all(near(t, tmax) for t in tmaxes):
        return False, "Tmax values %s, wanted %g" % (tmaxes, tmax)
    ys = sorted(float(line["y1"]) for line in lines)
    steps = [round(b - a, 2) for a, b in zip(ys, ys[1:])]
    if not all(near(s, spacing, absolute=0.2) for s in steps):
        return False, "layer spacings %s, wanted %g" % (steps, spacing)
    return True, "%d layers, %g lb/ft, %g apart" % (count, tmax, spacing)


def _crest_load(sd, pressure, width=None):
    rows = sd.get("dloads") or []
    values = [p.get("Normal") for row in rows for p in row]
    if not values:
        return False, "no distributed load"
    if not any(near(v, pressure, rel=0.02) for v in values if v is not None):
        return False, "load pressures %s, wanted %g" % (values, pressure)
    if width is not None:
        spans = [abs(row[-1]["X"] - row[0]["X"]) for row in rows if len(row) >= 2]
        if not any(near(s, width, rel=0.10) for s in spans):
            return False, "load widths %s, wanted %g" % (spans, width)
    return True, "%g over %s" % (pressure, "%g" % width if width else "the crest")


# --- editing ---------------------------------------------------------------- #
add(Scenario(
    "w1_modify", "edit", REINF,
    ["Change the slope face to 2:1 and rerun the search.",
     "Add a distributed load of 500 psf on the crest from x = 60 to x = 90 and "
     "rerun.",
     "Extend all the reinforcement lines 5 ft to the right and rerun."],
    discipline(completions=6) + [
        S.geometry_source_kept(),
        S.edited("the face is 2:1 after turn 1",
                 lambda before, after: slope_face_is(24.0, 2.0)(after), turn=1),
        S.edited("a second 500 psf block after turn 2",
                 lambda before, after: _second_load(after, 500.0, 60.0, 90.0),
                 turn=2),
        S.edited("every reinforcement line 5 ft longer after turn 3",
                 lambda before, after: _tails_extended(before, after, 5.0),
                 turn=3),
        S.per_turn_fs(),
        S.method_is("spencer"),
        S.claims_grounded(),
        S.reports_warnings(),
    ],
    timeout_s=1200, save_each_turn=True,
    what="W-1's three cumulative edits, unchanged"))


def _second_load(sd, pressure, x1, x2):
    rows = sd.get("dloads") or []
    if len(rows) < 2:
        return False, "%d distributed load block(s), wanted 2" % len(rows)
    for row in rows:
        if (near(row[0]["X"], x1, absolute=0.1)
                and near(row[-1]["X"], x2, absolute=0.1)
                and near(row[0].get("Normal"), pressure)):
            return True, "%g psf from x = %g to %g" % (pressure, x1, x2)
    return False, "no block runs %g to %g at %g psf" % (x1, x2, pressure)


def _tails_extended(before, after, by):
    old = reinforcement(before or {})
    new = reinforcement(after)
    if len(old) != len(new):
        return False, "%d line(s) became %d" % (len(old), len(new))
    moved = [round(float(b["x2"]) - float(a["x2"]), 3) for a, b in zip(old, new)]
    if not all(near(m, by, absolute=0.01) for m in moved):
        return False, "x2 moved by %s, wanted %g each" % (moved, by)
    return True, "every x2 moved %g ft" % by


add(Scenario(
    "edit_geom_polygon", "edit", POLY,
    ["Raise the crest of this slope by 5 ft, keeping the face angle the same, "
     "and rerun the search."],
    discipline(completions=5) + [
        S.geometry_source_kept(),
        S.edited("the crest is 5 ft higher",
                 lambda before, after: _crest_raised(before, after, 5.0)),
        S.fs_matches(),
        S.claims_grounded(),
        S.reports_warnings(),
    ],
    what="the same kind of edit on a polygon-native model"))


def _crest_raised(before, after, by):
    old = max(y for _x, y in ground(before or {})) if before else None
    new = max(y for _x, y in ground(after))
    if old is None:
        return False, "the starting model could not be read"
    return near(new - old, by, absolute=0.05), \
        "crest went from %g to %g (wanted +%g)" % (old, new, by)


add(Scenario(
    "edit_loads", "edit", SURCHARGE,
    ["Add a second distributed load of 300 psf on the crest from x = 40 to "
     "x = 55, then rerun the search and tell me the new factor of safety."],
    discipline(completions=5) + [
        S.edited("a second 300 psf block from x = 40 to 55",
                 lambda before, after: _second_load(after, 300.0, 40.0, 55.0)),
        S.fs_matches(),
        S.method_is(),
        S.claims_grounded(),
    ],
    what="a load added beside one that is already there"))

add(Scenario(
    "edit_materials", "edit", LAYERS,
    ["Raise the cohesion of soil 3 to 600 psf and rerun the search. How much "
     "does the factor of safety move?"],
    discipline(completions=5) + [
        S.edited("soil 3 carries c = 600 and nothing else changed",
                 lambda before, after: _only_c_changed(before, after, "soil 3",
                                                       600.0)),
        S.fs_matches(),
        S.states_the_measurement(),
        S.claims_grounded(),
    ],
    what="one material property, and the change it is worth"))


def _only_c_changed(before, after, name, value):
    got = material(after, name=name).get("c")
    if not near(got, value, rel=0.01):
        return False, "%s c = %s, wanted %g" % (name, got, value)
    for old, new in zip((before or {}).get("materials") or [],
                        after.get("materials") or []):
        if str(old.get("name")) == name:
            continue
        for field in ("gamma", "phi", "c"):
            if not near(old.get(field), new.get(field), absolute=1e-6):
                return False, "%s %s changed too (%s -> %s)" % (
                    old.get("name"), field, old.get(field), new.get(field))
    return True, "%s c = %g, the other seven untouched" % (name, value)


add(Scenario(
    "edit_reinforcement", "edit", GEOGRID,
    ["Add a sixth reinforcement layer 2 ft above the top one — same length, same "
     "Tmax and same pullout lengths as the others — and rerun the search."],
    discipline(completions=5) + [
        S.edited("six layers, the new one 2 ft above the old top",
                 lambda before, after: _layer_added(before, after, 2.0)),
        S.fs_matches(),
        S.claims_grounded(),
    ],
    what="a reinforcement layer added to a set of five"))


def _layer_added(before, after, spacing):
    old = reinforcement(before or {})
    new = reinforcement(after)
    if len(new) != len(old) + 1:
        return False, "%d line(s) became %d, wanted %d" % (
            len(old), len(new), len(old) + 1)
    top_old = max(float(line["y1"]) for line in old)
    top_new = max(float(line["y1"]) for line in new)
    if not near(top_new - top_old, spacing, absolute=0.05):
        return False, "the new layer is at %g, %g above the old top %g" % (
            top_new, top_new - top_old, top_old)
    added = max(new, key=lambda line: float(line["y1"]))
    ref = max(old, key=lambda line: float(line["y1"]))
    for field in ("t_max", "lp1", "lp2"):
        if not near(added.get(field), ref.get(field), rel=0.001):
            return False, "the new layer's %s is %s, the others carry %s" % (
                field, added.get(field), ref.get(field))
    span_new = abs(float(added["x2"]) - float(added["x1"]))
    span_ref = abs(float(ref["x2"]) - float(ref["x1"]))
    if not near(span_new, span_ref, rel=0.02):
        return False, "the new layer is %g long, the others %g" % (span_new,
                                                                   span_ref)
    return True, "six layers, the new one at y = %g" % top_new


add(Scenario(
    "edit_piles", "edit", PILES,
    ["Move the pile wall 5 ft to the right and take its tip 3 ft deeper, then "
     "rerun the search."],
    discipline(completions=5) + [
        S.edited("the wall moved 5 ft right and 3 ft deeper",
                 lambda before, after: _pile_moved(before, after, 5.0, 3.0)),
        S.fs_matches(),
        S.claims_grounded(),
    ],
    what="a pile wall relocated and lengthened"))


def _pile_moved(before, after, dx, deeper):
    old = ((before or {}).get("pile_lines") or [{}])[0]
    new = (after.get("pile_lines") or [{}])[0]
    if not old or not new:
        return False, "the model carries no pile line"
    checks = [("x1", float(old.get("x1", 0)) + dx),
              ("x2", float(old.get("x2", 0)) + dx),
              ("y2", float(old.get("y2", 0)) - deeper)]
    for field, want in checks:
        if not near(new.get(field), want, absolute=0.05):
            return False, "%s = %s, wanted %g" % (field, new.get(field), want)
    return True, "wall at x = %g, tip at %g" % (new.get("x1"), new.get("y2"))


# --- seepage ---------------------------------------------------------------- #
add(Scenario(
    "seep_dam_run", "seep", SEEP_DAM,
    ["Run a steady-state seepage analysis on this dam and tell me the total flow "
     "rate through it."],
    discipline(completions=4) + [
        S.helper_used("run_seep"),
        S.lock_holds(kind="seep", tol=0.05),
        S.flow_matches(tol=0.05),
    ],
    what="the seepage engine, on a zoned dam with its boundary set already on"))

add(Scenario(
    "seep_bc_edit", "seep", SEEP_DAM,
    ["Raise the reservoir 2 m — the upstream specified head goes from 18 to 20 — "
     "and rerun the seepage analysis. How much does the flow rate change?"],
    discipline(completions=6) + [
        S.model_matches("the upstream head is 20 m",
                        lambda sd: _head_is(sd, 20.0)),
        S.helper_used("run_seep"),
        S.flow_matches(tol=0.05),
        S.states_the_measurement(),
    ],
    what="a boundary condition edited, and the change measured"))


def _head_is(sd, want):
    heads = [row.get("head") for row
             in ((sd.get("seepage_bc") or {}).get("specified_heads") or [])]
    ok = any(near(h, want, absolute=0.05) for h in heads if h is not None)
    return ok, "specified heads %s (wanted one at %g)" % (heads, want)


# --- finite element --------------------------------------------------------- #
add(Scenario(
    "fem_ssrm", "fem", SSRM,
    ["Run the strength reduction analysis on this model and report the factor of "
     "safety."],
    discipline(completions=4) + [
        S.helper_used("run_fem"),
        S.warns_about(r"minute|slow|take a while|several minutes",
                      label="says an SSRM run costs minutes"),
        S.ssrm_matches(tol=0.03),
        S.claims_grounded(),
    ],
    timeout_s=1800,
    what="the finite element engine, on the smallest model that has a lock"))

add(Scenario(
    "w1_elastic_fem", "fem", REINF,
    ["Suggest values of Young's modulus and Poisson's ratio for these materials "
     "so I can run a finite element analysis, and explain your choice.",
     "Enter them, build a quadratic mesh at 2 ft, and run the strength reduction "
     "analysis."],
    discipline(completions=6) + [
        S.helper_used("suggest_elastic"),
        S.helper_used("run_fem"),
        S.model_matches("both materials carry a stiffness",
                        lambda sd: _stiffness_entered(sd)),
        S.ssrm_matches(tol=0.03),
        S.claims_grounded(),
    ],
    timeout_s=1800,
    what="W-1's stiffness-then-SSRM pair, unchanged"))


def _stiffness_entered(sd):
    bad = [m.get("name") for m in sd.get("materials") or []
           if not (m.get("E") or 0) > 0 or not (m.get("nu") or 0) > 0]
    if bad:
        return False, "no E / nu on %s" % ", ".join(str(b) for b in bad)
    return True, "E and nu on every material"


# --- parametric ------------------------------------------------------------- #
add(Scenario(
    "w1_sweep_builtin", "sweep", REINF,
    ["Sweep the geogrid Tmax from 500 to 3000 lb/ft in 6 steps with a search at "
     "each step and plot FS against Tmax."],
    discipline(completions=5) + [
        S.helper_any("sensitivity", "design_sweep", "parametric_design"),
        S.claims_grounded(),
        S.mechanism_claims_tested(),
        S.model_unchanged(),
    ],
    timeout_s=1200,
    what="W-1's parametric sweep, unchanged"))

add(Scenario(
    "w1_sweep_adhoc", "sweep", REINF,
    ["Run the analysis with 2, 3, 4, 5 and 6 geogrid layers (removing the top "
     "layers first), searching each time, and tabulate FS against the number of "
     "layers."],
    discipline(completions=6) + [
        S.claims_grounded(),
        S.method_is("spencer"),
        S.mechanism_claims_tested(),
        S.model_unchanged(),
    ],
    timeout_s=1200,
    what="W-1's hand-written sweep, unchanged"))

add(Scenario(
    "sweep_design", "sweep", DESIGN,
    ["What cohesion would this slope need to reach a factor of safety of 1.5?"],
    discipline(completions=5) + [
        S.helper_any("design_sweep", "parametric_design"),
        S.crossing_verified(
            r"(?:cohesion|c)\s*(?:of|=|is|:)?\s*(?:about\s*)?"
            r"([\d,]{2,7}(?:\.\d+)?)\s*psf",
            lambda sd, v: sd["materials"][0].__setitem__("c", v),
            target_fs=1.5, tol=0.03,
            label="the cohesion it reports really gives FS = 1.5"),
        S.claims_grounded(),
    ],
    what="design mode — the value that meets a target"))

add(Scenario(
    "back_analysis", "sweep", HASSIOTIS,
    ["This slope failed. Taking the cohesion as measured, what friction angle "
     "would put it exactly at a factor of safety of 1.0?"],
    discipline(completions=5) + [
        S.helper_any("parametric_back_analysis", "design_sweep",
                     "parametric_design"),
        S.crossing_verified(
            r"(?:phi|friction angle|φ)[^0-9\n]{0,24}(\d{1,2}(?:\.\d+)?)",
            lambda sd, v: sd["materials"][0].__setitem__("phi", v),
            target_fs=1.0, tol=0.03,
            label="the friction angle it reports really gives FS = 1.0"),
        S.claims_grounded(),
    ],
    what="back-analysis — the value that makes the slope limiting"))

# Six completions and twenty-one minutes of solving, measured: a tornado over
# eight materials searches on every leg, and the cap is set to what the work
# actually costs rather than to what a one-material sweep costs.
add(Scenario(
    "tornado", "sweep", LAYERS,
    ["Which inputs does the factor of safety of this model depend on most? Give "
     "me a tornado plot."],
    discipline(completions=6) + [
        S.helper_used("parametric_sweep"),
        S.model_unchanged(),
        S.claims_grounded(),
    ],
    timeout_s=2700,
    what="multi-parameter sensitivity in one plot"))

add(Scenario(
    "reliability_beta", "reliability", RELIABILITY,
    ["Run a reliability analysis on this model. Give me the reliability index "
     "and the probability of failure from the Taylor series method and from "
     "Monte Carlo, and tell me why they differ."],
    discipline(completions=6) + [
        S.helper_used("reliability_taylor"),
        S.helper_used("reliability_mc"),
        S.stated_values_grounded(r"(?:reliability index|beta|β)", "beta"),
        S.stated_values_grounded(r"(?:probability of failure|P_?f)",
                                 "probability of failure"),
        S.warns_about(r"fixed surface|holds the surface|does not re-?search|"
                      r"conditional on",
                      label="says Monte Carlo holds the surface fixed"),
        S.mechanism_claims_tested(),
    ],
    timeout_s=1200,
    what="two reliability engines, and the difference between them"))


# --- questions that change nothing ------------------------------------------ #
add(Scenario(
    "w1_conceptual", "conceptual", REINF,
    ["How does a reliability analysis work in XSLOPE?",
     "How do I decide standard deviations for a reliability analysis if I only "
     "have a few tests?"],
    discipline(completions=4) + [
        S.model_unchanged(),
        S.cites_docs(2),
        S.cites_corpus(),
        S.warns_about(r"lognormal|beta_ln",
                      label="names the index the engine actually reports"),
    ],
    save_after=False,
    what="W-1's two conceptual turns, unchanged"))

add(Scenario(
    "conceptual_rapid", "conceptual", RAPID,
    ["How does the rapid drawdown analysis in XSLOPE work, and what does it need "
     "in the input file?"],
    discipline(completions=3) + [
        S.model_unchanged(),
        S.cites_docs(1),
        S.cites_corpus(),
    ],
    save_after=False,
    what="a formulation question with a worked example behind it"))

add(Scenario(
    "conceptual_search", "conceptual", LAYERS,
    ["How does the automated search for the critical circle work, and when does "
     "it settle on the wrong surface?"],
    discipline(completions=3) + [
        S.model_unchanged(),
        S.cites_docs(1),
        S.cites_corpus(),
    ],
    save_after=False,
    what="a question about the search, answered from the pages"))


# --- diagnosis -------------------------------------------------------------- #
add(Scenario(
    "w1_diagnose", "diagnose", REINF,
    ["This model gives a factor of safety below 1. Can you find what is wrong?"],
    discipline(completions=14) + [
        S.faults_named(),
        S.no_false_accusation(r"sliver", r"degenerate", r"should be redrawn",
                              r"the shell polygon (?:is|should)"),
        S.varied_inputs(3),
        S.states_the_measurement(),
        S.mechanism_claims_tested(),
        S.claims_grounded(),
    ],
    faults=[F.phi_typo("base", 3, 37),
            F.load_decimal_slip(10.0, "240", "2400"),
            F.max_depth_typo(-100, -10)],
    timeout_s=1200,
    what="W-1's broken model — the session the audit says got worse"))

add(Scenario(
    "diagnose_layers", "diagnose", LAYERS,
    ["This model should give a Spencer factor of safety near 1.19 and it "
     "doesn't. Go over the inputs and tell me everything that looks wrong."],
    discipline(completions=14) + [
        S.faults_named(),
        S.no_false_accusation(r"the piezometric line is (?:wrong|missing)",
                              r"the profile lines? (?:are|is) (?:wrong|broken)"),
        S.varied_inputs(2),
        S.states_the_measurement(),
        S.mechanism_claims_tested(),
        S.claims_grounded(),
    ],
    faults=[F.gamma_typo("soil 3", 11.5, 115),
            F.circle_off_slope(0, 400.0)],
    timeout_s=1200,
    what="a decimal slip that moves the answer, beside one that does not"))

add(Scenario(
    "diagnose_dam", "diagnose", DAM,
    ["This dam comes back at a factor of safety around 3.5. The report it came "
     "from says 1.8. Find what is wrong with it."],
    discipline(completions=14) + [
        S.faults_named(),
        S.no_false_accusation(r"the shell (?:material )?is (?:wrong|missing)"),
        S.varied_inputs(2),
        S.states_the_measurement(),
        S.mechanism_claims_tested(),
        S.claims_grounded(),
    ],
    faults=[F.water_deleted(), F.phi_typo("Core", 2, 26)],
    timeout_s=1200,
    what="a deleted water definition, hiding a second fault under it"))


# --- the two behaviors that are neither ------------------------------------- #
def _units_criterion():
    """MPa handed over on a model written in kPa.

    Either answer is acceptable engineering: convert (50 MPa = 50,000 kPa in the
    model's own unit) or refuse to guess and say why. Writing 50 into a kPa field
    is the one that is not.
    """
    from .core import Criterion

    def check(ctx):
        sd = ctx.after() or load_model(ctx.start_model)
        sci = (material(sd, 0).get("hb_sci") if sd else None) or 0.0
        converted = near(sci, 50_000.0, rel=0.02)
        warned = re.search(r"MPa[^\n]{0,120}kPa|kPa[^\n]{0,120}MPa|50[,\s]?000",
                           ctx.prose, re.I)
        if converted:
            return True, "wrote sigma_ci = %g kPa" % sci
        if warned:
            return True, "did not convert but raised it: %r" % warned.group(0)[:48]
        return False, ("sigma_ci = %s on a kPa model, and the answer never "
                       "mentions the unit" % sci)
    return Criterion("catches MPa on a kPa model", check, kind="truth")


add(Scenario(
    "ambiguous_edit", "judgment", DAM,
    ["Move the water table up."],
    discipline(completions=3) + [
        S.asks_rather_than_guesses(),
        S.model_unchanged(),
    ],
    save_after=False,
    what="a request with no number in it"))

add(Scenario(
    "units_trap", "judgment", ROCK,
    ["Switch the rock to the Hoek-Brown option with sigma_ci = 50 MPa, GSI = 60, "
     "mi = 12 and D = 0, then rerun the search."],
    discipline(completions=6) + [
        _units_criterion(),
        S.claims_grounded(),
    ],
    timeout_s=2700,
    what="a strength handed over in MPa on a model that is written in kPa"))


# --- the report ------------------------------------------------------------- #
add(Scenario(
    "w1_report", "report", REINF,
    ["Run Spencer with a search.",
     "Generate the analysis report for this model."],
    discipline(completions=5) + [
        S.fs_matches(method="spencer", on="start"),
        S.wrote_file(".docx", contains=["Spencer", "1.587",
                                        lambda ctx: _input_digest(ctx)]),
        S.claims_grounded(),
    ],
    settings={"report/finalize": False}, timeout_s=1200,
    what="W-1's report pair, unchanged — the .docx must carry the run"))


def _input_digest(ctx):
    from .core import digest
    return digest(ctx.start_model)
