"""The task menu — round two of the corpus sweep.

Round one asked one question of every workbook: run it as the file declares it.
That measured the assistant across the INPUT columns, and it measured a single
verb. A model can be run correctly and still be described wrongly, edited on the
wrong geometry source, explained with numbers nobody computed, or reported on in
a document that never carries the run.

This module is the second axis: a menu of eight things to ask about a file, each
with its prompt and its scorers, drawn at random per file from the ones that file
can actually support — a tension-crack question is never asked of a model with no
tension crack — and a stratified sample of the corpus to ask them about, so a run
of sixty files still touches every input column rather than sixty embankments.

Two things are deliberately NOT random. The sample and the task draw are seeded,
so ``--seed 1`` names the same run twice and ``--replay`` re-scores exactly what
was played. And the SPEC inside a task — which material an edit changes, which
fault is planted, which three values a sweep walks — is seeded on the file and
the task alone, never on the run seed: a replay rebuilds the scenario from the
session's own record with no saved state to go stale, and the same file asked the
same question twice is asked it the same way.

A ninth and tenth family live here too, and they open no workbook at all:
building a model from a tutorial's problem drawing, and building one from the
prose and tables a tutorial page describes it in. The truth for both is the
completed workbook that tutorial ships, compared facet by facet — zone areas and
the surface trace rather than vertex lists, since two correct sections can be
drawn with different numbers of points.
"""

from __future__ import annotations

import glob
import json
import math
import os
import random
import re

from . import faults as F
from . import scorers as S
from .core import (Criterion, Scenario, claimed_fs, load_model, numbers_in,
                   repo, solve, solve_variant, strip_code)

# --------------------------------------------------------------------------- #
# Reading a model
# --------------------------------------------------------------------------- #
_UNIT_LENGTH = {"si": "m", "imperial": "ft"}


def _units(sd):
    return "si" if str((sd or {}).get("unit_system") or "").lower() == "si" \
        else "imperial"


def _length_unit(sd):
    return _UNIT_LENGTH[_units(sd)]


def _materials(sd):
    return list((sd or {}).get("materials") or [])


def _number(value, default=None):
    try:
        got = float(value)
    except (TypeError, ValueError):
        return default
    return default if got != got else got


def _ground(sd):
    """The ground surface as ``[(x, y)]``, whichever source built it."""
    gs = (sd or {}).get("ground_surface")
    if gs is None:
        return []
    try:
        return [(float(x), float(y)) for x, y in gs.coords]
    except Exception:
        try:
            return [(float(x), float(y)) for x, y in gs]
        except Exception:
            return []


def _face(sd):
    """``(rise, run_per_rise, x_toe, x_crest)`` of the steepest descending run."""
    pts = _ground(sd)
    if len(pts) < 2:
        return (0.0, 0.0, 0.0, 0.0)
    best = max(range(len(pts) - 1), key=lambda i: abs(pts[i + 1][1] - pts[i][1]))
    (x1, y1), (x2, y2) = pts[best], pts[best + 1]
    rise, run = abs(y2 - y1), abs(x2 - x1)
    lo, hi = (x1, x2) if y1 < y2 else (x2, x1)
    return (rise, (run / rise) if rise else 0.0, lo, hi)


def _one_clear_face(sd):
    """Whether "the slope face" names one segment and not several.

    A dam section has two faces and a benched cut has four, and "flatten the face
    to 3:1" is then a question about which one — an ambiguity the scorer would
    charge to the assistant. The edit is only offered where one descending run
    carries most of the height and no other run comes close.
    """
    pts = _ground(sd)
    if len(pts) < 3:
        return False
    drops = sorted((abs(b[1] - a[1]) for a, b in zip(pts, pts[1:])), reverse=True)
    if drops[0] < 1.0:
        return False
    total = max(y for _x, y in pts) - min(y for _x, y in pts)
    return drops[0] >= 0.6 * total and (len(drops) < 2 or drops[1] <= 0.3 * drops[0])


def _dload_values(sd):
    return [_number(p.get("Normal")) for row in ((sd or {}).get("dloads") or [])
            for p in row if _number(p.get("Normal"))]


def _piezo_span(sd):
    ys = [_number(p[1] if isinstance(p, (list, tuple)) else p.get("y"))
          for p in ((sd or {}).get("piezo_line") or [])]
    return [y for y in ys if y is not None]


def _line_length(line):
    return math.hypot(_number(line.get("x2"), 0.0) - _number(line.get("x1"), 0.0),
                      _number(line.get("y2"), 0.0) - _number(line.get("y1"), 0.0))


def _tidy(value):
    """A number written the way a person would type it into the prompt."""
    if abs(value - round(value)) < 1e-6:
        return "%d" % round(value)
    return ("%.2f" % value).rstrip("0").rstrip(".")


def _spec_rng(case_name, task_name):
    """The RNG a task's spec is drawn from — the file and the task, nothing else.

    Not the run seed. What a task ASKS is a property of the file it is asked
    about, so a replay can rebuild it from the session's name alone, and two runs
    that both drew ``edit_rerun`` on one file asked the same edit of it.
    """
    return random.Random("%s|%s" % (case_name, task_name))


# --------------------------------------------------------------------------- #
# The edits
# --------------------------------------------------------------------------- #
class Edit:
    """One change to ask for, with the words for it and the test that it landed.

    ``expect(before, after)`` is handed two models read back from disk and returns
    ``(ok, reason)``. ``facets`` names the fields it is allowed to move, so
    everything else moving is a separate, failing criterion.
    """

    def __init__(self, name, prompt, expect, facets, apply=None):
        self.name = name
        self.prompt = prompt
        self.expect = expect
        self.facets = set(facets)
        self.apply = apply             # the same edit, made here, for a re-run


def _edit_material_phi(sd, rng):
    pool = [m for m in _materials(sd)
            if (_number(m.get("phi")) or 0) >= 5.0
            and str(m.get("option") or "mc").strip().lower() in ("", "mc")]
    if not pool:
        return None
    mat = rng.choice(pool)
    name = str(mat.get("name") or "").strip()
    old = _number(mat.get("phi"))
    new = round(old - 5.0, 1)

    def expect(before, after):
        got = next((m for m in _materials(after)
                    if str(m.get("name") or "").strip() == name), None)
        if got is None:
            return False, "material %r is no longer in the table" % name
        if abs((_number(got.get("phi")) or -999) - new) > 0.05:
            return False, "%s phi = %s, wanted %s" % (name, got.get("phi"), new)
        others = [m for m in _materials(after)
                  if str(m.get("name") or "").strip() != name]
        was = [m for m in _materials(before)
               if str(m.get("name") or "").strip() != name]
        if S._canon(others) != S._canon(was):
            return False, "another material moved as well"
        return True, "%s phi %s -> %s, nothing else in the table" % (name, old, new)

    return Edit(
        "material_phi",
        "Change the friction angle of material '%s' from %s to %s degrees, then "
        "rerun the analysis as the file declares it and tell me the new factor "
        "of safety." % (name, _tidy(old), _tidy(new)),
        expect, {"materials"},
        apply=lambda model: _set_material(model, name, "phi", new))


def _edit_material_c(sd, rng):
    pool = [m for m in _materials(sd)
            if (_number(m.get("c")) or 0) > 0
            and str(m.get("option") or "mc").strip().lower() in ("", "mc")]
    if not pool:
        return None
    mat = rng.choice(pool)
    name = str(mat.get("name") or "").strip()
    old = _number(mat.get("c"))
    new = round(old * 1.25, 1)

    def expect(before, after):
        got = next((m for m in _materials(after)
                    if str(m.get("name") or "").strip() == name), None)
        if got is None:
            return False, "material %r is no longer in the table" % name
        if abs((_number(got.get("c")) or -1) - new) > max(0.05, 0.002 * new):
            return False, "%s c = %s, wanted %s" % (name, got.get("c"), new)
        return True, "%s c %s -> %s" % (name, old, new)

    return Edit(
        "material_c",
        "Change the cohesion of material '%s' from %s to %s, then rerun the "
        "analysis as the file declares it and tell me the new factor of safety."
        % (name, _tidy(old), _tidy(new)),
        expect, {"materials"},
        apply=lambda model: _set_material(model, name, "c", new))


def _set_material(sd, name, field, value):
    for mat in _materials(sd):
        if str(mat.get("name") or "").strip() == name:
            mat[field] = float(value)
            return
    raise KeyError("no material named %r" % name)


def _edit_dload(sd, rng):
    values = _dload_values(sd)
    # Only a uniform strip. "Change the distributed load to 360" names one number
    # on a model that has one; on a model whose load varies point to point it
    # names nothing, and the answer would be graded on an ambiguity.
    if not values or max(values) - min(values) > 1e-6:
        return None
    old = values[0]
    new = round(old * 1.5, 1)

    def expect(before, after):
        got = _dload_values(after)
        if not got:
            return False, "the distributed load is gone"
        if abs(max(got) - new) > max(0.05, 0.002 * new):
            return False, "load values %s, wanted %s" % (got[:4], new)
        rows_before = [len(r) for r in (before or {}).get("dloads") or []]
        rows_after = [len(r) for r in after.get("dloads") or []]
        if rows_before != rows_after:
            return False, "the load's extent changed too (%s -> %s)" % (
                rows_before, rows_after)
        return True, "the load reads %s (was %s), over the same extent" % (new, old)

    return Edit(
        "dload_value",
        "The distributed load on this model is %s. Change it to %s over the same "
        "extent, then rerun the analysis as the file declares it and tell me the "
        "new factor of safety." % (_tidy(old), _tidy(new)),
        expect, {"dloads"},
        apply=lambda model: _set_dloads(model, new))


def _set_dloads(sd, value):
    for row in sd.get("dloads") or []:
        for point in row:
            if _number(point.get("Normal")):
                point["Normal"] = float(value)


def _edit_piezo(sd, rng):
    ys = _piezo_span(sd)
    if len(ys) < 2:
        return None
    drop = 5.0 if _units(sd) == "imperial" else 2.0
    unit = _length_unit(sd)

    def expect(before, after):
        was, now = _piezo_span(before), _piezo_span(after)
        if len(was) != len(now):
            return False, "the piezometric line has %d point(s), was %d" % (
                len(now), len(was))
        off = [abs((a - b) - drop) for a, b in zip(was, now)]
        if max(off or [9]) > 0.05:
            return False, "the line moved by %s, wanted %g everywhere" % (
                [round(a - b, 2) for a, b in zip(was, now)][:4], drop)
        return True, "every piezometric point is %g %s lower" % (drop, unit)

    return Edit(
        "water_level",
        "Lower the piezometric line by %g %s — every point on it down by the same "
        "amount, the x values unchanged — then rerun the analysis as the file "
        "declares it and tell me the new factor of safety." % (drop, unit),
        expect, {"piezo_line"},
        apply=lambda model: _lower_piezo(model, drop))


def _lower_piezo(sd, drop):
    line = sd.get("piezo_line") or []
    for i, point in enumerate(line):
        if isinstance(point, dict):
            point["y"] = float(point["y"]) - drop
        else:
            line[i] = (float(point[0]), float(point[1]) - drop)


def _edit_reinforcement(sd, rng):
    lines = list(sd.get("reinforcement_lines") or [])
    if not lines:
        return None
    extra = 5.0 if _units(sd) == "imperial" else 2.0
    unit = _length_unit(sd)

    def expect(before, after):
        was = list((before or {}).get("reinforcement_lines") or [])
        now = list(after.get("reinforcement_lines") or [])
        if len(was) != len(now):
            return False, "%d reinforcement line(s), was %d" % (len(now), len(was))
        off = [abs((_line_length(b) - _line_length(a)) - extra)
               for a, b in zip(was, now)]
        if max(off or [9]) > 0.2:
            return False, "lengths grew by %s, wanted %g on every line" % (
                [round(_line_length(b) - _line_length(a), 2)
                 for a, b in zip(was, now)][:4], extra)
        return True, "every reinforcement line is %g %s longer" % (extra, unit)

    return Edit(
        "reinforcement_length",
        "Extend every reinforcement line %g %s further into the slope — move the "
        "embedded end back along the line, leave the face end where it is — then "
        "rerun the analysis as the file declares it and tell me the new factor of "
        "safety." % (extra, unit),
        expect, {"reinforcement_lines"},
        apply=lambda model: _extend_reinforcement(model, extra))


def _extend_reinforcement(sd, extra):
    for line in sd.get("reinforcement_lines") or []:
        length = _line_length(line)
        if length <= 0:
            continue
        dx = (_number(line.get("x2"), 0.0) - _number(line.get("x1"), 0.0)) / length
        dy = (_number(line.get("y2"), 0.0) - _number(line.get("y1"), 0.0)) / length
        line["x2"] = _number(line.get("x2"), 0.0) + dx * extra
        line["y2"] = _number(line.get("y2"), 0.0) + dy * extra


def _edit_pile_spacing(sd, rng):
    piles = list(sd.get("pile_lines") or [])
    spacings = [_number(p.get("S")) for p in piles if _number(p.get("S"))]
    if not spacings or max(spacings) - min(spacings) > 1e-6:
        return None
    old = spacings[0]
    new = round(old * 0.75, 2)

    def expect(before, after):
        got = [_number(p.get("S")) for p in after.get("pile_lines") or []]
        if not got:
            return False, "the pile row is gone"
        if any(g is None or abs(g - new) > 0.02 for g in got):
            return False, "pile spacings %s, wanted %s on every row" % (got, new)
        return True, "pile spacing %s -> %s" % (old, new)

    return Edit(
        "pile_spacing",
        "The piles are spaced %s apart. Tighten the spacing to %s on every pile "
        "row, then rerun the analysis as the file declares it and tell me the new "
        "factor of safety." % (_tidy(old), _tidy(new)),
        expect, {"pile_lines"},
        apply=lambda model: [p.__setitem__("S", new)
                             for p in model.get("pile_lines") or []])


def _edit_slope_angle(sd, rng):
    if not _one_clear_face(sd):
        return None
    rise, ratio, _toe, _crest = _face(sd)
    if rise <= 0 or ratio <= 0:
        return None
    old = round(ratio, 2)
    new = round(old + 0.5, 2)
    source = "profile lines" if sd.get("profile_lines") else "polygons"

    def expect(before, after):
        got_rise, got_ratio, _lo, _hi = _face(after)
        if abs(got_ratio - new) > 0.1:
            return False, "the face is %.2f:1, wanted %.2f:1" % (got_ratio, new)
        if abs(got_rise - rise) > max(0.2, 0.03 * rise):
            return False, "the face height moved (%.2f, was %.2f)" % (got_rise,
                                                                      rise)
        return True, "the face reads %.2f:1 at the same height" % got_ratio

    return Edit(
        "slope_angle",
        "Flatten the slope face from %s:1 to %s:1 (horizontal:vertical), keeping "
        "the toe and the crest elevation where they are, then rerun the analysis "
        "as the file declares it and tell me the new factor of safety."
        % (_tidy(old), _tidy(new)),
        expect, {"profile_lines", "polygons", "ground_surface"} if source ==
        "profile lines" else {"polygons", "ground_surface"})
    # No independent ``apply``: re-drawing a face here would be a second
    # implementation of the thing under test, and the two could disagree while
    # both were right. This edit is measured on the section the session saved.


def edit_answer_independent(case, spec, tol=0.01):
    """The number the answer reports matches the same edit made HERE.

    The other half of "make this change and re-run". Re-solving the workbook the
    session saved says the answer is consistent with what it wrote down; this
    says the answer is right about the change that was ASKED for — a session that
    edited the wrong material and reported a number consistent with that passes
    the first and fails this.

    Not every edit has an independent form. Re-drawing a slope face here would be
    a second implementation of the thing being measured, so the geometry edit is
    left to the saved-model comparison alone and this criterion is not asked of it.
    """
    def check(ctx):
        claims = claimed_fs(ctx.final_prose) or claimed_fs(ctx.prose)
        if not claims:
            return False, "no factor of safety stated in the answer"
        run = solve_variant(case.path, spec.apply, "edit:%s" % spec.name,
                            method=case.method, search=True)
        if run.get("FS") is None:
            return False, "the same edit made here does not solve: %s" \
                          % run.get("error")
        best = min(claims, key=lambda v: abs(v - run["FS"]))
        return abs(best - run["FS"]) <= tol, \
            "stated %s vs the same edit made here %s" % (best,
                                                         round(run["FS"], 4))
    return Criterion("the answer matches the edit made independently", check,
                     kind="truth")


#: Every edit, rarest first — the draw prefers one the file alone can support
#: over one every model could be asked, so a run spreads across the edit kinds.
_EDITS = (_edit_pile_spacing, _edit_reinforcement, _edit_piezo, _edit_dload,
          _edit_slope_angle, _edit_material_phi, _edit_material_c)


def edit_for(case):
    """The edit this file will be asked for, or ``None`` if none of them fit."""
    rng = _spec_rng(case.name, "edit_rerun")
    built = [make(case.model, rng) for make in _EDITS]
    offers = [e for e in built if e is not None]
    if not offers:
        return None
    # Rarest-first, then a seeded pick among the two rarest that fit: a file with
    # piles is usually asked about its piles, but not always, or the pile edit
    # would be the only one ever measured on a pile model.
    return rng.choice(offers[:2])


# --------------------------------------------------------------------------- #
# The sweeps
# --------------------------------------------------------------------------- #
class Sweep:
    """One parameter to walk over three values, with the words and the truth."""

    def __init__(self, name, prompt, values, apply, label):
        self.name = name
        self.prompt = prompt
        self.values = list(values)
        self.apply = apply             # apply(model, value)
        self.label = label


def sweep_for(case):
    """The parameter this file will be swept on, or ``None``."""
    sd = case.model or {}
    rng = _spec_rng(case.name, "small_sweep")
    mats = [m for m in _materials(sd)
            if str(m.get("option") or "mc").strip().lower() in ("", "mc")]
    with_c = [m for m in mats if (_number(m.get("c")) or 0) > 0]
    with_phi = [m for m in mats if (_number(m.get("phi")) or 0) >= 5]
    offers = []
    if with_phi:
        mat = rng.choice(with_phi)
        name = str(mat.get("name") or "").strip()
        base = _number(mat.get("phi"))
        values = [round(base * f, 1) for f in (0.8, 1.0, 1.2)]
        offers.append(Sweep(
            "phi", "Sweep the friction angle of material '%s' over %s, %s and %s "
                   "degrees and tabulate the factor of safety at each — the "
                   "analysis the file declares, once per value."
                   % (name, _tidy(values[0]), _tidy(values[1]), _tidy(values[2])),
            values,
            lambda model, v, n=name: _set_material(model, n, "phi", v),
            "phi of %s" % name))
    if with_c:
        mat = rng.choice(with_c)
        name = str(mat.get("name") or "").strip()
        base = _number(mat.get("c"))
        values = [round(base * f, 1) for f in (0.75, 1.0, 1.25)]
        offers.append(Sweep(
            "c", "Sweep the cohesion of material '%s' over %s, %s and %s and "
                 "tabulate the factor of safety at each — the analysis the file "
                 "declares, once per value."
                 % (name, _tidy(values[0]), _tidy(values[1]), _tidy(values[2])),
            values,
            lambda model, v, n=name: _set_material(model, n, "c", v),
            "c of %s" % name))
    if not offers:
        return None
    return rng.choice(offers)


# --------------------------------------------------------------------------- #
# The planted faults
# --------------------------------------------------------------------------- #
def fault_for(case):
    """One fault this file can carry, and the words for the question about it."""
    sd = case.model or {}
    rng = _spec_rng(case.name, "find_fault")
    offers = []
    mats = _materials(sd)
    typo = [m for m in mats if (_number(m.get("phi")) or 0) >= 20
            and str(m.get("option") or "mc").strip().lower() in ("", "mc")]
    if typo:
        mat = rng.choice(typo)
        right = int(round(_number(mat.get("phi"))))
        offers.append(F.phi_typo(str(mat.get("name") or "").strip(),
                                 right % 10 or 3, right))
    heavy = [m for m in mats if (_number(m.get("gamma")) or 0) > 0]
    if heavy:
        mat = rng.choice(heavy)
        right = _number(mat.get("gamma"))
        offers.append(F.gamma_typo(str(mat.get("name") or "").strip(),
                                   round(right / 10.0, 2), right))
    if _dload_values(sd):
        old = max(_dload_values(sd))
        offers.append(F.load_decimal_slip(10.0, _tidy(old * 10), _tidy(old)))
    if sd.get("piezo_line"):
        offers.append(F.water_deleted())
    if sd.get("circles"):
        offers.append(F.circle_off_slope(0, 400.0))
    if not offers:
        return None
    return rng.choice(offers)


# --------------------------------------------------------------------------- #
# describe — the facts a file carries, and the ones it does not
# --------------------------------------------------------------------------- #
DESCRIBE_PROMPT = ("Describe this model: geometry, materials and their strengths, "
                   "water, loads, reinforcement or piles, and the failure surface "
                   "it declares.")


def _facts(sd):
    """``[(label, [numbers any one of which proves it was read])]`` for a model."""
    out = []
    mats = _materials(sd)
    if mats:
        out.append(("the number of materials", [float(len(mats))]))
    for mat in mats[:6]:
        name = str(mat.get("name") or "?").strip()
        values = [v for v in (_number(mat.get("gamma")), _number(mat.get("c")),
                              _number(mat.get("phi"))) if v]
        if values:
            out.append(("material %r strengths" % name, values))
    ys = _piezo_span(sd)
    if ys:
        out.append(("the water elevation", [max(ys), min(ys)]))
    loads = _dload_values(sd)
    if loads:
        out.append(("the distributed load", [max(loads)]))
    if sd.get("reinforcement_lines"):
        out.append(("the reinforcement count",
                    [float(len(sd["reinforcement_lines"]))]))
    if sd.get("pile_lines"):
        out.append(("the pile rows", [float(len(sd["pile_lines"]))]))
    if _number(sd.get("tcrack_depth")):
        out.append(("the tension crack depth", [_number(sd["tcrack_depth"])]))
    if _number(sd.get("k_seismic")):
        out.append(("the seismic coefficient", [_number(sd["k_seismic"])]))
    if sd.get("circles"):
        out.append(("the starting circles", [float(len(sd["circles"]))]))
    return out


def describes_the_facts(case, allow_missing=1):
    """Every fact the file carries is in the description, with the file's numbers.

    Each fact passes when ONE of its numbers appears in the prose at the precision
    the answer wrote it to — a material is read correctly when its unit weight,
    its cohesion or its friction angle is quoted, not when all three are. One
    fact may be missed: an answer that describes a model well and rounds a single
    number away is not the failure this is looking for.
    """
    def check(ctx):
        prose = strip_code(ctx.prose)
        pool = numbers_in(prose)
        wanted = _facts(case.model)
        if not wanted:
            return True, "the file carries no fact worth quoting"
        missed = [label for label, values in wanted
                  if not any(_appears(v, pool) for v in values)]
        if len(missed) > allow_missing:
            return False, "%d of %d fact(s) not in the description: %s" % (
                len(missed), len(wanted), "; ".join(missed[:4]))
        return True, "%d of %d fact(s) quoted with the file's own numbers" % (
            len(wanted) - len(missed), len(wanted))
    return Criterion("describes what the file holds", check, kind="truth")


def _appears(value, pool):
    """Whether ``value`` is one of ``pool``, at the precision each was written."""
    from .core import rounds_to
    return any(rounds_to(seen, value) or rounds_to(value, seen) for seen in pool)


#: Features an answer can claim, and the words that claim them. A model without
#: one of these must not have it described — the invention this catches is the
#: assistant filling a gap with what a slope usually has.
_FEATURES = (
    ("reinforcement", r"\b(reinforcement|geogrid|geotextile|soil nail|tieback|"
                      r"anchor)s?\b", lambda sd: bool(sd.get("reinforcement_lines"))),
    ("piles", r"\bpiles?\b|\bpile row\b", lambda sd: bool(sd.get("pile_lines"))),
    ("a tension crack", r"\btension crack\b",
     lambda sd: bool(_number(sd.get("tcrack_depth")))),
    ("seismic loading", r"\bseismic\b|\bpseudo-?static\b",
     lambda sd: bool(_number(sd.get("k_seismic")))),
    ("a distributed load", r"\b(distributed load|surcharge)\b",
     lambda sd: bool(_dload_values(sd))),
    ("a piezometric line", r"\bpiezometric\b",
     lambda sd: bool(sd.get("piezo_line"))),
    ("line loads", r"\bline loads?\b", lambda sd: bool(sd.get("line_loads"))),
)
#: A feature word inside a negative sentence is a correct statement of ABSENCE.
_NEGATED = re.compile(r"\b(no|not|none|never|without|neither|nor|absent|"
                      r"lacks?|lacking|free of|zero|isn't|aren't|doesn't|"
                      r"there is no|there are no)\b", re.I)


def invents_nothing(case):
    """No feature the file lacks is described as present."""
    def check(ctx):
        # A markdown heading is a LABEL for the section under it, not a claim
        # about the model: "### Loads, reinforcement, piles" over a paragraph
        # saying there are none is a correct answer, and reading the heading as
        # an assertion failed every well-organized description in the round.
        prose = re.sub(r"^[ \t]{0,3}#{1,6}[^\n]*$", "", strip_code(ctx.prose),
                       flags=re.M)
        invented = []
        for label, pattern, present in _FEATURES:
            if present(case.model or {}):
                continue
            for hit in re.finditer(pattern, prose, re.I):
                if not _NEGATED.search(_sentence(prose, hit.start())):
                    invented.append(label)
                    break
        if invented:
            return False, "the file has no %s, and the answer says it does" % (
                ", ".join(sorted(set(invented))[:3]))
        return True, "nothing the file lacks was described"
    return Criterion("invents no input the file lacks", check, kind="truth")


def _sentence(text, index):
    """The sentence around ``index`` — where a negation would have to sit."""
    start = max(text.rfind(".", 0, index), text.rfind("\n", 0, index),
                text.rfind(";", 0, index)) + 1
    stop = min((p for p in (text.find(".", index), text.find("\n", index),
                            len(text)) if p != -1), default=len(text))
    return text[start:stop]


# --------------------------------------------------------------------------- #
# find_fault — what else was accused
# --------------------------------------------------------------------------- #
#: The input fields a diagnosis can accuse, and the words that name each one.
_FIELDS = (
    ("friction angle", r"\bphi\b|φ|friction angle"),
    ("unit weight", r"unit weight|\bgamma\b|γ"),
    ("cohesion", r"\bcohesion\b|\bc'?\s*=|\bc\s+value"),
    ("the load", r"distributed load|surcharge|\bdload"),
    ("the water", r"piezometric|water table|phreatic|pore[- ]?(?:water )?pressure"),
    ("the base of the model", r"max_?depth|bottom of the model"),
    ("a starting circle", r"starting circle|circle cent"),
    ("the geometry", r"profile line|polygon|ground surface"),
    ("the reinforcement", r"reinforcement|geogrid|\bt_?max\b"),
    ("the piles", r"\bpiles?\b|pile spacing"),
    ("the mesh", r"\bmesh\b|element size"),
)
#: The words that turn a mention of a field into an ACCUSATION of it.
_ACCUSES = re.compile(r"\b(wrong|incorrect|implausib|not credible|typo|mis-?typ|"
                      r"transpos|should be|error|suspect|off by|too (?:high|low|"
                      r"large|small)|missing|deleted|dropped)\b", re.I)
#: Which field each planted fault is about, keyed by the start of its id.
_FAULT_FIELD = {"phi_": "friction angle", "gamma_": "unit weight",
                "load_x": "the load", "water_deleted": "the water",
                "circle_off_slope": "a starting circle",
                "max_depth": "the base of the model"}


def _accused(prose):
    """The input fields the answer says are wrong."""
    found = []
    for line in re.split(r"(?<=[.\n])\s+", strip_code(prose)):
        if not _ACCUSES.search(line):
            continue
        for label, pattern in _FIELDS:
            if re.search(pattern, line, re.I) and label not in found:
                found.append(label)
    return found


def accuses_little_else(fault):
    """The planted field is accused, and at most one other input is.

    A diagnosis that lists every input in the model as suspect has found nothing —
    it has spread the answer wide enough to contain it. One further concern is
    ordinary engineering; a spray of five is the failure the ranked-table
    discipline in the brief exists to prevent.
    """
    planted = next((v for k, v in _FAULT_FIELD.items() if fault.id.startswith(k)),
                   None)

    def check(ctx):
        accused = _accused(ctx.prose)
        others = [a for a in accused if a != planted]
        if len(others) > 1:
            return False, "accused %d input(s) besides %s: %s" % (
                len(others), planted or "the planted one", ", ".join(others[:4]))
        return True, "accused %s%s" % (
            ", ".join(accused) or "nothing",
            " (one beyond the planted fault)" if others else "")
    return Criterion("accuses at most one input beyond the fault", check,
                     kind="truth")


# --------------------------------------------------------------------------- #
# explain_result — which slices carry the surface
# --------------------------------------------------------------------------- #
EXPLAIN_PROMPT = ("Why is the critical surface where it is, and which slices "
                  "contribute the most driving and resisting moment?")


#: ``(driving order, resisting order, slice count)`` per (file digest, method).
#: Ranking costs a solve, and a replay re-scores the same file; without this the
#: search on a slow section is paid again every time the scorer is fixed.
_RANKINGS = {}


def _ranks_searched(case):
    """Whether the file's own published answer is a SEARCH or one surface.

    Several benchmark files publish a ``single_circle`` value: the circle IS the
    problem, and running an automated search on one of them answers a different
    question — on VP30a it wanders to a sliver at FS 0.157 where the declared
    circle gives 1.19, and takes three minutes to do it. The ranking follows what
    the file declares, the same way the sweep's ground truth does.
    """
    types = {t.get("type") for t in case.tags}
    if types & {"circular_search", "noncircular_search"}:
        return True
    if types & {"single_circle", "single_noncirc"}:
        return False
    return False


def _moment_ranking(case):
    """Per-slice driving and resisting moment, computed here from the slice table.

    The solver returns a factor of safety and a slice table, not a moment column,
    so the ranking is made from the columns it DOES return, by the ordinary
    method's own statement of them: driving is ``W sin α`` and resisting is
    ``c·Δl + (W cos α − u·Δl) tan φ``, both times the radius, which is common to
    every slice and drops out of a ranking. Returns
    ``(driving_order, resisting_order, slice_count)`` — slice numbers, largest
    first — or ``(None, None, reason)``.
    """
    from .core import digest

    key = "%s|%s" % (digest(case.path), case.method)
    if key in _RANKINGS:
        return _RANKINGS[key]
    run = solve_variant(case.path, None, "", method=case.method,
                        search=_ranks_searched(case), want_slices=True)
    df = run.get("slice_df")
    if df is None or run.get("FS") is None:
        out = (None, None, (run.get("error") or "no slice table"))
        _RANKINGS[key] = out
        return out
    import numpy as np

    alpha = np.radians(df["alpha"].to_numpy(dtype=float))
    w = df["w"].to_numpy(dtype=float)
    dl = df["dl"].to_numpy(dtype=float)
    c = df["c"].to_numpy(dtype=float)
    phi = np.radians(df["phi"].to_numpy(dtype=float))
    u = df["u"].to_numpy(dtype=float)
    driving = w * np.sin(alpha)
    resisting = c * dl + np.maximum(w * np.cos(alpha) - u * dl, 0.0) * np.tan(phi)
    numbers = df["slice #"].to_numpy()
    order_d = [int(numbers[i]) for i in np.argsort(-driving)]
    order_r = [int(numbers[i]) for i in np.argsort(-resisting)]
    _RANKINGS[key] = (order_d, order_r, len(numbers))
    return _RANKINGS[key]


_SLICE_REF = re.compile(r"\bslices?\s*#?\s*(\d+)(?:\s*(?:-|–|to|through|and)\s*"
                        r"#?\s*(\d+))?", re.I)


def explains_the_slices(case):
    """The slices the answer names really are the ones carrying the surface.

    Every slice number the answer puts forward is checked against a ranking made
    here from the solver's own slice table. The band is the top three, widened to
    a tenth of the slices where the model is cut finer than thirty — a 200-slice
    surface has no meaningful "third largest" and scoring one there would measure
    the discretization.

    An answer that names ZONES rather than slice numbers is judged the same way,
    on the material the leading slices sit in.
    """
    def check(ctx):
        order_d, order_r, count = _moment_ranking(case)
        if order_d is None:
            return False, "no independent run to rank the slices against: %s" % count
        band = max(3, int(0.1 * count))
        top = set(order_d[:band]) | set(order_r[:band])
        named = []
        for hit in _SLICE_REF.finditer(strip_code(ctx.prose)):
            lo = int(hit.group(1))
            hi = int(hit.group(2)) if hit.group(2) else lo
            named += list(range(min(lo, hi), min(max(lo, hi), count) + 1))
        if not named:
            return False, ("the answer names no slice; the leading ones are "
                           "driving %s, resisting %s"
                           % (order_d[:3], order_r[:3]))
        hits = [n for n in named if n in top]
        if not hits:
            return False, ("named slice(s) %s; the top %d are driving %s, "
                           "resisting %s" % (sorted(set(named))[:6], band,
                                             order_d[:band], order_r[:band]))
        return True, ("named %d slice(s), %d of them in the leading %d "
                      "(driving %s, resisting %s)"
                      % (len(set(named)), len(set(hits)), band, order_d[:3],
                         order_r[:3]))
    return Criterion("the slices it names are the leading ones", check,
                     kind="truth")


# --------------------------------------------------------------------------- #
# compare_methods
# --------------------------------------------------------------------------- #
#: Three methods that answer on any model: the two that are always defensible and
#: the one the file itself declares, so a comparison always includes its own.
_COMPARE = ("bishop", "spencer", "janbu")


def compare_methods_of(case):
    methods = [case.method] + [m for m in _COMPARE if m != case.method]
    return methods[:3]


def _published(case, method):
    """The docs' own published answer for one method on this file, or ``None``.

    A file can publish both a circular and a non-circular answer — LEM-5's weak
    layer does — and they are answers to two different questions. The tag whose
    surface family is the one the model itself declares is the one that belongs
    beside a run made as the file declares it.
    """
    circular = not (not (case.model or {}).get("circular")
                    and (case.model or {}).get("non_circ"))
    family = ("circular_search", "single_circle") if circular \
        else ("noncircular_search", "single_noncirc")
    hits = [t for t in case.tags
            if str(t.get("method") or "").lower() == method
            and t.get("expected_fs") is not None]
    for pool in ([t for t in hits if t.get("type") in family], hits):
        for tag in pool:
            # A tag may carry the key with no value; the tag runner's own default
            # applies then, not ``None``.
            return float(tag["expected_fs"]), float(tag.get("tolerance") or 0.005)
    return None


_METHOD_WORD = {"mprice": r"m-?p(?:rice)?|morgenstern[ -]price", "oms": r"oms|ordinary",
                "bishop": r"bishop", "spencer": r"spencer", "janbu": r"janbu",
                "corps": r"corps", "lowe": r"lowe"}


def methods_tabulated(case, methods):
    """Each method's number in the answer matches an independent run of it.

    Where the docs publish a value for that method on this file, that published
    number IS the independent run — the regression suite makes it on every pass —
    and it is used rather than paying for a second search.
    """
    def check(ctx):
        prose = strip_code(ctx.final_prose or ctx.prose)
        rows, bad = [], []
        for method in methods:
            stated = _stated_for(prose, method)
            known = _published(case, method)
            if known is not None:
                want, tol = known
                label = "published"
            else:
                run = solve(case.path, method=method, search=True)
                if run.get("FS") is None:
                    bad.append("%s: no independent answer (%s)"
                               % (method, run.get("error")))
                    continue
                want, tol, label = run["FS"], 0.01, "re-run"
            if not stated:
                bad.append("%s: the answer states no value (%s %s)"
                           % (method, label, round(want, 4)))
                continue
            best = min(stated, key=lambda v: abs(v - want))
            rows.append("%s %s/%s" % (method, best, round(want, 4)))
            if abs(best - want) > max(tol, 0.005):
                bad.append("%s: stated %s, %s %s" % (method, best, label,
                                                     round(want, 4)))
        if bad:
            return False, "; ".join(bad[:3])
        if len(rows) < len(methods):
            return False, "only %d of %d methods answered" % (len(rows),
                                                              len(methods))
        return True, "stated/independent " + ", ".join(rows)
    return Criterion("every method's value matches an independent run", check,
                     kind="truth")


def _stated_for(prose, method):
    """The numbers the answer puts against one method — its table row, or the
    numbers in the sentence that names it."""
    word = _METHOD_WORD.get(method, re.escape(method))
    out = []
    for line in prose.splitlines():
        if not re.search(word, line, re.I):
            continue
        if line.strip().startswith("|"):
            out += [v for v in numbers_in(line) if 0.05 < v < 50]
        else:
            out += [v for v in numbers_in(line) if 0.05 < v < 50]
    return out


# --------------------------------------------------------------------------- #
# small_sweep
# --------------------------------------------------------------------------- #
def sweep_reproduced(case, spec, at_least=2):
    """Two of the three rows the answer tabulates match a run made here.

    Two rather than three because a sweep helper and a bare re-solve can settle on
    different critical circles at the ends of a range, and the question this asks
    is whether the table was computed rather than composed.
    """
    def check(ctx):
        prose = strip_code(ctx.final_prose or ctx.prose)
        stated = claimed_fs(prose) + _all_table_numbers(prose)
        if not stated:
            return False, "the answer tabulates nothing"
        rows, matched = [], 0
        for value in spec.values:
            # The prompt asks for "the analysis the file declares, once per
            # value", and the file declares a surface but not whether to search
            # it. Both readings are made here, and a row counts as reproduced
            # when it matches either — scoring the declared surface against a
            # search measured the ambiguity rather than the table.
            got = None
            for search in (False, True):
                run = solve_variant(case.path,
                                    lambda model, v=value: spec.apply(model, v),
                                    "%s=%s" % (spec.name, value),
                                    method=case.method, search=search)
                if run.get("FS") is None:
                    got = got or run
                    continue
                best = min(stated, key=lambda v: abs(v - run["FS"]))
                got = run
                if abs(best - run["FS"]) <= 0.01:
                    matched += 1
                    rows.append("%s -> %s (stated %s, %s)" % (
                        value, round(run["FS"], 4), best,
                        "as declared" if not search else "searched"))
                    break
            else:
                if got is None or got.get("FS") is None:
                    rows.append("%s: no run (%s)"
                                % (value, (got or {}).get("error")))
                else:
                    best = min(stated, key=lambda v: abs(v - got["FS"]))
                    rows.append("%s -> %s (stated %s)"
                                % (value, round(got["FS"], 4), best))
        if matched < at_least:
            return False, "%d of %d rows reproduced: %s" % (
                matched, len(spec.values), "; ".join(rows[:3]))
        return True, "%d of %d rows reproduced: %s" % (
            matched, len(spec.values), "; ".join(rows))
    return Criterion("the swept values were computed", check, kind="truth")


def _all_table_numbers(prose):
    return [v for line in prose.splitlines() if line.strip().startswith("|")
            for v in numbers_in(line) if 0.05 < v < 50]


# --------------------------------------------------------------------------- #
# The tasks
# --------------------------------------------------------------------------- #
class Task:
    """One thing to ask about a file: whether it fits, and what it becomes."""

    def __init__(self, name, fit, make, completions, cost, what):
        self.name = name
        self.fit = fit                 # fit(case) -> truthy when it can be asked
        self.make = make               # make(case) -> Scenario
        self.completions = completions
        self.cost = cost               # rough completions, for the estimate
        self.what = what

    def eligible(self, case):
        try:
            return bool(self.fit(case))
        except Exception:
            return False

    def __repr__(self):
        return "<Task %s>" % self.name


def _lem(case):
    return case.kind == "lem" and bool((case.model or {}).get("circles")
                                       or (case.model or {}).get("non_circ"))


def _lem_solvable(case):
    """A stability model this harness can re-solve on its own.

    Everything that compares a reported number with an independent run needs to
    make that run HERE, through ``run_lem_analysis``. A model whose materials read
    their pore pressures from a seepage solution has no field until that solve is
    staged, and driving the stability entry point straight is refused by the input
    checks — which is a fact about the run the scorer is making, not about the
    answer. ``run_declared`` still asks those files everything it asked in round
    one: its ground truth goes through the tag runner, which stages the flow.
    """
    return _lem(case) and not case.stages_seep


def _named(case, task):
    return "%s__%s" % (case.name, task)


def _describe(case):
    return Scenario(
        _named(case, "describe"), case.primary, case.path, [DESCRIBE_PROMPT],
        [describes_the_facts(case), invents_nothing(case), S.model_unchanged(),
         S.no_latex(), S.no_exploration(), S.completions_at_most(3)],
        what="describe — %s" % ", ".join(case.classes),
        timeout_s=min(case.timeout_s, 600), save_after=False)


def _run_declared(case):
    """Round one's question, asked again — the same prompt and the same criteria.

    Renamed to carry the task suffix every other scenario in a round-two run
    carries, so the scorecard can group by task and a replay can read the task
    back off the name. What is asked and what is scored are unchanged, which is
    what keeps the two rounds' numbers comparable.
    """
    from . import corpus

    built = corpus.scenario_for(case)
    built.name = _named(case, "run_declared")
    return built


def _edit_rerun(case):
    spec = edit_for(case)
    if spec is None:
        return None
    criteria = [S.edited("the edit asked for was made", spec.expect),
                S.nothing_else_changed(*spec.facets),
                S.no_geometry_warning(),
                S.geometry_source_kept(),
                S.fs_matches(tol=0.01, method=case.method, search=True,
                             on="saved"),
                S.claims_grounded(),
                S.no_exploration(),
                S.completions_at_most(6)]
    if spec.apply is not None:
        criteria.insert(4, edit_answer_independent(case, spec))
    return Scenario(
        _named(case, "edit_rerun"), case.primary, case.path, [spec.prompt],
        criteria, what="edit_rerun (%s) — %s" % (spec.name,
                                                 ", ".join(case.classes)),
        timeout_s=max(case.timeout_s, 900))


def _explain_result(case):
    from . import corpus

    return Scenario(
        _named(case, "explain_result"), case.primary, case.path,
        [corpus.PROMPT, EXPLAIN_PROMPT],
        [explains_the_slices(case),
         S.numbers_grounded_near(r"driving|resisting|moment", "the moments"),
         corpus.reported_value_matches(case),
         S.no_exploration(), S.completions_at_most(5)],
        what="explain_result — %s" % ", ".join(case.classes),
        timeout_s=max(case.timeout_s, 900), save_after=False)


def _compare_methods(case):
    methods = compare_methods_of(case)
    prompt = ("Run %s and %s on this model — each with its own search, as each "
              "method's own critical surface — and tabulate the factor of safety "
              "for each." % (", ".join(m.upper() for m in methods[:-1]),
                             methods[-1].upper()))
    return Scenario(
        _named(case, "compare_methods"), case.primary, case.path, [prompt],
        [methods_tabulated(case, methods), S.claims_grounded(),
         S.no_exploration(), S.completions_at_most(6)],
        what="compare_methods (%s) — %s" % ("/".join(methods),
                                            ", ".join(case.classes)),
        timeout_s=max(case.timeout_s, 1800), save_after=False)


def _small_sweep(case):
    spec = sweep_for(case)
    if spec is None:
        return None
    return Scenario(
        _named(case, "small_sweep"), case.primary, case.path, [spec.prompt],
        [sweep_reproduced(case, spec), S.claims_grounded(),
         S.helper_any("sensitivity", "design_sweep", "run_lem"),
         S.no_exploration(), S.completions_at_most(6)],
        what="small_sweep (%s) — %s" % (spec.label, ", ".join(case.classes)),
        timeout_s=max(case.timeout_s, 1800), save_after=False)


def _find_fault(case):
    fault = fault_for(case)
    if fault is None:
        return None
    return Scenario(
        _named(case, "find_fault"), case.primary, case.path,
        ["Something in this model is wrong — the answer it gives is not the one "
         "this section should give. Go over the inputs and tell me what you find."],
        [S.faults_named(),
         accuses_little_else(fault),
         invents_nothing(case),
         S.varied_inputs(2),
         S.states_the_measurement(),
         S.claims_grounded(),
         S.no_exploration(),
         S.completions_at_most(14)],
        faults=[fault],
        what="find_fault (%s) — %s" % (fault.id, ", ".join(case.classes)),
        timeout_s=max(case.timeout_s, 1200))


def _input_digest(ctx):
    """The hash of the model the session opened — what the report stamps on it."""
    from .core import digest
    return digest(ctx.start_model)


def _report(case):
    return Scenario(
        _named(case, "report"), case.primary, case.path,
        ["Run the analysis as the file declares it, then generate the analysis "
         "report for this model."],
        [S.wrote_file(".docx", contains=["Results", _input_digest]),
         S.claims_grounded(), S.no_exploration(), S.completions_at_most(5)],
        settings={"report/finalize": False},
        what="report — %s" % ", ".join(case.classes),
        timeout_s=max(case.timeout_s, 1800))


#: The menu. ``cost`` is the completions a turn of this kind has been costing,
#: and is what the per-file estimate is built from.
TASKS = [
    Task("describe", lambda c: True, _describe, 3, 2,
         "read the model back in words, with its own numbers"),
    Task("run_declared", lambda c: True, _run_declared, 4, 3,
         "run it as the file declares it"),
    Task("edit_rerun",
         lambda c: _lem_solvable(c) and edit_for(c) is not None,
         _edit_rerun, 6, 5, "one named change, made and re-solved"),
    Task("explain_result", _lem_solvable, _explain_result, 5, 6,
         "which slices carry the surface"),
    Task("compare_methods", _lem_solvable, _compare_methods, 6, 5,
         "three methods, each on its own search"),
    Task("small_sweep",
         lambda c: _lem_solvable(c) and sweep_for(c) is not None,
         _small_sweep, 6, 6, "one parameter over three values"),
    Task("find_fault", lambda c: _lem(c) and fault_for(c) is not None,
         _find_fault, 14, 10, "a planted fault, found or missed"),
    Task("report", _lem, _report, 5, 4, "the .docx, with the run in it"),
]
TASKS_BY_NAME = {t.name: t for t in TASKS}
#: ``run_declared`` is always eligible, so a file drawn into the sample always has
#: at least the round-one question asked of it and the two rounds stay comparable.
ALWAYS = "run_declared"


def eligible_tasks(case):
    return [t for t in TASKS if t.eligible(case)]


def draw_tasks(case, k, rng):
    """``k`` tasks for one file — ``run_declared`` always among them."""
    fits = eligible_tasks(case)
    always = TASKS_BY_NAME[ALWAYS]
    rest = [t for t in fits if t.name != ALWAYS]
    rng.shuffle(rest)
    chosen = ([always] if always in fits else []) + rest
    return chosen[:max(1, int(k))]


def scenario_for(case, task):
    """The scenario for one (file, task) pair, or ``None`` where it does not fit."""
    return task.make(case)


# --------------------------------------------------------------------------- #
# The stratified sample
# --------------------------------------------------------------------------- #
def swept_in(run_dir):
    """The case names a finished run already holds a session for.

    Read from that run's own ``scorecard.json`` where it has one, and from the
    directory names otherwise — a run stopped at its budget has the directories
    and may not have the card.
    """
    names = set()
    card = os.path.join(run_dir, "scorecard.json")
    if os.path.exists(card):
        try:
            with open(card, encoding="utf-8") as fh:
                stored = json.load(fh)
            names = {str(r.get("scenario") or "") for r in stored.get("results") or []}
        except Exception:
            names = set()
    if not names and os.path.isdir(run_dir):
        names = {name for name in os.listdir(run_dir)
                 if os.path.exists(os.path.join(run_dir, name, "session.json"))}
    # A round-two name carries its task; a round-one name is the file alone.
    return {n.split("__", 1)[0] for n in names if n}


def sample(every, n, seed, exclude=(), min_per_class=2):
    """``n`` files, drawn so every input class the corpus can fill is filled.

    Two passes. The first walks the classes rarest-first and takes files until
    each has ``min_per_class`` in the sample — counting the ones already drawn for
    another class, since a file with piles also carries profile lines and there is
    no reason to buy that column twice. The second fills the rest at random.

    ``n`` is a floor, not a cap. A class with two files in the corpus and none in
    the sample is a column the run says nothing about, and quietly dropping one to
    respect a round number would put the arithmetic above the measurement — so a
    sample too small to be stratified comes back larger than asked, and says so.
    """
    from .corpus import CLASS_ORDER

    rng = random.Random(seed)
    skip = set(exclude or ())
    pool = [c for c in every if c.loads and c.name not in skip]
    order = list(CLASS_ORDER)
    order += sorted({cls for c in pool for cls in c.classes} - set(order))

    chosen, taken, counts = [], set(), {}

    def take(case):
        chosen.append(case)
        taken.add(case.name)
        for cls in case.classes:
            counts[cls] = counts.get(cls, 0) + 1

    for label in order:
        available = [c for c in pool if label in c.classes and c.name not in taken]
        need = min_per_class - counts.get(label, 0)
        if need > 0 and available:
            for case in rng.sample(available, min(need, len(available))):
                take(case)
    quota = len(chosen)
    rest = [c for c in pool if c.name not in taken]
    rng.shuffle(rest)
    for case in rest[:max(0, int(n) - len(chosen))]:
        take(case)
    return chosen, {"pool": len(pool), "quota": quota, "asked": int(n),
                    "drawn": len(chosen), "min_per_class": min_per_class}


def plan(cases, tasks_per_file, seed):
    """``[(case, task)]`` for a sample — the seeded draw, in a stable order."""
    rng = random.Random("%s|tasks" % seed)
    out = []
    for case in cases:
        for task in draw_tasks(case, tasks_per_file, rng):
            out.append((case, task))
    return out


# --------------------------------------------------------------------------- #
# Building a model — from a drawing, and from a description
# --------------------------------------------------------------------------- #
#: Where the tutorials keep their problem drawings, and the pattern that marks one.
SKETCH_DIR = "docs/tutorials/images"
SKETCH_GLOB = "*_problem_sketch.png"
#: The headings a tutorial page states its problem under. A page's model is
#: described under exactly one of these, ahead of the build instructions.
_PROBLEM_HEADINGS = (r"## The problem", r"## The slope", r"## The dam",
                     r"## The section", r"## Part A\b")
#: The glance box that names a page's files, and a markdown link inside it.
_GLANCE = re.compile(r'<div class="tgm-model"[^>]*>(.*?)</div>', re.S)
_LINK = re.compile(r"\[([^\]]+\.xlsx)\]\(([^)]+\.xlsx)\)")


class Build:
    """One model to build from nothing, and the finished workbook it is judged on.

    ``truth`` is the completed model the tutorial ships. It is never opened by the
    session — it is the answer key, loaded here.
    """

    def __init__(self, key, kind, page, prompt, truth, image=None):
        self.key = key
        self.kind = kind               # "drawing" or "description"
        self.page = page
        self.prompt = prompt
        self.path = truth
        self.image = image
        self.model = load_model(truth)
        self.name = "build_%s__%s" % ("img" if kind == "drawing" else "txt", key)
        self.classes = ["build (%s)" % kind]
        self.primary = self.classes[0]
        self.kind_label = "build"

    @property
    def method(self):
        from .core import declared_method
        return declared_method(self.model)


def _pages():
    return sorted(glob.glob(repo("docs", "tutorials", "*.md")))


def _completed_model(page):
    """The workbook a tutorial's glance box calls its completed model.

    A page names two files: the STARTER it begins from and the COMPLETED model it
    ends on. The completed one is the answer key; the starter is a partial model
    and scoring a from-scratch build against it would mark every input the page
    adds as wrong. A ``_start`` link is skipped for that reason, and a page that
    names only a starter has no answer key and is left out.

    A page that names SEVERAL completed models is left out too. LEM-7 works two
    slopes and FEM-3 works a pile row and then a wall; "build this model" has no
    single answer there, and scoring a build of one against the other would
    measure the pairing rather than the assistant.
    """
    text = open(page, encoding="utf-8").read()
    found = []
    for block in _GLANCE.findall(text):
        for _label, target in _LINK.findall(block):
            if os.path.basename(target).lower().endswith("_start.xlsx"):
                continue
            path = os.path.normpath(os.path.join(os.path.dirname(page), target))
            if os.path.exists(path) and path not in found:
                found.append(path)
    return found[0] if len(found) == 1 else None


def _problem_text(page):
    """The section a tutorial states its problem in, as it is written."""
    text = open(page, encoding="utf-8").read()
    for heading in _PROBLEM_HEADINGS:
        found = re.search(r"^%s.*$" % heading, text, re.M)
        if not found:
            continue
        body = text[found.end():]
        stop = re.search(r"^## ", body, re.M)
        body = body[:stop.start()] if stop else body
        # Images are the page's figures, not the problem; a link to one reaches
        # the assistant as a path it cannot open.
        body = re.sub(r"!\[[^\]]*\]\([^)]*\)(\{[^}]*\})?", "", body)
        body = re.sub(r"\n{3,}", "\n\n", body).strip()
        if len(body) > 400:
            return body
    return None


def _sketch_for(page):
    """The problem drawing a tutorial page shows, if it shows one."""
    text = open(page, encoding="utf-8").read()
    found = re.search(r"([a-z0-9_]*_problem_sketch\.png)", text)
    if not found:
        return None
    path = repo(SKETCH_DIR, found.group(1))
    return path if os.path.exists(path) else None


#: The answer key for a drawing on a page whose glance box names SEVERAL completed
#: models. ``_completed_model`` cannot choose between them and the page's prose
#: states more than one problem, but the DRAWING is of one slope; these pairings
#: say which. Keyed by the drawing, so a page showing two of them (FEM-3 draws a
#: pile row and then a wall) contributes a build case for each.
DRAWING_KEYS = {
    "lem07_problem_sketch.png":
        ("docs", "lem", "files", "xslope_baker_clay.xlsx"),
    "lem10_problem_sketch.png":
        ("docs", "lem", "files", "xslope_mult_min_KEY.xlsx"),
    "fem03_piles_problem_sketch.png":
        ("docs", "lem", "files", "xslope_piles.xlsx"),
    "fem03_wall_problem_sketch.png":
        ("docs", "tutorials", "files", "xslope_pile_wall.xlsx"),
}


def _sketches_on(page):
    """Every problem drawing a page shows, in the order it shows them."""
    text = open(page, encoding="utf-8").read()
    out = []
    for name in re.findall(r"([a-z0-9_]*_problem_sketch\.png)", text):
        path = repo(SKETCH_DIR, name)
        if os.path.exists(path) and path not in out:
            out.append(path)
    return out


def _keyed_drawings(page):
    """``[(drawing, answer key)]`` for this page, from :data:`DRAWING_KEYS`."""
    out = []
    for sketch in _sketches_on(page):
        parts = DRAWING_KEYS.get(os.path.basename(sketch))
        if not parts:
            continue
        truth = repo(*parts)
        if os.path.exists(truth):
            out.append((sketch, truth))
    return out


def _drawing_key(page_key, sketch, count):
    """The case name for one drawing — the page, unless the page shows several."""
    if count < 2:
        return page_key
    return re.sub(r"_problem_sketch\.png$", "", os.path.basename(sketch))


def _unit_words(sd):
    return ("SI (m, kN/m3, kPa)" if _units(sd) == "si"
            else "US customary (ft, psf, pcf)")


def builds():
    """Every build case the tutorials support — drawings first, then descriptions.

    A page contributes a drawing case when it shows a problem sketch and an answer
    key can be named for it, and a description case when it states its problem in
    prose long enough to build from and names a single completed model. Most pages
    give both, and they are two different questions about the same slope. A page
    naming several completed models has no single answer for its prose, so it
    contributes drawings only — one per :data:`DRAWING_KEYS` pairing.
    """
    out, drawn = [], set()
    for page in _pages():
        key = os.path.splitext(os.path.basename(page))[0]
        pairs = _keyed_drawings(page)
        truth = _completed_model(page)
        if truth and not pairs:
            sketch = _sketch_for(page)
            if sketch:
                pairs = [(sketch, truth)]
        for sketch, answer in pairs:
            # Two pages can show one drawing of one model — COMBO-1 opens on
            # SEEP-2's dam. Building it twice would buy the same measurement twice.
            if (sketch, answer) in drawn:
                continue
            model = load_model(answer)
            if model is None:
                continue
            drawn.add((sketch, answer))
            method = (str(model.get("lem_method") or "") or "spencer").lower()
            out.append(Build(
                _drawing_key(key, sketch, len(pairs)), "drawing", page,
                "Build this model from the drawing. Use the dimensions and "
                "properties shown. Unit system: %s. Add a starting circle where "
                "the drawing suggests one and run %s with a search."
                % (_unit_words(model), method.upper()),
                answer, image=sketch))
        if not truth:
            continue
        model = load_model(truth)
        if model is None:
            continue
        method = (str(model.get("lem_method") or "") or "spencer").lower()
        problem = _problem_text(page)
        if problem:
            out.append(Build(
                key, "description", page,
                "Build this model in XSLOPE. Unit system: %s. Add a starting "
                "circle and run %s with a search.\n\n%s"
                % (_unit_words(model), method.upper(), problem),
                truth))
    return out


_BY_NAME = {}


def builds_by_name():
    """:func:`builds`, indexed and memoized — replay looks a name up per session,
    and rebuilding the family would reload thirty workbooks each time."""
    if not _BY_NAME:
        _BY_NAME.update({b.name: b for b in builds()})
    return _BY_NAME


# -- what a built model is judged on ---------------------------------------- #
def _zone_areas(sd):
    """``[area]`` per material zone, largest first — the geometry, comparably.

    Vertex lists are the wrong thing to compare: two correct sections can be drawn
    with different numbers of points, and one drawn with an extra collinear vertex
    would fail a vertex check while being the same slope. Areas per zone and the
    surface trace are the same section however it was entered.
    """
    out = []
    for entry in (sd or {}).get("polygons") or []:
        # The loader hands each zone back as {"polygon": <Polygon>, "mat_id": …},
        # not as the polygon itself. Reading .area off the dict raised, the
        # fallback then iterated the dict's KEYS, and every model measured zero
        # zones — so the zone half of the section check passed on everything.
        poly = entry.get("polygon") if isinstance(entry, dict) else entry
        try:
            out.append(abs(float(poly.area)))
            continue
        except Exception:
            pass
        try:
            pts = [(float(x), float(y)) for x, y in poly]
            out.append(abs(sum(pts[i][0] * pts[i - 1][1] - pts[i - 1][0] * pts[i][1]
                               for i in range(len(pts)))) / 2.0)
        except Exception:
            pass
    return sorted(out, reverse=True)


def _trace(sd, samples=25):
    """The ground surface sampled at even x — one number per station."""
    pts = _ground(sd)
    if len(pts) < 2:
        return []
    pts = sorted(pts)
    x0, x1 = pts[0][0], pts[-1][0]
    if x1 <= x0:
        return []
    out = []
    for i in range(samples):
        x = x0 + (x1 - x0) * i / (samples - 1.0)
        for (ax, ay), (bx, by) in zip(pts, pts[1:]):
            if ax <= x <= bx and bx > ax:
                out.append(ay + (by - ay) * (x - ax) / (bx - ax))
                break
        else:
            out.append(pts[-1][1])
    return out


def geometry_like(build, rel=0.05):
    """The section built is the section the tutorial ships — by area and by trace.

    Two measurements, both scale-free. The zone areas say the right amount of each
    material is in the section; the surface trace, sampled at even stations and
    compared against the height of the slope, says the ground goes where it should.
    """
    def check(ctx):
        got = ctx.after() or load_model(ctx.saved() or "")
        if got is None:
            return False, "the session saved no workbook"
        want = build.model
        wanted, made = _zone_areas(want), _zone_areas(got)
        if len(made) != len(wanted):
            return False, "%d zone(s), the model ships %d" % (len(made),
                                                              len(wanted))
        off = [abs(a - b) / max(1.0, b) for a, b in zip(made, wanted)]
        if off and max(off) > rel:
            return False, "zone areas %s vs %s (worst %.0f%% off)" % (
                [round(a) for a in made], [round(b) for b in wanted],
                100 * max(off))
        want_trace, got_trace = _trace(want), _trace(got)
        if len(want_trace) != len(got_trace) or not want_trace:
            return False, "the ground surface could not be compared"
        height = max(want_trace) - min(want_trace) or 1.0
        # A drawing gives dimensions, not an elevation datum, and neither does a
        # description: a section entered with its toe at el 10 rather than el 0 is
        # the same section. The best constant offset is removed before the residual
        # is measured, so a translation costs nothing and a shape error still shows.
        # Without this the check reported the datum — every "off by 10.00" was a
        # section that matched the shipped one exactly, ten feet higher.
        offset = sum(a - b for a, b in zip(want_trace, got_trace)) / len(want_trace)
        worst = max(abs(a - b - offset) for a, b in zip(want_trace, got_trace))
        if worst > rel * height:
            return False, "the ground surface is off by %.2f (%.0f%% of the rise) " \
                          "after the %.2f datum shift is removed" \
                          % (worst, 100 * worst / height, offset)
        return True, "%d zone(s) within %.0f%%, surface within %.2f" % (
            len(made), 100 * (max(off) if off else 0), worst)
    return Criterion("the section matches the shipped model", check, kind="truth")


def materials_like(build, rel=0.05):
    """Every material in the shipped model has a counterpart in the built one."""
    def check(ctx):
        got = ctx.after() or load_model(ctx.saved() or "")
        if got is None:
            return False, "the session saved no workbook"
        wanted = _materials(build.model)
        made = _materials(got)
        if len(made) != len(wanted):
            return False, "%d material(s), the model ships %d" % (len(made),
                                                                  len(wanted))
        unmatched = list(made)
        for want in wanted:
            hit = next((m for m in unmatched if _same_material(m, want, rel)), None)
            if hit is None:
                return False, "nothing matches %s (γ %s, c %s, φ %s)" % (
                    want.get("name"), want.get("gamma"), want.get("c"),
                    want.get("phi"))
            unmatched.remove(hit)
        return True, "%d material(s) match within %.0f%%" % (len(wanted),
                                                             100 * rel)
    return Criterion("the materials match the shipped model", check, kind="truth")


def _same_material(got, want, rel):
    """Whether one built material is the one the answer key states.

    Two things the plain c/phi comparison got wrong. A key whose strength model is
    NOT Mohr-Coulomb states its envelope in fields c and phi do not carry — a
    Hoek-Brown row reads c = 0, phi = 0 whatever the rock is — so the model itself
    is what is compared there. And a key that states no strength at all (a
    seepage-only row: no option, gamma, c and phi all blank) says nothing to match
    against, so anything sound matches it.
    """
    want_option = str(want.get("option") or "").strip().lower()
    if want_option not in ("", "mc"):
        if str(got.get("option") or "").strip().lower() != want_option:
            return False
        fields = (("gamma", None),)
    else:
        stated = [f for f in ("gamma", "c", "phi")
                  if _number(want.get(f)) not in (None, 0.0)]
        if not stated:
            return True
        fields = (("gamma", None), ("c", 1.0), ("phi", 0.6))
    for field, absolute in fields:
        a, b = _number(got.get(field)), _number(want.get(field))
        if a is None and b is None:
            continue
        if a is None or b is None:
            return False
        tol = absolute if absolute is not None else 0.0
        if abs(a - b) > max(tol, rel * max(1.0, abs(b))):
            return False
    return True


def features_like(build):
    """The counted inputs — loads, reinforcement, piles, circles — are all there."""
    def check(ctx):
        got = ctx.after() or load_model(ctx.saved() or "")
        if got is None:
            return False, "the session saved no workbook"
        rows, bad = [], []
        for label, key in (("distributed loads", "dloads"),
                           ("line loads", "line_loads"),
                           ("reinforcement lines", "reinforcement_lines"),
                           ("pile rows", "pile_lines"),
                           ("starting circles", "circles")):
            want = len(build.model.get(key) or [])
            made = len(got.get(key) or [])
            if not want and not made:
                continue
            rows.append("%s %d/%d" % (label, made, want))
            if key == "circles":
                if want and not made:
                    bad.append("no starting circle was added")
                continue
            if made != want:
                bad.append("%d %s, the model ships %d" % (made, label, want))
        if bad:
            return False, "; ".join(bad[:3])
        return True, ", ".join(rows) or "the model carries no counted input"
    return Criterion("loads, reinforcement and piles are all there", check,
                     kind="truth")


def water_like(build):
    """A model whose answer key defines water has water — however it is entered."""
    def check(ctx):
        want = build.model
        defines = bool(want.get("piezo_line")
                       or (want.get("seepage_bc") or {}).get("specified_heads")
                       or any(str(m.get("u") or "none").lower() != "none"
                              for m in _materials(want)))
        if not defines:
            return True, "the shipped model defines no water"
        got = ctx.after() or load_model(ctx.saved() or "")
        if got is None:
            return False, "the session saved no workbook"
        made = bool(got.get("piezo_line")
                    or (got.get("seepage_bc") or {}).get("specified_heads")
                    or any(str(m.get("u") or "none").lower() != "none"
                           for m in _materials(got)))
        if not made:
            return False, "the shipped model defines water and the built one does not"
        ys = _piezo_span(got)
        return True, "water is defined%s" % (
            " (piezometric line at %.1f–%.1f)" % (min(ys), max(ys)) if ys else "")
    return Criterion("the water definition is there", check, kind="truth")


def build_fs_close(build, rel=0.10):
    """The built model answers near the value the tutorial publishes for its own.

    The tolerance is wide on purpose. A drawing under-determines a section — the
    band a reinforced face is drawn over, the exact station a load ends at — so a
    build that is right in every input a reader could have taken off the drawing
    still lands a few percent away, and a tight band here would measure the
    drawing.
    """
    def check(ctx):
        from .core import lock

        published = lock(build.path, kind="circular_search", method=build.method)
        if published is None:
            published = lock(build.path, kind="circular_search")
        if published is None:
            run = solve(build.path, method=build.method, search=True)
            published = run.get("FS")
        if published is None:
            return True, "the shipped model publishes no value to compare with"
        path = ctx.saved() or ctx.start_model
        if not path:
            return False, "the session saved no workbook"
        made = solve(path, method=build.method, search=True)
        if made.get("FS") is None:
            return False, "the built model does not solve: %s" % made.get("error")
        off = abs(made["FS"] - published) / max(1e-9, abs(published))
        return off <= rel, "built model solves at %s; the shipped one at %s " \
                           "(%.1f%% apart)" % (round(made["FS"], 4),
                                               round(published, 4), 100 * off)
    return Criterion("within %d%% of the shipped model's answer" % int(rel * 100),
                     check, kind="truth")


def build_scenario(build):
    """The scenario that builds one model and scores it against the shipped one."""
    turns = [(build.prompt, build.image)] if build.image else [build.prompt]
    return Scenario(
        build.name, build.primary, None, turns,
        [geometry_like(build), materials_like(build), features_like(build),
         water_like(build), build_fs_close(build),
         S.preflight_clean(), S.claims_grounded(),
         S.no_latex(), S.no_exploration(), S.completions_at_most(8)],
        what="build from a %s — %s" % (build.kind,
                                       os.path.basename(build.page)),
        timeout_s=1800)


# --------------------------------------------------------------------------- #
# Replay
# --------------------------------------------------------------------------- #
def resolve(name, session):
    """Rebuild the scenario a recorded session played, from its name.

    Nothing about a scenario is stored in the run beyond the name it was played
    under and the workbook it opened, and nothing needs to be: the file names the
    case, the suffix names the task, and the spec inside the task is seeded on
    those two. A replay of a run therefore re-scores exactly what was played, and
    a scorer fixed afterwards re-grades it without buying it again.
    """
    from . import corpus

    if name.startswith("build_"):
        build = builds_by_name().get(name)
        return build_scenario(build) if build else None
    case_name, _sep, task_name = name.rpartition("__")
    task = TASKS_BY_NAME.get(task_name)
    path = session.get("lock_model") or session.get("start_model")
    if task is None or not path:
        return None
    case = corpus.Case(path)
    return scenario_for(case, task)


def row_for(name, session=None):
    """The case-shaped record the scorecard groups a finished scenario by."""
    from . import corpus

    if name.startswith("build_"):
        build = builds_by_name().get(name)
        if build is None:
            return None
        return Row(name, build.classes, build.primary, "build", build.path,
                   "build_%s" % build.kind)
    case_name, _sep, task_name = name.rpartition("__")
    path = (session or {}).get("lock_model")
    if not path:
        return None
    case = corpus.Case(path)
    return Row(name, case.classes, case.primary, case.kind, case.path,
               task_name)


class Row:
    """What :func:`tools.assistant_scenarios.corpus.render` reads off a case.

    A scenario in a sampled run is a (file, task) pair, so the scorecard cannot
    index its cases by file name the way round one does. This is the same shape
    the class renderer reads — name, classes, primary class, engine, path — keyed
    by the SCENARIO name, with the task carried alongside for the task table.
    """

    def __init__(self, name, classes, primary, kind, path, task):
        self.name = name
        self.classes = list(classes)
        self.primary = primary
        self.kind = kind
        self.path = path
        self.task = task


# --------------------------------------------------------------------------- #
# The scorecard
# --------------------------------------------------------------------------- #
def render(results, meta):
    """The corpus scorecard, with a table by TASK added beside the class one."""
    from . import corpus

    text = corpus.render(results, meta)
    index = {r.name: r for r in meta.get("cases") or []}
    tasks = []
    for row in results:
        task = getattr(index.get(row["scenario"]), "task", None) \
            or row["scenario"].rpartition("__")[2]
        tasks.append((task, row))
    names = []
    for task, _row in tasks:
        if task not in names:
            names.append(task)
    lines = ["", "## By task", "",
             "| Task | Run | Pass | Fail | Failing files |",
             "| --- | ---: | ---: | ---: | :--- |"]
    for task in names:
        rows = [r for t, r in tasks if t == task]
        bad = [r["scenario"] for r in rows if not r["pass"]]
        lines.append("| `%s` | %d | %d | %d | %s |" % (
            task, len(rows), len(rows) - len(bad), len(bad),
            ", ".join("`%s`" % n for n in bad[:5]) + (" …" if len(bad) > 5 else "")
            or "—"))
    draw = meta.get("draw") or {}
    if draw:
        lines += ["", "The sample: %d file(s) drawn from %d eligible (%d of them "
                      "to fill the class quota of %d each), %d scenario(s) "
                      "played, seed `%s`."
                  % (draw.get("drawn", 0), draw.get("pool", 0),
                     draw.get("quota", 0), draw.get("min_per_class", 2),
                     len(results), meta.get("seed"))]
    block = "\n".join(lines) + "\n"
    marker = "\n## By criterion"
    if marker in text:
        return text.replace(marker, block + marker, 1)
    return text + block

