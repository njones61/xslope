"""Interslice tension: one measure, one threshold, one consequence.

Spencer, Morgenstern-Price, the Corps of Engineers method and Lowe & Karafiath
all carry interslice forces, and all four used to judge the tension in them
differently. Morgenstern-Price refused a solution above 30% of its interior
boundaries in tension, force_equilibrium (Corps and Lowe) above 50%, and
Spencer's path refused nothing and reported the same state as a warning. Three
severities on one quantity, so a surface could be answered by one method and
declined by another for reasons that had nothing to do with the surface.

WHY THE CONSEQUENCE IS A REPORT

Soil carries no tension, so any tensile interslice resultant is a departure from
the Mohr-Coulomb model the method was solved with, and the classical remedy is a
tension crack. What the measurements say is that the EXTENT of the tension does
not tell a right answer from a wrong one, so no threshold on it can be a
refusal:

  * the tension-crack sample run with its crack removed — the very case a crack
    exists to handle — reports the LOWEST reading of any model here, 10%, and its
    factor of safety is 8% high;
  * GS-2.26's toe plane reports 69-71% and every method returns 1.35199772881,
    which is SLOPE/W's own answer on the same plane;
  * the tutorial embankment's phi = 0 circles sweep 26-36% across a smooth
    family and every one of them returns Bishop's exact phi = 0 moment factor of
    safety, which no interslice assumption can move.

The case that most needs an engineer's attention sits at the bottom of the
measure and the top of it is all correct answers. So the reading is reported —
which is what SLOPE/W does, and what Duncan & Wright ask the engineer to check —
and refusal is left to the quantities that contradict the strength model instead
of describing the internal forces: a non-positive factor of safety, and base
tension over more than half the slices (solve.MAX_BASE_TENSION_EXTENT).

WHAT IS CHECKED

* GS-2.26: Spencer, both Morgenstern-Price functions, Corps and Lowe all accept
  the plane, agree on the factor of safety to nine significant figures, and all
  four report the same interslice state;
* the counter-case, the tension-crack sample with its crack removed: all four
  report the crest tension, none refuses it, and the reading is below every
  reading in the GS-2.26 leg — the measurement that says no extent threshold can
  separate the two;
* no refusal anywhere names interslice tension, over every surface these legs
  touch;
* the phi = 0 family, where the exact answer is known: 15 circles whose readings
  straddle the old 30% line, with every rigorous answer equal to Bishop's to
  1e-6;
* both facings: a mirrored model reports the same state, since the three paths
  store Z in two different sign conventions;
* three mutations — the 30% refusal restored, the 50% refusal restored, and the
  magnitude floor that used to hide the counter-case from Spencer.

Run directly:  PYTHONPATH=. python3 test/interslice_admissibility_check.py
"""

import os
import re
import warnings

import numpy as np

warnings.filterwarnings('ignore')

from shapely import affinity

from xslope import solve
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

#: Baker's planar benchmark. On a single plane every slice has the same base
#: inclination, so all four side-force methods pose the same equilibrium; the
#: interslice tension belongs to the 76 degree face, not to any assumption.
PLANE_MODEL = os.path.join(REPO, 'docs', 'verification', 'files', 'geostudio',
                           'gs2_26.xlsx')
PLANE_SLICES = 50
PLANE_FS = 1.351997728  # SLOPE/W's own solve of the identical toe plane: 1.352

#: The counter-case: a model whose crest tension a tension crack exists to
#: remove, run with the crack and without it.
CRACK_MODEL = os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_tension_KEY.xlsx')
CRACK_SLICES = 40

#: The tutorial embankment: one phi = 0 material, so moment equilibrium fixes the
#: factor of safety whatever the interslice forces do, and Bishop's answer is the
#: exact one to compare every method against.
PHI0_MODEL = os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_simple_embankment.xlsx')
PHI0_SLICES = 40
PHI0_CIRCLES = [{'Xo': xo, 'Yo': yo, 'R': yo}
                for xo in (2.0, 4.0, 5.28, 6.0, 8.0)
                for yo in (39.0, 42.0, 48.0)]

#: The methods that carry interslice forces and therefore report this measure.
SIDE_FORCE_METHODS = ('spencer', 'mprice', 'corps', 'lowe')

NOTE = 'interslice tension'


def _solve(method, df):
    """(ok, result_or_message) for one method on a copy of the slice frame."""
    d = df.copy()
    if method == 'mprice_constant':
        return d, solve.mprice(d, f_type='constant')
    return d, getattr(solve, method)(d)


def _note_of(result):
    for w in result.get('warnings') or []:
        if w.startswith(NOTE):
            return w
    return None


def _extent(note):
    """The percentage the note quotes, as a float."""
    m = re.search(r'(\d+(?:\.\d+)?)%', note)
    if not m:
        raise AssertionError(f"no percentage in the note: {note!r}")
    return float(m.group(1))


def _plane_slices():
    sd = load_slope_data(PLANE_MODEL)
    ok, res = generate_slices(sd, non_circ=sd['non_circ'], num_slices=PLANE_SLICES)
    if not ok:
        raise AssertionError(f"could not slice {PLANE_MODEL}: {res}")
    return res[0]


def _crack_slices(crack):
    sd = load_slope_data(CRACK_MODEL)
    if not crack:
        sd['tcrack_depth'] = 0.0
        sd['tcrack_water'] = 0.0
    ok, res = generate_slices(sd, circle=sd['circles'][0], num_slices=CRACK_SLICES)
    if not ok:
        raise AssertionError(f"could not slice {CRACK_MODEL}: {res}")
    return res[0]


def _phi0_slices(spec, model=None):
    sd = model if isinstance(model, dict) else load_slope_data(PHI0_MODEL)
    circle = dict(spec)
    circle['Depth'] = circle['Yo'] - circle['R']
    ok, res = generate_slices(sd, circle=circle, num_slices=PHI0_SLICES,
                              check_inputs=False)
    if not ok:
        raise AssertionError(f"could not slice {spec}: {res}")
    return res[0]


# ---------------------------------------------------------------------------


def leg_plane_is_answered_by_every_method():
    """GS-2.26: one equilibrium, one answer, one reported state."""
    fails = []
    df = _plane_slices()
    answers = {}
    notes = {}
    for method in ('spencer', 'janbu', 'mprice', 'mprice_constant', 'corps', 'lowe'):
        _d, (ok, r) = _solve(method, df)
        if not ok:
            fails.append(f"{method} declines the plane: {r}")
            continue
        answers[method] = float(r['FS'])
        if method != 'janbu':
            notes[method] = _note_of(r)
    for method, fs in answers.items():
        if abs(fs - PLANE_FS) / PLANE_FS > 1e-8:
            fails.append(f"{method} returns {fs:.10f} on the plane, against "
                         f"SLOPE/W's {PLANE_FS} on the same plane")
    if answers:
        spread = max(answers.values()) - min(answers.values())
        if spread > 1e-6:
            fails.append(f"the methods disagree on the plane by {spread:.2e}; on a "
                         f"single plane they pose the same equilibrium")
        print(f"  plane        {len(answers)} methods, "
              f"FS {min(answers.values()):.9f}-{max(answers.values()):.9f}")
    for method, note in notes.items():
        if note is None:
            fails.append(f"{method} accepts the plane without reporting its "
                         f"interslice tension")
    if all(notes.values()):
        extents = {m: _extent(n) for m, n in notes.items()}
        if max(extents.values()) - min(extents.values()) > 5.0:
            fails.append(f"the methods read different extents on the same plane: "
                         f"{extents}")
        print(f"  plane state  {min(extents.values()):.0f}-{max(extents.values()):.0f}% "
              f"of interior boundaries in tension, reported by all four")
    return fails


def leg_counter_case_is_reported_not_refused():
    """A tension-cracked crest with no crack modeled: what the reading is worth."""
    fails = []
    with_crack = _crack_slices(True)
    without = _crack_slices(False)
    readings = {}
    for method in SIDE_FORCE_METHODS:
        _dc, (ok_c, rc) = _solve(method, with_crack)
        _dn, (ok_n, rn) = _solve(method, without)
        if not ok_c or not ok_n:
            fails.append(f"{method} could not solve the tension-crack sample: "
                         f"{rc if not ok_c else rn}")
            continue
        note = _note_of(rn)
        if note is None:
            fails.append(f"{method} accepts the uncracked crest without reporting "
                         f"the interslice tension a crack would remove")
            continue
        if _note_of(rc) is not None:
            fails.append(f"{method} reports interslice tension on the CRACKED "
                         f"model, where the crack has removed the tension zone")
        rise = 100 * (rn['FS'] - rc['FS']) / rc['FS']
        readings[method] = (_extent(note), rise)
        print(f"  counter-case {method:<8} {rc['FS']:.4f} -> {rn['FS']:.4f} "
              f"({rise:+.1f}%) at {_extent(note):.0f}% interslice tension")
    if readings:
        worst = max(r[0] for r in readings.values())
        rises = [r[1] for r in readings.values()]
        if min(rises) < 4.0:
            fails.append(f"removing the crack no longer moves the factor of safety "
                         f"({min(rises):.1f}%); the counter-case has gone stale")
        if worst >= 30.0:
            fails.append(f"the uncracked crest now reads {worst:.0f}%, at or above "
                         f"the 30% line the old guard drew — the measurement that "
                         f"says a threshold cannot separate it no longer holds")
    return fails


def leg_no_refusal_names_interslice_tension():
    """Interslice tension is reported. It never refuses a solution."""
    fails = []
    seen = 0
    frames = [_plane_slices(), _crack_slices(False), _crack_slices(True)]
    model = load_slope_data(PHI0_MODEL)
    frames += [_phi0_slices(spec, model) for spec in PHI0_CIRCLES]
    for df in frames:
        for method in ('spencer', 'mprice', 'mprice_constant', 'corps', 'lowe'):
            _d, (ok, r) = _solve(method, df)
            seen += 1
            if not ok and 'interslice tension' in str(r):
                fails.append(f"{method} refused a solution on interslice tension: {r}")
    if not fails:
        print(f"  refusals     {seen} solves, none refused on interslice tension")
    return fails


def leg_phi0_family_shows_the_measure_is_uninformative():
    """At phi = 0 the answer is known, so the measure can be scored against it."""
    fails = []
    model = load_slope_data(PHI0_MODEL)
    extents = []
    worst = 0.0
    for spec in PHI0_CIRCLES:
        df = _phi0_slices(spec, model)
        ok_b, rb = solve.bishop(df.copy())
        if not ok_b:
            continue
        for method in SIDE_FORCE_METHODS:
            _d, (ok, r) = _solve(method, df)
            if not ok:
                continue
            note = _note_of(r)
            if note is not None:
                extents.append(_extent(note))
            if method in ('spencer', 'mprice'):
                err = abs(r['FS'] - rb['FS']) / rb['FS']
                worst = max(worst, err)
                if err > 1e-6:
                    fails.append(
                        f"{method} on {spec} returns {r['FS']:.8f} against Bishop's "
                        f"exact phi = 0 moment factor {rb['FS']:.8f}; the rigorous "
                        f"methods must reproduce it whatever the interslice forces do")
    if extents:
        lo, hi = min(extents), max(extents)
        print(f"  phi = 0      {len(extents)} readings spanning {lo:.0f}-{hi:.0f}% "
              f"interslice tension, every answer within {worst:.1e} of Bishop's exact")
        if hi < 30.0 or lo > 30.0:
            fails.append(f"the phi = 0 family no longer straddles the old 30% line "
                         f"({lo:.0f}-{hi:.0f}%), so it no longer measures what a "
                         f"refusal there would have cost")
    else:
        fails.append("the phi = 0 family reported no interslice tension at all")
    return fails


def leg_both_facings_read_the_same_state():
    """Three paths store Z in two sign conventions. The measure resolves both."""
    fails = []
    sd = load_slope_data(PLANE_MODEL)
    axis = float(max(x for x, _ in sd['ground_surface'].coords))
    flip = lambda g: affinity.scale(g, xfact=-1.0, yfact=1.0, origin=(axis, 0.0))
    mirrored = dict(
        sd,
        ground_surface=flip(sd['ground_surface']),
        domain_polygon=flip(sd['domain_polygon']),
        polygons=[dict(p, polygon=flip(p['polygon'])) for p in sd['polygons']],
        profile_lines=[dict(p, coords=[(2 * axis - x, y) for x, y in p['coords']])
                       for p in sd['profile_lines']],
        non_circ=[dict(p, X=2 * axis - p['X']) for p in reversed(sd['non_circ'])])
    ok, res = generate_slices(mirrored, non_circ=mirrored['non_circ'],
                              num_slices=PLANE_SLICES)
    if not ok:
        fails.append(f"could not slice the mirrored plane: {res}")
        return fails
    left = _plane_slices()
    right = res[0]
    for method in SIDE_FORCE_METHODS:
        _dl, (ok_l, rl) = _solve(method, left)
        _dr, (ok_r, rr) = _solve(method, right)
        if not ok_l or not ok_r:
            fails.append(f"{method} answers one facing and not the other: "
                         f"{rl if not ok_l else rr}")
            continue
        if abs(rl['FS'] - rr['FS']) / rl['FS'] > 1e-6:
            fails.append(f"{method} is not mirror-symmetric: {rl['FS']:.6f} vs "
                         f"{rr['FS']:.6f}")
        nl, nr = _note_of(rl), _note_of(rr)
        if (nl is None) != (nr is None):
            fails.append(f"{method} reports interslice tension on one facing only "
                         f"({nl!r} vs {nr!r}) — the sign convention is leaking")
        elif nl is not None and abs(_extent(nl) - _extent(nr)) > 5.0:
            fails.append(f"{method} reads {_extent(nl):.0f}% one way and "
                         f"{_extent(nr):.0f}% mirrored")
    if not fails:
        print("  facings      mirrored plane: same factor, same reported state")
    return fails


# ---------------------------------------------------------------------------


def _mutation(label, patch, restore, leg, fails):
    patch()
    try:
        caught = leg()
    finally:
        restore()
    if not caught:
        fails.append(f"{label}: the leg passed with the defect in place")
    else:
        print(f"  mutation     {label} -> caught ({len(caught)} failure(s))")


def _refusing(name, limit, message):
    """Re-impose an extent refusal on one method, the way it used to be."""
    original = getattr(solve, name)

    def wrapped(slice_df, *a, **k):
        ok, r = original(slice_df, *a, **k)
        if ok and 'z' in slice_df.columns:
            Z = np.asarray(slice_df['z'].values, dtype=float)[1:]
            if Z.size and float((Z < 0).mean()) > limit:
                return False, message
        return ok, r

    return original, wrapped


def leg_mutations():
    """Each severity this rule replaced must make one of the legs fail."""
    fails = []

    orig_mp, refusing_mp = _refusing(
        'mprice', 0.30, solve.MPRICE_INADMISSIBLE + "(0% of base normals in tension, 33% interslice tension)")
    _mutation("Morgenstern-Price's 30% refusal",
              lambda: setattr(solve, 'mprice', refusing_mp),
              lambda: setattr(solve, 'mprice', orig_mp),
              leg_no_refusal_names_interslice_tension, fails)

    orig_corps, refusing_corps = _refusing(
        'corps', 0.50, "force_equilibrium: inadmissible solution "
                       "(0% of base normals in tension, 71% interslice tension)")
    _mutation("force_equilibrium's 50% refusal",
              lambda: setattr(solve, 'corps', refusing_corps),
              lambda: setattr(solve, 'corps', orig_corps),
              leg_no_refusal_names_interslice_tension, fails)

    original_note = solve.interslice_tension_note

    def floored(Z):
        """The old screen: silent unless the tension reached a tenth of the
        peak compression, which is how the counter-case went unreported."""
        note = original_note(Z)
        if note is None:
            return None
        z = np.asarray(Z, dtype=float)[1:-1]
        if -float(z.min()) <= 0.10 * float(z.max()):
            return None
        return note

    _mutation("the magnitude floor that hid the counter-case",
              lambda: setattr(solve, 'interslice_tension_note', floored),
              lambda: setattr(solve, 'interslice_tension_note', original_note),
              leg_counter_case_is_reported_not_refused, fails)
    return fails


LEGS = [
    ("every method answers the plane, and reports the same state",
     leg_plane_is_answered_by_every_method),
    ("the counter-case is reported, not refused",
     leg_counter_case_is_reported_not_refused),
    ("no refusal names interslice tension",
     leg_no_refusal_names_interslice_tension),
    ("at phi = 0 the measure is scored against the exact answer",
     leg_phi0_family_shows_the_measure_is_uninformative),
    ("both facings read the same state", leg_both_facings_read_the_same_state),
    ("mutations", leg_mutations),
]


def run():
    failures = []
    for label, fn in LEGS:
        print(f"{label}:")
        failures += fn()
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILED:")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print("\nPASS: interslice tension is one measure, reported identically by every "
          "method that carries interslice forces, and it refuses nothing.")


if __name__ == "__main__":
    main()
