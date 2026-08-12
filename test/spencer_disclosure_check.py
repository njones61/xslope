"""Surfaces a method admits no admissible solution on — classified, and disclosed.

Spencer's method assumes every interslice force acts at the SAME inclination
theta, which turns equilibrium into two equations in (F, theta). On some trial
surfaces that system has no root the method would accept. No starting guess and
no iteration count reaches one, and until this check's subject shipped, two
things followed from that.

First, the failure read as the solver's: "did not converge within the maximum
number of iterations" says the iteration ran out, when what happened is about the
method. Second — the one that reaches a printed factor of safety — the SEARCH
scored the unanswered trial exactly as it scores a circle that misses the model:
fs_fail, dropped, never mentioned. So a search could report its minimum as
converged while surfaces LOWER than it went unanswered, and nothing said so.

THE MECHANISM, AND WHY IT HAS A CLOSED FORM

At phi = 0, on a circular surface carrying no loads, the system decouples. With
n_i = W_i sin(a_i) - c_i dl_i / F,

    R1 = sum_i n_i sec(a_i - theta)
    R2 = (Xo sin(theta) - Yo cos(theta)) R1 + R sum_i n_i

because on a circle x_b = Xo + R sin(a) and y_b = Yo - R cos(a). Setting both to
zero forces sum_i n_i = 0 — so F is the MOMENT factor of safety, sum(c dl) /
sum(W sin a), which is Bishop's phi = 0 answer term for term — and leaves theta
to satisfy the scalar

    h(theta) = sum_i n_i sec(a_i - theta) = 0

inside the pole-free band max(a) - pi/2 < theta < min(a) + pi/2. h runs to
+/-infinity at both ends of that band, so a sign change means a root exists and
no sign change means there is none. That is a decision, not a tolerance.

THE BAND IS THE WHOLE CLAIM

Outside that band, roots of h are COMMON — and they are not solutions. Past a
pole a slice's m_alpha has changed sign, so its base normal has reversed; the
solver's own theta bounds exclude them, and where its unbounded scipy stage lands
on one it refuses it as anomalous base tension. Circle (8, 36, R=36) below is the
pinned case: no root inside the band, roots outside it, and the shipped cascade
reaches one of them and refuses it. So the verdict is scoped — no admissible
solution — and a message claiming the surface admits no solution AT ALL is a
defect this check fails on.

WHAT IS CHECKED

* the closed form against the SHIPPED equations, on the solver's own assembled
  arrays (``residual_hook``), not a second transcription of them;
* the tutorial embankment's critical circle, where h stays negative across the
  whole band: classified as admitting no admissible solution and refused BEFORE
  the retry cascade, since no retry reaches a root that is not there;
* circle (8, 36, 36), where the same is true INSIDE the band while roots exist
  outside it: the message must say so and must not claim absolute insolubility;
* INERTNESS, measured rather than asserted: every refused circle is re-run with
  the existence test disabled, so the full retry cascade runs, and none of them
  may reach a solution the solver accepts;
* both facings: the same family mirrored, where alpha and x_b flip together, must
  get the same verdicts;
* the general case, where phi > 0 and loads reach every slice so the closed form
  does not apply: there the verdict is a measured residual floor, and it must
  quote what it measured and decline on a surface the method goes on to solve;
* the search's disclosure line on that model (N > 0), and its absence on a model
  whose every trial solves, where the output must read exactly as it always has;
* the breakdown is only printed from messages this package has READ:
  Morgenstern-Price's two are mapped and pinned by fixtures here, an unmapped
  message suppresses the taxonomy entirely, and neither may ever be reported as
  an iteration failing to converge;
* four mutations: the silent skip restored, a classifier that calls the refused
  circle an iteration failure, a verdict that claims absolute insolubility, and
  the old catch-all that folded unmapped messages into "not_converged".

Run directly:  PYTHONPATH=. python3 test/spencer_disclosure_check.py
"""

import contextlib
import io
import os
import warnings

import numpy as np
from shapely import affinity

warnings.filterwarnings('ignore')

from xslope import search as xsearch
from xslope import solve
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

#: The tutorial embankment: one phi = 0 material, a circular search, and the
#: surfaces this whole check is about. Its critical circle is refused.
MODEL = os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_simple_embankment.xlsx')

#: A model whose Spencer search answers on every admissible trial: the control
#: for "clean searches print exactly what they always printed".
CLEAN_MODEL = os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_submerged.xlsx')

#: The circle the search reports as critical on MODEL.
CRITICAL = {'Xo': 5.28, 'Yo': 33.52, 'R': 33.52}

#: A circle from the model's own circles sheet, which IS solvable.
SOLVABLE = {'Xo': 10.0, 'Yo': 40.0, 'R': 40.0}

#: The pinned band case: no root inside the admissible band, 33 roots outside it
#: (all at F = F_m = 1.235922, theta -88.3 to +110.7 deg), and the shipped
#: cascade reaches one and refuses it as 25.8x cohesive capacity. A verdict of
#: "admits no solution" would be false here; "no ADMISSIBLE solution" is the fact.
BAND_CASE = {'Xo': 8.0, 'Yo': 36.0, 'R': 36.0}

#: The general case: phi = 37 deg with reinforcement on every slice, so the
#: closed form does not apply and the verdict rests on the measured residual
#: floor (2.2e-3, and the same at either sweep density).
GENERAL_MODEL = os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_reinforce.xlsx')
GENERAL_INSOLUBLE = {'Xo': -8.0, 'Yo': 48.0, 'R': 45.0}

#: Morgenstern-Price's two failure messages, pinned to circles on MODEL that
#: produce them. Neither is an iteration failure and neither may be reported as
#: one; the first is M-P's own no-admissible-root condition.
MP_NO_CROSSING = {'Xo': 4.0, 'Yo': 30.0, 'R': 30.0}
MP_INADMISSIBLE = {'Xo': 4.0, 'Yo': 42.0, 'R': 42.0}

NUM_SLICES = 40

#: The circles the inertness leg drives, spread over centre and radius.
FAMILY = [{'Xo': xo, 'Yo': 33.52 + dr, 'R': 33.52 + dr}
          for xo in (5.28, 8.0, 10.0, 12.0) for dr in (0.0, 3.0, 6.0)]


def _circle(spec):
    c = dict(spec)
    c['Depth'] = c['Yo'] - c['R']
    return c


def _slices(model, spec, num_slices=NUM_SLICES):
    sd = model if isinstance(model, dict) else load_slope_data(model)
    ok, res = generate_slices(sd, circle=_circle(spec), num_slices=num_slices,
                              check_inputs=False)
    if not ok:
        raise AssertionError(f"could not slice {spec}: {res}")
    return res[0]


def _mirror(model):
    """The same model reflected in x, and the reflection of a circle spec.

    A genuine mirror, not a facing flag: the ground surface, the domain, the
    material polygon and the profile line are all reflected, so the slicer sees a
    right-facing slope and Spencer reaches its equations with alpha negated AND
    x_b mirrored. That pair is why the circle identity the existence test checks
    survives the reflection, which is what this makes testable.
    """
    sd = load_slope_data(model)
    axis = float(max(x for x, _ in sd['ground_surface'].coords))
    flip = lambda g: affinity.scale(g, xfact=-1.0, yfact=1.0, origin=(axis, 0.0))
    md = dict(sd,
              ground_surface=flip(sd['ground_surface']),
              domain_polygon=flip(sd['domain_polygon']),
              polygons=[dict(p, polygon=flip(p['polygon'])) for p in sd['polygons']],
              profile_lines=[dict(p, coords=[(2 * axis - x, y) for x, y in p['coords']])
                             for p in sd['profile_lines']],
              circles=[dict(c, Xo=2 * axis - c['Xo']) for c in (sd['circles'] or [])])
    return md, (lambda spec: dict(spec, Xo=2 * axis - spec['Xo']))


def _closed_form(df):
    """(n, alpha, F_m) — the decoupled problem, from the slice table alone."""
    alpha = np.radians(df['alpha'].values)
    c = df['c'].values.astype(float)
    dl = df['dl'].values.astype(float)
    W = df['w'].values.astype(float)
    F_m = float(np.sum(c * dl) / np.sum(W * np.sin(alpha)))
    return W * np.sin(alpha) - c * dl / F_m, alpha, F_m


def _h_over_band(n, alpha, count=2001):
    lo = float(np.max(alpha)) - np.pi / 2
    hi = float(np.min(alpha)) + np.pi / 2
    th = np.linspace(lo + 1e-6 * (hi - lo), hi - 1e-6 * (hi - lo), count)
    return th, np.array([float(np.sum(n / np.cos(alpha - t))) for t in th])


def _quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


#: Claims that would be false on BAND_CASE, and must never appear in a verdict.
ABSOLUTE_CLAIMS = ("no root at any interslice inclination",
                   "admits no solution on this surface")


def leg_closed_form_is_the_shipped_residual():
    """h(theta) must BE R1, on the arrays the solver assembled for itself."""
    fails = []
    for label, spec in (('refused', CRITICAL), ('solvable', SOLVABLE)):
        df = _slices(MODEL, spec)
        box = {}
        solve.spencer(df, residual_hook=lambda f: box.setdefault('f', f))
        if 'f' not in box:
            fails.append(f"{label}: spencer never handed out its residual function")
            continue
        n, alpha, F_m = _closed_form(df)
        # sum n = 0 is the moment equation at F_m: the decoupling itself.
        scale = float(np.sum(np.abs(n)))
        if abs(float(np.sum(n))) > 1e-9 * scale:
            fails.append(f"{label}: sum n = {np.sum(n):.3e} is not zero at F_m")
        th, _ = _h_over_band(n, alpha, count=9)
        worst = 0.0
        for t in th:
            h = float(np.sum(n / np.cos(alpha - t)))
            R1, _R2, _Q, _y = box['f'](F_m, t)
            worst = max(worst, abs(h - R1) / max(abs(R1), 1.0))
        if worst > 1e-11:
            fails.append(f"{label}: closed form differs from the shipped R1 by "
                         f"{worst:.2e} (relative)")
        print(f"  closed form  {label:<10} F_m={F_m:.6f}  max |h - R1|/|R1| = {worst:.1e}")
    return fails


def leg_critical_circle_has_no_admissible_root():
    """The reported critical circle of MODEL has no root in the band, and says so."""
    fails = []
    df = _slices(MODEL, CRITICAL)
    n, alpha, F_m = _closed_form(df)
    th, h = _h_over_band(n, alpha)
    if h.min() < 0 < h.max():
        fails.append("closed form finds a sign change in the band: this circle is "
                     "solvable, and the check's premise is gone")
    ok, msg = solve.spencer(df)
    if ok:
        fails.append(f"spencer solved the refused circle: FS={msg.get('FS')}")
        return fails
    kind = solve.failure_kind(msg)
    if kind != 'no_admissible_solution':
        fails.append(f"classified {kind!r}, expected 'no_admissible_solution': {msg}")
    if 'phi = 0' not in str(msg):
        fails.append(f"the message does not say which test decided it: {msg}")
    print(f"  no adm root  h in [{h.min():.3g}, {h.max():.3g}] over the band; "
          f"F_m={F_m:.4f}")
    print(f"               {msg}")
    return fails


def leg_band_case_does_not_overclaim():
    """(8, 36, 36): no root in the band, roots outside — and the verdict says so.

    This is the case that makes an absolute claim false. The shipped cascade,
    given the chance, reaches one of the out-of-band roots and refuses it on base
    tension; the early verdict must be compatible with that, not contradict it.
    """
    fails = []
    df = _slices(MODEL, BAND_CASE)
    n, alpha, F_m = _closed_form(df)
    _th, h = _h_over_band(n, alpha)
    if h.min() < 0 < h.max():
        fails.append("(8, 36, 36) now has a root inside the band; the pinned case "
                     "no longer pins anything")
        return fails
    ok, msg = solve.spencer(df)
    if ok:
        fails.append(f"spencer solved the pinned band case: FS={msg.get('FS')}")
        return fails
    text = str(msg)
    if solve.failure_kind(msg) != 'no_admissible_solution':
        fails.append(f"pinned band case classified {solve.failure_kind(msg)!r}: {msg}")
    for claim in ABSOLUTE_CLAIMS:
        if claim in text:
            fails.append(f"the verdict claims {claim!r}, which is false here: "
                         f"roots exist outside the admissible band")
    if 'outside that band' not in text:
        fails.append(f"the verdict does not report the out-of-band roots: {msg}")
    # And the cascade's own answer on the same surface, with the test disabled.
    ok2, msg2 = solve.spencer(df, existence_test=False)
    if ok2:
        fails.append(f"with the existence test off the cascade ACCEPTS a solution "
                     f"(FS={msg2.get('FS')}): the early refusal is not inert")
    else:
        print(f"  band case    {solve.failure_kind(msg)} / cascade says "
              f"{solve.failure_kind(msg2)}")
        print(f"               {msg2}")
    return fails


def leg_early_refusal_is_inert():
    """Every refused circle, re-run with the FULL cascade, stays unsolved.

    The existence test is disabled with ``existence_test=False``, so each refused
    surface gets Newton, the Bishop-seeded restarts and all 66 scipy starts —
    the same work the shipped solver did before this round. None of them may
    reach a solution the solver accepts. That, not agreement with the check's own
    copy of the closed form, is the proof the early exit refuses nothing.
    """
    fails = []
    refused = solved = 0
    for spec in FAMILY + [BAND_CASE]:
        try:
            df = _slices(MODEL, spec)
        except AssertionError:
            continue
        ok, msg = solve.spencer(df)
        if ok:
            solved += 1
            continue
        if solve.failure_kind(msg) != 'no_admissible_solution':
            continue
        refused += 1
        ok2, msg2 = solve.spencer(_slices(MODEL, spec), existence_test=False)
        if ok2:
            fails.append(f"{spec}: refused by the existence test, but the full "
                         f"cascade accepts FS={msg2.get('FS'):.4f}")
    if refused < 5:
        fails.append(f"only {refused} circles were refused; the leg proves little")
    print(f"  inertness    {refused} refused circles re-run with the cascade: "
          f"none reached an accepted solution ({solved} solved normally)")
    return fails


def leg_both_facings_agree():
    """A genuine mirror gets the same verdicts, circle for circle.

    The existence test reaches a right-facing model with alpha negated and x_b
    mirrored, which together preserve the circle identity it checks — so the test
    fires on both facings, and this is the leg that says so rather than a comment
    claiming it declines them.
    """
    fails = []
    md, mspec = _mirror(MODEL)
    fired = compared = 0
    for spec in FAMILY:
        try:
            dl_ = _slices(MODEL, spec)
            dr_ = _slices(md, mspec(spec))
        except AssertionError:
            continue
        # The mirror must be faithful before its verdict means anything: Bishop
        # is facing-blind, so the two must agree to slicing noise.
        okl, resl = solve.bishop(dl_.copy())
        okr, resr = solve.bishop(dr_.copy())
        if not (okl and okr):
            fails.append(f"{spec}: Bishop failed on one facing, so the mirror is "
                         f"not comparable")
            continue
        if abs(resl['FS'] - resr['FS']) > 1e-3 * resl['FS']:
            fails.append(f"{spec}: mirror is not faithful — Bishop {resl['FS']:.6f} "
                         f"vs {resr['FS']:.6f}")
            continue
        compared += 1
        sl, ml = solve.spencer(dl_)
        sr, mr = solve.spencer(dr_)
        if sl != sr:
            fails.append(f"{spec}: solved on one facing and not the other "
                         f"({sl} vs {sr})")
            continue
        if not sl:
            kl, kr = solve.failure_kind(ml), solve.failure_kind(mr)
            if kl != kr:
                fails.append(f"{spec}: classified {kl!r} left-facing, {kr!r} "
                             f"right-facing")
            if kl == 'no_admissible_solution':
                fired += 1
        elif abs(ml['FS'] - mr['FS']) > 1e-3 * ml['FS']:
            fails.append(f"{spec}: Spencer {ml['FS']:.6f} vs mirrored "
                         f"{mr['FS']:.6f}")
    if compared < 6:
        fails.append(f"only {compared} mirrored pairs compared; the leg proves little")
    if not fired:
        fails.append("the existence test never fired on the mirrored family, so "
                     "the right-facing path is still untested")
    print(f"  facings      {compared} mirrored pairs agree; the existence test "
          f"fired on {fired} of them")
    return fails


def leg_general_case_is_measured_not_assumed():
    """With phi > 0 and loads on every slice, the verdict rests on measurement."""
    fails = []
    df = _slices(GENERAL_MODEL, GENERAL_INSOLUBLE)
    if float(np.abs(df['phi'].values).max()) <= 0:
        fails.append("the general fixture has phi = 0 everywhere, so it is not "
                     "the general case at all")
    ok, msg = solve.spencer(df)
    if ok:
        fails.append(f"the general fixture now solves (FS={msg.get('FS'):.4f}); "
                     f"it is no longer a no-solution fixture")
        return fails
    kind = solve.failure_kind(msg)
    if kind != 'no_admissible_solution':
        fails.append(f"general case classified {kind!r}, expected "
                     f"'no_admissible_solution': {msg}")
    elif 'times the' not in str(msg):
        fails.append(f"the message does not quote the floor it measured: {msg}")
    else:
        print(f"  general      {msg}")

    sd = load_slope_data(GENERAL_MODEL)
    solvable = dict(sd['circles'][0])
    df_ok = _slices(GENERAL_MODEL, {'Xo': solvable['Xo'], 'Yo': solvable['Yo'],
                                    'R': solvable['R']})
    ok2, msg2 = solve.spencer(df_ok)
    if not ok2:
        fails.append(f"the model's own starting circle no longer solves: {msg2}")
    else:
        print(f"  general      the model's starting circle still solves: "
              f"FS={msg2['FS']:.4f}")
    return fails


def leg_moment_fs_is_recorded_and_exact():
    """The moment factor of safety the ranking uses is Bishop's, to rounding."""
    fails = []
    for label, spec in (('refused', CRITICAL), ('solvable', SOLVABLE)):
        df = _slices(MODEL, spec)
        solve.spencer(df)
        got = df.attrs.get('moment_fs')
        ok_b, res_b = solve.bishop(df.copy())
        if not ok_b:
            fails.append(f"{label}: Bishop could not solve it, so there is nothing "
                         f"to compare the recorded moment FS with")
            continue
        if got is None:
            fails.append(f"{label}: spencer recorded no moment_fs for the search to read")
            continue
        rel = abs(got - res_b['FS']) / res_b['FS']
        if rel > 1e-9:
            fails.append(f"{label}: recorded moment FS {got:.10f} differs from "
                         f"Bishop's {res_b['FS']:.10f} ({rel:.1e})")
        print(f"  moment FS    {label:<10} recorded {got:.10f}  Bishop {res_b['FS']:.10f}")
    return fails


def _search(model, method='spencer'):
    sd = load_slope_data(model)
    out = {}
    (fs_cache, converged, _p, _c), text = _quiet(
        xsearch.circular_search, sd, method, num_slices=NUM_SLICES,
        diagnostic=False, unsolved_out=out)
    return fs_cache, converged, out, text


def _breakdown_adds_up(out):
    return (out['no_admissible_solution'] + out['not_converged']
            + out['inadmissible'] + out['unclassified']) == out['unsolved']


def leg_search_discloses():
    """The search says how many trial surfaces its method could not answer on."""
    fails = []
    fs_cache, converged, out, text = _search(MODEL)
    if not out.get('unsolved'):
        fails.append("the Spencer search on the tutorial embankment reported no "
                     "unsolved trials, so nothing was disclosed")
        return fails
    if out['unsolved'] > out['attempted']:
        fails.append(f"{out['unsolved']} unsolved of {out['attempted']} attempted")
    if not _breakdown_adds_up(out):
        fails.append(f"the breakdown does not add up: {out}")
    if not out['no_admissible_solution']:
        fails.append("none of the unsolved trials was classified as admitting no "
                     "admissible solution, though that is this model's class")
    if out['reported_moment_fs'] is None:
        fails.append("no moment measure on the reported minimum, so nothing could "
                     "be ranked against it")
    elif not out['lower_by_moment']:
        fails.append("no unsolved trial ranks below the reported minimum, though "
                     "the refused region contains it")
    line = [l for l in text.splitlines() if 'could not solve' in l]
    if len(line) != 1:
        fails.append(f"expected exactly one disclosure line, found {len(line)}")
    else:
        for want in ('Spencer could not solve', 'trial surfaces',
                     'rank lower than the reported minimum by the moment measure'):
            if want not in line[0]:
                fails.append(f"the disclosure line does not say {want!r}: {line[0]}")
        for claim in ABSOLUTE_CLAIMS:
            if claim in line[0]:
                fails.append(f"the disclosure line claims {claim!r}: {line[0]}")
        print(f"  disclosure   {line[0].strip()}")
    if not converged:
        fails.append("the search no longer converges on this model")
    print(f"  reported     FS={fs_cache[0]['FS']:.6f} converged={converged}")
    return fails


def leg_breakdown_only_from_read_messages():
    """A taxonomy is printed only where the messages behind it have been read.

    Morgenstern-Price fails on this model two ways, neither of them an iteration
    failure: no F_f/F_m crossing anywhere in its lambda range (its own
    no-admissible-root condition) and its admissibility refusal. Both are pinned
    to circles here, so the mapping is measured rather than assumed. An unmapped
    message must suppress the breakdown entirely rather than fall into a class.
    """
    fails = []
    for label, spec, want in (('no-crossing', MP_NO_CROSSING, 'no_admissible_solution'),
                              ('inadmissible', MP_INADMISSIBLE, 'inadmissible')):
        df = _slices(MODEL, spec)
        ok, msg = solve.mprice(df)
        if ok:
            fails.append(f"M-P fixture {label} now solves (FS={msg.get('FS'):.4f}); "
                         f"the mapping is no longer pinned by it")
            continue
        got = solve.failure_kind(msg)
        if got != want:
            fails.append(f"M-P {label} classified {got!r}, expected {want!r}: {msg}")
        else:
            print(f"  M-P mapping  {label:<12} -> {got}")

    fs_cache, _conv, out, text = _search(MODEL, 'mprice')
    if not out['unsolved']:
        fails.append("the M-P search reported no unsolved trials, so the line "
                     "this leg is about was never printed")
        return fails
    if not _breakdown_adds_up(out):
        fails.append(f"the M-P breakdown does not add up: {out}")
    if out['not_converged']:
        fails.append(f"{out['not_converged']} M-P trials reported as iteration "
                     f"failures, though M-P's two messages are a no-root condition "
                     f"and an admissibility refusal")
    line = [l for l in text.splitlines() if 'could not solve' in l]
    if len(line) != 1:
        fails.append(f"expected one M-P disclosure line, found {len(line)}")
    elif 'failed to converge' in line[0]:
        fails.append(f"the M-P line reports convergence failures: {line[0]}")
    else:
        print(f"  M-P line     {line[0].strip()}")

    # An unmapped message must leave the taxonomy off the line entirely.
    tally = xsearch.UnsolvedTrials('bishop')
    tally.record(1.0, 2.0, 3.0, False, None, "Some solver said something new.")
    if tally.unclassified != 1:
        fails.append("an unmapped message was folded into a class instead of "
                     "being counted as unclassified")
    if '(' in tally.sentence():
        fails.append(f"an unmapped message still produced a breakdown: "
                     f"{tally.sentence()}")
    else:
        print(f"  unmapped     {tally.sentence().strip()}")
    return fails


def leg_clean_search_is_silent():
    """A search whose every trial solved prints what it always printed."""
    fails = []
    fs_cache, converged, out, text = _search(CLEAN_MODEL)
    if out.get('unsolved'):
        fails.append(f"{CLEAN_MODEL} is no longer a clean control: {out}")
    noisy = [l for l in text.splitlines()
             if 'could not solve' in l or 'unsolved' in l]
    if noisy:
        fails.append(f"a clean search printed the disclosure anyway: {noisy}")
    if not converged:
        fails.append("the control search no longer converges")
    print(f"  clean        {os.path.basename(CLEAN_MODEL)}: "
          f"{out['attempted']} trials, none unsolved, FS={fs_cache[0]['FS']:.6f}")
    return fails


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


def leg_mutations():
    """Every defect this check exists for must make it fail."""
    fails = []

    # 1. The silent skip: the search stops counting what it could not solve.
    original = xsearch.UnsolvedTrials.record
    _mutation("silent skip",
              lambda: setattr(xsearch.UnsolvedTrials, 'record',
                              lambda self, *a, **k: None),
              lambda: setattr(xsearch.UnsolvedTrials, 'record', original),
              leg_search_discloses, fails)

    # 2. The classifier calls the refused circle an iteration failure.
    original_kind = solve.failure_kind
    _mutation("misclassification",
              lambda: setattr(solve, 'failure_kind', lambda msg: 'not_converged'),
              lambda: setattr(solve, 'failure_kind', original_kind),
              leg_critical_circle_has_no_admissible_root, fails)

    # 3. The verdict claims the surface admits no solution AT ALL — false on
    #    BAND_CASE, where roots exist outside the admissible band.
    original_spencer = solve.spencer

    def overclaiming(df, *a, **k):
        ok, msg = original_spencer(df, *a, **k)
        if not ok and solve.failure_kind(msg) == 'no_admissible_solution':
            return ok, ("Spencer's parallel interslice force assumption admits no "
                        "solution on this surface: no root at any interslice "
                        "inclination.")
        return ok, msg

    _mutation("absolute claim",
              lambda: setattr(solve, 'spencer', overclaiming),
              lambda: setattr(solve, 'spencer', original_spencer),
              leg_band_case_does_not_overclaim, fails)

    # 4. The old catch-all: every unmapped message becomes a convergence failure.
    _mutation("catch-all class",
              lambda: setattr(solve, 'failure_kind',
                              lambda msg: original_kind(msg) or 'not_converged'),
              lambda: setattr(solve, 'failure_kind', original_kind),
              leg_breakdown_only_from_read_messages, fails)
    return fails


LEGS = [
    ("the closed form is the shipped residual", leg_closed_form_is_the_shipped_residual),
    ("the critical circle has no admissible root", leg_critical_circle_has_no_admissible_root),
    ("the pinned band case does not overclaim", leg_band_case_does_not_overclaim),
    ("the early refusal is inert under the full cascade", leg_early_refusal_is_inert),
    ("both facings agree", leg_both_facings_agree),
    ("the general case is measured", leg_general_case_is_measured_not_assumed),
    ("the moment factor of safety is recorded", leg_moment_fs_is_recorded_and_exact),
    ("the search discloses what it could not solve", leg_search_discloses),
    ("the breakdown comes only from read messages", leg_breakdown_only_from_read_messages),
    ("a clean search stays silent", leg_clean_search_is_silent),
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
    print("\nPASS: surfaces with no admissible solution are classified as that and "
          "no more, refused without iterating, and disclosed by the search that "
          "met them.")


if __name__ == "__main__":
    main()
