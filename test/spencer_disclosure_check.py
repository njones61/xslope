"""Spencer's insoluble surfaces: classified as such, and disclosed by the search.

Spencer's method assumes every interslice force acts at the SAME inclination
theta, which turns equilibrium into two equations in (F, theta). On some trial
surfaces that system HAS NO ROOT. No starting guess and no iteration count
reaches an answer that is not there, and until this check's subject shipped, two
things followed from that.

First, the failure read as the solver's: "did not converge within the maximum
number of iterations" says the iteration ran out, when what happened is that the
method admits no answer on this surface. Second — the one that reaches a printed
factor of safety — the SEARCH scored the unsolvable trial exactly as it scores a
circle that misses the model: fs_fail, dropped, never mentioned. So a search
could report its minimum as converged while surfaces LOWER than it went
unanswered, and nothing in the output said so.

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

WHAT IS CHECKED

* the closed form against the SHIPPED equations, on the solver's own assembled
  arrays (``residual_hook``), not a second transcription of them;
* the tutorial embankment's critical circle (Xo=5.28, Yo=33.52, R=33.52), where
  h stays negative across the whole band: classified as admitting no solution,
  and rejected BEFORE the retry cascade, since a rootless system cannot be
  solved by retrying;
* inertness of that early exit — over a family of circles on the same model, the
  closed form's verdict and the solver's outcome agree case for case, so no
  surface that could have been solved is refused;
* the search's disclosure line on that model (N > 0), and its absence on a model
  whose every trial solves, where the output must read exactly as it always has;
* two mutations: the silent skip restored, and a classifier that calls the
  insoluble circle an iteration failure. Both must make this check fail.

Run directly:  PYTHONPATH=. python3 test/spencer_disclosure_check.py
"""

import contextlib
import io
import os
import warnings

import numpy as np

warnings.filterwarnings('ignore')

from xslope import search as xsearch
from xslope import solve
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

#: The tutorial embankment: one phi = 0 material, a circular search, and the
#: surfaces this whole check is about. Its critical circle is insoluble.
MODEL = os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_simple_embankment.xlsx')

#: A model whose Spencer search answers on every admissible trial: the control
#: for "clean searches print exactly what they always printed".
CLEAN_MODEL = os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_submerged.xlsx')

#: The circle the search reports as critical on MODEL, and the one the
#: diagnostic pinned the no-solution verdict to.
CRITICAL = {'Xo': 5.28, 'Yo': 33.52, 'R': 33.52}

#: A circle from the model's own circles sheet, which IS solvable.
SOLVABLE = {'Xo': 10.0, 'Yo': 40.0, 'R': 40.0}

NUM_SLICES = 40


def _circle(spec):
    c = dict(spec)
    c['Depth'] = c['Yo'] - c['R']
    return c


def _slices(model, spec, num_slices=NUM_SLICES):
    sd = load_slope_data(model)
    ok, res = generate_slices(sd, circle=_circle(spec), num_slices=num_slices,
                              check_inputs=False)
    if not ok:
        raise AssertionError(f"could not slice {spec}: {res}")
    return res[0]


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


def leg_closed_form_is_the_shipped_residual():
    """h(theta) must BE R1, on the arrays the solver assembled for itself."""
    fails = []
    for label, spec in (('insoluble', CRITICAL), ('solvable', SOLVABLE)):
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


def leg_critical_circle_admits_no_solution():
    """The reported critical circle of MODEL has no root, and is said to have none."""
    fails = []
    df = _slices(MODEL, CRITICAL)
    n, alpha, F_m = _closed_form(df)
    th, h = _h_over_band(n, alpha)
    if h.min() < 0 < h.max():
        fails.append("closed form finds a sign change: this circle is solvable, "
                     "and the check's premise is gone")
    ok, msg = solve.spencer(df)
    if ok:
        fails.append(f"spencer solved the insoluble circle: FS={msg.get('FS')}")
        return fails
    kind = solve.spencer_failure_kind(msg)
    if kind != 'no_solution':
        fails.append(f"classified {kind!r}, expected 'no_solution': {msg}")
    if 'phi = 0' not in str(msg):
        fails.append(f"the message does not say which test decided it: {msg}")
    print(f"  no solution  h in [{h.min():.3g}, {h.max():.3g}] over the band; "
          f"F_m={F_m:.4f}")
    print(f"               {msg}")
    return fails


def leg_moment_fs_is_recorded_and_exact():
    """The moment factor of safety the ranking uses is Bishop's, to rounding."""
    fails = []
    for label, spec in (('insoluble', CRITICAL), ('solvable', SOLVABLE)):
        df = _slices(MODEL, spec)
        solve.spencer(df)
        got = df.attrs.get('moment_fs')
        ok_b, res_b = solve.bishop(df)
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


def leg_early_exit_is_inert():
    """Over a family of circles: closed-form verdict and solver outcome agree.

    The early exit refuses a surface before the retry cascade runs. That is only
    safe if the verdict never lands on a surface the cascade would have solved,
    so both are run over a spread of circles and compared case for case.
    """
    fails = []
    tried = solved = refused = 0
    for xo in (5.28, 8.0, 10.0, 12.0):
        for r_extra in (0.0, 3.0, 6.0):
            spec = {'Xo': xo, 'Yo': 33.52 + r_extra, 'R': 33.52 + r_extra}
            try:
                df = _slices(MODEL, spec)
            except AssertionError:
                continue
            n, alpha, _F = _closed_form(df)
            _th, h = _h_over_band(n, alpha)
            has_root = bool(h.min() < 0 < h.max())
            ok, msg = solve.spencer(df)
            tried += 1
            solved += bool(ok)
            refused += (not ok) and solve.spencer_failure_kind(msg) == 'no_solution'
            if has_root and not ok and solve.spencer_failure_kind(msg) == 'no_solution':
                fails.append(f"{spec}: h changes sign (a root exists) yet the "
                             f"surface was refused as insoluble")
            if (not has_root) and ok:
                fails.append(f"{spec}: no root in the band, yet spencer returned "
                             f"FS={msg.get('FS'):.4f}")
    if tried < 6:
        fails.append(f"only {tried} circles in the family sliced; the leg proves little")
    print(f"  inertness    {tried} circles: {solved} solved, {refused} refused as "
          f"insoluble, 0 disagreements with the closed form")
    return fails


def _search(model, method='spencer'):
    sd = load_slope_data(model)
    out = {}
    (fs_cache, converged, _p, _c), text = _quiet(
        xsearch.circular_search, sd, method, num_slices=NUM_SLICES,
        diagnostic=False, unsolved_out=out)
    return fs_cache, converged, out, text


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
    if out['no_solution'] + out['not_converged'] + out['inadmissible'] != out['unsolved']:
        fails.append(f"the breakdown does not add up: {out}")
    if not out['no_solution']:
        fails.append("none of the unsolved trials was classified insoluble, though "
                     "this model's failures are the phi = 0 no-root class")
    if out['reported_moment_fs'] is None:
        fails.append("no moment measure on the reported minimum, so nothing could "
                     "be ranked against it")
    elif not out['lower_by_moment']:
        fails.append("no unsolved trial ranks below the reported minimum, though "
                     "the insoluble region contains it")
    line = [l for l in text.splitlines() if 'could not solve' in l]
    if len(line) != 1:
        fails.append(f"expected exactly one disclosure line, found {len(line)}")
    else:
        for want in ('Spencer could not solve', 'trial surfaces',
                     'rank lower than the reported minimum by the moment measure'):
            if want not in line[0]:
                fails.append(f"the disclosure line does not say {want!r}: {line[0]}")
        print(f"  disclosure   {line[0].strip()}")
    if not converged:
        fails.append("the search no longer converges on this model")
    print(f"  reported     FS={fs_cache[0]['FS']:.6f} converged={converged}")
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


def leg_mutations():
    """Both defects this check exists for must make it fail."""
    fails = []

    # 1. The silent skip: the search stops counting what it could not solve.
    original = xsearch.UnsolvedTrials.record
    xsearch.UnsolvedTrials.record = lambda self, *a, **k: None
    try:
        caught = leg_search_discloses()
    finally:
        xsearch.UnsolvedTrials.record = original
    if not caught:
        fails.append("silent-skip mutation: the disclosure leg passed with the "
                     "search counting nothing")
    else:
        print(f"  mutation     silent skip -> caught ({len(caught)} failure(s))")

    # 2. The classifier calls the insoluble circle an iteration failure.
    original_kind = solve.spencer_failure_kind
    solve.spencer_failure_kind = lambda msg: 'not_converged'
    try:
        caught = leg_critical_circle_admits_no_solution()
    finally:
        solve.spencer_failure_kind = original_kind
    if not caught:
        fails.append("classifier mutation: the insoluble circle passed while being "
                     "reported as an iteration failure")
    else:
        print(f"  mutation     misclassification -> caught ({len(caught)} failure(s))")
    return fails


LEGS = [
    ("the closed form is the shipped residual", leg_closed_form_is_the_shipped_residual),
    ("the critical circle admits no solution", leg_critical_circle_admits_no_solution),
    ("the moment factor of safety is recorded", leg_moment_fs_is_recorded_and_exact),
    ("the early exit is inert", leg_early_exit_is_inert),
    ("the search discloses what it could not solve", leg_search_discloses),
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
    print("\nPASS: Spencer's insoluble surfaces are classified as insoluble, "
          "refused without iterating, and disclosed by the search that met them.")


if __name__ == "__main__":
    main()
