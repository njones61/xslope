"""The force-equilibrium closure has several roots; only one is a body of soil.

Corps of Engineers and Lowe-Karafiath both solve the same residual: march the
slices left to right at a trial factor of safety and drive the side force left
over past the last slice, Z[n+1], to zero. That march divides by a per-slice
determinant which reduces to

    det_i = cos(alpha_i - theta_{i+1}) + tan(phi_i) * sin(alpha_i - theta_{i+1}) / F

and vanishes at F = -tan(phi_i) * tan(alpha_i - theta_{i+1}). Every slice whose
base is inclined more than ninety degrees away from its own side force puts a
pole into the residual, and each pole separates one root from the next. On the
two surfaces pinned here the residual carries fourteen and nineteen roots between
F = 0.05 and F = 10.

The solver used to take whichever root a secant from F = 1.5 happened to reach.
On the Corps critical surface of the pile model that was F = 1.061, a root whose
side forces run to 2.6 times the weight of the sliding mass, while Spencer,
Bishop, Morgenstern-Price and Lowe all read about 6.55 on the identical slices;
on the right-facing reinforced model it was F = 0.696 against about 2.14. Both
numbers were published. An interslice-tension cap had been excluding some of
those roots and not others, which is why the values published before that cap
was lifted (1.369 and 1.091) were of the same class rather than better.

WHAT DECIDES WHICH ROOT

Three measures, each calibrated against every Corps and Lowe solution the
shipped corpus produces on its own surfaces:

  * the base factor above stays clear of zero (solve.FE_MIN_BASE_FACTOR).
    Answers that agree with the moment methods keep it above 0.2 on every slice,
    while the roots wedged between two poles sit at 0.003 to 0.022 — the
    residual sweeps the whole real line between consecutive poles, so a root
    always exists there, and it is always this spurious branch.
  * the largest side force stays within the weight of the sliding mass
    (solve.FE_MAX_Z_OVER_W). Accepted corpus solutions run 0.11 to 0.61 of it.
  * base tension does not saturate (solve.MAX_BASE_TENSION_EXTENT), the same
    extent test applied to the reported solution.

Where more than one root passes, the one nearest the moment-equilibrium answer
on the same slices is reported — Bishop on a circular surface, Spencer on a
non-circular one. Where none passes, the method says so; no other method's
factor of safety is ever substituted.

A per-boundary form of the second measure was tried and is deliberately absent:
comparing each |Z_j| against the soil column standing above that boundary,
accepted corpus solutions reach 11.5 times the column on a full-height boundary
and 58 times when the short boundaries at crest and toe are counted, so no bar
on it separates a right answer from a wrong one.

WHAT IS CHECKED

* the two surfaces the defect was found on carry fourteen and nineteen roots,
  and exactly one of each survives the three measures;
* the reported root is not the low branch, and its side forces stay within the
  mass;
* the tie-break, on a closure built to have two admissible roots either side of
  the moment answer: the nearer is reported, not the smaller;
* the roots passed over are named in results['warnings'] with the measure each
  failed;
* on each model's own Corps critical surface, Corps reads within 5% of Spencer
  on the identical slices — the agreement every other model in the corpus shows;
* the two no-root failures classify as no_admissible_solution and the base
  tension refusal as inadmissible, so a search reports them by class;
* a mutation: the unbracketed secant restored, which puts the low branch back.

Run directly:  PYTHONPATH=. python3 test/force_closure_root_check.py
"""

import os
import re
import warnings

import numpy as np

warnings.filterwarnings('ignore')

from scipy.optimize import newton

from xslope import solve
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

NUM_SLICES = 40


def _circle(xo, yo, r):
    return {'Xo': xo, 'Yo': yo, 'R': r, 'Depth': yo - r}


#: The two surfaces the unbracketed secant reported a wrong root on. These are
#: the circles its own search settled on, kept as geometry so the solver is
#: checked without re-running a search.
PATHOLOGICAL = {
    'piles': dict(
        model=os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_piles.xlsx'),
        circle=_circle(35.0, 43.75, 46.674805),
        secant_fs=1.0609,          # what the secant reported here
        roots=13,                  # bracketed roots between FE_FS_MIN and FE_FS_MAX
    ),
    'reinforce_rface': dict(
        model=os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_reinforce_rface.xlsx'),
        circle=_circle(86.235687, 33.280951, 37.503698),
        secant_fs=0.6958,
        roots=16,
    ),
}

#: Each model's Corps critical surface under the bracketed root finder. Corps is
#: a force-equilibrium method and reads a few percent above Spencer; what is
#: pinned is that it reads the same branch, not the same number.
CRITICAL = {
    'piles': dict(
        model=os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_piles.xlsx'),
        circle=_circle(-7.633636475, 51.459037781, 48.145221978),
    ),
    'reinforce_rface': dict(
        model=os.path.join(REPO, 'docs', 'lem', 'files', 'xslope_reinforce_rface.xlsx'),
        circle=_circle(107.303637695, 48.281365967, 48.830256438),
    ),
}

#: How far a force-equilibrium factor of safety may stand from Spencer's on the
#: same slices before it is a different branch rather than a different side-force
#: assumption. Corps runs above Spencer by construction — the ground-parallel
#: convention is the least conservative of the family — and across the corpus
#: that gap is a few percent.
SPENCER_TOLERANCE = 0.05

ROOT_NOTE = 'force closure has other roots'


def _slices(model, circle):
    sd = load_slope_data(model)
    ok, res = generate_slices(sd, circle=dict(circle), num_slices=NUM_SLICES)
    if not ok:
        raise AssertionError(f"could not slice {os.path.basename(model)}: {res}")
    return res[0]


def _corps_closure(df):
    """The closure, its arrays and its poles for the Corps convention on `df`.

    Rebuilt the way `solve.corps` builds it so the check reads the same residual
    the solver roots, without reaching into the solver's internals twice.
    """
    n = len(df)
    x_l, y_lt = df['x_l'].values, df['y_lt'].values
    x_r, y_rt = df['x_r'].values, df['y_rt'].values
    slope_top = (y_rt - y_lt) / (x_r - x_l)
    theta_list = np.zeros(n + 1)
    for j in range(n + 1):
        s = slope_top[0] if j == 0 else (slope_top[-1] if j == n
                                         else 0.5 * (slope_top[j - 1] + slope_top[j]))
        theta_list[j] = np.degrees(np.arctan(s))
    right_facing = bool(df.attrs.get('right_facing',
                                     df['y_lb'].iat[0] > df['y_rb'].iat[-1]))
    if right_facing:
        theta_list = -theta_list
    return solve._force_closure(df.copy(), theta_list, right_facing)


def _brackets(residual, poles):
    """Every root the solver's own scan brackets, in order."""
    from scipy.optimize import brentq
    edge = 1e-7 * (solve.FE_FS_MAX - solve.FE_FS_MIN)
    inside = poles[(poles > solve.FE_FS_MIN + edge) & (poles < solve.FE_FS_MAX - edge)]
    grid = np.unique(np.concatenate([
        np.linspace(solve.FE_FS_MIN, solve.FE_FS_MAX, solve.FE_FS_SCAN),
        inside - edge, inside + edge]))
    at_pole = np.isin(grid, inside - edge)
    vals = np.array([residual(float(F)) for F in grid])
    out = []
    for i in range(grid.size - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b < 0.0 and not at_pole[i]:
            try:
                out.append(float(brentq(residual, grid[i], grid[i + 1],
                                        xtol=1e-9, maxiter=100)))
            except Exception:
                pass
    return out


def _classify(residual, Z, N, det, poles, W_sum):
    """(admissible, rejected) over the bracketed roots, by the shipped measures."""
    admissible, rejected = [], []
    for r in _brackets(residual, poles):
        residual(r)
        if float(np.abs(det).min()) < solve.FE_MIN_BASE_FACTOR:
            rejected.append(r)
        elif float(np.abs(Z).max()) > solve.FE_MAX_Z_OVER_W * W_sum:
            rejected.append(r)
        elif float(np.mean(N < 0)) > solve.MAX_BASE_TENSION_EXTENT:
            rejected.append(r)
        else:
            admissible.append(r)
    return admissible, rejected


def leg_the_closure_has_several_roots():
    """Each pathological surface carries many roots and few admissible ones."""
    fails = []
    for name, spec in PATHOLOGICAL.items():
        df = _slices(spec['model'], spec['circle'])
        residual, N, Z, det, poles = _corps_closure(df)
        W_sum = float(np.abs(df['w'].values).sum())
        admissible, rejected = _classify(residual, Z, N, det, poles, W_sum)
        total = len(admissible) + len(rejected)
        if total < spec['roots']:
            fails.append(f"{name}: the closure brackets {total} root(s), fewer than "
                         f"the {spec['roots']} this surface carried, so it no longer "
                         f"exercises the selection")
        if len(admissible) != 1:
            fails.append(f"{name}: {len(admissible)} admissible root(s), expected 1 "
                         f"({[round(r, 4) for r in admissible]})")
        else:
            print(f"  {name:16s} {total} roots, 1 admissible "
                  f"{[round(r, 4) for r in admissible]}")
    return fails


def leg_the_low_branch_is_not_reported():
    """The root the secant used to return is refused, on its own measures."""
    fails = []
    for name, spec in PATHOLOGICAL.items():
        df = _slices(spec['model'], spec['circle'])
        ok, res = solve.corps(df)
        if not ok:
            fails.append(f"{name}: corps returned no answer ({res})")
            continue
        if abs(res['FS'] - spec['secant_fs']) < 1e-2:
            fails.append(f"{name}: corps still reports the secant's branch "
                         f"{res['FS']:.4f}")
        Z = np.asarray(df['z'].values, dtype=float)
        W_sum = float(np.abs(df['w'].values).sum())
        if float(np.abs(Z).max()) > solve.FE_MAX_Z_OVER_W * W_sum:
            fails.append(f"{name}: the reported root carries side forces of "
                         f"{np.abs(Z).max() / W_sum:.1f} times the weight of the mass")
        else:
            print(f"  {name:16s} corps {res['FS']:.4f} (was {spec['secant_fs']}), "
                  f"max|Z| = {np.abs(Z).max() / W_sum:.2f} of the mass")
    return fails


def leg_the_moment_answer_breaks_the_tie():
    """Two admissible roots: the reported one is the nearer, not the smaller.

    Every surface in the corpus leaves exactly one root standing once the three
    measures are applied, so the tie-break is exercised on a closure built for
    the purpose: two roots either side of the moment answer, the nearer one the
    larger, with the arrays it fills held at values every measure accepts.
    """
    spec = CRITICAL['reinforce_rface']
    df = _slices(spec['model'], spec['circle'])
    ok_b, res_b = solve.bishop(df.copy())
    if not ok_b:
        return ["reinforce_rface: Bishop gives no moment reference here"]
    ref = res_b['FS']
    lo, hi = ref - 0.28, ref + 0.12          # the nearer root is the larger one
    n = len(df)
    W_sum = float(np.abs(df['w'].values).sum())
    Z = np.zeros(n + 1)
    N = np.zeros(n)
    det = np.zeros(n)

    def residual(F):
        Z[:] = 0.01 * W_sum
        Z[0] = 0.0
        N[:] = 1.0
        det[:] = 1.0
        return (F - lo) * (F - hi)

    got, warns = solve._force_closure_root(residual, Z, N, det,
                                           np.array([]), df.copy())
    fails = []
    if got is None:
        return [f"the synthetic closure was refused: {warns}"]
    if abs(got - hi) > 1e-6:
        fails.append(f"the tie-break returned {got:.4f}, not the root nearest "
                     f"Bishop's {ref:.4f} ({hi:.4f}); the smaller root is {lo:.4f}")
    if not any(f"{lo:.3f}" in w for w in warns):
        fails.append(f"the root passed over ({lo:.3f}) is not named in {warns}")
    if not fails:
        print(f"  synthetic       roots {lo:.4f} / {hi:.4f}, Bishop {ref:.4f} "
              f"-> {got:.4f}")
    return fails


def leg_discarded_roots_are_reported():
    """Every root passed over is named, with the measure it failed."""
    fails = []
    for name, spec in PATHOLOGICAL.items():
        df = _slices(spec['model'], spec['circle'])
        ok, res = solve.corps(df)
        if not ok:
            fails.append(f"{name}: corps returned no answer ({res})")
            continue
        note = next((w for w in res.get('warnings') or []
                     if w.startswith(ROOT_NOTE)), None)
        if note is None:
            fails.append(f"{name}: no warning names the roots that were passed over")
            continue
        named = re.findall(r'FS=([0-9.]+) \(([^)]+)\)', note)
        if not named:
            fails.append(f"{name}: the warning names no root: {note!r}")
            continue
        if any(why.strip() == '' for _, why in named):
            fails.append(f"{name}: a root is named without a reason: {note!r}")
        if any(abs(float(fs) - res['FS']) < 1e-9 for fs, _ in named):
            fails.append(f"{name}: the reported root appears in its own list of "
                         f"discarded roots")
        print(f"  {name:16s} {len(named)} root(s) named: "
              + ", ".join(f"{fs} ({why})" for fs, why in named[:3])
              + (" ..." if len(named) > 3 else ""))
    return fails


def leg_corps_tracks_spencer_on_its_critical_surface():
    """On each model's own Corps critical surface, Corps reads Spencer's branch."""
    fails = []
    for name, spec in CRITICAL.items():
        df = _slices(spec['model'], spec['circle'])
        ok_c, res_c = solve.corps(df.copy())
        ok_s, res_s = solve.spencer(df.copy())
        if not (ok_c and ok_s):
            fails.append(f"{name}: corps {'ok' if ok_c else 'failed'}, "
                         f"spencer {'ok' if ok_s else 'failed'}")
            continue
        gap = res_c['FS'] / res_s['FS'] - 1.0
        if abs(gap) > SPENCER_TOLERANCE:
            fails.append(f"{name}: corps {res_c['FS']:.4f} is {100 * gap:+.1f}% from "
                         f"Spencer's {res_s['FS']:.4f} on the same slices")
        else:
            print(f"  {name:16s} corps {res_c['FS']:.4f} vs spencer "
                  f"{res_s['FS']:.4f} ({100 * gap:+.1f}%)")
    return fails


def leg_failures_are_classified():
    """A search must be able to say which kind of failure it skipped."""
    fails = []
    cases = (
        (solve.FORCE_EQ_NO_ROOT + "0.05 to 10", 'no_admissible_solution'),
        (solve.FORCE_EQ_NO_ADMISSIBLE_ROOT + "; rejected FS=0.500 (base factor 0.001)",
         'no_admissible_solution'),
        (solve.FORCE_EQ_INADMISSIBLE + "(80% of base normals in tension)",
         'inadmissible'),
    )
    for message, expected in cases:
        got = solve.failure_kind(message)
        if got != expected:
            fails.append(f"{message[:48]!r} classified {got!r}, expected {expected!r}")
    if not fails:
        print(f"  {len(cases)} force-equilibrium failure messages classified")
    return fails


def _mutation(label, apply, restore, leg, fails):
    apply()
    try:
        caught = leg()
    finally:
        restore()
    if not caught:
        fails.append(f"{label}: the leg passed with the defect in place")
    else:
        print(f"  mutation     {label} -> caught ({len(caught)} failure(s))")


def leg_mutations():
    """The unbracketed secant must put the low branch back."""
    fails = []
    original = solve._force_closure_root

    def secant(residual, Z, N, det, poles, slice_df):
        """The root finder as it stood: one secant from F = 1.5, no bracket."""
        try:
            return float(newton(residual, 1.5, tol=1e-6, maxiter=50)), []
        except Exception as e:                              # pragma: no cover
            return None, str(e)

    _mutation("the unbracketed secant",
              lambda: setattr(solve, '_force_closure_root', secant),
              lambda: setattr(solve, '_force_closure_root', original),
              leg_the_low_branch_is_not_reported, fails)
    return fails


LEGS = [
    ("the closure carries several roots", leg_the_closure_has_several_roots),
    ("the low branch is not reported", leg_the_low_branch_is_not_reported),
    ("the moment answer breaks the tie", leg_the_moment_answer_breaks_the_tie),
    ("the discarded roots are reported", leg_discarded_roots_are_reported),
    ("corps reads Spencer's branch on its critical surface",
     leg_corps_tracks_spencer_on_its_critical_surface),
    ("failures are classified", leg_failures_are_classified),
    ("mutations", leg_mutations),
]


def run():
    failures = []
    for label, fn in LEGS:
        print(f"[{label}]")
        try:
            failures.extend(fn())
        except Exception as e:
            failures.append(f"{label}: raised {type(e).__name__}: {e}")
    return failures


if __name__ == '__main__':
    import sys
    fails = run()
    if fails:
        print("\nFAILURES:")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("\nforce-closure root selection: all legs pass")
