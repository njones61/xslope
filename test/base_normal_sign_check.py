"""The simplified methods had no test on the sign of the base normal.

Bishop's and Janbu's methods take the base normal from vertical equilibrium of
each slice with the interslice shear neglected, which puts the factor of safety
inside the answer:

    N' = [ W - (c dl - u dl tan(phi)) sin(alpha) / F ] / m_alpha,
    m_alpha = cos(alpha) + sin(alpha) tan(phi) / F

m_alpha vanishes at tan(alpha) = -F / tan(phi). Past that base inclination it is
negative, N' has reversed, and the iteration goes on converging: it lands on the
far side of a singularity and reports a number. Nothing tested for it.

WHAT IT COST

`noncircular_search` on `xslope_noncircular.xlsx`, run with the base-angle cap
opened to 75 degrees, reaches a surface whose entry ramp is a single straight
segment at alpha = -59.96 degrees in a phi = 33 material. Four of its forty-one
slices sit past the singularity: m_alpha falls to -0.2507 and the least base
normal is -963.6. Janbu reports 0.8152 there. On the identical slices Spencer
reads 4.3825, Morgenstern-Price 3.4116, Lowe & Karafiath 1.6107 and Corps of
Engineers 1.4036 — Janbu is a factor of five below the complete-equilibrium
answer, not a few percent below it.

The default 65-degree `max_base_angle` cap kept the search off that basin, but
only incidentally: the winning surface's own steepest base is 59.96 degrees, well
inside the cap, and a 66-degree cap is enough to reach it. The cap is a geometric
admissibility rule and stays one; it is not this test and was never meant to be.

WHAT DECIDES IT NOW

`solve._base_normal_admissibility`, applied by the Ordinary, Bishop and Janbu
methods to their own converged answer:

  * m_alpha stays clear of zero on every slice, at `solve.FE_MIN_BASE_FACTOR` —
    the same bar the force-equilibrium closure keeps its own base factor above,
    and the same quantity: that closure's determinant IS this m_alpha under a
    different interslice assumption. Measured over every Bishop and Janbu
    solution the corpus produces on its own surfaces (570 solves), the nearest
    an accepted answer comes to the bar is 0.2439, five times it. The Ordinary
    method's normal does not divide by m_alpha and is not tested on it.
  * base tension past `solve.MAX_BASE_TENSION_EXTENT` of the slices refuses the
    answer, and below that extent is reported on the results dict — the same
    treatment `force_equilibrium` and `mprice` already give it. No accepted
    corpus answer carries any.

WHAT IS CHECKED

* the surface the defect was found on: Janbu refused with a reason that names
  m_alpha and the count of slices past it, where it used to report 0.8152, and
  the four complete-equilibrium methods still answer on the same slices;
* the refusal classifies as `inadmissible`, so a search reports it by class;
* every shipped sample's own surfaces still solve under all three methods, and
  the closest any of them comes to the bar is far above it;
* base tension below the extent is reported and not refused, on a purpose-built
  slice set;
* two mutations: the m_alpha test removed, which puts 0.8152 back, and the
  extent test inverted.

Run directly:  PYTHONPATH=. python3 test/base_normal_sign_check.py
"""

import glob
import os
import warnings

import numpy as np

warnings.filterwarnings('ignore')

from xslope import solve
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SAMPLES = os.path.join(REPO, 'docs', 'lem', 'files')

NUM_SLICES = 40

#: The surface `noncircular_search` reaches on `xslope_noncircular.xlsx` once the
#: base-angle cap is opened past 65 degrees, kept as geometry so the solver is
#: checked without re-running a search. Its entry ramp is one straight segment at
#: alpha = -59.96 degrees in a phi = 33 material.
DEFECT_MODEL = os.path.join(SAMPLES, 'xslope_noncircular.xlsx')
DEFECT_POINTS = [
    (0.0, 0.0),
    (3.3533597789999994, -5.8),
    (11.68815851153653, -5.8),
    (24.22341144119898, 8.074470480399661),
]

#: What Janbu used to report there, and what the four complete-equilibrium
#: methods read on the identical slices.
DEFECT_JANBU = 0.8152
DEFECT_COMPANIONS = {'corps': 1.4036, 'lowe': 1.6107,
                     'spencer': 4.3825, 'mprice': 3.4116}


def _poly(points):
    out = [{'X': float(x), 'Y': float(y), 'Movement': 'Horiz'} for x, y in points]
    out[0]['Movement'] = 'Free'
    out[-1]['Movement'] = 'Free'
    return out


def _defect_slices():
    sd = load_slope_data(DEFECT_MODEL)
    ok, res = generate_slices(sd, non_circ=_poly(DEFECT_POINTS),
                              num_slices=NUM_SLICES, debug=False)
    if not ok:
        raise AssertionError(f"could not slice the defect surface: {res}")
    return res[0]


def _m_alpha(df, F):
    a = np.radians(df['alpha'].values)
    tp = np.tan(np.radians(df['phi'].values))
    return np.cos(a) + np.sin(a) * tp / F


def leg_the_singular_surface_is_refused():
    """Janbu must not report a factor of safety from past its own singularity."""
    fails = []
    df = _defect_slices()
    ok, res = solve.janbu(df.copy())
    if ok:
        fails.append(f"janbu reported {res['FS']:.4f} on a surface whose m_alpha "
                     f"reaches {_m_alpha(df, res.get('FS_base', res['FS'])).min():.4f}")
    elif 'm_alpha' not in str(res):
        fails.append(f"janbu refused but the reason does not name m_alpha: {res}")
    elif solve.failure_kind(str(res)) != 'inadmissible':
        fails.append(f"the refusal classifies as "
                     f"{solve.failure_kind(str(res))!r}, expected 'inadmissible'")
    else:
        print(f"  janbu    refused: {str(res).split('(', 1)[1].rstrip(')')}")

    # The surface itself is still solvable — this is a defect in one method, not
    # a slice set nothing can answer.
    for method, expected in DEFECT_COMPANIONS.items():
        ok_c, res_c = getattr(solve, method)(df.copy())
        if not ok_c:
            fails.append(f"{method}: no longer answers on the defect surface "
                         f"({res_c}); the case no longer isolates Janbu")
        elif abs(res_c['FS'] - expected) > 5e-3:
            fails.append(f"{method}: {res_c['FS']:.4f} against the {expected:.4f} "
                         f"this surface is pinned at")
    if not fails:
        print("  companions " + ", ".join(
            f"{m} {v:.3f}" for m, v in sorted(DEFECT_COMPANIONS.items()))
            + f" — against Janbu's old {DEFECT_JANBU}")
    return fails


def leg_the_corpus_stands_clear_of_the_bar():
    """Every sample's own surfaces solve, and not by a narrow margin."""
    fails = []
    worst = (np.inf, None)
    solves = 0
    for path in sorted(glob.glob(os.path.join(SAMPLES, '*.xlsx'))):
        try:
            sd = load_slope_data(path)
        except Exception:
            continue
        surfaces = []
        if sd.get('non_circ'):
            surfaces.append(('nc', dict(non_circ=sd['non_circ'])))
        for i, circle in enumerate((sd.get('circles') or [])[:2]):
            surfaces.append((f'c{i}', dict(circle=circle)))
        for tag, kwargs in surfaces:
            try:
                ok, res = generate_slices(sd, num_slices=NUM_SLICES, debug=False,
                                          **kwargs)
            except Exception:
                continue
            if not ok:
                continue
            df = res[0]
            for method in ('bishop', 'janbu'):
                ok_m, res_m = getattr(solve, method)(df.copy())
                if not ok_m:
                    if 'inadmissible solution' in str(res_m):
                        fails.append(f"{os.path.basename(path)}|{tag}: {method} "
                                     f"refuses a shipped sample's own surface "
                                     f"({res_m})")
                    continue
                solves += 1
                F = res_m.get('FS_base', res_m['FS'])
                m = float(_m_alpha(df, F).min())
                if m < worst[0]:
                    worst = (m, f"{os.path.basename(path)}|{tag} {method}")
    if solves < 30:
        fails.append(f"only {solves} sample solves reached; this leg no longer "
                     f"measures the margin")
    elif worst[0] < 2.0 * solve.FE_MIN_BASE_FACTOR:
        fails.append(f"the corpus now comes within {worst[0]:.4f} of the "
                     f"{solve.FE_MIN_BASE_FACTOR} bar ({worst[1]}), so the bar is "
                     f"no longer clear of the answers it must not refuse")
    else:
        print(f"  {solves} sample solves, nearest m_alpha {worst[0]:.4f} "
              f"({worst[1]}) against a bar of {solve.FE_MIN_BASE_FACTOR}")
    return fails


def leg_base_tension_below_the_extent_is_reported():
    """A few base normals in tension are a note on the answer, not a refusal."""
    fails = []
    n = 20
    N = np.full(n, 100.0)
    N[:3] = -5.0                                    # 15%, under the extent
    reason, warns = solve._base_normal_admissibility(
        "Bishop's Simplified Method", N)
    if reason is not None:
        fails.append(f"15% base tension was refused: {reason}")
    elif not warns or 'base tension' not in warns[0]:
        fails.append(f"15% base tension was not reported: {warns}")
    else:
        print(f"  under the extent: {warns[0]}")

    N_bad = np.full(n, 100.0)
    N_bad[:15] = -5.0                               # 75%, past the extent
    reason, warns = solve._base_normal_admissibility(
        "Bishop's Simplified Method", N_bad)
    if reason is None:
        fails.append("75% base tension was accepted")
    elif solve.failure_kind(reason) != 'inadmissible':
        fails.append(f"the extent refusal classifies as "
                     f"{solve.failure_kind(reason)!r}")
    else:
        print(f"  past the extent: {reason}")
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
        print(f"  mutation  {label} -> caught ({len(caught)} failure(s))")


def leg_mutations():
    """Removing either test must put its own failure back."""
    fails = []
    original = solve._base_normal_admissibility

    def no_m_alpha_test(method_name, N_eff, m_alpha=None):
        """The engine as it stood: no test on the sign of the base normal."""
        return None, []

    def inverted_extent(method_name, N_eff, m_alpha=None):
        """The extent test inverted: tension refused only when there is none."""
        N_eff = np.asarray(N_eff, dtype=float)
        if N_eff.size and float((N_eff < 0).mean()) <= solve.MAX_BASE_TENSION_EXTENT:
            return f"{method_name}{solve.SIMPLIFIED_INADMISSIBLE}(inverted)", []
        return None, []

    _mutation("the m_alpha test removed",
              lambda: setattr(solve, '_base_normal_admissibility', no_m_alpha_test),
              lambda: setattr(solve, '_base_normal_admissibility', original),
              leg_the_singular_surface_is_refused, fails)
    _mutation("the extent test inverted",
              lambda: setattr(solve, '_base_normal_admissibility', inverted_extent),
              lambda: setattr(solve, '_base_normal_admissibility', original),
              leg_base_tension_below_the_extent_is_reported, fails)
    return fails


LEGS = [
    ("the singular surface is refused", leg_the_singular_surface_is_refused),
    ("the corpus stands clear of the bar", leg_the_corpus_stands_clear_of_the_bar),
    ("base tension below the extent is reported",
     leg_base_tension_below_the_extent_is_reported),
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
    print("\nbase-normal sign: all legs pass")
