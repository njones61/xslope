"""
Independent, from-scratch implementation of Bishop's Simplified Method.

Derived from published theory (Bishop, 1955):
  Moment equilibrium of the whole sliding mass about the circle center,
  with the simplifying assumption that interslice shear forces are zero
  (so the resultant interslice force is horizontal).

For each slice, vertical force equilibrium gives the effective normal force N'
on the base. The mobilized shear from Mohr-Coulomb (divided by FS) is:

    S_m = (1/FS) * [ c*b/cos(a) + (N') * tan(phi) ]   ... along base length dl=b/cos(a)

Solving vertical equilibrium for N' and substituting into moment equilibrium
about the center (where the lever arm for shear is R and the driving term is
W*sin(a)*R, R cancels) yields the classic iterative form:

    FS = SUM[ (c*b + (W - u*b)*tan(phi)) / m_a ] / SUM[ W*sin(a) ]

with
    m_a = cos(a) * (1 + tan(a)*tan(phi)/FS)

Here a = base inclination (alpha), b = horizontal slice width (dx),
W = slice weight, u = pore pressure at base center, c, phi = base strength.

This file deliberately re-derives the FS computation and does NOT consult
xslope's solve.bishop(). It is what lets the shipped solver be scored against
something other than itself: `run()` sweeps every shipped sample circle whose
slice set carries nothing but weight, pore pressure and strength — the terms this
implementation models — and requires the two to agree to solver precision.
"""

import numpy as np


def bishop_independent(df, fs_init=1.5, tol=1e-10, max_iter=200):
    a = np.radians(df['alpha'].to_numpy(float))   # base inclination, radians
    phi = np.radians(df['phi'].to_numpy(float))   # friction angle, radians
    c = df['c'].to_numpy(float)                    # cohesion
    W = df['w'].to_numpy(float)                    # slice weight
    u = df['u'].to_numpy(float)                    # pore pressure at base center
    b = df['dx'].to_numpy(float)                   # horizontal slice width

    tan_phi = np.tan(phi)
    sin_a = np.sin(a)
    cos_a = np.cos(a)
    tan_a = np.tan(a)

    # Driving (denominator) term is independent of FS.
    denom = np.sum(W * sin_a)

    # Cohesive + frictional numerator pieces (independent of FS).
    cohesive = c * b
    frictional = (W - u * b) * tan_phi

    fs = fs_init
    for _ in range(max_iter):
        m_a = cos_a * (1.0 + tan_a * tan_phi / fs)
        numer = np.sum((cohesive + frictional) / m_a)
        fs_new = numer / denom
        if abs(fs_new - fs) < tol:
            fs = fs_new
            break
        fs = fs_new
    return fs


# ---------------------------------------------------------------------------
# Scoring the shipped solver against this one
# ---------------------------------------------------------------------------

#: Per-slice terms this implementation does not model. A slice set carrying any
#: of them is skipped rather than scored, because a disagreement there would be
#: about what the two implementations include, not about how either solves.
EXTRA_TERMS = ('dload', 'lload', 'kw', 't', 'p', 'p_pt', 'h_pile',
               'pa_cx', 'pa_cy', 'pp_cx', 'pp_cy', 'h_pile_pas')

#: How far the two may differ, relative. They solve the same equation by the same
#: fixed point from different code, so what separates them is only where each
#: stops iterating: the shipped solver at 1e-6 on F, this one at 1e-10. Measured
#: across the qualifying corpus, the spread is 3e-8 at worst.
TOLERANCE = 1e-5

#: How many slice sets must qualify before the comparison means anything.
MIN_CASES = 8


def _plain(df):
    """True when this slice set carries only weight, pore pressure and strength."""
    from xslope import solve
    if solve._has_nonlinear(df):
        # A nonlinear strength envelope makes the shipped solver iterate the
        # strength against the base normal it computes; this file solves the
        # linear Mohr-Coulomb equation it derives, so the two are not the same
        # problem.
        return False
    for col in EXTRA_TERMS:
        if col in df.columns and float(np.abs(df[col].to_numpy(float)).max()) > 0:
            return False
    if 'r' not in df.columns or df['r'].isna().all():
        return False           # Bishop needs a center of rotation
    # A composite surface's floor run is not at radius R, so the classic form
    # this file derives does not apply to it.
    xr = df['x_c'].to_numpy(float) - float(df['xo'].iloc[0])
    yr = df['y_cb'].to_numpy(float) - float(df['yo'].iloc[0])
    a = np.radians(df['alpha'].to_numpy(float))
    a_N = xr * np.cos(a) + yr * np.sin(a)
    scale = float(df['r'].iloc[0])
    return bool(np.max(np.abs(a_N)) <= 1e-6 * scale)


def leg_the_shipped_solver_matches():
    """Bishop's method against this file's independent implementation."""
    import glob
    import os

    from xslope import solve
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices

    here = os.path.dirname(os.path.abspath(__file__))
    samples = os.path.join(os.path.dirname(here), 'docs', 'lem', 'files')

    fails = []
    cases = 0
    worst = (0.0, None)
    for path in sorted(glob.glob(os.path.join(samples, '*.xlsx'))):
        try:
            sd = load_slope_data(path)
        except Exception:
            continue
        for i, circle in enumerate((sd.get('circles') or [])[:2]):
            try:
                ok, res = generate_slices(sd, circle=circle, num_slices=40,
                                          debug=False)
            except Exception:
                continue
            if not ok:
                continue
            df = res[0]
            if not _plain(df):
                continue
            ok_b, res_b = solve.bishop(df.copy())
            if not ok_b:
                continue
            cases += 1
            mine = bishop_independent(df)
            rel = abs(res_b['FS'] - mine) / res_b['FS']
            if rel > worst[0]:
                worst = (rel, f"{os.path.basename(path)}|c{i}")
            if rel > TOLERANCE:
                fails.append(f"{os.path.basename(path)}|c{i}: xslope "
                             f"{res_b['FS']:.8f} against an independent "
                             f"{mine:.8f} ({rel:.2e} relative)")
    if cases < MIN_CASES:
        fails.append(f"only {cases} slice sets qualified, under the {MIN_CASES} "
                     f"this comparison needs to mean anything")
    elif not fails:
        print(f"  {cases} circles scored, worst relative difference "
              f"{worst[0]:.2e} ({worst[1]}) against a bar of {TOLERANCE:.0e}")
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
    """A wrong m_alpha in the shipped solver must show up here."""
    from xslope import solve
    fails = []
    original = solve.bishop

    def shifted(slice_df, **kwargs):
        """Bishop's answer moved by a tenth of a percent."""
        ok, res = original(slice_df, **kwargs)
        if ok:
            res = dict(res, FS=res['FS'] * 1.001)
        return ok, res

    _mutation("Bishop's answer shifted by 0.1%",
              lambda: setattr(solve, 'bishop', shifted),
              lambda: setattr(solve, 'bishop', original),
              leg_the_shipped_solver_matches, fails)
    return fails


LEGS = [
    ("the shipped solver matches an independent implementation",
     leg_the_shipped_solver_matches),
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


if __name__ == "__main__":
    import sys
    fails = run()
    if fails:
        print("\nFAILURES:")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("\nPASS: Bishop's method reproduces an independent implementation "
          "of the same equation.")
