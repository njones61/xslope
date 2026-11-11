"""Spencer on a sliver: the residual bar scales with the mass.

Spencer's Newton iteration accepts a solution when both residuals fall below
`tol` = 1e-4: R1, the sum of the interslice resultants (a force), and R2, their
moment. Both are in the model's own units, so the bar is a different precision on
every model. On an ordinary slope, whose vertical forces sum to thousands of lb
or kN per unit width, it is a part in ten million or finer. On LEM-10's shallow
critical surface, a 0.33 ft sliver on the sand face weighing 0.009 lb/ft, it was
1% of the weight, and the iteration stopped 2.4% short of the root: the shallow
search reported 1.268 where the infinite-slope limit tan(phi)/tan(beta) is
1.299, and Bishop, OMS, Janbu and Morgenstern-Price all return 1.299 on the same
slices. The fault showed when Spencer's starting inclination moved from a fixed
8 degrees to the slip surface's weight-averaged base inclination: from 8 degrees
Newton had far to travel and happened to land near the root, while from the
sliver's own inclination it met the loose bar almost at once.

Below SPENCER_TOL_REF (1 force unit for R1, 1 force-length unit for R2) the bar
now shrinks in proportion to the mass; at and above it the bar is `tol`, so a
model of real size solves exactly as before.

WHAT IS CHECKED

* the shallow search on the tutorial's model, from LEM-10's embankment circle
  (Xo 16.875, Yo 30, Depth 0) and from the CE 544 exercise's (Xo 11.25, Yo 40,
  Depth 0), returns the infinite-slope limit, 1.299, to 0.002, and never falls
  more than 0.001 below it (with the fixed bar they read 1.2985 and 1.268);
* on each search's critical sliver Spencer agrees with Bishop to 0.001;
* a mutation that restores the fixed bar (SPENCER_TOL_REF = 1e-30, so both
  factors are 1) brings 1.268 back on the exercise's circle, so the first leg
  can fail;
* on the tutorial's deep circle, a model of real size, Spencer returns the same
  factor of safety and inclination to the last bit with the scaling and without
  it.

Run directly:  PYTHONPATH=. python3 test/spencer_sliver_check.py
"""

import contextlib
import io
import math
import os
import warnings

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.path.join(HERE, '..', 'docs', 'lem', 'files', 'xslope_mult_min_KEY.xlsx')
SHALLOW = ({'Xo': 16.875, 'Yo': 30.0, 'Depth': 0.0, 'R': 30.0},   # LEM-10's embankment circle
           {'Xo': 11.25, 'Yo': 40.0, 'Depth': 0.0, 'R': 40.0})    # the CE 544 two-minima exercise's
LIMIT = math.tan(math.radians(30.0)) / (15.0 / 33.75)        # tan(phi) / tan(beta), the 2.25:1 sand face


def _quiet(fn, *a, **k):
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter('ignore')
        return fn(*a, **k)


def _shallow_search(sd, circle):
    from xslope.search import run_lem_analysis
    return _quiet(run_lem_analysis, dict(sd, circles=[dict(circle)], circular=True), 'spencer',
                  analysis='auto_search', surface='circular', num_slices=40)


def run():
    from xslope import solve
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices

    failures = []
    sd = _quiet(load_slope_data, MODEL)

    for circle in SHALLOW:
        tag = f"shallow search from ({circle['Xo']:g}, {circle['Yo']:g})"
        b = _shallow_search(sd, circle)
        F = b['results']['FS']
        if abs(F - LIMIT) > 0.002 or F < LIMIT - 0.001:
            failures.append(f"{tag}: Spencer {F:.4f}, the infinite-slope limit is {LIMIT:.4f}")
        df = b['slice_df']
        ok_s, r_s = _quiet(solve.spencer, df.copy())
        ok_b, r_b = _quiet(solve.bishop, df.copy())
        if not (ok_s and ok_b):
            failures.append(f"{tag}: Spencer or Bishop did not solve the critical sliver")
        elif abs(r_s['FS'] - r_b['FS']) > 0.001:
            failures.append(f"{tag}: on the critical sliver Spencer {r_s['FS']:.4f}, Bishop {r_b['FS']:.4f}")

    saved = solve.SPENCER_TOL_REF
    try:
        solve.SPENCER_TOL_REF = 1e-30          # the fixed bar, as before the fix
        F_m = _shallow_search(sd, SHALLOW[1])['results']['FS']
        if F_m > LIMIT - 0.01:
            failures.append(f"mutation: the fixed bar gave {F_m:.4f} from the exercise's circle, so the "
                            f"shallow-search leg cannot tell the fault from the fix")
    finally:
        solve.SPENCER_TOL_REF = saved

    ok, res = _quiet(generate_slices, sd, circle=sd['circles'][0], num_slices=40)
    if not ok:
        failures.append(f"deep circle: slicing failed ({res})")
    else:
        deep = res[0]
        ok_a, r_a = _quiet(solve.spencer, deep.copy())
        try:
            solve.SPENCER_TOL_REF = 1e-30
            ok_o, r_o = _quiet(solve.spencer, deep.copy())
        finally:
            solve.SPENCER_TOL_REF = saved
        if not (ok_a and ok_o):
            failures.append("deep circle: Spencer did not solve it")
        elif r_a['FS'] != r_o['FS'] or r_a['theta'] != r_o['theta']:
            failures.append(f"deep circle: {r_a['FS']!r} with the scaled bar, {r_o['FS']!r} without; "
                            f"a model of real size must solve exactly as before")
    return failures


if __name__ == '__main__':
    f = run()
    print("PASS" if not f else "FAIL\n  " + "\n  ".join(f))
