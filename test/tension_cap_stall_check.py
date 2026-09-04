"""A tension cap that cannot be enforced must not decide a factor of safety.

A material with no tensile strength has no admissible state at a sharp
free-surface tip. At the toe of an embankment the elastic solution puts the
soil there in tension, and the discretization cannot open the crack that would
relieve it: the Gauss point relaxes for a while, stops, and then flows forever
at constant stress. Nothing about the slope is happening — on the model below
the displacement field is frozen to seven digits over 33,000 iterations and
7,871 of its 7,881 nodes reach equilibrium at 1e-8 — but the per-node
out-of-balance rate stays above the force tolerance at the ten nodes around the
tip, and the Newton corrector has no equilibrium to find. Left to decide the
trial, four elements out of 3,923 walk the bisection down from 1.230 to 0.696.

The model is the seepage dam of Tutorial COMBO-1 (`xslope_johnson_res.xlsx`) at
the mesh, bracket and tolerance the page's own tag pins, with `t_cut = 0`
declared on its three materials in memory. That cap is very nearly inert on this
section: the shell's Mohr-Coulomb apex already limits tension to 143 psf, the
same run on a coarser mesh moves by 0.002, and Spencer reads 1.248 on the same
geometry. So the capped answer has to be the blank-cap answer, and this check
pins it at 1% — which is a fifth of the distance to the 0.696 the defect gave.

Two legs beside it say the release is confined to the case it exists for: the
same model with the cap left blank must release nothing at all, and the run must
report the release rather than making it silently.

Cost: about five minutes for the capped bisection, a minute for the blank one.
A build without the release takes half an hour to fail, because every trial then
spends its whole 50,000-iteration ceiling creeping at the toe.

Run directly:  PYTHONPATH=. python3 test/tension_cap_stall_check.py
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# The COMBO-1 tag, read from the page rather than restated here, so the check
# runs the configuration the lock beside it was measured on.
PAGE = 'docs/tutorials/combo01_seepage_stability.md'
BENCHMARK = 'COMBO-1-ssrm'

#: The blank-cap factor of safety on this configuration — the COMBO-1 lock — and
#: the fraction the capped run may differ from it by.
BLANK_CAP_FS = 1.23047
FS_TOL = 0.01

#: What the defect produced, quoted so a failure can be read at a glance.
DEFECT_FS = 0.696


def _case(capped):
    """The COMBO-1 case, with ``t_cut`` set to 0 on every material or left as the
    file has it. Built through the suite's own ``build_fem_ssrm_case`` so the
    mesh, the seepage staging and the solver options are the tag's, and through a
    patched loader so the workbook itself is never touched."""
    import run_tests
    import xslope.fileio as fileio

    tests = run_tests.parse_test_tags(run_tests._repo(PAGE))
    tag = [t for t in tests
           if t.get('benchmark') == BENCHMARK and t.get('type') == 'fem_ssrm']
    if len(tag) != 1:
        raise AssertionError(f"{BENCHMARK} matched {len(tag)} tags on {PAGE}")

    orig = fileio.load_slope_data

    def patched(*a, **k):
        sd = orig(*a, **k)
        if capped:
            for m in sd['materials']:
                m['t_cut'] = 0.0
        return sd

    fileio.load_slope_data = patched
    try:
        return run_tests.build_fem_ssrm_case(tag[0])
    finally:
        fileio.load_slope_data = orig


def run():
    from xslope.fem import solve_ssrm

    failures = []

    fem_data, kwargs, f_min, f_max, tol = _case(capped=True)
    res = solve_ssrm(fem_data, F_min=f_min, F_max=f_max, tolerance=tol,
                     debug_level=0, capture_failure_state=False, **kwargs)
    fs = res.get('FS')
    released = res.get('tension_cap_released', 0)
    print(f"  t_cut = 0:     FS = {fs}, {released} element(s) released")
    if fs is None:
        failures.append(f"the capped bisection returned no factor of safety: "
                        f"{res.get('error')}")
    elif abs(fs - BLANK_CAP_FS) > FS_TOL * BLANK_CAP_FS:
        failures.append(
            f"with t_cut = 0 the dam reads FS = {fs:.5f} against the blank-cap "
            f"{BLANK_CAP_FS} — more than {FS_TOL:.0%} away. A cap this section's "
            f"own Mohr-Coulomb apex almost entirely covers cannot move the "
            f"answer; {DEFECT_FS} is what the toe-tip creep gave before the "
            f"stalled cap was released (see _release_stalled_tension_caps).")
    if not released:
        failures.append(
            "the capped run released no tension cap at all, so either the toe "
            "tip is no longer stalling or the release is not reached — the "
            "factor of safety above is not evidence either way")

    fem_data, kwargs, f_min, f_max, tol = _case(capped=False)
    res = solve_ssrm(fem_data, F_min=f_min, F_max=f_max, tolerance=tol,
                     debug_level=0, capture_failure_state=False, **kwargs)
    blank_released = res.get('tension_cap_released', 0)
    print(f"  t_cut blank:   FS = {res.get('FS')}, "
          f"{blank_released} element(s) released")
    if blank_released:
        failures.append(
            f"the blank-cap run released {blank_released} element(s); a model "
            "that declares no cap has none to release, so the reading is firing "
            "on something other than a stalled cap")
    if res.get('FS') is None or abs(res['FS'] - BLANK_CAP_FS) > FS_TOL * BLANK_CAP_FS:
        failures.append(
            f"the blank-cap run reads FS = {res.get('FS')} against its own lock "
            f"{BLANK_CAP_FS}; the comparison above has no baseline")

    return failures


def main():
    print("Tension cap release — the COMBO-1 dam with t_cut = 0")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nThe cap the toe tip cannot carry is released, and the dam reads the "
          "factor of safety it reads without the cap.")


if __name__ == '__main__':
    main()
