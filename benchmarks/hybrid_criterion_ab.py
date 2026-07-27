"""A/B driver for the SSRM's opt-in HYBRID failure criterion.

Runs one benchmark under both `failure_criterion='non_convergence'` (today's
default) and `failure_criterion='hybrid'` on the SAME mesh and the same solver
options, and prints FS_old vs FS_new together with the per-trial verdict table the
hybrid records (F, converged, stable, verdict, u_ratio, growth, exit reason).

WHY IT EXISTS. Non-convergence is a statement about the solver, not about the slope.
The hybrid keeps non-convergence as the trigger but requires displacement evidence
before a non-converged trial counts as a failed slope (see
`xslope.fem.classify_nonconvergence`). Whether that changes any locked factor of
safety is an empirical question; this driver is how it is asked, one benchmark at a
time, without touching the default.

CASES (--case):

  griffiths1     Griffiths & Lane (1999) Example 1, the flagship SSRM benchmark
                 whose displacement-vs-F sweep (docs/fem/images/griffiths1_sweep.png)
                 is the design template for the criterion. Both criteria should land
                 in the published 1.36-1.40 band.

  rs2_62c        RS2-62c as it stands in the corpus, with the vendor tensile
                 strengths in place. Its verdicts sit on genuine displacement
                 runaway, so both criteria should agree near 0.769.

  rs2_62c_nocut  RS2-62c with the per-material t_cut caps STRIPPED IN MEMORY
                 (tension_cutoff_by_material={}) — the pre-fix configuration in which
                 Mohr-Coulomb hands the cap soil fictitious tension. Under the legacy
                 criterion this model's verdicts wandered because trials were stuck
                 rather than failing; the hybrid should recover the equilibrium
                 boundary of that (unphysical but self-consistent) uncapped physics.
                 No corpus file is edited: the strip is a kwarg on this run only.

The mesh, refinement and solver options come from the SAME test tag the regression
suite locks, via run_tests.build_fem_ssrm_case — this driver cannot benchmark a
different model than the lock.

Run from the repo root:

    PYTHONPATH=. python3 benchmarks/hybrid_criterion_ab.py --case griffiths1
    PYTHONPATH=. python3 benchmarks/hybrid_criterion_ab.py --case rs2_62c
    PYTHONPATH=. python3 benchmarks/hybrid_criterion_ab.py --case rs2_62c_nocut
    PYTHONPATH=. python3 benchmarks/hybrid_criterion_ab.py --case all --json out.json

This is a REPORTING driver, not a guard: it always exits 0 unless a solve raises.
The corpus-wide A/B is a separate, owner-gated decision.
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from run_tests import parse_test_tags, build_fem_ssrm_case   # noqa: E402
from xslope.fem import solve_ssrm                            # noqa: E402

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'docs')

# name -> (docs markdown file, benchmark tag, kwarg overrides, tag overrides)
CASES = {
    'griffiths1': (
        os.path.join(DOCS, 'verification', 'ssrm.md'), 'SSRM-1', {}, {},
    ),
    'rs2_62c': (
        os.path.join(DOCS, 'verification', 'rs2.md'), 'RS2-62c', {}, {},
    ),
    'rs2_62c_nocut': (
        os.path.join(DOCS, 'verification', 'rs2.md'), 'RS2-62c',
        # Explicit {} DISABLES the file-carried t_cut column (see solve_ssrm's
        # template-defaults resolution) — the pre-fix, no-tensile-cap physics.
        {'tension_cutoff_by_material': {}},
        # The uncapped model equilibrates far above the capped one, so give the
        # bisection a bracket that reaches it, and a smaller budget: the point here
        # is which trials are STUCK versus FAILING, and that is settled long before
        # 40k iterations (the stuck states are frozen from ~10k onward).
        {'f_min': 0.5, 'f_max': 2.0, 'max_iter': '10000', 'tolerance': 0.05},
    ),
}


def _find_tag(md_path, benchmark):
    for t in parse_test_tags(md_path):
        if t.get('benchmark') == benchmark and t.get('type') == 'fem_ssrm':
            return dict(t)
    raise SystemExit(f"no fem_ssrm tag with benchmark={benchmark!r} in {md_path}")


def _fmt(v, spec='.3f'):
    return '   -  ' if v is None else format(v, spec)


def _trial_table(trials):
    lines = ["    F      role    conv  stable  verdict        u/u_el   growth   exit          iters",
             "    " + "-" * 84]
    for t in trials:
        lines.append(
            f"    {t['F']:<6.3f} {t['role']:<7s} {str(t['converged']):<5s} "
            f"{str(t['stable']):<7s} {str(t['verdict'] or '-'):<14s} "
            f"{_fmt(t['u_ratio'], '6.2f')}  {_fmt(t['growth'], '+7.3f')}  "
            f"{str(t['exit_reason'] or '-'):<13s} {t['iterations']}")
    return "\n".join(lines)


def run_case(name, max_iter_override=None, verbose=0):
    md_path, benchmark, kw_over, tag_over = CASES[name]
    test = _find_tag(md_path, benchmark)
    test.update(tag_over)
    if max_iter_override:
        test['max_iter'] = str(max_iter_override)

    print(f"\n{'=' * 88}\nCASE {name}  (tag benchmark={benchmark}, file={test['file']})")
    print(f"  mesh: {test.get('element_type')} target_size={test.get('target_size')} "
          f"refine={test.get('refine_factor', '-')}/{test.get('refine_features', '-')}")
    print(f"  bracket [{test.get('f_min')}, {test.get('f_max')}] "
          f"tol={test.get('tolerance')} max_iter={test.get('max_iter')}"
          + (f"  overrides={kw_over}" if kw_over else ""))
    if test.get('expected_fs') is not None:
        print(f"  locked FS = {test['expected_fs']}")

    fem_data, kwargs, f_min, f_max, ssrm_tol = build_fem_ssrm_case(test)
    kwargs.update(kw_over)
    # The at-failure capture is a post-processing extra that does not touch FS or the
    # bracket; turn it off so the A/B pays only for the bisection.
    kwargs['capture_failure_state'] = False

    out = {'case': name, 'benchmark': benchmark, 'file': test['file'],
           'expected_fs': test.get('expected_fs'), 'runs': {}}

    for criterion in ('non_convergence', 'hybrid'):
        t0 = time.perf_counter()
        res = solve_ssrm(fem_data, F_min=f_min, F_max=f_max, tolerance=ssrm_tol,
                         debug_level=verbose, failure_criterion=criterion, **kwargs)
        dt = time.perf_counter() - t0
        trials = res.get('trials', [])
        print(f"\n  --- {criterion} ---")
        if res.get('FS') is None:
            print(f"    FS = (none): {res.get('error', 'unknown')}")
        else:
            print(f"    FS = {res['FS']:.4f}   bracket {res['final_interval']}   "
                  f"{len(trials)} trials, {dt:.0f}s")
        if trials:
            print(_trial_table(trials))
        out['runs'][criterion] = {
            'FS': res.get('FS'),
            'final_interval': list(res['final_interval']) if res.get('final_interval') else None,
            'error': res.get('error'),
            'elapsed_s': dt,
            'trials': trials,
        }

    fs_old = out['runs']['non_convergence']['FS']
    fs_new = out['runs']['hybrid']['FS']
    n_rescued = sum(1 for t in out['runs']['hybrid']['trials']
                    if t.get('verdict') == 'STABLE_STUCK')
    n_amb = sum(1 for t in out['runs']['hybrid']['trials']
                if t.get('verdict') == 'AMBIGUOUS')
    out['delta'] = None if (fs_old is None or fs_new is None) else fs_new - fs_old
    out['n_stable_stuck'] = n_rescued
    out['n_ambiguous'] = n_amb
    print(f"\n  SUMMARY {name}: FS_old = {_fmt(fs_old, '.4f')}   FS_new = {_fmt(fs_new, '.4f')}   "
          f"delta = {_fmt(out['delta'], '+.4f')}   "
          f"(hybrid: {n_rescued} STABLE_STUCK, {n_amb} AMBIGUOUS)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--case', default='griffiths1',
                    choices=list(CASES) + ['all'])
    ap.add_argument('--max-iter', type=int, default=None,
                    help='override the tag iteration ceiling (both criteria)')
    ap.add_argument('--json', default=None, help='write the full report to this path')
    ap.add_argument('-v', '--verbose', action='count', default=0,
                    help='solve_ssrm debug_level (repeat for more)')
    args = ap.parse_args()

    names = list(CASES) if args.case == 'all' else [args.case]
    report = [run_case(n, args.max_iter, args.verbose) for n in names]

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(report, fh, indent=2)
        print(f"\nwrote {args.json}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
