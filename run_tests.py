#!/usr/bin/env python3
"""
Regression test suite for xslope.

Scans docs/{lem,fem,seep}/samples.md and docs/seep/seep_slope.md for test
tags of the form:

    <!-- test: file=files/foo.xlsx, type=circular_search, method=spencer, expected_fs=1.234, num_slices=30 -->
    <!-- test: file=files/foo.xlsx, type=fem_ssrm, expected_fs=1.38, element_type=quad8, target_size=3.5, tolerance=0.025 -->
    <!-- test: file=files/foo.xlsx, type=seep, expected_flowrate=40.062, tolerance=0.05 -->
    <!-- test: file=files/foo.xlsx, type=seep, expected_flowrate=28.6, element_type=tri6, target_size=2.0, tolerance=0.01 -->

An optional benchmark=<ID> key (e.g. benchmark=SSRM-2) marks a verification
benchmark; these can be excluded with --skip-benchmarks for quicker routine
runs.

Runs each test and compares results against expected values.
As new samples are added with test tags, they automatically become part of the suite.

Usage:
    python run_tests.py              # run all tests
    python run_tests.py --lem        # run only LEM tests
    python run_tests.py --fem        # run only FEM tests
    python run_tests.py --seep       # run only seepage tests
    python run_tests.py --tolerance 0.02  # custom FS tolerance (default 0.01)
    python run_tests.py --skip-benchmarks # exclude verification benchmarks (faster)
    python run_tests.py --verbose    # print details for passing tests too
"""

import argparse
import os
import re
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend — no plot windows


# Short method names used in compact fs_<method> test tags -> solver names.
FS_METHOD_NAMES = {
    'oms': 'oms', 'bishop': 'bishop', 'janbu': 'janbu',
    'corps': 'corps_engineers', 'lowe': 'lowe_karafiath', 'spencer': 'spencer',
}


def parse_test_tags(md_path):
    """Parse <!-- test: ... --> tags from a markdown file.

    A tag may name a single method the classic way (method=spencer,
    expected_fs=1.23) or carry per-method values (fs_oms=0.94, fs_bishop=0.99,
    ...). The latter form is expanded into one test per method, so every method
    becomes an independently-checked regression.
    """
    tests = []
    md_dir = Path(md_path).parent

    with open(md_path, 'r') as f:
        content = f.read()

    pattern = r'<!--\s*test:\s*(.*?)\s*-->'
    for match in re.finditer(pattern, content):
        tag_str = match.group(1)
        params = {}
        for pair in tag_str.split(','):
            key, _, value = pair.partition('=')
            params[key.strip()] = value.strip()

        # Resolve file path relative to the markdown file's directory
        if 'file' in params:
            params['file'] = str(md_dir / params['file'])

        # Convert numeric fields
        for key in ['expected_fs', 'expected_flowrate', 'expected_beta', 'tolerance', 'target_size', 'f_min', 'f_max']:
            if key in params:
                params[key] = float(params[key])
        if 'num_slices' in params:
            params['num_slices'] = int(params['num_slices'])

        params['source'] = str(md_path)

        # Expand a compact multi-method tag (fs_oms=..., fs_bishop=..., ...) into
        # one test per method; otherwise keep the single tag as-is.
        fs_keys = [k for k in params if k.startswith('fs_')]
        if fs_keys:
            shared = {k: v for k, v in params.items() if not k.startswith('fs_')}
            for k in fs_keys:
                method = FS_METHOD_NAMES.get(k[3:], k[3:])
                t = dict(shared)
                t['method'] = method
                t['expected_fs'] = float(params[k])
                tests.append(t)
        else:
            tests.append(params)

    return tests


def run_lem_test(test):
    """Run a single LEM test (single_circle, circular_search, or noncircular_search)."""
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope.solve import solve_selected
    from xslope.search import circular_search, noncircular_search

    file_path = test['file']
    test_type = test['type']
    method = test['method']
    num_slices = test.get('num_slices', 30)

    slope_data = load_slope_data(file_path)

    if test_type == 'single_circle':
        circle = slope_data['circles'][0]
        success, result = generate_slices(slope_data, circle=circle, num_slices=num_slices)
        if not success:
            return None, f"generate_slices failed: {result}"
        slice_df, failure_surface = result
        solver_result = solve_selected(method, slice_df)
        if isinstance(solver_result, str):
            return None, f"solve failed: {solver_result}"
        return solver_result['FS'], None

    elif test_type == 'circular_search':
        fs_cache, converged, search_path, circle_cache = circular_search(
            slope_data, method, num_slices=num_slices
        )
        if not fs_cache or fs_cache[0]['FS'] >= 9999:
            return None, "circular_search found no valid surface"
        return fs_cache[0]['FS'], None

    elif test_type == 'noncircular_search':
        fs_cache, converged, search_path = noncircular_search(
            slope_data, method, num_slices=num_slices
        )
        if not fs_cache or fs_cache[0]['FS'] >= 9999:
            return None, "noncircular_search found no valid surface"
        return fs_cache[0]['FS'], None

    else:
        return None, f"Unknown LEM test type: {test_type}"


def run_fem_test(test):
    """Run a single FEM SSRM test."""
    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data, solve_ssrm
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons, extract_constraint_line_geometry

    file_path = test['file']
    element_type = test.get('element_type', 'tri6')
    target_size = test.get('target_size')  # default computed from domain extent below
    ssrm_tolerance = test.get('tolerance', 0.05)

    slope_data = load_slope_data(file_path)

    # Seepage-coupled models (material u option = 'seep') must run on the SAME
    # mesh as the stored seepage solution: seep_u is nodal, and build_fem_data
    # silently zeroes the pore pressures if the node count differs.
    uses_seep = any(str(m.get('u', '')).strip().lower() == 'seep'
                    for m in slope_data.get('materials', []))
    if uses_seep and slope_data.get('mesh') is not None:
        mesh = slope_data['mesh']
    else:
        # Default mesh size scales with the domain (like the seepage path) rather
        # than a fixed value, so a wide domain does not blow up to a huge mesh.
        if target_size is None:
            x_coords = [x for x, _ in slope_data['ground_surface'].coords]
            target_size = (max(x_coords) - min(x_coords)) / 100
        # constraint lines include BOTH reinforcement and pile axes — pile
        # beam elements require their axis lines in the mesh (a
        # reinforcement-only extraction silently drops the piles)
        constraint_lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
        polygons = get_material_polygons(slope_data, reinf_lines=constraint_lines)
        mesh = build_mesh_from_polygons(
            polygons, target_size=target_size, element_type=element_type,
            lines=constraint_lines
        )

    fem_data = build_fem_data(slope_data, mesh)
    f_min = test.get('f_min', 0.5)
    f_max = test.get('f_max', 3.0)
    kwargs = {}
    if 'criterion' in test:
        kwargs['failure_criterion'] = test['criterion']
    if 'max_iter' in test:
        kwargs['max_iterations'] = int(test['max_iter'])
    if test.get('cutoff', '').lower() in ('true', '1', 'yes'):
        kwargs['tension_cutoff'] = True
    if 'char_x' in test and 'char_y' in test:
        kwargs['char_point'] = (float(test['char_x']), float(test['char_y']))
    result = solve_ssrm(fem_data, F_min=f_min, F_max=f_max, tolerance=ssrm_tolerance,
                        debug_level=0, **kwargs)

    if result.get('converged', False):
        return result['FS'], None
    else:
        return None, f"SSRM failed: {result.get('error', 'Unknown error')}"


def run_seep_test(test):
    """Run a single seepage test."""
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import build_seep_data, run_seepage_analysis

    file_path = test['file']
    slope_data = load_slope_data(file_path)

    polygons = get_material_polygons(slope_data)
    element_type = test.get('element_type', 'tri3')
    target_size = test.get('target_size')
    if target_size is None:
        x_coords = [x for x, _ in slope_data['ground_surface'].coords]
        target_size = (max(x_coords) - min(x_coords)) / 120
    mesh = build_mesh_from_polygons(polygons, target_size, element_type)

    seep_data = build_seep_data(mesh, slope_data)
    solution = run_seepage_analysis(seep_data, tol=1e-4)

    if not solution.get('converged', True):
        return None, "seepage solution did not converge (flowrate unreliable)"
    return solution['flowrate'], None


def run_reliability_test(test):
    """Run a single reliability analysis, returning the lognormal reliability index beta."""
    from xslope.fileio import load_slope_data
    from xslope.advanced import reliability as reliability_analysis

    file_path = test['file']
    method = test.get('method', 'spencer')
    circular = str(test.get('circular', 'true')).lower() not in ('false', '0', 'no')

    slope_data = load_slope_data(file_path)
    success, result = reliability_analysis(slope_data, method, circular=circular, debug_level=0)
    if not success:
        return None, f"reliability failed: {result}"
    return result['beta_ln'], None


def run_test(test):
    """Run a single test and return (computed_value, error_msg)."""
    test_type = test.get('type', '')
    if test_type == 'fem_ssrm':
        return run_fem_test(test)
    elif test_type == 'seep':
        return run_seep_test(test)
    elif test_type == 'reliability':
        return run_reliability_test(test)
    else:
        return run_lem_test(test)


def main():
    parser = argparse.ArgumentParser(description='xslope regression test suite')
    parser.add_argument('--lem', action='store_true', help='Run only LEM tests')
    parser.add_argument('--fem', action='store_true', help='Run only FEM tests')
    parser.add_argument('--seep', action='store_true', help='Run only seepage tests')
    parser.add_argument('--tolerance', type=float, default=0.01,
                        help='Default tolerance for FS comparison (default: 0.01)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print details for passing tests')
    parser.add_argument('--skip-benchmarks', action='store_true',
                        help='Skip verification benchmark tests (annotations '
                             'carrying a benchmark=<ID> tag), e.g. the slow '
                             'SSRM dam cases')
    parser.add_argument('--quick', action='store_true',
                        help='For LEM problems that list several methods, check '
                             'only one method per problem (prefers Spencer) so '
                             'routine runs stay fast')
    args = parser.parse_args()

    # If no specific flags, run all
    run_all = not args.lem and not args.fem and not args.seep
    run_lem = args.lem or run_all
    run_fem = args.fem or run_all
    run_seep = args.seep or run_all

    # Discover tests from markdown files
    tests = []
    if run_lem:
        lem_samples = Path('docs/lem/samples.md')
        if lem_samples.exists():
            tests.extend(parse_test_tags(lem_samples))
    if run_fem:
        fem_samples = Path('docs/fem/samples.md')
        if fem_samples.exists():
            tests.extend(parse_test_tags(fem_samples))
    if run_seep:
        seep_samples = Path('docs/seep/samples.md')
        if seep_samples.exists():
            tests.extend(parse_test_tags(seep_samples))

    # docs/seep/seep_slope.md mixes seepage, LEM, and FEM tests (the combined
    # Johnson Reservoir example) — always scan it and route each test by type.
    seep_slope = Path('docs/seep/seep_slope.md')
    if seep_slope.exists():
        for t in parse_test_tags(seep_slope):
            ttype = t.get('type', '')
            if ttype == 'fem_ssrm':
                if run_fem:
                    tests.append(t)
            elif ttype == 'seep':
                if run_seep:
                    tests.append(t)
            elif run_lem:
                tests.append(t)

    # Private test problems (CE 544 homework/exam keys) live in a separate
    # repo kept out of the public xslope tree. Look for it as a sibling
    # directory or at $XSLOPE_PRIVATE_TESTS; scan every markdown file there for
    # test tags and route by type. Silently skipped when the repo is absent
    # (public CI, other clones), so the public suite is unaffected.
    private_dir = os.environ.get('XSLOPE_PRIVATE_TESTS') or str(
        Path(__file__).resolve().parent.parent / 'xslope_private_tests')
    private_path = Path(private_dir)
    if private_path.is_dir():
        n_priv = 0
        for md in sorted(private_path.glob('*.md')):
            if md.name.lower() == 'readme.md':
                continue
            for t in parse_test_tags(md):
                ttype = t.get('type', '')
                if ttype == 'fem_ssrm':
                    if run_fem:
                        tests.append(t); n_priv += 1
                elif ttype == 'seep':
                    if run_seep:
                        tests.append(t); n_priv += 1
                elif run_lem:
                    tests.append(t); n_priv += 1
        if n_priv:
            print(f"Including {n_priv} private tests from {private_path.name}/")

    if args.skip_benchmarks:
        n_before = len(tests)
        tests = [t for t in tests if 'benchmark' not in t]
        n_skipped = n_before - len(tests)
        if n_skipped:
            print(f"Skipping {n_skipped} benchmark tests (--skip-benchmarks)")

    if args.quick:
        # Collapse each problem's per-method expansion to a single check
        # (prefer Spencer). Singleton tests (FEM, seep, single-method LEM) pass through.
        from collections import OrderedDict
        groups = OrderedDict()
        for t in tests:
            groups.setdefault((t.get('source'), t.get('file'), t.get('type')), []).append(t)
        tests = [next((t for t in grp if t.get('method') == 'spencer'), grp[0])
                 for grp in groups.values()]
        print(f"--quick: checking one method per problem ({len(tests)} tests)")

    if not tests:
        print("No test tags found in documentation files.")
        sys.exit(1)

    print(f"Found {len(tests)} tests\n")
    print(f"{'#':<4} {'File':<45} {'Type':<20} {'Method':<10} {'Expected':>10}  {'Computed':>10}  {'Status'}")
    print("-" * 120)

    passed = 0
    failed = 0
    errors = 0
    total_time = 0

    for i, test in enumerate(tests, 1):
        file_name = Path(test['file']).name
        test_type = test.get('type', '?')
        method = test.get('method', '-')
        is_seep = test_type == 'seep'

        # Get expected value and tolerance
        if is_seep:
            expected = test.get('expected_flowrate')
            tol = test.get('tolerance', 0.05) * abs(expected) if expected else 0
            label = 'flowrate'
        elif test_type == 'reliability':
            expected = test.get('expected_beta')
            tol = test.get('tolerance', 0.02)
            label = 'beta'
        else:
            expected = test.get('expected_fs')
            tol = test.get('tolerance', args.tolerance)
            label = 'FS'

        t0 = time.time()
        try:
            computed, error_msg = run_test(test)
        except Exception as e:
            computed = None
            error_msg = str(e)
        elapsed = time.time() - t0
        total_time += elapsed

        if error_msg:
            status = f"ERROR: {error_msg}"
            errors += 1
            comp_str = "    --    "
        elif expected is not None and abs(computed - expected) <= tol:
            status = f"PASS ({elapsed:.1f}s)"
            passed += 1
            comp_str = f"{computed:10.3f}"
        else:
            diff = computed - expected if computed is not None else float('nan')
            status = f"FAIL (diff={diff:+.4f}, {elapsed:.1f}s)"
            failed += 1
            comp_str = f"{computed:10.3f}" if computed is not None else "    --    "

        exp_str = f"{expected:10.3f}" if expected is not None else "    --    "
        print(f"{i:<4} {file_name:<45} {test_type:<20} {method:<10} {exp_str}  {comp_str}  {status}")

    print("-" * 120)
    print(f"\nResults: {passed} passed, {failed} failed, {errors} errors out of {len(tests)} tests")
    print(f"Total time: {total_time:.1f}s")

    if failed > 0 or errors > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
