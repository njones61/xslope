#!/usr/bin/env python3
"""
Regression test suite for xslope.

Scans docs/{lem,fem,seep}/samples.md for test tags of the form:

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
import re
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # non-interactive backend — no plot windows


def parse_test_tags(md_path):
    """Parse <!-- test: ... --> tags from a markdown file."""
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
        for key in ['expected_fs', 'expected_flowrate', 'tolerance', 'target_size', 'f_min', 'f_max']:
            if key in params:
                params[key] = float(params[key])
        if 'num_slices' in params:
            params['num_slices'] = int(params['num_slices'])

        params['source'] = str(md_path)
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
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons, extract_reinforcement_line_geometry

    file_path = test['file']
    element_type = test.get('element_type', 'tri6')
    target_size = test.get('target_size', 2.0)
    ssrm_tolerance = test.get('tolerance', 0.05)

    slope_data = load_slope_data(file_path)

    reinf_geom = extract_reinforcement_line_geometry(slope_data)
    polygons = get_material_polygons(slope_data, reinf_lines=reinf_geom)
    mesh = build_mesh_from_polygons(
        polygons, target_size=target_size, element_type=element_type, lines=reinf_geom
    )

    fem_data = build_fem_data(slope_data, mesh)
    f_min = test.get('f_min', 0.5)
    f_max = test.get('f_max', 3.0)
    result = solve_ssrm(fem_data, F_min=f_min, F_max=f_max, tolerance=ssrm_tolerance, debug_level=0)

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

    return solution['flowrate'], None


def run_test(test):
    """Run a single test and return (computed_value, error_msg)."""
    test_type = test.get('type', '')
    if test_type == 'fem_ssrm':
        return run_fem_test(test)
    elif test_type == 'seep':
        return run_seep_test(test)
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

    if args.skip_benchmarks:
        n_before = len(tests)
        tests = [t for t in tests if 'benchmark' not in t]
        n_skipped = n_before - len(tests)
        if n_skipped:
            print(f"Skipping {n_skipped} benchmark tests (--skip-benchmarks)")

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
