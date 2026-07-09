#!/usr/bin/env python3
"""
Regression test suite for xslope.

Scans docs/{lem,fem,seep}/samples.md and docs/seep/seep_slope.md for test
tags of the form:

    <!-- test: file=files/foo.xlsx, type=circular_search, method=spencer, expected_fs=1.234, num_slices=30 -->
    <!-- test: file=files/foo.xlsx, type=fem_ssrm, expected_fs=1.38, element_type=quad8, target_size=3.5, tolerance=0.025 -->
    <!-- test: file=files/foo.xlsx, type=seep, expected_flowrate=40.062, tolerance=0.05 -->
    <!-- test: file=files/foo.xlsx, type=seep, expected_flowrate=28.6, element_type=tri6, target_size=2.0, tolerance=0.01 -->
    <!-- test: file=files/foo.xlsx, type=seep_elements, expected_flowrate=40.062, target_size=1.5, tolerance=0.05 -->
    <!-- test: file=files/foo.xlsx, type=fem_elements, expected_fs=1.36, target_size=3.5, tolerance=0.04, f_min=1.0, f_max=1.8, max_iter=4000, benchmark=SSRM-elements -->

The seep_elements / fem_elements types solve ONE problem with every supported
element type (seep: tri3/tri6/quad4/quad8/quad9; FEM: the quadratic tri6/quad8/
quad9) and check each converges and matches the expected flowrate/FS within
tolerance — coverage that a solver/mesh change hasn't broken an element type.
element_types= overrides the default set.

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
    python run_tests.py --roundtrip  # run only the Excel save/load round-trip tests
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


# === Excel round-trip regression (save_slope_data_to_xlsx) ===
# load -> save into a blank template -> reload must reproduce every input
# category. The file set spans circular & non-circular surfaces, profile &
# polygon geometry, reinforcement, piles, distributed loads (both sets), a
# second piezo line, reliability sigmas, and seepage BCs (both sets).
ROUNDTRIP_TEMPLATE = 'docs/inputs/input_template.xlsx'   # editable master (docs link)
BUNDLED_TEMPLATE = 'xslope/resources/input_template.xlsx'  # copy shipped in the wheel
# The /xslope skill body: docs master + the copy shipped in the wheel (used by the
# Studio assistant on pip installs, where docs/ is absent). Must stay byte-identical.
SKILL_MASTER = 'docs/usage/claude/xslope.md'
BUNDLED_SKILL = 'xslope/resources/xslope_skill.md'
ROUNDTRIP_FILES = [
    'docs/inputs/slope/xslope_simple1.xlsx',
    'docs/inputs/slope/xslope_dam.xlsx',
    'docs/inputs/slope/xslope_rapid.xlsx',
    'docs/inputs/slope/xslope_reliability.xlsx',
    'docs/inputs/slope/xslope_rface.xlsx',
    'docs/fem/files/xslope_reinforce_fem.xlsx',
    'docs/fem/files/xslope_piles_fem.xlsx',
    'docs/fem/files/xslope_noncircular_fem.xlsx',
    'docs/fem/files/xslope_griffiths1_load.xlsx',
    'docs/seep/files/xslope_earth_dam1.xlsx',
    'docs/inputs/seep/xslope_earth_dam_bc2.xlsx',
    'docs/seep/files/xslope_levee_poly.xlsx',
    'docs/inputs/seep/xslope_lost_lake.xlsx',
]
# Structured-DXF round-trip files (export_dxf -> read_dxf_layers -> default wizard
# mapping -> build_from_dxf_mapping). Each entry is (file, kind) where kind drives
# the per-feature assertion: 'profile', 'polygon', or 'reinforce'. Needs ezdxf +
# PySide6 (build_from_dxf_mapping lives in studio); skipped when either is absent.
DXF_FILES = [
    ('docs/inputs/slope/xslope_simple1.xlsx', 'profile'),
    ('docs/inputs/slope/xslope_reinf.xlsx', 'reinforce'),
    ('docs/seep/files/xslope_levee_poly.xlsx', 'polygon'),
]
# Source (non-derived) keys that must survive a round-trip. Derived geometry
# (ground_surface, domain_polygon, and polygons-from-profile) is recomputed by
# the loader and so is compared separately (polygon geometry only).
ROUNDTRIP_KEYS = [
    'gamma_water', 'tcrack_depth', 'tcrack_water', 'k_seismic', 'max_depth',
    'profile_lines', 'materials', 'piezo_line', 'piezo_line2', 'circles',
    'non_circ', 'dloads', 'dloads2', 'reinforcement_lines', 'pile_lines',
    'seepage_bc', 'seepage_bc2',
]


def _roundtrip_eq(a, b):
    """Scalar equality with float tolerance; strings compared verbatim (so the
    loader's empty-cell 'nan' sentinel matches itself)."""
    import math
    if str(a) == str(b):
        return True
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-6, abs_tol=1e-6)
    except (TypeError, ValueError):
        return False


def _roundtrip_diff(a, b, path=''):
    """Recursively compare two slope_data values; return a list of mismatch paths."""
    out = []
    if isinstance(a, dict):
        if not isinstance(b, dict):
            return [f"{path}: dict vs {type(b).__name__}"]
        for k in a:
            if k not in b:
                out.append(f"{path}.{k}: missing")
            else:
                out += _roundtrip_diff(a[k], b[k], f"{path}.{k}")
        return out
    if isinstance(a, (list, tuple)):
        if not isinstance(b, (list, tuple)) or len(a) != len(b):
            blen = len(b) if hasattr(b, '__len__') else '?'
            return [f"{path}: len {len(a)} vs {blen}"]
        for i, (x, y) in enumerate(zip(a, b)):
            out += _roundtrip_diff(x, y, f"{path}[{i}]")
        return out
    if not _roundtrip_eq(a, b):
        out.append(f"{path}: {a!r} != {b!r}")
    return out


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
        for key in ['expected_fs', 'expected_flowrate', 'expected_beta', 'tolerance', 'target_size', 'f_min', 'f_max', 'beta']:
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
                t = dict(shared)
                t['method'] = k[3:]      # fs_<method> tag -> solver function name
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
    # `rapid=true` runs the 3-stage rapid-drawdown analysis (needs the seep solution
    # + dloads/seep-bc stage-2 data in the input). Works with any test type.
    rapid = str(test.get('rapid', 'false')).strip().lower() in ('true', '1', 'yes')

    slope_data = load_slope_data(file_path)

    if test_type == 'single_circle':
        circle = slope_data['circles'][0]
        success, result = generate_slices(slope_data, circle=circle, num_slices=num_slices)
        if not success:
            return None, f"generate_slices failed: {result}"
        slice_df, failure_surface = result
        solver_result = solve_selected(method, slice_df, rapid=rapid)
        if isinstance(solver_result, str):
            return None, f"solve failed: {solver_result}"
        # Regression guard: rapid drawdown must write the winning stage's stresses
        # back to the caller's slice_df (else plots get n_eff=0 — see
        # advanced.rapid_drawdown). The FS can be correct while n_eff is stale.
        if rapid and not bool((slice_df['n_eff'].abs() > 0).any()):
            return None, "rapid drawdown left n_eff all-zero in the slice_df (not written back)"
        return solver_result['FS'], None

    elif test_type == 'circular_search':
        fs_cache, converged, search_path, circle_cache = circular_search(
            slope_data, method, num_slices=num_slices, rapid=rapid
        )
        if not fs_cache or fs_cache[0]['FS'] >= 9999:
            return None, "circular_search found no valid surface"
        if rapid:
            crit = fs_cache[0].get('slices')
            if crit is not None and not bool((crit['n_eff'].abs() > 0).any()):
                return None, "rapid search: critical-surface n_eff all-zero (not written back)"
        return fs_cache[0]['FS'], None

    elif test_type == 'noncircular_search':
        fs_cache, converged, search_path = noncircular_search(
            slope_data, method, num_slices=num_slices, rapid=rapid
        )
        if not fs_cache or fs_cache[0]['FS'] >= 9999:
            return None, "noncircular_search found no valid surface"
        return fs_cache[0]['FS'], None

    else:
        return None, f"Unknown LEM test type: {test_type}"


def run_design_test(test):
    """Run a design-sweep regression: edit a profile point to set the slope-face
    angle, rebuild the polygon geometry, then run a circular search and return FS.

    This guards the bug found in the design driver (main_design.py), where editing
    profile_lines and rebuilding only the ground surface left the slice weights —
    and therefore the factor of safety — pinned to the ORIGINAL geometry. Slice
    weights are computed from slope_data['polygons'], so the polygons must be
    regenerated from the edited profile. The expected FS reflects the NEW slope
    angle; if the polygon resync regresses, the search returns the stale
    base-geometry FS and this test fails.
    """
    import math
    from shapely.geometry import Polygon
    from xslope.fileio import load_slope_data, build_ground_surface_from_polygons
    from xslope.mesh import build_polygons
    from xslope.search import circular_search

    file_path = test['file']
    method = test.get('method', 'bishop')
    beta = float(test['beta'])
    toe_index = int(test.get('toe_index', 1))
    slope_index = int(test.get('slope_index', 2))
    num_slices = test.get('num_slices', 30)

    slope_data = load_slope_data(file_path)

    # Move the slope-top point so the slope face makes angle `beta`:
    #   tan(beta) = (y_top - y_toe) / (x_top - x_toe)
    profile = slope_data['profile_lines'][0]['coords']
    x_toe, y_toe = profile[toe_index]
    _x_top, y_top = profile[slope_index]
    x_top_new = x_toe + (y_top - y_toe) / math.tan(math.radians(beta))
    slope_data['profile_lines'][0]['coords'][slope_index] = (x_top_new, y_top)

    # Rebuild material polygons from the edited profile, then re-derive the ground
    # surface / domain polygon from them — the resync that the original bug omitted.
    polys = [
        {'polygon': Polygon(p['coords']), 'mat_id': p['mat_id']}
        for p in build_polygons(slope_data={'profile_lines': slope_data['profile_lines'],
                                             'max_depth': slope_data.get('max_depth')})
    ]
    slope_data['polygons'] = polys
    ground_surface, domain_polygon = build_ground_surface_from_polygons(polys)
    slope_data['ground_surface'] = ground_surface
    slope_data['domain_polygon'] = domain_polygon

    fs_cache, converged, search_path, circle_cache = circular_search(
        slope_data, method, num_slices=num_slices)
    if not fs_cache or fs_cache[0]['FS'] >= 9999:
        return None, "design_search found no valid surface"
    return fs_cache[0]['FS'], None


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


def run_seep_elements_test(test):
    """Element-type coverage: solve ONE seepage problem with every supported
    element type (tri3, tri6, quad4, quad8, quad9) and check each converges and
    reproduces the expected flowrate within tolerance. The same physical problem
    must give the same flowrate on every element type, so a single expected value
    covers all — a change to the seepage solver or mesh that breaks (or shifts) an
    element type is caught. Returns (0.0, None) if all pass, else (None, message
    naming the failing types). element_types= overrides the default set."""
    import io
    import contextlib
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import build_seep_data, run_seepage_analysis

    slope_data = load_slope_data(test['file'])
    polygons = get_material_polygons(slope_data)
    element_types = [s.strip() for s in
                     test.get('element_types', 'tri3,tri6,quad4,quad8,quad9').split(',')]
    target_size = test.get('target_size')
    if target_size is None:
        x_coords = [x for x, _ in slope_data['ground_surface'].coords]
        target_size = (max(x_coords) - min(x_coords)) / 120
    expected = test['expected_flowrate']
    tol = test.get('tolerance', 0.05) * abs(expected)

    problems = []
    for et in element_types:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mesh = build_mesh_from_polygons(polygons, target_size, et)
                seep_data = build_seep_data(mesh, slope_data)
                solution = run_seepage_analysis(seep_data, tol=1e-4)
        except Exception as e:
            problems.append(f"{et}: {type(e).__name__}: {e}")
            continue
        if not solution.get('converged', True):
            problems.append(f"{et}: did not converge")
        elif abs(solution['flowrate'] - expected) > tol:
            problems.append(f"{et}: flowrate {solution['flowrate']:.3f} vs {expected:.3f}")
    if problems:
        return None, "; ".join(problems)
    return 0.0, None


def run_fem_elements_test(test):
    """Element-type coverage: solve ONE FEM (SSRM) problem with every QUADRATIC
    element type (tri6, quad8, quad9) and check each converges to the expected
    factor of safety within tolerance. FEM uses quadratic elements only (linear
    tri3/quad4 lock and overestimate FS). Returns (0.0, None) if all pass, else
    (None, message naming the failing types). element_types= overrides the set."""
    import io
    import contextlib
    from xslope.fileio import load_slope_data
    from xslope.fem import build_fem_data, solve_ssrm
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry)

    slope_data = load_slope_data(test['file'])
    element_types = [s.strip() for s in
                     test.get('element_types', 'tri6,quad8,quad9').split(',')]
    target_size = test.get('target_size')
    if target_size is None:
        x_coords = [x for x, _ in slope_data['ground_surface'].coords]
        target_size = (max(x_coords) - min(x_coords)) / 100
    expected = test['expected_fs']
    tol = test.get('tolerance', 0.03)
    kwargs = {}
    if 'max_iter' in test:
        kwargs['max_iterations'] = int(test['max_iter'])
    # ssrm_tol is the bisection precision (how tightly FS is bracketed); keep it
    # well below the cross-element comparison tolerance.
    ssrm_tol = float(test.get('ssrm_tol', 0.01))

    constraint_lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
    polygons = get_material_polygons(slope_data, reinf_lines=constraint_lines)
    problems = []
    for et in element_types:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                mesh = build_mesh_from_polygons(
                    polygons, target_size=target_size, element_type=et,
                    lines=constraint_lines)
                fem_data = build_fem_data(slope_data, mesh)
                result = solve_ssrm(fem_data, F_min=test.get('f_min', 0.5),
                                    F_max=test.get('f_max', 3.0), tolerance=ssrm_tol,
                                    debug_level=0, **kwargs)
        except Exception as e:
            problems.append(f"{et}: {type(e).__name__}: {e}")
            continue
        if not result.get('converged', False):
            problems.append(f"{et}: SSRM failed ({result.get('error', '?')})")
        elif abs(result['FS'] - expected) > tol:
            problems.append(f"{et}: FS {result['FS']:.3f} vs {expected:.3f}")
    if problems:
        return None, "; ".join(problems)
    return 0.0, None


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


def run_fem_reliability_test(test):
    """FEM reliability regression: run reliability_fem on a mesh and return the
    lognormal reliability index beta_ln (compared to expected_beta). Guards the
    whole TSPM-over-SSRM pipeline — F_MLV, the F+/F- perturbation solves, and the
    beta combination — on a coarse mesh so it stays affordable. ssrm_tol sets the
    bisection precision; the tag's `tolerance` is the beta comparison tolerance."""
    import io
    import contextlib
    from xslope.fileio import load_slope_data
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry)
    from xslope.advanced import reliability_fem

    slope_data = load_slope_data(test['file'])
    element_type = test.get('element_type', 'tri6')
    target_size = test.get('target_size')
    if target_size is None:
        x_coords = [x for x, _ in slope_data['ground_surface'].coords]
        target_size = (max(x_coords) - min(x_coords)) / 100
    constraint_lines, _n_reinf, _n_pile = extract_constraint_line_geometry(slope_data)
    polygons = get_material_polygons(slope_data, reinf_lines=constraint_lines)
    mesh = build_mesh_from_polygons(polygons, target_size=target_size,
                                    element_type=element_type, lines=constraint_lines)
    with contextlib.redirect_stdout(io.StringIO()):
        success, result = reliability_fem(
            slope_data, mesh=mesh, F_min=test.get('f_min', 0.5),
            F_max=test.get('f_max', 3.0), tolerance=float(test.get('ssrm_tol', 0.001)),
            failure_criterion=test.get('criterion', 'non_convergence'))
    if not success:
        return None, f"reliability_fem failed: {result}"
    return result['beta_ln'], None


def run_roundtrip_test(test):
    """Verify save_slope_data_to_xlsx round-trips a file: load -> save into a
    blank template -> reload must reproduce every input category.

    Returns (0.0, None) on success, or (None, message) listing the mismatches.
    """
    import tempfile
    from xslope.fileio import load_slope_data, save_slope_data_to_xlsx

    d1 = load_slope_data(test['file'])
    tmp = tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False).name
    try:
        save_slope_data_to_xlsx(d1, tmp, template=test['template'])
        d2 = load_slope_data(tmp)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)

    mismatches = []
    for k in ROUNDTRIP_KEYS:
        mismatches += _roundtrip_diff(d1.get(k), d2.get(k), k)

    # Polygon-sourced geometry has no profile_lines to compare, so check the
    # reconstructed polygon zones directly (geometry + material assignment).
    if not d1.get('profile_lines'):
        p1 = d1.get('polygons') or []
        p2 = d2.get('polygons') or []
        if len(p1) != len(p2):
            mismatches.append(f"polygons: len {len(p1)} vs {len(p2)}")
        else:
            for i, (a, b) in enumerate(zip(p1, p2)):
                if a.get('mat_id') != b.get('mat_id'):
                    mismatches.append(f"polygons[{i}].mat_id")
                if not a['polygon'].equals(b['polygon']):
                    mismatches.append(f"polygons[{i}].geom")

    if mismatches:
        return None, "round-trip mismatch: " + "; ".join(mismatches[:5])
    return 0.0, None


def run_template_sync_test(test):
    """Verify a copy shipped in the wheel (xslope/resources) is byte-identical to
    its editable docs master, so the two can't silently drift when the master is
    tweaked. test['master']/test['copy'] name the pair (defaults: input template).

    Returns (0.0, None) if identical, else (None, message).
    """
    import filecmp
    master = test.get('master', ROUNDTRIP_TEMPLATE)
    copy = test.get('copy', BUNDLED_TEMPLATE)
    if not os.path.exists(master):
        return None, f"master file missing: {master}"
    if not os.path.exists(copy):
        return None, f"packaged copy missing: {copy}"
    if not filecmp.cmp(master, copy, shallow=False):
        return None, (f"packaged copy differs from master — run: "
                      f"cp {master} {copy}")
    return 0.0, None


def _default_dxf_mapping(layers):
    """The per-layer mapping the import wizard seeds by default: target from
    ``suggest_dxf_target`` (xslope's own export layer names + geometry kind), and
    the material name from the layer name (PROFILE_ stripped). Mirrors
    studio.dialogs.DxfImportDialog so the suite exercises the same defaults."""
    from xslope.cad import suggest_dxf_target
    out = {}
    for name, geom in layers.items():
        target = suggest_dxf_target(name, geom)
        mat = (name[8:] if name.upper().startswith('PROFILE_') else name)
        out[name] = {'target': target,
                     'material': mat if target in ('material_zone', 'profile') else None}
    return out


def run_dxf_roundtrip_test(test):
    """Round-trip a model through the structured DXF export/import path: load ->
    export_dxf -> read_dxf_layers -> default wizard mapping -> build_from_dxf_mapping
    must recover the geometry (feature counts + circles). Returns (0.0, None) or
    (None, message)."""
    import tempfile
    # Ensure a (headless) QApplication exists — ProjectDocument is a QObject.
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])
    from xslope.fileio import load_slope_data
    from xslope.cad import export_dxf, read_dxf_layers
    from studio.document import ProjectDocument

    sd = load_slope_data(test['file'])
    tmp = tempfile.NamedTemporaryFile(suffix='.dxf', delete=False).name
    try:
        export_dxf(sd, tmp)
        layers, _ = read_dxf_layers(tmp)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)

    mapping = _default_dxf_mapping(layers)
    doc = ProjectDocument()
    doc.build_from_dxf_mapping(layers, mapping)
    out = doc.slope_data

    problems = []
    kind = test.get('kind')
    if kind == 'profile':
        n0, n1 = len(sd.get('profile_lines') or []), len(out.get('profile_lines') or [])
        if n1 != n0:
            problems.append(f"profile_lines {n1} vs {n0}")
        if out.get('ground_surface') is None:
            problems.append("ground_surface missing")
    elif kind == 'polygon':
        n0, n1 = len(sd.get('polygons') or []), len(out.get('polygons') or [])
        if n1 != n0:
            problems.append(f"polygons {n1} vs {n0}")
        if not out.get('materials'):
            problems.append("no materials")
    elif kind == 'reinforce':
        n0 = len(sd.get('reinforcement_lines') or sd.get('reinforce_lines') or [])
        n1 = len(out.get('reinforcement_lines') or [])
        if n1 != n0:                          # rejoin must not split into per-segment lines
            problems.append(f"reinforcement {n1} vs {n0}")
    # Circles are recovered (radius fit from the arc) wherever present.
    if sd.get('circles'):
        n0, n1 = len(sd['circles']), len(out.get('circles') or [])
        if n1 != n0:
            problems.append(f"circles {n1} vs {n0}")

    if problems:
        return None, "DXF round-trip: " + "; ".join(problems[:5])
    return 0.0, None


def run_vg_kr_test(test):
    """Unit check for the van Genuchten relative-permeability function and the
    kr-model dispatch (xslope.seep). Verifies kr_vg_vec against an independent
    scalar evaluation of the Mualem-vG formula across pressure heads and (alpha, n),
    the saturation / floor / monotonicity behavior, and that the dispatcher reduces
    EXACTLY to the linear-front kr when no vG model is active. Returns (0.0, None)
    on success, else (None, message). No mesh / gmsh needed."""
    import numpy as np
    from xslope.seep import (kr_vg_vec, kr_relative_vec, kr_relative,
                             kr_frontal_vec, KR_LF, KR_VG)
    KR_MIN = 1e-4

    def vg_scalar(p, a, n):
        # Independent scalar reference for Mualem-van Genuchten kr.
        n = max(n, 1.0 + 1e-6)
        m = 1.0 - 1.0 / n
        if p >= 0.0:
            return 1.0
        se = (1.0 + (a * abs(p)) ** n) ** (-m)
        kr = (se ** 0.5) * (1.0 - (1.0 - se ** (1.0 / m)) ** m) ** 2
        return min(max(kr, KR_MIN), 1.0)

    problems = []
    heads = np.array([-200., -50., -10., -3., -1., -0.3, -0.05, 0.0, 2.0])
    for a, n in [(0.075, 1.89), (2.286, 1.89), (1.097, 1.56), (0.5, 1.1), (3.0, 4.0)]:
        got = kr_vg_vec(heads, a, n)
        ref = np.array([vg_scalar(float(p), a, n) for p in heads])
        if not np.allclose(got, ref, atol=1e-12, rtol=1e-9):
            problems.append(f"kr_vg_vec != Mualem-vG for (a={a}, n={n}): max|d|={np.max(np.abs(got-ref)):.2e}")
        if got[7] != 1.0 or got[8] != 1.0:
            problems.append(f"kr not 1.0 at saturation for (a={a}, n={n})")
        unsat = got[:7]
        if np.any(unsat < KR_MIN - 1e-15) or np.any(unsat > 1.0 + 1e-12):
            problems.append(f"kr out of [{KR_MIN},1] for (a={a}, n={n})")
        if np.any(np.diff(unsat) < -1e-12):     # drier (more negative p) -> lower kr
            problems.append(f"kr not monotonic in p for (a={a}, n={n})")

    # Dispatch must reduce EXACTLY to the linear front when no vG model is active.
    P = np.array([[-2.0, -1.0, 0.5], [-5.0, -0.2, 1.0]])
    kr0 = np.array([0.001, 0.01]); h0 = np.array([-1.0, -2.0])
    base = kr_frontal_vec(P, kr0[:, None], h0[:, None])
    if not np.array_equal(kr_relative_vec(P, kr0[:, None], h0[:, None]), base):
        problems.append("kr_relative_vec(model=None) != kr_frontal_vec (lf not bit-identical)")
    alllf = kr_relative_vec(P, kr0[:, None], h0[:, None],
                            vg_a=np.zeros((2, 1)), vg_n=np.zeros((2, 1)),
                            model=np.array([[KR_LF], [KR_LF]]))
    if not np.array_equal(alllf, base):
        problems.append("kr_relative_vec(all-lf) != kr_frontal_vec")
    mixed = kr_relative_vec(P, kr0[:, None], h0[:, None],
                            vg_a=np.array([[0.075], [0.0]]), vg_n=np.array([[1.89], [0.0]]),
                            model=np.array([[KR_VG], [KR_LF]]))
    if not np.array_equal(mixed[1], base[1]):
        problems.append("mixed dispatch perturbed the lf row")
    if not np.allclose(mixed[0], kr_vg_vec(P[0], 0.075, 1.89)):
        problems.append("mixed dispatch wrong on the vG row")

    # Scalar kr_relative agrees with the vector path and dispatches on model.
    if abs(kr_relative(-1.5, 0.01, -1.0, model=KR_LF) - float(kr_frontal_vec(np.array(-1.5), 0.01, -1.0))) > 1e-12:
        problems.append("scalar kr_relative(lf) != kr_frontal")
    if abs(kr_relative(-1.5, 0.0, 0.0, vg_a=0.075, vg_n=1.89, model=KR_VG) - vg_scalar(-1.5, 0.075, 1.89)) > 1e-12:
        problems.append("scalar kr_relative(vg) != Mualem-vG")

    if problems:
        return None, "; ".join(problems[:5])
    return 0.0, None


def run_mesh_conform_test(test):
    """Guard the conforming-edge preprocessing in build_mesh_from_polygons: a vertex
    in the interior of a neighbour's shared edge (a T-junction) must be inserted so
    the interface meshes without a slit. Builds the earth-dam shell/core case with a
    T-junction at (55,18) on the core-top interface and checks (a) the preprocessor
    inserts it into the core edge and (b) every interface node is shared by elements
    of BOTH materials. Returns (0.0, None) on success, else (None, message)."""
    import numpy as np
    from collections import defaultdict
    from xslope.mesh import build_mesh_from_polygons, make_polygons_conforming

    core = [(46.5, 0.0), (63.5, 0.0), (59.0, 18.0), (51.0, 18.0)]
    shell = [(0.0, 0.0), (51.0, 22.0), (59.0, 22.0), (110.0, 0.0), (63.5, 0.0),
             (59.0, 18.0), (55.0, 18.0), (51.0, 18.0), (46.5, 0.0)]   # T-junction (55,18)

    problems = []
    pc = [list(shell), list(core)]
    make_polygons_conforming(pc)
    if (55.0, 18.0) not in [tuple(p) for p in pc[1]]:
        problems.append("preprocessor did not insert the T-junction vertex (55,18) into the core edge")

    mesh = build_mesh_from_polygons([{'coords': shell, 'mat_id': 0},
                                     {'coords': core, 'mat_id': 1}], 3.0, 'tri3')
    nodes = np.asarray(mesh['nodes']); el = np.asarray(mesh['elements'])
    et = np.asarray(mesh['element_types']); em = np.asarray(mesh['element_materials'])
    node_mats = defaultdict(set)
    for ei in range(len(el)):
        for nd in el[ei][:int(et[ei])]:
            node_mats[int(nd)].add(int(em[ei]))
    nonconf = sum(1 for i, (x, y) in enumerate(nodes)
                  if abs(y - 18.0) < 1e-6 and 51.0 + 1e-6 < x < 59.0 - 1e-6
                  and len(node_mats.get(i, set())) < 2)
    if nonconf:
        problems.append(f"{nonconf} non-conforming interface node(s) on the core-top edge")

    if problems:
        return None, "; ".join(problems[:5])
    return 0.0, None


def run_test(test):
    """Run a single test and return (computed_value, error_msg)."""
    test_type = test.get('type', '')
    if test_type == 'mesh_conform':
        return run_mesh_conform_test(test)
    if test_type == 'vg_kr':
        return run_vg_kr_test(test)
    if test_type == 'roundtrip':
        return run_roundtrip_test(test)
    if test_type == 'dxf':
        return run_dxf_roundtrip_test(test)
    if test_type == 'template_sync':
        return run_template_sync_test(test)
    if test_type == 'fem_ssrm':
        return run_fem_test(test)
    if test_type == 'seep_elements':
        return run_seep_elements_test(test)
    if test_type == 'fem_elements':
        return run_fem_elements_test(test)
    if test_type == 'fem_reliability':
        return run_fem_reliability_test(test)
    elif test_type == 'seep':
        return run_seep_test(test)
    elif test_type == 'reliability':
        return run_reliability_test(test)
    elif test_type == 'design_search':
        return run_design_test(test)
    else:
        return run_lem_test(test)


def main():
    parser = argparse.ArgumentParser(description='xslope regression test suite')
    parser.add_argument('--lem', action='store_true', help='Run only LEM tests')
    parser.add_argument('--fem', action='store_true', help='Run only FEM tests')
    parser.add_argument('--seep', action='store_true', help='Run only seepage tests')
    parser.add_argument('--roundtrip', action='store_true',
                        help='Run only the Excel save/load round-trip tests')
    parser.add_argument('--dxf', action='store_true',
                        help='Run only the structured DXF export/import round-trip '
                             'tests (needs ezdxf + PySide6)')
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
    run_all = not (args.lem or args.fem or args.seep or args.roundtrip or args.dxf)
    run_lem = args.lem or run_all
    run_fem = args.fem or run_all
    run_seep = args.seep or run_all
    run_roundtrip = args.roundtrip or run_all
    run_dxf = args.dxf or run_all

    # Discover tests from markdown files
    tests = []
    if run_lem:
        lem_samples = Path('docs/lem/samples.md')
        if lem_samples.exists():
            tests.extend(parse_test_tags(lem_samples))
        # The design example carries a design_search regression (guards the
        # profile-edit / polygon-rebuild bug found in main_design.py).
        lem_design = Path('docs/lem/design.md')
        if lem_design.exists():
            tests.extend(parse_test_tags(lem_design))
    if run_fem:
        fem_samples = Path('docs/fem/samples.md')
        if fem_samples.exists():
            tests.extend(parse_test_tags(fem_samples))
    if run_seep:
        # Fast, file-less unit check for the van Genuchten kr function + dispatch.
        tests.append({'type': 'vg_kr', 'file': 'kr_vg_vec (unit)', 'method': '-',
                      'source': 'vg_kr'})
        # Mesh conforming-edge (T-junction) regression.
        tests.append({'type': 'mesh_conform', 'file': 'conforming edges (mesh)',
                      'method': '-', 'source': 'mesh_conform'})
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
    private_dir = os.environ.get('XSLOPE_PRIVATE_TESTS')
    if not private_dir:
        siblings = Path(__file__).resolve().parent.parent
        # 'xslope_private_tests' is the pre-rename name, still used by old clones.
        for name in ('xslope_private', 'xslope_private_tests'):
            if (siblings / name).is_dir():
                private_dir = str(siblings / name)
                break
    private_path = Path(private_dir or '')
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

    # Excel round-trip tests (save_slope_data_to_xlsx). Built from a curated file
    # list rather than markdown tags, since they check load/save fidelity, not FS.
    if run_roundtrip:
        # Guard that the packaged copies (shipped in the wheel) haven't drifted
        # from their editable docs masters.
        tests.append({'type': 'template_sync', 'file': BUNDLED_TEMPLATE,
                      'method': '-', 'source': 'template'})
        tests.append({'type': 'template_sync', 'file': BUNDLED_SKILL,
                      'master': SKILL_MASTER, 'copy': BUNDLED_SKILL,
                      'method': '-', 'source': 'skill'})
        if Path(ROUNDTRIP_TEMPLATE).exists():
            n_rt = 0
            for fp in ROUNDTRIP_FILES:
                if Path(fp).exists():
                    tests.append({'type': 'roundtrip', 'file': fp,
                                  'template': ROUNDTRIP_TEMPLATE, 'method': '-',
                                  'source': 'roundtrip'})
                    n_rt += 1
            if n_rt and not run_all:
                print(f"Including {n_rt} Excel round-trip tests")
        else:
            print(f"Skipping round-trip tests (template not found: {ROUNDTRIP_TEMPLATE})")

    # Structured DXF export/import round-trip tests. These touch the studio layer
    # (build_from_dxf_mapping) and ezdxf, so they're skipped cleanly when either
    # dependency is absent — the engine-only suite stays runnable without the GUI.
    if run_dxf:
        try:
            import ezdxf            # noqa: F401
            import PySide6          # noqa: F401
            dxf_ok = True
        except ImportError:
            dxf_ok = False
            print("Skipping DXF round-trip tests (ezdxf/PySide6 not installed)")
        if dxf_ok:
            n_dxf = 0
            for fp, kind in DXF_FILES:
                if Path(fp).exists():
                    tests.append({'type': 'dxf', 'file': fp, 'kind': kind,
                                  'method': '-', 'source': 'dxf'})
                    n_dxf += 1
            if n_dxf and not run_all:
                print(f"Including {n_dxf} DXF round-trip tests")

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
        elif test_type in ('reliability', 'fem_reliability'):
            expected = test.get('expected_beta')
            tol = test.get('tolerance', 0.02)
            label = 'beta'
        elif test_type in ('roundtrip', 'template_sync', 'dxf', 'vg_kr',
                            'mesh_conform', 'seep_elements', 'fem_elements'):
            expected = 0.0          # these return 0.0 on success (pass/fail tests)
            tol = 0.0
            label = 'mismatch'
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
