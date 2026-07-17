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
# M-P with f(x)==1 must reproduce Spencer exactly, on both slope facings (the S3b gate).
AXIAL_MIRROR_LEFT = 'docs/inputs/slope/xslope_nail_axial.xlsx'
AXIAL_MIRROR_RIGHT = 'docs/inputs/slope/xslope_nail_axial_rface.xlsx'
MP_SPENCER_LEFT = 'docs/inputs/slope/xslope_simple1.xlsx'
MP_SPENCER_RIGHT = 'docs/inputs/slope/xslope_rface.xlsx'
# Stage-1 FS 1.91 unweakened; halving c/phi/d/psi drops it below 1, which must be refused.
DRAWDOWN_GUARD_FILE = 'docs/inputs/slope/xslope_rapid.xlsx'
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
    'line_loads', 'seepage_bc', 'seepage_bc2',
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
        if 'file2' in params:
            params['file2'] = str(md_dir / params['file2'])

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
    # `composite=true` lets a circle be truncated at the bottom of the model and run
    # along it (slice.CompositeSurface). Only for problems whose base is a real
    # impenetrable boundary; off by default.
    composite = str(test.get('composite', 'false')).strip().lower() in ('true', '1', 'yes')
    # `seed=grid` runs the circular search from an automatic grid-and-tangent sweep
    # instead of (only) the circles sheet — the global-search mode.
    seed = str(test.get('seed', 'circles')).strip().lower()

    slope_data = load_slope_data(file_path)

    if test_type == 'single_circle':
        circle = slope_data['circles'][0]
        success, result = generate_slices(slope_data, circle=circle, num_slices=num_slices,
                                          composite=composite)
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
            slope_data, method, num_slices=num_slices, rapid=rapid, composite=composite,
            seed=seed
        )
        if not fs_cache or fs_cache[0]['FS'] >= 9999:
            return None, "circular_search found no valid surface"
        if rapid:
            crit = fs_cache[0].get('slices')
            if crit is not None and not bool((crit['n_eff'].abs() > 0).any()):
                return None, "rapid search: critical-surface n_eff all-zero (not written back)"
        return fs_cache[0]['FS'], None

    elif test_type == 'single_noncirc':
        # Evaluate the file's specified non-circular surface as-is (no search) —
        # the "predefined slip surface" form of the verification problems.
        success, result = generate_slices(slope_data, non_circ=slope_data['non_circ'],
                                          num_slices=num_slices)
        if not success:
            return None, f"generate_slices failed: {result}"
        slice_df, failure_surface = result
        solver_result = solve_selected(method, slice_df, rapid=rapid)
        if isinstance(solver_result, str):
            return None, f"solve failed: {solver_result}"
        return solver_result['FS'], None

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
    from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                             extract_constraint_line_geometry, extract_point_constraints)

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
            lines=constraint_lines,
            point_constraints=extract_point_constraints(slope_data)
        )

    fem_data = build_fem_data(slope_data, mesh)
    f_min = test.get('f_min', 0.5)
    f_max = test.get('f_max', 3.0)
    kwargs = {}
    if 'criterion' in test:
        kwargs['failure_criterion'] = test['criterion']
    if 'max_iter' in test:
        kwargs['max_iterations'] = int(test['max_iter'])
    # TEMP (#50): the Dawson per-node criterion needs ~3x the iterations the old rate-based
    # test did (displacement settles ~5k, the per-node max only ~11k). Let the experiment
    # raise the floor without editing 64 tags.
    import os as _os
    _floor = int(_os.environ.get('XSLOPE_MIN_MAX_ITER', '0'))
    if _floor:
        kwargs['max_iterations'] = max(kwargs.get('max_iterations', 0), _floor)
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
    solution = run_seepage_analysis(seep_data, tol=1e-4,
                                    max_iter=int(test.get('max_iter', 400)))

    if not solution.get('converged', True):
        return None, "seepage solution did not converge (flowrate unreliable)"
    return solution['flowrate'], None


def run_seep_head_test(test):
    """Check solved total head at named points (groundwater-corpus tags).

    Tag keys: file, points="x:y:h;x:y:h;...", tolerance (m, absolute),
    optional target_size / element_type. Meshes and solves live like
    run_seep_test, then interpolates head at each point by inverse-distance
    over the four nearest nodes. Pass/fail: returns 0.0 on success."""
    import numpy as np
    from xslope.fileio import load_slope_data
    from xslope.mesh import get_material_polygons, build_mesh_from_polygons
    from xslope.seep import build_seep_data, run_seepage_analysis

    slope_data = load_slope_data(test['file'])
    polygons = get_material_polygons(slope_data)
    target_size = test.get('target_size')
    if target_size is None:
        xs = [x for x, _ in slope_data['ground_surface'].coords]
        target_size = (max(xs) - min(xs)) / 120
    mesh = build_mesh_from_polygons(polygons, float(target_size),
                                    test.get('element_type', 'tri3'))
    seep_data = build_seep_data(mesh, slope_data)
    solution = run_seepage_analysis(seep_data, tol=1e-5,
                                    max_iter=int(test.get('max_iter', 400)))
    if not solution.get('converged', True):
        return None, "seepage solution did not converge"

    nodes = seep_data['nodes']
    h = np.asarray(solution['head'])

    def head_at(xq, yq):
        d2 = (nodes[:, 0] - xq) ** 2 + (nodes[:, 1] - yq) ** 2
        idx = np.argsort(d2)[:4]
        w = 1.0 / np.maximum(d2[idx], 1e-12)
        return float(np.sum(w * h[idx]) / np.sum(w))

    tol = float(test.get('tolerance', 0.01))
    errs = []
    for triplet in str(test['points']).split(';'):
        xs_, ys_, hs_ = (float(v) for v in triplet.split(':'))
        got = head_at(xs_, ys_)
        if abs(got - hs_) > tol:
            errs.append(f"({xs_:g},{ys_:g}): expected {hs_:.3f}, got {got:.3f}")
    if errs:
        return None, "head mismatch: " + "; ".join(errs)
    return 0.0, None


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


def run_sensitivity_test(test):
    """Run a sensitivity-sweep regression (docs/lem/sensitivity.md tags).

    Tag keys: file, param, method, num_slices, n, rel_range, fs_base, fs_low,
    fs_high, tolerance. Runs sensitivity() with a searched sweep and checks
    the base row and the two range endpoints. Returns (fs_base, None) on pass
    so the table shows the base FS; a mismatch anywhere returns an error."""
    import io, contextlib
    from xslope.fileio import load_slope_data
    from xslope.sensitivity import sensitivity

    slope_data = load_slope_data(test['file'])
    tol = float(test.get('tolerance', 0.01))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ok, res = sensitivity(slope_data, param=test['param'],
                              rel_range=float(test.get('rel_range', 0.5)),
                              n=int(test.get('n', 3)),
                              methods=(test['method'],),
                              search=str(test.get('search', 'true')).lower()
                                     not in ('false', '0', 'no'),
                              num_slices=int(test.get('num_slices', 30)))
    if not ok:
        return None, f"sensitivity failed: {res}"
    df = res['df']
    swept = df.loc[~df['is_base']].sort_values('value')
    checks = [('expected_base', float(df.loc[df['is_base'], 'fs'].iloc[0])),
              ('expected_low', float(swept['fs'].iloc[0])),
              ('expected_high', float(swept['fs'].iloc[-1]))]
    if not bool(swept['success'].all()):
        return None, f"sweep had failed points: {swept.loc[~swept['success'], 'msg'].tolist()}"
    for key, got in checks:
        if key in test and abs(got - float(test[key])) > tol:
            return None, f"{key}: expected {float(test[key]):.3f}, got {got:.3f}"
    return checks[0][1], None


def run_reliability_test(test):
    """Run a single reliability analysis, returning the lognormal reliability index beta."""
    from xslope.fileio import load_slope_data
    from xslope.advanced import reliability as reliability_analysis

    file_path = test['file']
    method = test.get('method', 'spencer')
    circular = str(test.get('circular', 'true')).lower() not in ('false', '0', 'no')

    do_search = str(test.get('search', 'true')).lower() not in ('false', '0', 'no')

    slope_data = load_slope_data(file_path)
    success, result = reliability_analysis(slope_data, method, circular=circular, debug_level=0,
                                           search=do_search)
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


def run_gsat_pair_test(test):
    """The gamma_sat sidecar equivalence gate (plan §1.2 S1): a model hand-zoned
    into moist/saturated polygons at the water table must give the same FS as the
    single-material gamma/gamma_sat + piezo-sidecar formulation. Discretization
    differs (the two force different slice boundaries), so the comparison runs at
    300 slices where both converge to the same answer (verified 5e-8 at build
    time); tolerance 1e-5."""
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope import solve

    vals = []
    for f in (test['file'], test['file2']):
        sd = load_slope_data(f)
        ok, result = generate_slices(sd, circle=sd['circles'][0], num_slices=300)
        if not ok:
            return None, f"slice generation failed for {f}: {result}"
        df, _ = result
        s_ok, r = solve.spencer(df)
        if not s_ok:
            return None, f"spencer failed for {f}: {r}"
        vals.append(r['FS'])
    d = abs(vals[0] - vals[1])
    if d > 1e-5:
        return None, (f"zoned vs sidecar FS differ by {d:.2e} "
                      f"({vals[0]:.6f} vs {vals[1]:.6f})")
    return 0.0, None    # pass/fail test: 0.0 = the two formulations agree


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
        save_slope_data_to_xlsx(d1, tmp, template=test.get('template', ROUNDTRIP_TEMPLATE))
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


# === Editor round-trip guard (studio.editors) ===
# Opening each CategoryEditor on a fully-populated record and applying with NO
# changes must leave its inputs deep-equal — the editor module's contract that "a
# round-trip through an editor never drops data". This is the guard that would
# have caught the materials-combo corruption: a Field whose `choices` lagged the
# loader's accepted set left the combo unselectable for a valid value, so OK
# silently rewrote option/u/unsat to the combo's first item. The Excel --roundtrip
# test never routes through the editors, so it can't see this. Needs PySide6 (the
# studio layer), so it's registered only when that imports — engine-only installs
# skip it cleanly, like the DXF tests.

# slope_data keys each editor round-trips (its "managed" inputs). Keys the editors
# recompute from these (ground_surface, domain_polygon, tcrack_surface,
# reinforce_lines, circular, has_seepage_bc2, ...) are derived, not source, so are
# excluded — only the source inputs must survive a no-op round-trip.
_EDITOR_MANAGED_KEYS = {
    "global": ["gamma_water", "tcrack_depth", "tcrack_water", "k_seismic"],
    "materials": ["materials"],
    "circles": ["circles"],
    "non_circ": ["non_circ"],
    "piezo": ["piezo_line", "piezo_line2"],
    "dloads": ["dloads", "dloads2"],
    "seep_bc": ["seepage_bc", "seepage_bc2"],
    "piles": ["pile_lines"],
    "reinforce": ["reinforcement_lines"],
    "line_loads": ["line_loads"],
    "profile": ["profile_lines"],
    "polygons": ["polygons"],
}


def _editor_full_material(name, option, u, unsat):
    """A material with EVERY loader-produced key set to a distinct non-default
    value, so a dropped key is caught; option/u/unsat carry the enum value under
    test. Together the fixture's rows exercise every accepted option (mc/cp/pow/hb),
    u (none/piezo/seep/ru) and unsat (lf/vg/gard) value."""
    return {
        "name": name, "gamma": 120.0, "gamma_sat": 125.0, "option": option,
        "c": 100.0, "phi": 30.0, "cp": 0.5, "r_elev": 12.0, "d": 3.0, "psi": 5.0,
        "pow_a": 1.1, "pow_b": 0.9, "pow_c": 2.0, "pow_d": 4.0,
        "u": u, "ru": 0.35,
        "sigma_gamma": 1.0, "sigma_c": 2.0, "sigma_phi": 3.0, "sigma_cp": 0.1,
        "sigma_d": 0.2, "sigma_psi": 0.3,
        "k1": 1.0, "k2": 2.0, "alpha": 0.4, "unsat": unsat,
        "kr0": 0.6, "h0": 7.0, "vg_a": 0.05, "vg_n": 1.4, "E": 5000.0, "nu": 0.3,
        "hb_sci": 50.0, "hb_gsi": 60.0, "hb_mi": 10.0, "hb_d": 0.0,
    }


def _editor_fixture():
    """A fully-populated slope_data spanning every editable category, canonical so
    each editor's apply-time transform (circle R/Depth, pile theta, geometry
    resync, ...) is idempotent on its managed keys."""
    import math
    from shapely.geometry import Polygon

    materials = [
        _editor_full_material("m-mc-none-lf",    "mc",  "none",  "lf"),
        _editor_full_material("m-cp-piezo-vg",   "cp",  "piezo", "vg"),
        _editor_full_material("m-pow-seep-gard", "pow", "seep",  "gard"),
        _editor_full_material("m-hb-ru-lf",      "hb",  "ru",    "lf"),
        # A seep-only material with a BLANK strength option ('' — valid per the
        # loader). Locks the MaterialsEditor combo's empty entry: without it the
        # combo would normalize blank -> 'mc' and drop the seep-only classification.
        _editor_full_material("m-blank-seep-vg", "",    "seep",  "vg"),
    ]

    def pile(x1, y1, x2, y2, appl="active"):
        p = {"label": "P", "x1": x1, "y1": y1, "x2": x2, "y2": y2, "H": 10.0,
             "theta_p": math.degrees(math.atan2(x2 - x1, -(y2 - y1))),
             "D_pile": 2.0, "S": 6.0, "E": 3000.0, "I": 1.5, "area": 4.0,
             "V_cap": 50.0, "M_cap": 200.0, "fixity": "fixed", "appl": appl}
        return p

    return {
        "gamma_water": 62.4, "tcrack_depth": 0.0, "tcrack_water": 0.0,
        "k_seismic": 0.15, "max_depth": 0.0,
        "profile_lines": [{"coords": [(0.0, 0.0), (20.0, 20.0), (100.0, 20.0)],
                           "mat_id": 0}],
        "polygons": [{"polygon": Polygon([(0.0, 0.0), (20.0, 20.0), (100.0, 20.0),
                                          (100.0, 0.0)]), "mat_id": 0}],
        "ground_surface": None, "domain_polygon": None, "tcrack_surface": None,
        "materials": materials,
        "piezo_line": [(0.0, 5.0), (100.0, 5.0)],
        "piezo_line2": [(0.0, 3.0), (100.0, 3.0)],
        "circular": True,
        "circles": [{"Xo": 10.0, "Yo": 40.0, "Option": "Depth", "Depth": 5.0,
                     "Xi": 0.0, "Yi": 0.0, "R": 35.0}],
        "non_circ": [{"X": 0.0, "Y": 10.0, "Movement": "Free"},
                     {"X": 50.0, "Y": 2.0, "Movement": "Horiz"},
                     {"X": 100.0, "Y": 10.0, "Movement": "Fixed"}],
        "dloads": [[{"X": 0.0, "Y": 20.0, "Normal": 100.0},
                    {"X": 30.0, "Y": 20.0, "Normal": 100.0}]],
        "dloads2": [[{"X": 0.0, "Y": 20.0, "Normal": 50.0},
                     {"X": 30.0, "Y": 20.0, "Normal": 50.0}]],
        # One row per support type (blank/geosynthetic/nail/tieback/anchor) so the
        # editor exercises every Type value plus both Dir (tangent/axial) and Appl
        # (active/passive) enums and the tend1/tend2/spacing numerics.
        "reinforcement_lines": [
            {"x1": 0.0, "y1": 5.0, "x2": 40.0, "y2": 5.0, "t_max": 1000.0,
             "t_res": 800.0, "lp1": 2.0, "lp2": 3.0, "E": 2000.0, "area": 1.2,
             "type": "", "dir": "tangent", "appl": "active",
             "tend1": 0.0, "tend2": 0.0, "spacing": 1.0},
            {"x1": 0.0, "y1": 6.0, "x2": 38.0, "y2": 6.0, "t_max": 900.0,
             "t_res": 700.0, "lp1": 1.5, "lp2": 2.5, "E": 1800.0, "area": 1.1,
             "type": "geosynthetic", "dir": "tangent", "appl": "active",
             "tend1": 5.0, "tend2": 6.0, "spacing": 1.0},
            {"x1": 0.0, "y1": 7.0, "x2": 36.0, "y2": 7.0, "t_max": 800.0,
             "t_res": 600.0, "lp1": 1.0, "lp2": 2.0, "E": 1600.0, "area": 1.0,
             "type": "nail", "dir": "axial", "appl": "passive",
             "tend1": 10.0, "tend2": 12.0, "spacing": 1.5},
            {"x1": 0.0, "y1": 8.0, "x2": 34.0, "y2": 8.0, "t_max": 700.0,
             "t_res": 500.0, "lp1": 0.5, "lp2": 1.5, "E": 1400.0, "area": 0.9,
             "type": "tieback", "dir": "axial", "appl": "active",
             "tend1": 15.0, "tend2": 18.0, "spacing": 2.0},
            {"x1": 0.0, "y1": 9.0, "x2": 32.0, "y2": 9.0, "t_max": 600.0,
             "t_res": 400.0, "lp1": 0.5, "lp2": 1.0, "E": 1200.0, "area": 0.8,
             "type": "anchor", "dir": "axial", "appl": "passive",
             "tend1": 20.0, "tend2": 22.0, "spacing": 2.5},
        ],
        "pile_lines": [pile(20.0, 20.0, 20.0, 0.0, "passive"),
                       pile(35.0, 20.0, 35.0, 2.0, "active")],
        # ground_surface is None in this fixture, so LineLoadsEditor.apply performs
        # no snapping and the record must round-trip byte-for-byte.
        "line_loads": [{"x": 30.0, "y": 20.0, "P": 500.0, "angle": -90.0,
                        "label": "Facing plate"}],
        "seepage_bc": {
            "specified_heads": [
                {"head": 18.0, "coords": [(0.0, 0.0), (0.0, 18.0)]},
                {"head": 5.0, "coords": [(100.0, 0.0), (100.0, 5.0)]}],
            # Two flux BCs (v15) — the SeepBcEditor now edits fluxes in the master
            # list alongside heads/exit; multiple fluxes lock their list ordering.
            "specified_fluxes": [
                {"flux": 1.5, "coords": [(40.0, 20.0), (60.0, 20.0)]},
                {"flux": -0.75, "coords": [(70.0, 20.0), (85.0, 20.0)]}],
            "exit_face": [(60.0, 20.0), (100.0, 5.0)],
        },
        "seepage_bc2": {
            "specified_heads": [{"head": 10.0, "coords": [(0.0, 0.0), (0.0, 10.0)]}],
            "specified_fluxes": [{"flux": 2.25, "coords": [(30.0, 20.0), (45.0, 20.0)]}],
            "exit_face": [(60.0, 20.0), (100.0, 5.0)],
        },
        "has_seepage_bc2": True,
        "mesh": None,
    }


def _editor_norm(key, val):
    """Normalize a managed value for deep comparison: material-zone polygons hold
    shapely objects rebuilt fresh on apply, so compare their exterior coordinate
    rings + mat_id rather than object identity."""
    if key == "polygons":
        out = []
        for p in (val or []):
            poly = p.get("polygon")
            coords = list(poly.exterior.coords) if poly is not None else []
            out.append({"coords": coords, "mat_id": p.get("mat_id")})
        return out
    return val


def run_editor_roundtrip_test(test):
    """Open each studio CategoryEditor on a fully-populated record, apply with NO
    changes, and assert its managed keys are deep-equal in -> out.

    Returns (0.0, None) on success, or (None, message) naming the dropped/corrupted
    fields. Fast (one offscreen QApplication, no file or solver work).
    """
    import copy
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    from PySide6.QtWidgets import QApplication
    from studio.editors import CATEGORY_EDITORS

    QApplication.instance() or QApplication([])
    problems = []
    for key, editor in CATEGORY_EDITORS.items():
        keys = _EDITOR_MANAGED_KEYS.get(key, [])
        sd = _editor_fixture()
        before = {k: _editor_norm(k, copy.deepcopy(sd.get(k))) for k in keys}
        dlg = editor.build(sd, None)
        editor.apply(sd, dlg)
        for k in keys:
            problems += _roundtrip_diff(before[k], _editor_norm(k, sd.get(k)),
                                        f"{key}:{k}")
    if problems:
        return None, "editor round-trip dropped/corrupted data: " + "; ".join(problems[:6])
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


def run_deps_declared_test(test):
    """Every third-party module imported at module scope in xslope/ must be a declared
    runtime dependency in pyproject.toml.

    A module-scope import that isn't declared makes `pip install xslope` install a
    package that cannot be imported — while every test here still passes, because the
    suite runs against the source tree in a dev environment that happens to have the
    module. `twine check` doesn't catch it either; it only validates metadata. That is
    exactly how lxml (fileio, since 2026-06-08) and tabulate (solve, advanced) shipped
    undeclared through releases 0.1.50-0.1.54: six weeks of broken clean installs.

    Import-time is the criterion, not usage. A module imported lazily inside a function
    may legitimately be an optional extra (ezdxf in cad.py is, and __init__ imports only
    _version, so a plain `import xslope` never pulls it). One imported at module scope
    cannot be optional for any module a user imports.

    Returns (0.0, None) if clean, else (None, message).
    """
    import ast as _ast
    root = Path(__file__).parent
    pyproject = root / "pyproject.toml"
    if not pyproject.exists():
        return None, "pyproject.toml missing"

    # Declared runtime deps: the [project] dependencies array. Parsed by hand so the
    # check runs on 3.9/3.10 without tomllib.
    txt = pyproject.read_text()
    m = re.search(r"^dependencies\s*=\s*\[(.*?)\]", txt, re.S | re.M)
    if not m:
        return None, "could not find [project] dependencies in pyproject.toml"
    declared = {re.split(r"[<>=!\[; ]", s)[0].strip().lower()
                for s in re.findall(r'"([^"]+)"', m.group(1))}

    # Optional extras are allowed at module scope only in modules the package never
    # imports on its own; those are opt-in imports by the user.
    extras = set()
    for em in re.finditer(r'^\w[\w-]*\s*=\s*\[(.*?)\]', txt, re.S | re.M):
        for s in re.findall(r'"([^"]+)"', em.group(1)):
            extras.add(re.split(r"[<>=!\[; ]", s)[0].strip().lower())
    OPTIONAL_OK = {"ezdxf", "gmsh", "pyside6", "litellm", "keyring", "pyobjc-framework-cocoa"}

    stdlib = getattr(sys, "stdlib_module_names", frozenset())
    offenders = {}
    for py in sorted((root / "xslope").glob("*.py")):
        try:
            tree = _ast.parse(py.read_text())
        except SyntaxError as e:
            return None, f"could not parse {py.name}: {e}"
        for node in tree.body:                       # module scope only
            mods = []
            if isinstance(node, _ast.Import):
                mods = [a.name.split('.')[0] for a in node.names]
            elif isinstance(node, _ast.ImportFrom):
                if node.level:                        # relative -> in-package
                    continue
                if node.module:
                    mods = [node.module.split('.')[0]]
            for mod in mods:
                low = mod.lower()
                if mod in stdlib or mod == "xslope" or low in declared:
                    continue
                if low in OPTIONAL_OK:
                    continue
                offenders.setdefault(mod, []).append(py.name)

    if offenders:
        parts = [f"{mod} ({', '.join(sorted(set(files)))})"
                 for mod, files in sorted(offenders.items())]
        return None, ("imported at module scope but not a declared runtime dependency: "
                      + "; ".join(parts)
                      + " — add to [project] dependencies in pyproject.toml, or make the "
                        "import lazy if it is genuinely optional")
    return 0.0, None


def run_drawdown_tauff_test(test):
    """The Stage-2 undrained strength pipeline, checked against the worked example in
    Duncan, Wright & Brandon, *Soil Strength and Slope Stability*, 2nd ed., Table 9.2.

    An infinite slope, so the whole tau_ff chain -- eqs (4) K1, (5) interpolation,
    (6) Kf, (7)/(8) negative-stress fallback -- is exercised with hand-computable
    inputs and no slice machinery at all.

    gamma = 125 pcf, c' = 0, phi' = 40 deg, d = 2000 psf, psi = 20 deg, beta = 18.4 deg,
    fully submerged under 100 ft of water before a complete drawdown.

    The book carries K1 = 2.0 where the exact value is 2.0304, so its tau_ff (1585 psf
    at z = 5 ft) is a rounded figure. We assert against the exact recomputation and
    record the book's rounded values alongside; substituting K1 = 2.0 and Kf = 4.6
    reproduces 1585 to the digit.

    Also asserts the invariant that motivated the Stage-1 FS >= 1 guard: K1 rises
    monotonically as FS1 falls and equals Kf exactly at FS1 = 1, so K1 <= Kf there and
    eq (5) interpolates rather than extrapolating.

    Returns (0.0, None) on success, else (None, message) — a pass/fail test.
    """
    import numpy as np

    gamma, gamma_w = 125.0, 62.4
    c1, phi1, d_val, psi_val = 0.0, 40.0, 2000.0, 20.0
    # The 3H:1V slope is arctan(1/3) = 18.4349 deg, but the book carries 18.4 deg through
    # eqs (9.19)-(9.42). Use its rounded angle so the printed values are reproducible.
    beta = 18.4
    pr = np.radians(phi1)
    cos_phi, sin_phi = np.cos(pr), np.sin(pr)

    def tau_ff_of(sigma_fc, tau_fc):
        """Mirror of advanced.rapid_drawdown's Stage-2 strength block."""
        t_k1 = d_val + sigma_fc * np.tan(np.radians(psi_val))          # eq (9.26)
        t_kf = c1 + sigma_fc * np.tan(pr)                              # eq (9.28)
        kf_first = sigma_fc - c1 * cos_phi
        if abs(cos_phi) < 1e-12 or abs(kf_first) < 1e-12:
            return t_k1, t_kf, None, None, min(t_k1, t_kf)
        s3_k1 = sigma_fc + tau_fc * (sin_phi - 1) / cos_phi            # eq (9.14)
        s3_kf = kf_first * (1 - sin_phi) / (cos_phi ** 2)              # eq (9.15)
        if s3_k1 <= 0 or s3_kf <= 0:
            return t_k1, t_kf, None, None, min(t_k1, t_kf)
        K1 = (sigma_fc + tau_fc * (sin_phi + 1) / cos_phi) / s3_k1     # eq (9.10)
        Kf = ((sigma_fc + c1 * cos_phi) * (1 + sin_phi)) / (kf_first * (1 - sin_phi))  # eq (9.13)
        t_ff = ((Kf - K1) * t_k1 + (K1 - 1) * t_kf) / (Kf - 1)         # eq (9.11)
        return t_k1, t_kf, K1, Kf, t_ff

    problems = []
    br = np.radians(beta)

    # Book Table 9.2, per depth: sigma'_fc, tau_fc, tau_ff(Kc=1), tau_ff(Kc=Kf), drained
    book = {5: (282.0, 94.0, 2103.0, 237.0, 472.0),
            30: (1691.0, 562.0, 2615.0, 1419.0, 2833.0)}

    for z, (b_sfc, b_tfc, b_k1, b_kf, b_drained) in book.items():
        # Stage 1: submerged infinite slope. With no flow, sigma'_fc = gamma' z cos^2(beta)
        # (eq 9.23) and tau_fc = (gamma - gamma_w) z sin(beta) cos(beta) (eq 9.24).
        sigma_fc = (gamma - gamma_w) * z * np.cos(br) ** 2
        tau_fc = (gamma - gamma_w) * z * np.sin(br) * np.cos(br)
        if abs(sigma_fc - b_sfc) > 1.0:
            problems.append(f"z={z}: sigma'_fc {sigma_fc:.1f} vs book {b_sfc}")
        if abs(tau_fc - b_tfc) > 1.0:
            problems.append(f"z={z}: tau_fc {tau_fc:.1f} vs book {b_tfc}")

        t_k1, t_kf, K1, Kf, t_ff = tau_ff_of(sigma_fc, tau_fc)
        if abs(t_k1 - b_k1) > 1.0:
            problems.append(f"z={z}: tau_ff(Kc=1) {t_k1:.1f} vs book {b_k1}")
        if abs(t_kf - b_kf) > 1.0:
            problems.append(f"z={z}: tau_ff(Kc=Kf) {t_kf:.1f} vs book {b_kf}")
        if K1 is None or not (1.0 <= K1 <= Kf):
            problems.append(f"z={z}: K1={K1} outside [1, Kf={Kf}] — eq (5) would extrapolate")

        # Book's rounded K1=2.0, Kf=4.6 must reproduce its printed tau_ff exactly.
        rounded = ((4.6 - 2.0) * t_k1 + (2.0 - 1) * t_kf) / (4.6 - 1)
        book_tff = {5: 1585.0, 30: 2283.0}[z]
        if abs(rounded - book_tff) > 1.0:
            problems.append(f"z={z}: tau_ff at book's rounded K1/Kf {rounded:.1f} "
                            f"vs book {book_tff}")

        # Stage 3: total drawdown, zero pore pressure -> sigma' = gamma z cos^2(beta).
        drained = c1 + gamma * z * np.cos(br) ** 2 * np.tan(pr)        # eqs (9.37), (9.39)
        if abs(drained - b_drained) > 1.0:
            problems.append(f"z={z}: drained strength {drained:.1f} vs book {b_drained}")

    # K1 == Kf exactly at Stage-1 FS == 1, and K1 > Kf below it. This is why
    # advanced.rapid_drawdown refuses a Stage-1 FS < 1 rather than extrapolating.
    sigma_fc = 282.0
    tau_fail = c1 + sigma_fc * np.tan(pr)          # tau_fc when FS1 == 1
    prev_K1 = None
    for FS1 in (3.0, 2.0, 1.5, 1.0, 0.95, 0.9):
        _, _, K1, Kf, _ = tau_ff_of(sigma_fc, tau_fail / FS1)
        if K1 is None:
            continue
        if prev_K1 is not None and not K1 > prev_K1:
            problems.append(f"K1 not monotonically increasing as FS1 falls (FS1={FS1})")
        prev_K1 = K1
        if abs(FS1 - 1.0) < 1e-12 and abs(K1 - Kf) > 1e-6:
            problems.append(f"K1={K1:.6f} != Kf={Kf:.6f} at FS1 = 1")
        if FS1 < 1.0 and not K1 > Kf:
            problems.append(f"FS1={FS1} < 1 but K1={K1:.4f} <= Kf={Kf:.4f}")
        if FS1 > 1.0 and not K1 < Kf:
            problems.append(f"FS1={FS1} > 1 but K1={K1:.4f} >= Kf={Kf:.4f}")

    if problems:
        return None, "Duncan/Wright/Brandon Table 9.2: " + "; ".join(problems)
    return 0.0, None


def run_drawdown_guard_test(test):
    """`rapid_drawdown` must REFUSE a slope whose Stage-1 (full pool) FS < 1 rather
    than silently extrapolating eq (5) past the Kc=Kf envelope.

    Below FS1 = 1 the mobilized shear tau_fc exceeds the failure envelope, so K1 > Kf
    on every low-K slice and tau_ff falls below the physical floor -- eventually
    negative, where `max(0, tau_ff)` laundered it into a zero-strength slice. A search
    would then rank that trial surface as catastrophic on a fictitious FS ~ 0. The
    negative-stress fallback does not catch it: sigma'_3c stays positive.

    Returns (0.0, None) on success, else (None, message).
    """
    import copy
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope.advanced import rapid_drawdown

    path = test['file']
    if not Path(path).exists():
        return None, f"missing fixture {path}"
    base = load_slope_data(path)
    problems = []

    def run(scale):
        data = copy.deepcopy(base)
        for m in data['materials']:
            for key in ('c', 'phi', 'd', 'psi'):
                m[key] = m.get(key, 0) * scale
        ok, payload = generate_slices(data, circle=data['circles'][0],
                                      num_slices=40, debug=False)
        if not ok:
            return None
        return rapid_drawdown(payload[0], 'bishop', debug_level=0)

    # Unweakened: Stage 1 is stable, analysis proceeds.
    got = run(1.0)
    if got is None or not got[0]:
        problems.append(f"unweakened model should succeed, got {got}")
    elif got[1]['stage1_FS'] < 1.0:
        problems.append("fixture no longer has Stage-1 FS >= 1; pick another")

    # Weakened until Stage 1 fails: must be refused, not answered.
    got = run(0.5)
    if got is None:
        problems.append("weakened model failed to generate slices")
    elif got[0]:
        problems.append(f"weakened model (Stage-1 FS < 1) returned FS={got[1]['FS']:.4f} "
                        "instead of being refused")
    elif 'stable before drawdown' not in str(got[1]):
        problems.append(f"refused, but with an unexpected message: {got[1]}")

    if problems:
        return None, "; ".join(problems)
    return 0.0, None


def run_axial_mirror_test(test):
    """Axial (nail) reinforcement must be facing-invariant: a nailed wall and its
    exact left/right mirror must give the same FS in every method, on both a
    circular and a planar surface. Locks the right-facing axial sign conventions
    (vertical force component, facing detection against a vertical wall face,
    spencer/mprice full-vector negation) fixed while building Slide VP47/VP48.

    corps/lowe carry a small inherent directional asymmetry (their interslice
    inclinations are read from the geometry left-to-right), so they get a looser
    tolerance. Pass/fail test: returns (0.0, None) on success.
    """
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope import solve

    lp, rp = test['file'], test['file2']
    for fp in (lp, rp):
        if not os.path.exists(fp):
            return None, f"input missing: {fp}"
    left = load_slope_data(lp)
    right = load_slope_data(rp)

    problems = []
    for kind in ('circle', 'plane'):
        if kind == 'circle':
            ok_l, pl = generate_slices(left, circle=left['circles'][0], num_slices=50)
            ok_r, pr = generate_slices(right, circle=right['circles'][0], num_slices=50)
            methods = ('oms', 'bishop', 'janbu', 'spencer', 'mprice', 'corps', 'lowe')
        else:
            ok_l, pl = generate_slices(left, non_circ=left['non_circ'], num_slices=50)
            ok_r, pr = generate_slices(right, non_circ=right['non_circ'], num_slices=50)
            methods = ('janbu', 'spencer', 'mprice', 'corps', 'lowe')
        if not (ok_l and ok_r):
            problems.append(f"{kind}: slice generation failed")
            continue
        df_l, df_r = pl[0], pr[0]
        for m in methods:
            fn = getattr(solve, m)
            ok1, r1 = fn(df_l.copy())
            ok2, r2 = fn(df_r.copy())
            if ok1 != ok2:
                # both-fail is acceptable (method admissibility, mirror-symmetric);
                # one-sided failure is a facing bug
                problems.append(f"{kind}/{m}: one facing failed "
                                f"(L={'ok' if ok1 else 'fail'}, R={'ok' if ok2 else 'fail'})")
                continue
            if not ok1:
                continue
            tol = 0.01 if m in ('corps', 'lowe') else 1e-3
            rel = abs(r1['FS'] - r2['FS']) / max(abs(r1['FS']), 1e-9)
            if rel > tol:
                problems.append(f"{kind}/{m}: L={r1['FS']:.4f} R={r2['FS']:.4f} "
                                f"rel diff {rel:.1e} > {tol:.0e}")
    if problems:
        return None, "; ".join(problems)
    return 0.0, None


def run_mp_spencer_test(test):
    """Morgenstern-Price with f(x) == 1 must reproduce Spencer exactly, on BOTH
    slope facings. Spencer is the f == 1 special case of M-P, so any divergence
    means the M-P march (or its right-facing mirror in `_mp_extract`) is wrong.

    This is the S3b gate. It is asserted here because the right-facing mirror is
    otherwise only exercised implicitly, and a stale docstring once claimed the
    facing was unsupported when it had in fact been implemented.

    Returns (0.0, None) on success, else (None, message) — a pass/fail test, like
    the round-trip checks.
    """
    from xslope.fileio import load_slope_data
    from xslope.slice import generate_slices
    from xslope import solve

    path = test['file']
    if not os.path.exists(path):
        return None, f"input missing: {path}"

    data = load_slope_data(path)
    if not data.get('circles'):
        return None, f"{path} has no circle to analyze"

    ok, payload = generate_slices(data, circle=data['circles'][0],
                                  num_slices=test.get('num_slices', 40), debug=False)
    if not ok:
        return None, f"slice generation failed: {payload}"
    slice_df = payload[0]

    facing = 'right' if slice_df['y_ct'].values[0] > slice_df['y_ct'].values[-1] else 'left'
    expected = test.get('facing')
    if expected and expected != facing:
        return None, f"{path} is {facing}-facing, expected {expected}-facing"

    ok_s, res_s = solve.spencer(slice_df)
    if not ok_s:
        return None, f"spencer failed on {facing}-facing: {res_s}"
    ok_m, res_m = solve.mprice(slice_df, f_type='constant')
    if not ok_m:
        return None, f"mprice(f=constant) failed on {facing}-facing: {res_m}"

    diff = abs(res_s['FS'] - res_m['FS'])
    tol = test.get('tolerance', 1e-8)
    if diff > tol:
        return None, (f"{facing}-facing: mprice(f=1) FS={res_m['FS']:.8f} != "
                      f"spencer FS={res_s['FS']:.8f} (diff {diff:.2e} > {tol:.0e})")
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


# Every XML tag path (and attribute) that GeoStudio itself writes, among those
# export_gsz emits. Harvested from Seequent-authored .gsz files; the FILES are
# their copyrighted Materials and are not in this repo, but the tag vocabulary is
# a fact about the format, not a copy of their content.
#
# This guards the failure that a round-trip test CANNOT catch: a reader and writer
# that share the same wrong idea of the schema agree perfectly with each other and
# produce a file GeoStudio rejects. (That is exactly what happened — the first
# exporter put PiezometricSurfaces under <Geometry> instead of the StabilityItem,
# and omitted <Lines> entirely, and the round-trip test passed regardless.)
GSZ_SCHEMA_PATHS = {
    "/GSIData",
    "/GSIData/Analyses",
    "/GSIData/Analyses/Analysis",
    "/GSIData/Analyses/Analysis/GeometryId",
    "/GSIData/Analyses/Analysis/ID",
    "/GSIData/Analyses/Analysis/Kind",
    "/GSIData/Analyses/Analysis/Method",
    "/GSIData/Analyses/Analysis/Name",
    "/GSIData/Analyses@Len",
    "/GSIData/Contexts",
    "/GSIData/Contexts/Context",
    "/GSIData/Contexts/Context/AnalysisID",
    "/GSIData/Contexts/Context/GeometryUsesMaterials",
    "/GSIData/Contexts/Context/GeometryUsesMaterials/GeometryUsesMaterial",
    "/GSIData/Contexts/Context/GeometryUsesMaterials/GeometryUsesMaterial@Entry",
    "/GSIData/Contexts/Context/GeometryUsesMaterials/GeometryUsesMaterial@ID",
    "/GSIData/Contexts/Context/GeometryUsesMaterials@Len",
    "/GSIData/Contexts/Context/IsDefined",
    "/GSIData/Contexts@Len",
    "/GSIData/Coordinates",
    "/GSIData/Coordinates/EngCoords",
    "/GSIData/Coordinates/EngCoords@HorzScale",
    "/GSIData/Coordinates/EngCoords@LockScales",
    "/GSIData/Coordinates/EngCoords@MaxSnapDist",
    "/GSIData/Coordinates/EngCoords@UnitSystem",
    "/GSIData/Coordinates/EngCoords@VertScale",
    "/GSIData/Coordinates/EngCoords@XPageLeft",
    "/GSIData/Coordinates/EngCoords@XPageOrg",
    "/GSIData/Coordinates/EngCoords@XPageRight",
    "/GSIData/Coordinates/EngCoords@YPageBottom",
    "/GSIData/Coordinates/EngCoords@YPageOrg",
    "/GSIData/Coordinates/EngCoords@YPageTop",
    "/GSIData/Coordinates/PageCoords",
    "/GSIData/Coordinates/PageCoords@PageHeight",
    "/GSIData/Coordinates/PageCoords@PageWidth",
    "/GSIData/Coordinates/PageCoords@PageXOrg",
    "/GSIData/Coordinates/PageCoords@PageYOrg",
    "/GSIData/Coordinates/PageCoords@Units",
    "/GSIData/FileInfo",
    "/GSIData/FileInfo@Comments",
    "/GSIData/FileInfo@Title",
    "/GSIData/Geometries",
    "/GSIData/Geometries/Geometry",
    "/GSIData/Geometries/Geometry/Lines",
    "/GSIData/Geometries/Geometry/Lines/Line",
    "/GSIData/Geometries/Geometry/Lines/Line/ID",
    "/GSIData/Geometries/Geometry/Lines/Line/PointID1",
    "/GSIData/Geometries/Geometry/Lines/Line/PointID2",
    "/GSIData/Geometries/Geometry/Lines@Len",
    "/GSIData/Geometries/Geometry/Name",
    "/GSIData/Geometries/Geometry/Points",
    "/GSIData/Geometries/Geometry/Points/Point",
    "/GSIData/Geometries/Geometry/Points/Point@ID",
    "/GSIData/Geometries/Geometry/Points/Point@Pinned",
    "/GSIData/Geometries/Geometry/Points/Point@X",
    "/GSIData/Geometries/Geometry/Points/Point@Y",
    "/GSIData/Geometries/Geometry/Points@Len",
    "/GSIData/Geometries/Geometry/Regions",
    "/GSIData/Geometries/Geometry/Regions/Region",
    "/GSIData/Geometries/Geometry/Regions/Region/ID",
    "/GSIData/Geometries/Geometry/Regions/Region/Mesh",
    "/GSIData/Geometries/Geometry/Regions/Region/Mesh@Pattern",
    "/GSIData/Geometries/Geometry/Regions/Region/PointIDs",
    "/GSIData/Geometries/Geometry/Regions@Len",
    "/GSIData/Geometries/Geometry/Window",
    "/GSIData/Geometries/Geometry/Window/Base",
    "/GSIData/Geometries/Geometry/Window/Base@X",
    "/GSIData/Geometries/Geometry/Window/Base@Y",
    "/GSIData/Geometries/Geometry/Window/Zoom",
    "/GSIData/Geometries@Len",
    "/GSIData/Materials",
    "/GSIData/Materials/Material",
    "/GSIData/Materials/Material/Color",
    "/GSIData/Materials/Material/ID",
    "/GSIData/Materials/Material/Name",
    "/GSIData/Materials/Material/SlopeModel",
    "/GSIData/Materials/Material/StressStrain",
    "/GSIData/Materials/Material/StressStrain/CohesionPrime",
    "/GSIData/Materials/Material/StressStrain/PhiPrime",
    "/GSIData/Materials/Material/StressStrain/UnitWeight",
    "/GSIData/Materials/Material/StressModel",
    "/GSIData/Materials/Material/StressStrain/JointEffectiveCohesion",
    "/GSIData/Materials/Material/StressStrain/JointEffectiveCohesion@Missing",
    "/GSIData/Materials/Material/StressStrain/OCRatio",
    "/GSIData/Materials/Material/StressStrain/OCRatio@Missing",
    "/GSIData/Materials/Material/StressStrain/TensileStrength",
    "/GSIData/Materials/Material/StressStrain/TensileStrength@Missing",
    "/GSIData/Materials/Material/StressStrain/UseTensileStrength",
    "/GSIData/Materials@Len",
    "/GSIData/StabilityItems",
    "/GSIData/StabilityItems/StabilityItem",
    "/GSIData/StabilityItems/StabilityItem/AnalysisID",
    "/GSIData/StabilityItems/StabilityItem/Entry",
    # The analysis's LOCAL point list. The piezometric surface, tension crack, surcharges
    # and reinforcement all index INTO this, by Number — not into the geometry <Points>.
    "/GSIData/StabilityItems/StabilityItem/Entry/DataPoints",
    "/GSIData/StabilityItems/StabilityItem/Entry/DataPoints/DataPoint",
    "/GSIData/StabilityItems/StabilityItem/Entry/DataPoints/DataPoint@Number",
    "/GSIData/StabilityItems/StabilityItem/Entry/DataPoints/DataPoint@X",
    "/GSIData/StabilityItems/StabilityItem/Entry/DataPoints/DataPoint@Y",
    "/GSIData/StabilityItems/StabilityItem/Entry/DataPoints@Len",
    "/GSIData/StabilityItems/StabilityItem/Entry/MaterialUsesPiezs",
    "/GSIData/StabilityItems/StabilityItem/Entry/MaterialUsesPiezs/MaterialUsesPiez",
    "/GSIData/StabilityItems/StabilityItem/Entry/MaterialUsesPiezs/MaterialUsesPiez@ID",
    "/GSIData/StabilityItems/StabilityItem/Entry/MaterialUsesPiezs/MaterialUsesPiez@UsesID",
    "/GSIData/StabilityItems/StabilityItem/Entry/MaterialUsesPiezs@Len",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces/PiezometricSurface",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces/PiezometricSurface/CapSuction",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces/PiezometricSurface/DataPoints",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces/PiezometricSurface/DataPoints/DataPoint",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces/PiezometricSurface/DataPoints@Len",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces/PiezometricSurface/ID",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces/PiezometricSurface/MaxSuction",
    "/GSIData/StabilityItems/StabilityItem/Entry/PiezometricSurfaces@Len",
    "/GSIData/StabilityItems/StabilityItem/Entry/Seismic",
    "/GSIData/StabilityItems/StabilityItem/Entry/Seismic@Horizontal",
    "/GSIData/StabilityItems/StabilityItem/Entry/Seismic@Vertical",
    "/GSIData/StabilityItems@Len",
    "/GSIData/WaterItems",
    "/GSIData/WaterItems/WaterItem",
    "/GSIData/WaterItems/WaterItem/AnalysisID",
    "/GSIData/WaterItems/WaterItem/Entry",
    "/GSIData/WaterItems/WaterItem/Entry/ResultInputInfo",
    "/GSIData/WaterItems/WaterItem/Entry/ResultInputInfo/Option",
    "/GSIData/WaterItems/WaterItem/Entry/UnitWaterWeight",
    "/GSIData/WaterItems@Len",
    "/GSIData@AppVersion",
    "/GSIData@Version",
}


# The OTHER half of the conformance check, and the one that matters more.
#
# GSZ_SCHEMA_PATHS above says "every tag we write must be one GeoStudio writes". That
# can only catch a tag we invent. It CANNOT catch one we leave out -- and an omission is
# the more dangerous bug, because the file still opens and still draws. We shipped a
# .gsz with no <ComputedPhysics>, so GeoStudio did not know the analysis solved slope
# stability: the geometry rendered, the materials arrived named and correctly coloured,
# and their strengths were simply unreachable -- the Properties button greyed out, with
# nothing said. Only opening it in GeoStudio revealed that.
#
# These are the tag paths GeoStudio writes in EVERY file of a real corpus, minus the ones
# xslope consciously does not write (search definitions, solved results, the saved window
# scroll state, SEEP/W physics -- see tools/gsz_export_diff.py, which regenerates this).
# Like the vocabulary above, a list of tag names is a fact about the format and commits
# cleanly even though the vendor files themselves never can.
GSZ_REQUIRED_PATHS = {
    "/GSIData",
    "/GSIData/Analyses",
    "/GSIData/Analyses/Analysis",
    "/GSIData/Analyses/Analysis/ComputedPhysics",
    "/GSIData/Analyses/Analysis/ComputedPhysics@SlopeStability",
    "/GSIData/Analyses/Analysis/GeometryId",
    "/GSIData/Analyses/Analysis/ID",
    "/GSIData/Analyses/Analysis/IterationControls",
    "/GSIData/Analyses/Analysis/IterationControls/IterationControl",
    "/GSIData/Analyses/Analysis/IterationControls/IterationControl/Entry",
    "/GSIData/Analyses/Analysis/IterationControls/IterationControl/Entry@MaxIterations",
    "/GSIData/Analyses/Analysis/IterationControls/IterationControl/Entry@MaxReviewIterations",
    "/GSIData/Analyses/Analysis/IterationControls/IterationControl/Key",
    "/GSIData/Analyses/Analysis/IterationControls@Len",
    "/GSIData/Analyses/Analysis/Kind",
    "/GSIData/Analyses/Analysis/Method",
    "/GSIData/Analyses/Analysis/Name",
    "/GSIData/Analyses/Analysis/TimeIncrements",
    "/GSIData/Analyses/Analysis/TimeIncrements/Duration",
    "/GSIData/Analyses/Analysis/TimeIncrements/Duration@Missing",
    "/GSIData/Analyses/Analysis/TimeIncrements/TimeSteps",
    "/GSIData/Analyses/Analysis/TimeIncrements/TimeSteps/TimeStep",
    "/GSIData/Analyses/Analysis/TimeIncrements/TimeSteps/TimeStep@Save",
    "/GSIData/Analyses/Analysis/TimeIncrements/TimeSteps@Len",
    "/GSIData/Analyses@Len",
    "/GSIData/Contexts",
    "/GSIData/Contexts/Context",
    "/GSIData/Contexts/Context/AnalysisID",
    "/GSIData/Contexts/Context/GeometryUsesMaterials",
    "/GSIData/Contexts/Context/GeometryUsesMaterials/GeometryUsesMaterial",
    "/GSIData/Contexts/Context/GeometryUsesMaterials/GeometryUsesMaterial@Entry",
    "/GSIData/Contexts/Context/GeometryUsesMaterials/GeometryUsesMaterial@ID",
    "/GSIData/Contexts/Context/GeometryUsesMaterials@Len",
    "/GSIData/Contexts/Context/IsDefined",
    "/GSIData/Contexts@Len",
    "/GSIData/Coordinates",
    "/GSIData/Coordinates/EngCoords",
    "/GSIData/Coordinates/EngCoords@HorzScale",
    "/GSIData/Coordinates/EngCoords@LockScales",
    "/GSIData/Coordinates/EngCoords@MaxSnapDist",
    "/GSIData/Coordinates/EngCoords@UnitSystem",
    "/GSIData/Coordinates/EngCoords@VertScale",
    "/GSIData/Coordinates/EngCoords@XPageLeft",
    "/GSIData/Coordinates/EngCoords@XPageOrg",
    "/GSIData/Coordinates/EngCoords@XPageRight",
    "/GSIData/Coordinates/EngCoords@YPageBottom",
    "/GSIData/Coordinates/EngCoords@YPageOrg",
    "/GSIData/Coordinates/EngCoords@YPageTop",
    "/GSIData/Coordinates/PageCoords",
    "/GSIData/Coordinates/PageCoords@PageHeight",
    "/GSIData/Coordinates/PageCoords@PageWidth",
    "/GSIData/Coordinates/PageCoords@PageXOrg",
    "/GSIData/Coordinates/PageCoords@PageYOrg",
    "/GSIData/Coordinates/PageCoords@Units",
    "/GSIData/FileInfo",
    "/GSIData/FileInfo@AppVersion",
    "/GSIData/FileInfo@Author",
    "/GSIData/FileInfo@Comments",
    "/GSIData/FileInfo@Date",
    "/GSIData/FileInfo@FileVersion",
    "/GSIData/FileInfo@LastAuthor",
    "/GSIData/FileInfo@RevNumber",
    "/GSIData/FileInfo@Time",
    "/GSIData/Geometries",
    "/GSIData/Geometries/Geometry",
    "/GSIData/Geometries/Geometry/Lines",
    "/GSIData/Geometries/Geometry/Lines/Line",
    "/GSIData/Geometries/Geometry/Lines/Line/ID",
    "/GSIData/Geometries/Geometry/Lines/Line/PointID1",
    "/GSIData/Geometries/Geometry/Lines/Line/PointID2",
    "/GSIData/Geometries/Geometry/Lines@Len",
    "/GSIData/Geometries/Geometry/Name",
    "/GSIData/Geometries/Geometry/Points",
    "/GSIData/Geometries/Geometry/Points/Point",
    "/GSIData/Geometries/Geometry/Points/Point@ID",
    "/GSIData/Geometries/Geometry/Points/Point@X",
    "/GSIData/Geometries/Geometry/Points/Point@Y",
    "/GSIData/Geometries/Geometry/Points@Len",
    "/GSIData/Geometries/Geometry/Regions",
    "/GSIData/Geometries/Geometry/Regions/Region",
    "/GSIData/Geometries/Geometry/Regions/Region/ID",
    "/GSIData/Geometries/Geometry/Regions/Region/PointIDs",
    "/GSIData/Geometries/Geometry/Regions@Len",
    "/GSIData/Geometries@Len",
    "/GSIData/Materials",
    "/GSIData/Materials/Material",
    "/GSIData/Materials/Material/Color",
    "/GSIData/Materials/Material/ID",
    "/GSIData/Materials/Material/Name",
    "/GSIData/Materials/Material/SlopeModel",
    "/GSIData/Materials/Material/StressStrain",
    "/GSIData/Materials/Material/StressStrain/CohesionPrime",
    "/GSIData/Materials/Material/StressStrain/DisturbanceFactor",
    "/GSIData/Materials/Material/StressStrain/DisturbanceFactor@Missing",
    "/GSIData/Materials/Material/StressStrain/GeologicalStrengthIndex",
    "/GSIData/Materials/Material/StressStrain/GeologicalStrengthIndex@Missing",
    "/GSIData/Materials/Material/StressStrain/IntactRockParam",
    "/GSIData/Materials/Material/StressStrain/IntactRockParam@Missing",
    "/GSIData/Materials/Material/StressStrain/MaxConfiningStress",
    "/GSIData/Materials/Material/StressStrain/MaxConfiningStress@Missing",
    "/GSIData/Materials/Material/StressStrain/UnitWeight",
    "/GSIData/Materials@Len",
    "/GSIData/StabilityItems",
    "/GSIData/StabilityItems/StabilityItem",
    "/GSIData/StabilityItems/StabilityItem/AnalysisID",
    "/GSIData/StabilityItems/StabilityItem/Entry",
    "/GSIData/StabilityItems/StabilityItem/Entry/Seismic",
    "/GSIData/StabilityItems/StabilityItem/Entry/Seismic@Horizontal",
    "/GSIData/StabilityItems/StabilityItem/Entry/Seismic@Vertical",
    "/GSIData/StabilityItems@Len",
    "/GSIData/WaterItems",
    "/GSIData/WaterItems/WaterItem",
    "/GSIData/WaterItems/WaterItem/AnalysisID",
    "/GSIData/WaterItems/WaterItem/Entry",
    "/GSIData/WaterItems/WaterItem/Entry/UnitWaterWeight",
    "/GSIData/WaterItems/WaterItem/Entry/WaterBulkMod",
    "/GSIData/WaterItems@Len",
    "/GSIData@AppVersion",
    "/GSIData@Version",
}

# Every REQUIRED path is by construction a path GeoStudio writes, so the vocabulary
# must contain it. Keeping the two in sync here means a new required path can never
# be rejected by the allow-list it belongs to.
GSZ_SCHEMA_PATHS |= GSZ_REQUIRED_PATHS



def _write_synthetic_gsz(path):
    """Author a minimal GeoStudio .gsz — a two-analysis, two-material model over one
    region — to exercise the importer.

    The fixture is written here rather than shipped because GeoStudio's own sample and
    verification files are Seequent's copyrighted Materials and must not be committed to
    this repository. Authoring the XML also pins the exact schema the importer relies on
    -- and every fact below was learned by scoring the importer against SLOPE/W's own
    answers over a corpus of real files (tools/gsz_corpus.py). Each was silent when
    wrong, and each cost real accuracy:

      * material assignment lives in <Contexts>, PER ANALYSIS, not on the region;
      * an undrained material keeps its strength in <Cohesion>, not <CohesionPrime>;
      * <PiezometricSurfaces> hangs off StabilityItem/Entry, and indexes that analysis's
        LOCAL <DataPoints> -- not the shared geometry <Points> table (10% of FS);
      * a tension crack is switched off by DROPPING <TensionOption>; GeoStudio keeps the
        crack's geometry, so the points alone mean nothing (2% of FS);
      * a <Surcharge>'s <Pressure> is a UNIT WEIGHT, and the load is the weight of the
        fill between the drawn line and the ground -- not a uniform pressure (4% of FS);
      * reinforcement carries no direction: it is implied by <Type>;
      * a stability analysis whose water <Option> is "Parent" takes a whole pore-pressure
        FIELD from a SEEP/W run, on a binary PLY mesh -- and SLOPE/W raises the reservoir
        standing on the slope from that same field, storing no water surface at all
        (13% of FS, and most of it in the missing reservoir rather than the pressures).
    """
    import zipfile
    xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<GSIData Version="11.11" AppVersion="25.2.1.4">
  <FileInfo Title="synthetic" />
  <Coordinates><EngCoords UnitSystem="Metric" /></Coordinates>
  <Analyses Len="8">
    <Analysis><ID>1</ID><Name>dry</Name><Kind>SLOPE/W</Kind>
      <Method>Morgenstern-Price</Method><GeometryId>1</GeometryId></Analysis>
    <Analysis><ID>2</ID><Name>wet</Name><Kind>SLOPE/W</Kind>
      <Method>Spencer</Method><GeometryId>1</GeometryId></Analysis>
    <Analysis><ID>3</ID><Name>loaded</Name><Kind>SLOPE/W</Kind>
      <Method>Spencer</Method><GeometryId>1</GeometryId></Analysis>
    <Analysis><ID>4</ID><Name>crack off</Name><Kind>SLOPE/W</Kind>
      <Method>Spencer</Method><GeometryId>1</GeometryId></Analysis>
    <Analysis><ID>5</ID><Name>nailed</Name><Kind>SLOPE/W</Kind>
      <Method>Spencer</Method><GeometryId>1</GeometryId></Analysis>
    <Analysis><ID>6</ID><Name>seepage</Name><Kind>SEEP/W</Kind>
      <Method>SteadyState</Method><GeometryId>1</GeometryId></Analysis>
    <Analysis><ID>7</ID><Name>seep stability</Name><Kind>SLOPE/W</Kind>
      <ParentID>6</ParentID><Method>Spencer</Method><GeometryId>1</GeometryId></Analysis>
    <Analysis><ID>8</ID><Name>tall</Name><Kind>SLOPE/W</Kind>
      <Method>Morgenstern-Price</Method><GeometryId>2</GeometryId></Analysis>
  </Analyses>
  <Geometries Len="2">
    <Geometry><Name>2D</Name>
      <Points Len="6">
        <Point ID="1" X="0"  Y="0" /><Point ID="2" X="60" Y="0" />
        <Point ID="3" X="60" Y="30" /><Point ID="4" X="40" Y="30" />
        <Point ID="5" X="20" Y="10" /><Point ID="6" X="0"  Y="10" />
      </Points>
      <Regions Len="1">
        <Region><ID>1</ID><PointIDs>1,2,3,4,5,6</PointIDs></Region>
      </Regions>
    </Geometry>
    <!-- A SECOND geometry in the same file — the "same embankment at another height"
         shape (Borges & Cardoso Case 3). It REUSES point IDs 1-6 and region ID 1, but a
         taller slope: crest at y=45, not y=30. GeoStudio gives the geometry element no
         <ID>; the analysis names it by POSITION via <GeometryId>. Because the IDs collide,
         a single merged point table resolves every analysis to whichever geometry was
         parsed last — here geometry 2 — so analysis 1 would silently import THIS slope. -->
    <Geometry><Name>2D tall</Name>
      <Points Len="6">
        <Point ID="1" X="0"  Y="0" /><Point ID="2" X="60" Y="0" />
        <Point ID="3" X="60" Y="45" /><Point ID="4" X="40" Y="45" />
        <Point ID="5" X="20" Y="15" /><Point ID="6" X="0"  Y="15" />
      </Points>
      <Regions Len="1">
        <Region><ID>1</ID><PointIDs>1,2,3,4,5,6</PointIDs></Region>
      </Regions>
    </Geometry>
  </Geometries>
  <Materials Len="3">
    <Material><ID>1</ID><Name>strong</Name><Color>RGB=(211,201,137)</Color>
      <SlopeModel>MohrCoulomb</SlopeModel>
      <StressStrain><UnitWeight>20</UnitWeight><CohesionPrime>25</CohesionPrime>
        <PhiPrime>30</PhiPrime></StressStrain></Material>
    <Material><ID>2</ID><Name>weak</Name><Color>RGB=(0,128,255)</Color>
      <SlopeModel>MohrCoulomb</SlopeModel>
      <StressStrain><UnitWeight>18</UnitWeight><CohesionPrime>5</CohesionPrime>
        <PhiPrime>20</PhiPrime></StressStrain></Material>
    <Material><ID>3</ID><Name>undrained</Name><SlopeModel>UndrainedPhiZero</SlopeModel>
      <StressStrain><UnitWeight>19</UnitWeight><Cohesion>40</Cohesion></StressStrain></Material>
  </Materials>
  <Reinforcements Len="1">
    <Reinforcement><ID>1</ID><Name>Soil Nails</Name><Type>Nail</Type>
      <Spacing>2</Spacing><Tensile>100</Tensile><PlateCapacity>60</PlateCapacity>
      <PulloutResistance>10</PulloutResistance></Reinforcement>
  </Reinforcements>
  <Contexts Len="7">
    <Context><AnalysisID>1</AnalysisID><GeometryUsesMaterials Len="1">
        <GeometryUsesMaterial ID="Regions-1" Entry="1" /></GeometryUsesMaterials></Context>
    <Context><AnalysisID>2</AnalysisID><GeometryUsesMaterials Len="1">
        <GeometryUsesMaterial ID="Regions-1" Entry="2" /></GeometryUsesMaterials></Context>
    <Context><AnalysisID>3</AnalysisID><GeometryUsesMaterials Len="1">
        <GeometryUsesMaterial ID="Regions-1" Entry="3" /></GeometryUsesMaterials></Context>
    <Context><AnalysisID>4</AnalysisID><GeometryUsesMaterials Len="1">
        <GeometryUsesMaterial ID="Regions-1" Entry="1" /></GeometryUsesMaterials></Context>
    <Context><AnalysisID>5</AnalysisID><GeometryUsesMaterials Len="1">
        <GeometryUsesMaterial ID="Regions-1" Entry="1" /></GeometryUsesMaterials></Context>
    <Context><AnalysisID>7</AnalysisID><GeometryUsesMaterials Len="1">
        <GeometryUsesMaterial ID="Regions-1" Entry="1" /></GeometryUsesMaterials></Context>
    <Context><AnalysisID>8</AnalysisID><GeometryUsesMaterials Len="1">
        <GeometryUsesMaterial ID="Regions-1" Entry="1" /></GeometryUsesMaterials></Context>
  </Contexts>
  <WaterItems Len="7">
    <WaterItem><AnalysisID>1</AnalysisID>
      <Entry><UnitWaterWeight>9.807</UnitWaterWeight></Entry></WaterItem>
    <WaterItem><AnalysisID>2</AnalysisID>
      <Entry><ResultInputInfo><Option>PiezoSurface</Option></ResultInputInfo>
        <UnitWaterWeight>9.807</UnitWaterWeight></Entry></WaterItem>
    <WaterItem><AnalysisID>3</AnalysisID>
      <Entry><UnitWaterWeight>9.807</UnitWaterWeight></Entry></WaterItem>
    <WaterItem><AnalysisID>4</AnalysisID>
      <Entry><UnitWaterWeight>9.807</UnitWaterWeight></Entry></WaterItem>
    <WaterItem><AnalysisID>5</AnalysisID>
      <Entry><UnitWaterWeight>9.807</UnitWaterWeight></Entry></WaterItem>
    <WaterItem><AnalysisID>7</AnalysisID>
      <Entry><ResultInputInfo><Option>Parent</Option></ResultInputInfo>
        <UnitWaterWeight>9.807</UnitWaterWeight></Entry></WaterItem>
    <WaterItem><AnalysisID>8</AnalysisID>
      <Entry><UnitWaterWeight>9.807</UnitWaterWeight></Entry></WaterItem>
  </WaterItems>
  <StabilityItems Len="7">
    <StabilityItem><AnalysisID>1</AnalysisID>
      <Entry><Seismic Horizontal="" Vertical="" /></Entry></StabilityItem>

    <StabilityItem><AnalysisID>2</AnalysisID>
      <Entry>
        <DataPoints Len="2">
          <DataPoint Number="1" X="0"  Y="8" />
          <DataPoint Number="2" X="60" Y="24" />
        </DataPoints>
        <Seismic Horizontal="0.15" Vertical="" />
        <PiezometricSurfaces Len="1">
          <PiezometricSurface><ID>1</ID>
            <DataPoints Len="2"><DataPoint>1</DataPoint><DataPoint>2</DataPoint></DataPoints>
            <CapSuction>false</CapSuction><MaxSuction>0</MaxSuction>
          </PiezometricSurface>
        </PiezometricSurfaces>
        <MaterialUsesPiezs Len="1"><MaterialUsesPiez ID="2" UsesID="1" /></MaterialUsesPiezs>
      </Entry></StabilityItem>

    <StabilityItem><AnalysisID>3</AnalysisID>
      <Entry>
        <DataPoints Len="4">
          <DataPoint Number="1" X="40" Y="26" />
          <DataPoint Number="2" X="60" Y="26" />
          <DataPoint Number="3" X="42" Y="31" />
          <DataPoint Number="4" X="58" Y="33" />
        </DataPoints>
        <Seismic Horizontal="" Vertical="" />
        <TensionCrack>
          <TensionOption>Surface</TensionOption>
          <PctFilledWithWater>0.5</PctFilledWithWater>
          <DataPoints Len="2"><DataPoint>1</DataPoint><DataPoint>2</DataPoint></DataPoints>
        </TensionCrack>
        <Surcharges Len="1">
          <Surcharge><ID>1</ID>
            <DataPoints Len="2"><DataPoint>3</DataPoint><DataPoint>4</DataPoint></DataPoints>
            <Pressure>20</Pressure>
          </Surcharge>
        </Surcharges>
        <SomeFutureGeoStudioFeature>whatever</SomeFutureGeoStudioFeature>
      </Entry></StabilityItem>

    <StabilityItem><AnalysisID>4</AnalysisID>
      <Entry>
        <DataPoints Len="2">
          <DataPoint Number="1" X="40" Y="26" />
          <DataPoint Number="2" X="60" Y="26" />
        </DataPoints>
        <Seismic Horizontal="" Vertical="" />
        <TensionCrack>
          <PctFilledWithWater>1</PctFilledWithWater>
          <DataPoints Len="2"><DataPoint>1</DataPoint><DataPoint>2</DataPoint></DataPoints>
        </TensionCrack>
      </Entry></StabilityItem>

    <StabilityItem><AnalysisID>5</AnalysisID>
      <Entry>
        <DataPoints Len="2">
          <DataPoint Number="1" X="40" Y="30" />
          <DataPoint Number="2" X="50" Y="20" />
        </DataPoints>
        <Seismic Horizontal="" Vertical="" />
        <ReinforcementLines Len="1">
          <ReinforcementLine ID="1" Reinforcement="1" Point1Id="1" Point2Id="2" />
        </ReinforcementLines>
      </Entry></StabilityItem>

    <StabilityItem><AnalysisID>7</AnalysisID>
      <Entry><Seismic Horizontal="" Vertical="" /></Entry></StabilityItem>

    <StabilityItem><AnalysisID>8</AnalysisID>
      <Entry><Seismic Horizontal="" Vertical="" /></Entry></StabilityItem>
  </StabilityItems>
</GSIData>
"""
    nodes, elements, u = _synthetic_seep_field()
    with zipfile.ZipFile(path, 'w') as zf:
        zf.writestr("synthetic.xml", xml)
        # A SEEP/W-fed stability analysis keeps the TRANSFERRED field in its own result
        # folder, on the mesh SEEP/W solved on. Both are written the way GeoStudio writes
        # them -- a binary PLY, and a node.csv numbered from 1.
        zf.writestr("seep stability/Mesh.ply", _synthetic_ply(nodes, elements))
        zf.writestr("seep stability/000/node.csv",
                    "Node,PoreWaterPressure,XTotalStress\n" +
                    "".join(f"{i + 1},{v}\n" for i, v in enumerate(u)))


# The synthetic seepage model: a flat water table at elevation 15 over the fixture's
# slope, which rises from y=10 at x=20 to y=30 at x=40. So the reservoir stands 5 deep on
# the flat ground left of x=20 and tapers to nothing at x=25, where the face reaches
# elevation 15. Every number the importer must produce is hand-checkable from that.
_SEEP_WATER_LEVEL = 15.0
_SEEP_GAMMA_W = 9.807
_SEEP_COLUMNS = [(0.0, 10.0), (10.0, 10.0), (20.0, 10.0), (25.0, 15.0),
                 (30.0, 20.0), (40.0, 30.0), (50.0, 30.0), (60.0, 30.0)]
_SEEP_ROWS = 5


def _synthetic_seep_field():
    """Nodes, elements and nodal pore pressure for the fixture's SEEP/W analysis.

    The mesh CONFORMS to the slope: its top row of nodes is the ground surface itself,
    so a point sampled just below the ground always lands inside an element. (A mesh that
    merely covered the bounding box would leave the face outside it, and the reservoir
    would silently not be found.)
    """
    nodes, u = [], []
    for x, yg in _SEEP_COLUMNS:
        for r in range(_SEEP_ROWS):
            y = yg * r / (_SEEP_ROWS - 1)
            nodes.append((x, y))
            u.append(_SEEP_GAMMA_W * (_SEEP_WATER_LEVEL - y))    # hydrostatic; <0 above

    def nid(c, r):
        return c * _SEEP_ROWS + r + 1                            # PLY ids are 1-based

    elements = []                                                # (shape_kind, [ids])
    for c in range(len(_SEEP_COLUMNS) - 1):
        for r in range(_SEEP_ROWS - 1):
            quad = [nid(c, r), nid(c + 1, r), nid(c + 1, r + 1), nid(c, r + 1)]
            if c == 0 and r == 0:
                # One quad split into two triangles, so the fixture carries BOTH element
                # types and the importer's element_types has to come out right.
                elements.append((2, [quad[0], quad[1], quad[2]]))
                elements.append((2, [quad[0], quad[2], quad[3]]))
            else:
                elements.append((3, quad))

    # GeoStudio also writes the region's edges and corners into the same block. They are
    # not domain elements; importing them would leave zero-area cells in the mesh.
    elements.append((1, [nid(0, 0), nid(1, 0)]))                 # a line entity
    elements.append((0, [nid(0, 0)]))                            # a point entity
    return nodes, elements, u


def _synthetic_ply(nodes, elements):
    """Write nodes/elements as the binary PLY GeoStudio writes.

    The layout is copied from a real GeoStudio mesh, including the ``version`` block that
    precedes the nodes: a reader that skips it reads every coordinate 4 bytes out of step
    and produces plausible garbage. Writing it here is what pins that.
    """
    import struct
    header = (
        "ply\n"
        "format binary_little_endian 1.0\n"
        "element version 1\n"
        "property ushort major\n"
        "property ushort minor\n"
        f"element node {len(nodes)}\n"
        "property double x\n"
        "property double y\n"
        "property double z\n"
        f"element element {len(elements)}\n"
        "property uchar shape_kind\n"
        "property uchar integration_order\n"
        "property uchar category\n"
        "property uint owner\n"
        "property list uchar uint id\n"
        "element element_line_association 0\n"
        "property uint owner\n"
        "property uint id\n"
        "end_header\n"
    ).encode("ascii")

    body = bytearray(struct.pack("<HH", 1, 0))                   # version major/minor
    for x, y in nodes:
        body += struct.pack("<ddd", x, y, 0.0)
    for shape_kind, ids in elements:
        body += struct.pack("<BBBI", shape_kind, 2, 0, 1)
        body += struct.pack("<B", len(ids))
        body += struct.pack(f"<{len(ids)}I", *ids)
    return bytes(header + body)


def run_gsz_import_test(test):
    """Import a synthetic GeoStudio .gsz and check the model that comes out.

    Guards the three schema facts the importer depends on, each of which is easy to
    get wrong and silent when wrong:
      - material assignment is per ANALYSIS (<Contexts>), not per region, so the two
        analyses must yield different materials over the same geometry;
      - a piezometric surface indexes the shared <Points> table by ID;
      - the seismic coefficient rides on the analysis's <StabilityItem>.
    Also checks that an unsolved file reports the missing failure surface as a caveat
    rather than importing a wrong one. Returns (0.0, None) or (None, message).
    """
    import tempfile
    from xslope.geostudio import read_gsz, list_analyses, gsz_to_slope_data, gsz_style

    problems = []
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "synthetic.gsz")
        _write_synthetic_gsz(path)
        gsz = read_gsz(path)

        analyses = list_analyses(gsz)
        if len(analyses) != 8:
            return None, f"GeoStudio import: {len(analyses)} analyses, expected 8"

        sd1, cav1 = gsz_to_slope_data(gsz, 1)
        sd2, cav2 = gsz_to_slope_data(gsz, 2)
        sd3, cav3 = gsz_to_slope_data(gsz, 3)
        sd4, cav4 = gsz_to_slope_data(gsz, 4)
        sd5, cav5 = gsz_to_slope_data(gsz, 5)
        sd7, cav7 = gsz_to_slope_data(gsz, 7)
        sd8, cav8 = gsz_to_slope_data(gsz, 8)

        # --- a stability analysis fed by a parent SEEP/W run -------------------------
        # Its water <Option> is "Parent": pore pressure is a whole finite element FIELD,
        # which SLOPE/W stores on a binary PLY mesh in this analysis's own result folder.
        mesh = sd7.get('mesh')
        n_nodes = len(_SEEP_COLUMNS) * _SEEP_ROWS                    # 8 x 5 = 40
        n_elems = (len(_SEEP_COLUMNS) - 1) * (_SEEP_ROWS - 1) + 1    # 28 quads, one split
        if mesh is None:
            problems.append("the SEEP/W pore-pressure field did not import: no mesh")
        else:
            if len(mesh['nodes']) != n_nodes:
                problems.append(f"SEEP/W mesh has {len(mesh['nodes'])} nodes, expected "
                                f"{n_nodes} — the PLY was misread")
            # The line and point entities in the same PLY block are NOT domain elements.
            if len(mesh['elements']) != n_elems:
                problems.append(f"SEEP/W mesh has {len(mesh['elements'])} elements, "
                                f"expected {n_elems} — GeoStudio writes region edges and "
                                f"corners into the same block, and they are not elements")
            types = sorted(set(int(t) for t in mesh['element_types']))
            if types != [3, 4]:
                problems.append(f"SEEP/W element types {types}, expected [3, 4]")
            # The coordinates themselves: a reader that skips the PLY's <version> block
            # is 4 bytes out of step and returns plausible garbage.
            xs = [float(p[0]) for p in mesh['nodes']]
            ys = [float(p[1]) for p in mesh['nodes']]
            if (min(xs), max(xs), min(ys), max(ys)) != (0.0, 60.0, 0.0, 30.0):
                problems.append(f"SEEP/W mesh spans x {min(xs)}..{max(xs)}, y {min(ys)}.."
                                f"{max(ys)}, expected 0..60 and 0..30 — the binary PLY "
                                f"was decoded at the wrong offset")

        u = sd7.get('seep_u')
        if u is None or len(u) != n_nodes:
            problems.append(f"seep_u is {None if u is None else len(u)}, expected "
                            f"{n_nodes} nodal pore pressures")
        if [m['u'] for m in sd7['materials']] != ['seep']:
            problems.append(f"materials take u={[m['u'] for m in sd7['materials']]}, "
                            f"expected ['seep'] — a Parent water option is a seepage field")

        # THE RESERVOIR. GeoStudio stores no ponded-water object and, with a SEEP/W
        # field, no water surface either — yet SLOPE/W still puts the weight of the water
        # on the submerged face. It can only be deriving it from the head field, and it
        # is: water stands to y + u/gamma_w at the ground. Here that is a flat table at
        # elevation 15 over ground at elevation 10, so 5 deep (49.035) out to x=20, then
        # tapering to zero at x=25 where the slope face reaches 15. Missing this cost 13%
        # of the factor of safety on GeoStudio's own rapid-drawdown example — MORE than
        # the pore pressures did.
        dl7 = sd7['dloads']
        if len(dl7) != 1:
            problems.append(f"the reservoir implied by the SEEP/W field gave {len(dl7)} "
                            f"dload blocks, expected 1")
        else:
            deep = max(p['Normal'] for p in dl7[0])
            xmax = max(p['X'] for p in dl7[0])
            if abs(deep - _SEEP_GAMMA_W * 5.0) > 0.05:
                problems.append(f"ponded water from the SEEP/W field is {deep:.2f} at its "
                                f"deepest, expected {_SEEP_GAMMA_W * 5.0:.2f} (5 deep)")
            if abs(xmax - 25.0) > 0.3:
                problems.append(f"ponded water from the SEEP/W field reaches x={xmax:.2f}, "
                                f"expected 25 — where the face rises to the water level")
            if abs(min(p['Normal'] for p in dl7[0])) > 1e-6:
                problems.append("ponded water must taper to zero pressure at the "
                                "waterline, not stop at full depth")

        # An UNDRAINED material keeps its strength in <Cohesion>, not <CohesionPrime>.
        # Reading only the drained field gave c = 0, phi = 0 — a soil with no strength
        # and a meaningless factor of safety, with nothing said about it.
        m3 = sd3['materials'][0]
        if (m3['name'], m3['c'], m3['phi']) != ('undrained', 40.0, 0.0):
            problems.append(f"undrained material -> c={m3['c']}, phi={m3['phi']}, "
                            f"expected c=40 (from <Cohesion>), phi=0")

        # Tension crack: depth below ground, and water DEPTH in the crack (not a %).
        # The crest is at y=30 over x=40..60 and the crack line sits at y=26, so the
        # crack is 4 deep and — at PctFilledWithWater=0.5 — holds 2 of water.
        if abs(sd3['tcrack_depth'] - 4.0) > 1e-6:
            problems.append(f"tcrack_depth {sd3['tcrack_depth']}, expected 4")
        if abs(sd3['tcrack_water'] - 2.0) > 1e-6:
            problems.append(f"tcrack_water {sd3['tcrack_water']}, expected 2 "
                            f"(a DEPTH, = 50% of the crack depth, not a fraction)")

        # A crack with NO <TensionOption> is switched OFF — GeoStudio keeps its geometry
        # anyway. Analysis 4 has the very same crack line as analysis 3 and no option,
        # so it must import NO crack. Reading the points as a live crack put a 2 m
        # water-filled crack into models that had none, and cost 2% of FS in c'=0 soil.
        if sd4['tcrack_depth'] != 0.0 or sd4['tcrack_water'] != 0.0:
            problems.append(f"a tension crack with no <TensionOption> is switched OFF, "
                            f"but it imported at depth {sd4['tcrack_depth']}")

        # A <Surcharge> is a WEDGE OF FILL between the drawn line and the ground, and
        # <Pressure> is its unit weight — so the load varies with depth. Here a 20 kN/m3
        # fill runs from 1 deep at x=42 to 3 deep at x=58 over flat ground at y=30, so
        # the pressure must run 20 -> 60, NOT a constant 20.
        dl = sd3['dloads']
        if len(dl) != 1:
            problems.append(f"surcharge -> {len(dl)} dload blocks, expected 1")
        else:
            got = [(round(p['X'], 3), round(p['Normal'], 3)) for p in dl[0]]
            if got != [(42.0, 20.0), (58.0, 60.0)]:
                problems.append(
                    f"surcharge -> {got}, expected [(42.0, 20.0), (58.0, 60.0)]: "
                    f"<Pressure> is a UNIT WEIGHT and the load is the weight of the "
                    f"fill above the ground, not a uniform pressure")

        # The piezometric surface indexes the analysis's LOCAL <DataPoints>, not the
        # shared geometry <Points>. Read against the geometry it silently produced a
        # water table that doubled back on itself, and 10% of FS.
        if sd2['piezo_line'] != [(0.0, 8.0), (60.0, 24.0)]:
            problems.append(f"piezo line {sd2['piezo_line']}, expected "
                            f"[(0.0, 8.0), (60.0, 24.0)] — a piezometric surface indexes "
                            f"the analysis's LOCAL DataPoints, not the geometry points")

        # Reinforcement. GeoStudio quotes per-element values with an out-of-plane
        # spacing; xslope's lines are per unit width. Tensile 100 / spacing 2 = 50;
        # plate 60 / 2 = 30 at the end that is at the face; pullout 10 / 2 = 5 per unit
        # length, so the bond length that develops full capacity is 50 / 5 = 10.
        r5 = (sd5.get('reinforcement_lines') or [None])[0]
        if r5 is None:
            problems.append("reinforcement line was not imported")
        else:
            want = {'t_max': 50.0, 'lp1': 10.0, 'lp2': 10.0, 'tend1': 30.0, 'tend2': 0.0,
                    'spacing': 1.0, 'type': 'nail', 'dir': 'axial', 'appl': 'active'}
            bad = {k: r5.get(k) for k, v in want.items()
                   if not (abs(r5.get(k) - v) < 1e-9 if isinstance(v, float)
                           else r5.get(k) == v)}
            if bad:
                problems.append(f"reinforcement mapped to {bad}, expected "
                                f"{ {k: want[k] for k in bad} }")

        # Material colours come across, so an imported model still looks like the user's.
        if gsz_style(gsz, 1) != {'materials': {'0': {'color': '#d3c989'}}}:
            problems.append(f"material colour not imported: {gsz_style(gsz, 1)}")

        # FAIL LOUD: an element we do not recognise must be reported, never ignored.
        if not any('SomeFutureGeoStudioFeature' in c for c in cav3):
            problems.append("an unrecognised GeoStudio element was silently ignored")

        # Per-analysis material assignment: same geometry, different soil.
        m1 = sd1['materials'][0]
        m2 = sd2['materials'][0]
        if (m1['name'], m1['c'], m1['phi'], m1['gamma']) != ('strong', 25.0, 30.0, 20.0):
            problems.append(f"analysis 1 material {m1['name']}/{m1['c']}/{m1['phi']}")
        if (m2['name'], m2['c'], m2['phi'], m2['gamma']) != ('weak', 5.0, 20.0, 18.0):
            problems.append(f"analysis 2 material {m2['name']}/{m2['c']}/{m2['phi']}")

        # Geometry: the region's 6 points become one zone.
        if len(sd1['polygons']) != 1:
            problems.append(f"{len(sd1['polygons'])} polygons, expected 1")
        elif len(sd1['polygons'][0]['polygon'].exterior.coords) != 7:   # ring closes
            problems.append("zone ring did not round-trip")
        if not sd1['ground_surface'] or sd1['ground_surface'].is_empty:
            problems.append("no ground surface derived")

        # PER-ANALYSIS GEOMETRY. This file carries TWO geometries that reuse the same
        # point IDs (1-6) and region ID (1); each analysis names its own by <GeometryId>.
        # Analyses 1-7 bind geometry 1 (crest at y=30); analysis 8 binds geometry 2 (crest
        # at y=45). If GeometryId is ignored and the point tables are merged, both resolve
        # to whichever geometry was parsed LAST (geometry 2), so analysis 1 silently
        # imports the taller slope — distinct heights, distinct factors of safety, no word
        # said. This was the Borges & Cardoso Case 3 defect (7 m vs 8.75 m embankment).
        def _crest_y(sd):
            return max(xy[1] for pg in sd['polygons']
                       for xy in pg['polygon'].exterior.coords)
        if abs(_crest_y(sd1) - 30.0) > 1e-6:
            problems.append(f"analysis 1 resolved to crest y={_crest_y(sd1)}, expected 30 "
                            f"— its <GeometryId> was ignored, so it took the last geometry "
                            f"parsed instead of its own")
        if abs(_crest_y(sd8) - 45.0) > 1e-6:
            problems.append(f"analysis 8 (GeometryId 2) resolved to crest y={_crest_y(sd8)}, "
                            f"expected 45 — its own geometry, which reuses geometry 1's "
                            f"point IDs, was not selected")

        # Pore pressure: analysis 1 dry, analysis 2 on the piezo surface (points 7-8).
        if sd1['piezo_line'] or sd1['materials'][0]['u'] != 'none':
            problems.append("analysis 1 should be dry")
        if sd2['piezo_line'] != [(0.0, 8.0), (60.0, 24.0)]:
            problems.append(f"piezo line {sd2['piezo_line']}")
        if sd2['materials'][0]['u'] != 'piezo':
            problems.append(f"analysis 2 u={sd2['materials'][0]['u']}, expected piezo")

        # Seismic coefficient rides on the analysis.
        if sd1['k_seismic'] != 0.0 or sd2['k_seismic'] != 0.15:
            problems.append(f"k_seismic {sd1['k_seismic']}/{sd2['k_seismic']}")

        # Unsolved file: no surface invented, and the gap is reported.
        if sd1['circles'] or sd1['non_circ']:
            problems.append("invented a failure surface for an unsolved file")
        if not any('no failure surface' in c for c in cav1):
            problems.append("missing failure surface not reported as a caveat")

        # A SEEP/W field does not fit in a spreadsheet, so import_gsz writes it BESIDE
        # the .xlsx as _mesh.json + _seep.csv, and load_slope_data reads it back. If it
        # did not, the .xlsx would say every material takes its pore pressure from a
        # seepage solution and there would be none -- a dry model wearing a wet one's
        # clothes. Prove the whole product path: .gsz -> .xlsx (+ sidecars) -> reload.
        from xslope.geostudio import import_gsz
        from xslope.fileio import load_slope_data, default_template_path
        xlsx = os.path.join(td, "seep_roundtrip.xlsx")
        rcav = import_gsz(path, default_template_path(), xlsx, analysis_id=7)
        side = {os.path.basename(p) for p in os.listdir(td)}
        missing_side = {"seep_roundtrip_mesh.json", "seep_roundtrip_seep.csv"} - side
        if missing_side:
            problems.append(f"import_gsz did not write the SEEP/W sidecar(s) "
                            f"{sorted(missing_side)} — the .xlsx alone cannot carry a "
                            f"pore-pressure field, so the model would reload dry")
        try:
            reload = load_slope_data(xlsx)
        except Exception as e:
            # load_slope_data waives the failure-surface requirement only when a mesh is
            # present; if the sidecar was not written, it raises here instead. Report that
            # as the sidecar loss it is, rather than letting it read as an unrelated crash.
            reload = None
            problems.append(f"a SEEP/W field imported to .xlsx would not reload ({e}) — "
                            f"most likely its mesh/seep sidecar was not written beside it")
        if reload is not None:
            if reload.get("mesh") is None or reload.get("seep_u") is None:
                problems.append("a SEEP/W field imported to .xlsx did not reload — the "
                                "mesh or the pore pressure was lost between save and open")
            elif len(reload["seep_u"]) != n_nodes:
                problems.append(f"reloaded seep_u has {len(reload['seep_u'])} values, "
                                f"expected {n_nodes}")
            elif [m['u'] for m in reload['materials']] != ['seep']:
                problems.append(f"reloaded material u={[m['u'] for m in reload['materials']]}"
                                f", expected ['seep'] — the field reloaded but nothing reads it")
            else:
                # export_seep_u writes head+u only. import_seep_solution must load u as
                # real and the flow-net columns as NaN (not zeros), so a flow-net plot
                # fails visibly rather than drawing a field that was never solved.
                from xslope.seep import build_seep_data, import_seep_solution
                import numpy as _np
                sday = build_seep_data(reload["mesh"], reload, seep_bc=1)
                sol = import_seep_solution(sday, os.path.join(td, "seep_roundtrip_seep.csv"))
                if _np.any(_np.isnan(sol["u"])):
                    problems.append("imported seepage field lost its pore pressure to NaN")
                if not _np.all(_np.isnan(sol["phi"])):
                    problems.append("a field with no flow net came back with fabricated "
                                    "phi instead of NaN — a flow-net plot would draw a lie")

        # Export round-trip: slope_data -> .gsz -> slope_data must preserve the
        # geometry, the materials, and the piezo line.
        from xslope.geostudio import export_gsz
        out = os.path.join(td, "exported.gsz")
        export_gsz(sd2, out, analysis_name="rt")

        # SCHEMA CONFORMANCE — the check the round-trip cannot make. Every tag path we
        # write must be one GeoStudio itself writes; otherwise the file round-trips
        # through our own reader perfectly and GeoStudio still refuses it.
        import zipfile as _zip, xml.etree.ElementTree as _ET
        zf = _zip.ZipFile(out)
        root = _ET.fromstring(zf.read(zf.namelist()[0]))
        emitted = set()

        def _walk(e, p=""):
            q = f"{p}/{e.tag}"
            emitted.add(q)
            for a in e.attrib:
                emitted.add(f"{q}@{a}")
            for c in e:
                _walk(c, q)
        _walk(root)
        unknown = sorted(emitted - GSZ_SCHEMA_PATHS)
        if unknown:
            problems.append("export writes tag(s) GeoStudio does not use: "
                            + ", ".join(unknown[:3]))

        # THE OTHER DIRECTION, and the one that actually bites. The check above can only
        # catch a tag we write that GeoStudio does not; it is structurally incapable of
        # catching one we FAIL to write. Everything GeoStudio needs and we omitted got
        # through it — and the file then opens, draws, and looks fine while quietly
        # missing its physics. GeoStudio was the only thing that could tell us.
        #
        # So: every path GeoStudio writes in EVERY vendor file must be one we write too,
        # unless it is on the list of things we consciously do not write (a search
        # definition, solved results, the window scroll state, SEEP/W). Anything else
        # missing is a defect. This is what would have caught <ComputedPhysics> — without
        # it GeoStudio does not know the analysis solves slope stability, so the materials
        # arrive named and coloured with their strengths unreachable and nothing said.
        missing = sorted(GSZ_REQUIRED_PATHS - emitted)
        if missing:
            problems.append("export OMITS path(s) GeoStudio always writes: "
                            + ", ".join(missing[:3]))
        # And the structural facts a wrong-schema writer gets wrong:
        if not root.findall("./Geometries/Geometry/Lines/Line"):
            problems.append("export wrote no <Lines> — GeoStudio draws regions from them")
        if not root.findall("./StabilityItems/StabilityItem/Entry/"
                            "PiezometricSurfaces/PiezometricSurface"):
            problems.append("export put the piezometric surface outside the StabilityItem")
        if root.find("./Coordinates/EngCoords") is None:
            problems.append("export declared no unit system / view extent")
        # GeoStudio wants a piezometric surface's point IDs to ascend along the
        # polyline. Reusing a coincident region vertex hands it a low ID out of
        # sequence and it reports the surface as corrupt — so the surface gets its
        # own points, in order.
        pz = [int(d.text) for d in root.findall(
            "./StabilityItems/StabilityItem/Entry/PiezometricSurfaces/"
            "PiezometricSurface/DataPoints/DataPoint")]
        if pz != sorted(pz):
            problems.append(f"piezo surface point IDs not ascending: {pz}")
        # Materials must be visually distinguishable, not all on GeoStudio's default.
        cols = [m.findtext("Color") for m in root.findall("./Materials/Material")]
        if len(cols) > 1 and len(set(cols)) != len(cols):
            problems.append("materials exported without distinct colours")

        back, _ = gsz_to_slope_data(read_gsz(out), 1)
        b = back['materials'][0]
        if (b['name'], b['c'], b['phi'], b['gamma']) != ('weak', 5.0, 20.0, 18.0):
            problems.append(f"export round-trip material {b['name']}/{b['c']}/{b['phi']}")
        a0 = round(sd2['polygons'][0]['polygon'].area, 6)
        a1 = round(back['polygons'][0]['polygon'].area, 6)
        if a0 != a1:
            problems.append(f"export round-trip zone area {a1} vs {a0}")
        if back['piezo_line'] != sd2['piezo_line']:
            problems.append(f"export round-trip piezo {back['piezo_line']}")
        if back['k_seismic'] != sd2['k_seismic']:
            problems.append(f"export round-trip k_seismic {back['k_seismic']}")

    if problems:
        return None, "GeoStudio import: " + "; ".join(problems[:5])
    return 0.0, None


# A minimal Slide2 model in the sectioned ASCII grammar Slide2 writes. Authored here
# rather than shipped because Rocscience's own tutorial files are their copyrighted
# material and are not in this repository (the same reason the .gsz fixture above is
# synthesised). Two materials over a two-layer slope, a water table the upper layer
# draws from, and one specified circle — enough to exercise the cell-union geometry,
# the shared-strength/per-scenario-flag material merge, the 'hu' water mapping, and a
# round-trip through the .xlsx writer and load_slope_data.
_SYNTH_SLI = """
model description:
  version: 9.027
  methods: 6
  units: metric
  seismic: 0
  seismicv: 0
  water: hu
  direction: right to left
  gammaw: 9.81
  nummaterials: 2
  transient: no
  design_selection: 0

material types:
  soil1 = type: 0 water: 1 wtable: 1 c: 5 phi: 30 uw: 19 hutype: 0 withru: 0
  soil2 = type: 0 water: 1 wtable: 1 c: 25 phi: 32 uw: 20 hutype: 0 withru: 0

vertices:
  1 x: 0  y: 0
  2 x: 40  y: 0
  3 x: 40  y: 6
  4 x: 0  y: 6
  5 x: 40  y: 18
  6 x: 20  y: 18
  7 x: 0  y: 9
  8 x: 40  y: 9

cells:
  1  vertices: [1,2,3] material: soil2
  2  vertices: [1,3,4] material: soil2
  3  vertices: [4,3,5] material: soil1
  4  vertices: [4,5,6] material: soil1

water table:
  from_boreholes: 0  vertices: [7,8]

exterior:
  1  vertices: [1,2,5,6,4]

centers:
  1 x: 12 y: 30 r: 26 unique_id: {AAAAAAAA-0000-0000-0000-000000000001}

grids:

surfaces:

anchors:
"""


def _write_synthetic_slim(path):
    """Author a minimal Slide2 .slim — a ZIP of one .sli plus empty sidecars."""
    import zipfile
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("synthetic.sli", _SYNTH_SLI)
        zf.writestr("synthetic.slv", "")     # view options — ignored by the reader
        zf.writestr("synthetic.sltm", "")    # theme — ignored


def run_slide2_import_test(test):
    """Import a synthetic Slide2 .slim and check the model that comes out.

    Guards the schema facts the importer depends on, each silent when wrong:
      - triangular cells are unioned per material into one zone each, so two
        materials yield two polygons over the shared vertex table;
      - c/phi/uw come across, and only Mohr-Coulomb (type 0) maps cleanly;
      - a 'hu' water table becomes a piezo line, and only materials with wtable = 1
        draw pore pressure from it;
      - a specified circle is imported (a search would not be);
      - the whole model round-trips through the .xlsx writer and load_slope_data.
    """
    import tempfile
    from xslope.slide2 import (read_slim, read_slide2, list_scenarios,
                               slide2_to_slope_data, import_slmd)
    from xslope.fileio import (save_slope_data_to_xlsx, load_slope_data,
                               default_template_path)

    problems = []
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "synthetic.slim")
        _write_synthetic_slim(path)

        d = read_slim(path)
        scens = list_scenarios(d)
        if len(scens) != 1:
            return None, f"Slide2 import: {len(scens)} scenarios, expected 1"

        sd, caveats = slide2_to_slope_data(d)

        # Two materials, two zones, strengths intact.
        mats = sd["materials"]
        if len(mats) != 2:
            problems.append(f"{len(mats)} materials, expected 2")
        else:
            if (round(mats[0]["c"], 3), round(mats[0]["phi"], 3),
                    round(mats[0]["gamma"], 3)) != (5.0, 30.0, 19.0):
                problems.append(f"soil1 came across as c={mats[0]['c']} phi="
                                f"{mats[0]['phi']} gamma={mats[0]['gamma']}, "
                                f"expected 5/30/19")
            if (round(mats[1]["c"], 3), round(mats[1]["phi"], 3),
                    round(mats[1]["gamma"], 3)) != (25.0, 32.0, 20.0):
                problems.append(f"soil2 came across as c={mats[1]['c']} phi="
                                f"{mats[1]['phi']} gamma={mats[1]['gamma']}, "
                                f"expected 25/32/20")
        if len(sd["polygons"]) != 2:
            problems.append(f"{len(sd['polygons'])} polygons, expected 2 (the "
                            f"triangular cells did not union per material)")

        # 'hu' water table -> piezo line; both materials have wtable = 1.
        if len(sd["piezo_line"]) != 2:
            problems.append(f"piezo_line has {len(sd['piezo_line'])} points, "
                            f"expected 2 (the water table did not import)")
        if [m["u"] for m in mats] != ["piezo", "piezo"]:
            problems.append(f"materials take u={[m['u'] for m in mats]}, expected "
                            f"['piezo', 'piezo'] — a 'hu' water table is a piezo line")

        # A specified circle is imported; a search would not be.
        if len(sd["circles"]) != 1:
            problems.append(f"{len(sd['circles'])} circles, expected 1")
        elif round(sd["circles"][0]["R"], 3) != 26.0:
            problems.append(f"circle R={sd['circles'][0]['R']}, expected 26")

        # A search-only model must report its missing surface, never invent one.
        d2 = read_slim(path)
        d2["scenarios"][0]["sections"]["centers"] = []
        sd_ns, cav_ns = slide2_to_slope_data(d2)
        if sd_ns["circles"] or sd_ns["non_circ"]:
            problems.append("a search-only scenario imported a surface it should not")
        if not any("no failure surface" in c for c in cav_ns):
            problems.append("a search-only scenario did not report its missing surface")

        # read_slide2 dispatches on extension and agrees with read_slim.
        if len(list_scenarios(read_slide2(path))) != 1:
            problems.append("read_slide2 did not dispatch a .slim correctly")

        # Round-trip: the model must survive the .xlsx writer AND reload — solving the
        # raw dict is not enough (load_slope_data is the real consumer).
        xlsx = os.path.join(td, "synthetic.xlsx")
        rcav = import_slmd(path, default_template_path(), xlsx)
        if not os.path.exists(xlsx):
            problems.append("import_slmd did not write the .xlsx")
        else:
            reloaded = load_slope_data(xlsx)
            if len(reloaded["materials"]) != 2:
                problems.append(f"reloaded model has {len(reloaded['materials'])} "
                                f"materials, expected 2")
            if not reloaded.get("circular") or len(reloaded["circles"]) != 1:
                problems.append("the circle did not survive the .xlsx round-trip")
            if len(reloaded.get("piezo_line") or []) != 2:
                problems.append("the piezo line did not survive the .xlsx round-trip")
        if not isinstance(rcav, list):
            problems.append("import_slmd did not return a caveat list")

    if problems:
        return None, "Slide2 import: " + "; ".join(problems[:5])
    return 0.0, None


# --------------------------------------------------------------------------------------
# RS2 (.fez) import test
# --------------------------------------------------------------------------------------

# A hand-authored RS2 model, synthesised so the test carries NO Rocscience file. The
# real .fea (megabytes of finite-element mesh) is reduced to just the sections the
# importer reads: the model description with its SSR settings and the version-on-the-
# next-line dialect, two Mohr-Coulomb materials, a two-zone geometry (an external
# boundary cut by one material boundary), the coarse material-mesh seed triangles that
# label the zones, and a piezometric line the two materials draw from. Enough to
# exercise the boundary polygonisation, the seed-triangle material assignment, the
# next-line value dialect, the water mapping, and a round-trip through the .xlsx writer
# and load_slope_data.
_SYNTH_FEA = """model description:
version:
11.028
title:
analysis:  solid
strength_reduction_analysis: ON
auto_SRF: ON
change_in_SRF: 0.2
initial_SRF: 1
final_SRF: 2
maxiter_SRF: 500
tolerance_SRF: 0.001
delta_FS: 0.01
tensilestrength_SRF: 1
RFCunits: Metric kPa
stages: 1

material types:
material 1: soilA
 solid properties:
  rhoS: 2 rhoF: 1 porosity: 0.5
 Elastic Properties: LinearElastic
  nu: 0.3 E: 50000
 Plasticity Specifications: MohrCoulomb
  C: 5 phi: 30 dil: 0 T: 5 Cr: 5 phir: 30 Tr: 5 Apply_SSR: 1
material 2: soilB
 solid properties:
  rhoS: 2 rhoF: 1 porosity: 0.5
 Elastic Properties: LinearElastic
  nu: 0.3 E: 50000
 Plasticity Specifications: MohrCoulomb
  C: 25 phi: 32 dil: 0 T: 25 Cr: 25 phir: 32 Tr: 25 Apply_SSR: 1

material properties:
soilA
1 19 0 50000 8333 0.3 0 20000 20000 20000 0.2 0.2 0.2 1 0 5 0 30 30 5 5 100000
soilB
1 20 0 50000 8333 0.3 0 20000 20000 20000 0.2 0.2 0.2 1 0 25 0 32 32 25 25 100000

material piezos:
 1
 1

groundwater setup:
  gw_type: "Static Analysis"

piezos:
1
1
2
0 6
40 6

materials mesh start:
  num elements: 2
  num valid stages: 1
  element 0 start:
    P1: 20, 10
    P2: 21, 10
    P3: 20, 11
    materials start:
      stage 1: 1
    materials end:
  element 0 end:
  element 1 start:
    P1: 20, 2
    P2: 21, 2
    P3: 20, 3
    materials start:
      stage 1: 2
    materials end:
  element 1 end:
materials mesh end:

v6 geometry start:
  num boundaries: 2
  boundary 1 start:
    type: "external"
    guid: "{00000000-0000-0000-0000-000000000001}"
    vertices start:
      dpoint array start:
        num points: 5
        0: 0, 0
        1: 40, 0
        2: 40, 18
        3: 20, 18
        4: 0, 9
      dpoint array end:
    vertices end:
  boundary 1 end:
  boundary 2 start:
    type: "material"
    guid: "{00000000-0000-0000-0000-000000000002}"
    vertices start:
      dpoint array start:
        num points: 2
        0: 0, 4.5
        1: 40, 4.5
      dpoint array end:
    vertices end:
  boundary 2 end:
v6 geometry end:
"""


def _write_synthetic_fez(path):
    """Author a minimal RS2 .fez — a ZIP of one .fea plus empty sidecars."""
    import zipfile
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("synthetic.fea", _SYNTH_FEA)
        zf.writestr("synthetic.p2m", "")       # binary mesh sidecar — ignored
        zf.writestr("synthetic.config", "")    # view config — ignored


def run_rs2_import_test(test):
    """Import a synthetic RS2 .fez and check the model that comes out.

    Guards the schema facts the importer depends on, each silent when wrong:
      - the 'version:' then '11.028' next-line dialect parses;
      - two Mohr-Coulomb materials come across with c/phi (from 'material types:')
        and unit weight (from the 'material properties:' numeric row);
      - an external boundary cut by one material boundary polygonises into two zones,
        each labelled by its material-mesh seed triangle;
      - a piezometric line becomes the piezo line and the materials draw from it;
      - the SSR settings are surfaced as metadata, never turned into an LEM search;
      - NO failure surface is imported (an SSR analysis has none), and that is said;
      - the model round-trips through the .xlsx writer and load_slope_data.
    """
    import tempfile
    from xslope.rs2 import read_fez, fez_to_slope_data, import_fez
    from xslope.fileio import (save_slope_data_to_xlsx, load_slope_data,
                               default_template_path)

    problems = []
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "synthetic.fez")
        _write_synthetic_fez(path)

        d = read_fez(path)
        if d["version"] != "11.028":
            problems.append(f"version parsed as {d['version']!r}, expected '11.028' "
                            f"(the next-line value dialect did not parse)")
        if d["srf"].get("strength_reduction_analysis") != "ON":
            problems.append("SSR settings did not parse from the model description")

        sd, caveats = fez_to_slope_data(d)

        # Two materials, two zones, strengths intact.
        mats = sd["materials"]
        if len(mats) != 2:
            problems.append(f"{len(mats)} materials, expected 2")
        else:
            if (round(mats[0]["c"], 3), round(mats[0]["phi"], 3),
                    round(mats[0]["gamma"], 3)) != (5.0, 30.0, 19.0):
                problems.append(f"soilA came across as c={mats[0]['c']} phi="
                                f"{mats[0]['phi']} gamma={mats[0]['gamma']}, "
                                f"expected 5/30/19")
            if (round(mats[1]["c"], 3), round(mats[1]["phi"], 3),
                    round(mats[1]["gamma"], 3)) != (25.0, 32.0, 20.0):
                problems.append(f"soilB came across as c={mats[1]['c']} phi="
                                f"{mats[1]['phi']} gamma={mats[1]['gamma']}, "
                                f"expected 25/32/20")
        if len(sd["polygons"]) != 2:
            problems.append(f"{len(sd['polygons'])} polygons, expected 2 (the external "
                            f"boundary did not polygonise into two zones)")

        # A piezometric line -> piezo line; both materials draw from it.
        if len(sd["piezo_line"]) != 2:
            problems.append(f"piezo_line has {len(sd['piezo_line'])} points, "
                            f"expected 2 (the water table did not import)")
        if [m["u"] for m in mats] != ["piezo", "piezo"]:
            problems.append(f"materials take u={[m['u'] for m in mats]}, expected "
                            f"['piezo', 'piezo']")

        # An RS2 model imports NO failure surface, and must say so.
        if sd["circles"] or sd["non_circ"]:
            problems.append("an RS2 model imported a failure surface it should not have")
        if not any("no failure surface" in c for c in caveats):
            problems.append("the missing failure surface was not reported")
        if not any("Shear-Strength-Reduction" in c for c in caveats):
            problems.append("the SSR analysis metadata was not reported as a caveat")

        # Round-trip: the geometry must survive the .xlsx writer AND reload. An RS2
        # import carries no surface, so give it one (as a user must) before saving —
        # load_slope_data rejects a surface-less file, which is the real consumer.
        sd["circular"] = True
        sd["circles"] = [{"Xo": 20.0, "Yo": 30.0, "R": 26.0, "Depth": 4.0}]
        xlsx = os.path.join(td, "synthetic.xlsx")
        save_slope_data_to_xlsx(sd, xlsx, template=default_template_path())
        reloaded = load_slope_data(xlsx)
        if len(reloaded["materials"]) != 2:
            problems.append(f"reloaded model has {len(reloaded['materials'])} "
                            f"materials, expected 2")
        if len(reloaded["polygons"]) != 2:
            problems.append("the two zones did not survive the .xlsx round-trip")
        if len(reloaded.get("piezo_line") or []) != 2:
            problems.append("the piezo line did not survive the .xlsx round-trip")

        # import_fez writes an .xlsx and returns the caveat list (the surface-less file
        # it writes is intentionally incomplete — an SSR model has no LEM surface).
        out2 = os.path.join(td, "direct.xlsx")
        rcav = import_fez(path, default_template_path(), out2)
        if not os.path.exists(out2):
            problems.append("import_fez did not write the .xlsx")
        if not isinstance(rcav, list):
            problems.append("import_fez did not return a caveat list")

    if problems:
        return None, "RS2 import: " + "; ".join(problems[:5])
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
    if test_type == 'editor_roundtrip':
        return run_editor_roundtrip_test(test)
    if test_type == 'dxf':
        return run_dxf_roundtrip_test(test)
    if test_type == 'gsz':
        return run_gsz_import_test(test)
    if test_type == 'slide2':
        return run_slide2_import_test(test)
    if test_type == 'rs2':
        return run_rs2_import_test(test)
    if test_type == 'template_sync':
        return run_template_sync_test(test)
    if test_type == 'deps_declared':
        return run_deps_declared_test(test)
    if test_type == 'gsat_pair':
        return run_gsat_pair_test(test)
    if test_type == 'axial_mirror':
        return run_axial_mirror_test(test)
    if test_type == 'mp_spencer':
        return run_mp_spencer_test(test)
    if test_type == 'drawdown_tauff':
        return run_drawdown_tauff_test(test)
    if test_type == 'drawdown_guard':
        return run_drawdown_guard_test(test)
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
    elif test_type == 'seep_head':
        return run_seep_head_test(test)
    elif test_type == 'reliability':
        return run_reliability_test(test)
    elif test_type == 'design_search':
        return run_design_test(test)
    elif test_type == 'sensitivity':
        return run_sensitivity_test(test)
    else:
        return run_lem_test(test)


def _expected_and_tol(test, default_tolerance):
    """Expected value and tolerance for a test dict (shared by both runners)."""
    test_type = test.get('type', '?')
    if test_type == 'seep':
        expected = test.get('expected_flowrate')
        tol = test.get('tolerance', 0.05) * abs(expected) if expected else 0
    elif test_type in ('reliability', 'fem_reliability'):
        expected = test.get('expected_beta')
        tol = test.get('tolerance', 0.02)
    elif test_type == 'sensitivity':
        # the runner checks base/low/high internally; the framework-level
        # comparison re-checks the base row
        expected = float(test['expected_base']) if 'expected_base' in test else None
        tol = float(test.get('tolerance', 0.01))
    elif test_type in ('roundtrip', 'editor_roundtrip', 'template_sync', 'deps_declared', 'dxf', 'gsz', 'slide2', 'rs2', 'vg_kr',
                       'mesh_conform', 'seep_elements', 'fem_elements',
                       'mp_spencer', 'axial_mirror', 'drawdown_tauff', 'drawdown_guard',
                       'gsat_pair', 'seep_head'):
        expected = 0.0          # these return 0.0 on success (pass/fail tests)
        tol = 0.0
    else:
        expected = test.get('expected_fs')
        tol = test.get('tolerance', default_tolerance)
    return expected, tol


# Rough per-type cost ranks so the parallel scheduler starts the slow tests
# first (wall time is otherwise dominated by an FEM case landing last).
_COST_RANK = {'fem_reliability': 6, 'fem_ssrm': 5, 'fem_elements': 5,
              'reliability': 4, 'seep_elements': 3, 'seep': 3,
              'noncircular_search': 2, 'circular_search': 2}


def _parallel_worker(item):
    """Run one test in a worker process, stdout/stderr suppressed."""
    i, test = item
    import io as _io
    import contextlib as _ctx
    import warnings as _w
    _w.filterwarnings('ignore')
    t0 = time.time()
    buf = _io.StringIO()
    try:
        with _ctx.redirect_stdout(buf), _ctx.redirect_stderr(buf):
            computed, error_msg = run_test(test)
    except Exception as e:
        computed, error_msg = None, str(e)
    return i, computed, error_msg, time.time() - t0


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
    parser.add_argument('--gsz', action='store_true',
                        help='Run only the GeoStudio (.gsz) import test')
    parser.add_argument('--slide2', action='store_true',
                        help='Run only the Slide2 (.slim/.slmd) import test')
    parser.add_argument('--rs2', action='store_true',
                        help='Run only the RS2 (.fez) import test')
    parser.add_argument('--tolerance', type=float, default=0.01,
                        help='Default tolerance for FS comparison (default: 0.01)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print details for passing tests')
    parser.add_argument('--benchmark', default=None,
                        help='Run only tests whose benchmark id matches the '
                             'given comma-separated list; each entry matches '
                             'exactly or as a family prefix (VP28 runs VP28a, '
                             'VP28a-beta, ... but not VP280). For tag-only or '
                             'builder-only changes, run just the new tags '
                             'instead of the whole suite.')
    parser.add_argument('--skip-benchmarks', action='store_true',
                        help='Skip verification benchmark tests (annotations '
                             'carrying a benchmark=<ID> tag), e.g. the slow '
                             'SSRM dam cases')
    parser.add_argument('--jobs', default='auto',
                        help="Parallel worker processes: an integer or 'auto' "
                             "(cpu_count - 2, the default). Use --jobs 1 for "
                             "sequential runs (in-order rows, full solver "
                             "output — the debugging mode).")
    parser.add_argument('--quick', action='store_true',
                        help='For LEM problems that list several methods, check '
                             'only one method per problem (prefers Spencer) so '
                             'routine runs stay fast')
    args = parser.parse_args()

    # If no specific flags, run all
    run_all = not (args.lem or args.fem or args.seep or args.roundtrip or args.dxf
                   or args.gsz or args.slide2 or args.rs2)
    run_lem = args.lem or run_all
    run_fem = args.fem or run_all
    run_seep = args.seep or run_all
    run_roundtrip = args.roundtrip or run_all
    run_dxf = args.dxf or run_all
    run_gsz = args.gsz or run_all
    run_slide2 = args.slide2 or run_all
    run_rs2 = args.rs2 or run_all

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
        lem_sens = Path('docs/lem/sensitivity.md')
        if lem_sens.exists():
            tests.extend(parse_test_tags(lem_sens))

    if run_fem:
        fem_samples = Path('docs/fem/samples.md')
        if fem_samples.exists():
            tests.extend(parse_test_tags(fem_samples))
    # Verification corpus pages (docs/verification/*.md): tags are routed by
    # type so LEM, SSRM, and seepage assertions can live on the same page.
    for verification_md in sorted(Path('docs/verification').glob('*.md')):
        for t in parse_test_tags(verification_md):
            ttype = t.get('type', '')
            if ttype in ('fem_ssrm', 'fem_elements', 'fem_reliability'):
                if run_fem:
                    tests.append(t)
            elif ttype in ('seep', 'seep_elements', 'seep_head'):
                if run_seep:
                    tests.append(t)
            elif run_lem:
                tests.append(t)

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

    # Private test problems (CE 544 homework/exam keys, plus synthetic variants)
    # live in a separate repo kept out of the public xslope tree. Look for it as a
    # sibling directory or at $XSLOPE_PRIVATE_TESTS; scan its markdown files for
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
        # Fixtures live under tests/. Scope the walk there rather than over the whole
        # repo: sibling trees (plans/, ref_docs/) hold prose that quotes the
        # `<!-- test: ... -->` syntax as documentation, and parsing those yields
        # garbage entries with no 'file' key.
        scan_root = private_path / 'tests'
        if not scan_root.is_dir():
            scan_root = private_path
        for md in sorted(scan_root.rglob('*.md')):
            if md.name.lower() == 'readme.md':
                continue
            if 'ref_docs' in md.relative_to(private_path).parts:
                continue
            for t in parse_test_tags(md):
                if not t.get('file'):
                    continue    # a prose example of the tag syntax, not a fixture
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

    # Spencer is the f(x) == 1 special case of Morgenstern-Price, so the two must
    # agree to machine precision. Asserted on BOTH facings: the right-facing mirror
    # in _mp_extract is otherwise only exercised implicitly.
    if run_lem:
        for _fp, _facing in ((MP_SPENCER_LEFT, 'left'), (MP_SPENCER_RIGHT, 'right')):
            if Path(_fp).exists():
                tests.append({'type': 'mp_spencer', 'file': _fp, 'facing': _facing,
                              'method': 'mprice=spencer', 'source': f'{_facing}-facing'})

        # Axial reinforcement facing-invariance: nailed wall vs its exact mirror
        if Path(AXIAL_MIRROR_LEFT).exists() and Path(AXIAL_MIRROR_RIGHT).exists():
            tests.append({'type': 'axial_mirror', 'file': AXIAL_MIRROR_LEFT,
                          'file2': AXIAL_MIRROR_RIGHT, 'method': 'all-methods L==R',
                          'source': 'axial_mirror'})

        # Rapid-drawdown Stage-2 strength pipeline, against the worked example in
        # Duncan/Wright/Brandon Table 9.2. File-less: an infinite slope, hand-computable.
        tests.append({'type': 'drawdown_tauff', 'file': 'DWB Table 9.2 (unit)',
                      'method': 'tau_ff', 'source': 'drawdown_tauff'})
        # And the guard that stops eq (5) extrapolating past Kc=Kf when Stage-1 FS < 1.
        if Path(DRAWDOWN_GUARD_FILE).exists():
            tests.append({'type': 'drawdown_guard', 'file': DRAWDOWN_GUARD_FILE,
                          'method': 'stage1 FS>=1', 'source': 'drawdown_guard'})

    # Excel round-trip tests (save_slope_data_to_xlsx). Built from a curated file
    # list rather than markdown tags, since they check load/save fidelity, not FS.
    if run_roundtrip:
        # Guard that the packaged copies (shipped in the wheel) haven't drifted
        # from their editable docs masters.
        tests.append({'type': 'template_sync', 'file': BUNDLED_TEMPLATE,
                      'method': '-', 'source': 'template'})
        # Guard that every module-scope third-party import is a declared runtime dep,
        # so `pip install xslope` yields a package that actually imports.
        tests.append({'type': 'deps_declared', 'file': 'pyproject.toml',
                      'method': '-', 'source': 'packaging'})
        tests.append({'type': 'template_sync', 'file': BUNDLED_SKILL,
                      'master': SKILL_MASTER, 'copy': BUNDLED_SKILL,
                      'method': '-', 'source': 'skill'})
        # Editor no-drop guard (studio.editors). Touches the studio/Qt layer, so
        # it's skipped cleanly when PySide6 is absent (engine-only installs), like
        # the DXF round-trip tests.
        try:
            import PySide6          # noqa: F401
            tests.append({'type': 'editor_roundtrip', 'file': '(studio editors)',
                          'method': '-', 'source': 'editor_roundtrip'})
            if not run_all:
                print("Including 1 editor round-trip guard")
        except ImportError:
            print("Skipping editor round-trip guard (PySide6 not installed)")
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

    # The GeoStudio import test authors its own .gsz fixture (GeoStudio's own sample
    # files are Seequent's copyrighted Materials and are not in this repository), so it
    # needs no input file and no GUI — engine-only, always runnable.
    if run_gsz:
        tests.append({'type': 'gsz', 'file': '(synthetic .gsz)',
                      'method': '-', 'source': 'gsz'})
        if not run_all:
            print("Including 1 GeoStudio import test")

    # The Slide2 import test authors its own .slim fixture (Rocscience's own tutorial
    # files are their copyrighted material and are not in this repository), so it needs
    # no input file and no GUI — engine-only, always runnable.
    if run_slide2:
        tests.append({'type': 'slide2', 'file': '(synthetic .slim)',
                      'method': '-', 'source': 'slide2'})
        if not run_all:
            print("Including 1 Slide2 import test")

    # The RS2 import test authors its own .fez fixture (Rocscience's own verification
    # files are their copyrighted material and are not in this repository), so it needs
    # no input file and no GUI — engine-only, always runnable.
    if run_rs2:
        tests.append({'type': 'rs2', 'file': '(synthetic .fez)',
                      'method': '-', 'source': 'rs2'})
        if not run_all:
            print("Including 1 RS2 import test")

    if args.skip_benchmarks:
        n_before = len(tests)
        tests = [t for t in tests if 'benchmark' not in t]
        n_skipped = n_before - len(tests)
        if n_skipped:
            print(f"Skipping {n_skipped} benchmark tests (--skip-benchmarks)")

    if args.benchmark:
        wanted = [w.strip().lower() for w in args.benchmark.split(',') if w.strip()]

        def _matches(bench_id):
            b = str(bench_id).lower()
            for w in wanted:
                # exact id, or a family prefix followed by a non-digit
                # (VP28 matches VP28a and VP28-beta but not VP280; VP3 does
                # not match VP33)
                if b == w or (b.startswith(w) and not b[len(w):][:1].isdigit()):
                    return True
            return False

        tests = [t for t in tests if 'benchmark' in t and _matches(t['benchmark'])]
        print(f"Running {len(tests)} tests matching --benchmark {args.benchmark}")
        if not tests:
            print("No benchmark tags matched.")
            return 1

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
    wall_t0 = time.time()

    def report(i, test, computed, error_msg, elapsed):
        nonlocal passed, failed, errors
        file_name = Path(test['file']).name
        test_type = test.get('type', '?')
        method = test.get('method', '-')
        expected, tol = _expected_and_tol(test, args.tolerance)
        if error_msg:
            status = f"ERROR: {error_msg}"
            errors += 1
            comp_str = "    --    "
        elif expected is None:
            # A tag with no expected value for its type (e.g. a typo'd key or a
            # type missing from the pass/fail list above) is a broken tag, not
            # a pass — surface it instead of crashing the summary.
            status = "ERROR: tag has no expected value for this test type"
            errors += 1
            comp_str = f"{computed:10.3f}" if computed is not None else "    --    "
        elif abs(computed - expected) <= tol:
            status = f"PASS ({elapsed:.1f}s)"
            passed += 1
            comp_str = f"{computed:10.3f}"
        else:
            diff = computed - expected if computed is not None else float('nan')
            status = f"FAIL (diff={diff:+.4f}, {elapsed:.1f}s)"
            failed += 1
            comp_str = f"{computed:10.3f}" if computed is not None else "    --    "
        exp_str = f"{expected:10.3f}" if expected is not None else "    --    "
        print(f"{i:<4} {file_name:<45} {test_type:<20} {method:<10} {exp_str}  {comp_str}  {status}",
              flush=True)

    jobs = os.cpu_count() - 2 if str(args.jobs) == 'auto' else int(args.jobs)
    jobs = max(1, min(jobs, len(tests) or 1))

    if jobs > 1:
        # Keep numpy's BLAS single-threaded inside each worker so N processes
        # don't oversubscribe the CPU (inherited by spawned children).
        for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                   'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
            os.environ.setdefault(_v, '1')
        from concurrent.futures import ProcessPoolExecutor, as_completed
        indexed = list(enumerate(tests, 1))
        # slowest types first: wall time is set by the last straggler
        indexed.sort(key=lambda it: -(_COST_RANK.get(it[1].get('type'), 1)
                                      + (2 if it[1].get('rapid') else 0)))
        print(f"Running {len(tests)} tests on {jobs} workers "
              f"(rows print as they finish)")
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            futures = [ex.submit(_parallel_worker, item) for item in indexed]
            for fut in as_completed(futures):
                i, computed, error_msg, elapsed = fut.result()
                total_time += elapsed
                report(i, tests[i - 1], computed, error_msg, elapsed)
    else:
        for i, test in enumerate(tests, 1):
            t0 = time.time()
            try:
                computed, error_msg = run_test(test)
            except Exception as e:
                computed = None
                error_msg = str(e)
            elapsed = time.time() - t0
            total_time += elapsed
            report(i, test, computed, error_msg, elapsed)

    print("-" * 120)
    print(f"\nResults: {passed} passed, {failed} failed, {errors} errors out of {len(tests)} tests")
    if jobs > 1:
        print(f"Total time: {time.time() - wall_t0:.1f}s wall on {jobs} workers "
              f"({total_time:.1f}s of test time)")
    else:
        print(f"Total time: {total_time:.1f}s")

    if failed > 0 or errors > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
