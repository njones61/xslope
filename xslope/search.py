# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import time

import numpy as np
from shapely.geometry import LineString, Point

from . import solve
from .advanced import rapid_drawdown, validate_rapid_drawdown
from .slice import generate_slices, get_y_from_intersection

# A valid limit-equilibrium failure surface must have a meaningful net gravitational
# driving force. A surface under flat ground (e.g. a circle in a reservoir bottom)
# has |sum W sin(alpha)| ~ 0 — no failure mechanism — yet with high pore pressure the
# factor of safety can blow up to a huge negative/near-zero value and win the search
# minimum. Reject surfaces whose |sum W sin(alpha)| is below this fraction of the
# total weight. Valid criticals run 0.2-0.6; degenerate flat circles run ~5e-5.
MIN_DRIVING_FRAC = 0.01

# A surface with more than this fraction of slices in base tension (negative
# effective normal) is rejected as a non-failure surface. Valid criticals run
# 0-10%; the degenerate high-pore-pressure circles run 40-80%.
MAX_BASE_TENSION_FRAC = 0.25


def _net_driving_too_small(df_slices):
    """True if the surface has negligible net gravitational driving force (a flat,
    non-failure surface) and should be rejected as a search candidate."""
    # Degenerate circles (e.g. grazing a stepped wall face) can yield an empty
    # slice table with no columns: reject those outright.
    if df_slices is None or len(df_slices) == 0 or 'w' not in df_slices.columns:
        return True
    W = df_slices['w'].values
    total = float(np.abs(W).sum())
    if total <= 0:
        return True
    driving = abs(float((W * np.sin(np.radians(df_slices['alpha'].values))).sum()))
    return driving < MIN_DRIVING_FRAC * total


MIN_THICKNESS_FRAC = 0.05


def _too_thin(df_slices, h_ground, frac):
    """True if the sliding mass is thinner than ``frac`` of the slope height.

    GRID MODE ONLY (frac=0 disables). A global sweep of a benched geometry finds
    surficial skins hugging any steep face — on a c~0 face they undercut every
    real mechanism. They are usually a raveling mode, not the slope failure the
    search is after, and Slide2 suppresses them the same way (its 'minimum
    depth' filter). But thin criticals can also be THE answer (a submerged c=0
    slope fails as an infinite-slope skin at zero depth), which is why this must
    never apply to the default seeded search: there the user's circles define
    the mechanism of interest, thin or not."""
    if frac <= 0:
        return False
    if df_slices is None or len(df_slices) == 0 or 'y_ct' not in df_slices.columns:
        return True
    thick = float((df_slices['y_ct'] - df_slices['y_cb']).max())
    return thick < frac * h_ground


def _too_shallow(df_slices, min_slip_depth):
    """True if the deepest point of the slip surface lies less than
    ``min_slip_depth`` below the ground surface — an ABSOLUTE minimum-slip-depth
    filter (Slide2's 'minimum depth'), and the search-side twin of solve_fem's
    ``min_slip_depth`` in SSRM. ``None`` or a non-positive value disables it.

    Unlike ``_too_thin`` (a fraction of slope height, grid-mode only), this
    applies to EVERY trial when set, so a caller can suppress surficial
    cohesionless skins (FS = tan phi / tan beta) on demand while leaving the
    default seeded search untouched."""
    if min_slip_depth is None or min_slip_depth <= 0:
        return False
    if df_slices is None or len(df_slices) == 0 or 'y_ct' not in df_slices.columns:
        return True
    max_depth = float((df_slices['y_ct'] - df_slices['y_cb']).max())
    return max_depth < min_slip_depth


def _base_tension_too_extensive(df_slices):
    """True if a large fraction of slices carry a negative effective base normal
    (base in tension) — a non-physical surface, e.g. a circle through a
    high-pore-pressure zone where the moment methods then return a negative/near-zero
    factor of safety. Valid criticals run ~0%; the degenerate cases run 60-80%. A
    few negative base normals are normal and are NOT rejected (extent only)."""
    if 'n_eff' not in df_slices.columns:
        return False
    N = df_slices['n_eff'].values
    return len(N) > 0 and float((N < 0).mean()) > MAX_BASE_TENSION_FRAC


class AnalysisCancelled(Exception):
    """Raised by a search/reliability run when its ``cancel_check`` asks it to stop.

    Cooperative cancellation: long loops call ``cancel_check()`` at iteration
    boundaries and raise this if it returns True, so a caller (e.g. the GUI's
    worker thread) can abort cleanly without killing the thread mid-computation."""


def _check_cancel(cancel_check):
    if cancel_check is not None and cancel_check():
        raise AnalysisCancelled()


def _grid_seed_circles(slope_data, method_name, num_slices=20, fs_fail=9999,
                       rapid=False, composite=False, nx=10, ny=5, n_tangents=6,
                       keep=4, diagnostic=False, cancel_check=None, circle_cache=None,
                       min_slip_depth=None):
    """Coarse grid-and-tangent sweep that SEEDS the adaptive circular search.

    The adaptive 9-point search is a local optimizer: it refines whatever
    neighborhood its starting circles put it in, and a single seed in the wrong
    family converges to a local minimum with no warning (a shallow circle in a
    fill can read 20%+ high while the true critical rides the base of the
    foundation). This sweep is the global stage that protects against that.

    It is Slide2's grid-and-tangent model. Centers are laid out on an nx-by-ny
    grid over a box derived from the slope geometry — spanning the slope
    horizontally (extended 0.75H past the toe and crest) and sitting 0.25H to
    2H above the crest, which brackets where the center of a circle that cuts
    the slope can be. For each center, one trial circle is made TANGENT to each
    of ``n_tangents`` elevations between the bottom of the model and mid-slope.
    Sweeping tangent elevations, not depths-below-center, is what finds
    competing families: the shallow face circle and the deep base-tangent
    circle live at similar centers but different tangent lines, and a local
    depth optimizer started in one family never sees the other.

    Every trial is solved (at a reduced slice count) with the same degenerate-
    surface guards as the search proper. The best circle from each tangent
    family within 25% of the global best FS is returned as a seed, deepest
    families first among ties, capped at ``keep``. The caller merges these with
    any user-entered circles and runs the normal adaptive search from all of
    them.
    """
    solver = getattr(solve, method_name)
    coords = list(slope_data['ground_surface'].coords)
    ys_g = [y for _, y in coords]
    y_hi, y_lo = max(ys_g), min(ys_g)
    H = y_hi - y_lo
    if H <= 0:
        return []
    right_facing = coords[0][1] > coords[-1][1]
    tol_y = 1e-6 * max(1.0, H)
    his = [x for x, y in coords if y > y_hi - tol_y]
    los = [x for x, y in coords if y < y_lo + tol_y]
    x_crest = max(his) if right_facing else min(his)
    x_toe = min(los) if right_facing else max(los)

    # The box must scale with the FULL model depth, not just the slope height: a
    # circle tangent to the base of a deep foundation spans the whole model, and
    # its center sits several slope-heights above the crest (VP75's critical
    # centers run to ~2.5x the crest-to-floor depth above the toe).
    depth_floor = slope_data['domain_polygon'].bounds[1]
    Hm = max(H, y_hi - depth_floor)
    xs = np.linspace(min(x_toe, x_crest) - 0.75 * Hm, max(x_toe, x_crest) + 0.75 * Hm, nx)
    ys = np.linspace(y_hi + 0.25 * Hm, y_hi + 2.5 * Hm, ny)
    y_mid = 0.5 * (y_lo + y_hi)
    tangents = list(np.linspace(depth_floor, y_mid, n_tangents))
    if composite:
        # one below-floor family: its circles truncate at the base and run along it
        tangents.insert(0, depth_floor - 0.5 * H)

    if rapid:
        from .advanced import rapid_drawdown

    sweep_slices = min(num_slices, 20)
    best_by_tan = {}
    n_valid = 0
    for yt in tangents:
        for yc in ys:
            _check_cancel(cancel_check)
            for xc in xs:
                R = yc - yt
                if R <= 0:
                    continue
                circle = {'Xo': float(xc), 'Yo': float(yc), 'Depth': float(yt), 'R': float(R)}
                success, result = generate_slices(slope_data, circle=circle,
                                                  num_slices=sweep_slices,
                                                  composite=composite)
                if not success:
                    continue
                df_slices, failure_surface = result
                if (_net_driving_too_small(df_slices) or _too_thin(df_slices, H, MIN_THICKNESS_FRAC)
                        or _too_shallow(df_slices, min_slip_depth)):
                    continue
                if rapid:
                    ok, solver_result = rapid_drawdown(df_slices, method_name, debug_level=0)
                else:
                    ok, solver_result = solver(df_slices)
                if not ok:
                    continue
                FS = solver_result['FS']
                if not (0 < FS < fs_fail) or _base_tension_too_extensive(df_slices):
                    continue
                n_valid += 1
                if circle_cache is not None:
                    circle_cache.append({"Xo": circle['Xo'], "Yo": circle['Yo'],
                                         "Depth": circle['Depth'], "R": circle['R'],
                                         "FS": FS, "failure_surface": failure_surface})
                if yt not in best_by_tan or FS < best_by_tan[yt][0]:
                    best_by_tan[yt] = (FS, circle)

    if not best_by_tan:
        return []
    fs_best = min(fs for fs, _ in best_by_tan.values())
    ranked = sorted(best_by_tan.items(), key=lambda kv: kv[1][0])
    seeds = [circ for yt, (fs, circ) in ranked if fs <= 1.25 * fs_best][:keep]
    print(f"[grid seed] swept {len(tangents) * nx * ny} circles ({n_valid} valid), "
          f"best FS={fs_best:.4f}; seeding {len(seeds)} tangent "
          f"famil{'y' if len(seeds) == 1 else 'ies'}")
    if diagnostic:
        for fs, c in [v for _, v in ranked]:
            tag = 'seed' if c in seeds else '    '
            print(f"  [{tag}] tangent y={c['Depth']:8.2f}  FS={fs:7.4f}  "
                  f"center=({c['Xo']:.1f}, {c['Yo']:.1f})")
    return seeds


def circular_search(slope_data, method_name, rapid=False, tol=1e-2, fs_tol=5e-4, max_iter=50,
                    shrink_factor=0.5, fs_fail=9999, min_grid_frac=0.03, depth_tol_frac=0.03,
                    diagnostic=False, num_slices=40, cancel_check=None, composite=False,
                    seed='circles', min_slip_depth=None):
    """
    Global 9-point circular search with adaptive grid refinement.

    Convergence is driven by the factor of safety, not the grid spacing: the
    center grid keeps refining until the best FS stops improving — specifically
    until two successive refinement levels each gain less than ``fs_tol`` (so the
    reported FS is stable to ~3 decimal places). FS convergence is only accepted
    once the center grid has refined below ``min_grid_frac`` of the domain height,
    so a coarse-grid FS plateau cannot be mistaken for the true minimum. This grid
    gate is deliberately loose (3% of the domain height): near the critical circle
    the FS surface is flat, so once the FS has plateaued the center is already
    located well enough — tightening the gate only spends extra refinement levels
    that don't move the reported FS. ``tol`` is a geometric backstop on the grid
    spacing and ``max_iter`` caps the count. Keying the stop on FS (not an absolute
    length like ``grid < 0.01``, which means different things on a 20 ft vs a
    500 ft slope) makes it scale-invariant.

    ``seed`` selects the starting circles. ``'circles'`` (default) refines from
    the circles sheet, exactly as before. ``'grid'`` first runs a coarse
    grid-and-tangent sweep over a box derived from the slope geometry and seeds
    the refinement with the best circle from each competing tangent family (plus
    any user circles) — see _grid_seed_circles. Use it when a single starting
    circle risks trapping the search in a local minimum, e.g. an embankment on a
    layered foundation where a shallow fill circle and a deep foundation circle
    compete.

    ``composite`` (default False) lets trial circles dip below the bottom of the
    model, where they are truncated into a composite surface that runs along the
    base (see slice.CompositeSurface). Turn it on when the base is a real
    impenetrable boundary with a weak seam or soft layer on it — the critical
    mechanism there rides the base, and a search clamped to the floor can only ever
    reach a circle tangent to it. Leave it off when the floor is just ``max_depth``,
    a search bound rather than bedrock.

    ``min_slip_depth`` (default None = off) is an absolute minimum-slip-depth filter:
    trial surfaces whose deepest point is less than this far below the ground surface
    are rejected, so a shallow surficial skin (FS = tan phi / tan beta on a
    cohesionless face) cannot win. This is Slide2's "minimum depth" filter and the
    search-side twin of solve_fem's ``min_slip_depth`` in SSRM — set the same value in
    both to keep an LEM/SSRM comparison on the same mechanism. Off by default, so the
    reported surface is the true minimum unless the caller opts in.

    Returns:
        list of dict: sorted fs_cache by FS
        bool: convergence flag
        list of dict: search path
        list of dict: circle_cache - all circles tested during search
    """

    if seed != 'grid' and not slope_data.get('circles'):
        raise ValueError(
            "Circular search requires at least one circle defined in the input. "
            "Add circle data to the 'circles' sheet in the input template, "
            "or run with seed='grid' to generate starting circles automatically.")

    if rapid:
        validate_rapid_drawdown(slope_data)

    solver = getattr(solve, method_name)
    circle_cache = []  # Store ALL circles tested for plotting

    start_time = time.time()  # Start timing

    ground_surface = slope_data['ground_surface']
    ground_surface = LineString([(x, y) for x, y in ground_surface.coords])
    y_max = max(y for _, y in ground_surface.coords)
    h_ground = y_max - min(y for _, y in ground_surface.coords)   # slope height, for the skin-slide guard
    # The deepest allowable failure-surface elevation is the bottom of the domain
    # polygon (not max_depth, which is only a profile-sheet input). Circles below an
    # irregular bottom are rejected by the containment check in generate_slices.
    depth_floor = slope_data['domain_polygon'].bounds[1]
    y_min = depth_floor
    delta_y = y_max - y_min

    # composite=True lets trial circles dip BELOW the floor, where generate_slices
    # truncates them into a composite surface that runs along the floor (see
    # slice.CompositeSurface). This is what a weak seam or a soft layer sitting on
    # bedrock actually fails along, and a search clamped to the floor can never find
    # it: the best it can do is a circle tangent to the base. Off by default,
    # because on a profile-line model the floor is `max_depth` — a search bound the
    # user picked, not a real impenetrable boundary, and truncating there would be
    # meaningless. The lower bound is generous; a circle that goes too deep flattens
    # out along the floor and its FS climbs again, so the optimizer turns back.
    search_floor = depth_floor - delta_y if composite else depth_floor
    min_grid = delta_y * min_grid_frac   # required center-grid resolution before FS convergence

    # The skin-slide filter runs only in grid mode: a global sweep must not report
    # a raveling sliver as the critical, but the seeded search must stay able to
    # follow a user's circles to a genuinely thin critical (e.g. a submerged c=0
    # slope, whose infinite-slope skin IS the answer). See _too_thin.
    thin_frac = MIN_THICKNESS_FRAC if seed == 'grid' else 0.0

    if seed == 'grid':
        # Global stage: sweep a geometry-derived center grid against tangent
        # elevations and seed the adaptive search with the best circle from each
        # competing family (see _grid_seed_circles). User circles, if any, stay in
        # as extra seeds. This is the protection against the local-minimum trap —
        # a single user seed in the wrong family otherwise converges 20%+ high
        # with no warning.
        _check_cancel(cancel_check)
        circles = _grid_seed_circles(
            slope_data, method_name, num_slices=num_slices, fs_fail=fs_fail,
            rapid=rapid, composite=composite, diagnostic=diagnostic,
            cancel_check=cancel_check, circle_cache=circle_cache,
            min_slip_depth=min_slip_depth)
        circles = circles + list(slope_data.get('circles') or [])
        if not circles:
            return [], False, [], circle_cache
    else:
        circles = slope_data['circles']

    def optimize_depth(x, y, depth_guess, depth_step_init, depth_shrink_factor, tol_frac, fs_fail, circle_cache, diagnostic=False):
        depth_step = min(10.0, depth_step_init)
        best_depth = max(depth_guess, search_floor)
        best_fs = fs_fail
        # These three are only rebound inside the loop below, so they must be seeded
        # here: a degenerate depth_step (<= depth_tol on entry) skips the loop
        # entirely, and the return would otherwise raise UnboundLocalError instead of
        # reporting an honest fs_fail.
        best_df = None
        best_surface = None
        best_solver_result = None
        depth_tol = depth_step * tol_frac
        iterations = 0

        while depth_step > depth_tol:
            depths = [
                max(best_depth - depth_step, search_floor),
                best_depth,
                best_depth + depth_step
            ]
            fs_results = []
            for d in depths:
                test_circle = {'Xo': x, 'Yo': y, 'Depth': d, 'R': y - d}
                success, result = generate_slices(slope_data, circle=test_circle,
                                                  num_slices=num_slices, composite=composite)
                if not success:
                    FS = fs_fail
                    df_slices = None
                    failure_surface = None
                    solver_result = None
                else:
                    df_slices, failure_surface = result
                    if (_net_driving_too_small(df_slices) or _too_thin(df_slices, h_ground, thin_frac)
                            or _too_shallow(df_slices, min_slip_depth)):
                        FS = fs_fail  # degenerate (near-zero driving; grid mode also drops skins)
                        solver_result = None
                    else:
                        if rapid:
                            solver_success, solver_result = rapid_drawdown(df_slices, method_name, debug_level=0)
                        else:
                            solver_success, solver_result = solver(df_slices)
                        FS = solver_result['FS'] if solver_success else fs_fail
                        if solver_success and _base_tension_too_extensive(df_slices):
                            FS = fs_fail  # degenerate surface (base mostly in tension)
                        if not (FS > 0):  # reject non-positive / NaN factor of safety
                            FS = fs_fail
                fs_results.append((FS, d, df_slices, failure_surface, solver_result))

                # Add to circle_cache for plotting all tested circles
                if FS != fs_fail:
                    circle_cache.append({
                        "Xo": x,
                        "Yo": y,
                        "Depth": d,
                        "R": y - d,
                        "FS": FS,
                        "failure_surface": failure_surface
                    })

            fs_results.sort(key=lambda t: t[0])
            best_fs, best_depth, best_df, best_surface, best_solver_result = fs_results[0]

            if all(FS == fs_fail for FS, *_ in fs_results):
                if diagnostic:
                    print(f"[❌ all fail] x={x:.2f}, y={y:.2f}")
                return best_depth, fs_fail, None, None, None

            if diagnostic:
                print(f"[✓ depth-opt] x={x:.2f}, y={y:.2f}, depth={best_depth:.2f}, FS={best_fs:.4f}, step={depth_step:.2f}")

            depth_step *= depth_shrink_factor
            iterations += 1
            if iterations > 50:
                if diagnostic:
                    print(f"[⚠️ warning] depth iterations exceeded at (x={x:.2f}, y={y:.2f})")
                break

        return best_depth, best_fs, best_df, best_surface, best_solver_result

    def evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=False, fs_cache=None, circle_cache=None):
        if fs_cache is None:
            fs_cache = {}

        Xs = [x0 - grid_size, x0, x0 + grid_size]
        Ys = [y0 - grid_size, y0, y0 + grid_size]
        points = [(x, y) for y in Ys for x in Xs]

        for i, (x, y) in enumerate(points):
            if (x, y) in fs_cache:
                result = fs_cache[(x, y)]
                if diagnostic:
                    print(f"[cache hit] grid pt {i + 1}/9 at (x={x:.2f}, y={y:.2f}) → FS={result['FS']:.4f}")
                continue

            depth_step_init = grid_size * 0.75
            d, FS, df_slices, failure_surface, solver_result = optimize_depth(
                x, y, depth_guess, depth_step_init, depth_shrink_factor=0.25, tol_frac=0.01, fs_fail=fs_fail,
                circle_cache=circle_cache, diagnostic=diagnostic
            )

            fs_cache[(x, y)] = {
                "Xo": x,
                "Yo": y,
                "Depth": d,
                "FS": FS,
                "slices": df_slices,
                "failure_surface": failure_surface,
                "solver_result": solver_result
            }

            if diagnostic:
                print(f"[grid pt {i + 1}/9] x={x:.2f}, y={y:.2f} → FS={FS:.4f} at d={d:.2f}")

        sorted_fs = sorted(fs_cache.items(), key=lambda item: item[1]['FS'])
        best_point = sorted_fs[0][1]
        best_index = list(fs_cache.keys()).index((best_point['Xo'], best_point['Yo']))

        if diagnostic:
            print(f"[★ grid best {best_index + 1}/9] FS={best_point['FS']:.4f} at (x={best_point['Xo']:.2f}, y={best_point['Yo']:.2f})")

        return fs_cache, best_point

    # === Step 1: Evaluate starting circles ===
    all_starts = []
    fs_cache = {}  # Shared cache for all starting circles
    for i, start_circle in enumerate(circles):
        _check_cancel(cancel_check)
        x0 = start_circle['Xo']
        y0 = start_circle['Yo']
        r0 = y0 - start_circle['Depth']
        if diagnostic:
            print(f"\n[⏱ starting circle {i+1}] x={x0:.2f}, y={y0:.2f}, r={r0:.2f}")
        grid_size = r0 * 0.15
        depth_guess = start_circle['Depth']
        fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=diagnostic, fs_cache=fs_cache, circle_cache=circle_cache)
        all_starts.append((start_circle, best_point))

    all_starts.sort(key=lambda t: t[1]['FS'])

    def refine(start_circle, best_start, launch_label=""):
        """Adaptive 9-point refinement from one start. Returns (fs, converged, path)."""
        nonlocal fs_cache
        x0 = best_start['Xo']
        y0 = best_start['Yo']
        depth_guess = best_start['Depth']
        grid_size = (y0 - depth_guess) * 0.15
        best_fs = best_start['FS']

        # Include initial jump from the seed circle to the best point on its grid
        path = [
            {"x": start_circle['Xo'], "y": start_circle['Yo'], "FS": None},
            {"x": x0, "y": y0, "FS": best_fs}
        ]
        converged = False

        if diagnostic:
            print(f"\n[✅ launch grid{launch_label}] Starting refinement from FS={best_fs:.4f} at ({x0:.2f}, {y0:.2f})")

        fs_level_start = best_fs   # best FS at the start of the current grid resolution
        small_levels = 0           # consecutive refinements that barely improved FS

        for iteration in range(max_iter):
            _check_cancel(cancel_check)
            print(f"[🔁 iteration {iteration+1}{launch_label}] center=({x0:.2f}, {y0:.2f}), FS={best_fs:.4f}, grid={grid_size:.4f}")
            fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=diagnostic, fs_cache=fs_cache, circle_cache=circle_cache)

            if best_point['FS'] < best_fs:
                # Found a better center at this resolution — move there and keep exploring.
                best_fs = best_point['FS']
                x0 = best_point['Xo']
                y0 = best_point['Yo']
                depth_guess = best_point['Depth']
                path.append({"x": x0, "y": y0, "FS": best_fs})
                continue

            # No improvement at this resolution. Refine the grid, and converge once the
            # FS gain over a refinement level has been negligible twice in a row (so the
            # third decimal of FS is stable). The grid floor `tol` is only a backstop.
            level_gain = fs_level_start - best_fs
            if level_gain < fs_tol and grid_size < min_grid:
                small_levels += 1
                if small_levels >= 2:
                    converged = True
                    elapsed = time.time() - start_time
                    print(f"[✅ converged] Iter={iteration+1}, FS={best_fs:.4f} (ΔFS<{fs_tol}) at (x={x0:.2f}, y={y0:.2f}, depth={depth_guess:.2f}), elapsed time={elapsed:.2f} seconds")
                    break
            else:
                small_levels = 0
            fs_level_start = best_fs
            grid_size *= shrink_factor

            if grid_size < tol:
                converged = True
                elapsed = time.time() - start_time
                print(f"[✅ converged: grid floor] Iter={iteration+1}, FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f}, depth={depth_guess:.2f}), elapsed time={elapsed:.2f} seconds")
                break

        if not converged and diagnostic:
            print(f"\n[❌ max iterations reached] FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f})")
        return best_fs, converged, path

    # Grid seeding refines EVERY surviving family, not just the best start —
    # competing families can differ by only a few percent at the coarse stage yet
    # refine to very different minima, and launching only the best start is
    # exactly the local-minimum trap the grid exists to remove. Seeded mode
    # ('circles') keeps the original single-launch behavior unchanged.
    if seed == 'grid':
        fs_cut = all_starts[0][1]['FS'] * 1.15
        launches = [t for t in all_starts if t[1]['FS'] <= fs_cut][:4]
    else:
        launches = all_starts[:1]

    search_path = []
    best_fs, converged = None, False
    for i, (sc, bs) in enumerate(launches):
        label = f" (family {i+1}/{len(launches)})" if len(launches) > 1 else ""
        fs_i, conv_i, path_i = refine(sc, bs, launch_label=label)
        search_path.extend(path_i)
        if best_fs is None or fs_i < best_fs:
            best_fs, converged = fs_i, conv_i

    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    return sorted_fs_cache, converged, search_path, circle_cache

def noncircular_search(slope_data, method_name, rapid=False, diagnostic=True, movement_distance=4.0, shrink_factor=0.8, fs_tol=0.001, max_iter=100, move_tol=0.1, num_slices=30, max_base_angle=65.0, cancel_check=None, min_slip_depth=None):
    """
    Non-circular search using the specified solver.

    A geometric admissibility guard (`max_base_angle`) rejects trial surfaces with
    an over-steep base segment. Without it the coordinate-descent search drives the
    points into a near-vertical base running up to the toe, which the rigorous
    methods (Spencer especially) score as a spurious low minimum. Capping the base
    inclination keeps the search on physically realistic surfaces. (The convergence
    criterion itself is kept absolute, but the guard removes the degeneracy that
    blocked refining it.)

    Parameters:
    -----------
    data : dict
        Input data dictionary containing all necessary parameters
    method_name : str
        The method name to use (e.g., 'lowe', 'spencer')
    diagnostic : bool
        If True, print diagnostic information during search
    movement_distance : float
        Initial distance to move points in each iteration
    shrink_factor : float
        Factor to reduce movement_distance by when no improvement is found
    fs_tol : float
        Factor of safety convergence tolerance
    max_iter : int
        Maximum number of iterations
    move_tol : float
        Minimum movement distance for convergence (AND logic with fs_tol)
    max_base_angle : float
        Maximum allowed base inclination (degrees) for any slice. Trial surfaces
        with a steeper base are rejected as inadmissible. Default 65 deg, the
        active-wedge angle (45 + phi/2) for phi ~ 40 deg, near the steep end of
        realistic soils; steeper bases are geometrically unrealistic slip surfaces.

    Returns:
    --------
    tuple : (fs_cache, converged, search_path)
        fs_cache : dict of all evaluated surfaces and their FS values
        converged : bool indicating if search converged
        search_path : list of surfaces evaluated during search
    """
    if not slope_data.get('non_circ'):
        raise ValueError(
            "Non-circular search requires a non-circular surface defined in the "
            "input. Add surface point data to the 'circles' sheet (non-circular "
            "section) in the input template.")

    if rapid:
        validate_rapid_drawdown(slope_data)

    # Get the solver function from solve module
    solver = getattr(solve, method_name)
    def move_point(points, i, dx, dy, movement_type, ground_surface, depth_floor):
        """Move a point while respecting constraints"""
        # Get current point
        point = points[i]
        
        # Calculate new position
        new_x = point[0] + dx
        new_y = point[1] + dy
        
        # For endpoints, ensure they stay on ground surface
        if i == 0 or i == len(points)-1:
            # Create vertical line at new_x
            vertical_line = LineString([(new_x, 0), (new_x, 1000)])  # Arbitrary high y value
            intersection = ground_surface.intersection(vertical_line)
            y = get_y_from_intersection(intersection)
            if y is None:
                return False
            new_y = y
        else:
            # For middle points, ensure they stay below ground surface but above the
            # domain floor (surfaces leaving an irregular bottom are rejected later
            # by the containment check in generate_slices).
            if new_y > ground_surface.interpolate(ground_surface.project(Point(new_x, new_y))).y:
                return False
            if new_y < depth_floor:
                return False
        
        # Check x-ordering constraints
        if i > 0 and new_x <= points[i-1][0]:  # Don't move past left neighbor
            return False
        if i < len(points)-1 and new_x >= points[i+1][0]:  # Don't move past right neighbor
            return False
        
        # Update point
        points[i] = [new_x, new_y]
        return True

    def evaluate_surface(points, distance, fs_cache=None):
        """Evaluate factor of safety for current surface configuration"""
        if fs_cache is None:
            fs_cache = {}
            
        # Create non_circ format from points
        non_circ = [{'X': x, 'Y': y, 'Movement': movements[i]} for i, (x, y) in enumerate(points)]

        # Geometric admissibility on the surface itself: slicing breaks at
        # every vertex, so a steep segment with real width gets judged by the
        # per-slice alpha guard below - but generate_slices merges slices
        # thinner than its min_width (1e-2), so a near-vertical segment
        # narrower than that is absorbed into its neighbour and its angle
        # diluted out of the alpha table. Reject exactly those hidden risers.
        seg = np.diff(np.array([[p_['X'], p_['Y']] for p_ in non_circ]), axis=0)
        seg_ang = np.degrees(np.arctan2(np.abs(seg[:, 1]),
                                        np.maximum(np.abs(seg[:, 0]), 1e-300)))
        hidden_riser = (np.abs(seg[:, 0]) < 1e-2) & (seg_ang > max_base_angle)
        if hidden_riser.any():
            return float('inf'), None, None, None, fs_cache

        # Generate slices and compute FS
        success, result = generate_slices(slope_data, non_circ=non_circ, num_slices=num_slices)
        if not success:
            return float('inf'), None, None, None, fs_cache
            
        df_slices, failure_surface = result

        # Admissibility guard: reject surfaces with an over-steep base segment. A
        # near-vertical base (typically running up to the toe) is geometrically
        # unrealistic and is the spurious low-FS minimum the coordinate-descent
        # search would otherwise slide into (worst for Spencer/Lowe).
        if np.abs(df_slices['alpha'].values).max() > max_base_angle:
            return float('inf'), None, None, None, fs_cache

        # Reject surfaces with negligible net driving force (flat / non-failure).
        if _net_driving_too_small(df_slices):
            return float('inf'), None, None, None, fs_cache

        # Optional minimum-slip-depth filter: drop surficial skins on demand.
        if _too_shallow(df_slices, min_slip_depth):
            return float('inf'), None, None, None, fs_cache

        if rapid:
            solver_success, solver_result = rapid_drawdown(df_slices, method_name, debug_level=0)
        else:
            solver_success, solver_result = solver(df_slices)
        FS = solver_result['FS'] if solver_success else float('inf')
        if solver_success and _base_tension_too_extensive(df_slices):
            FS = float('inf')  # degenerate surface (base mostly in tension)
        if not (FS > 0):  # reject non-positive / NaN factor of safety
            FS = float('inf')
        
        # Cache result
        key = tuple(map(tuple, points))
        fs_cache[key] = {
            'points': points.copy(),
            'FS': FS,
            'slices': df_slices,
            'failure_surface': failure_surface,
            'solver_result': solver_result
        }
        
        return FS, df_slices, failure_surface, solver_result, fs_cache

    # Get initial surface from non_circ data
    non_circ = slope_data['non_circ']
    points = np.array([[p['X'], p['Y']] for p in non_circ])
    movements = [p['Movement'] for p in non_circ]
    ground_surface = slope_data['ground_surface']
    
    # Initialize cache and search path
    fs_cache = {}
    search_path = []
    
    # Evaluate initial surface
    FS, df_slices, failure_surface, solver_result, fs_cache = evaluate_surface(
        points, movement_distance, fs_cache)
    
    # A dead starting surface (slice generation or the solver failed on it)
    # leaves failure_surface=None; the search loop can't recover from that, so
    # report "no valid surface" instead of crashing on best_surface.coords.
    if failure_surface is None or not np.isfinite(FS):
        print("[❌ noncircular_search] the starting surface is not viable "
              "(slice generation or the solver failed on it) — adjust the "
              "non_circ starting points.")
        return [], False, []

    # Initialize best surface with initial evaluation
    best_points = points.copy()
    best_fs = FS
    best_df = df_slices
    best_surface = failure_surface
    best_solver_result = solver_result
    
    # Track convergence
    converged = False
    start_time = time.time()
    prev_fs = best_fs
    
    if diagnostic:
        print(f"\n[✅ starting search] Initial FS={best_fs:.4f}\n")
        print("Initial failure surface:")
        for i, point in enumerate(points):
            print(f"Point {i}: ({point[0]:.2f}, {point[1]:.2f})")
        print("\nGround surface:")
        for i, point in enumerate(ground_surface.coords):
            print(f"Point {i}: ({point[0]:.2f}, {point[1]:.2f})")
    
    # Main search loop
    for iteration in range(max_iter):
        _check_cancel(cancel_check)
        improved = False
        
        if diagnostic:
            print(f"\nIteration {iteration + 1}")
            print("Current surface points:")
            for i, point in enumerate(best_surface.coords):
                print(f"Point {i}: ({point[0]:.2f}, {point[1]:.2f})")
        
        # Try moving each point
        for i in range(len(points)):
            # Try both positive and negative directions
            for direction in [-1, 1]:
                test_points = points.copy()
                
                # Get movement direction based on point type
                if i == 0 or i == len(points)-1:  # End points
                    dx = direction * movement_distance
                    dy = 0  # y will be determined by ground surface
                elif movements[i] == 'Horiz':
                    dx = direction * movement_distance
                    dy = 0
                elif movements[i] == 'Free':
                    # For free points, move perpendicular to tangent
                    dx_tangent = points[i+1][0] - points[i-1][0] if i > 0 and i < len(points)-1 else 1
                    dy_tangent = points[i+1][1] - points[i-1][1] if i > 0 and i < len(points)-1 else 0
                    length = np.sqrt(dx_tangent**2 + dy_tangent**2)
                    if length > 0:
                        dx = -dy_tangent/length * direction * movement_distance
                        dy = dx_tangent/length * direction * movement_distance
                    else:
                        dx = direction * movement_distance
                        dy = 0
                else:  # Fixed
                    continue
                
                # Try to move the point
                if move_point(test_points, i, dx, dy, movements[i], ground_surface, slope_data['domain_polygon'].bounds[1]):
                    # Evaluate new surface
                    FS, df_slices, failure_surface, solver_result, fs_cache = evaluate_surface(
                        test_points, movement_distance, fs_cache)
                    
                    if FS < best_fs:
                        best_fs = FS
                        best_points = test_points.copy()
                        best_df = df_slices
                        best_surface = failure_surface
                        best_solver_result = solver_result
                        improved = True
                        if diagnostic:
                            print(f"[✓ improved] iter={iteration}, point={i}, FS={FS:.4f}")
        
        # print iteration results
        print(f"iteration {iteration+1} FS={best_fs:.4f}")
        
        # Check convergence based on FS change and movement distance (AND logic)
        fs_change = abs(best_fs - prev_fs)
        if fs_change < fs_tol and movement_distance < move_tol:
            converged = True
            if diagnostic:
                print(f"[✓ converged] FS change {fs_change:.6f} < tolerance {fs_tol} and movement_distance {movement_distance:.4f} < move_tol {move_tol}")
            break
        prev_fs = best_fs
        
        if not improved or fs_change < fs_tol:
            movement_distance *= shrink_factor
            if True:
                print(f"[↘️ shrinking] movement_distance={movement_distance:.4f}")
        
        points = best_points.copy()
        
    end_time = time.time()
    elapsed = end_time - start_time
    
    if converged:
        print(f"\n[✅ converged] Iter={iteration+1}, FS={best_fs:.4f}, elapsed time={elapsed:.2f} seconds")
    else:
        print(f"\n[❌ max iterations reached] FS={best_fs:.4f}, elapsed time={elapsed:.2f} seconds")
    
    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    return sorted_fs_cache, converged, search_path
