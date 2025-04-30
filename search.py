import numpy as np
from slice import generate_slices
from shapely.geometry import LineString
import time

def circular_search(data, solver, tol=1e-2, max_iter=50, shrink_factor=0.5,
                    fs_fail=9999, depth_tol_frac=0.03, diagnostic=False):
    """
    Global 9-point circular search with adaptive grid refinement.

    Returns:
        list of dict: sorted fs_cache by FS
        bool: convergence flag
        list of dict: search path
    """

    start_time = time.time()  # Start timing

    ground_surface = data['ground_surface']
    ground_surface = LineString([(x, y) for x, y in ground_surface.coords])
    y_max = max(y for _, y in ground_surface.coords)
    y_min = data['max_depth']
    delta_y = y_max - y_min
    tol = delta_y * depth_tol_frac

    circles = data['circles']
    max_depth = data['max_depth']

    def optimize_depth(x, y, depth_guess, depth_step_init, depth_shrink_factor, tol_frac, fs_fail, diagnostic=False):
        depth_step = min(10.0, depth_step_init)
        best_depth = max(depth_guess, max_depth)
        best_fs = fs_fail
        best_result = None
        depth_tol = depth_step * tol_frac
        iterations = 0

        while depth_step > depth_tol:
            depths = [
                max(best_depth - depth_step, max_depth),
                best_depth,
                best_depth + depth_step
            ]
            fs_results = []
            for d in depths:
                test_circle = {'Xo': x, 'Yo': y, 'Depth': d, 'R': y - d}
                success, result = generate_slices(data, circle=test_circle)
                if not success:
                    FS = fs_fail
                    df_slices = None
                    failure_surface = None
                    solver_result = None
                else:
                    df_slices, failure_surface = result
                    solver_success, solver_result = solver(df_slices)
                    FS = solver_result['FS'] if solver_success else fs_fail
                fs_results.append((FS, d, df_slices, failure_surface, solver_result))

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

    def evaluate_grid(x0, y0, grid_size, depth_guess, data, diagnostic=False, fs_cache=None):
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
                diagnostic=diagnostic
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
    for i, start_circle in enumerate(circles):
        x0 = start_circle['Xo']
        y0 = start_circle['Yo']
        r0 = y0 - start_circle['Depth']
        if diagnostic:
            print(f"\n[⏱ starting circle {i+1}] x={x0:.2f}, y={y0:.2f}, r={r0:.2f}")
        grid_size = r0 * 0.15
        depth_guess = start_circle['Depth']
        fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, data, diagnostic=diagnostic)
        all_starts.append((start_circle, best_point, fs_cache))

    all_starts.sort(key=lambda t: t[1]['FS'])
    start_circle, best_start, fs_cache = all_starts[0]
    x0 = best_start['Xo']
    y0 = best_start['Yo']
    depth_guess = best_start['Depth']
    grid_size = (y0 - depth_guess) * 0.15
    best_fs = best_start['FS']

    # Include initial jump from user-defined circle to best point on its grid
    search_path = [
        {"x": start_circle['Xo'], "y": start_circle['Yo'], "FS": None},
        {"x": x0, "y": y0, "FS": best_fs}
    ]
    converged = False

    if diagnostic:
        print(f"\n[✅ launch grid] Starting refinement from FS={best_fs:.4f} at ({x0:.2f}, {y0:.2f})")

    for iteration in range(max_iter):
        print(f"[🔁 iteration {iteration+1}] center=({x0:.2f}, {y0:.2f}), FS={best_fs:.4f}, grid={grid_size:.4f}")
        fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, data, diagnostic=diagnostic, fs_cache=fs_cache)

        if best_point['FS'] < best_fs:
            best_fs = best_point['FS']
            x0 = best_point['Xo']
            y0 = best_point['Yo']
            depth_guess = best_point['Depth']
            search_path.append({"x": x0, "y": y0, "FS": best_fs})
        else:
            grid_size *= shrink_factor

        if grid_size < tol:
            converged = True
            end_time = time.time()
            elapsed = end_time - start_time
            print(f"[✅ converged] Iter={iteration+1}, FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f}), elapsed time={elapsed:.2f} seconds")
            break

    if not converged and diagnostic:
        print(f"\n[❌ max iterations reached] FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f})")

    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    return sorted_fs_cache, converged, search_path
