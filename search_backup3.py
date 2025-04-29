import numpy as np
from slice import generate_slices

def circular_search(data, solver, tol=1e-2, max_iter=50, shrink_factor=0.5, fs_fail=9999, depth_tol_frac=0.01, diagnostic=False):
    """
    Global 9-point circular search using multiple starting circles.

    Returns:
        list of dict: fs_cache sorted by ascending FS
        bool: convergence status
        list of dict: search path points (x, y, FS)
    """

    ground_surface = data['ground_surface']
    circles = data['circles']
    max_depth = data['max_depth']

    def optimize_depth(x, y, r0, depth_guess, data):
        # Start with smaller, capped depth step
        depth_step = min(10.0, r0 * 0.1)  # ⬅️ change from r0 * 0.25
        best_depth = max(depth_guess, data['max_depth'])
        best_fs = fs_fail
        best_result = None
        depth_tol = r0 * depth_tol_frac
        iterations = 0

        while depth_step > depth_tol:
            depths = [
                max(best_depth - depth_step, data['max_depth']),
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
                    if not solver_success:
                        FS = fs_fail
                        solver_result = None
                    else:
                        FS = solver_result['FS']
                fs_results.append((FS, d, df_slices, failure_surface, solver_result))

            fs_results.sort(key=lambda t: t[0])
            best_fs, best_depth, best_df, best_surface, best_solver_result = fs_results[0]

            if all(FS == fs_fail for FS, *_ in fs_results):
                if diagnostic:
                    print(f"[optimize_depth] All FS failed at (x={x:.2f}, y={y:.2f})")
                return best_depth, fs_fail, None, None, None

            if diagnostic:
                print(
                    f"[optimize_depth] x={x:.2f}, y={y:.2f}, depth={best_depth:.2f}, best FS={best_fs:.4f}, depth_step={depth_step:.4f}")

            depth_step *= shrink_factor
            iterations += 1

            if iterations > 50:
                if diagnostic:
                    print(f"[optimize_depth] Max depth iterations exceeded at (x={x:.2f}, y={y:.2f})")
                break

        return best_depth, best_fs, best_df, best_surface, best_solver_result

    def evaluate_grid(x0, y0, r0, prev_depths={}):
        Xs = [x0 - r0 * 0.25, x0, x0 + r0 * 0.25]
        Ys = [y0 - r0 * 0.25, y0, y0 + r0 * 0.25]
        points = [(x, y) for y in Ys for x in Xs]

        fs_cache = {}
        for x, y in points:
            depth_guess = prev_depths.get((x, y), y - r0)
            depth, FS, df_slices, failure_surface, solver_result = optimize_depth(x, y, r0, depth_guess, data)
            fs_cache[(x, y)] = {
                "Xo": x,
                "Yo": y,
                "Depth": depth,
                "FS": FS,
                "slices": df_slices,
                "failure_surface": failure_surface,
                "solver_result": solver_result
            }
            if diagnostic:
                print(f"[evaluate_grid] (x={x:.2f}, y={y:.2f}) FS={FS:.4f}, Depth={depth:.2f}")

        sorted_fs = sorted(fs_cache.items(), key=lambda item: item[1]['FS'])
        best_point = sorted_fs[0][1]

        return fs_cache, best_point

    # === Step 1: Evaluate grids for all starting circles ===
    all_starts = []
    for start_circle in circles:
        x0 = start_circle['Xo']
        y0 = start_circle['Yo']
        r0 = y0 - start_circle['Depth']
        if diagnostic:
            print(f"\n[START CIRCLE] (x0={x0:.2f}, y0={y0:.2f}, r0={r0:.2f})")
        fs_cache, best_point = evaluate_grid(x0, y0, r0)
        all_starts.append((best_point, fs_cache))

    # === Step 2: Pick the best starting location ===
    all_starts = sorted(all_starts, key=lambda t: t[0]['FS'])
    best_start, initial_cache = all_starts[0]

    x0 = best_start['Xo']
    y0 = best_start['Yo']
    r0 = y0 - best_start['Depth']

    fs_cache = initial_cache
    best_fs = best_start['FS']
    best_circle = best_start
    converged = False
    prev_depths = {(x, y): data['Depth'] for (x, y), data in initial_cache.items()}

    search_path = [{"x": x0, "y": y0, "FS": best_fs}]

    # === Step 3: Refine search iteratively ===
    for iteration in range(max_iter):
        if diagnostic:
            print(f"\n[ITERATION {iteration+1}] Center=({x0:.2f}, {y0:.2f}), Radius={r0:.2f}")

        grid_cache, best_grid_point = evaluate_grid(x0, y0, r0, prev_depths)

        for key, value in grid_cache.items():
            fs_cache[key] = value  # update cache with new points

        if best_grid_point['FS'] < best_fs:
            best_fs = best_grid_point['FS']
            best_circle = best_grid_point
            prev_depths = {(p['Xo'], p['Yo']): p['Depth'] for p in grid_cache.values()}
            search_path.append({"x": best_circle['Xo'], "y": best_circle['Yo'], "FS": best_fs})

        sorted_keys = sorted(grid_cache.keys(), key=lambda k: grid_cache[k]['FS'])
        center_index = 4  # center point of 3x3 grid

        min_index = sorted_keys.index((x0, y0))
        if min_index == center_index:
            r0 *= shrink_factor
        else:
            x0 = best_grid_point['Xo']
            y0 = best_grid_point['Yo']
            r0 = y0 - best_grid_point['Depth']

        if r0 * 0.25 < tol:
            converged = True
            break

    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    return sorted_fs_cache, converged, search_path