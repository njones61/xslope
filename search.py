import numpy as np

def circular_search(data, solver, tol=1e-2, max_iter=50, shrink_factor=0.5, fs_fail=9999, depth_tol_frac=0.01):
    """
    Global 9-point circular search using multiple starting circles.

    Returns:
        list of dict: fs_cache sorted by ascending FS
        bool: convergence status
        list of dict: search path points (x, y, FS)
    """
    from slice import generate_slices, build_ground_surface

    ground_surface = build_ground_surface(data['profile_lines'])
    circles = data['circles']
    max_depth = data['max_depth']

    # === Helper Functions ===
    def optimize_depth(x, y, r0, depth_guess):
        depth_step = r0 * 0.25
        best_depth = max(depth_guess, max_depth)
        best_fs = fs_fail
        best_result = None
        depth_tol = r0 * depth_tol_frac

        while depth_step > depth_tol:
            depths = [
                max(best_depth - depth_step, max_depth),
                best_depth,
                best_depth + depth_step
            ]
            fs_results = []
            for d in depths:
                test_circle = {'Xo': x, 'Yo': y, 'Depth': d}
                success, result = generate_slices(data, ground_surface=ground_surface, circle=test_circle)
                if not success:
                    FS = fs_fail
                    df_slices = None
                    failure_surface = None
                else:
                    df_slices, failure_surface = result
                    try:
                        solver_success, solver_result = solver(df_slices)
                        if not solver_success:
                            FS = fs_fail
                        else:
                            FS = solver_result['FS']
                    except:
                        FS = fs_fail

                fs_results.append((FS, d, df_slices, failure_surface))

            fs_results.sort(key=lambda t: t[0])
            best_fs, best_depth, best_df, best_surface = fs_results[0]
            best_result = (best_depth, best_fs, best_df, best_surface)
            depth_step *= shrink_factor

        return best_result

    def evaluate_grid(x0, y0, r0, prev_depths={}):
        Xs = [x0 - r0 * 0.25, x0, x0 + r0 * 0.25]
        Ys = [y0 - r0 * 0.25, y0, y0 + r0 * 0.25]
        points = [(x, y) for y in Ys for x in Xs]

        fs_cache = {}
        for x, y in points:
            depth_guess = prev_depths.get((x, y), y - r0)
            depth, FS, df_slices, failure_surface = optimize_depth(x, y, r0, depth_guess)
            fs_cache[(x, y)] = {
                "Xo": x,
                "Yo": y,
                "Depth": depth,
                "FS": FS,
                "slices": df_slices,
                "failure_surface": failure_surface
            }

        sorted_fs = sorted(fs_cache.items(), key=lambda item: item[1]['FS'])
        best_point = sorted_fs[0][1]

        return fs_cache, best_point

    # === Step 1: Evaluate grids for all starting circles ===
    all_starts = []
    for start_circle in circles:
        x0 = start_circle['Xo']
        y0 = start_circle['Yo']
        r0 = y0 - start_circle['Depth']
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