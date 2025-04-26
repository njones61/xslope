import numpy as np

def circular_search(data, solver, circle, tol=1e-2, max_iter=50, shrink_factor=0.5, fs_fail=9999, depth_tol_frac=0.01):
    """
    Adaptive 9-point circular search with dynamic slicing, depth optimization, and convergence tracking.

    Parameters:
        data (dict): Global data structure from load_globals()
        solver (callable): Solver function (e.g., oms, bishop, spencer, janbu_corrected)
        circle (dict): Starting circle {'Xo', 'Yo', 'Depth'}
        tol (float): Convergence tolerance
        max_iter (int): Maximum number of iterations
        shrink_factor (float): Grid shrink factor when center is minimum
        fs_fail (float): FS assigned if solver fails
        depth_tol_frac (float): Depth search convergence fraction

    Returns:
        list of dict: fs_cache sorted by ascending FS
        bool: convergence status
    """
    from slice import generate_slices, build_ground_surface

    x0 = circle['Xo']
    y0 = circle['Yo']
    r0 = y0 - circle['Depth']

    max_depth = data["max_depth"]

    ground_surface = build_ground_surface(data['profile_lines'])

    grid_size = r0 * 0.25
    best_fs = np.inf
    best_circle = None
    fs_cache = {}
    converged = False
    depth_tol = r0 * depth_tol_frac
    prev_depths = {}
    search_path = []

    def optimize_depth(x, y, r0, depth_guess):
        depth_step = r0 * 0.25
        best_depth = max(depth_guess, max_depth)
        best_fs = fs_fail
        best_result = None

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

    for iteration in range(max_iter):
        Xs = [x0 - grid_size, x0, x0 + grid_size]
        Ys = [y0 - grid_size, y0, y0 + grid_size]
        points = [(x, y) for y in Ys for x in Xs]

        results = {}
        for x, y in points:
            key = (x, y)
            if key in fs_cache:
                results[key] = fs_cache[key]
            else:
                depth_guess = prev_depths.get((x, y), y - r0)
                depth, FS, df_slices, failure_surface = optimize_depth(x, y, r0, depth_guess)
                result = {
                    "Xo": x,
                    "Yo": y,
                    "Depth": depth,
                    "FS": FS,
                    "slices": df_slices,
                    "failure_surface": failure_surface
                }
                results[key] = result
                fs_cache[key] = result

        sorted_results = sorted(results.items(), key=lambda item: item[1]['FS'])
        (min_x, min_y), min_data = sorted_results[0]
        min_index = points.index((min_x, min_y))

        if min_data['FS'] < best_fs:
            best_fs = min_data['FS']
            best_circle = min_data
            prev_depths = {(x, y): data['Depth'] for (x, y), data in results.items()}
            search_path.append({"x": min_x, "y": min_y, "FS": best_fs})

        if min_index == 4:
            grid_size *= shrink_factor
        else:
            x0, y0 = min_x, min_y
            r0 = y0 - min_data['Depth']

        if grid_size < tol:
            converged = True
            break

    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    return sorted_fs_cache, converged