import numpy as np

def circular_search(df, solver, circle, max_depth, tol=1e-2, max_iter=50, shrink_factor=0.5, fs_fail=9999, depth_tol_frac=0.01):
    """
    Adaptive 9-point circular search with depth optimization and max depth limit.

    Parameters:
        df (pd.DataFrame): Slice geometry
        solver (callable): Function that takes a circle dict {'Xo', 'Yo', 'Depth'} and returns FS
        circle (dict): Starting circle with keys 'Xo', 'Yo', 'Depth'
        max_depth (float): Minimum allowable circle base depth (Yo - R)
        tol (float): Grid tolerance for convergence
        max_iter (int): Maximum iterations
        shrink_factor (float): Shrinks grid size when center point is min
        fs_fail (float): FS value to assign on solver failure
        depth_tol_frac (float): Depth optimization convergence fraction

    Returns:
        dict: Best circle found
        list of dict: Grid points with FS values
        bool: Converged
        list of dict: Search path with min FS progression
    """
    x0 = circle['Xo']
    y0 = circle['Yo']
    r0 = y0 - circle['Depth']

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

        while depth_step > depth_tol:
            depths = [
                max(best_depth - depth_step, max_depth),
                best_depth,
                best_depth + depth_step
            ]
            fs_results = []
            for d in depths:
                test_circle = {'Xo': x, 'Yo': y, 'Depth': d}
                try:
                    FS, *_ = solver(df, test_circle)
                except:
                    FS = fs_fail
                fs_results.append((FS, d))

            fs_results.sort(key=lambda t: t[0])
            best_fs, best_depth = fs_results[0]
            depth_step *= shrink_factor

        return best_depth, best_fs

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
                depth, FS = optimize_depth(x, y, r0, depth_guess)
                result = {"Xo": x, "Yo": y, "Depth": depth, "FS": FS}
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

    grid_points = [{"x": x, "y": y, "FS": data["FS"]} for (x, y), data in fs_cache.items()]
    return best_circle, grid_points, converged, search_path