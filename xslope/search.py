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
from .advanced import rapid_drawdown
from .slice import generate_slices, get_y_from_intersection

def circular_search(slope_data, method_name, rapid=False, tol=1e-2, max_iter=50, shrink_factor=0.5,
                    fs_fail=9999, depth_tol_frac=0.03, diagnostic=False):
    """
    Global 9-point circular search with adaptive grid refinement.

    Returns:
        list of dict: sorted fs_cache by FS
        bool: convergence flag
        list of dict: search path
    """

    solver = getattr(solve, method_name)

    start_time = time.time()  # Start timing

    ground_surface = slope_data['ground_surface']
    ground_surface = LineString([(x, y) for x, y in ground_surface.coords])
    y_max = max(y for _, y in ground_surface.coords)
    y_min = slope_data['max_depth']
    delta_y = y_max - y_min
    tol = delta_y * depth_tol_frac

    circles = slope_data['circles']
    max_depth = slope_data['max_depth']

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
                success, result = generate_slices(slope_data, circle=test_circle)
                if not success:
                    FS = fs_fail
                    df_slices = None
                    failure_surface = None
                    solver_result = None
                else:
                    df_slices, failure_surface = result
                    if rapid:
                        solver_success, solver_result = rapid_drawdown(df_slices, method_name)
                    else:
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

    def evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=False, fs_cache=None):
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
        fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=diagnostic)
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
        fs_cache, best_point = evaluate_grid(x0, y0, grid_size, depth_guess, slope_data, diagnostic=diagnostic, fs_cache=fs_cache)

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
            print(f"[✅ converged] Iter={iteration+1}, FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f}, depth={depth_guess:.2f}), elapsed time={elapsed:.2f} seconds")
            break

    if not converged and diagnostic:
        print(f"\n[❌ max iterations reached] FS={best_fs:.4f} at (x={x0:.2f}, y={y0:.2f})")

    sorted_fs_cache = sorted(fs_cache.values(), key=lambda d: d['FS'])
    return sorted_fs_cache, converged, search_path

def noncircular_search(slope_data, method_name, rapid=False, diagnostic=True, movement_distance=4.0, shrink_factor=0.8, fs_tol=0.001, max_iter=100, move_tol=0.1):
    """
    Non-circular search using the specified solver.
    
    Parameters:
    -----------
    data : dict
        Input data dictionary containing all necessary parameters
    method_name : str
        The method name to use (e.g., 'lowe_karafiath', 'spencer')
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
        
    Returns:
    --------
    tuple : (fs_cache, converged, search_path)
        fs_cache : dict of all evaluated surfaces and their FS values
        converged : bool indicating if search converged
        search_path : list of surfaces evaluated during search
    """
    # Get the solver function from solve module
    solver = getattr(solve, method_name)
    def move_point(points, i, dx, dy, movement_type, ground_surface, max_depth):
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
            # For middle points, ensure they stay below ground surface but above max_depth
            if new_y > ground_surface.interpolate(ground_surface.project(Point(new_x, new_y))).y:
                return False
            if new_y < max_depth:
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
        
        # Generate slices and compute FS
        success, result = generate_slices(slope_data, non_circ=non_circ)
        if not success:
            return float('inf'), None, None, None, fs_cache
            
        df_slices, failure_surface = result
        if rapid:
            solver_success, solver_result = rapid_drawdown(df_slices, method_name)
        else:
            solver_success, solver_result = solver(df_slices)
        FS = solver_result['FS'] if solver_success else float('inf')
        
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
                if move_point(test_points, i, dx, dy, movements[i], ground_surface, slope_data['max_depth']):
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
