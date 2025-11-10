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

from math import sin, cos, tan, radians, atan, atan2, degrees, sqrt

import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint, GeometryCollection

from .mesh import find_element_containing_point, interpolate_at_point

def get_circular_y_coordinates(x_coords, Xo, Yo, R):
    """
    Calculate y-coordinates on a circular failure surface for given x-coordinates.
    
    Parameters:
        x_coords (array-like): X-coordinates to evaluate
        Xo, Yo (float): Center coordinates of the circle
        R (float): Radius of the circle
        
    Returns:
        numpy.ndarray: Y-coordinates on the circle (bottom half)
    """
    x_coords = np.asarray(x_coords)
    # Calculate y-coordinates for the bottom half of the circle
    # y = Yo - sqrt(R^2 - (x - Xo)^2)
    dx_squared = (x_coords - Xo) ** 2
    # Handle cases where x is outside the circle
    valid_mask = dx_squared <= R ** 2
    y_coords = np.full_like(x_coords, np.nan)
    y_coords[valid_mask] = Yo - np.sqrt(R ** 2 - dx_squared[valid_mask])
    return y_coords


def get_circular_intersection_points(ground_surface, Xo, Yo, R, x_min, x_max):
    """
    Find intersection points between a circular failure surface and ground surface.
    
    Parameters:
        ground_surface (LineString): Ground surface geometry
        Xo, Yo, R (float): Circle parameters
        x_min, x_max (float): X-range to search
        
    Returns:
        tuple: (x_min, x_max, y_left, y_right, success)
    """
    # Create a dense set of points on the circle for intersection testing
    x_test = np.linspace(x_min, x_max, 1000)
    y_circle = get_circular_y_coordinates(x_test, Xo, Yo, R)
    
    # Create circle line for intersection
    valid_mask = ~np.isnan(y_circle)
    if not np.any(valid_mask):
        return None, None, None, None, False
    
    circle_coords = list(zip(x_test[valid_mask], y_circle[valid_mask]))
    circle_line = LineString(circle_coords)
    
    # Find intersections
    intersections = circle_line.intersection(ground_surface)
    
    if isinstance(intersections, Point):
        points = [intersections]
    elif isinstance(intersections, MultiPoint):
        points = list(intersections.geoms)
    elif isinstance(intersections, GeometryCollection):
        points = [g for g in intersections.geoms if isinstance(g, Point)]
    else:
        points = []
    
    if len(points) < 2:
        return None, None, None, None, False
    
    # Sort by x and take the two endpoints
    points = sorted(points, key=lambda p: p.x)
    x_min, x_max = points[0].x, points[-1].x
    y_left, y_right = points[0].y, points[-1].y
    
    return x_min, x_max, y_left, y_right, True


def get_ground_surface_y_coordinates(x_coords, ground_surface):
    """
    Get y-coordinates on the ground surface for given x-coordinates using interpolation.
    
    Parameters:
        x_coords (array-like): X-coordinates to evaluate
        ground_surface (LineString): Ground surface geometry
        
    Returns:
        numpy.ndarray: Y-coordinates on the ground surface
    """
    x_coords = np.asarray(x_coords)
    ground_coords = np.array(ground_surface.coords)
    ground_x = ground_coords[:, 0]
    ground_y = ground_coords[:, 1]
    
    # Sort by x to ensure proper interpolation
    sort_idx = np.argsort(ground_x)
    ground_x = ground_x[sort_idx]
    ground_y = ground_y[sort_idx]
    
    # Interpolate y-coordinates
    y_coords = np.interp(x_coords, ground_x, ground_y, left=np.nan, right=np.nan)
    return y_coords


def get_profile_layer_y_coordinates(x_coords, profile_lines):
    """
    Get y-coordinates for each profile layer at given x-coordinates.
    
    Parameters:
        x_coords (array-like): X-coordinates to evaluate
        profile_lines (list): List of profile layer lines
        
    Returns:
        list: List of arrays, each containing y-coordinates for a profile layer
    """
    x_coords = np.asarray(x_coords)
    layer_y_coords = []
    
    for line in profile_lines:
        line_coords = np.array(line)
        line_x = line_coords[:, 0]
        line_y = line_coords[:, 1]
        
        # Sort by x to ensure proper interpolation
        sort_idx = np.argsort(line_x)
        line_x = line_x[sort_idx]
        line_y = line_y[sort_idx]
        
        # Interpolate y-coordinates
        y_coords = np.interp(x_coords, line_x, line_y, left=np.nan, right=np.nan)
        layer_y_coords.append(y_coords)
    
    return layer_y_coords


def get_piezometric_y_coordinates(x_coords, piezo_line):
    """
    Get y-coordinates on the piezometric surface for given x-coordinates.
    
    Parameters:
        x_coords (array-like): X-coordinates to evaluate
        piezo_line (list): Piezometric line coordinates
        
    Returns:
        numpy.ndarray: Y-coordinates on the piezometric surface
    """
    if not piezo_line:
        return np.full_like(x_coords, np.nan)
    
    x_coords = np.asarray(x_coords)
    piezo_coords = np.array(piezo_line)
    piezo_x = piezo_coords[:, 0]
    piezo_y = piezo_coords[:, 1]
    
    # Sort by x to ensure proper interpolation
    sort_idx = np.argsort(piezo_x)
    piezo_x = piezo_x[sort_idx]
    piezo_y = piezo_y[sort_idx]
    
    # Interpolate y-coordinates
    y_coords = np.interp(x_coords, piezo_x, piezo_y, left=np.nan, right=np.nan)
    return y_coords


def circle_polyline_intersections(Xo, Yo, R, polyline):
    """
    Find intersection points between the bottom half of a circle and a polyline (LineString).
    Returns a list of shapely Point objects.
    """
    intersections = []
    coords = list(polyline.coords)
    for i in range(len(coords) - 1):
        x1, y1 = coords[i]
        x2, y2 = coords[i+1]
        dx = x2 - x1
        dy = y2 - y1

        # Quadratic coefficients for t
        a = dx**2 + dy**2
        b = 2 * (dx * (x1 - Xo) + dy * (y1 - Yo))
        c = (x1 - Xo)**2 + (y1 - Yo)**2 - R**2

        discriminant = b**2 - 4*a*c
        if discriminant < 0:
            continue  # No intersection

        sqrt_disc = np.sqrt(discriminant)
        for sign in [-1, 1]:
            t = (-b + sign * sqrt_disc) / (2 * a)
            if 0 <= t <= 1:
                xi = x1 + t * dx
                yi = y1 + t * dy
                if yi < Yo:  # Only keep points below the center (bottom half)
                    intersections.append(Point(xi, yi))
    return intersections


def get_sorted_intersections(failure_surface, ground_surface, circle_params=None):
    """
    Find and sort the intersection points between the failure and ground surfaces,
    pruning extras if the circle exits and re-enters the ground beyond the toe.
    If circle_params is provided, use analytic intersection.
    Returns:
        success (bool), msg (str), points (list of shapely Point)
    """
    if circle_params is not None:
        Xo, Yo, R = circle_params['Xo'], circle_params['Yo'], circle_params['R']
        points = circle_polyline_intersections(Xo, Yo, R, ground_surface)
    else:
        intersections = failure_surface.intersection(ground_surface)
        if isinstance(intersections, MultiPoint):
            points = list(intersections.geoms)
        elif isinstance(intersections, Point):
            points = [intersections]
        elif isinstance(intersections, GeometryCollection):
            points = [g for g in intersections.geoms if isinstance(g, Point)]
        else:
            points = []

    # need at least two
    if len(points) < 2:
        return False, f"Expected at least 2 intersection points, but got {len(points)}.", None

    # sort by x
    points = sorted(points, key=lambda p: p.x)

    # if exactly two, we're done
    if len(points) == 2:
        left, right = points[0], points[1]
        tol = 1e-6
        if abs(left.y - right.y) < tol:
            return False, "Rejected: left-most and right-most intersection points have the same y-value (flat arc).", None
        return True, "", points

    # more than two: decide facing
    y_first, y_last = points[0].y, points[-1].y
    if y_first > y_last:
        # right-facing: keep first two
        pruned = points[:2]
    else:
        # left-facing: keep last two
        pruned = points[-2:]

    # sort those two again by x (just in case)
    pruned = sorted(pruned, key=lambda p: p.x)
    left, right = pruned[0], pruned[1]
    tol = 1e-6
    if abs(left.y - right.y) < tol:
        return False, "Rejected: left-most and right-most intersection points have the same y-value (flat arc).", None
    return True, "", pruned


def adjust_ground_for_tcrack(ground_surface, x_center, tcrack_depth, right_facing):
    # helper function to adjust the ground surface for tension crack
    if tcrack_depth <= 0:
        return ground_surface

    new_coords = []
    for x, y in ground_surface.coords:
        if right_facing and x < x_center:
            new_coords.append((x, y - tcrack_depth))
        elif not right_facing and x > x_center:
            new_coords.append((x, y - tcrack_depth))
        else:
            new_coords.append((x, y))
    return LineString(new_coords)


def generate_failure_surface(ground_surface, circular, circle=None, non_circ=None, tcrack_depth=0):
    """
    Generates a failure surface based on either a circular or non-circular definition.

    Parameters:
        ground_surface (LineString): The ground surface geometry.
        circular (bool): Whether to use circular failure surface.
        circle (dict, optional): Dictionary with keys 'Xo', 'Yo', 'Depth', and 'R'.
        non_circ (list, optional): List of dicts with keys 'X', 'Y', and 'Movement'.
        tcrack_depth (float, optional): Tension crack depth.

    Returns:
        tuple: (success, result)
            - If success is True:
                result = (x_min, x_max, y_left, y_right, clipped_surface)
            - If success is False:
                result = error message string
    """
    # --- Step 1: Build failure surface ---
    if circular and circle:
        Xo, Yo, depth, R = circle['Xo'], circle['Yo'], circle['Depth'], circle['R']
        theta_range = np.linspace(np.pi, 2 * np.pi, 100)
        arc = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
        failure_coords = arc
        failure_surface = LineString(arc)
    elif non_circ:
        failure_coords = [(pt['X'], pt['Y']) for pt in non_circ]
        failure_surface = LineString(failure_coords)
    else:
        return False, "Either a circular or non-circular failure surface must be provided."

    # --- Step 2: Intersect with original ground surface to determine slope facing ---
    if circular and circle:
        success, msg, points = get_sorted_intersections(failure_surface, ground_surface, circle_params=circle)
    else:
        success, msg, points = get_sorted_intersections(failure_surface, ground_surface)
    if not success:
        return False, msg

    x_min, x_max = points[0].x, points[1].x
    y_left, y_right = points[0].y, points[1].y
    right_facing = y_left > y_right
    x_center = 0.5 * (x_min + x_max)

    # --- Step 3: If tension crack exists, adjust surface and re-intersect ---
    if tcrack_depth > 0:
        modified_surface = adjust_ground_for_tcrack(ground_surface, x_center, tcrack_depth, right_facing)
        if circular and circle:
            success, msg, points = get_sorted_intersections(failure_surface, modified_surface, circle_params=circle)
        else:
            success, msg, points = get_sorted_intersections(failure_surface, modified_surface)
        if not success:
            return False, msg
        x_min, x_max = points[0].x, points[1].x
        y_left, y_right = points[0].y, points[1].y

    # --- Step 4: Clip the failure surface between intersection x-range ---
    # Filter coordinates within the x-range
    filtered_coords = [pt for pt in failure_coords if x_min <= pt[0] <= x_max]
    
    # Add the exact intersection points if they're not already in the filtered list
    left_intersection = (x_min, y_left)
    right_intersection = (x_max, y_right)
    
    # Check if intersection points are already in the filtered list (with tolerance)
    tol = 1e-6
    has_left = any(abs(pt[0] - x_min) < tol and abs(pt[1] - y_left) < tol for pt in filtered_coords)
    has_right = any(abs(pt[0] - x_max) < tol and abs(pt[1] - y_right) < tol for pt in filtered_coords)
    
    if not has_left:
        filtered_coords.insert(0, left_intersection)
    if not has_right:
        filtered_coords.append(right_intersection)
    
    # Sort by x-coordinate to ensure proper ordering
    filtered_coords.sort(key=lambda pt: pt[0])
    
    clipped_surface = LineString(filtered_coords)

    return True, (x_min, x_max, y_left, y_right, clipped_surface)


def get_y_from_intersection(geom):
    """
    Extracts the maximum Y-coordinate from a geometric intersection result.

    This function handles different geometric types resulting from intersections,
    including Point, MultiPoint, LineString, and GeometryCollection. If the input
    geometry is not one of these or is empty, the function returns None.

    Parameters:
        geom (shapely.geometry.base.BaseGeometry): The geometry object from which
            to extract the Y-coordinate(s).

    Returns:
        float or None: The maximum Y-coordinate found in the geometry, or None if not found.
    """
    if isinstance(geom, Point):
        return geom.y
    elif isinstance(geom, MultiPoint):
        return max(pt.y for pt in geom.geoms)
    elif isinstance(geom, LineString):
        return max(y for _, y in geom.coords) if geom.coords else None
    elif isinstance(geom, GeometryCollection):
        pts = [g for g in geom.geoms if isinstance(g, Point)]
        return max(pt.y for pt in pts) if pts else None
    return None


def calc_dload_resultant(x_l, y_lt, x_r, y_rt, qL, qR, dl):
    """
    Compute:
      - D    : total resultant force from a trapezoidal load varying
               linearly from intensity qL at (x_l,y_lt) to qR at (x_r,y_rt).
      - d_x  : x‐coordinate of the resultant's centroid on the top edge
      - d_y  : y‐coordinate of the resultant's centroid on the top edge

    Parameters
    ----------
    x_l, y_lt : float
        Coordinates of the left‐end of the top edge.
    x_r, y_rt : float
        Coordinates of the right‐end of the top edge.
    qL : float
        Load intensity (force per unit length) at (x_l, y_lt).
    qR : float
        Load intensity (force per unit length) at (x_r, y_rt).
    dl : float
        Actual length along the inclined surface.

    Returns
    -------
    D    : float
           Total resultant (area of trapezoid) = ½ (qL + qR) * dl
    d_x  : float
           Global x‐coordinate of the centroid of that trapezoid
    d_y  : float
           Global y‐coordinate of the centroid (lies on the line segment
           between (x_l,y_lt) and (x_r,y_rt))

    Notes
    -----
    1.  If x_r == x_l (zero‐width slice), this will return D=0 and place
        the "centroid" at (x_l, y_lt).
    2.  For a nonzero‐width trapezoid, the horizontal centroid‐offset from
        x_l is:
             x_offset = (x_r – x_l) * ( qL + 2 qR ) / [3 (qL + qR) ]
        provided (qL + qR) ≠ 0.  If qL + qR ≈ 0, it simply places the
        centroid at the midpoint in x.
    3.  The vertical coordinate d_y is found by linear‐interpolation:
          t = x_offset / (x_r – x_l)
          d_y = y_lt + t ·(y_rt – y_lt)

    """
    dx = x_r - x_l

    # 1) Total resultant force (area under trapezoid) using actual length
    D = 0.5 * (qL + qR) * dl

    # 2) Horizontal centroid offset from left end
    sum_q = qL + qR
    if abs(sum_q) < 1e-12:
        # nearly zero trapezoid => centroid at geometric midpoint
        x_offset = dx * 0.5
    else:
        x_offset = dx * (qL + 2.0 * qR) / (3.0 * sum_q)

    # 3) Global x‐coordinate of centroid
    d_x = x_l + x_offset

    # 4) Corresponding y‐coordinate by linear interpolation along top edge
    t = x_offset / dx
    d_y = y_lt + t * (y_rt - y_lt)

    return D, d_x, d_y


def generate_slices(slope_data, circle=None, non_circ=None, num_slices=40, debug=True):

    """
    Generates vertical slices between the ground surface and a failure surface for slope stability analysis.

    This function supports both circular and non-circular failure surfaces and computes
    geometric and mechanical properties for each slice, including weight, base geometry,
    water pressures, distributed loads, and reinforcement effects.

    Parameters:
        data (dict): Dictionary containing all input data
        circle (dict, optional): Dictionary with keys 'Xo', 'Yo', 'Depth', and 'R' defining the circular failure surface.
        non_circ (list, optional): List of dicts defining a non-circular failure surface with keys 'X', 'Y', and 'Movement'.
        num_slices (int, optional): Desired number of slices to generate (default is 40).
        debug (bool, optional): Whether to print debug information (default is True).

    Returns:
        tuple:
            - pd.DataFrame: Slice table where each row includes geometry, strength, and external force values.
            - shapely.geometry.LineString: The clipped failure surface between the ground surface intersections.

    Notes:
        - Supports Method A interpretation of reinforcement: T reduces driving forces.
        - Handles pore pressure and distributed loads using linear interpolation at slice centers.
        - Automatically includes all geometry breakpoints in slice generation.
        - Must specify exactly one of 'circle' or 'non_circ'.
    """

    # Unpack data
    profile_lines = slope_data["profile_lines"]
    ground_surface = slope_data["ground_surface"]
    materials = slope_data["materials"]
    piezo_line = slope_data["piezo_line"]
    piezo_line2 = slope_data.get("piezo_line2", [])  # Second piezometric line
    gamma_w = slope_data["gamma_water"]
    tcrack_depth = slope_data["tcrack_depth"]
    tcrack_water = slope_data["tcrack_water"]
    k_seismic = slope_data['k_seismic']
    dloads = slope_data["dloads"]
    dloads2 = slope_data.get("dloads2", [])
    max_depth = slope_data["max_depth"]

    # Determine failure surface type
    if circle is not None:
        circular = True
        Xo, Yo, depth, R = circle['Xo'], circle['Yo'], circle['Depth'], circle['R']
    else:
        circular = False

    # Prepare reinforcement lines data
    reinf_lines_data = []
    if slope_data.get("reinforce_lines"):
        for line in slope_data["reinforce_lines"]:
            xs = [pt["X"] for pt in line]
            ts = [pt["T"] for pt in line]
            geom = LineString([(pt["X"], pt["Y"]) for pt in line])
            reinf_lines_data.append({"xs": xs, "ts": ts, "geom": geom})

    ground_surface = LineString([(x, y) for x, y in ground_surface.coords])

    # Generate failure surface
    success, result = generate_failure_surface(ground_surface, circular, circle=circle, non_circ=non_circ, tcrack_depth=tcrack_depth)
    if success:
        x_min, x_max, y_left, y_right, clipped_surface = result
    else:
        return False, "Failed to generate surface:" + result

    # Determine if the failure surface is right-facing
    right_facing = y_left > y_right

    # === BEGIN : Find set of points that should correspond to slice boundaries. ===

    # Find set of points that are on the profile lines if the points are above the failure surface.
    fixed_xs = set()
    
    # Vectorized approach for profile line points
    for line in profile_lines:
        line_coords = np.array(line)
        x_coords = line_coords[:, 0]
        y_coords = line_coords[:, 1]
        
        # Filter points within x-range
        mask = (x_coords >= x_min) & (x_coords <= x_max)
        x_filtered = x_coords[mask]
        y_filtered = y_coords[mask]
        
        if len(x_filtered) == 0:
            continue
            
        # Check if points are above the failure surface
        if circular:
            # Use parametric equation for circular failure surface
            failure_y = get_circular_y_coordinates(x_filtered, Xo, Yo, R)
            above_mask = y_filtered > failure_y
        else:
            # For non-circular, use geometric intersection (slower but necessary)
            above_mask = np.zeros(len(x_filtered), dtype=bool)
            for i, (x, y) in enumerate(zip(x_filtered, y_filtered)):
                vertical_line = LineString([(x, -1e6), (x, 1e6)])
                failure_y = get_y_from_intersection(clipped_surface.intersection(vertical_line))
                if failure_y is not None and y > failure_y:
                    above_mask[i] = True
        
        # Add points that are above the failure surface
        fixed_xs.update(x_filtered[above_mask])
    
    fixed_xs.update([x_min, x_max])

    # Add transition points from dloads
    if dloads:
        fixed_xs.update(
            pt['X'] for line in dloads for pt in line
            if x_min <= pt['X'] <= x_max
        )

    # Add transition points from dloads2
    if dloads2:
        fixed_xs.update(
            pt['X'] for line in dloads2 for pt in line
            if x_min <= pt['X'] <= x_max
        )

    # Add transition points from non_circ
    if non_circ:
        fixed_xs.update(
            pt['X'] for pt in non_circ
            if x_min <= pt['X'] <= x_max
        )

    # Find intersections with profile lines and failure surface
    if circular:
        # For circular failure surfaces, we can use a more efficient approach
        # by creating a dense circle representation and finding intersections
        for line in profile_lines:
            line_geom = LineString(line)
            # Create a dense circle representation for intersection
            theta_range = np.linspace(np.pi, 2 * np.pi, 200)
            circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
            circle_line = LineString(circle_coords)
            
            intersection = line_geom.intersection(circle_line)
            if not intersection.is_empty:
                if hasattr(intersection, 'x'):
                    fixed_xs.add(intersection.x)
                elif hasattr(intersection, 'geoms'):
                    for geom in intersection.geoms:
                        if hasattr(geom, 'x'):
                            fixed_xs.add(geom.x)
    else:
        # For non-circular failure surfaces, use the original approach
        for i in range(len(profile_lines)):
            intersection = LineString(profile_lines[i]).intersection(clipped_surface)
            if not intersection.is_empty:
                if hasattr(intersection, 'x'):
                    fixed_xs.add(intersection.x)
                elif hasattr(intersection, 'geoms'):
                    for geom in intersection.geoms:
                        if hasattr(geom, 'x'):
                            fixed_xs.add(geom.x)

    # Find intersections with piezometric lines
    if piezo_line:
        piezo_geom1 = LineString(piezo_line)
        if circular:
            # Use dense circle representation for intersection
            theta_range = np.linspace(np.pi, 2 * np.pi, 200)
            circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
            circle_line = LineString(circle_coords)
            intersection1 = piezo_geom1.intersection(circle_line)
        else:
            intersection1 = piezo_geom1.intersection(clipped_surface)
            
        if not intersection1.is_empty:
            if hasattr(intersection1, 'x'):
                # Single point intersection
                if x_min <= intersection1.x <= x_max:
                    fixed_xs.add(intersection1.x)
            elif hasattr(intersection1, 'geoms'):
                # Multiple points or line intersection
                for geom in intersection1.geoms:
                    if hasattr(geom, 'x') and x_min <= geom.x <= x_max:
                        fixed_xs.add(geom.x)

    if piezo_line2:
        piezo_geom2 = LineString(piezo_line2)
        if circular:
            theta_range = np.linspace(np.pi, 2 * np.pi, 200)
            circle_coords = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
            circle_line = LineString(circle_coords)
            intersection2 = piezo_geom2.intersection(circle_line)
        else:
            intersection2 = piezo_geom2.intersection(clipped_surface)
            
        if not intersection2.is_empty:
            if hasattr(intersection2, 'x'):
                # Single point intersection
                if x_min <= intersection2.x <= x_max:
                    fixed_xs.add(intersection2.x)
            elif hasattr(intersection2, 'geoms'):
                # Multiple points or line intersection
                for geom in intersection2.geoms:
                    if hasattr(geom, 'x') and x_min <= geom.x <= x_max:
                        fixed_xs.add(geom.x)

    # Remove duplicate points that are very close to each other
    tolerance = 1e-6
    cleaned_xs = []
    for x in sorted(fixed_xs):
        if not cleaned_xs or abs(x - cleaned_xs[-1]) > tolerance:
            cleaned_xs.append(x)
    
    fixed_xs = cleaned_xs

    # Generate slice boundaries
    segment_lengths = [fixed_xs[i + 1] - fixed_xs[i] for i in range(len(fixed_xs) - 1)]
    total_length = sum(segment_lengths)
    all_xs = [fixed_xs[0]]
    for i in range(len(fixed_xs) - 1):
        x_start = fixed_xs[i]
        x_end = fixed_xs[i + 1]
        segment_length = x_end - x_start
        n_subdiv = max(1, int(round((segment_length / total_length) * num_slices)))
        xs = np.linspace(x_start, x_end, n_subdiv + 1).tolist()
        all_xs.extend(xs[1:])

    # === END : Find set of points that should correspond to slice boundaries. ===

    # Remove thin slices (width < 1e-2), including at the ends
    min_width = 1e-2
    cleaned_xs = [all_xs[0]]
    for x in all_xs[1:]:
        if abs(x - cleaned_xs[-1]) >= min_width:
            cleaned_xs.append(x)
    
    # If the last slice is thin, merge it with the previous one
    if len(cleaned_xs) > 2 and abs(cleaned_xs[-1] - cleaned_xs[-2]) < min_width:
        cleaned_xs.pop(-2)
    
    all_xs = cleaned_xs

    # Pre-compute all y-coordinates for efficiency
    slice_x_coords = np.array(all_xs)
    slice_centers = (slice_x_coords[:-1] + slice_x_coords[1:]) / 2
    
    # Get failure surface y-coordinates
    if circular:
        # Use parametric equations for circular failure surface
        y_lb_all = get_circular_y_coordinates(slice_x_coords[:-1], Xo, Yo, R)
        y_rb_all = get_circular_y_coordinates(slice_x_coords[1:], Xo, Yo, R)
        y_cb_all = get_circular_y_coordinates(slice_centers, Xo, Yo, R)
    else:
        # For non-circular, we need to use geometric intersections
        y_lb_all = np.array([get_y_from_intersection(clipped_surface.intersection(LineString([(x, -1e6), (x, 1e6)]))) for x in slice_x_coords[:-1]])
        y_rb_all = np.array([get_y_from_intersection(clipped_surface.intersection(LineString([(x, -1e6), (x, 1e6)]))) for x in slice_x_coords[1:]])
        y_cb_all = np.array([get_y_from_intersection(clipped_surface.intersection(LineString([(x, -1e6), (x, 1e6)]))) for x in slice_centers])

    # Get ground surface y-coordinates
    y_lt_all = get_ground_surface_y_coordinates(slice_x_coords[:-1], ground_surface)
    y_rt_all = get_ground_surface_y_coordinates(slice_x_coords[1:], ground_surface)
    y_ct_all = get_ground_surface_y_coordinates(slice_centers, ground_surface)

    # Get profile layer y-coordinates
    profile_y_coords = get_profile_layer_y_coordinates(slice_centers, profile_lines)

    # Get piezometric y-coordinates
    piezo_y_all = get_piezometric_y_coordinates(slice_centers, piezo_line)
    piezo_y2_all = get_piezometric_y_coordinates(slice_centers, piezo_line2)

    # Interpolation functions for distributed loads
    dload_interp_funcs = []
    if dloads:
        for line in dloads:
            xs = [pt['X'] for pt in line]
            normals = [pt['Normal'] for pt in line]
            dload_interp_funcs.append(lambda x, xs=xs, normals=normals: np.interp(x, xs, normals, left=0, right=0))

    # Interpolation functions for second set of distributed loads
    dload2_interp_funcs = []
    if dloads2:
        for line in dloads2:
            xs = [pt['X'] for pt in line]
            normals = [pt['Normal'] for pt in line]
            dload2_interp_funcs.append(lambda x, xs=xs, normals=normals: np.interp(x, xs, normals, left=0, right=0))

    # Generate slices
    slices = []
    for i in range(len(all_xs) - 1):
        x_l, x_r = slice_x_coords[i], slice_x_coords[i + 1]
        x_c = slice_centers[i]
        dx = x_r - x_l

        # Get pre-computed y-coordinates
        y_lb = y_lb_all[i]
        y_rb = y_rb_all[i]
        y_cb = y_cb_all[i]
        y_lt = y_lt_all[i]
        y_rt = y_rt_all[i]
        y_ct = y_ct_all[i]

        # Skip if any coordinates are invalid
        if any(np.isnan([y_lb, y_rb, y_cb, y_lt, y_rt, y_ct])):
            continue

        # Calculate beta (slope angle of the top edge) in degrees
        beta = degrees(atan2(y_rt - y_lt, x_r - x_l))
        if right_facing:
            beta = -beta

        # Calculate dl for the top surface (for distributed loads)
        dl_top = sqrt((x_r - x_l)**2 + (y_rt - y_lt)**2)

        # Calculate layer heights using pre-computed profile coordinates
        heights = []
        soil_weight = 0
        base_material_idx = None
        sum_gam_h_y = 0  # for calculating center of gravity of slice
        sum_gam_h = 0    # ditto
        
        for mat_index, layer_y in enumerate(profile_y_coords):
            layer_top_y = layer_y[i]
            
            # Bottom: highest of all other profile lines at x, or failure surface
            layer_bot_y = y_cb  # Start with failure surface as default bottom
            for j in range(mat_index + 1, len(profile_y_coords)):
                next_y = profile_y_coords[j][i]
                if not np.isnan(next_y) and next_y > layer_bot_y:
                    # Take the highest of the lower profile lines
                    layer_bot_y = next_y

            if np.isnan(layer_top_y) or np.isnan(layer_bot_y):
                h = 0
            else:
                overlap_top = min(y_ct, layer_top_y)
                overlap_bot = max(y_cb, layer_bot_y)
                h = max(0, overlap_top - overlap_bot)
                sum_gam_h_y += h * materials[mat_index]['gamma'] * (overlap_top + overlap_bot) / 2
                sum_gam_h += h * materials[mat_index]['gamma']

            heights.append(h)
            soil_weight += h * materials[mat_index]['gamma'] * dx

            if h > 0:
                base_material_idx = mat_index

        # Center of gravity
        y_cg = (sum_gam_h_y) / sum_gam_h if sum_gam_h > 0 else None

        # Distributed load
        qC = sum(func(x_c) for func in dload_interp_funcs) if dload_interp_funcs else 0   # intensity at center
        if qC > 0: # We need to check qC to distinguish between a linear ramp up (down) and the case where the load starts or ends on one of the sides
            qL = sum(func(x_l) for func in dload_interp_funcs) if dload_interp_funcs else 0   # intensity at left‐top corner
            qR = sum(func(x_r) for func in dload_interp_funcs) if dload_interp_funcs else 0   # intensity at right‐top corner
        else:
            qL = 0
            qR = 0
        dload, d_x, d_y = calc_dload_resultant(x_l, y_lt, x_r, y_rt, qL, qR, dl_top)

        # Second distributed load
        qC2 = sum(func(x_c) for func in dload2_interp_funcs) if dload2_interp_funcs else 0   # intensity at center
        if qC2 > 0: # We need to check qC2 to distinguish between a linear ramp up (down) and the case where the load starts or ends on one of the sides
            qL2 = sum(func(x_l) for func in dload2_interp_funcs) if dload2_interp_funcs else 0   # intensity at left‐top corner
            qR2 = sum(func(x_r) for func in dload2_interp_funcs) if dload2_interp_funcs else 0   # intensity at right‐top corner
        else:
            qL2 = 0
            qR2 = 0
        dload2, d_x2, d_y2 = calc_dload_resultant(x_l, y_lt, x_r, y_rt, qL2, qR2, dl_top)

        # Seismic force
        kw = k_seismic * soil_weight

        # === BEGIN : "Tension crack water force" ===

        # By default, zero out t and its line‐of‐action:
        t_force = 0.0
        y_t_loc  = 0.0

        # Only nonzero for the appropriate end‐slice:
        if tcrack_water is not None and tcrack_water > 0:
            # Horizontal resultant of water in tension crack (triangular distribution):
            #    t = (1/2) * γ_w * (d_tc)^2
            # Here, gamma_w is the unit weight of the crack‐water (y_w),
            # and tcrack_water is the depth of water in the crack (d_tc).
            t_force = 0.5 * gamma_w * (tcrack_water ** 2)

            if right_facing:
                # Right‐facing slope → water pushes on left side of the first slice (i == 0)
                if i == 0:
                    # line of action is d_tc/3 above the bottom left corner y_lb
                    t_force = - t_force  # negative because it acts to the right on free body diagram
                    y_t_loc = y_lb + (tcrack_water / 3.0)
                else:
                    # other slices = no tension‐crack force
                    t_force = 0.0
                    y_t_loc = 0.0

            else:
                # Left‐facing slope → water pushes on right side of the last slice (i == n-1)
                if i == (len(all_xs) - 2):  # last slice index = (number_of_slices − 1)
                    # line of action is d_tc/3 above the bottom right corner y_rb
                    y_t_loc = y_rb + (tcrack_water / 3.0)
                else:
                    t_force = 0.0
                    y_t_loc = y_rb
        # === END: "Tension crack water force" ===

        # === BEGIN : "Reinforcement lines" ===

        # 1) Build this slice's base as a LineString from (x_l, y_lb) to (x_r, y_rb):
        slice_base = LineString([(x_l, y_lb), (x_r, y_rb)])

        # 2) For each reinforcement line, check a single‐point intersection:
        p_sum = 0.0
        for rl in reinf_lines_data:
            intersec = slice_base.intersection(rl["geom"])
            if intersec.is_empty:
                continue

            # Since we guarantee only one intersection point, it must be a Point:
            if isinstance(intersec, Point):
                xi = intersec.x
                # interpolated T at xi
                t_i = np.interp(xi, rl["xs"], rl["ts"], left=0.0, right=0.0)
                p_sum += t_i
            else:
                # (In the extremely unlikely case that intersection is not a Point,
                #  skip it. Our assumption is only one Point per slice-base.)
                continue

        # Now p_sum is the TOTAL T‐pull acting at this slice's base.
        # === END: "Tension crack water force" ===

        # Process piezometric line and pore pressures using pre-computed coordinates
        piezo_y = piezo_y_all[i]
        piezo_y2 = piezo_y2_all[i]
        
        hw = 0
        hw2 = 0
        u = 0
        u2 = 0
        # Determine pore pressure method from material property
        mat_u = materials[base_material_idx]['u'] if base_material_idx is not None else 'none'
        if mat_u == 'none':
            u = 0
            u2 = 0
        elif mat_u == 'piezo':
            if not np.isnan(piezo_y) and piezo_y > y_cb:
                hw = piezo_y - y_cb
            if not np.isnan(piezo_y2) and piezo_y2 > y_cb:
                hw2 = piezo_y2 - y_cb
            u = hw * gamma_w if not np.isnan(piezo_y) else 0
            u2 = hw2 * gamma_w if not np.isnan(piezo_y2) else 0
        elif mat_u == 'seep':
            # Seepage-based pore pressure calculation using mesh interpolation
            if 'seep_mesh' in data and 'seep_u' in data:
                seep_mesh = data['seep_mesh']
                seep_u = data['seep_u']
                
                # Interpolate pore pressure at the slice center base point
                point = (x_c, y_cb)
                u = interpolate_at_point(
                    seep_mesh['nodes'], 
                    seep_mesh['elements'], 
                    seep_mesh['element_types'], 
                    seep_u, 
                    point
                )
            else:
                u = 0
                
            # Check for second seepage solution (rapid drawdown)
            if 'seep_u2' in data:
                seep_mesh = data['seep_mesh']
                seep_u2 = data['seep_u2']
                
                # Interpolate pore pressure at the slice center base point
                point = (x_c, y_cb)
                u2 = interpolate_at_point(
                    seep_mesh['nodes'], 
                    seep_mesh['elements'], 
                    seep_mesh['element_types'], 
                    seep_u2, 
                    point
                )
            else:
                u2 = 0
        else:
            u = 0
            u2 = 0

        # Calculate alpha (slope angle of the failure surface) more efficiently
        delta = 0.01
        if circular:
            # For circular failure surface, use parametric equation for derivative
            # The slope at any point on the circle is: dy/dx = (x - Xo) / sqrt(R^2 - (x - Xo)^2)
            dx_circle = x_c - Xo
            if abs(dx_circle) < R:  # Check if point is on the circle
                alpha = degrees(atan(dx_circle / sqrt(R**2 - dx_circle**2)))
            else:
                # Fallback to numerical method
                y1 = get_circular_y_coordinates([x_c - delta], Xo, Yo, R)[0]
                y2 = get_circular_y_coordinates([x_c + delta], Xo, Yo, R)[0]
                alpha = degrees(atan2(y2 - y1, 2 * delta))
        else:
            # For non-circular failure surface, use geometric intersection
            failure_line = clipped_surface
            y1 = get_y_from_intersection(
                failure_line.intersection(LineString([(x_c - delta, -1e6), (x_c - delta, 1e6)])))
            y2 = get_y_from_intersection(
                failure_line.intersection(LineString([(x_c + delta, -1e6), (x_c + delta, 1e6)])))
            if y1 is not None and y2 is not None:
                alpha = degrees(atan2(y2 - y1, 2 * delta))
            else:
                alpha = 0

        if right_facing:
            alpha = -alpha
        dl = dx / cos(radians(alpha))

        if base_material_idx is None:
            phi = 0
            c = 0
            c1 = 0      # not used in rapid drawdown, but must be defined
            phi1 = 0    # not used in rapid drawdown, but must be defined
            d = 0       # not used in rapid drawdown, but must be defined
            psi = 0     # not used in rapid drawdown, but must be defined
        else:
            if materials[base_material_idx]['option'] == 'mc':
                c = materials[base_material_idx]['c']
                phi = materials[base_material_idx]['phi']
                c1 = c       # make a copy for use in rapid drawdown
                phi1 = phi   # make a copy for use in rapid drawdown
                d = materials[base_material_idx]['d']
                psi = materials[base_material_idx]['psi']
            else:
                c = (materials[base_material_idx]['r_elev'] - y_cb) * materials[base_material_idx]['cp']
                phi = 0
                c1 = 0      # not used in rapid drawdown, but must be defined
                phi1 = 0    # not used in rapid drawdown, but must be defined
                d = 0       # not used in rapid drawdown, but must be defined
                psi = 0     # not used in rapid drawdown, but must be defined

        # Prepare slice data with conditional circle parameters
        slice_data = {
            'slice #': i + 1, # Slice numbering starts at 1
            'x_l': x_l,     # left x-coordinate of the slice
            'y_lb': y_lb,   # left y-coordinate of the slice base
            'y_lt': y_lt,   # left y-coordinate of the slice top
            'x_r': x_r,     # right x-coordinate of the slice
            'y_rb': y_rb,   # right y-coordinate of the slice base
            'y_rt': y_rt,   # right y-coordinate of the slice top
            'x_c': x_c,     # center x-coordinate of the slice
            'y_cb': y_cb,   # center y-coordinate of the slice base
            'y_ct': y_ct,   # center y-coordinate of the slice top
            'y_cg': y_cg,   # center of gravity y-coordinate of the slice
            'dx': dx,       # width of the slice
            'alpha': alpha,  # slope angle of the bottom of the slice in degrees
            'dl': dl,        # length of the slice along the failure surface
            **{f'h{j+1}': h for j, h in enumerate(heights)},  # heights of each layer in the slice
            'w': soil_weight,  # weight of the slice
            'qL': qL,  # distributed load intensity at left edge
            'qR': qR,  # distributed load intensity at right edge
            'dload': dload,     # distributed load resultant (area of trapezoid)
            'd_x': d_x, # dist load resultant x-coordinate (point d)
            'd_y': d_y, # dist load resultant y-coordinate (point d)
            'qL2': qL2,  # second distributed load intensity at left edge
            'qR2': qR2,  # second distributed load intensity at right edge
            'dload2': dload2,     # second distributed load resultant (area of trapezoid)
            'd_x2': d_x2, # second dist load resultant x-coordinate (point d)
            'd_y2': d_y2, # second dist load resultant y-coordinate (point d)
            'beta': beta, # slope angle of the top edge in degrees
            'kw': kw,   # seismic force
            't': t_force,  # tension crack water force
            'y_t': y_t_loc,  # y-coordinate of the tension crack water force line of action
            'p': p_sum,   # sum of reinforcement line T values that intersect base of slice.
            'n_eff': 0, # Placeholder for effective normal force
            'z': 0,     # Placeholder for interslice side forces
            'theta': 0, # Placeholder for interslice angles
            'piezo_y': piezo_y,  # y-coordinate of the piezometric surface at x_c
            'piezo_y2': piezo_y2,  # y-coordinate of the piezometric surface at x_c for second piezometric line (rapid drawdown)
            'hw': hw,   # height of water at x_c
            'u': u,     # pore pressure at x_c
            'hw2': hw2, # height of water at x_c for second piezometric line (rapid drawdown)
            'u2': u2,   # pore pressure at x_c for second piezometric line (rapid drawdown)
            'mat': base_material_idx + 1 if base_material_idx is not None else None,  # index of the base material (1-indexed)
            'c': c,      # cohesion of the base material
            'phi': phi,  # friction angle of the base material in degrees
            'c1': c1,    # cohesion of the base material for rapid drawdown
            'phi1': phi1,  # friction angle of the base material for rapid drawdown
            'd': d,       # d cohesion of the base material for rapid drawdown
            'psi': psi,   # psi friction angle of the base material for rapid drawdown
        }
        
        # Add circle parameters only for circular failure surfaces
        if circular:
            slice_data.update({
                'r': R,      # radius of the circular failure surface
                'xo': Xo,    # x-coordinate of the center of the circular failure surface
                'yo': Yo,    # y-coordinate of the center of the circular failure surface
            })
        else:
            slice_data.update({
                'r': None,   # not applicable for non-circular failure surface
                'xo': None,  # not applicable for non-circular failure surface
                'yo': None,  # not applicable for non-circular failure surface
            })
        slices.append(slice_data)

    df = pd.DataFrame(slices)

    # Slice data were built by iterating from left to right. Flip the order slice data for right-facing slopes.
    # Slice 1 should be at the bottom and slice n at the top. This makes the slice data consistent with the
    # sign convention for alpha and the free-body diagram used to calculate forces.
    # if right_facing:
    #     df = df.iloc[::-1].reset_index(drop=True)

    return True, (df, clipped_surface)