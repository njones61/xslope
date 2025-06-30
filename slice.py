import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint, GeometryCollection
from math import sin, tan, atan, atan2, degrees, sqrt, cos, radians


def get_sorted_intersections(failure_surface, ground_surface):
    """
    Find and sort the intersection points between the failure and ground surfaces,
    pruning extras if the circle exits and re-enters the ground beyond the toe.

    Returns:
        success (bool), msg (str), points (list of shapely Point)
    """
    # get all the intersection geometries
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


def generate_slices(data, circle=None, non_circ=None, num_slices=40, debug=True):

    """
    Generates vertical slices between the ground surface and a failure surface for slope stability analysis.

    This function supports both circular and non-circular failure surfaces and computes
    geometric and mechanical properties for each slice, including weight, base geometry,
    water pressures, distributed loads, and reinforcement effects.

    Parameters:
        profile_lines (list): List of profile layers, each a list of (x, y) tuples.
        materials (list): List of material property dictionaries corresponding to each layer.
        ground_surface (shapely.geometry.LineString): LineString representing the top surface of the slope.
        circle (dict, optional): Dictionary with keys 'Xo', 'Yo', and 'Depth' defining the circular failure surface.
        non_circ (list, optional): List of dicts defining a non-circular failure surface with keys 'X', 'Y', and 'Movement'.
        num_slices (int, optional): Desired number of slices to generate (default is 20).
        gamma_w (float, optional): Unit weight of water (default is 62.4).
        piezo_line (list, optional): List of (x, y) tuples defining the piezometric surface.
        dloads (list, optional): List of distributed load lines, each a list of dicts with 'X', 'Y', and 'Normal'.
        dloads2 (list, optional): Second list of distributed load lines, each a list of dicts with 'X', 'Y', and 'Normal'.
        reinforce_lines (list, optional): List of reinforcement lines, each a list of dicts with 'X', 'Y', 'FL', and 'FT'.

    Returns:
        tuple:
            - pd.DataFrame: Slice table where each row includes geometry, strength, and external force values.
            - shapely.geometry.LineString: The clipped failure surface between the ground surface intersections.

    Notes:
        - Supports Method A interpretation of reinforcement: FL reduces driving forces, FT adds to normal force.
        - Handles pore pressure and distributed loads using linear interpolation at slice centers.
        - Automatically includes all geometry breakpoints in slice generation.
        - Must specify exactly one of 'circle' or 'non_circ'.
    """

    # Unpack data
    profile_lines = data["profile_lines"]
    ground_surface = data["ground_surface"]
    materials = data["materials"]
    piezo_line = data["piezo_line"]
    piezo_line2 = data.get("piezo_line2", [])  # Second piezometric line
    gamma_w = data["gamma_water"]
    tcrack_depth = data["tcrack_depth"]
    tcrack_water = data["tcrack_water"]
    k_seismic = data['k_seismic']
    if circle is not None:
        circular = True
        Xo, Yo, depth, R = circle['Xo'], circle['Yo'], circle['Depth'], circle['R']
    else:
        circular = False
    dloads = data["dloads"]
    dloads2 = data.get("dloads2", [])  # Second set of distributed loads
    max_depth = data["max_depth"]

    reinf_lines_data = []
    if data.get("reinforce_lines"):
        for line in data["reinforce_lines"]:
            xs = [pt["X"] for pt in line]
            fls = [pt["FL"] for pt in line]
            geom = LineString([(pt["X"], pt["Y"]) for pt in line])
            reinf_lines_data.append({"xs": xs, "fls": fls, "geom": geom})

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
    for line in profile_lines:
        for x, y in line:
            if x_min <= x <= x_max:
                # Check if this point is above the failure surface
                if non_circ:
                    # For non-circular failure surface, check intersection with vertical line
                    vertical_line = LineString([(x, -1e6), (x, 1e6)])
                    failure_y = get_y_from_intersection(clipped_surface.intersection(vertical_line))
                else:
                    # For circular failure surface, calculate y-coordinate
                    try:
                        failure_y = Yo - sqrt(R**2 - (x - Xo)**2)
                    except ValueError:
                        continue
                
                # Only add the point if it's above the failure surface
                if failure_y is not None and y > failure_y:
                    fixed_xs.add(x)
    
    fixed_xs.update([x_min, x_max])

    # Add transition points from dloads.
    if dloads:
        fixed_xs.update(
            pt['X'] for line in dloads for pt in line
            if x_min <= pt['X'] <= x_max
        )

    # Add transition points from dloads2.
    if dloads2:
        fixed_xs.update(
            pt['X'] for line in dloads2 for pt in line
            if x_min <= pt['X'] <= x_max
        )

    # Add transition points from non_circ.
    if non_circ:
        fixed_xs.update(
            pt['X'] for pt in non_circ
            if x_min <= pt['X'] <= x_max
        )

    # Find points associated with intersections of the profile lines and the failure surface.
    for i in range(len(profile_lines)):
        intersection = LineString(profile_lines[i]).intersection(clipped_surface)
        if not intersection.is_empty:
            if hasattr(intersection, 'x'):
                # Single point intersection
                fixed_xs.add(intersection.x)
            elif hasattr(intersection, 'geoms'):
                # Multiple points or line intersection
                for geom in intersection.geoms:
                    if hasattr(geom, 'x'):
                        fixed_xs.add(geom.x)

    fixed_xs = sorted(fixed_xs)

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
        x_l, x_r = all_xs[i], all_xs[i + 1]
        x_c = (x_l + x_r) / 2
        dx = x_r - x_l

        # Find the y-coordinates of the failure surface at the left, right, and center of the slice
        if non_circ:
            failure_line = clipped_surface
            y_cb = get_y_from_intersection(failure_line.intersection(LineString([(x_c, -1e6), (x_c, 1e6)])))
            y_lb = get_y_from_intersection(failure_line.intersection(LineString([(x_l, -1e6), (x_l, 1e6)])))
            y_rb = get_y_from_intersection(failure_line.intersection(LineString([(x_r, -1e6), (x_r, 1e6)])))
        else:
            try:
                y_cb = Yo - sqrt(R**2 - (x_c - Xo)**2)
                y_lb = Yo - sqrt(R**2 - (x_l - Xo)**2)
                y_rb = Yo - sqrt(R**2 - (x_r - Xo)**2)
            except ValueError:
                continue

        # Find the y-coordinates of the ground surface at the left, right, and center of the slice
        vertical_l = LineString([(x_l, -1e6), (x_l, 1e6)])
        vertical_r = LineString([(x_r, -1e6), (x_r, 1e6)])
        vertical_c = LineString([(x_c, -1e6), (x_c, 1e6)])

        y_lt = get_y_from_intersection(ground_surface.intersection(vertical_l))
        y_rt = get_y_from_intersection(ground_surface.intersection(vertical_r))
        y_ct = get_y_from_intersection(ground_surface.intersection(vertical_c))

        # Calculate beta (slope angle of the top edge) in degrees
        beta = degrees(atan2(y_rt - y_lt, x_r - x_l))
        if right_facing:
            beta = -beta

        # Calculate dl for the top surface (for distributed loads)
        dl_top = sqrt((x_r - x_l)**2 + (y_rt - y_lt)**2)

        heights = []
        soil_weight = 0
        base_material_idx = None

        # Make a vertical line at the center of the slice
        vertical = LineString([(x_c, -1e6), (x_c, 1e6)])

        sum_gam_h_y = 0 # for calculating center of gravity of slice
        sum_gam_h = 0   # ditto
        for mat_index, line in enumerate(profile_lines):
            layer_line = LineString(line)
            layer_top_y = get_y_from_intersection(layer_line.intersection(vertical))

            # Bottom: highest of all other profile lines at x, or failure surface
            layer_bot_y = y_cb  # Start with failure surface as default bottom
            for j in range(mat_index + 1, len(profile_lines)):
                # Check each lower profile line
                next_line = LineString(profile_lines[j])
                next_y = get_y_from_intersection(next_line.intersection(vertical))
                if next_y is not None and next_y > layer_bot_y:
                    # Take the highest of the lower profile lines
                    layer_bot_y = next_y

            if layer_top_y is None or layer_bot_y is None:
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
                # interpolated FL at xi
                fl_i = np.interp(xi, rl["xs"], rl["fls"], left=0.0, right=0.0)
                p_sum += fl_i
            else:
                # (In the extremely unlikely case that intersection is not a Point,
                #  skip it. Our assumption is only one Point per slice-base.)
                continue

        # Now p_sum is the TOTAL FL‐pull acting at this slice's base.
        # === END: "Tension crack water force" ===

        # --Process piezometric line and pore pressures---
        hw = 0
        hw2 = 0
        piezo_y = None
        piezo_y2 = None
        if piezo_line:
            piezo_geom1 = LineString(piezo_line)
            piezo_geom2 = LineString(piezo_line2)
            if circle:
                piezo_vertical = LineString([(x_c, y_cb - 2 * R), (x_c, y_ct + 2 * R)])
            else: # non-circular
                span = abs(y_ct - y_cb)
                piezo_vertical = LineString([(x_c, y_cb - 0.5 * span), (x_c, y_ct + 0.5 * span)])
            piezo_y = get_y_from_intersection(piezo_geom1.intersection(piezo_vertical))
            piezo_y2 = get_y_from_intersection(piezo_geom2.intersection(piezo_vertical))
            if piezo_y is not None and piezo_y > y_cb:
                hw = piezo_y - y_cb
            if piezo_y2 is not None and piezo_y2 > y_cb:
                hw2 = piezo_y2 - y_cb
        u = hw * gamma_w if piezo_y is not None else 0
        u2 = hw2 * gamma_w if piezo_y2 is not None else 0

        delta = 0.01
        if not non_circ:
            p1 = Point(x_c - delta, Yo - sqrt(R ** 2 - (x_c - delta - Xo) ** 2))
            p2 = Point(x_c + delta, Yo - sqrt(R ** 2 - (x_c + delta - Xo) ** 2))
        else:
            failure_line = clipped_surface
            y1 = get_y_from_intersection(
                failure_line.intersection(LineString([(x_c - delta, -1e6), (x_c - delta, 1e6)])))
            y2 = get_y_from_intersection(
                failure_line.intersection(LineString([(x_c + delta, -1e6), (x_c + delta, 1e6)])))
            p1 = Point(x_c - delta, y1)
            p2 = Point(x_c + delta, y2)

        alpha = degrees(atan2(p2.y - p1.y, p2.x - p1.x))
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
            'p': p_sum,   # sum of reinforcement line FL values that intersect base of slice.
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
            'r': R,      # radius of the circular failure surface
            'xo': Xo,    # x-coordinate of the center of the circular failure surface
            'yo': Yo,    # y-coordinate of the center of the circular failure surface
        }
        slices.append(slice_data)

    df = pd.DataFrame(slices)

    # Slice data were built by iterating from left to right. Flip the order slice data for right-facing slopes.
    # Slice 1 should be at the bottom and slice n at the top. This makes the slice data consistent with the
    # sign convention for alpha and the free-body diagram used to calculate forces.
    # if right_facing:
    #     df = df.iloc[::-1].reset_index(drop=True)

    return True, (df, clipped_surface)