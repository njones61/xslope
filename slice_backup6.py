import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint, GeometryCollection
from math import atan2, degrees, sqrt, cos, radians

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

def generate_slices(profile_lines, materials, ground_surface, *,
                    circle=None, non_circ=None, num_slices=20,
                    gamma_w=62.4, piezo_line=None, dloads=None,
                    reinforce_lines=None):
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

    if ground_surface is None or ground_surface.is_empty:
        return pd.DataFrame(), LineString([])

    ground_surface = LineString([(x, y) for x, y in ground_surface.coords])

    # Generate failure surface (either circular or non-circular)
    if non_circ:
        failure_coords = [(pt['X'], pt['Y']) for pt in non_circ]
        failure_surface = LineString(failure_coords)
    else:
        Xo, Yo, depth, R = circle['Xo'], circle['Yo'], circle['Depth'], circle['R']
        theta_range = np.linspace(np.pi, 2 * np.pi, 1000)
        arc = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
        failure_coords = arc
        failure_surface = LineString(arc)

    intersections = failure_surface.intersection(ground_surface)

    x_min = x_max = None
    y_left = y_right = None

    if isinstance(intersections, MultiPoint):
        points = sorted(intersections.geoms, key=lambda p: p.x)
        x_min, x_max = points[0].x, points[-1].x
        y_left, y_right = points[0].y, points[-1].y
    elif isinstance(intersections, Point):
        x_min = x_max = intersections.x
        y_left = y_right = intersections.y
    elif isinstance(intersections, GeometryCollection):
        points = [g for g in intersections.geoms if isinstance(g, Point)]
        if points:
            points = sorted(points, key=lambda p: p.x)
            x_min, x_max = points[0].x, points[-1].x
            y_left, y_right = points[0].y, points[-1].y
        else:
            return pd.DataFrame(), LineString([])
    else:
        return pd.DataFrame(), LineString([])

    right_facing = y_left > y_right

    clipped_surface = LineString([pt for pt in failure_coords if x_min <= pt[0] <= x_max])

    # Build fixed_xs set
    fixed_xs = set(
        x for line in profile_lines for x, _ in line
        if x_min <= x <= x_max
    )
    fixed_xs.update([x_min, x_max])

    if dloads:
        fixed_xs.update(
            pt['X'] for line in dloads for pt in line
            if x_min <= pt['X'] <= x_max
        )

    if non_circ:
        fixed_xs.update(
            pt['X'] for pt in non_circ
            if x_min <= pt['X'] <= x_max
        )

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

    # Interpolation functions for distributed loads
    dload_interp_funcs = []
    if dloads:
        for line in dloads:
            xs = [pt['X'] for pt in line]
            normals = [pt['Normal'] for pt in line]
            dload_interp_funcs.append(lambda x, xs=xs, normals=normals: np.interp(x, xs, normals, left=0, right=0))

    # Interpolation functions for reinforcement
    reinforce_interp_FL = []
    reinforce_interp_FT = []
    if reinforce_lines:
        for line in reinforce_lines:
            xs = [pt['X'] for pt in line]
            fls = [pt['FL'] for pt in line]
            fts = [pt['FT'] for pt in line]
            reinforce_interp_FL.append(lambda x, xs=xs, fls=fls: np.interp(x, xs, fls, left=0, right=0))
            reinforce_interp_FT.append(lambda x, xs=xs, fts=fts: np.interp(x, xs, fts, left=0, right=0))

    slices = []
    for i in range(len(all_xs) - 1):
        x_l, x_r = all_xs[i], all_xs[i + 1]
        x_c = (x_l + x_r) / 2
        dx = x_r - x_l

        if non_circ:
            failure_line = failure_surface
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

        vertical_l = LineString([(x_l, ground_surface.bounds[1] - 10), (x_l, ground_surface.bounds[3] + 10)])
        vertical_r = LineString([(x_r, ground_surface.bounds[1] - 10), (x_r, ground_surface.bounds[3] + 10)])
        vertical_c = LineString([(x_c, ground_surface.bounds[1] - 10), (x_c, ground_surface.bounds[3] + 10)])

        y_lt = get_y_from_intersection(ground_surface.intersection(vertical_l))
        y_rt = get_y_from_intersection(ground_surface.intersection(vertical_r))
        y_ct = get_y_from_intersection(ground_surface.intersection(vertical_c))

        heights = []
        total_weight = 0
        base_material_idx = None

        vertical = LineString([(x_c, ground_surface.bounds[1] - 10), (x_c, ground_surface.bounds[3] + 10)])

        for mat_index, line in enumerate(profile_lines):
            layer_line = LineString(line)
            layer_top_y = get_y_from_intersection(layer_line.intersection(vertical))

            if mat_index + 1 < len(profile_lines):
                next_line = LineString(profile_lines[mat_index + 1])
                layer_bot_y = get_y_from_intersection(next_line.intersection(vertical))
            else:
                layer_bot_y = y_cb

            if layer_top_y is None or layer_bot_y is None:
                h = 0
            else:
                overlap_top = min(y_ct, layer_top_y)
                overlap_bot = max(y_cb, layer_bot_y)
                h = max(0, overlap_top - overlap_bot)

            heights.append(h)
            total_weight += h * materials[mat_index]['gamma'] * dx
            if base_material_idx is None and h > 0:
                base_material_idx = mat_index

        # Distributed load
        dload_normal = sum(func(x_c) for func in dload_interp_funcs) if dload_interp_funcs else 0
        dload = dload_normal * dx
        total_weight += dload

        # Interpolate reinforcement at x_c
        shear_reinf = sum(func(x_c) for func in reinforce_interp_FL) if reinforce_interp_FL else 0
        normal_reinf = sum(func(x_c) for func in reinforce_interp_FT) if reinforce_interp_FT else 0

        hw = 0
        piezo_y = None
        if piezo_line:
            piezo_geom = LineString(piezo_line)
            if circle:
                piezo_vertical = LineString([(x_c, y_cb - 2 * R), (x_c, y_ct + 2 * R)])
            else: # non-circular
                span = abs(y_ct - y_cb)
                piezo_vertical = LineString([(x_c, y_cb - 0.5 * span), (x_c, y_ct + 0.5 * span)])
            piezo_y = get_y_from_intersection(piezo_geom.intersection(piezo_vertical))
            if piezo_y is not None and piezo_y > y_cb:
                hw = piezo_y - y_cb

        delta = 0.01
        if not non_circ:
            p1 = Point(x_c - delta, Yo - sqrt(R ** 2 - (x_c - delta - Xo) ** 2))
            p2 = Point(x_c + delta, Yo - sqrt(R ** 2 - (x_c + delta - Xo) ** 2))
        else:
            failure_line = failure_surface
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
        else:
            if materials[base_material_idx]['option'] == 'mc':
                c = materials[base_material_idx]['c']
                phi = materials[base_material_idx]['phi']
            else:
                c = (materials[base_material_idx]['r_elev'] - y_cb) * materials[base_material_idx]['cp']
                phi = 0

        slice_data = {
            'slice #': i + 1,
            'x_l': x_l,
            'y_lb': y_lb,
            'y_lt': y_lt,
            'x_r': x_r,
            'y_rb': y_rb,
            'y_rt': y_rt,
            'x_c': x_c,
            'y_cb': y_cb,
            'y_ct': y_ct,
            'dx': dx,
            'alpha': alpha,
            'dl': dl,
            **{f'h{j+1}': h for j, h in enumerate(heights)},
            'dload': dload,
            'w': total_weight,
            'shear_reinf': shear_reinf,
            'normal_reinf': normal_reinf,
            'piezo_y': piezo_y,
            'hw': hw,
            'u': hw * gamma_w if piezo_y is not None else 0,
            'mat': base_material_idx + 1 if base_material_idx is not None else None,
            'phi': phi,
            'c': c
        }
        slices.append(slice_data)

    df = pd.DataFrame(slices)
    return df, clipped_surface