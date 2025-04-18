import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint, GeometryCollection
from math import atan2, degrees, sqrt, cos, radians

def get_y_from_intersection(geom):
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

def generate_slices(profile_lines, materials, circle, surface_polyline, num_slices=20, gamma_w=62.4, piezo_line=None, dloads=None):
    Xo, Yo, depth = circle['Xo'], circle['Yo'], circle['Depth']
    R = Yo - depth

    theta_range = np.linspace(np.pi, 2 * np.pi, 1000)
    arc = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
    arc_line = LineString([(x, y) for x, y in arc])

    if surface_polyline is None or surface_polyline.is_empty:
        return pd.DataFrame(), LineString([])

    surface_polyline = LineString([(x, y) for x, y in surface_polyline.coords])

    intersections = arc_line.intersection(surface_polyline)

    intersections = arc_line.intersection(surface_polyline)

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

    if y_left > y_right:
        right_facing = True
    else:
        right_facing = False

    clipped_arc = LineString([pt for pt in arc if x_min <= pt[0] <= x_max])

    # DETERMINE SLICE LOCATIONS
    # Start with all profile x-points that intersect the arc
    fixed_xs = set(
        x for line in profile_lines for x, _ in line
        if x_min <= x <= x_max
    )
    # Explicitly include x_min and x_max
    fixed_xs.update([x_min, x_max])

    # Add the points on the distributed loads lines
    fixed_xs.update(
        pt['X'] for line in dloads for pt in line
        if x_min <= pt['X'] <= x_max
    )

    # Convert to a sorted list
    fixed_xs = sorted(fixed_xs)

    # Compute total arc span and how to divide num_slices proportionally
    segment_lengths = [fixed_xs[i + 1] - fixed_xs[i] for i in range(len(fixed_xs) - 1)]
    total_length = sum(segment_lengths)
    all_xs = [fixed_xs[0]]
    for i in range(len(fixed_xs) - 1):
        x_start = fixed_xs[i]
        x_end = fixed_xs[i + 1]
        segment_length = x_end - x_start

        # Proportional allocation
        n_subdiv = max(1, int(round((segment_length / total_length) * num_slices)))
        xs = np.linspace(x_start, x_end, n_subdiv + 1).tolist()
        all_xs.extend(xs[1:])  # skip duplicate start

    # Preprocess distributed loads for center-line interpolation
    dload_interp_funcs = []
    if dloads:
        for line in dloads:
            xs = [pt['X'] for pt in line]
            normals = [pt['Normal'] for pt in line]
            dload_interp_funcs.append(lambda x, xs=xs, normals=normals: np.interp(x, xs, normals, left=0, right=0))

    slices = []
    for i in range(len(all_xs) - 1):
        x_l, x_r = all_xs[i], all_xs[i + 1]
        x_c = (x_l + x_r) / 2
        dx = x_r - x_l

        try:
            y_cb = Yo - sqrt(R**2 - (x_c - Xo)**2)
            y_lb = Yo - sqrt(R**2 - (x_l - Xo)**2)
            y_rb = Yo - sqrt(R**2 - (x_r - Xo)**2)
        except ValueError:
            continue

        # Use vertical intersection for top of slice
        vertical_l = LineString([(x_l, surface_polyline.bounds[1] - 10), (x_l, surface_polyline.bounds[3] + 10)])
        vertical_r = LineString([(x_r, surface_polyline.bounds[1] - 10), (x_r, surface_polyline.bounds[3] + 10)])
        vertical_c = LineString([(x_c, surface_polyline.bounds[1] - 10), (x_c, surface_polyline.bounds[3] + 10)])

        y_lt = get_y_from_intersection(surface_polyline.intersection(vertical_l))
        y_rt = get_y_from_intersection(surface_polyline.intersection(vertical_r))
        y_ct = get_y_from_intersection(surface_polyline.intersection(vertical_c))

        heights = []
        total_weight = 0
        base_material_idx = None

        vertical = LineString([(x_c, surface_polyline.bounds[1] - 10), (x_c, surface_polyline.bounds[3] + 10)])

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

        # Interpolate distributed load at x_c and compute total dload
        dload_normal = 0
        if dload_interp_funcs:
            for func in dload_interp_funcs:
                dload_normal += func(x_c)

        dload = dload_normal * dx
        total_weight += dload

        hw = 0
        piezo_y = None
        if piezo_line:
            piezo_geom = LineString(piezo_line)
            piezo_vertical = LineString([(x_c, y_cb - 2 * R), (x_c, y_ct + 2 * R)])
            piezo_y = get_y_from_intersection(piezo_geom.intersection(piezo_vertical))
            if piezo_y is not None and piezo_y > y_cb:
                hw = piezo_y - y_cb
        delta = 0.01
        p1 = Point(x_c - delta, Yo - sqrt(R**2 - (x_c - delta - Xo)**2))
        p2 = Point(x_c + delta, Yo - sqrt(R**2 - (x_c + delta - Xo)**2))
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
            'piezo_y': piezo_y,
            'hw': hw,
            'u': hw * gamma_w if piezo_y is not None else 0,
            'mat': base_material_idx + 1,
            'phi': phi,
            'c': c
        }
        slices.append(slice_data)

    df = pd.DataFrame(slices)



    return df, clipped_arc