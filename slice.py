import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint, GeometryCollection
from math import atan2, degrees, sqrt

def generate_slices(profile_lines, materials, circle, piezo_line=None, num_slices=20):
    Xo, Yo, depth = circle['Xo'], circle['Yo'], circle['Depth']
    R = Yo - depth

    theta_range = np.linspace(np.pi, 2 * np.pi, 1000)
    arc = [(Xo + R * np.cos(t), Yo - R * np.sin(t)) for t in theta_range]
    arc_line = LineString(arc)

    # Step 1: Collect all points and map highest y for each x
    all_points = sorted(set(pt for line in profile_lines for pt in line))
    top_candidates = {}
    for x, y in all_points:
        if x not in top_candidates or y > top_candidates[x]:
            top_candidates[x] = y

    # Step 2: Validate each top candidate
    top_surface_points = []
    for x, y in sorted(top_candidates.items()):
        keep = True
        for other_line in profile_lines:
            line = LineString(other_line)
            if line.length == 0:
                continue
            proj = line.project(Point(x, 0))
            if proj == 0 or proj == line.length:
                continue  # avoid edge extrapolation
            ipt = line.interpolate(proj)
            if ipt.y > y + 1e-6:
                keep = False
                break
        if keep:
            top_surface_points.append((x, y))

    if len(top_surface_points) < 2:
        return pd.DataFrame(), LineString([])

    surface_polyline = LineString(top_surface_points)

    # Intersect arc with surface
    intersections = arc_line.intersection(surface_polyline)
    if isinstance(intersections, MultiPoint):
        x_coords = sorted(pt.x for pt in intersections)
        x_min, x_max = x_coords[0], x_coords[-1]
    elif isinstance(intersections, Point):
        x_min = x_max = intersections.x
    elif isinstance(intersections, GeometryCollection):
        pts = [g for g in intersections.geoms if isinstance(g, Point)]
        if pts:
            x_coords = sorted(pt.x for pt in pts)
            x_min, x_max = x_coords[0], x_coords[-1]
        else:
            return pd.DataFrame(), LineString([])
    else:
        return pd.DataFrame(), LineString([])

    clipped_arc = LineString([pt for pt in arc if x_min <= pt[0] <= x_max])
    uniform_xs = np.linspace(x_min, x_max, num_slices + 1)
    slices = []

    for i in range(len(uniform_xs) - 1):
        x_l, x_r = uniform_xs[i], uniform_xs[i + 1]
        x_c = (x_l + x_r) / 2
        dx = x_r - x_l

        try:
            y_cb = Yo - sqrt(R**2 - (x_c - Xo)**2)
            y_lb = Yo - sqrt(R**2 - (x_l - Xo)**2)
            y_rb = Yo - sqrt(R**2 - (x_r - Xo)**2)
        except ValueError:
            continue

        y_lt = surface_polyline.interpolate(surface_polyline.project(Point(x_l, 0))).y
        y_rt = surface_polyline.interpolate(surface_polyline.project(Point(x_r, 0))).y
        y_ct = surface_polyline.interpolate(surface_polyline.project(Point(x_c, 0))).y

        slice_line = LineString([(x_c, y_ct), (x_c, y_cb)])

        heights = []
        total_weight = 0
        base_material_idx = None

        for mat_index, line in enumerate(profile_lines):
            layer_line = LineString(line)
            layer_top_y = layer_line.interpolate(layer_line.project(Point(x_c, 0))).y
            if mat_index + 1 < len(profile_lines):
                next_line = LineString(profile_lines[mat_index + 1])
                layer_bot_y = next_line.interpolate(next_line.project(Point(x_c, 0))).y
            else:
                layer_bot_y = y_cb
            top = min(y_ct, layer_top_y)
            bot = max(y_cb, layer_bot_y)
            h = max(0, top - bot)
            heights.append(h)
            total_weight += h * materials[mat_index]['gamma']
            if base_material_idx is None and h > 0:
                base_material_idx = mat_index

        hw = 0
        piezo_y = None
        if piezo_line:
            piezo_geom = LineString(piezo_line)
            vert_line = LineString([(x_c, y_cb - 2 * R), (x_c, y_ct + 2 * R)])
            intersection = piezo_geom.intersection(vert_line)
            if isinstance(intersection, Point):
                piezo_y = intersection.y
            elif isinstance(intersection, MultiPoint):
                piezo_y = max(pt.y for pt in intersection.geoms)
            elif isinstance(intersection, LineString) and intersection.coords:
                piezo_y = max(y for _, y in intersection.coords)
            if piezo_y is not None and piezo_y > y_cb:
                hw = piezo_y - y_cb

        delta = 0.01
        p1 = Point(x_c - delta, Yo - sqrt(R**2 - (x_c - delta - Xo)**2))
        p2 = Point(x_c + delta, Yo - sqrt(R**2 - (x_c + delta - Xo)**2))
        alpha = degrees(atan2(p2.y - p1.y, p2.x - p1.x))

        phi = materials[base_material_idx]['phi'] if base_material_idx is not None else 0
        c = materials[base_material_idx]['c'] if base_material_idx is not None else 0

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
            **{f'h{j+1}': h for j, h in enumerate(heights)},
            'w': total_weight,
            'piezo_y': piezo_y,
            'hw': hw,
            'alpha': alpha,
            'phi': phi,
            'c': c
        }
        slices.append(slice_data)

    df = pd.DataFrame(slices)
    return df, clipped_arc