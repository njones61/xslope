import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint
from math import atan2, degrees, sqrt

def generate_slices(profile_lines, materials, circle, piezo_line=None, num_slices=20):
    Xo, Yo, depth = circle['Xo'], circle['Yo'], circle['Depth']
    R = Yo - depth  # Radius of circle

    # Define the bottom arc only
    theta_range = np.linspace(np.pi, 2 * np.pi, 1000)
    arc = [(Xo + R * np.cos(t), Yo - R * np.sin(t)) for t in theta_range]
    arc_line = LineString(arc)

    x_vals = [pt[0] for pt in arc]
    x_min, x_max = min(x_vals), max(x_vals)

    profile_xs = sorted(set(x for line in profile_lines for x, _ in line if x_min <= x <= x_max))
    uniform_xs = np.linspace(x_min, x_max, num_slices)
    all_xs = sorted(set(np.concatenate([uniform_xs, profile_xs])))

    # Construct the uppermost profile line
    top_profile_points = {}
    for line in profile_lines:
        for x, y in line:
            if x_min <= x <= x_max:
                if x not in top_profile_points or y > top_profile_points[x]:
                    top_profile_points[x] = y
    top_profile = sorted(top_profile_points.items())
    top_profile_line = LineString(top_profile)

    slices = []
    for i in range(len(all_xs) - 1):
        x1, x2 = all_xs[i], all_xs[i + 1]
        xc = (x1 + x2) / 2
        dx = x2 - x1

        try:
            base_y = Yo - sqrt(R**2 - (xc - Xo)**2)
        except ValueError:
            continue

        top_y = top_profile_line.interpolate(top_profile_line.project(Point(xc, 0))).y
        slice_line = LineString([(xc, top_y), (xc, base_y)])

        heights = []
        total_weight = 0
        base_material_idx = None

        for mat_index, line in enumerate(profile_lines):
            layer_line = LineString(line)
            layer_top_y = layer_line.interpolate(layer_line.project(Point(xc, 0))).y

            if mat_index + 1 < len(profile_lines):
                next_line = LineString(profile_lines[mat_index + 1])
                layer_bot_y = next_line.interpolate(next_line.project(Point(xc, 0))).y
            else:
                layer_bot_y = base_y

            # Calculate thickness within this layer
            top = min(top_y, layer_top_y)
            bot = max(base_y, layer_bot_y)
            h = max(0, top - bot)
            heights.append(h)
            total_weight += h * materials[mat_index]['gamma']

            if base_material_idx is None and h > 0:
                base_material_idx = mat_index

        hw = 0
        py = None
        if piezo_line:
            piezo_geom = LineString(piezo_line)
            vert_line = LineString([(xc, base_y - 2 * R), (xc, top_y + 2 * R)])
            intersection = piezo_geom.intersection(vert_line)

            piezo_y = None
            if isinstance(intersection, Point):
                piezo_y = intersection.y
            elif isinstance(intersection, MultiPoint):
                piezo_y = max(pt.y for pt in intersection.geoms)
            elif isinstance(intersection, LineString):
                if intersection.coords:
                    piezo_y = max(y for _, y in intersection.coords)

            if piezo_y is not None:
                py = piezo_y
                if piezo_y > base_y:
                    hw = piezo_y - base_y

        delta = 0.01
        p1 = Point(xc - delta, Yo - sqrt(R**2 - (xc - delta - Xo)**2))
        p2 = Point(xc + delta, Yo - sqrt(R**2 - (xc + delta - Xo)**2))
        dx_slope = p2.x - p1.x
        dy_slope = p2.y - p1.y
        alpha = degrees(atan2(dy_slope, dx_slope))

        phi = materials[base_material_idx]['phi'] if base_material_idx is not None else 0
        c = materials[base_material_idx]['c'] if base_material_idx is not None else 0

        slice_data = {
            'slice #': i + 1,
            'dx': dx,
            **{f'h{j+1}': h for j, h in enumerate(heights)},
            'w': total_weight,
            'alpha': alpha,
            'hw': hw,
            'phi': phi,
            'c': c,
            'xc': xc,
            'base_y': base_y,
            'top_y': top_y,
            'xb': xc,
            'yb': base_y,
            'py': py
        }
        slices.append(slice_data)

    df = pd.DataFrame(slices)
    return df
