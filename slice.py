import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point
from math import atan2, degrees, sqrt

def generate_slices(profile_lines, materials, circle, num_slices=20):
    Xo, Yo, depth = circle['Xo'], circle['Yo'], circle['Depth']
    R = Yo - depth  # Corrected: Radius is Yo - Depth

    # Define the circle arc from left to right
    theta_range = np.linspace(np.pi, 2 * np.pi, 1000)
    arc = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
    arc_line = LineString(arc)

    # Get x-range of the arc
    x_vals = [pt[0] for pt in arc]
    x_min, x_max = min(x_vals), max(x_vals)

    # Slice boundaries: include all profile x-values within circle + evenly spaced
    profile_xs = sorted(set(x for line in profile_lines for x, _ in line if x_min <= x <= x_max))
    uniform_xs = np.linspace(x_min, x_max, num_slices)
    all_xs = sorted(set(np.concatenate([uniform_xs, profile_xs])))

    slices = []
    for i in range(len(all_xs) - 1):
        x1, x2 = all_xs[i], all_xs[i + 1]
        xc = (x1 + x2) / 2
        dx = x2 - x1

        # Compute base point on arc
        base_y = None
        for t in theta_range:
            x = Xo + R * np.cos(t)
            if abs(x - xc) < dx / 2:
                base_y = Yo + R * np.sin(t)
                break
        if base_y is None:
            continue

        # Build vertical line through slice center
        slice_line = LineString([(xc, base_y), (xc, base_y + 1000)])

        # Determine height segments through soil layers
        heights = []
        total_weight = 0
        hw = 0  # Water height
        base_material_idx = None
        current_y = base_y

        for mat_index, line in enumerate(profile_lines):
            layer_line = LineString(line)
            int_pt = slice_line.intersection(layer_line)
            if int_pt.is_empty:
                continue
            if isinstance(int_pt, Point):
                top_y = int_pt.y
                h = top_y - current_y
                if h < 0:
                    h = 0
                heights.append(h)
                total_weight += h * materials[mat_index]['gamma']
                current_y = top_y
                if base_material_idx is None:
                    base_material_idx = mat_index
            else:
                heights.append(0)

        # Compute slope angle alpha at base center
        delta = 0.01
        p1 = Point(xc - delta, sqrt(R ** 2 - (xc - delta - Xo) ** 2) + Yo)
        p2 = Point(xc + delta, sqrt(R ** 2 - (xc + delta - Xo) ** 2) + Yo)
        dx_slope = p2.x - p1.x
        dy_slope = p2.y - p1.y
        alpha = degrees(atan2(dy_slope, dx_slope))

        # Get material at base
        phi = materials[base_material_idx]['phi'] if base_material_idx is not None else 0
        c = materials[base_material_idx]['c'] if base_material_idx is not None else 0

        slice_data = {
            'slice #': i + 1,
            'dx': dx,
            **{f'h{j + 1}': h for j, h in enumerate(heights)},
            'w': total_weight,
            'alpha': alpha,
            'hw': hw,
            'phi': phi,
            'c': c
        }
        slices.append(slice_data)

    df = pd.DataFrame(slices)
    return df