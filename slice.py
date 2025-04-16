import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point, MultiPoint, GeometryCollection
from math import atan2, degrees, sqrt

def generate_slices(profile_lines, materials, circle, piezo_line=None, surface_polyline=None, num_slices=20):
    Xo, Yo, depth = circle['Xo'], circle['Yo'], circle['Depth']
    R = Yo - depth

    theta_range = np.linspace(np.pi, 2 * np.pi, 1000)
    arc = [(Xo + R * np.cos(t), Yo + R * np.sin(t)) for t in theta_range]
    arc_line = LineString([(x, y) for x, y in arc])  # force 2D

    for pt in arc:
      print(pt)

    if surface_polyline is None or surface_polyline.is_empty:
        return pd.DataFrame(), LineString([])

    # Ensure surface is 2D
    surface_polyline = LineString([(x, y) for x, y in surface_polyline.coords])

    # Additional diagnostics
    print("arc type:", type(arc_line))
    print("arc has_z:", arc_line.has_z)
    print("surface type:", type(surface_polyline))
    print("surface has_z:", surface_polyline.has_z)

    # Debug: bounds and intersections
    intersections = arc_line.intersection(surface_polyline)
    print("arc bounds:", arc_line.bounds)
    print("surface bounds:", surface_polyline.bounds)
    print("intersection type:", type(intersections))
    print("intersection result:", intersections)

    # Diagnostic plot of arc and surface
    import matplotlib.pyplot as plt
    arc_x, arc_y = zip(*arc)
    surf_x, surf_y = zip(*surface_polyline.coords)
    plt.figure(figsize=(10, 5))
    plt.plot(arc_x, arc_y, label='Arc')
    plt.plot(surf_x, surf_y, label='Surface')
    plt.title("Arc and Surface Intersection Diagnostic")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()

    if isinstance(intersections, MultiPoint):
        x_coords = sorted(pt.x for pt in intersections.geoms)
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

        vertical_l = LineString([(x_l, surface_polyline.bounds[1] - 10), (x_l, surface_polyline.bounds[3] + 10)])
        vertical_r = LineString([(x_r, surface_polyline.bounds[1] - 10), (x_r, surface_polyline.bounds[3] + 10)])
        vertical_c = LineString([(x_c, surface_polyline.bounds[1] - 10), (x_c, surface_polyline.bounds[3] + 10)])

        inter_l = surface_polyline.intersection(vertical_l)
        inter_r = surface_polyline.intersection(vertical_r)
        inter_c = surface_polyline.intersection(vertical_c)

        y_lt = inter_l.y if isinstance(inter_l, Point) else inter_l.geoms[0].y
        y_rt = inter_r.y if isinstance(inter_r, Point) else inter_r.geoms[0].y
        y_ct = inter_c.y if isinstance(inter_c, Point) else inter_c.geoms[0].y

        slice_line = LineString([(x_c, y_ct), (x_c, y_cb)])

        heights = []
        total_weight = 0
        base_material_idx = None

        for mat_index, line in enumerate(profile_lines):
            layer_line = LineString(line)
            proj_pt = Point(x_c, 0)
            layer_top_y = layer_line.interpolate(layer_line.project(proj_pt)).y

            if mat_index + 1 < len(profile_lines):
                next_line = LineString(profile_lines[mat_index + 1])
                layer_bot_y = next_line.interpolate(next_line.project(proj_pt)).y
            else:
                layer_bot_y = y_cb

            overlap_top = min(y_ct, layer_top_y)
            overlap_bot = max(y_cb, layer_bot_y)
            h = max(0, overlap_top - overlap_bot)

            heights.append(h)
        if i in [5, 6]:
            print(f"Slice {i+1}: x_c = {x_c:.2f}, y_ct = {y_ct:.2f}, y_cb = {y_cb:.2f}")
            print(f"  Material {mat_index+1}: h = {h:.2f}")
            total_weight += h * materials[mat_index]['gamma']
            if base_material_idx is None and h > 0:
                base_material_idx = mat_index
            total_weight += h * materials[mat_index]['gamma']
            if base_material_idx is None and h > 0:
                base_material_idx = mat_index
            total_weight += h * materials[mat_index]['gamma']
            if base_material_idx is None and h > 0:
                base_material_idx = mat_index
            total_weight += h * materials[mat_index]['gamma']
            if base_material_idx is None and h > 0:
                base_material_idx = mat_index
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
