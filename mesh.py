import numpy as np
import gmsh
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee

def build_tri_mesh(points, target_size):

    gmsh.initialize()
    gmsh.model.add("paving_mesh")

    # Add points to Gmsh with target size
    point_tags = [gmsh.model.geo.addPoint(x, y, 0, target_size) for x, y in points]

    # Create lines between consecutive points (wrap last to first)
    line_tags = [
        gmsh.model.geo.addLine(point_tags[i], point_tags[(i + 1) % len(point_tags)])
        for i in range(len(point_tags))
    ]

    # Create loop and surface
    loop = gmsh.model.geo.addCurveLoop(line_tags)
    gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # Extract nodes and triangles
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(node_coords).reshape(-1, 3)[:, :2]
    elements = gmsh.model.mesh.getElementsByType(2)[1].reshape(-1, 3) - 1  # zero-based

    gmsh.finalize()
    return nodes, elements

def build_tri_mesh_with_regions(polygons, region_ids, target_size):
    """
    Parameters:
        polygons    : List of lists of (x, y) tuples
        region_ids  : List of material IDs (one per polygon)
        target_size : Desired triangle size

    Returns:
        nodes       : np.ndarray of node coordinates (n_nodes, 2)
        elements    : np.ndarray of triangle vertex indices (n_elements, 3)
        mat_ids     : np.ndarray of material ID for each triangle (n_elements,)
    """
    import gmsh
    import numpy as np

    gmsh.initialize()
    gmsh.model.add("multi_region_mesh")
    point_map = {}  # maps (x, y) to Gmsh point tag
    point_tag_counter = 1

    def add_point(x, y):
        nonlocal point_tag_counter
        key = (x, y)
        if key not in point_map:
            tag = gmsh.model.geo.addPoint(x, y, 0, target_size)
            point_map[key] = tag
        return point_map[key]

    surface_tags = []
    for poly_pts, region_id in zip(polygons, region_ids):
        pt_tags = [add_point(x, y) for x, y in poly_pts]
        line_tags = [gmsh.model.geo.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                     for i in range(len(pt_tags))]
        loop = gmsh.model.geo.addCurveLoop(line_tags)
        surface = gmsh.model.geo.addPlaneSurface([loop])
        surface_tags.append(surface)
        gmsh.model.addPhysicalGroup(2, [surface], region_id)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # Get nodes
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)[:, :2]

    # Get triangle elements and material IDs
    elements = []
    mat_ids = []

    for region_id in region_ids:
        phys_group = gmsh.model.getEntitiesForPhysicalGroup(2, region_id)
        for surface in phys_group:
            elem_type, elem_tags, node_tags = gmsh.model.mesh.getElements(2, surface)
            tri_conn = np.array(node_tags[0]).reshape(-1, 3) - 1  # zero-based
            elements.append(tri_conn)
            mat_ids.extend([region_id] * len(tri_conn))

    elements = np.vstack(elements)
    mat_ids = np.array(mat_ids)

    gmsh.finalize()
    return nodes, elements, mat_ids


def reorder_mesh(nodes, elements):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    n_nodes = nodes.shape[0]
    rows, cols = [], []

    for tri in elements:
        for i in range(3):
            for j in range(i + 1, 3):
                rows.append(tri[i])
                cols.append(tri[j])
                rows.append(tri[j])
                cols.append(tri[i])

    adj = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_nodes, n_nodes))

    def bandwidth(el):
        if len(el) == 0:
            return 0
        return max(abs(int(i) - int(j)) for tri in el for i in tri for j in tri)

    bw_before = bandwidth(elements)

    # RCM reordering
    perm = reverse_cuthill_mckee(adj.tocsr())
    inverse_perm = np.zeros_like(perm)
    inverse_perm[perm] = np.arange(len(perm))

    new_nodes = nodes[perm]
    new_elements = np.array([[inverse_perm[i] for i in tri] for tri in elements])

    # Bandwidth after
    bw_after = bandwidth(new_elements)

    print(f"Bandwidth before reordering: {bw_before}")
    print(f"Bandwidth after reordering:  {bw_after}")

    return new_nodes, new_elements


# Simple plot showing material regions
def plot_mesh_with_materials(nodes, elements, mat_ids, figsize=(14, 4), pad_frac=0.03):
    import matplotlib.pyplot as plt
    from matplotlib.cm import tab10

    fig, ax = plt.subplots(figsize=figsize)
    for tri, mid in zip(elements, mat_ids):
        pts = nodes[np.append(tri, tri[0])]
        xs, ys = pts[:, 0], pts[:, 1]
        ax.fill(xs, ys, color=tab10(mid % 10), edgecolor='k', alpha=0.6, linewidth=0.5)
    ax.set_aspect('equal')
    ax.set_title("Finite Element Mesh with Material IDs")

    # Add cushion
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    x_pad = (x_max - x_min) * pad_frac
    y_pad = (y_max - y_min) * pad_frac
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    plt.tight_layout()
    plt.show()

def build_polygons(profile_lines, max_depth=None):
    import numpy as np

    if not profile_lines or len(profile_lines) < 2:
        raise ValueError("Need at least 2 profile lines to create material zones")

    def get_avg_y(line):
        return sum(y for _, y in line) / len(line)
    # Top to bottom
    sorted_lines = sorted(profile_lines, key=get_avg_y, reverse=True)
    n = len(sorted_lines)

    def highest_lower_y(x, lower_lines):
        ys = []
        for line in lower_lines:
            xs, ys_line = zip(*line)
            if min(xs) <= x <= max(xs):
                y = np.interp(x, xs, ys_line)
                ys.append(y)
        return max(ys) if ys else None

    polygons = []
    for i, top in enumerate(sorted_lines):
        lower_lines = sorted_lines[i+1:] if i+1 < n else []
        is_last = (i == n-1)
        xs_top, ys_top = zip(*top)
        xs_top = list(xs_top)
        ys_top = list(ys_top)
        # Only endpoints
        left_x, left_y = xs_top[0], ys_top[0]
        right_x, right_y = xs_top[-1], ys_top[-1]
        # Project endpoints down
        if not is_last:
            left_y_low = highest_lower_y(left_x, lower_lines)
            right_y_low = highest_lower_y(right_x, lower_lines)
        else:
            left_y_low = right_y_low = max_depth
        # Skip if no lower boundary
        if left_y_low is None or right_y_low is None:
            continue
        # Build polygon: walk top line left to right, project right endpoint down, walk lower boundary right to left, project left endpoint up
        poly = []
        # Top line
        for x, y in zip(xs_top, ys_top):
            poly.append((x, y))
        # Project right endpoint down if not coincident
        if abs(right_y - right_y_low) > 1e-8:
            poly.append((right_x, right_y_low))
        # Walk lower boundary right to left
        if not is_last:
            # Find which lower line is highest at each endpoint
            right_y_low_line = None
            left_y_low_line = None
            for line in lower_lines:
                xs, ys_line = zip(*line)
                if min(xs) <= right_x <= max(xs) and np.isclose(np.interp(right_x, xs, ys_line), right_y_low):
                    right_y_low_line = line
                if min(xs) <= left_x <= max(xs) and np.isclose(np.interp(left_x, xs, ys_line), left_y_low):
                    left_y_low_line = line
            # Walk along the lower line from right_x to left_x
            if right_y_low_line is not None:
                xs, ys_line = zip(*right_y_low_line)
                # Get all xs between left_x and right_x (inclusive), in reverse order
                xs_between = [x for x in xs if left_x < x < right_x]
                xs_between = sorted(xs_between, reverse=True)
                for x in xs_between:
                    y = np.interp(x, xs, ys_line)
                    poly.append((x, y))
            # Project left endpoint up if not coincident
            if abs(left_y - left_y_low) > 1e-8:
                poly.append((left_x, left_y_low))
        else:
            # For the last region, walk along the base (max_depth)
            if right_x != left_x:
                xs_between = [x for x in np.linspace(right_x, left_x, num=10)][1:-1]
                for x in xs_between:
                    poly.append((x, max_depth))
            if abs(left_y - left_y_low) > 1e-8:
                poly.append((left_x, left_y_low))
        # Close polygon
        if poly[0] != poly[-1]:
            poly.append(poly[0])
        polygons.append(poly)
    return polygons


def print_polygon_summary(polygons):
    """
    Prints a summary of the generated polygons for diagnostic purposes.
    
    Parameters:
        polygons: List of polygon coordinate lists
    """
    print("=== POLYGON SUMMARY ===")
    print(f"Number of material zones: {len(polygons)}")
    print()
    
    for i, polygon in enumerate(polygons):
        print(f"Material Zone {i+1} (Material ID: {i}):")
        print(f"  Number of vertices: {len(polygon)}")
        
        # Calculate area (simple shoelace formula)
        area = 0
        for j in range(len(polygon) - 1):
            x1, y1 = polygon[j]
            x2, y2 = polygon[j + 1]
            area += (x2 - x1) * (y2 + y1) / 2
        area = abs(area)
        
        print(f"  Approximate area: {area:.2f} square units")
        
        # Print bounding box
        xs = [x for x, y in polygon]
        ys = [y for x, y in polygon]
        print(f"  Bounding box: x=[{min(xs):.2f}, {max(xs):.2f}], y=[{min(ys):.2f}, {max(ys):.2f}]")
        print()


def plot_polygons(polygons, title="Material Zone Polygons"):
    """
    Plots the generated polygons to visualize the material zones.
    
    Parameters:
        polygons: List of polygon coordinate lists
        title: Plot title
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Use a colormap for different materials
    colors = plt.cm.tab10(np.linspace(0, 1, len(polygons)))
    
    for i, polygon in enumerate(polygons):
        # Extract x and y coordinates
        xs = [x for x, y in polygon]
        ys = [y for x, y in polygon]
        
        # Plot filled polygon
        ax.fill(xs, ys, color=colors[i], alpha=0.6, 
                label=f'Material {i}')
        
        # Plot polygon boundary
        ax.plot(xs, ys, 'k-', linewidth=1)
    
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()