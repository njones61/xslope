import numpy as np
import gmsh
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from plot import get_material_color

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

def build_mesh_with_regions_OLD(polygons, region_ids, target_size, element_type='tri', debug=False):
    """
    Build a finite element mesh with material regions using Gmsh.
    
    Parameters:
        polygons     : List of lists of (x, y) tuples defining material boundaries
        region_ids   : List of material IDs (one per polygon)
        target_size  : Desired element size
        element_type : 'tri' for triangles or 'quad' for quadrilaterals (may include a few triangles)

    Returns:
        nodes        : np.ndarray of node coordinates (n_nodes, 2)
        elements     : list of element vertex indices (each element is a list of 3 or 4 ints)
        mat_ids      : list of material ID for each element
    """
    import gmsh
    import numpy as np

    if element_type not in ['tri', 'quad']:
        raise ValueError("element_type must be 'tri' or 'quad'")

    gmsh.initialize()
    gmsh.model.add("multi_region_mesh")
    point_map = {}  # maps (x, y) to Gmsh point tag

    def add_point(x, y):
        key = (x, y)
        if key not in point_map:
            tag = gmsh.model.geo.addPoint(x, y, 0, target_size)
            point_map[key] = tag
        return point_map[key]

    # Create geometry for all regions
    surface_tags = []
    region_to_physical = {region_id: idx + 1 for idx, region_id in enumerate(region_ids)}
    for idx, (poly_pts, region_id) in enumerate(zip(polygons, region_ids)):
        poly_pts_clean = remove_duplicate_endpoint(list(poly_pts))  # make a copy
        pt_tags = [add_point(x, y) for x, y in poly_pts_clean]
        line_tags = [gmsh.model.geo.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                     for i in range(len(pt_tags))]
        loop = gmsh.model.geo.addCurveLoop(line_tags)
        surface = gmsh.model.geo.addPlaneSurface([loop])
        surface_tags.append(surface)
        
        # Add physical group for this region
        physical_tag = region_to_physical[region_id]
        gmsh.model.addPhysicalGroup(2, [surface], physical_tag)

    # Synchronize geometry
    gmsh.model.geo.synchronize()
    
    # Set mesh algorithm and recombination options BEFORE generating mesh

    if element_type == 'quad':
        # Set global options for quad meshing
        gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for quads
        gmsh.option.setNumber("Mesh.RecombineAll", 1)  # Recombine triangles into quads
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # Simple recombination
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # All quads
        
        # Set recombination for each surface
        for surface in surface_to_region.keys():
            gmsh.model.mesh.setRecombine(2, surface)

    else:
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for triangles
    
    # Generate mesh
    gmsh.model.mesh.generate(2)

    # Get nodes
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)[:, :2]
    
    # Create node tag to index mapping
    node_tag_to_index = {tag: i for i, tag in enumerate(node_tags)}

    elements = []
    mat_ids = []

    # Extract elements for each region
    for region_id in region_ids:
        try:
            # Get surfaces for this physical group
            physical_tag = region_to_physical[region_id]
            phys_group_surfaces = gmsh.model.getEntitiesForPhysicalGroup(2, physical_tag)
            
            for surface in phys_group_surfaces:
                # Get all elements for this surface
                elem_types, elem_tags_list, node_tags_list = gmsh.model.mesh.getElements(2, surface)
                
                for elem_type, elem_tags, node_tags in zip(elem_types, elem_tags_list, node_tags_list):
                    if elem_type == 2:  # 3-node triangle
                        elements_array = np.array(node_tags).reshape(-1, 3)
                        for element in elements_array:
                            idxs = [node_tag_to_index[tag] for tag in element]
                            elements.append(idxs)
                            mat_ids.append(region_id)
                    elif elem_type == 3:  # 4-node quadrilateral
                        elements_array = np.array(node_tags).reshape(-1, 4)
                        for element in elements_array:
                            idxs = [node_tag_to_index[tag] for tag in element]
                            elements.append(idxs)
                            mat_ids.append(region_id)
                    elif elem_type == 9:  # 6-node triangle (second-order)
                        elements_array = np.array(node_tags).reshape(-1, 6)
                        for element in elements_array:
                            # Take only the first 3 nodes (corner nodes)
                            idxs = [node_tag_to_index[tag] for tag in element[:3]]
                            elements.append(idxs)
                            mat_ids.append(region_id)
                    elif elem_type == 10:  # 9-node quadrilateral (second-order)
                        elements_array = np.array(node_tags).reshape(-1, 9)
                        for element in elements_array:
                            # Take only the first 4 nodes (corner nodes)
                            idxs = [node_tag_to_index[tag] for tag in element[:4]]
                            elements.append(idxs)
                            mat_ids.append(region_id)
        except Exception as e:
            print(f"Warning: Could not extract elements for region {region_id}: {e}")
            continue
    
    gmsh.finalize()
    return nodes, elements, mat_ids

def build_mesh_with_regions(polygons, region_ids, target_size, element_type='tri', debug=False):
    """
    Build a finite element mesh with material regions using Gmsh.
    Alternative approach that doesn't use physical groups to avoid warnings.
    
    Parameters:
        polygons     : List of lists of (x, y) tuples defining material boundaries
        region_ids   : List of material IDs (one per polygon)
        target_size  : Desired element size
        element_type : 'tri' for triangles or 'quad' for quadrilaterals (may include a few triangles)

    Returns:
        nodes        : np.ndarray of node coordinates (n_nodes, 2)
        elements     : list of element vertex indices (each element is a list of 3 or 4 ints)
        mat_ids      : list of material ID for each element
    """
    import gmsh
    import numpy as np

    if element_type not in ['tri', 'quad']:
        raise ValueError("element_type must be 'tri' or 'quad'")

    # Adjust target_size for quads to compensate for recombination creating finer meshes
    if element_type == 'quad':
        # Increase target_size by ~1.4 (sqrt(2)) to get similar element count as triangles
        # You can adjust this factor based on your needs
        adjusted_target_size = target_size * 1.4
    else:
        adjusted_target_size = target_size

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 4)  # Reduce verbosity
    gmsh.model.add("multi_region_mesh")
    point_map = {}  # maps (x, y) to Gmsh point tag

    def add_point(x, y):
        key = (x, y)
        if key not in point_map:
            tag = gmsh.model.geo.addPoint(x, y, 0, adjusted_target_size)
            point_map[key] = tag
        return point_map[key]

    # Create geometry for all regions and track surface -> region mapping
    surface_to_region = {}
    
    for idx, (poly_pts, region_id) in enumerate(zip(polygons, region_ids)):
        poly_pts_clean = remove_duplicate_endpoint(list(poly_pts))  # make a copy
        pt_tags = [add_point(x, y) for x, y in poly_pts_clean]
        line_tags = [gmsh.model.geo.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                     for i in range(len(pt_tags))]
        loop = gmsh.model.geo.addCurveLoop(line_tags)
        surface = gmsh.model.geo.addPlaneSurface([loop])
        
        # Map surface tag directly to region ID (no physical groups needed)
        surface_to_region[surface] = region_id

    # Synchronize geometry
    gmsh.model.geo.synchronize()
    
    # Set mesh algorithm and recombination options BEFORE generating mesh
    if element_type == 'quad':
        # Set global options for quad meshing
        gmsh.option.setNumber("Mesh.Algorithm", 8)  # Frontal-Delaunay for quads
        gmsh.option.setNumber("Mesh.RecombineAll", 1)  # Recombine triangles into quads
        gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)  # Simple recombination
        gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)  # All quads
        
        # Alternative: Use more aggressive recombination to reduce element count
        # gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)  # Blossom recombination
        # gmsh.option.setNumber("Mesh.RecombineOptimize", 1)  # Optimize recombination
        
        # Set recombination for each surface
        for surface in surface_to_region.keys():
            gmsh.model.mesh.setRecombine(2, surface)
    else:
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Frontal-Delaunay for triangles
    
    # Generate mesh
    gmsh.model.mesh.generate(2)

    # Get nodes
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = np.array(coords).reshape(-1, 3)[:, :2]
    
    # Create node tag to index mapping
    node_tag_to_index = {tag: i for i, tag in enumerate(node_tags)}

    elements = []
    mat_ids = []

    # Extract elements directly from surfaces (no physical groups)
    for surface, region_id in surface_to_region.items():
        try:
            # Get all elements for this surface directly
            elem_types, elem_tags_list, node_tags_list = gmsh.model.mesh.getElements(2, surface)
            
            for elem_type, elem_tags, node_tags in zip(elem_types, elem_tags_list, node_tags_list):
                if elem_type == 2:  # 3-node triangle
                    elements_array = np.array(node_tags).reshape(-1, 3)
                    for element in elements_array:
                        idxs = [node_tag_to_index[tag] for tag in element]
                        elements.append(idxs)
                        mat_ids.append(region_id)
                elif elem_type == 3:  # 4-node quadrilateral
                    elements_array = np.array(node_tags).reshape(-1, 4)
                    for element in elements_array:
                        idxs = [node_tag_to_index[tag] for tag in element]
                        elements.append(idxs)
                        mat_ids.append(region_id)
                elif elem_type == 9:  # 6-node triangle (second-order)
                    elements_array = np.array(node_tags).reshape(-1, 6)
                    for element in elements_array:
                        # Take only the first 3 nodes (corner nodes)
                        idxs = [node_tag_to_index[tag] for tag in element[:3]]
                        elements.append(idxs)
                        mat_ids.append(region_id)
                elif elem_type == 10:  # 9-node quadrilateral (second-order)
                    elements_array = np.array(node_tags).reshape(-1, 9)
                    for element in elements_array:
                        # Take only the first 4 nodes (corner nodes)
                        idxs = [node_tag_to_index[tag] for tag in element[:4]]
                        elements.append(idxs)
                        mat_ids.append(region_id)
        except Exception as e:
            print(f"Warning: Could not extract elements for surface {surface} (region {region_id}): {e}")
            continue
    
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
def plot_mesh_with_materials(nodes, elements, mat_ids, materials=None, figsize=(14, 6), pad_frac=0.05):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.collections import PolyCollection
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Group elements by material ID for efficient plotting
    material_elements = {}
    for element, mid in zip(elements, mat_ids):
        if mid not in material_elements:
            material_elements[mid] = []
        # Get element vertices
        pts = nodes[element]
        material_elements[mid].append(pts)
    
    # Plot each material's elements as a single collection
    legend_elements = []
    for mid, elements_list in material_elements.items():
        # Create polygon collection for this material
        poly_collection = PolyCollection(elements_list, 
                                       facecolor=get_material_color(mid),
                                       edgecolor='k',
                                       alpha=0.4,
                                       linewidth=0.5)
        ax.add_collection(poly_collection)
        
        # Add to legend
        if materials and mid < len(materials) and materials[mid].get('name'):
            label = materials[mid]['name']
        else:
            label = f'Material {mid}'
        
        legend_elements.append(Patch(facecolor=get_material_color(mid), 
                                   edgecolor='k', 
                                   alpha=0.4, 
                                   label=label))
    
    ax.set_aspect('equal')
    ax.set_title("Finite Element Mesh with Material Regions (Triangles and Quads)")
    
    # Add legend if we have materials
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=min(len(legend_elements), 4))

    # Add cushion
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    x_pad = (x_max - x_min) * pad_frac
    y_pad = (y_max - y_min) * pad_frac
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    
    # Add extra cushion for legend space
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    plt.tight_layout()
    plt.show()

def build_polygons(profile_lines, max_depth=None):
    """
    For each endpoint of each profile line (except the lowest):
      - For all lower profiles, if the x is within the lower profile's x-range, compute y at that x.
      - Find the lower profile with the highest y at that x (i.e., the one directly below the endpoint).
      - If the endpoint (x, y) is coincident (within tolerance) with that lower profile at that x, but the lower profile does not have a vertex at that x, insert a new point at (x, y) at the correct position.
      - If not coincident, project vertically to that lower profile (x, y_lower(x)), and if the lower profile does not have a vertex at that x, insert it at the correct position.
    Then build polygons as before, using the (possibly augmented) profile lines.
    The cleaning step remains as a safeguard.
    """
    import numpy as np
    import copy

    if not profile_lines or len(profile_lines) < 2:
        raise ValueError("Need at least 2 profile lines to create material zones")

    def get_avg_y(line):
        return sum(y for _, y in line) / len(line)

    # Sort profile lines from top to bottom by average y
    sorted_lines = sorted(profile_lines, key=get_avg_y, reverse=True)
    n = len(sorted_lines)
    # Deep copy so we can insert points
    lines = [list(line) for line in copy.deepcopy(sorted_lines)]
    tol = 1e-8

    for i in range(n - 1):
        top = lines[i]
        for endpoint in [0, -1]:  # left and right
            x_top, y_top = top[endpoint]
            # Find the highest lower profile at this x
            best_j = None
            best_y = -np.inf
            for j in range(i + 1, n):
                lower = lines[j]
                xs_lower = np.array([x for x, y in lower])
                ys_lower = np.array([y for x, y in lower])
                if xs_lower[0] - tol <= x_top <= xs_lower[-1] + tol:
                    y_proj = np.interp(x_top, xs_lower, ys_lower)
                    if y_proj > best_y:
                        best_y = y_proj
                        best_j = j
            if best_j is not None:
                lower = lines[best_j]
                xs_lower = np.array([x for x, y in lower])
                ys_lower = np.array([y for x, y in lower])
                y_proj = np.interp(x_top, xs_lower, ys_lower)
                # Check if lower profile already has a point at this x (within tol)
                found = False
                for (x_l, y_l) in lower:
                    if abs(x_l - x_top) < tol:
                        found = True
                        break
                if abs(y_proj - y_top) < tol:
                    # Coincident: insert (x_top, y_top) if not present
                    if not found:
                        insert_idx = np.searchsorted(xs_lower, x_top)
                        lower.insert(insert_idx, (round(x_top, 6), round(y_top, 6)))
                else:
                    # Not coincident: insert (x_top, y_proj) if not present
                    if not found:
                        insert_idx = np.searchsorted(xs_lower, x_top)
                        lower.insert(insert_idx, (round(x_top, 6), round(y_proj, 6)))

    def clean_polygon(poly, tol=1e-8):
        # Remove consecutive duplicate points (except for closing point)
        if not poly:
            return poly
        cleaned = [poly[0]]
        for pt in poly[1:]:
            if abs(pt[0] - cleaned[-1][0]) > tol or abs(pt[1] - cleaned[-1][1]) > tol:
                cleaned.append(pt)
        # Ensure closed
        if abs(cleaned[0][0] - cleaned[-1][0]) > tol or abs(cleaned[0][1] - cleaned[-1][1]) > tol:
            cleaned.append(cleaned[0])
        return cleaned

    # Now build polygons as before
    polygons = []
    for i, top_line in enumerate(lines):
        xs_top, ys_top = zip(*top_line)
        xs_top = np.array(xs_top)
        ys_top = np.array(ys_top)
        left_x, left_y = xs_top[0], ys_top[0]
        right_x, right_y = xs_top[-1], ys_top[-1]

        if i < n - 1:
            lower_line = lines[i + 1]
            xs_bot, ys_bot = zip(*lower_line)
            xs_bot = np.array(xs_bot)
            ys_bot = np.array(ys_bot)
            # Project left and right endpoints vertically to lower profile
            left_y_bot = np.interp(left_x, xs_bot, ys_bot)
            right_y_bot = np.interp(right_x, xs_bot, ys_bot)
            # Find all lower profile points between left_x and right_x (exclusive)
            mask = (xs_bot > left_x) & (xs_bot < right_x)
            xs_bot_in = xs_bot[mask]
            ys_bot_in = ys_bot[mask]
            # Build bottom boundary: right projection, lower profile points (right to left), left projection
            bottom = []
            bottom.append((right_x, right_y_bot))
            for x, y in zip(xs_bot_in[::-1], ys_bot_in[::-1]):
                bottom.append((x, y))
            bottom.append((left_x, left_y_bot))
        else:
            # For the lowest polygon, bottom is at max_depth
            bottom = []
            bottom.append((right_x, max_depth))
            for x in xs_top[::-1][1:-1]:
                bottom.append((x, max_depth))
            bottom.append((left_x, max_depth))

        # Build polygon: top left-to-right, bottom right-to-left
        poly = []
        for x, y in zip(xs_top, ys_top):
            poly.append((round(x, 6), round(y, 6)))
        for x, y in bottom:
            poly.append((round(x, 6), round(y, 6)))
        # Clean up polygon (should rarely do anything)
        poly = clean_polygon(poly)
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
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 8))
    for i, polygon in enumerate(polygons):
        xs = [x for x, y in polygon]
        ys = [y for x, y in polygon]
        ax.fill(xs, ys, color=get_material_color(i), alpha=0.6, label=f'Material {i}')
        ax.plot(xs, ys, color=get_material_color(i), linewidth=1)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()

def plot_polygons_separately(polygons, title_prefix='Material Zone'):
    """
    Plot each polygon in a separate matplotlib frame (subplot), with vertices as round dots.
    """
    import matplotlib.pyplot as plt
    from plot import get_material_color
    n = len(polygons)
    fig, axes = plt.subplots(n, 1, figsize=(8, 3 * n), squeeze=False)
    for i, polygon in enumerate(polygons):
        xs = [x for x, y in polygon]
        ys = [y for x, y in polygon]
        ax = axes[i, 0]
        ax.fill(xs, ys, color=get_material_color(i), alpha=0.6, label=f'Material {i}')
        ax.plot(xs, ys, color=get_material_color(i), linewidth=1)
        ax.scatter(xs, ys, color='k', s=30, marker='o', zorder=3, label='Vertices')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'{title_prefix} {i}')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        ax.legend()
    plt.tight_layout()
    plt.show()

def remove_duplicate_endpoint(poly, tol=1e-8):
    if len(poly) > 1 and abs(poly[0][0] - poly[-1][0]) < tol and abs(poly[0][1] - poly[-1][1]) < tol:
        return poly[:-1]
    return poly