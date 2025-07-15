import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.ticker import MaxNLocator
import numpy as np


def plot_seep_data(seep_data, figsize=(14, 6), show_nodes=False, show_bc=False, material_table=False, label_elements=False, label_nodes=False, alpha=0.4):
    """
    Plots a mesh colored by material zone.
    Supports both triangular and quadrilateral elements.
    
    Args:
        seep_data: Dictionary containing seepage data from import_seep2d
        show_nodes: If True, plot node points
        show_bc: If True, plot boundary condition nodes
        material_table: If True, show material table
        label_elements: If True, label each element with its number at its centroid
        label_nodes: If True, label each node with its number just above and to the right
    """

    from matplotlib.patches import Polygon

    # Extract data from seep_data
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)  # New field for element types
    bc_type = seep_data["bc_type"]

    fig, ax = plt.subplots(figsize=figsize)
    materials = np.unique(element_materials)
    
    # Import get_material_color to ensure consistent colors with plot_mesh
    from plot import get_material_color
    mat_to_color = {mat: get_material_color(mat) for mat in materials}

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        
        if element_type == 3:
            # Triangle: use first 3 nodes (4th node is repeated)
            polygon_coords = nodes[element_nodes[:3]]
        else:
            # Quad: use all 4 nodes
            polygon_coords = nodes[element_nodes]
            
        color = mat_to_color[element_materials[idx]]
        
        # Create polygon patch
        polygon = Polygon(polygon_coords, edgecolor='k', facecolor=color, linewidth=0.5, alpha=alpha)
        ax.add_patch(polygon)

        # Label element number at centroid if requested
        if label_elements:
            centroid = np.mean(polygon_coords, axis=0)
            ax.text(centroid[0], centroid[1], str(idx+1),
                    ha='center', va='center', fontsize=6, color='black', alpha=0.4,
                    zorder=10)

    if show_nodes:
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=2)

    # Label node numbers if requested
    if label_nodes:
        for i, (x, y) in enumerate(nodes):
            ax.text(x + 0.5, y + 0.5, str(i+1), fontsize=6, color='blue', alpha=0.7,
                    ha='left', va='bottom', zorder=11)

    # Get material names if available
    material_names = seep_data.get("material_names", [])
    
    legend_handles = []
    for mat in materials:
        # Use material name if available, otherwise use "Material {mat}"
        if material_names and mat <= len(material_names):
            label = material_names[mat - 1]  # Convert to 0-based index
        else:
            label = f"Material {mat}"
        
        legend_handles.append(
            plt.Line2D([0], [0], color=mat_to_color[mat], lw=4, label=label)
        )

    if show_bc:
        bc1 = nodes[bc_type == 1]
        bc2 = nodes[bc_type == 2]
        if len(bc1) > 0:
            h1, = ax.plot(bc1[:, 0], bc1[:, 1], 'ro', label="Fixed Head (bc_type=1)")
            legend_handles.append(h1)
        if len(bc2) > 0:
            h2, = ax.plot(bc2[:, 0], bc2[:, 1], 'bs', label="Exit Face (bc_type=2)")
            legend_handles.append(h2)

    # Single combined legend outside the plot
    ax.legend(
        handles=legend_handles,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.1),
        ncol=3,  # or more, depending on how many items you have
        frameon=False
    )
    ax.set_aspect("equal")
    
    # Count element types for title
    num_triangles = np.sum(element_types == 3)
    num_quads = np.sum(element_types == 4)
    if num_triangles > 0 and num_quads > 0:
        title = f"SEEP2D Mesh with Material Zones ({num_triangles} triangles, {num_quads} quads)"
    elif num_quads > 0:
        title = f"SEEP2D Mesh with Material Zones ({num_quads} quadrilaterals)"
    else:
        title = f"SEEP2D Mesh with Material Zones ({num_triangles} triangles)"
    
    # Place the table in the upper left
    if material_table:
        plot_seep_material_table(ax, seep_data, xloc=0.3, yloc=1.1)  # upper left
    
    ax.set_title(title)
    # plt.subplots_adjust(bottom=0.2)  # Add vertical cushion
    plt.tight_layout()
    plt.show()


def plot_seep_solution(seep_data, solution, figsize=(14, 6), levels=20, base_mat=1, fill_contours=True, phreatic=True, alpha=0.4, pad_frac=0.05):
    """
    Plots head contours and optionally overlays flowlines (phi) based on flow function.
    Fixed version that properly handles mesh aspect ratio and doesn't clip the plot.
    Supports both triangular and quadrilateral elements.

    Arguments:
        seep_data: Dictionary containing seepage data from import_seep2d
        solution: Dictionary containing solution results from run_analysis
        levels: number of head contour levels
        base_mat: material ID (1-based) used to compute k for flow function
        fill_contours: bool, if True shows filled contours, if False only black solid lines
        phreatic: bool, if True plots phreatic surface (pressure head = 0) as thick red line
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    from matplotlib.ticker import MaxNLocator
    from matplotlib.patches import Polygon
    import numpy as np

    # Extract data from seep_data and solution
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)  # New field for element types
    k1_by_mat = seep_data.get("k1_by_mat")  # Use .get() in case it's not present
    head = solution["head"]
    phi = solution.get("phi")
    flowrate = solution.get("flowrate")


    # Use constrained_layout for best layout
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    # Separate triangles and quads
    triangle_mask = element_types == 3
    quad_mask = element_types == 4
    
    triangle_elements = elements[triangle_mask]
    quad_elements = elements[quad_mask]
    
    print(f"Plotting {len(triangle_elements)} triangles and {len(quad_elements)} quadrilaterals")

    # Plot material zones first (if element_materials provided)
    if element_materials is not None:
        materials = np.unique(element_materials)
        
        # Import get_material_color to ensure consistent colors with plot_mesh
        from plot import get_material_color
        mat_to_color = {mat: get_material_color(mat) for mat in materials}

        # Plot triangles
        for idx, element_nodes in enumerate(triangle_elements):
            element_type = element_types[np.where(triangle_mask)[0][idx]]
            if element_type == 3:
                # Triangle: use first 3 nodes (4th node is repeated)
                polygon = nodes[element_nodes[:3]]
            else:
                # Quad: use all 4 nodes
                polygon = nodes[element_nodes]
            color = mat_to_color[element_materials[np.where(triangle_mask)[0][idx]]]
            ax.fill(*zip(*polygon), edgecolor='none', facecolor=color, alpha=alpha)
        
        # Plot quads
        for idx, element_nodes in enumerate(quad_elements):
            element_type = element_types[np.where(quad_mask)[0][idx]]
            if element_type == 3:
                # Triangle: use first 3 nodes (4th node is repeated)
                polygon_coords = nodes[element_nodes[:3]]
            else:
                # Quad: use all 4 nodes
                polygon_coords = nodes[element_nodes]
            color = mat_to_color[element_materials[np.where(quad_mask)[0][idx]]]
            polygon = Polygon(polygon_coords, edgecolor='none', facecolor=color, alpha=alpha)
            ax.add_patch(polygon)

    vmin = np.min(head)
    vmax = np.max(head)
    hdrop = vmax - vmin
    contour_levels = np.linspace(vmin, vmax, levels)

    # For contouring, we need to create triangulations
    all_triangles_for_contouring = []
    for tri_nodes in triangle_elements:
        all_triangles_for_contouring.append(tri_nodes[:3])
    for quad_nodes in quad_elements:
        tri1 = [quad_nodes[0], quad_nodes[1], quad_nodes[2]]
        tri2 = [quad_nodes[0], quad_nodes[2], quad_nodes[3]]
        all_triangles_for_contouring.extend([tri1, tri2])
    triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], all_triangles_for_contouring)

    # Filled contours (only if fill_contours=True)
    if fill_contours:
        contourf = ax.tricontourf(triang, head, levels=contour_levels, cmap="Spectral_r", vmin=vmin, vmax=vmax, alpha=0.5)
        cbar = plt.colorbar(contourf, ax=ax, label="Total Head", shrink=0.8, pad=0.02)
        cbar.locator = MaxNLocator(nbins=10, steps=[1, 2, 5])
        cbar.update_ticks()

    # Solid lines for head contours
    ax.tricontour(triang, head, levels=contour_levels, colors="k", linewidths=0.5)

    # Phreatic surface (pressure head = 0)
    if phreatic:
        elevation = nodes[:, 1]  # y-coordinate is elevation
        pressure_head = head - elevation
        ax.tricontour(triang, pressure_head, levels=[0], colors="red", linewidths=2.0)

    # Overlay flowlines if phi is available
    if phi is not None and flowrate is not None and k1_by_mat is not None:
        if base_mat > len(k1_by_mat):
            print(f"Warning: base_mat={base_mat} is larger than number of materials ({len(k1_by_mat)}). Using material 1.")
            base_mat = 1
        elif base_mat < 1:
            print(f"Warning: base_mat={base_mat} is less than 1. Using material 1.")
            base_mat = 1
        base_k = k1_by_mat[base_mat - 1]
        ne = levels - 1
        nf = (flowrate * ne) / (base_k * hdrop)
        phi_levels = round(nf) + 1
        print(f"Computed nf: {nf:.2f}, using {phi_levels} φ contours (flowrate={flowrate:.3f}, base k={base_k}, head drop={hdrop:.3f})")
        phi_contours = np.linspace(np.min(phi), np.max(phi), phi_levels)
        ax.tricontour(triang, phi, levels=phi_contours, colors="blue", linewidths=0.7, linestyles="solid")

    # Plot the mesh boundary
    try:
        boundary = get_ordered_mesh_boundary(nodes, elements, element_types)
        ax.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.0, label="Mesh Boundary")
    except Exception as e:
        print(f"Warning: Could not plot mesh boundary: {e}")

    # Add cushion around the mesh
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    x_pad = (x_max - x_min) * pad_frac
    y_pad = (y_max - y_min) * pad_frac
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    title = "Flow Net: Head Contours"
    if phi is not None:
        title += " and Flowlines"
    if phreatic:
        title += " with Phreatic Surface"
    if flowrate is not None:
        title += f" — Total Flowrate: {flowrate:.3f}"
    ax.set_title(title)

    # Set equal aspect ratio AFTER setting limits
    ax.set_aspect("equal")

    # Remove tight_layout and subplots_adjust for best constrained layout
    # plt.tight_layout()
    # plt.subplots_adjust(top=0.78)
    plt.show()


def plot_seep_material_table(ax, seep_data, xloc=0.6, yloc=0.7):
    """
    Adds a seepage material properties table to the plot.

    Parameters:
        ax: matplotlib Axes object
        seep_data: Dictionary containing seepage data with material properties
        xloc: x-location of table (0-1)
        yloc: y-location of table (0-1)

    Returns:
        None
    """
    # Extract material properties from seep_data
    k1_by_mat = seep_data.get("k1_by_mat")
    k2_by_mat = seep_data.get("k2_by_mat")
    angle_by_mat = seep_data.get("angle_by_mat")
    kr0_by_mat = seep_data.get("kr0_by_mat")
    h0_by_mat = seep_data.get("h0_by_mat")
    material_names = seep_data.get("material_names", [])
    
    if k1_by_mat is None or len(k1_by_mat) == 0:
        return

    # Column headers for seepage properties
    col_labels = ["Mat", "Name", "k₁", "k₂", "Angle", "kr₀", "h₀"]

    # Build table rows
    table_data = []
    for idx in range(len(k1_by_mat)):
        k1 = k1_by_mat[idx]
        k2 = k2_by_mat[idx] if k2_by_mat is not None else 0.0
        angle = angle_by_mat[idx] if angle_by_mat is not None else 0.0
        kr0 = kr0_by_mat[idx] if kr0_by_mat is not None else 0.0
        h0 = h0_by_mat[idx] if h0_by_mat is not None else 0.0
        
        # Get material name, use default if not available
        material_name = material_names[idx] if idx < len(material_names) else f"Material {idx+1}"
        
        # Format values with appropriate precision
        row = [
            idx + 1,  # Material number (1-based)
            material_name,  # Material name
            f"{k1:.3f}",  # k1 in scientific notation
            f"{k2:.3f}",  # k2 in scientific notation
            f"{angle:.1f}",  # angle in degrees
            f"{kr0:.4f}",  # kr0
            f"{h0:.2f}"   # h0
        ]
        table_data.append(row)

    # Add the table
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, 0.45, 0.25])  # Increased width to accommodate name column
    table.auto_set_font_size(False)
    table.set_fontsize(8)


def get_ordered_mesh_boundary(nodes, elements, element_types=None):
    """
    Extracts the outer boundary of the mesh and returns it as an ordered array of points.
    Supports both triangular and quadrilateral elements.

    Returns:
        np.ndarray of shape (N, 2): boundary coordinates in order (closed loop)
    """
    import numpy as np
    from collections import defaultdict, deque

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    # Step 1: Count all edges
    edge_count = defaultdict(int)
    edge_to_nodes = {}

    for i, element_nodes in enumerate(elements):
        element_type = element_types[i]
        
        if element_type == 3:
            # Triangle: 3 edges
            for j in range(3):
                a, b = sorted((element_nodes[j], element_nodes[(j + 1) % 3]))
                edge_count[(a, b)] += 1
                edge_to_nodes[(a, b)] = (element_nodes[j], element_nodes[(j + 1) % 3])  # preserve direction
        elif element_type == 4:
            # Quadrilateral: 4 edges
            for j in range(4):
                a, b = sorted((element_nodes[j], element_nodes[(j + 1) % 4]))
                edge_count[(a, b)] += 1
                edge_to_nodes[(a, b)] = (element_nodes[j], element_nodes[(j + 1) % 4])  # preserve direction

    # Step 2: Keep only boundary edges (appear once)
    boundary_edges = [edge_to_nodes[e] for e, count in edge_count.items() if count == 1]

    if not boundary_edges:
        raise ValueError("No boundary edges found.")

    # Step 3: Build adjacency for boundary walk
    adj = defaultdict(list)
    for a, b in boundary_edges:
        adj[a].append(b)
        adj[b].append(a)

    # Step 4: Walk all boundary segments
    all_boundary_nodes = []
    remaining_edges = set(boundary_edges)
    
    while remaining_edges:
        # Start a new boundary segment
        start_edge = remaining_edges.pop()
        start_node = start_edge[0]
        current_node = start_edge[1]
        
        segment = [start_node, current_node]
        remaining_edges.discard((current_node, start_node))  # Remove reverse edge if present
        
        # Walk this segment until we can't continue
        while True:
            # Find next edge from current node
            next_edge = None
            for edge in remaining_edges:
                if edge[0] == current_node:
                    next_edge = edge
                    break
                elif edge[1] == current_node:
                    next_edge = (edge[1], edge[0])  # Reverse the edge
                    break
            
            if next_edge is None:
                break
                
            next_node = next_edge[1]
            segment.append(next_node)
            remaining_edges.discard(next_edge)
            remaining_edges.discard((next_node, current_node))  # Remove reverse edge if present
            current_node = next_node
            
            # Check if we've closed the loop
            if current_node == start_node:
                break
        
        all_boundary_nodes.extend(segment)
    
    # If we have multiple segments, we need to handle them properly
    # For now, just return the first complete segment
    if all_boundary_nodes:
        # Ensure the boundary is closed
        if all_boundary_nodes[0] != all_boundary_nodes[-1]:
            all_boundary_nodes.append(all_boundary_nodes[0])
        return nodes[all_boundary_nodes]
    else:
        raise ValueError("No boundary nodes found.")