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
    from .plot import get_material_color
    mat_to_color = {mat: get_material_color(mat) for mat in materials}

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        color = mat_to_color[element_materials[idx]]
        
        if element_type == 3:  # Linear triangle
            polygon_coords = nodes[element_nodes[:3]]
            polygon = Polygon(polygon_coords, edgecolor='k', facecolor=color, linewidth=0.5, alpha=alpha)
            ax.add_patch(polygon)
            
        elif element_type == 6:  # Quadratic triangle - subdivide into 4 sub-triangles
            # Corner nodes
            n0, n1, n2 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]]
            # Midpoint nodes - standard GMSH pattern: n3=edge 0-1, n4=edge 1-2, n5=edge 2-0
            n3, n4, n5 = nodes[element_nodes[3]], nodes[element_nodes[4]], nodes[element_nodes[5]]
            
            # Create 4 sub-triangles with standard GMSH connectivity
            sub_triangles = [
                [n0, n3, n5],  # Corner triangle at node 0 (uses midpoints 0-1 and 2-0)
                [n3, n1, n4],  # Corner triangle at node 1 (uses midpoints 0-1 and 1-2)
                [n5, n4, n2],  # Corner triangle at node 2 (uses midpoints 2-0 and 1-2)
                [n3, n4, n5]   # Center triangle (connects all midpoints)
            ]
            
            # Add all sub-triangles without internal edges
            for sub_tri in sub_triangles:
                polygon = Polygon(sub_tri, edgecolor='none', facecolor=color, alpha=alpha)
                ax.add_patch(polygon)
            
            # Add outer boundary of the tri6 element
            outer_boundary = [n0, n1, n2, n0]  # Close the triangle
            ax.plot([p[0] for p in outer_boundary], [p[1] for p in outer_boundary], 
                   'k-', linewidth=0.5)
                
        elif element_type == 4:  # Linear quadrilateral
            polygon_coords = nodes[element_nodes[:4]]
            polygon = Polygon(polygon_coords, edgecolor='k', facecolor=color, linewidth=0.5, alpha=alpha)
            ax.add_patch(polygon)
            
        elif element_type == 8:  # Quadratic quadrilateral - subdivide into 4 sub-quads
            # Corner nodes
            n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
            # Midpoint nodes
            n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
            
            # Calculate center point (average of all 8 nodes)
            center = ((n0[0] + n1[0] + n2[0] + n3[0] + n4[0] + n5[0] + n6[0] + n7[0]) / 8,
                     (n0[1] + n1[1] + n2[1] + n3[1] + n4[1] + n5[1] + n6[1] + n7[1]) / 8)
            
            # Create 4 sub-quadrilaterals
            sub_quads = [
                [n0, n4, center, n7],  # Sub-quad at corner 0
                [n4, n1, n5, center],  # Sub-quad at corner 1
                [center, n5, n2, n6],  # Sub-quad at corner 2
                [n7, center, n6, n3]   # Sub-quad at corner 3
            ]
            
            # Add all sub-quads without internal edges
            for sub_quad in sub_quads:
                polygon = Polygon(sub_quad, edgecolor='none', facecolor=color, alpha=alpha)
                ax.add_patch(polygon)
            
            # Add outer boundary of the quad8 element
            outer_boundary = [n0, n1, n2, n3, n0]  # Close the quadrilateral
            ax.plot([p[0] for p in outer_boundary], [p[1] for p in outer_boundary], 
                   'k-', linewidth=0.5)
                
        elif element_type == 9:  # 9-node quadrilateral - subdivide using actual center node
            # Corner nodes
            n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
            # Midpoint nodes
            n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
            # Center node
            center = nodes[element_nodes[8]]
            
            # Create 4 sub-quadrilaterals using the actual center node
            sub_quads = [
                [n0, n4, center, n7],  # Sub-quad at corner 0
                [n4, n1, n5, center],  # Sub-quad at corner 1
                [center, n5, n2, n6],  # Sub-quad at corner 2
                [n7, center, n6, n3]   # Sub-quad at corner 3
            ]
            
            # Add all sub-quads without internal edges
            for sub_quad in sub_quads:
                polygon = Polygon(sub_quad, edgecolor='none', facecolor=color, alpha=alpha)
                ax.add_patch(polygon)
            
            # Add outer boundary of the quad9 element
            outer_boundary = [n0, n1, n2, n3, n0]  # Close the quadrilateral
            ax.plot([p[0] for p in outer_boundary], [p[1] for p in outer_boundary], 
                   'k-', linewidth=0.5)

        # Label element number at centroid if requested
        if label_elements:
            # Calculate centroid based on element type
            if element_type in [3, 4]:
                # For linear elements, use the polygon_coords
                if element_type == 3:
                    element_coords = nodes[element_nodes[:3]]
                else:
                    element_coords = nodes[element_nodes[:4]]
            else:
                # For quadratic elements, use all nodes to calculate centroid
                if element_type == 6:
                    element_coords = nodes[element_nodes[:6]]
                elif element_type == 8:
                    element_coords = nodes[element_nodes[:8]]
                else:  # element_type == 9
                    element_coords = nodes[element_nodes[:9]]
            
            centroid = np.mean(element_coords, axis=0)
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


def plot_seep_solution(seep_data, solution, figsize=(14, 6), levels=20, base_mat=1, fill_contours=True, phreatic=True, alpha=0.4, pad_frac=0.05, show_mesh=True):
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
        show_mesh: bool, if True overlays element edges in light gray
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

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)
    
    # Count element types
    tri3_count = np.sum(element_types == 3)
    tri6_count = np.sum(element_types == 6) 
    quad4_count = np.sum(element_types == 4)
    quad8_count = np.sum(element_types == 8)
    quad9_count = np.sum(element_types == 9)
    
    print(f"Plotting {tri3_count} linear triangles, {tri6_count} quadratic triangles, "
          f"{quad4_count} linear quads, {quad8_count} 8-node quads, {quad9_count} 9-node quads")

    # Plot material zones first (if element_materials provided)
    if element_materials is not None:
        materials = np.unique(element_materials)
        
        # Import get_material_color to ensure consistent colors with plot_mesh
        from .plot import get_material_color
        mat_to_color = {mat: get_material_color(mat) for mat in materials}

        # Plot all elements with proper subdivision for quadratic elements
        for idx, element_nodes in enumerate(elements):
            element_type = element_types[idx]
            color = mat_to_color[element_materials[idx]]
            
            if element_type == 3:  # Linear triangle
                polygon = nodes[element_nodes[:3]]
                ax.fill(*zip(*polygon), edgecolor='none', facecolor=color, alpha=alpha)
                
            elif element_type == 6:  # Quadratic triangle - subdivide into 4 sub-triangles
                # Corner nodes
                n0, n1, n2 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]]
                # Midpoint nodes - standard GMSH pattern: n3=edge 0-1, n4=edge 1-2, n5=edge 2-0
                n3, n4, n5 = nodes[element_nodes[3]], nodes[element_nodes[4]], nodes[element_nodes[5]]
                
                # Create 4 sub-triangles with standard GMSH connectivity
                sub_triangles = [
                    [n0, n3, n5],  # Corner triangle at node 0 (uses midpoints 0-1 and 2-0)
                    [n3, n1, n4],  # Corner triangle at node 1 (uses midpoints 0-1 and 1-2)
                    [n5, n4, n2],  # Corner triangle at node 2 (uses midpoints 2-0 and 1-2)
                    [n3, n4, n5]   # Center triangle (connects all midpoints)
                ]
                
                # Plot all sub-triangles
                for sub_tri in sub_triangles:
                    ax.fill(*zip(*sub_tri), edgecolor='none', facecolor=color, alpha=alpha)
                    
            elif element_type == 4:  # Linear quadrilateral
                polygon = nodes[element_nodes[:4]]
                ax.fill(*zip(*polygon), edgecolor='none', facecolor=color, alpha=alpha)
                
            elif element_type == 8:  # Quadratic quadrilateral - subdivide into 4 sub-quads
                # Corner nodes
                n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
                # Midpoint nodes
                n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
                
                # Calculate center point (average of all 8 nodes)
                center = ((n0[0] + n1[0] + n2[0] + n3[0] + n4[0] + n5[0] + n6[0] + n7[0]) / 8,
                         (n0[1] + n1[1] + n2[1] + n3[1] + n4[1] + n5[1] + n6[1] + n7[1]) / 8)
                
                # Create 4 sub-quadrilaterals
                sub_quads = [
                    [n0, n4, center, n7],  # Sub-quad at corner 0
                    [n4, n1, n5, center],  # Sub-quad at corner 1
                    [center, n5, n2, n6],  # Sub-quad at corner 2
                    [n7, center, n6, n3]   # Sub-quad at corner 3
                ]
                
                # Plot all sub-quads
                for sub_quad in sub_quads:
                    ax.fill(*zip(*sub_quad), edgecolor='none', facecolor=color, alpha=alpha)
                    
            elif element_type == 9:  # 9-node quadrilateral - subdivide using actual center node
                # Corner nodes
                n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
                # Midpoint nodes
                n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
                # Center node
                center = nodes[element_nodes[8]]
                
                # Create 4 sub-quadrilaterals using the actual center node
                sub_quads = [
                    [n0, n4, center, n7],  # Sub-quad at corner 0
                    [n4, n1, n5, center],  # Sub-quad at corner 1
                    [center, n5, n2, n6],  # Sub-quad at corner 2
                    [n7, center, n6, n3]   # Sub-quad at corner 3
                ]
                
                # Plot all sub-quads
                for sub_quad in sub_quads:
                    ax.fill(*zip(*sub_quad), edgecolor='none', facecolor=color, alpha=alpha)

    vmin = np.min(head)
    vmax = np.max(head)
    hdrop = vmax - vmin
    contour_levels = np.linspace(vmin, vmax, levels)

    # For contouring, subdivide tri6 elements into 4 subtriangles
    all_triangles_for_contouring = []
    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        if element_type == 3:  # Linear triangular elements
            all_triangles_for_contouring.append(element_nodes[:3])
        elif element_type == 6:  # Quadratic triangular elements
            # Standard GMSH tri6 ordering: 3 = edge 0-1; 4 = edge 1-2; 5 = edge 2-0
            # Create 4 subtriangles: 0-3-5, 3-1-4, 5-4-2, 3-4-5
            subtriangles = [
                [element_nodes[0], element_nodes[3], element_nodes[5]],  # 0-3-5 (corner at 0)
                [element_nodes[3], element_nodes[1], element_nodes[4]],  # 3-1-4 (corner at 1)
                [element_nodes[5], element_nodes[4], element_nodes[2]],  # 5-4-2 (corner at 2)
                [element_nodes[3], element_nodes[4], element_nodes[5]]   # 3-4-5 (center)
            ]
            all_triangles_for_contouring.extend(subtriangles)
        elif element_type in [4, 8, 9]:  # Quadrilateral elements
            tri1 = [element_nodes[0], element_nodes[1], element_nodes[2]]
            tri2 = [element_nodes[0], element_nodes[2], element_nodes[3]]
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

    # Plot element edges if requested
    if show_mesh:
        # Draw all element edges
        for element, elem_type in zip(elements, element_types if element_types is not None else [3]*len(elements)):
            if elem_type == 3:
                # Triangle: connect nodes 0-1-2-0
                edge_nodes = [element[0], element[1], element[2], element[0]]
            elif elem_type == 4:
                # Quadrilateral: connect nodes 0-1-2-3-0
                edge_nodes = [element[0], element[1], element[2], element[3], element[0]]
            elif elem_type == 6:
                # 6-node triangle: only connect corner nodes 0-1-2-0
                edge_nodes = [element[0], element[1], element[2], element[0]]
            elif elem_type in [8, 9]:
                # Higher-order quads: only connect corner nodes 0-1-2-3-0
                edge_nodes = [element[0], element[1], element[2], element[3], element[0]]
            else:
                continue  # Skip unknown element types
                
            # Get coordinates of edge nodes
            edge_coords = nodes[edge_nodes]
            ax.plot(edge_coords[:, 0], edge_coords[:, 1], color="darkgray", linewidth=0.5, alpha=0.7)

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
        elif element_type == 6:
            # 6-node triangle: 3 edges (use only corner nodes 0,1,2)
            for j in range(3):
                a, b = sorted((element_nodes[j], element_nodes[(j + 1) % 3]))
                edge_count[(a, b)] += 1
                edge_to_nodes[(a, b)] = (element_nodes[j], element_nodes[(j + 1) % 3])  # preserve direction
        elif element_type in [8, 9]:
            # Higher-order quadrilaterals: 4 edges (use only corner nodes 0,1,2,3)
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