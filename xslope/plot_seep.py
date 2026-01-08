import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.ticker import MaxNLocator
import numpy as np


def plot_seep_data(seep_data, figsize=(14, 6), show_nodes=False, show_bc=False, label_elements=False, label_nodes=False, alpha=0.4, save_png=False, dpi=300):
    """
    Plots a mesh colored by material zone.
    Supports both triangular and quadrilateral elements.
    
    Args:
        seep_data: Dictionary containing seep data from import_seep2d
        show_nodes: If True, plot node points
        show_bc: If True, plot boundary condition nodes
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
            h1, = ax.plot(bc1[:, 0], bc1[:, 1], 'bs', label="Fixed Head (bc_type=1)")
            legend_handles.append(h1)
        if len(bc2) > 0:
            h2, = ax.plot(bc2[:, 0], bc2[:, 1], 'ro', label="Exit Face (bc_type=2)")
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

    # Add a bit of headroom so the mesh/BC markers don't touch the top border
    y0, y1 = ax.get_ylim()
    if y1 > y0:
        pad = 0.05 * (y1 - y0)
        ax.set_ylim(y0, y1 + pad)
    
    # Count element types for title
    num_triangles = np.sum(element_types == 3)
    num_quads = np.sum(element_types == 4)
    if num_triangles > 0 and num_quads > 0:
        title = f"Finite Element Mesh with Material Zones ({num_triangles} triangles, {num_quads} quads)"
    elif num_quads > 0:
        title = f"Finite Element Mesh with Material Zones ({num_quads} quadrilaterals)"
    else:
        title = f"Finite Element Mesh with Material Zones ({num_triangles} triangles)"
    
    ax.set_title(title)
    # plt.subplots_adjust(bottom=0.2)  # Add vertical cushion
    plt.tight_layout()
    
    if save_png:
        filename = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('(', '').replace(')', '') + '.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()


def plot_seep_solution(seep_data, solution, figsize=(14, 6), levels=20, base_mat=1, fill_contours=True, phreatic=True, alpha=0.4, pad_frac=0.05, mesh=True, variable="head", vectors=False, vector_scale=0.05, flowlines=True, save_png=False, dpi=300):
    """
    Plot seep analysis results including head contours, flowlines, and phreatic surface.
    
    This function visualizes the results of a seep analysis by plotting contours of various
    nodal variables (head, pore pressure, velocity magnitude, or gradient magnitude). When
    plotting head, flowlines are also overlaid. The plot properly handles mesh aspect ratios
    and supports both linear and quadratic triangular and quadrilateral elements.
    
    Parameters:
    -----------
    seep_data : dict
        Dictionary containing seep mesh data from import_seep2d. Required keys include:
        'nodes', 'elements', 'element_materials', 'element_types' (optional), and
        'k1_by_mat' (optional, for flowline calculation).
    solution : dict
        Dictionary containing solution results from run_seepage_analysis. Required keys include:
        'head' (array of total head values at nodes), 'velocity' (array of velocity vectors),
        'gradient' (array of hydraulic gradient vectors). Optional keys: 'phi' (stream function),
        'flowrate' (total flow rate), 'u' (pore pressure), 'v_mag' (velocity magnitude),
        'i_mag' (gradient magnitude).
    figsize : tuple of float, optional
        Figure size in inches (width, height). Default is (14, 6).
    levels : int, optional
        Number of contour levels to plot. Default is 20.
    base_mat : int, optional
        Material ID (1-based) used to compute hydraulic conductivity for flow function
        calculation. Default is 1. Only used when variable="head".
    fill_contours : bool, optional
        If True, shows filled contours with color map. If False, only black solid
        contour lines are shown. Default is True.
    phreatic : bool, optional
        If True, plots the phreatic surface (where pressure head = 0) as a thick red line.
        Default is True. Only plotted if pore pressure is negative somewhere in the domain.
    alpha : float, optional
        Transparency level (0-1) for material zone fill colors. Default is 0.4.
    pad_frac : float, optional
        Fraction of mesh extent to add as padding around the plot boundaries. Default is 0.05.
    mesh : bool, optional
        If True, overlays element edges in light gray. Default is True.
    variable : str, optional
        Nodal variable to contour. Options: "head" (default), "u" (pore pressure),
        "v_mag" (velocity magnitude), "i_mag" (gradient magnitude). When "head" is selected,
        flowlines can be overlaid if flowlines=True. Other variables do not include flowlines.
    vectors : bool, optional
        If True, plots velocity vectors as arrows at each node. Default is False.
    vector_scale : float, optional
        Scale factor for vector lengths. Maximum vector length will be x_range * vector_scale,
        where x_range is the x-extent of the mesh. Default is 0.05.
    flowlines : bool, optional
        If True and variable="head", overlays flowlines (stream function contours) on the plot.
        Default is True. Only applicable when variable="head".
    
    Returns:
    --------
    None
        Displays the plot using matplotlib.pyplot.show().
    
    Notes:
    ------
    - The function automatically subdivides quadratic elements (tri6, quad8, quad9) for
      proper visualization and contouring.
    - Flowlines are only plotted when variable="head" and if 'phi' and 'flowrate' are present
      in solution and 'k1_by_mat' is present in seep_data.
    - The plot includes a colorbar for contours when fill_contours=True.
    - The title includes flowrate information if available in the solution dictionary and
      variable="head".
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    from matplotlib.ticker import MaxNLocator
    from matplotlib.patches import Polygon
    import numpy as np

    # Validate variable parameter
    valid_variables = ["head", "u", "v_mag", "i_mag"]
    if variable not in valid_variables:
        raise ValueError(f"variable must be one of {valid_variables}, got '{variable}'")
    
    # Extract data from seep_data and solution
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)  # New field for element types
    k1_by_mat = seep_data.get("k1_by_mat")  # Use .get() in case it's not present
    
    # Extract the variable to plot
    if variable not in solution:
        raise ValueError(f"Variable '{variable}' not found in solution dictionary. Available keys: {list(solution.keys())}")
    contour_data = solution[variable]
    
    # Extract head and flowline-related data (only needed for head plots)
    head = solution["head"]
    phi = solution.get("phi")
    flowrate = solution.get("flowrate")
    
    # Determine if we should plot flowlines (only for head and if flowlines=True)
    plot_flowlines = (variable == "head" and flowlines)


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

    # Set up contour levels
    vmin = np.min(contour_data)
    vmax = np.max(contour_data)
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

    # Variable labels for colorbar and title
    variable_labels = {
        "head": "Total Head",
        "u": "Pore Pressure",
        "v_mag": "Velocity Magnitude",
        "i_mag": "Hydraulic Gradient Magnitude"
    }
    variable_label = variable_labels[variable]

    # Filled contours (only if fill_contours=True)
    if fill_contours:
        contourf = ax.tricontourf(triang, contour_data, levels=contour_levels, cmap="Spectral_r", vmin=vmin, vmax=vmax, alpha=0.5)
        cbar = plt.colorbar(contourf, ax=ax, label=variable_label, shrink=0.8, pad=0.02)
        cbar.locator = MaxNLocator(nbins=10, steps=[1, 2, 5])
        cbar.update_ticks()

    # Solid lines for contours
    ax.tricontour(triang, contour_data, levels=contour_levels, colors="k", linewidths=0.5)

    # Phreatic surface (pressure head = 0)
    # Check if phreatic surface exists (pore pressure must be negative somewhere)
    has_phreatic = False
    if phreatic:
        # Check if pore pressure goes negative (indicating a phreatic surface exists)
        u = solution.get("u")
        if u is not None and np.min(u) < 0:
            elevation = nodes[:, 1]  # y-coordinate is elevation
            pressure_head = head - elevation
            ax.tricontour(triang, pressure_head, levels=[0], colors="black", linewidths=2.0)
            has_phreatic = True

    # Overlay flowlines if variable is head and phi is available
    if plot_flowlines and phi is not None and flowrate is not None and k1_by_mat is not None:
        # Compute head drop for flowline calculation
        hdrop = vmax - vmin
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

    # Plot velocity vectors if requested
    if vectors:
        velocity = solution.get("velocity")
        if velocity is not None:
            # Calculate x_range for scaling
            x_min_vec = nodes[:, 0].min()
            x_max_vec = nodes[:, 0].max()
            x_range = x_max_vec - x_min_vec
            max_vector_length = x_range * vector_scale
            
            # Get velocity magnitude
            v_mag = solution.get("v_mag")
            if v_mag is None:
                # Calculate v_mag if not available
                v_mag = np.linalg.norm(velocity, axis=1)
            
            # Find maximum velocity magnitude
            max_v_mag = np.max(v_mag)
            
            # Scale vectors: if max_v_mag > 0, scale so max vector has length max_vector_length
            if max_v_mag > 0:
                scale_factor = max_vector_length / max_v_mag
                velocity_scaled = velocity * scale_factor
            else:
                velocity_scaled = velocity
            
            # Plot vectors using quiver
            ax.quiver(nodes[:, 0], nodes[:, 1], velocity_scaled[:, 0], velocity_scaled[:, 1],
                     angles='xy', scale_units='xy', scale=1, width=0.002, headwidth=2.5,
                     headlength=3, headaxislength=2.5, color='black', alpha=0.7)

    # Plot element edges if requested
    if mesh:
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

    # Build title based on variable
    if variable == "head":
        title = f"Flow Net: {variable_label} Contours"
        if plot_flowlines and phi is not None:
            title += " and Flowlines"
        if has_phreatic:
            title += " with Phreatic Surface"
        if flowrate is not None:
            title += f" — Total Flowrate: {flowrate:.3f}"
    else:
        title = f"{variable_label} Contours"
    ax.set_title(title)

    # Set equal aspect ratio AFTER setting limits
    ax.set_aspect("equal")

    # Remove tight_layout and subplots_adjust for best constrained layout
    # plt.tight_layout()
    # plt.subplots_adjust(top=0.78)
    
    if save_png:
        filename = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('—', '').replace('(', '').replace(')', '') + '.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()


    # plot_seep_material_table has been moved to xslope/plot.py


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