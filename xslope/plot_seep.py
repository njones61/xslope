import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
import numpy as np


def plot_seep_data(seep_data, figsize=(14, 6), show_nodes=False, show_bc=False, label_elements=False, label_nodes=False, alpha=0.4, save_png=False, save_dxf=False, dpi=300):
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

    # Extract data from seep_data
    nodes = seep_data["nodes"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    element_types = seep_data.get("element_types", None)
    bc_type = seep_data["bc_type"]

    fig, ax = plt.subplots(figsize=figsize)
    materials = np.unique(element_materials)

    # Import get_material_color to ensure consistent colors with plot_mesh
    from .plot import get_material_color
    mat_to_color = {mat: get_material_color(mat) for mat in materials}

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    # Batch polygons and edge lines by material for efficient rendering
    mat_fill_polys = {mat: [] for mat in materials}
    edge_segments = []

    for idx, element_nodes in enumerate(elements):
        element_type = element_types[idx]
        mat = element_materials[idx]

        if element_type == 3:  # Linear triangle
            mat_fill_polys[mat].append(nodes[element_nodes[:3]])

        elif element_type == 6:  # Quadratic triangle - subdivide into 4 sub-triangles
            n0, n1, n2 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]]
            n3, n4, n5 = nodes[element_nodes[3]], nodes[element_nodes[4]], nodes[element_nodes[5]]
            mat_fill_polys[mat].extend([
                np.array([n0, n3, n5]),
                np.array([n3, n1, n4]),
                np.array([n5, n4, n2]),
                np.array([n3, n4, n5]),
            ])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n0]])

        elif element_type == 4:  # Linear quadrilateral
            mat_fill_polys[mat].append(nodes[element_nodes[:4]])

        elif element_type == 8:  # Quadratic quadrilateral - subdivide into 4 sub-quads
            n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
            n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
            center = np.array([(n0[0]+n1[0]+n2[0]+n3[0]+n4[0]+n5[0]+n6[0]+n7[0]) / 8,
                               (n0[1]+n1[1]+n2[1]+n3[1]+n4[1]+n5[1]+n6[1]+n7[1]) / 8])
            mat_fill_polys[mat].extend([
                np.array([n0, n4, center, n7]),
                np.array([n4, n1, n5, center]),
                np.array([center, n5, n2, n6]),
                np.array([n7, center, n6, n3]),
            ])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n3], [n3, n0]])

        elif element_type == 9:  # 9-node quadrilateral
            n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
            n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
            center = nodes[element_nodes[8]]
            mat_fill_polys[mat].extend([
                np.array([n0, n4, center, n7]),
                np.array([n4, n1, n5, center]),
                np.array([center, n5, n2, n6]),
                np.array([n7, center, n6, n3]),
            ])
            edge_segments.extend([[n0, n1], [n1, n2], [n2, n3], [n3, n0]])

    # Render filled polygons as batched PatchCollections (one per material)
    for mat in materials:
        polys = mat_fill_polys[mat]
        if not polys:
            continue
        has_edge = any(element_types[i] in (3, 4) for i, m in enumerate(element_materials) if m == mat)
        has_no_edge = any(element_types[i] in (6, 8, 9) for i, m in enumerate(element_materials) if m == mat)
        color = mat_to_color[mat]

        if has_edge and not has_no_edge:
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='k', linewidth=0.5, alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        elif has_no_edge and not has_edge:
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='none', alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        else:
            linear_polys = []
            sub_polys = []
            sub_idx = 0
            for i in range(len(elements)):
                if element_materials[i] != mat:
                    continue
                et = element_types[i]
                if et in (3, 4):
                    linear_polys.append(polys[sub_idx]); sub_idx += 1
                elif et in (6, 8, 9):
                    sub_polys.extend(polys[sub_idx:sub_idx+4]); sub_idx += 4
            if linear_polys:
                pc = PatchCollection([Polygon(p) for p in linear_polys], facecolor=color, edgecolor='k', linewidth=0.5, alpha=alpha, gid='MESH_FILL')
                ax.add_collection(pc)
            if sub_polys:
                pc = PatchCollection([Polygon(p) for p in sub_polys], facecolor=color, edgecolor='none', alpha=alpha, gid='MESH_FILL')
                ax.add_collection(pc)

    # Render outer boundary edges of quadratic elements as a single LineCollection
    if edge_segments:
        lc = LineCollection(edge_segments, colors='k', linewidths=0.5, gid='MESH')
        ax.add_collection(lc)

    # Label element numbers at centroids if requested
    if label_elements:
        for idx, element_nodes in enumerate(elements):
            element_type = element_types[idx]
            if element_type == 3:
                element_coords = nodes[element_nodes[:3]]
            elif element_type == 4:
                element_coords = nodes[element_nodes[:4]]
            elif element_type == 6:
                element_coords = nodes[element_nodes[:6]]
            elif element_type == 8:
                element_coords = nodes[element_nodes[:8]]
            else:
                element_coords = nodes[element_nodes[:9]]
            centroid = np.mean(element_coords, axis=0)
            ax.text(centroid[0], centroid[1], str(idx+1),
                    ha='center', va='center', fontsize=6, color='black', alpha=0.4,
                    zorder=10)

    if show_nodes:
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=0.75, gid='MESH_NODES')

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
            h1, = ax.plot(bc1[:, 0], bc1[:, 1], 'bs', label="Fixed Head (bc_type=1)", gid='SEEP_FIXED_HEAD')
            legend_handles.append(h1)
        if len(bc2) > 0:
            h2, = ax.plot(bc2[:, 0], bc2[:, 1], 'ro', label="Exit Face (bc_type=2)", gid='SEEP_EXIT_FACE')
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
    
    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('(', '').replace(')', '')
    if save_png:
        plt.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')
    
    plt.show()


def plot_seep_solution(seep_data, solution, figsize=(14, 6), levels=20, base_mat=1, fill_contours=True, phreatic=True, alpha=0.4, pad_frac=0.05, mesh=True, variable="head", vectors=False, vector_scale=0.05, flowlines=True, save_png=False, save_dxf=False, dpi=300):
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
    k2_by_mat = seep_data.get("k2_by_mat")
    
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

        # Batch polygons by material for efficient rendering
        mat_fill_polys = {mat: [] for mat in materials}

        for idx, element_nodes in enumerate(elements):
            element_type = element_types[idx]
            mat = element_materials[idx]

            if element_type == 3:  # Linear triangle
                mat_fill_polys[mat].append(nodes[element_nodes[:3]])

            elif element_type == 6:  # Quadratic triangle - subdivide into 4 sub-triangles
                n0, n1, n2 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]]
                n3, n4, n5 = nodes[element_nodes[3]], nodes[element_nodes[4]], nodes[element_nodes[5]]
                mat_fill_polys[mat].extend([
                    np.array([n0, n3, n5]),
                    np.array([n3, n1, n4]),
                    np.array([n5, n4, n2]),
                    np.array([n3, n4, n5]),
                ])

            elif element_type == 4:  # Linear quadrilateral
                mat_fill_polys[mat].append(nodes[element_nodes[:4]])

            elif element_type == 8:  # Quadratic quadrilateral - subdivide into 4 sub-quads
                n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
                n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
                center = np.array([(n0[0]+n1[0]+n2[0]+n3[0]+n4[0]+n5[0]+n6[0]+n7[0]) / 8,
                                   (n0[1]+n1[1]+n2[1]+n3[1]+n4[1]+n5[1]+n6[1]+n7[1]) / 8])
                mat_fill_polys[mat].extend([
                    np.array([n0, n4, center, n7]),
                    np.array([n4, n1, n5, center]),
                    np.array([center, n5, n2, n6]),
                    np.array([n7, center, n6, n3]),
                ])

            elif element_type == 9:  # 9-node quadrilateral
                n0, n1, n2, n3 = nodes[element_nodes[0]], nodes[element_nodes[1]], nodes[element_nodes[2]], nodes[element_nodes[3]]
                n4, n5, n6, n7 = nodes[element_nodes[4]], nodes[element_nodes[5]], nodes[element_nodes[6]], nodes[element_nodes[7]]
                center = nodes[element_nodes[8]]
                mat_fill_polys[mat].extend([
                    np.array([n0, n4, center, n7]),
                    np.array([n4, n1, n5, center]),
                    np.array([center, n5, n2, n6]),
                    np.array([n7, center, n6, n3]),
                ])

        # Render as batched PatchCollections (one per material)
        for mat in materials:
            polys = mat_fill_polys[mat]
            if polys:
                patch_list = [Polygon(p) for p in polys]
                pc = PatchCollection(patch_list, facecolor=mat_to_color[mat], edgecolor='none', alpha=alpha, gid='ZONE_FILL')
                ax.add_collection(pc)

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
        contourf.set_gid('CONTOUR_FILL')
        cbar = plt.colorbar(contourf, ax=ax, label=variable_label, shrink=0.8, pad=0.02)
        cbar.locator = MaxNLocator(nbins=10, steps=[1, 2, 5])
        cbar.update_ticks()

    # Solid lines for contours
    _cs = ax.tricontour(triang, contour_data, levels=contour_levels, colors="k", linewidths=0.5)
    _cs.set_gid('HEAD_CONTOURS')

    # Phreatic surface (pressure head = 0)
    # Check if phreatic surface exists (pore pressure must be negative somewhere)
    has_phreatic = False
    if phreatic:
        # Check if pore pressure goes negative (indicating a phreatic surface exists)
        u = solution.get("u")
        if u is not None and np.min(u) < 0:
            elevation = nodes[:, 1]  # y-coordinate is elevation
            pressure_head = head - elevation
            _csp = ax.tricontour(triang, pressure_head, levels=[0], colors="black", linewidths=2.0)
            _csp.set_gid('PHREATIC')
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
        base_k1 = k1_by_mat[base_mat - 1]
        base_k2 = k2_by_mat[base_mat - 1] if k2_by_mat is not None else base_k1
        # For anisotropic media, the equivalent conductivity is sqrt(k1*k2)
        # This ensures the flow net cells have the correct aspect ratio sqrt(k1/k2)
        base_k = np.sqrt(base_k1 * base_k2)
        ne = levels - 1
        nf = (flowrate * ne) / (base_k * hdrop)
        phi_levels = max(round(nf) + 1, 2)
        print(f"Computed nf: {nf:.2f}, using {phi_levels} φ contours (flowrate={flowrate:.3f}, base k={base_k:.4f}, head drop={hdrop:.3f})")
        phi_contours = np.linspace(np.min(phi), np.max(phi), phi_levels)
        _csf = ax.tricontour(triang, phi, levels=phi_contours, colors="blue", linewidths=0.7, linestyles="solid")
        _csf.set_gid('FLOWLINES')

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
            _q = ax.quiver(nodes[:, 0], nodes[:, 1], velocity_scaled[:, 0], velocity_scaled[:, 1],
                     angles='xy', scale_units='xy', scale=1, width=0.002, headwidth=2.5,
                     headlength=3, headaxislength=2.5, color='black', alpha=0.7)
            _q.set_gid('VELOCITY')

    # Plot element edges if requested
    if mesh:
        mesh_segments = []
        for element, elem_type in zip(elements, element_types if element_types is not None else [3]*len(elements)):
            if elem_type in (3, 6):
                # Triangle: corner edges 0-1, 1-2, 2-0
                mesh_segments.extend([
                    [nodes[element[0]], nodes[element[1]]],
                    [nodes[element[1]], nodes[element[2]]],
                    [nodes[element[2]], nodes[element[0]]],
                ])
            elif elem_type in (4, 8, 9):
                # Quad: corner edges 0-1, 1-2, 2-3, 3-0
                mesh_segments.extend([
                    [nodes[element[0]], nodes[element[1]]],
                    [nodes[element[1]], nodes[element[2]]],
                    [nodes[element[2]], nodes[element[3]]],
                    [nodes[element[3]], nodes[element[0]]],
                ])
        if mesh_segments:
            lc = LineCollection(mesh_segments, colors='darkgray', linewidths=0.5, alpha=0.7, gid='MESH')
            ax.add_collection(lc)

    # Plot the mesh boundary
    try:
        boundary = get_ordered_mesh_boundary(nodes, elements, element_types)
        ax.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.0, label="Mesh Boundary", gid='MESH_BOUNDARY')
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
    
    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('—', '').replace('(', '').replace(')', '')
    if save_png:
        plt.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')
    
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