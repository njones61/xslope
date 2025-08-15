# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from matplotlib.collections import LineCollection
import warnings


def plot_fem_data(fem_data, figsize=(14, 6), show_nodes=False, show_bc=True, material_table=False, 
                  label_elements=False, label_nodes=False, alpha=0.4, bc_symbol_size=0.03):
    """
    Plots a FEM mesh colored by material zone with boundary conditions displayed.
    
    Args:
        fem_data: Dictionary containing FEM data from build_fem_data
        figsize: Figure size
        show_nodes: If True, plot node points
        show_bc: If True, plot boundary condition symbols
        material_table: If True, show material table
        label_elements: If True, label each element with its number at its centroid
        label_nodes: If True, label each node with its number just above and to the right
        alpha: Transparency for element faces
        bc_symbol_size: Size factor for boundary condition symbols (as fraction of mesh size)
    """
    
    # Extract data from fem_data
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_materials = fem_data["element_materials"]
    element_types = fem_data.get("element_types", None)
    bc_type = fem_data["bc_type"]
    bc_values = fem_data["bc_values"]
    
    fig, ax = plt.subplots(figsize=figsize)
    materials = np.unique(element_materials)
    
    # Import get_material_color to ensure consistent colors with plot_mesh
    from plot import get_material_color
    mat_to_color = {mat: get_material_color(mat) for mat in materials}
    
    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)
    
    # Plot mesh elements with material colors
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
    material_names = fem_data.get("material_names", [])
    
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

    # Plot boundary conditions
    if show_bc:
        _plot_boundary_conditions(ax, nodes, bc_type, bc_values, legend_handles, bc_symbol_size)

    # Single combined legend outside the plot
    ax.legend(
        handles=legend_handles,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.1),
        ncol=3,  # or more, depending on how many items you have
        frameon=False
    )
    # Adjust plot limits to accommodate force arrows
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    
    # Add extra space for force arrows if they exist
    force_nodes = np.where(bc_type == 4)[0]
    if len(force_nodes) > 0:
        # Find the extent of force arrows
        mesh_size = min(x_max - x_min, y_max - y_min)
        symbol_size = mesh_size * bc_symbol_size
        
        # Add padding for force arrows (they extend outward from nodes)
        y_padding = symbol_size * 4  # Extra space above for upward arrows
        x_padding = (x_max - x_min) * 0.05  # Standard padding
        y_padding_bottom = (y_max - y_min) * 0.05
    else:
        # Standard padding
        x_padding = (x_max - x_min) * 0.05
        y_padding = (y_max - y_min) * 0.05
        y_padding_bottom = y_padding
    
    ax.set_xlim(x_min - x_padding, x_max + x_padding)
    ax.set_ylim(y_min - y_padding_bottom, y_max + y_padding)
    ax.set_aspect("equal")
    
    # Count element types for title
    num_triangles = np.sum(element_types == 3)
    num_quads = np.sum(element_types == 4)
    if num_triangles > 0 and num_quads > 0:
        title = f"FEM Mesh with Material Zones ({num_triangles} triangles, {num_quads} quads)"
    elif num_quads > 0:
        title = f"FEM Mesh with Material Zones ({num_quads} quadrilaterals)"
    else:
        title = f"FEM Mesh with Material Zones ({num_triangles} triangles)"
    
    # Place the table in the upper left
    if material_table:
        _plot_fem_material_table(ax, fem_data, xloc=0.3, yloc=1.1)  # upper left
    
    ax.set_title(title)
    plt.tight_layout()
    plt.show()


def _plot_boundary_conditions(ax, nodes, bc_type, bc_values, legend_handles, bc_symbol_size=0.03):
    """
    Plot boundary condition symbols on the mesh.
    
    BC types:
    0 = free (do nothing)
    1 = fixed (small triangle below node)
    2 = x roller (shouldn't have any)
    3 = y roller (small circle + line, left/right sides)
    4 = specified force (vector arrow)
    """
    
    # Get mesh bounds for symbol sizing
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    mesh_size = min(x_max - x_min, y_max - y_min)
    symbol_size = mesh_size * bc_symbol_size  # Adjustable symbol size
    
    # Fixed boundary conditions (type 1) - triangle below node
    fixed_nodes = np.where(bc_type == 1)[0]
    if len(fixed_nodes) > 0:
        for node_idx in fixed_nodes:
            x, y = nodes[node_idx]
            # Create small isosceles triangle below the node
            triangle_height = symbol_size
            triangle_width = symbol_size * 0.8
            triangle = patches.Polygon([
                [x - triangle_width/2, y - triangle_height],
                [x + triangle_width/2, y - triangle_height],
                [x, y]
            ], closed=True, facecolor='none', edgecolor='red', linewidth=1.5)
            ax.add_patch(triangle)
        
        # Add to legend
        legend_handles.append(
            plt.Line2D([0], [0], marker='^', color='red', linestyle='None', 
                      markersize=8, label='Fixed (bc_type=1)')
        )
    
    # Y-roller boundary conditions (type 3) - circle + line on left/right sides
    y_roller_nodes = np.where(bc_type == 3)[0]
    if len(y_roller_nodes) > 0:
        for node_idx in y_roller_nodes:
            x, y = nodes[node_idx]
            
            # Determine if node is on left or right side of mesh
            is_left_side = x < (x_min + x_max) / 2
            
            circle_radius = symbol_size * 0.4
            
            if is_left_side:
                # Put roller symbol on the left of node (circle touching node)
                circle_center_x = x - circle_radius
                line_x = circle_center_x - circle_radius
            else:
                # Put roller symbol on the right of node (circle touching node)
                circle_center_x = x + circle_radius
                line_x = circle_center_x + circle_radius
            
            # Create small hollow circle
            circle = patches.Circle((circle_center_x, y), circle_radius, 
                                  facecolor='none', edgecolor='blue', linewidth=1)
            ax.add_patch(circle)
            
            # Create tangent line
            line_length = symbol_size
            ax.plot([line_x, line_x], [y - line_length/2, y + line_length/2], 
                   'b-', linewidth=1)
        
        # Add to legend
        legend_handles.append(
            plt.Line2D([0], [0], marker='o', color='blue', linestyle='None', 
                      markersize=6, markerfacecolor='none', markeredgewidth=1, label='Y-Roller (bc_type=3)')
        )
    
    # Specified force boundary conditions (type 4) - vector arrows
    force_nodes = np.where(bc_type == 4)[0]
    if len(force_nodes) > 0:
        # Find max force magnitude for scaling
        force_magnitudes = []
        for node_idx in force_nodes:
            fx, fy = bc_values[node_idx]
            force_magnitudes.append(np.sqrt(fx**2 + fy**2))
        
        if force_magnitudes:
            max_force = max(force_magnitudes)
            if max_force > 0:
                scale = symbol_size * 3 / max_force  # Scale arrows to reasonable size
                
                for node_idx in force_nodes:
                    x, y = nodes[node_idx]
                    fx, fy = bc_values[node_idx]
                    
                    # Scale force components
                    scaled_fx = fx * scale
                    scaled_fy = fy * scale
                    
                    # Draw arrow from force end to node (so arrow points to node)
                    ax.annotate('', xy=(x, y), xytext=(x - scaled_fx, y - scaled_fy),
                               arrowprops=dict(arrowstyle='->', color='green', lw=2))
        
        # Add to legend
        legend_handles.append(
            plt.Line2D([0], [0], marker='>', color='green', linestyle='-', 
                      markersize=8, label='Applied Force (bc_type=4)')
        )


def _plot_fem_material_table(ax, fem_data, xloc=0.6, yloc=0.7):
    """
    Adds a FEM material properties table to the plot.

    Parameters:
        ax: matplotlib Axes object
        fem_data: Dictionary containing FEM data with material properties
        xloc: x-location of table (0-1)
        yloc: y-location of table (0-1)

    Returns:
        None
    """
    # Extract material properties from fem_data
    c_by_mat = fem_data.get("c_by_mat")
    phi_by_mat = fem_data.get("phi_by_mat")
    E_by_mat = fem_data.get("E_by_mat")
    nu_by_mat = fem_data.get("nu_by_mat")
    gamma_by_mat = fem_data.get("gamma_by_mat")
    material_names = fem_data.get("material_names", [])
    
    if c_by_mat is None or len(c_by_mat) == 0:
        return

    # Column headers for FEM properties
    col_labels = ["Mat", "Name", "γ", "c", "φ", "E", "ν"]

    # Build table rows
    table_data = []
    for idx in range(len(c_by_mat)):
        c = c_by_mat[idx]
        phi = phi_by_mat[idx] if phi_by_mat is not None else 0.0
        E = E_by_mat[idx] if E_by_mat is not None else 0.0
        nu = nu_by_mat[idx] if nu_by_mat is not None else 0.0
        gamma = gamma_by_mat[idx] if gamma_by_mat is not None else 0.0
        
        # Get material name, use default if not available
        material_name = material_names[idx] if idx < len(material_names) else f"Material {idx+1}"
        
        # Format values with appropriate precision
        row = [
            idx + 1,  # Material number (1-based)
            material_name,  # Material name
            f"{gamma:.1f}",   # unit weight
            f"{c:.1f}",  # cohesion
            f"{phi:.1f}",  # friction angle
            f"{E:.0f}",  # Young's modulus
            f"{nu:.2f}"  # Poisson's ratio
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


def plot_fem_results(fem_data, solution, plot_type='displacement', deform_scale=None, 
                    show_mesh=True, show_reinforcement=True, figsize=(12, 8)):
    """
    Plot FEM results with various visualization options.
    
    Parameters:
        fem_data (dict): FEM data dictionary
        solution (dict): FEM solution dictionary
        plot_type (str): Type(s) of plot. Single type ('stress', 'displacement', 'deformation') 
            or comma-separated multiple types ('stress,deformation', 'displacement,deformation').
            Multiple types are stacked vertically in the order specified.
        deform_scale (float or None): Scale factor for deformed mesh visualization.
            If None, automatically calculates scale factor so max deformation is 10% of mesh size.
            If 1.0, shows actual displacements (may be too small or too large to see).
        show_mesh (bool): Whether to show mesh lines
        show_reinforcement (bool): Whether to show reinforcement elements
        figsize (tuple): Figure size
    
    Returns:
        matplotlib figure and axes (or list of axes for multiple plots)
    """
    
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))
    
    # Parse plot types (support comma-separated list)
    plot_types = [pt.strip().lower() for pt in plot_type.split(',')]
    valid_types = ['displacement', 'deformation', 'stress']
    
    # Validate plot types
    for pt in plot_types:
        if pt not in valid_types:
            raise ValueError(f"Unknown plot_type: '{pt}'. Valid types: {valid_types}")
    
    # Auto-calculate deformation scale if not provided and deformation plots are requested
    needs_deform_scale = any(pt in ['deformation', 'stress'] for pt in plot_types)
    if deform_scale is None and needs_deform_scale:
        # Extract displacement components
        u = displacements[0::2]  # x-displacements
        v = displacements[1::2]  # y-displacements
        max_disp_mag = np.max(np.sqrt(u**2 + v**2))
        
        if max_disp_mag > 0:
            # Calculate mesh dimensions
            mesh_x_size = np.max(nodes[:, 0]) - np.min(nodes[:, 0])
            mesh_y_size = np.max(nodes[:, 1]) - np.min(nodes[:, 1])
            
            # Scale so max deformation is 10% of smallest mesh dimension
            target_deform = min(mesh_x_size, mesh_y_size) * 0.1
            deform_scale = target_deform / max_disp_mag
            
            print(f"Auto-calculated deformation scale factor: {deform_scale:.6f}")
            print(f"  Max displacement: {max_disp_mag:.4f} units")
            print(f"  Target visualization: {target_deform:.4f} units (10% of mesh size)")
        else:
            deform_scale = 1.0
            print("No displacements detected, using scale factor 1.0")
    elif deform_scale is None:
        deform_scale = 1.0  # Default for non-deformation plots
    
    # Create subplots based on number of plot types
    n_plots = len(plot_types)
    if n_plots == 1:
        fig, ax = plt.subplots(figsize=figsize)
        axes = [ax]
    else:
        # For multiple plots, adjust height scaling and use tighter spacing
        height_factor = min(0.8, 1.2 / n_plots)  # Reduce height factor for more plots
        fig, axes = plt.subplots(n_plots, 1, figsize=(figsize[0], figsize[1] * n_plots * height_factor))
        if n_plots == 1:  # Handle case where subplots returns single axis for n=1
            axes = [axes]
        

    
    # Calculate overall mesh bounds for consistent axis limits
    nodes = fem_data["nodes"]
    x_min, x_max = np.min(nodes[:, 0]), np.max(nodes[:, 0])
    y_min, y_max = np.min(nodes[:, 1]), np.max(nodes[:, 1])
    
    # Add small margin
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05
    
    # Plot each type
    for i, pt in enumerate(plot_types):
        ax = axes[i]
        
        # Calculate colorbar parameters based on number of plots
        if n_plots == 1:
            cbar_shrink = 0.8
            cbar_labelpad = 20
        elif n_plots == 2:
            cbar_shrink = 0.7  # Slightly larger than before
            cbar_labelpad = 15
        else:  # 3 or more plots
            cbar_shrink = 0.5  # Slightly larger than before
            cbar_labelpad = 12
        
        if pt == 'displacement':
            plot_displacement_contours(ax, fem_data, solution, show_mesh, show_reinforcement, 
                                     cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad)
        elif pt == 'deformation':
            plot_deformed_mesh(ax, fem_data, solution, deform_scale, show_mesh, show_reinforcement,
                             cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad)
        elif pt == 'stress':
            plot_stress_contours(ax, fem_data, solution, show_mesh, show_reinforcement,
                               cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad)
        
        # Set consistent axis limits for all plots (including single plots)
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()
    
    # Return appropriate values
    if n_plots == 1:
        return fig, axes[0]
    else:
        return fig, axes


def plot_displacement_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True, 
                              cbar_shrink=0.8, cbar_labelpad=20):
    """
    Plot displacement magnitude contours.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))
    
    # Calculate displacement magnitudes
    u = displacements[0::2]  # x-displacements
    v = displacements[1::2]  # y-displacements
    disp_mag = np.sqrt(u**2 + v**2)
    
    # Create triangulation for contouring
    triangles = []
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type == 4:  # Quad - split into triangles
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
    
    if triangles:
        triangles = np.array(triangles)
        
        # Create contour plot
        tcf = ax.tricontourf(nodes[:, 0], nodes[:, 1], triangles, disp_mag, 
                           levels=20, cmap='viridis', alpha=0.8)
        
        # Colorbar
        cbar = plt.colorbar(tcf, ax=ax, shrink=cbar_shrink)
        cbar.set_label('Displacement Magnitude', rotation=270, labelpad=cbar_labelpad)
    
    # Plot mesh
    if show_mesh:
        plot_mesh_lines(ax, fem_data, color='black', alpha=0.3, linewidth=0.5)
    
    # Plot reinforcement
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution)
    
    ax.set_aspect('equal')
    ax.set_title('Displacement Magnitude Contours')
    ax.set_xlabel('x')
    ax.set_ylabel('y')


def plot_stress_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                        cbar_shrink=0.8, cbar_labelpad=20):
    """
    Plot von Mises stress contours.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    stresses = solution.get("stresses", np.zeros((len(elements), 4)))
    plastic_elements = solution.get("plastic_elements", np.zeros(len(elements), dtype=bool))
    
    # Extract von Mises stresses
    von_mises = stresses[:, 3]  # 4th column is von Mises stress
    
    # Create element patches with color based on stress
    patches_list = []
    stress_values = []
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            coords = nodes[elem[:3]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
        elif elem_type == 4:  # Quadrilateral
            coords = nodes[elem[:4]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
    
    if patches_list:
        from matplotlib.collections import PatchCollection
        
        # Create patch collection
        p = PatchCollection(patches_list, alpha=0.8, edgecolors='none')
        p.set_array(np.array(stress_values))
        p.set_cmap('plasma')
        ax.add_collection(p)
        
        # Colorbar
        cbar = plt.colorbar(p, ax=ax, shrink=cbar_shrink)
        cbar.set_label('von Mises Stress', rotation=270, labelpad=cbar_labelpad)
    
    # Highlight plastic elements with thick boundary
    if np.any(plastic_elements):
        for i, elem in enumerate(elements):
            if plastic_elements[i]:
                elem_type = element_types[i]
                if elem_type == 3:  # Triangle
                    coords = nodes[elem[:3]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
                elif elem_type == 4:  # Quadrilateral
                    coords = nodes[elem[:4]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
    
    # Plot mesh
    if show_mesh:
        plot_mesh_lines(ax, fem_data, color='gray', alpha=0.3, linewidth=0.3)
    
    # Plot reinforcement with force visualization
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_forces(ax, fem_data, solution)
    
    ax.set_aspect('equal')
    ax.set_title('von Mises Stress (Red lines = Plastic Elements)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')


def plot_deformed_mesh(ax, fem_data, solution, deform_scale=1.0, show_mesh=True, show_reinforcement=True, 
                       cbar_shrink=0.8, cbar_labelpad=20):
    """
    Plot deformed mesh overlay on original mesh.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))
    
    # Calculate deformed node positions
    u = displacements[0::2]
    v = displacements[1::2]
    nodes_deformed = nodes + deform_scale * np.column_stack([u, v])
    
    # Plot original mesh
    if show_mesh:
        plot_mesh_lines(ax, fem_data, color='lightgray', alpha=0.5, linewidth=1.0, label='Original')
    
    # Plot deformed mesh
    fem_data_deformed = fem_data.copy()
    fem_data_deformed["nodes"] = nodes_deformed
    plot_mesh_lines(ax, fem_data_deformed, color='blue', alpha=0.8, linewidth=1.5, label='Deformed')
    
    # Plot reinforcement in both original and deformed configurations
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution, color='gray', alpha=0.5, linewidth=2, label='Original Reinforcement')
        plot_reinforcement_lines(ax, fem_data_deformed, solution, color='red', alpha=0.8, linewidth=2, label='Deformed Reinforcement')
    
    # Add a dummy colorbar to maintain consistent spacing with other plots
    # This ensures the x-axis alignment is consistent across all subplots
    dummy_data = np.array([[0, 1]])
    dummy_im = ax.imshow(dummy_data, cmap='viridis', alpha=0)
    cbar = plt.colorbar(dummy_im, ax=ax, shrink=cbar_shrink)
    cbar.set_label('Deformation Scale', rotation=270, labelpad=cbar_labelpad, color='white')
    cbar.set_ticks([])  # Remove tick marks
    cbar.set_ticklabels([])  # Remove tick labels
    
    # Make the colorbar completely invisible by setting colors to background
    cbar.outline.set_color('white')  # Make the border invisible
    cbar.outline.set_linewidth(0)    # Remove the border line
    
    # Note: Axis limits will be set by the calling function for consistent multi-plot alignment
    # When used as a standalone plot, matplotlib will auto-scale appropriately
    ax.set_title(f'Mesh Deformation (Scale Factor = {deform_scale:.1f})')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if show_mesh or show_reinforcement:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=2)


def plot_mesh_lines(ax, fem_data, color='black', alpha=1.0, linewidth=1.0, label=None):
    """
    Plot mesh element boundaries.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    lines = []
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            # Add triangle edges
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type == 4:  # Quadrilateral
            # Add quad edges
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        else:
            continue
        
        for edge in edges:
            line_coords = nodes[[edge[0], edge[1]]]
            lines.append(line_coords)
    
    if lines:
        lc = LineCollection(lines, colors=color, alpha=alpha, linewidths=linewidth, label=label)
        ax.add_collection(lc)


def plot_reinforcement_lines(ax, fem_data, solution, color='red', alpha=1.0, linewidth=2, label=None):
    """
    Plot reinforcement elements as lines.
    """
    if 'elements_1d' not in fem_data:
        return
    
    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    element_types_1d = fem_data["element_types_1d"]
    
    lines = []
    for i, elem in enumerate(elements_1d):
        elem_type = element_types_1d[i]
        if elem_type >= 2:  # At least 2 nodes
            line_coords = nodes[elem[:2]]  # Use first two nodes for line
            lines.append(line_coords)
    
    if lines:
        lc = LineCollection(lines, colors=color, alpha=alpha, linewidths=linewidth, label=label)
        ax.add_collection(lc)


def plot_reinforcement_forces(ax, fem_data, solution):
    """
    Plot reinforcement elements with color based on force magnitude.
    """
    if 'elements_1d' not in fem_data:
        return
    
    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    element_types_1d = fem_data["element_types_1d"]
    forces_1d = solution.get("forces_1d", np.zeros(len(elements_1d)))
    t_allow = fem_data.get("t_allow_by_1d_elem", np.ones(len(elements_1d)))
    failed_1d = solution.get("failed_1d_elements", np.zeros(len(elements_1d), dtype=bool))
    
    lines = []
    force_ratios = []
    
    for i, elem in enumerate(elements_1d):
        elem_type = element_types_1d[i]
        if elem_type >= 2:  # At least 2 nodes
            line_coords = nodes[elem[:2]]
            lines.append(line_coords)
            
            # Compute force ratio (force / allowable)
            if t_allow[i] > 0:
                force_ratio = abs(forces_1d[i]) / t_allow[i]
            else:
                force_ratio = 0.0
            
            # Cap at 1.5 for color scaling
            force_ratios.append(min(force_ratio, 1.5))
    
    if lines:
        # Create line collection with colors based on force ratio
        lc = LineCollection(lines, linewidths=3, alpha=0.8)
        lc.set_array(np.array(force_ratios))
        lc.set_cmap('coolwarm')  # Blue = low force, Red = high force
        ax.add_collection(lc)
        
        # Colorbar for reinforcement forces
        cbar = plt.colorbar(lc, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label('Force Ratio (Force/Allowable)', rotation=270, labelpad=15, fontsize=10)
        
        # Mark failed elements with thick black outline
        if np.any(failed_1d):
            failed_lines = [lines[i] for i in range(len(lines)) if i < len(failed_1d) and failed_1d[i]]
            if failed_lines:
                lc_failed = LineCollection(failed_lines, colors='black', linewidths=5, alpha=0.6)
                ax.add_collection(lc_failed)


def plot_reinforcement_force_profiles(fem_data, solution, figsize=(12, 8)):
    """
    Plot force profiles along each reinforcement line.
    """
    if 'elements_1d' not in fem_data:
        print("No reinforcement elements found")
        return None, None
    
    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    element_materials_1d = fem_data["element_materials_1d"]
    forces_1d = solution.get("forces_1d", np.zeros(len(elements_1d)))
    t_allow = fem_data.get("t_allow_by_1d_elem", np.ones(len(elements_1d)))
    t_res = fem_data.get("t_res_by_1d_elem", np.zeros(len(elements_1d)))
    failed_1d = solution.get("failed_1d_elements", np.zeros(len(elements_1d), dtype=bool))
    
    # Group elements by reinforcement line (material ID)
    unique_lines = np.unique(element_materials_1d)
    n_lines = len(unique_lines)
    
    if n_lines == 0:
        print("No reinforcement lines found")
        return None, None
    
    # Create subplot layout
    if n_lines <= 3:
        fig, axes = plt.subplots(n_lines, 1, figsize=figsize, squeeze=False)
        axes = axes.flatten()
    else:
        rows = int(np.ceil(n_lines / 2))
        fig, axes = plt.subplots(rows, 2, figsize=figsize, squeeze=False)
        axes = axes.flatten()
    
    for line_idx, line_id in enumerate(unique_lines):
        ax = axes[line_idx]
        
        # Get elements for this line
        line_elements = np.where(element_materials_1d == line_id)[0]
        
        if len(line_elements) == 0:
            continue
        
        # Get element positions along the line
        positions = []
        forces = []
        t_allow_line = []
        t_res_line = []
        failed_line = []
        
        for elem_idx in line_elements:
            elem = elements_1d[elem_idx]
            # Use midpoint of element
            mid_point = 0.5 * (nodes[elem[0]] + nodes[elem[1]])
            # Distance along line (simplified - use x-coordinate)
            positions.append(mid_point[0])
            forces.append(forces_1d[elem_idx])
            t_allow_line.append(t_allow[elem_idx])
            t_res_line.append(t_res[elem_idx])
            failed_line.append(failed_1d[elem_idx])
        
        # Sort by position
        sorted_indices = np.argsort(positions)
        positions = np.array(positions)[sorted_indices]
        forces = np.array(forces)[sorted_indices]
        t_allow_line = np.array(t_allow_line)[sorted_indices]
        t_res_line = np.array(t_res_line)[sorted_indices]
        failed_line = np.array(failed_line)[sorted_indices]
        
        # Plot force profile
        ax.plot(positions, forces, 'b-o', linewidth=2, markersize=6, label='Tensile Force')
        ax.plot(positions, t_allow_line, 'g--', linewidth=1, label='Allowable Force')
        
        if np.any(t_res_line > 0):
            ax.plot(positions, t_res_line, 'orange', linestyle='--', linewidth=1, label='Residual Force')
        
        # Mark failed elements
        if np.any(failed_line):
            failed_positions = positions[failed_line]
            failed_forces = forces[failed_line]
            ax.scatter(failed_positions, failed_forces, color='red', s=100, marker='x', 
                      linewidth=3, label='Failed Elements', zorder=10)
        
        # Formatting
        ax.set_xlabel('Position along line')
        ax.set_ylabel('Force')
        ax.set_title(f'Reinforcement Line {line_id} Force Profile')
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Set y-limits to show all relevant values
        max_val = max(np.max(np.abs(forces)), np.max(t_allow_line))
        if max_val > 0:
            ax.set_ylim([-max_val * 0.1, max_val * 1.1])
    
    # Hide unused subplots
    for i in range(n_lines, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    return fig, axes


def plot_ssrm_convergence(ssrm_solution, figsize=(10, 6)):
    """
    Plot SSRM convergence history.
    """
    if 'F_history' not in ssrm_solution:
        print("No SSRM convergence history found")
        return None, None
    
    F_history = ssrm_solution['F_history']
    convergence_history = ssrm_solution['convergence_history']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    
    # Plot F vs iteration
    iterations = range(1, len(F_history) + 1)
    colors = ['green' if conv else 'red' for conv in convergence_history]
    
    ax1.scatter(iterations, F_history, c=colors, s=50, alpha=0.7)
    ax1.plot(iterations, F_history, 'k-', alpha=0.5)
    
    # Mark final FS
    if 'FS' in ssrm_solution and ssrm_solution['FS'] is not None:
        ax1.axhline(y=ssrm_solution['FS'], color='blue', linestyle='--', 
                   linewidth=2, label=f"FS = {ssrm_solution['FS']:.3f}")
        ax1.legend()
    
    ax1.set_xlabel('SSRM Iteration')
    ax1.set_ylabel('Reduction Factor F')
    ax1.set_title('SSRM Convergence History')
    ax1.grid(True, alpha=0.3)
    
    # Plot convergence status
    conv_status = [1 if conv else 0 for conv in convergence_history]
    ax2.bar(iterations, conv_status, color=colors, alpha=0.7, width=0.8)
    ax2.set_xlabel('SSRM Iteration')
    ax2.set_ylabel('Converged')
    ax2.set_title('Convergence Status (Green=Converged, Red=Failed)')
    ax2.set_ylim([0, 1.2])
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig, (ax1, ax2)