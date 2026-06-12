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

import warnings

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon


def _extract_uv(disp, fem_data):
    """Extract per-node u,v displacements from a mixed-DOF displacement vector."""
    dof_offset = fem_data.get("dof_offset", None)
    if dof_offset is not None:
        n_nodes = len(fem_data["nodes"])
        u = np.array([disp[dof_offset[i]] for i in range(n_nodes)])
        v = np.array([disp[dof_offset[i] + 1] for i in range(n_nodes)])
    else:
        u = disp[0::2]
        v = disp[1::2]
    return u, v


def plot_fem_data(fem_data, figsize=(14, 6), show_nodes=False, show_bc=True,
                  label_elements=False, label_nodes=False, alpha=0.4, bc_symbol_size=0.03, save_png=False, save_dxf=False, dpi=300):
    """
    Plots a FEM mesh colored by material zone with boundary conditions displayed.

    Args:
        fem_data: Dictionary containing FEM data from build_fem_data
        figsize: Figure size
        show_nodes: If True, plot node points
        show_bc: If True, plot boundary condition symbols
        label_elements: If True, label each element with its number at its centroid
        label_nodes: If True, label each node with its number just above and to the right
        alpha: Transparency for element faces
        bc_symbol_size: Size factor for boundary condition symbols (as fraction of mesh size)
    """
    from matplotlib.collections import PatchCollection

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
    from .plot import get_material_color
    mat_to_color = {mat: get_material_color(mat) for mat in materials}

    # If element_types is not provided, assume all triangles (backward compatibility)
    if element_types is None:
        element_types = np.full(len(elements), 3)

    # Batch polygons and edge lines by material for efficient rendering
    # Key: material -> list of polygon vertex arrays
    mat_fill_polys = {mat: [] for mat in materials}
    edge_segments = []  # for outer boundaries of quadratic elements

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
            # Outer boundary edges
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
            # All linear elements — draw with edges
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='k', linewidth=0.5, alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        elif has_no_edge and not has_edge:
            # All quadratic sub-polys — no edges on fills
            patch_list = [Polygon(p) for p in polys]
            pc = PatchCollection(patch_list, facecolor=color, edgecolor='none', alpha=alpha, gid='MESH_FILL')
            ax.add_collection(pc)
        else:
            # Mixed — separate linear (with edges) and quadratic sub-polys (no edges)
            linear_polys = []
            sub_polys = []
            sub_idx = 0
            for i in range(len(elements)):
                if element_materials[i] != mat:
                    continue
                et = element_types[i]
                if et == 3:
                    linear_polys.append(polys[sub_idx]); sub_idx += 1
                elif et == 4:
                    linear_polys.append(polys[sub_idx]); sub_idx += 1
                elif et == 6:
                    sub_polys.extend(polys[sub_idx:sub_idx+4]); sub_idx += 4
                elif et in (8, 9):
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
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=2, gid='MESH_NODES')

    # Label node numbers if requested
    if label_nodes:
        for i, (x, y) in enumerate(nodes):
            ax.text(x + 0.5, y + 0.5, str(i+1), fontsize=6, color='blue', alpha=0.7,
                    ha='left', va='bottom', zorder=11)

    # Get material names if available
    material_names = fem_data.get("material_names", [])

    legend_handles = []
    for mat in materials:
        if material_names and mat <= len(material_names):
            label = material_names[mat - 1]
        else:
            label = f"Material {mat}"
        legend_handles.append(
            plt.Line2D([0], [0], color=mat_to_color[mat], lw=4, label=label)
        )

    # Plot 1D elements (reinforcement truss + pile beam) using LineCollection
    elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(elements_1d), dtype=bool))
    n_reinf_plotted = 0
    n_pile_plotted = 0
    if len(elements_1d) > 0:
        reinf_segs = []
        pile_segs = []
        for elem_idx in range(len(elements_1d)):
            elem_nodes_1d = elements_1d[elem_idx]
            seg = [nodes[elem_nodes_1d[0]], nodes[elem_nodes_1d[1]]]
            if pile_elem_mask[elem_idx]:
                pile_segs.append(seg)
                n_pile_plotted += 1
            else:
                reinf_segs.append(seg)
                n_reinf_plotted += 1
        if reinf_segs:
            lc = LineCollection(reinf_segs, colors='red', linewidths=2.5, zorder=5, gid='REINFORCEMENT')
            ax.add_collection(lc)
            legend_handles.append(
                plt.Line2D([0], [0], color='red', lw=2.5, label=f'Reinforcement ({n_reinf_plotted} elements)')
            )
        if pile_segs:
            lc = LineCollection(pile_segs, colors='green', linewidths=3.5, zorder=5, gid='PILES')
            ax.add_collection(lc)
            legend_handles.append(
                plt.Line2D([0], [0], color='green', lw=3.5, label=f'Pile ({n_pile_plotted} elements)')
            )

    # Plot boundary conditions
    if show_bc:
        saved_roller_x = fem_data.get("roller_x_nodes", set())
        _plot_boundary_conditions(ax, nodes, bc_type, bc_values, legend_handles, bc_symbol_size, saved_roller_x)

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
    num_tri = np.sum((element_types == 3) | (element_types == 6))
    num_quad = np.sum((element_types == 4) | (element_types == 8) | (element_types == 9))
    num_1d = len(elements_1d)
    parts = []
    if num_tri > 0:
        parts.append(f"{num_tri} triangles")
    if num_quad > 0:
        parts.append(f"{num_quad} quads")
    if n_reinf_plotted > 0:
        parts.append(f"{n_reinf_plotted} reinforcement")
    if n_pile_plotted > 0:
        parts.append(f"{n_pile_plotted} pile")
    title = f"FEM Mesh with Material Zones ({', '.join(parts)})"
    
    ax.set_title(title)
    plt.tight_layout()
    
    base_name = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('(', '').replace(')', '')
    if save_png:
        plt.savefig(base_name + '.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        axes_to_dxf(ax, base_name + '.dxf')

    plt.show()


def _plot_boundary_conditions(ax, nodes, bc_type, bc_values, legend_handles, bc_symbol_size=0.03, saved_roller_x=None):
    """
    Plot boundary condition symbols on the mesh.

    BC types:
    0 = free (do nothing)
    1 = fixed (small triangle below node)
    2 = x roller (small circle + line, left/right sides)
    3 = y roller (shouldn't have any)
    4 = specified force (vector arrow)
    """
    from matplotlib.collections import PatchCollection

    # Get mesh bounds for symbol sizing
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    mesh_size = min(x_max - x_min, y_max - y_min)
    symbol_size = mesh_size * bc_symbol_size

    # Fixed boundary conditions (type 1) - triangle below node
    fixed_nodes = np.where(bc_type == 1)[0]
    if len(fixed_nodes) > 0:
        triangle_height = symbol_size
        triangle_width = symbol_size * 0.8
        tri_patches = []
        for node_idx in fixed_nodes:
            x, y = nodes[node_idx]
            tri_patches.append(patches.Polygon([
                [x - triangle_width/2, y - triangle_height],
                [x + triangle_width/2, y - triangle_height],
                [x, y]
            ], closed=True))
        pc = PatchCollection(tri_patches, facecolor='none', edgecolor='red', linewidth=1.5)
        ax.add_collection(pc)

        legend_handles.append(
            plt.Line2D([0], [0], marker='^', color='red', linestyle='None',
                      markersize=8, label='Fixed (bc_type=1)')
        )

    # X-roller boundary conditions (type 2) - circle + line on left/right sides
    x_roller_nodes = np.where(bc_type == 2)[0]
    if saved_roller_x:
        x_roller_nodes = np.unique(np.concatenate([x_roller_nodes, np.array(sorted(saved_roller_x), dtype=int)]))
    if len(x_roller_nodes) > 0:
        circle_radius = symbol_size * 0.4
        line_length = symbol_size
        x_mid = (x_min + x_max) / 2
        circle_patches = []
        roller_line_segs = []
        for node_idx in x_roller_nodes:
            x, y = nodes[node_idx]
            is_left_side = x < x_mid
            if is_left_side:
                cx = x - circle_radius
                lx = cx - circle_radius
            else:
                cx = x + circle_radius
                lx = cx + circle_radius
            circle_patches.append(patches.Circle((cx, y), circle_radius))
            roller_line_segs.append([[lx, y - line_length/2], [lx, y + line_length/2]])

        pc = PatchCollection(circle_patches, facecolor='none', edgecolor='blue', linewidth=1)
        ax.add_collection(pc)
        lc = LineCollection(roller_line_segs, colors='blue', linewidths=1)
        ax.add_collection(lc)

        legend_handles.append(
            plt.Line2D([0], [0], marker='o', color='blue', linestyle='None',
                      markersize=6, markerfacecolor='none', markeredgewidth=1, label='Y-Roller (bc_type=3)')
        )

    # Specified force boundary conditions (type 4) - vector arrows
    force_nodes = np.where(bc_type == 4)[0]
    if len(force_nodes) > 0:
        force_magnitudes = np.array([np.sqrt(bc_values[ni][0]**2 + bc_values[ni][1]**2) for ni in force_nodes])
        max_force = force_magnitudes.max() if len(force_magnitudes) > 0 else 0
        if max_force > 0:
            scale = symbol_size * 3 / max_force
            gap = symbol_size * 0.5
            arrow_segs = []
            tip_xs, tip_ys = [], []

            for node_idx in force_nodes:
                x, y = nodes[node_idx]
                fx, fy = bc_values[node_idx]
                scaled_fx = fx * scale
                scaled_fy = fy * scale
                mag = np.sqrt(scaled_fx**2 + scaled_fy**2)
                if mag == 0:
                    continue
                ux, uy = scaled_fx / mag, scaled_fy / mag
                tail_x = x - scaled_fx
                tail_y = y - scaled_fy
                tip_x = x - ux * gap
                tip_y = y - uy * gap
                arrow_segs.append([[tail_x, tail_y], [tip_x, tip_y]])
                tip_xs.append(tip_x)
                tip_ys.append(tip_y)

            if arrow_segs:
                lc = LineCollection(arrow_segs, colors='green', linewidths=1.2)
                ax.add_collection(lc)
                ax.plot(tip_xs, tip_ys, marker='v', color='green', markersize=3, linestyle='None')

        legend_handles.append(
            plt.Line2D([0], [0], marker=r'$\rightarrow$', color='green', linestyle='None',
                      markersize=12, label='Applied Force')
        )

def plot_fem_results(fem_data, solution, plot_type=['deformation', 'shear_strain', 'displace_vector'],
                    deform_percent=15, show_mesh=True, show_reinforcement=True, figsize=(12, 8), label_elements=False,
                    plot_nodes=False, plot_elements=False, plot_boundary=True, displacement_tolerance=0.5,
                    scale_vectors=False, save_png=False, save_dxf=False, dpi=300):
    """
    Plot FEM results with various visualization options.

    Parameters:
        fem_data: FEM data dictionary from build_fem_data
        solution: FEM solution dictionary from solve_fem
        plot_type: Comma-separated plot types. Valid types:
            'deformation' - deformed mesh overlay
            'displace_mag' - displacement magnitude contours
            'displace_vector' - displacement vectors at corner nodes
            'stress' - von Mises stress contours
            'strain' - equivalent strain contours
            'shear_strain' - viscoplastic max shear strain contours
            'yield' - Mohr-Coulomb yield function contours
        deform_percent: Target deformation as percentage of mesh height (default 15).
        show_mesh: Show mesh lines
        show_reinforcement: Show reinforcement elements
        figsize: Figure size (width, height)
        label_elements: Show element ID labels at centroids
        plot_nodes: For displace_vector, show dots at node locations
        plot_elements: For displace_vector, show all element edges
        plot_boundary: For displace_vector, show boundary edges only (default)
        displacement_tolerance: Fraction of max displacement below which vectors are hidden
        scale_vectors: For displace_vector, auto-scale vectors for visibility
        save_png: Save figure to PNG file
        dpi: Resolution for saved PNG
    """
    
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))
    
    # Accept a single string or a list of strings
    if isinstance(plot_type, str):
        plot_types = [plot_type.strip().lower()]
    else:
        plot_types = [pt.strip().lower() for pt in plot_type]
    valid_types = ['displace_mag', 'displace_vector', 'deformation', 'stress', 'strain', 'shear_strain', 'yield']
    
    # Validate plot types
    for pt in plot_types:
        if pt not in valid_types:
            raise ValueError(f"Unknown plot_type: '{pt}'. Valid types: {valid_types}")
    
    # Auto-calculate deformation scale so max displacement is deform_percent of mesh height
    # Use VP displacement if available (matches what plot_deformed_mesh will plot)
    disp_elastic = solution.get("displacements_elastic", None)
    disp_for_scale = displacements - disp_elastic if disp_elastic is not None else displacements
    u_arr, v_arr = _extract_uv(disp_for_scale, fem_data)
    max_disp = np.max(np.sqrt(u_arr**2 + v_arr**2))
    mesh_height = np.max(nodes[:, 1]) - np.min(nodes[:, 1])
    if max_disp > 1e-30:
        deform_scale = max(1.0, (mesh_height * deform_percent / 100) / max_disp)
    else:
        deform_scale = 1.0
    
    # Create subplots based on number of plot types
    # When the first plot is 'deformation', add a thin extra row for its legend
    n_plots = len(plot_types)
    has_deform_legend = n_plots > 1 and plot_types[0] == 'deformation'

    # Mesh bounds (used below for axis limits and single-plot figure sizing)
    nodes = fem_data["nodes"]
    x_min, x_max = np.min(nodes[:, 0]), np.max(nodes[:, 0])
    y_min, y_max = np.min(nodes[:, 1]), np.max(nodes[:, 1])
    x_margin = (x_max - x_min) * 0.05
    y_margin = (y_max - y_min) * 0.05

    if n_plots == 1:
        # With equal aspect, a wide/short slope only fills a thin band of a tall
        # figure, leaving the colorbar towering over the actual plot. Size the
        # figure height to the data aspect ratio so the image fills the figure
        # and the colorbar matches its height.
        data_w = (x_max - x_min) + 2 * x_margin
        data_h = (y_max - y_min) + 2 * y_margin
        if data_w > 0:
            single_height = figsize[0] * (data_h / data_w)
            # Clamp to a sensible range so very flat/steep slopes stay readable
            single_height = float(np.clip(single_height, 2.0, figsize[1]))
        else:
            single_height = figsize[1]
        fig, ax = plt.subplots(figsize=(figsize[0], single_height), layout='constrained')
        axes = [ax]
        legend_ax = None
    else:
        height_factor = min(0.8, 1.2 / n_plots)
        total_height = figsize[1] * n_plots * height_factor
        if has_deform_legend:
            # Add a thin row after the first plot for the legend
            height_ratios = [1] + [0.08] + [1] * (n_plots - 1)
            fig, all_axes = plt.subplots(n_plots + 1, 1,
                                         figsize=(figsize[0], total_height),
                                         layout='constrained',
                                         gridspec_kw={'height_ratios': height_ratios})
            axes = [all_axes[0]] + list(all_axes[2:])
            legend_ax = all_axes[1]
            legend_ax.set_axis_off()
        else:
            fig, axes = plt.subplots(n_plots, 1,
                                     figsize=(figsize[0], total_height),
                                     layout='constrained')
            legend_ax = None
        if not isinstance(axes, (list, np.ndarray)):
            axes = [axes]
        elif isinstance(axes, np.ndarray):
            axes = list(axes)


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
        
        if pt == 'displace_mag':
            plot_displacement_contours(ax, fem_data, solution, show_mesh, show_reinforcement, 
                                     cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'displace_vector':
            plot_displacement_vectors(ax, fem_data, solution, show_mesh, show_reinforcement, 
                                    cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements,
                                    plot_nodes=plot_nodes, plot_elements=plot_elements, plot_boundary=plot_boundary,
                                    displacement_tolerance=displacement_tolerance, scale_vectors=scale_vectors)
        elif pt == 'deformation':
            plot_deformed_mesh(ax, fem_data, solution, deform_scale, show_mesh, show_reinforcement,
                             cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'stress':
            plot_stress_contours(ax, fem_data, solution, show_mesh, show_reinforcement,
                               cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'strain':
            plot_strain_contours(ax, fem_data, solution, show_mesh, show_reinforcement,
                               cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'shear_strain':
            plot_shear_strain_contours(ax, fem_data, solution, show_mesh, show_reinforcement,
                                     cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        elif pt == 'yield':
            plot_yield_function_contours(ax, fem_data, solution, show_mesh, show_reinforcement,
                                        cbar_shrink=cbar_shrink, cbar_labelpad=cbar_labelpad, label_elements=label_elements)
        
        # Set consistent axis limits for all plots (including single plots)
        ax.set_xlim(x_min - x_margin, x_max + x_margin)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.set_aspect('equal')

    # Place deformation legend in the dedicated legend row
    if has_deform_legend and legend_ax is not None:
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            legend_ax.legend(handles, labels, loc='center', ncol=4, fontsize=10,
                             frameon=False)

    if save_png:
        plt.savefig('fem_results.png', dpi=dpi, bbox_inches='tight')
    if save_dxf:
        from .cad import axes_to_dxf
        # One DXF per panel (each plot type), since the figure is multi-panel.
        for i, pt in enumerate(plot_types):
            if i < len(axes):
                axes_to_dxf(axes[i], f'fem_results_{pt}.dxf')
    
    plt.show()
    
    # Return appropriate values
    if n_plots == 1:
        return fig, axes[0]
    else:
        return fig, axes


def plot_displacement_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                              cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot total displacement magnitude as filled contours using the viridis colormap.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))
    
    # Calculate displacement magnitudes
    u, v = _extract_uv(displacements, fem_data)
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
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
    
    if triangles:
        triangles = np.array(triangles)
        
        # Create contour plot
        tcf = ax.tricontourf(nodes[:, 0], nodes[:, 1], triangles, disp_mag, gid='DISPLACEMENT_CONTOURS',
                         
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
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)
    
    ax.set_aspect('equal')
    ax.set_title('Displacement Magnitude Contours')
    ax.set_xlabel('x')
    ax.set_ylabel('y')


def _get_mesh_boundary(fem_data):
    """
    Compute the boundary edges of the mesh.
    
    Returns:
        boundary_edges: List of (node1, node2) tuples representing boundary edges
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    # Count how many times each edge appears
    edge_count = {}
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        
        # Define edges for each element type
        if elem_type == 3:  # Triangle
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type == 4:  # Quadrilateral
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        else:
            continue
        
        # Count each edge (both directions)
        for edge in edges:
            # Normalize edge direction (smaller node first)
            normalized_edge = tuple(sorted(edge))
            edge_count[normalized_edge] = edge_count.get(normalized_edge, 0) + 1
    
    # Boundary edges appear only once
    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]
    
    return boundary_edges


def plot_displacement_vectors(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                             cbar_shrink=0.8, cbar_labelpad=20, label_elements=False,
                             plot_nodes=False, plot_elements=False, plot_boundary=True,
                             displacement_tolerance=1e-6, scale_vectors=True):
    """
    Plot displacement vectors at corner nodes of each element.

    If viscoplastic solution data is available (displacements_elastic in solution),
    plots VP displacement (total - elastic) to show the failure mechanism rather
    than the gravity settlement. Otherwise plots total displacement.

    Parameters:
        ax: Matplotlib axes
        fem_data: FEM data dictionary
        solution: FEM solution dictionary
        show_mesh: Show mesh lines or boundary
        show_reinforcement: Show reinforcement elements
        label_elements: Show element ID labels
        plot_nodes: If True, show dots at all node locations
        plot_elements: If True, show all element edges
        plot_boundary: If True, show only boundary edges (default)
        displacement_tolerance: Fraction of max displacement below which vectors are hidden
        scale_vectors: If True (default), auto-scale vectors for visibility
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))

    # Use VP displacement (total - elastic) if available, to show failure mechanism
    # This removes the gravity settlement and shows only plastic deformation
    disp_elastic = solution.get("displacements_elastic", None)
    if disp_elastic is not None:
        disp_vp = displacements - disp_elastic
        u, v = _extract_uv(disp_vp, fem_data)
    else:
        u, v = _extract_uv(displacements, fem_data)
    disp_mag = np.sqrt(u**2 + v**2)
    max_disp_mag = np.max(disp_mag)

    if max_disp_mag < 1e-30:
        print("Warning: No VP displacements to plot")
        return

    # Collect corner nodes only (avoid mid-side nodes of quad8/tri6)
    corner_nodes = set()
    for i, elem in enumerate(elements):
        et = element_types[i]
        if et == 8 or et == 9:
            for j in range(4):
                corner_nodes.add(elem[j])
        elif et == 6:
            for j in range(3):
                corner_nodes.add(elem[j])
        else:
            for j in range(et):
                corner_nodes.add(elem[j])

    # Absolute threshold
    abs_tol = displacement_tolerance * max_disp_mag

    # Plot boundary outline first
    if show_mesh:
        if plot_elements:
            plot_mesh_lines(ax, fem_data, color='lightgray', alpha=0.5, linewidth=0.5)
        elif plot_boundary:
            boundary_edges = _get_mesh_boundary(fem_data)
            for edge in boundary_edges:
                x_coords = [nodes[edge[0], 0], nodes[edge[1], 0]]
                y_coords = [nodes[edge[0], 1], nodes[edge[1], 1]]
                ax.plot(x_coords, y_coords, 'k-', alpha=0.7, linewidth=1.0)

    # Plot small vectors at corner nodes
    corner_list = sorted(corner_nodes)
    cx = nodes[corner_list, 0]
    cy = nodes[corner_list, 1]
    cu = u[corner_list]
    cv = v[corner_list]
    cmag = disp_mag[corner_list]

    mask = cmag > abs_tol

    if np.sum(mask) == 0:
        print("Warning: All displacements below tolerance")
        return

    _q = ax.quiver(cx[mask], cy[mask], cu[mask], cv[mask], gid='DISPLACE_VECTORS',
              angles='xy', color='black', alpha=0.7,
              scale=None, width=0.002, headwidth=3, headlength=4,
              headaxislength=3, pivot='tail')

    # Plot node dots if requested
    if plot_nodes:
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=1, alpha=0.4, gid='MESH_NODES')

    # Plot reinforcement
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution)

    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)

    # Dummy colorbar for axis alignment with other subplots
    dummy_data = np.array([[0, 1]])
    dummy_im = ax.imshow(dummy_data, cmap='viridis', alpha=0)
    cbar = plt.colorbar(dummy_im, ax=ax, shrink=cbar_shrink)
    cbar.set_label('', color='white')
    cbar.set_ticks([])
    cbar.set_ticklabels([])
    cbar.outline.set_color('white')
    cbar.outline.set_linewidth(0)

    F = solution.get("F", None)
    title = 'Viscoplastic Displacement Vectors' if disp_elastic is not None else 'Displacement Vectors'
    if F is not None:
        title += f'  F={F:.2f}'
    ax.set_title(title, fontsize=12, pad=15)


def plot_stress_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                        cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot von Mises stress contours with yielding elements highlighted.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    stresses = solution.get("stresses", np.zeros((len(elements), 4)))
    
    # Use yield function to determine plastic elements for consistency
    # If yield_function is available, use it; otherwise fall back to plastic_elements
    yield_function = solution.get("yield_function", None)
    if yield_function is not None:
        plastic_elements = yield_function > 0  # F > 0 means yielding
    else:
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
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            coords = nodes[elem[:3]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            coords = nodes[elem[:4]]
            patch = Polygon(coords, closed=True)
            patches_list.append(patch)
            stress_values.append(von_mises[i])
    
    if patches_list:
        from matplotlib.collections import PatchCollection
        
        # Create patch collection
        p = PatchCollection(patches_list, alpha=0.8, edgecolors='none', gid='STRESS_CONTOURS')
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
                elif elem_type == 6:  # 6-node triangle - use corner nodes
                    coords = nodes[elem[:3]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
                elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                    coords = nodes[elem[:4]]
                    coords = np.vstack([coords, coords[0]])  # Close the polygon
                    ax.plot(coords[:, 0], coords[:, 1], 'r-', linewidth=2, alpha=0.8)
    
    # Plot mesh
    if show_mesh:
        plot_mesh_lines(ax, fem_data, color='gray', alpha=0.3, linewidth=0.3)
    
    # Plot reinforcement with force visualization
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_forces(ax, fem_data, solution)
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)
    
    ax.set_aspect('equal')
    ax.set_title('von Mises Stress (Red outline = Yielding/Plastic Elements)')
    ax.set_xlabel('x')
    ax.set_ylabel('y')


def plot_deformed_mesh(ax, fem_data, solution, deform_scale=1.0, show_mesh=True, show_reinforcement=True, 
                       cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot deformed mesh overlay on original mesh.

    Shows the original mesh in light gray and the deformed mesh in blue.
    Uses VP displacement (total - elastic) when available to show the failure
    mechanism rather than gravity settlement. The deform_scale parameter
    amplifies the displacements for visibility.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))

    # Use VP displacement (total - elastic) if available, to show failure mechanism
    disp_elastic = solution.get("displacements_elastic", None)
    if disp_elastic is not None:
        disp = displacements - disp_elastic
    else:
        disp = displacements

    # Calculate deformed node positions
    u, v = _extract_uv(disp, fem_data)
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
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data_deformed)  # Label on deformed mesh
    
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
    F = solution.get("F", None)
    disp_label = 'Viscoplastic Deformation' if disp_elastic is not None else 'Mesh Deformation'
    scale_str = f'{deform_scale:.0f}' if deform_scale >= 10 else f'{deform_scale:.1f}'
    title = f'{disp_label} (Scale = {scale_str}x)'
    if F is not None:
        title += f'  F={F:.2f}'
    ax.set_title(title, fontsize=12, pad=15)
    ax.set_xlabel('x')
    ax.set_ylabel('y')


def _add_element_labels(ax, fem_data):
    """
    Add element ID labels at element centers.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        
        # Get element nodes for centroid calculation
        if elem_type == 3:  # Triangle
            elem_nodes = nodes[elem[:3]]
        elif elem_type == 4:  # Quad
            elem_nodes = nodes[elem[:4]]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            elem_nodes = nodes[elem[:3]]
        elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
            elem_nodes = nodes[elem[:4]]
        else:
            continue
            
        # Calculate centroid
        centroid = np.mean(elem_nodes, axis=0)
        
        # Add label (1-based indexing for display)
        ax.text(centroid[0], centroid[1], str(i+1),
                ha='center', va='center', fontsize=6, 
                color='darkblue', alpha=0.7, zorder=100)


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
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[0])]
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            edges = [(elem[0], elem[1]), (elem[1], elem[2]), (elem[2], elem[3]), (elem[3], elem[0])]
        else:
            continue
        
        for edge in edges:
            line_coords = nodes[[edge[0], edge[1]]]
            lines.append(line_coords)
    
    if lines:
        lc = LineCollection(lines, colors=color, alpha=alpha, linewidths=linewidth, label=label, gid='MESH')
        ax.add_collection(lc)


def plot_reinforcement_lines(ax, fem_data, solution, color='red', alpha=1.0, linewidth=2, label=None):
    """
    Plot reinforcement and pile elements as lines with distinct colors.
    """
    if 'elements_1d' not in fem_data:
        return

    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    element_types_1d = fem_data["element_types_1d"]
    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(elements_1d), dtype=bool))

    reinf_lines = []
    pile_lines = []
    for i, elem in enumerate(elements_1d):
        elem_type = element_types_1d[i]
        if elem_type >= 2:
            line_coords = nodes[elem[:2]]
            if pile_elem_mask[i]:
                pile_lines.append(line_coords)
            else:
                reinf_lines.append(line_coords)

    if reinf_lines:
        lc = LineCollection(reinf_lines, colors=color, alpha=alpha, linewidths=linewidth, label=label, gid='REINFORCEMENT')
        ax.add_collection(lc)
    if pile_lines:
        pile_label = label.replace('Reinforcement', 'Pile') if label and 'Reinforcement' in label else None
        lc = LineCollection(pile_lines, colors='green', alpha=alpha, linewidths=linewidth + 1, label=pile_label, gid='PILES')
        ax.add_collection(lc)


def plot_reinforcement_forces(ax, fem_data, solution):
    """
    Plot reinforcement elements colored by force level.

    Color scheme:
    - Blue to green to yellow to red: 0 to Tmax (tension force ramp)
    - Magenta: element has yielded and is at residual capacity Tres
    - White/open with dashed outline: element has pulled out (broken, T=0)
    - Gray: element carrying no tension (inactive or in compression)
    """
    if 'elements_1d' not in fem_data:
        return

    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.cm as cm

    nodes = fem_data["nodes"]
    elements_1d = fem_data["elements_1d"]
    forces_1d = solution.get("forces_1d", np.zeros(len(elements_1d)))
    t_allow = fem_data.get("t_allow_by_1d_elem", np.ones(len(elements_1d)))
    t_res = fem_data.get("t_res_by_1d_elem", np.zeros(len(elements_1d)))
    failed_1d = solution.get("failed_1d_elements", np.zeros(len(elements_1d), dtype=bool))

    # Find global Tmax (max of all t_allow values)
    t_max_global = t_allow.max() if len(t_allow) > 0 else 1.0

    # Custom colormap: blue -> white -> red (coolwarm style)
    force_cmap = LinearSegmentedColormap.from_list(
        'force_ramp', ['#2166ac', '#f7f7f7', '#d73027'], N=256)

    # Classify and draw each element
    normal_lines = []
    normal_colors = []
    tres_lines = []
    pullout_lines = []
    inactive_lines = []

    pile_elem_mask = fem_data.get("pile_elem_mask", np.zeros(len(elements_1d), dtype=bool))
    pile_force_lines = []
    pile_force_colors = []
    forces_pile_lateral = solution.get("forces_pile_lateral", np.array([]))

    # Build pile element index mapping: global 1d index -> pile force index
    pile_force_idx = 0

    for i in range(len(elements_1d)):
        elem = elements_1d[i]
        coords = nodes[elem[:2]]

        if pile_elem_mask[i]:
            # Pile element — color by lateral (shear) force
            if pile_force_idx < len(forces_pile_lateral):
                pile_force_lines.append(coords)
                pile_force_colors.append(abs(forces_pile_lateral[pile_force_idx]))
            pile_force_idx += 1
            continue

        force = forces_1d[i]
        is_failed = failed_1d[i]

        if is_failed and t_res[i] < 1e-6 and force < 1e-6:
            pullout_lines.append(coords)
        elif is_failed and t_res[i] > 1e-6:
            tres_lines.append(coords)
        elif force > 1e-6:
            ratio = min(force / t_max_global, 1.0) if t_max_global > 0 else 0.0
            normal_lines.append(coords)
            normal_colors.append(force_cmap(ratio))
        else:
            inactive_lines.append(coords)

    # Draw inactive elements (cyan, solid)
    if inactive_lines:
        lc_outline = LineCollection(inactive_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=3.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(inactive_lines, colors='#00CC00', linewidths=3, alpha=0.9, zorder=4)
        ax.add_collection(lc)
        ax.plot([], [], '-', color='#00CC00', linewidth=3, alpha=0.9, label='Inactive (no tension)')

    # Draw normal tension elements (force-colored)
    if normal_lines:
        lc_outline = LineCollection(normal_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=4.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(normal_lines, colors=normal_colors, linewidths=3, alpha=0.9, zorder=5)
        ax.add_collection(lc)

        # Add colorbar
        sm = cm.ScalarMappable(cmap=force_cmap, norm=plt.Normalize(0, t_max_global))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label('Reinforcement Force', rotation=270, labelpad=15, fontsize=10)

    # Draw elements at Tres (magenta)
    if tres_lines:
        lc_outline = LineCollection(tres_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=5.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(tres_lines, colors='magenta', linewidths=3, alpha=0.9, zorder=6)
        ax.add_collection(lc)
        ax.plot([], [], '-', color='magenta', linewidth=3, label='At residual (Tres)')

    # Draw pulled-out elements (orange, solid)
    if pullout_lines:
        lc_outline = LineCollection(pullout_lines, colors='black', linewidths=4.5, alpha=0.9, zorder=5.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(pullout_lines, colors='black', linewidths=3, alpha=0.9, zorder=6)
        ax.add_collection(lc)
        ax.plot([], [], '-', color='black', linewidth=3, alpha=0.9, label='Pulled out')

    # Draw pile elements colored by lateral (shear) force
    if pile_force_lines:
        from matplotlib.colors import Normalize
        max_lateral = max(pile_force_colors) if pile_force_colors else 1.0
        pile_cmap = plt.cm.Greens
        pile_norm = Normalize(vmin=0, vmax=max_lateral if max_lateral > 0 else 1.0)
        colors = [pile_cmap(pile_norm(v)) for v in pile_force_colors]
        lc_outline = LineCollection(pile_force_lines, colors='black', linewidths=5, alpha=0.9, zorder=5.9)
        ax.add_collection(lc_outline)
        lc = LineCollection(pile_force_lines, colors=colors, linewidths=3.5, alpha=0.9, zorder=6)
        ax.add_collection(lc)
        sm = cm.ScalarMappable(cmap=pile_cmap, norm=pile_norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label('Pile Shear Force', rotation=270, labelpad=15, fontsize=10)

    # Add legend if any special states exist
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(loc='lower right', fontsize=9, framealpha=0.9)


def plot_reinforcement_force_profiles(fem_data, solution, figsize=(12, 8), save_png=False, dpi=300):
    """
    Plot axial force profiles along each reinforcement line as subplots.
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
    
    if save_png:
        filename = 'plot_reinforcement_force_profiles.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    return fig, axes


def plot_ssrm_convergence(ssrm_solution, figsize=(10, 6), save_png=False, dpi=300):
    """
    Plot SSRM bisection convergence history showing F vs iteration and convergence status.
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
    
    if save_png:
        filename = 'plot_ssrm_convergence.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    return fig, (ax1, ax2)


def plot_strain_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                        cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot von Mises equivalent strain contours computed from total strains.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    strains = solution.get("strains", np.zeros((len(elements), 4)))
    
    if strains.shape[1] < 3:
        print("Warning: Strain data not available or incomplete")
        return
    
    # Calculate equivalent strain (von Mises equivalent strain)
    # For plane strain: equiv_strain = sqrt(2/3) * sqrt(eps_x^2 + eps_y^2 + eps_x*eps_y + 3/4*gamma_xy^2)
    eps_x = strains[:, 0]
    eps_y = strains[:, 1]
    gamma_xy = strains[:, 2]
    
    equiv_strain = np.sqrt((2/3) * (eps_x**2 + eps_y**2 + eps_x*eps_y + 0.75*gamma_xy**2))
    
    # Plot contours
    _plot_element_contours(ax, fem_data, equiv_strain, 'Equivalent Strain', 
                          show_mesh, show_reinforcement, cbar_shrink, cbar_labelpad, label_elements)


def plot_shear_strain_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True,
                              cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot viscoplastic max shear strain contours.

    Uses accumulated viscoplastic strains from the solution (vp_shear_strain key).
    Falls back to total shear strain if VP data is not available.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]

    vp_shear_strain = solution.get("vp_shear_strain", None)
    if vp_shear_strain is None:
        # Fallback to total shear strain if VP not available
        strains = solution.get("strains", np.zeros((len(elements), 4)))
        if strains.shape[1] >= 4:
            vp_shear_strain = strains[:, 3]
        else:
            print("Warning: Shear strain data not available")
            return

    _plot_nodal_contours(ax, fem_data, vp_shear_strain, 'VP Max Shear Strain',
                        False, False, cbar_shrink, cbar_labelpad,
                        colormap='coolwarm', label_elements=label_elements)

    # Draw reinforcement with force-based coloring
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_forces(ax, fem_data, solution)

    F = solution.get("F", None)
    title = 'Viscoplastic Shear Strain'
    if F is not None:
        title += f'  F={F:.2f}'
    ax.set_title(title, fontsize=12, pad=15)


def plot_yield_function_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True, 
                                cbar_shrink=0.8, cbar_labelpad=20, label_elements=False):
    """
    Plot yield function values (Mohr-Coulomb failure criterion).
    Positive values indicate yielding/failure, negative values indicate elastic state.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    yield_function = solution.get("yield_function", None)
    
    if yield_function is None:
        print("Warning: Yield function data not available in solution")
        # Create dummy data
        yield_function = np.zeros(len(elements))
    
    # Create custom colormap for yield function visualization
    # Strong blue for very negative (very safe), white near zero, red for positive (yielding)
    from matplotlib.colors import LinearSegmentedColormap
    
    # Define color transitions for yield function
    # F < 0: shades of blue/green (elastic/safe)
    # F = 0: white/light gray (critical)
    # F > 0: shades of red (yielding/plastic)
    colors_below = ['#0000FF', '#0066FF', '#00AAFF', '#00DDDD', '#CCCCCC']  # Blue to gray
    colors_above = ['#FFCCCC', '#FF9999', '#FF6666', '#FF3333', '#FF0000', '#CC0000']  # Light red to dark red
    
    # Create custom colormap with sharp transition at F=0
    n_bins = 256
    n_below = int(n_bins * 0.7)  # 70% for negative values
    n_above = n_bins - n_below   # 30% for positive values
    
    from matplotlib.colors import ListedColormap
    colors_below_interp = plt.cm.Blues_r(np.linspace(0.2, 0.9, n_below))
    colors_above_interp = plt.cm.Reds(np.linspace(0.3, 1.0, n_above))
    colors_all = np.vstack([colors_below_interp, colors_above_interp])
    cmap_yield = ListedColormap(colors_all)
    
    # Set visualization bounds - asymmetric to focus on near-yield region
    vmin = -200  # Cap negative values for better contrast
    vmax = 50    # Positive values are more important
    
    # Plot each element as a colored patch
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon
    patches_list = []
    values_list = []
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            coords = nodes[elem[:3]]
        elif elem_type == 4:  # Quad
            coords = nodes[elem[:4]]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            coords = nodes[elem[:3]]
        elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
            coords = nodes[elem[:4]]
        else:
            continue
            
        patch = Polygon(coords, closed=True)
        patches_list.append(patch)
        # Clip values for visualization
        values_list.append(np.clip(yield_function[i], vmin, vmax))
    
    if patches_list:
        p = PatchCollection(patches_list, alpha=0.9, edgecolors='gray', linewidths=0.3)
        p.set_array(np.array(values_list))
        p.set_cmap(cmap_yield)
        p.set_clim(vmin, vmax)
        ax.add_collection(p)
        
        # Add colorbar with custom ticks
        cbar = plt.colorbar(p, ax=ax, shrink=cbar_shrink)
        cbar.set_label('Yield Function F', rotation=270, labelpad=cbar_labelpad)
        
        # Set custom ticks to highlight key values
        tick_values = [-200, -100, -50, -20, -10, -5, 0, 5, 10, 20, 50]
        tick_labels = ['-200', '-100', '-50', '-20', '-10', '-5', '0', '5', '10', '20', '50']
        # Filter ticks to those within bounds
        valid_ticks = [(v, l) for v, l in zip(tick_values, tick_labels) if vmin <= v <= vmax]
        if valid_ticks:
            tick_values, tick_labels = zip(*valid_ticks)
            cbar.set_ticks(tick_values)
            cbar.set_ticklabels(tick_labels)
        
        # Add a line at F=0
        cbar.ax.axhline(y=0, color='black', linewidth=2)
    
    # Add yield function values as text on elements (if requested or for yielding elements)
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        
        # Get element centroid
        if elem_type == 3:  # Triangle
            elem_nodes = nodes[elem[:3]]
        elif elem_type == 4:  # Quad
            elem_nodes = nodes[elem[:4]]
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            elem_nodes = nodes[elem[:3]]
        elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
            elem_nodes = nodes[elem[:4]]
        else:
            continue
            
        centroid = np.mean(elem_nodes, axis=0)
        
        # Show values for elements that are close to yielding or already yielding
        # or if label_elements is True
        f_val = yield_function[i]
        
        if label_elements or f_val > -50:  # Show if requested or if close to yielding
            # Format the number based on magnitude
            if abs(f_val) < 10:
                text = f'{f_val:.1f}'
            else:
                text = f'{f_val:.0f}'
            
            # Choose text color based on value
            if f_val > 0:
                color = 'white'  # White on red background
                fontweight = 'bold'
            elif f_val > -10:
                color = 'black'  # Black on light background
                fontweight = 'normal'
            else:
                color = 'white'  # White on blue background
                fontweight = 'normal'
            
            # Only show for elements near yield or if explicitly requested
            if label_elements or f_val > -30:
                ax.text(centroid[0], centroid[1], text,
                       ha='center', va='center', fontsize=5,
                       color=color, fontweight=fontweight, alpha=0.8)
    
    # Highlight yielding elements with thick red border
    for i, elem in enumerate(elements):
        if yield_function[i] > 0:
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
            elif elem_type in [8, 9]:  # 8 or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
            else:
                continue
            
            # Close the polygon
            coords = np.vstack([coords, coords[0]])
            ax.plot(coords[:, 0], coords[:, 1], 'k-', linewidth=2.5, alpha=1.0)  # Black border for yielding elements
    
    # Add reinforcement if requested
    if show_reinforcement and 'elements_1d' in fem_data:
        plot_reinforcement_lines(ax, fem_data, solution)
    
    # Add title indicating yield state
    ax.set_title('Yield Function (Red: F>0 Yielding/Plastic, Blue: F<0 Elastic)', fontsize=12, pad=15)
    
    # Add statistics to the plot
    n_yielding = np.sum(yield_function > 0)
    n_total = len(yield_function)
    n_critical = np.sum((yield_function > -10) & (yield_function <= 0))  # Near yielding
    
    stats_text = f'Yielding: {n_yielding}/{n_total} elements\n'
    stats_text += f'Critical (F>-10): {n_critical} elements\n'
    stats_text += f'Max F: {np.max(yield_function):.1f}\n'
    stats_text += f'Min F: {np.min(yield_function):.1f}'
    
    ax.text(0.02, 0.98, stats_text,
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))


def _plot_element_contours(ax, fem_data, values, label, show_mesh=True, show_reinforcement=True,
                          cbar_shrink=0.8, cbar_labelpad=20, label_elements=False, colormap='viridis'):
    """
    Plot element-based scalar data as colored patches (one color per element).
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    # For element-based values, we need to interpolate to nodes or use a different approach
    # Let's use a simpler approach: plot each element as a colored patch
    
    # Create contour plot by directly coloring elements
    if np.max(values) > np.min(values):  # Only plot if there's variation
        # Normalize values for colormap
        vmin, vmax = np.min(values), np.max(values)
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.get_cmap(colormap)
        
        # Plot each element as colored patch
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            color = cmap(norm(values[i]))
            
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor=color, edgecolor='none', alpha=0.8)
                ax.add_patch(quad)
        
        # Create colorbar using a ScalarMappable
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=cbar_shrink, pad=0.05)
        cbar.set_label(label, rotation=270, labelpad=cbar_labelpad)
    else:
        # Uniform values - just color all elements the same
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, facecolor='lightblue', edgecolor='none', alpha=0.7)
                ax.add_patch(quad)
    
    # Overlay mesh if requested
    if show_mesh:
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(quad)
    
    # Add reinforcement if requested
    if show_reinforcement:
        elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
        if len(elements_1d) > 0:
            for elem in elements_1d:
                if len(elem) >= 2:
                    x_coords = [nodes[elem[0], 0], nodes[elem[1], 0]]
                    y_coords = [nodes[elem[0], 1], nodes[elem[1], 1]]
                    ax.plot(x_coords, y_coords, 'r-', linewidth=2, alpha=0.8)
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)
    
    ax.set_aspect('equal')


def _plot_nodal_contours(ax, fem_data, element_values, label, show_mesh=True, show_reinforcement=True,
                        cbar_shrink=0.8, cbar_labelpad=20, colormap='viridis', label_elements=False):
    """
    Plot smooth filled contours by averaging element values to nodes and triangulating.
    """
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    
    # Interpolate element values to nodes
    nodal_values = np.zeros(len(nodes))
    node_counts = np.zeros(len(nodes))  # For averaging
    
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        elem_nodes = elem[:elem_type] if elem_type <= len(elem) else elem
        
        # Add this element's value to all its nodes
        for node_id in elem_nodes:
            if node_id < len(nodes):
                nodal_values[node_id] += element_values[i]
                node_counts[node_id] += 1
    
    # Average values at nodes (avoid division by zero)
    valid_nodes = node_counts > 0
    nodal_values[valid_nodes] /= node_counts[valid_nodes]
    
    # Create triangulation for smooth contouring
    triangles = []
    for i, elem in enumerate(elements):
        elem_type = element_types[i]
        if elem_type == 3:  # Triangle
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type == 4:  # Quad - split into triangles
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
        elif elem_type == 6:  # 6-node triangle - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
        elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
            triangles.append([elem[0], elem[1], elem[2]])
            triangles.append([elem[0], elem[2], elem[3]])
    
    if not triangles:
        print("No valid elements for contouring")
        return
    
    import matplotlib.tri as tri
    triangles = np.array(triangles)
    
    # Create triangulation
    triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], triangles)
    
    # Create smooth contour plot
    if np.max(nodal_values) > np.min(nodal_values):  # Only plot if there's variation
        levels = np.linspace(np.min(nodal_values), np.max(nodal_values), 20)
        cs = ax.tricontourf(triang, nodal_values, levels=levels, cmap=colormap)
        # DXF layer named after the plotted quantity (e.g. "VP Max Shear Strain").
        cs.set_gid((label or 'CONTOURS').upper().replace(' ', '_') + '_CONTOURS')

        # Add colorbar
        cbar = plt.colorbar(cs, ax=ax, shrink=cbar_shrink, pad=0.02)
        cbar.set_label(label, rotation=270, labelpad=20)
    else:
        # Uniform values - just color all elements the same
        uniform_color = plt.get_cmap(colormap)(0.5)
        for triangle_nodes in triangles:
            coords = nodes[triangle_nodes]
            triangle = plt.Polygon(coords, facecolor=uniform_color, edgecolor='none', alpha=0.8)
            ax.add_patch(triangle)
    
    # Overlay mesh if requested
    if show_mesh:
        for i, elem in enumerate(elements):
            elem_type = element_types[i]
            if elem_type == 3:  # Triangle
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type == 4:  # Quad
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(quad)
            elif elem_type == 6:  # 6-node triangle - use corner nodes
                coords = nodes[elem[:3]]
                triangle = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(triangle)
            elif elem_type in [8, 9]:  # 8-node or 9-node quad - use corner nodes
                coords = nodes[elem[:4]]
                quad = plt.Polygon(coords, fill=False, edgecolor='black', linewidth=0.5, alpha=0.7)
                ax.add_patch(quad)
    
    # Add reinforcement if requested
    if show_reinforcement:
        elements_1d = fem_data.get("elements_1d", np.array([]).reshape(0, 3))
        if len(elements_1d) > 0:
            for elem in elements_1d:
                if len(elem) >= 2:
                    x_coords = [nodes[elem[0], 0], nodes[elem[1], 0]]
                    y_coords = [nodes[elem[0], 1], nodes[elem[1], 1]]
                    ax.plot(x_coords, y_coords, 'r-', linewidth=2, alpha=0.8)
    
    # Add element labels if requested
    if label_elements:
        _add_element_labels(ax, fem_data)
    
    ax.set_aspect('equal')