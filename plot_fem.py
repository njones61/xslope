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


def plot_fem_results(fem_data, solution, plot_type='displacement', deform_scale=1.0, 
                    show_mesh=True, show_reinforcement=True, figsize=(12, 8)):
    """
    Plot FEM results with various visualization options.
    
    Parameters:
        fem_data (dict): FEM data dictionary
        solution (dict): FEM solution dictionary
        plot_type (str): Type of plot - 'displacement', 'stress', or 'deformed_mesh'
        deform_scale (float): Scale factor for deformed mesh visualization
        show_mesh (bool): Whether to show mesh lines
        show_reinforcement (bool): Whether to show reinforcement elements
        figsize (tuple): Figure size
    
    Returns:
        matplotlib figure and axes
    """
    
    nodes = fem_data["nodes"]
    elements = fem_data["elements"]
    element_types = fem_data["element_types"]
    displacements = solution.get("displacements", np.zeros(2 * len(nodes)))
    
    if plot_type in ['displacement', 'deformed_mesh']:
        fig, ax = plt.subplots(figsize=figsize)
        
        if plot_type == 'displacement':
            plot_displacement_contours(ax, fem_data, solution, show_mesh, show_reinforcement)
        elif plot_type == 'deformed_mesh':
            plot_deformed_mesh(ax, fem_data, solution, deform_scale, show_mesh, show_reinforcement)
            
    elif plot_type == 'stress':
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(figsize[0], figsize[1]*1.2))
        plot_stress_contours(ax1, fem_data, solution, show_mesh, show_reinforcement)
        plot_deformed_mesh(ax2, fem_data, solution, deform_scale, show_mesh, show_reinforcement)
        ax2.set_title('Deformed Mesh')
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}")
    
    plt.tight_layout()
    return fig, ax if plot_type != 'stress' else (ax1, ax2)


def plot_displacement_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True):
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
        cbar = plt.colorbar(tcf, ax=ax, shrink=0.8)
        cbar.set_label('Displacement Magnitude', rotation=270, labelpad=20)
    
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


def plot_stress_contours(ax, fem_data, solution, show_mesh=True, show_reinforcement=True):
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
        cbar = plt.colorbar(p, ax=ax, shrink=0.8)
        cbar.set_label('von Mises Stress', rotation=270, labelpad=20)
    
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


def plot_deformed_mesh(ax, fem_data, solution, deform_scale=1.0, show_mesh=True, show_reinforcement=True):
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
    
    ax.set_aspect('equal')
    ax.set_title(f'Mesh Deformation (Scale Factor = {deform_scale:.1f})')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    if show_mesh or show_reinforcement:
        ax.legend()


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