import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.ticker import MaxNLocator
import numpy as np


def plot_seep_mesh(seep_data, show_nodes=False, show_bc=False):
    """
    Plots a mesh colored by material zone.
    
    Args:
        seep_data: Dictionary containing seepage data from import_seep2d
        show_nodes: If True, plot node points
        show_bc: If True, plot boundary condition nodes
    """
    import matplotlib.pyplot as plt

    # Extract data from seep_data
    coords = seep_data["coords"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    nbc = seep_data["nbc"]

    fig, ax = plt.subplots(figsize=(12, 5))
    materials = np.unique(element_materials)
    cmap = plt.get_cmap("tab10", len(materials))
    mat_to_color = {mat: cmap(i) for i, mat in enumerate(materials)}

    for idx, tri_nodes in enumerate(elements):
        polygon = coords[tri_nodes]
        color = mat_to_color[element_materials[idx]]
        ax.fill(*zip(*polygon), edgecolor='k', facecolor=color, linewidth=0.5)

    if show_nodes:
        ax.plot(coords[:, 0], coords[:, 1], 'k.', markersize=2)

    legend_handles = [
        plt.Line2D([0], [0], color=cmap(i), lw=4, label=f"Material {mat}")
        for i, mat in enumerate(materials)
    ]

    if show_bc:
        bc1 = coords[nbc == 1]
        bc2 = coords[nbc == 2]
        if len(bc1) > 0:
            h1, = ax.plot(bc1[:, 0], bc1[:, 1], 'ro', label="Fixed Head (nbc=1)")
            legend_handles.append(h1)
        if len(bc2) > 0:
            h2, = ax.plot(bc2[:, 0], bc2[:, 1], 'bs', label="Exit Face (nbc=2)")
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
    ax.set_title("SEEP2D Mesh with Material Zones")
    # plt.subplots_adjust(bottom=0.2)  # Add vertical cushion
    plt.tight_layout()
    plt.show()


def plot_seep_solution(seep_data, solution, levels=20, base_mat=None, fill_contours=True, phreatic=True):
    """
    Plots head contours and optionally overlays flowlines (phi) based on flow function.
    Fixed version that properly handles mesh aspect ratio and doesn't clip the plot.

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
    import numpy as np

    # Extract data from seep_data and solution
    coords = seep_data["coords"]
    elements = seep_data["elements"]
    element_materials = seep_data["element_materials"]
    k1_by_mat = seep_data["k1_by_mat"]
    head = solution["head"]
    phi = solution.get("phi")
    flowrate = solution.get("flowrate")

    # Calculate proper figure size based on mesh aspect ratio
    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()
    y_min, y_max = coords[:, 1].min(), coords[:, 1].max()

    x_range = x_max - x_min
    y_range = y_max - y_min
    mesh_aspect = x_range / y_range if y_range > 0 else 1.0

    # Set figure size to accommodate the mesh properly
    if mesh_aspect > 2.0:  # Wide mesh
        fig_width = 12
        fig_height = 12 / mesh_aspect
    elif mesh_aspect < 0.5:  # Tall mesh
        fig_height = 10
        fig_width = 10 * mesh_aspect
    else:  # Roughly square mesh
        fig_width = 10
        fig_height = 10 / mesh_aspect

    # Ensure minimum size
    fig_width = max(fig_width, 6)
    fig_height = max(fig_height, 4)

    print(f"Mesh bounds: x=[{x_min:.1f}, {x_max:.1f}], y=[{y_min:.1f}, {y_max:.1f}]")
    print(f"Mesh aspect ratio: {mesh_aspect:.2f}, Figure size: {fig_width:.1f} x {fig_height:.1f}")

    triang = tri.Triangulation(coords[:, 0], coords[:, 1], elements)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Plot material zones first (if element_materials provided)
    if element_materials is not None:
        materials = np.unique(element_materials)
        cmap = plt.get_cmap("tab10", len(materials))
        mat_to_color = {mat: cmap(i) for i, mat in enumerate(materials)}

        for idx, tri_nodes in enumerate(elements):
            polygon = coords[tri_nodes]
            color = mat_to_color[element_materials[idx]]
            ax.fill(*zip(*polygon), edgecolor='none', facecolor=color, alpha=0.5)

    vmin = np.min(head)
    vmax = np.max(head)
    hdrop = vmax - vmin
    contour_levels = np.linspace(vmin, vmax, levels)

    # Filled contours (only if fill_contours=True)
    if fill_contours:
        contourf = ax.tricontourf(triang, head, levels=contour_levels, cmap="Spectral_r", vmin=vmin, vmax=vmax,
                                  alpha=0.5)
        cbar = plt.colorbar(contourf, ax=ax, label="Total Head")
        cbar.locator = MaxNLocator(nbins=10, steps=[1, 2, 5])
        cbar.update_ticks()

    # Solid lines for head contours
    ax.tricontour(triang, head, levels=contour_levels, colors="k", linewidths=0.5)

    # Phreatic surface (pressure head = 0)
    if phreatic:
        elevation = coords[:, 1]  # y-coordinate is elevation
        pressure_head = head - elevation
        ax.tricontour(triang, pressure_head, levels=[0], colors="red", linewidths=2.0)

    # Overlay flowlines if phi is available
    if phi is not None and flowrate is not None and base_mat is not None and k1_by_mat is not None:
        # Materials are 1-based, so adjust index
        base_k = k1_by_mat[base_mat - 1]
        ne = levels - 1
        nf = (flowrate * ne) / (base_k * hdrop)
        phi_levels = round(nf) + 1
        print(f"Computed nf: {nf:.2f}, using {phi_levels} φ contours (base k={base_k}, head drop={hdrop:.3f})")
        phi_contours = np.linspace(np.min(phi), np.max(phi), phi_levels)
        ax.tricontour(triang, phi, levels=phi_contours, colors="blue", linewidths=0.7, linestyles="solid")

    # Plot the mesh boundary
    try:
        boundary = get_ordered_mesh_boundary(coords, elements)
        ax.plot(boundary[:, 0], boundary[:, 1], color="black", linewidth=1.0, label="Mesh Boundary")
    except Exception as e:
        print(f"Warning: Could not plot mesh boundary: {e}")

    # Set explicit axis limits to ensure full mesh is shown
    margin = 0.10  # 10% margin
    x_margin = x_range * margin
    y_margin = y_range * margin

    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)

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

    # Adjust layout to prevent clipping
    plt.tight_layout()
    plt.show()

def get_ordered_mesh_boundary(coords, elements):
    """
    Extracts the outer boundary of the mesh and returns it as an ordered array of points.

    Returns:
        np.ndarray of shape (N, 2): boundary coordinates in order (closed loop)
    """
    import numpy as np
    from collections import defaultdict, deque

    # Step 1: Count all edges
    edge_count = defaultdict(int)
    edge_to_nodes = {}

    for tri in elements:
        for i in range(3):
            a, b = sorted((tri[i], tri[(i + 1) % 3]))
            edge_count[(a, b)] += 1
            edge_to_nodes[(a, b)] = (tri[i], tri[(i + 1) % 3])  # preserve direction

    # Step 2: Keep only boundary edges (appear once)
    boundary_edges = [edge_to_nodes[e] for e, count in edge_count.items() if count == 1]

    if not boundary_edges:
        raise ValueError("No boundary edges found.")

    # Step 3: Build adjacency for boundary walk
    adj = defaultdict(list)
    for a, b in boundary_edges:
        adj[a].append(b)
        adj[b].append(a)

    # Step 4: Walk the boundary in order
    start = boundary_edges[0][0]
    boundary_loop = [start]
    visited = set([start])
    current = start

    while True:
        neighbors = [n for n in adj[current] if n not in visited]
        if not neighbors:
            break
        next_node = neighbors[0]
        boundary_loop.append(next_node)
        visited.add(next_node)
        current = next_node
        if current == start:
            break

    return coords[boundary_loop]