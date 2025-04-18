import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_slices(profile_lines, df, piezo_line=None, failure_surface=None, fs=None, dloads=None, max_depth=None):
    """
    Plots a slopetools cross-section with profile lines, slices, piezometric line, failure surface,
    distributed loads, and optional max depth and factor of safety.

    Parameters:
        profile_lines (list): A list of profile segments, each a list of (x, y) tuples.
        df (pd.DataFrame): DataFrame containing slice geometry with required columns:
                           ['x_l', 'x_r', 'y_lb', 'y_lt', 'y_rt', 'y_rb'].
        piezo_line (list, optional): A list of (x, y) tuples representing the piezometric line.
        failure_surface (shapely.geometry.LineString, optional): Failure surface to plot, assumed to have `.coords`.
        fs (float, optional): Factor of Safety to display in the plot title.
        dloads (list, optional): List of distributed load blocks, each a list of dictionaries with 'X', 'Y', and 'Normal'.
        max_depth (float, optional): Horizontal line to indicate the analysis depth limit.

    Notes:
        - Profile lines are labeled by material index.
        - Slices are drawn as filled red quadrilaterals based on the geometry in `df`.
        - Piezometric line is marked with a blue line and a downward triangle for water level.
        - Distributed loads are shown as downward-pointing purple arrows with a connecting top line.
        - The function uses matplotlib for plotting and shows the plot immediately.

    Returns:
        None
    """

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot profile lines
    for i, line in enumerate(profile_lines):
        xs, ys = zip(*line)
        ax.plot(xs, ys, label=f'Material {i+1}')

    # Plot max depth line
    if max_depth is not None:
        x_vals = [x for line in profile_lines for x, _ in line]
        x_min = min(x_vals)
        x_max = max(x_vals)

        # Solid line
        ax.hlines(max_depth, x_min, x_max, colors='black', linewidth=1.5, label='Max Depth')

        # Down-left diagonal hash marks
        spacing = 5
        length = 4
        angle_rad = np.radians(60)

        dx = length * np.cos(angle_rad)
        dy = length * np.sin(angle_rad)

        x_hashes = np.arange(x_min, x_max, spacing)[1:]  # skip first hash
        for x in x_hashes:
            ax.plot([x, x - dx], [max_depth, max_depth - dy], color='black', linewidth=1)

    # Plot the clipped failure surface if provided
    if failure_surface:
        x_clip, y_clip = zip(*failure_surface.coords)
        ax.plot(x_clip, y_clip, 'k-', linewidth=2, label="Failure Surface")

    # Plot slices using left/right geometry
    for _, row in df.iterrows():
        xs = [row['x_l'], row['x_l'], row['x_r'], row['x_r'], row['x_l']]
        ys = [row['y_lb'], row['y_lt'], row['y_rt'], row['y_rb'], row['y_lb']]
        ax.plot(xs, ys, 'r-')
        ax.fill(xs, ys, color='red', alpha=0.1)

    # Plot piezometric line
    if piezo_line:
        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, 'b-', label="Piezometric Line")
        # Add inverted triangle marker to indicate water level
        mid_index = len(piezo_xs) // 2  # use middle of piezo line
        marker_offset = 2  # adjust this for triangle size
        ax.plot(piezo_xs[mid_index], piezo_ys[mid_index] + marker_offset,
                marker='v', color='b', markersize=8)

    # Plot distributed loads as vertical arrows
    for line in dloads:
        xs = [pt['X'] for pt in line]
        ys = [pt['Y'] for pt in line]
        ns = [pt['Normal'] for pt in line]

        # Interpolation functions
        interp_y = lambda x: np.interp(x, xs, ys)
        interp_n = lambda x: np.interp(x, xs, ns)

        # Create a smooth sequence of x-values for arrows
        x_arrows = np.linspace(min(xs), max(xs), 15)

        # Track arrow tops for connection line
        top_xs = []
        top_ys = []

        for x in x_arrows:
            y = interp_y(x)
            n = interp_n(x)
            arrow_height = min(n / 100, 10)
            y_top = y + arrow_height + 2

            ax.arrow(x, y_top, 0, -arrow_height,
                     head_width=2, head_length=2, fc='purple', ec='purple')

            top_xs.append(x)
            top_ys.append(y_top)

        ax.plot(top_xs, top_ys, color='purple', linewidth=1.5)


    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    ax.grid(False) # change to True if you want grid lines

    if fs is not None:
        ax.set_title(f"Factor of Safety = {fs:.3f}")

    plt.tight_layout()
    plt.show()