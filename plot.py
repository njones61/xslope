import matplotlib.pyplot as plt
import numpy as np
from slice import generate_failure_surface


def plot_profile_lines(ax, profile_lines):
    for i, line in enumerate(profile_lines):
        xs, ys = zip(*line)
        ax.plot(xs, ys, label=f'Material {i+1}')

def plot_max_depth(ax, profile_lines, max_depth):
    if max_depth is None:
        return
    x_vals = [x for line in profile_lines for x, _ in line]
    x_min = min(x_vals)
    x_max = max(x_vals)
    ax.hlines(max_depth, x_min, x_max, colors='black', linewidth=1.5, label='Max Depth')

    spacing = 5
    length = 4
    angle_rad = np.radians(60)
    dx = length * np.cos(angle_rad)
    dy = length * np.sin(angle_rad)
    x_hashes = np.arange(x_min, x_max, spacing)[1:]
    for x in x_hashes:
        ax.plot([x, x - dx], [max_depth, max_depth - dy], color='black', linewidth=1)

def plot_failure_surface(ax, failure_surface):
    if failure_surface:
        x_clip, y_clip = zip(*failure_surface.coords)
        ax.plot(x_clip, y_clip, 'k-', linewidth=2, label="Failure Surface")

def plot_slices(ax, df):
    if df is not None:
        for _, row in df.iterrows():
            xs = [row['x_l'], row['x_l'], row['x_r'], row['x_r'], row['x_l']]
            ys = [row['y_lb'], row['y_lt'], row['y_rt'], row['y_rb'], row['y_lb']]
            ax.plot(xs, ys, 'r-')
            ax.fill(xs, ys, color='red', alpha=0.1)

def plot_piezo_line(ax, piezo_line):
    if piezo_line:
        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, 'b-', label="Piezometric Line")
        mid_index = len(piezo_xs) // 2
        marker_offset = 2
        ax.plot(piezo_xs[mid_index], piezo_ys[mid_index] + marker_offset,
                marker='v', color='b', markersize=8)

def plot_dloads(ax, dloads):
    for line in dloads:
        xs = [pt['X'] for pt in line]
        ys = [pt['Y'] for pt in line]
        ns = [pt['Normal'] for pt in line]

        interp_y = lambda x: np.interp(x, xs, ys)
        interp_n = lambda x: np.interp(x, xs, ns)

        x_arrows = np.linspace(min(xs), max(xs), 15)
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

def plot_circles(ax, data):
    """
    Plots starting circles with center markers and arrows.

    Parameters:
        ax (matplotlib axis): The plotting axis
        circles (list of dicts): List of circles with 'Xo', 'Yo', 'R'

    Returns:
        None
    """
    circles = data['circles']
    for circle in circles:
        Xo = circle['Xo']
        Yo = circle['Yo']
        R = circle['R']
        # theta = np.linspace(0, 2 * np.pi, 100)
        # x_circle = Xo + R * np.cos(theta)
        # y_circle = Yo + R * np.sin(theta)
        # ax.plot(x_circle, y_circle, 'r--', label='Circle')

        # Plot the portion of the circle in the slope
        ground_surface = data['ground_surface']
        success, result = generate_failure_surface(ground_surface, circular=True, circle=circle)
        if not success:
            continue  # or handle error
        # result = (x_min, x_max, y_left, y_right, clipped_surface)
        x_min, x_max, y_left, y_right, clipped_surface = result
        x_clip, y_clip = zip(*clipped_surface.coords)
        ax.plot(x_clip, y_clip, 'r--', label="Circle")

        # Center marker
        ax.plot(Xo, Yo, 'r+', markersize=10)

        # Arrow direction: point from center to midpoint of failure surface
        mid_idx = len(x_clip) // 2
        x_mid = x_clip[mid_idx]
        y_mid = y_clip[mid_idx]

        dx = x_mid - Xo
        dy = y_mid - Yo

        # Normalize direction vector
        length = np.hypot(dx, dy)
        if length != 0:
            dx /= length
            dy /= length

        # Shorten shaft length slightly
        shaft_length = R - 5

        ax.arrow(Xo, Yo, dx * shaft_length, dy * shaft_length,
                 head_width=5, head_length=5, fc='red', ec='red')

def plot_non_circ(ax, non_circ):
    xs, ys = zip(*non_circ)
    ax.plot(xs, ys, 'r--', label='Non-Circular Surface')

def plot_material_table(ax, materials, xloc=0.6, yloc=0.7):
    """
    Adds a material properties table to the plot.
    """
    if not materials:
        return

    # Check material options
    options = set(mat['option'] for mat in materials)

    # Decide column headers
    if options == {'mc'}:
        col_labels = ["Mat", "γ", "c", "φ"]
    elif options == {'cp'}:
        col_labels = ["Mat", "γ", "cp", "rₑ"]
    else:
        col_labels = ["Mat", "γ", "c / cp", "φ / rₑ"]

    # Build table rows
    table_data = []
    for idx, mat in enumerate(materials):
        gamma = mat['gamma']
        option = mat['option']
        if option == 'mc':
            c = mat['c']
            phi = mat['phi']
            row = [idx+1, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}"]
        elif option == 'cp':
            cp = mat['cp']
            r_elev = mat['r_elev']
            row = [idx+1, f"{gamma:.1f}", f"{cp:.2f}", f"{r_elev:.1f}"]
        else:
            row = [idx+1, f"{gamma:.1f}", "-", "-"]
        table_data.append(row)

    # Add the table
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, 0.2, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(8)

# ========== FOR PLOTTING INPUT DATA  =========


def plot_inputs(data, title="Slope Geometry and Inputs", width=12, height=6):
    """
    Simple clean slope plot with equal aspect ratio and normal legend.
    """
    fig, ax = plt.subplots(figsize=(width, height))

    # Plot contents
    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'])

    if data['circular']:
        plot_circles(ax, data)
    else:
        plot_non_circ(ax, data['non_circ'])

    plot_material_table(ax, data['materials'], xloc=0.62) # Adjust this so that it fits with the legend

    ax.set_aspect('equal')  # ✅ Equal aspect
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)

    # Normal legend inside plot
    ax.legend()

    # # Adjust dynamic height based on geometry
    # x_min, x_max = ax.get_xlim()
    # y_min, y_max = ax.get_ylim()
    # x_span = x_max - x_min
    # y_span = y_max - y_min
    #
    # aspect = y_span / x_span
    # dynamic_height = width * aspect + 1
    # fig.set_size_inches(width, dynamic_height)

    ax.set_title(title)

    plt.tight_layout()
    plt.show()

# ========== Main Plotting Function =========

def plot_slope(data, df=None, failure_surface=None, fs=None):
    fig, ax = plt.subplots(figsize=(10, 6))

    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_failure_surface(ax, failure_surface)
    plot_slices(ax, df)
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'])
    if df is None and data['circular']:
        plot_circles(ax, data['circles'])

    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    ax.grid(False)

    if fs is not None:
        ax.set_title(f"Factor of Safety = {fs:.3f}")

    plt.tight_layout()
    plt.show()

# ========== Functions for Search Results =========

def plot_failure_surfaces(ax, fs_cache):
    """
    Plots all failure surfaces.
    Critical surface (lowest FS) is plotted in red, others in gray.
    """
    for i, result in enumerate(fs_cache):
        surface = result['failure_surface']
        if surface is None or surface.is_empty:
            continue
        x, y = zip(*surface.coords)
        color = 'red' if i == 0 else 'gray'
        lw = 2 if i == 0 else 1
        ax.plot(x, y, color=color, linestyle='-', linewidth=lw, alpha=1.0 if i == 0 else 0.6)


def plot_circle_centers(ax, fs_cache):
    """
    Plots a dot at each circle center tried during the search.
    """
    for result in fs_cache:
        ax.plot(result['Xo'], result['Yo'], 'ko', markersize=3, alpha=0.6)


def plot_search_path(ax, search_path):
    """
    Plots small arrows showing the search refinement path.
    """
    if len(search_path) < 2:
        return  # need at least two points to draw an arrow

    for i in range(len(search_path) - 1):
        start = search_path[i]
        end = search_path[i + 1]
        dx = end['x'] - start['x']
        dy = end['y'] - start['y']
        ax.arrow(start['x'], start['y'], dx, dy,
                 head_width=1, head_length=2, fc='green', ec='green', length_includes_head=True)


def plot_circular_search_results(data, fs_cache, search_path=None, highlight_fs=True):
    """
    Main function to plot slope geometry and circular search results.

    Parameters:
        data (dict): Global data dictionary
        fs_cache (list of dict): Search results sorted by FS
        search_path (list of dict, optional): Path of search refinements
        highlight_fs (bool): Whether to label critical FS on plot
    """
    fig, ax = plt.subplots(figsize=(12, 7))

    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'])

    plot_failure_surfaces(ax, fs_cache)
    plot_circle_centers(ax, fs_cache)
    if search_path:
        plot_search_path(ax, search_path)

    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True)
    ax.legend()

    if highlight_fs and fs_cache:
        critical_fs = fs_cache[0]['FS']
        ax.set_title(f"Critical Factor of Safety = {critical_fs:.3f}")

    plt.tight_layout()
    plt.show()