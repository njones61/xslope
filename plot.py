import matplotlib.pyplot as plt
import numpy as np
from slice import generate_failure_surface, calculate_normal_stresses
from solve import extract_spencer_Q, compute_line_of_thrust_sweep

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

def plot_slices(ax, df, fill=True):
    if df is not None:
        for _, row in df.iterrows():
            if fill:
                xs = [row['x_l'], row['x_l'], row['x_r'], row['x_r'], row['x_l']]
                ys = [row['y_lb'], row['y_lt'], row['y_rt'], row['y_rb'], row['y_lb']]
                ax.plot(xs, ys, 'r-')
                ax.fill(xs, ys, color='red', alpha=0.1)
            else:
                ax.plot([row['x_l'], row['x_l']], [row['y_lb'], row['y_lt']], 'k-', linewidth=0.5)
                ax.plot([row['x_r'], row['x_r']], [row['y_rb'], row['y_rt']], 'k-', linewidth=0.5)

def plot_piezo_line(ax, piezo_line):
    if piezo_line:
        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, 'b-', label="Piezometric Line")
        mid_index = len(piezo_xs) // 2
        marker_offset = 2
        ax.plot(piezo_xs[mid_index], piezo_ys[mid_index] + marker_offset,
                marker='v', color='b', markersize=8)

def plot_tcrack_surface(ax, tcrack_surface):
    """
    Plots the tension crack surface as a thin dashed red line.

    Parameters:
        ax: matplotlib Axes object
        tcrack_surface: Shapely LineString

    Returns:
        None
    """
    if tcrack_surface is None:
        return

    x_vals, y_vals = tcrack_surface.xy
    ax.plot(x_vals, y_vals, linestyle='--', color='red', linewidth=1.0, label='Tension Crack Depth')

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

def plot_base_stresses(ax, df, FS, scale_frac=0.5, alpha=0.3):
    """
    Plots trapezoidal bars at the base of each slice representing normal stresses and pore pressures.

    Parameters:
        ax: matplotlib axis
        df: DataFrame with slice geometry (must contain x_l, y_lb, x_r, y_rb, y_ct, y_cb, u)
        FS: Factor of Safety (used for stress calculation)
        scale_frac: fraction of max slice height used to scale stress bar length

    Returns:
        None
    """
    import numpy as np

    N = calculate_normal_stresses(df, FS)
    u = df['u'].values
    N_eff = N - u

    heights = df['y_ct'] - df['y_cb']
    max_ht = heights.max() if not heights.empty else 1.0
    max_bar_len = max_ht * scale_frac

    max_stress = np.max(np.abs(N_eff)) if len(N_eff) > 0 else 1.0
    max_u = np.max(u) if len(u) > 0 else 1.0

    for i, (index, row) in enumerate(df.iterrows()):
        if i >= len(N_eff):
            break

        x1, y1 = row['x_l'], row['y_lb']
        x2, y2 = row['x_r'], row['y_rb']
        stress = N_eff[i]
        pore = u[i]

        dx = x2 - x1
        dy = y2 - y1
        length = np.hypot(dx, dy)
        if length == 0:
            continue

        nx = -dy / length
        ny = dx / length

        # --- Normal stress trapezoid ---
        bar_len = (abs(stress) / max_stress) * max_bar_len
        direction = -np.sign(stress)

        x1_top = x1 + direction * bar_len * nx
        y1_top = y1 + direction * bar_len * ny
        x2_top = x2 + direction * bar_len * nx
        y2_top = y2 + direction * bar_len * ny

        poly_x = [x1, x2, x2_top, x1_top]
        poly_y = [y1, y2, y2_top, y1_top]

        # ax.fill(poly_x, poly_y, color='red' if stress <= 0 else 'green', alpha=alpha, edgecolor='k', linewidth=0.5)
        ax.fill(poly_x, poly_y, facecolor='none', edgecolor='red' if stress <= 0 else 'limegreen', hatch='.....', linewidth=1)

        # --- Pore pressure trapezoid ---
        u_len = (pore / max_stress) * max_bar_len
        u_dir = -1  # always into the base

        ux1_top = x1 + u_dir * u_len * nx
        uy1_top = y1 + u_dir * u_len * ny
        ux2_top = x2 + u_dir * u_len * nx
        uy2_top = y2 + u_dir * u_len * ny

        poly_ux = [x1, x2, ux2_top, ux1_top]
        poly_uy = [y1, y2, uy2_top, uy1_top]

        ax.fill(poly_ux, poly_uy, color='blue', alpha=alpha, edgecolor='k', linewidth=1)
        #ax.fill(poly_ux, poly_uy, facecolor='none', edgecolor='blue', hatch='.....', linewidth=1)



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
    plot_tcrack_surface(ax, data['tcrack_surface'])

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

    ax.set_title(title)

    plt.tight_layout()
    plt.show()

# ========== Main Plotting Function =========

def plot_thrust_line(ax, df, y_thrust, color='magenta', label='Line of Thrust', **kwargs):
    """
    Plots the line of thrust over the slice cross section using interslice boundary locations.
    Assumes y_thrust is the distance above the base of the slice at each boundary.

    Parameters:
        ax (matplotlib.axes.Axes): Axis to draw the plot on
        df (pd.DataFrame): Slice data with 'x_l', 'x_r', 'y_lb', and 'y_rb' columns
        y_thrust (list or array): Distance above base (length = n+1)
        color (str): Line color
        label (str): Legend label for the line
        kwargs: Additional arguments passed to ax.plot()

    Returns:
        None
    """
    import numpy as np

    # X coordinates of slice boundaries
    x_vals = list(df['x_l']) + [df['x_r'].iloc[-1]]

    # Corresponding y base elevations at those boundaries
    y_base = list(df['y_lb']) + [df['y_rb'].iloc[-1]]

    # Convert to numpy arrays
    x_vals = np.array(x_vals)
    y_base = np.array(y_base)
    y_thrust = np.array(y_thrust)

    # Add base elevation to thrust height
    y_plot = y_base + y_thrust

    # Mask invalid values
    mask = np.isfinite(y_plot)
    ax.plot(x_vals[mask], y_plot[mask], color=color, label=label, **kwargs)

def print_thrust_line_debug(df, y_left, y_right, y_avg, digits=2):
    """
    Prints a formatted comparison table of left/right/average thrust line values.

    Parameters:
        df (pd.DataFrame): For determining slice count
        y_left (list): y-values from left-to-right sweep (length = n+1)
        y_right (list): y-values from right-to-left sweep (length = n+1)
        y_avg (list): averaged y-values (length = n+1)
        digits (int): Decimal precision for printing
    """
    fmt = f"{{:>{digits+5}.{digits}f}}"
    n = len(df)

    print(f"{'Slice':>5} | {'y_left':>{digits+7}} | {'y_right':>{digits+7}} | {'y_avg':>{digits+7}} | {'Δ = y_L - y_R':>{digits+9}}")
    print("-" * (6 + 4*(digits+7) + 1))

    for i in range(n + 1):
        yl = y_left[i]
        yr = y_right[i]
        ya = y_avg[i]

        if all(v is not None and np.isfinite(v) for v in [yl, yr, ya]):
            delta = yl - yr
            print(f"{i:5d} | {fmt.format(yl)} | {fmt.format(yr)} | {fmt.format(ya)} | {fmt.format(delta)}")
        else:
            print(f"{i:5d} | {'None':>{digits+7}} | {'None':>{digits+7}} | {'None':>{digits+7}} | {'':>{digits+9}}")

def plot_solution(data, df, failure_surface, results, width=12, height=7):
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)

    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_slices(ax, df, fill=False)
    plot_failure_surface(ax, failure_surface)
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'])
    plot_tcrack_surface(ax, data['tcrack_surface'])


    if results['method'] == 'spencer':
        FS = results['FS']
        theta = results['theta']
        Q = extract_spencer_Q(df, FS, theta, debug=True)
        y_left, y_right, y_avg = compute_line_of_thrust_sweep(df, Q, theta)
        plot_thrust_line(ax, df, y_avg, color='magenta', linestyle='-', linewidth=2)
        print_thrust_line_debug(df, y_left, y_right, y_avg)

    alpha = 0.3
    plot_base_stresses(ax, df, results['FS'], alpha=alpha)

    import matplotlib.patches as mpatches
    normal_patch = mpatches.Patch(facecolor='none', edgecolor='green', hatch='.....',  label="Eff Normal Stress (σ')")
    pore_patch = mpatches.Patch(color='blue', alpha=alpha, label='Pore Pressure (u)')

    # Add these to the legend
    # ax.legend(handles=ax.get_legend_handles_labels()[0] + [normal_patch, pore_patch])
    # Add legend below the plot
    ax.legend(
        handles=ax.get_legend_handles_labels()[0] + [normal_patch, pore_patch],
        loc='upper center',
        bbox_to_anchor=(0.5, -0.12),  # x=centered, y=slightly below axes
        ncol=2
    )

    # Add vertical space below for the legend
    plt.subplots_adjust(bottom=0.2)  # Increase if needed
    ax.set_aspect('equal')

    fs = results['FS']
    method = results['method']
    if method == 'oms':
        title = f'OMS: FS = {fs:.3f}'
    elif method == 'bishop':
        title = f'Bishop: FS = {fs:.3f}'
    elif method == 'spencer':
        theta = results['theta']
        title = f'Spencer: FS = {fs:.3f}, θ = {theta:.2f}°'
    elif method == 'janbu_corrected':
        fo = results['fo']
        title = f'Janbu-Corrected: FS = {fs:.3f}, fo = {fo:.2f}'
    ax.set_title(title)

    plt.tight_layout()
    plt.show()

# ========== Functions for Search Results =========

def plot_failure_surfaces(ax, fs_cache):
    """
    Plots all failure surfaces.
    Critical surface (lowest FS) is plotted in red, others in gray.
    Drawing order is reversed so the critical surface is on top.
    """
    for i, result in reversed(list(enumerate(fs_cache))):
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


def plot_circular_search_results(data, fs_cache, search_path=None, highlight_fs=True, width=12, height=7):
    """
    Main function to plot slope geometry and circular search results.

    Parameters:
        data (dict): Global data dictionary
        fs_cache (list of dict): Search results sorted by FS
        search_path (list of dict, optional): Path of search refinements
        highlight_fs (bool): Whether to label critical FS on plot
    """
    fig, ax = plt.subplots(figsize=(width, height))

    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'])
    plot_tcrack_surface(ax, data['tcrack_surface'])

    plot_failure_surfaces(ax, fs_cache)
    plot_circle_centers(ax, fs_cache)
    if search_path:
        plot_search_path(ax, search_path)

    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)
    ax.legend()

    if highlight_fs and fs_cache:
        critical_fs = fs_cache[0]['FS']
        ax.set_title(f"Critical Factor of Safety = {critical_fs:.3f}")

    plt.tight_layout()
    plt.show()