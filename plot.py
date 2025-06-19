import matplotlib.pyplot as plt
import numpy as np
from slice import generate_failure_surface
from solve import compute_line_of_thrust
from shapely.geometry import LineString
from matplotlib.lines import Line2D
from matplotlib.path import Path


def get_dload_legend_handler():
    """
    Creates and returns a custom legend entry for distributed loads.
    Returns a tuple of (handler_class, dummy_patch) for use in matplotlib legends.
    """
    # Create a line with built-in arrow marker
    dummy_line = Line2D([0.0, 1.0], [0, 0],  # Two points to define line
                       color='purple', 
                       alpha=0.7, 
                       linewidth=2,
                       marker='>',  # Built-in right arrow marker
                       markersize=6,  # Smaller marker size
                       markerfacecolor='purple',
                       markeredgecolor='purple',
                       drawstyle='steps-post',  # Draw line then marker
                       solid_capstyle='butt')
    
    return None, dummy_line


def plot_profile_lines(ax, profile_lines):
    """
    Plots the profile lines for each material in the slope.

    Parameters:
        ax: matplotlib Axes object
        profile_lines: List of line coordinates representing material boundaries

    Returns:
        None
    """
    for i, line in enumerate(profile_lines):
        xs, ys = zip(*line)
        ax.plot(xs, ys, label=f'Material {i+1}')

def plot_max_depth(ax, profile_lines, max_depth):
    """
    Plots a horizontal line representing the maximum depth limit with hash marks.

    Parameters:
        ax: matplotlib Axes object
        profile_lines: List of line coordinates representing material boundaries
        max_depth: Maximum allowed depth for analysis

    Returns:
        None
    """
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
    """
    Plots the failure surface as a black line.

    Parameters:
        ax: matplotlib Axes object
        failure_surface: Shapely LineString representing the failure surface

    Returns:
        None
    """
    if failure_surface:
        x_clip, y_clip = zip(*failure_surface.coords)
        ax.plot(x_clip, y_clip, 'k-', linewidth=2, label="Failure Surface")

def plot_slices(ax, df, fill=True):
    """
    Plots the slices used in the analysis.

    Parameters:
        ax: matplotlib Axes object
        df: DataFrame containing slice data
        fill: Boolean indicating whether to fill the slices with color

    Returns:
        None
    """
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

def plot_slice_numbers(ax, df):
    """
    Plots the slice number in the middle of each slice at the middle height.
    Numbers are 1-indexed.

    Parameters:
        ax: matplotlib Axes object
        df: DataFrame containing slice data

    Returns:
        None
    """
    if df is not None:
        for _, row in df.iterrows():
            # Calculate middle x-coordinate of the slice
            x_middle = row['x_c']
            
            # Calculate middle height of the slice
            y_middle = (row['y_cb'] + row['y_ct']) / 2
            
            # Plot the slice number (1-indexed)
            slice_number = int(row['slice #'])
            ax.text(x_middle, y_middle, str(slice_number), 
                   ha='center', va='center', fontsize=8, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.8))

def plot_piezo_line(ax, piezo_line):
    """
    Plots the piezometric line with a marker at its midpoint.

    Parameters:
        ax: matplotlib Axes object
        piezo_line: List of coordinates representing the piezometric line

    Returns:
        None
    """
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

def plot_dloads(ax, dloads, data=None, max_height_frac=0.3):
    """
    Plots distributed loads as arrows along the surface.

    Parameters:
        ax: matplotlib Axes object
        dloads: List of distributed load data points with X, Y, and Normal force components
        data: Dictionary containing plot data (needed for ground_surface to calculate slope height)
        max_height_frac: Maximum arrow height as fraction of slope height (default 0.3)

    Returns:
        None
    """
    if not dloads:
        return
    
    # Calculate slope height from ground surface if available
    slope_height = None
    if data and 'ground_surface' in data:
        ground_surface = data['ground_surface']
        if hasattr(ground_surface, 'coords'):
            y_vals = [y for _, y in ground_surface.coords]
            slope_height = max(y_vals) - min(y_vals)
    
    # If we can't get slope height, use a default scaling
    if slope_height is None:
        # Fallback: find max load value across all dloads for scaling
        max_load = 0
        for line in dloads:
            max_load = max(max_load, max(pt['Normal'] for pt in line))
        slope_height = max_load / 10  # Arbitrary scaling if no ground surface available
    
    # Find the maximum load value for scaling
    max_load = 0
    for line in dloads:
        max_load = max(max_load, max(pt['Normal'] for pt in line))
    
    # Calculate maximum arrow height
    max_arrow_height = slope_height * max_height_frac
    
    for line in dloads:
        if len(line) < 2:
            continue
            
        xs = [pt['X'] for pt in line]
        ys = [pt['Y'] for pt in line]
        ns = [pt['Normal'] for pt in line]
        
        # Process line segments
        for i in range(len(line) - 1):
            x1, y1, n1 = xs[i], ys[i], ns[i]
            x2, y2, n2 = xs[i+1], ys[i+1], ns[i+1]
            
            # Calculate segment direction (perpendicular to this segment)
            dx = x2 - x1
            dy = y2 - y1
            segment_length = np.sqrt(dx**2 + dy**2)
            
            if segment_length == 0:
                continue
                
            # Normalize the segment direction
            dx_norm = dx / segment_length
            dy_norm = dy / segment_length
            
            # Perpendicular direction (rotate 90 degrees CCW)
            perp_dx = -dy_norm
            perp_dy = dx_norm
            
            # Generate arrows along this segment
            num_arrows = max(3, min(10, int(segment_length / 5)))  # Adaptive number of arrows
            t_values = np.linspace(0, 1, num_arrows)
            
            # Store arrow top points for connecting line
            top_xs = []
            top_ys = []
            
            # Add start point if it's the first segment and load is zero
            if i == 0 and n1 == 0:
                top_xs.append(x1)
                top_ys.append(y1)
            
            for t in t_values:
                # Interpolate position along segment
                x = x1 + t * dx
                y = y1 + t * dy
                
                # Interpolate load value
                n = n1 + t * (n2 - n1)
                
                # Scale arrow height
                if max_load > 0:
                    arrow_height = (n / max_load) * max_arrow_height
                else:
                    arrow_height = 0
                
                # For very small arrows, just store surface point for connecting line
                if arrow_height < 0.5:
                    top_xs.append(x)
                    top_ys.append(y)
                    continue
                
                # Calculate arrow dimensions
                head_length = arrow_height * 0.2
                head_width = arrow_height * 0.15
                
                # Calculate arrow start point (above surface)
                arrow_start_x = x + perp_dx * arrow_height
                arrow_start_y = y + perp_dy * arrow_height
                
                # Store points for connecting line
                top_xs.append(arrow_start_x)
                top_ys.append(arrow_start_y)
                
                # Draw arrow - extend all the way to surface point
                ax.arrow(arrow_start_x, arrow_start_y, 
                        x - arrow_start_x, y - arrow_start_y,
                        head_width=head_width, head_length=head_length, 
                        fc='purple', ec='purple', alpha=0.7,
                        length_includes_head=True)
            
            # Add end point if it's the last segment and load is zero
            if i == len(line) - 2 and n2 == 0:
                top_xs.append(x2)
                top_ys.append(y2)
            
            # Draw connecting line at arrow tops
            if top_xs:
                ax.plot(top_xs, top_ys, color='purple', linewidth=1.5, alpha=0.8)
        
        # Draw the surface line itself
        ax.plot(xs, ys, color='purple', linewidth=1.5, alpha=0.8)

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
    """
    Plots a non-circular failure surface.

    Parameters:
        ax: matplotlib Axes object
        non_circ: List of coordinates representing the non-circular failure surface

    Returns:
        None
    """
    xs, ys = zip(*non_circ)
    ax.plot(xs, ys, 'r--', label='Non-Circular Surface')

def plot_material_table(ax, materials, xloc=0.6, yloc=0.7):
    """
    Adds a material properties table to the plot.

    Parameters:
        ax: matplotlib Axes object
        materials: List of material property dictionaries
        xloc: x-location of table (0-1)
        yloc: y-location of table (0-1)

    Returns:
        None
    """
    if not materials:
        return

    # Check material options
    options = set(mat['option'] for mat in materials)

    # Decide column headers - need 5 columns to match the data
    if options == {'mc'}:
        col_labels = ["Mat", "Name", "γ", "c", "φ"]
    elif options == {'cp'}:
        col_labels = ["Mat", "Name", "γ", "cp", "rₑ"]
    else:
        col_labels = ["Mat", "Name", "γ", "c / cp", "φ / rₑ"]

    # Build table rows
    table_data = []
    for idx, mat in enumerate(materials):
        name = mat['name']
        gamma = mat['gamma']
        option = mat['option']
        if option == 'mc':
            c = mat['c']
            phi = mat['phi']
            row = [idx+1, name, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}"]
        elif option == 'cp':
            cp = mat['cp']
            r_elev = mat['r_elev']
            row = [idx+1, name, f"{gamma:.1f}", f"{cp:.2f}", f"{r_elev:.1f}"]
        else:
            row = [idx+1, name, f"{gamma:.1f}", "-", "-"]
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

def plot_base_stresses(ax, df, scale_frac=0.5, alpha=0.3):
    """
    Plots the effective normal stresses and pore pressures on the base of each slice.

    Parameters:
        ax: matplotlib Axes object
        df: DataFrame containing slice data
        scale_frac: Scaling factor for stress visualization
        alpha: Transparency value for filled areas

    Returns:
        None
    """
    u = df['u'].values
    n_eff = df['n_eff'].values
    dl = df['dl'].values
    sigma_eff = n_eff / dl

    heights = df['y_ct'] - df['y_cb']
    max_ht = heights.max() if not heights.empty else 1.0
    max_bar_len = max_ht * scale_frac

    max_stress = np.max(np.abs(sigma_eff)) if len(sigma_eff) > 0 else 1.0
    max_u = np.max(u) if len(u) > 0 else 1.0

    for i, (index, row) in enumerate(df.iterrows()):
        if i >= len(sigma_eff):
            break

        x1, y1 = row['x_l'], row['y_lb']
        x2, y2 = row['x_r'], row['y_rb']

        stress = sigma_eff[i]
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

        ax.fill(poly_x, poly_y, facecolor='none', edgecolor='red' if stress <= 0 else 'limegreen', hatch='.....',
                linewidth=1)

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

def plot_thrust_line(ax, thrust_line: LineString,
                    color: str = 'red',
                    linestyle: str = '--',
                    linewidth: float = 1,
                    label: str = 'Line of Thrust'):
    """
    Plots the line of thrust on the slope.

    Parameters:
        ax: matplotlib Axes object
        thrust_line: Shapely LineString representing the line of thrust
        color: Color of the line
        linestyle: Style of the line
        linewidth: Width of the line
        label: Label for the line in the legend

    Returns:
        None
    """
    # extract x,y coords
    xs, ys = zip(*list(thrust_line.coords))
    ax.plot(xs, ys,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label)

def compute_ylim(data, df, scale_frac=0.5, pad_fraction=0.1):
    """
    Computes y‐axis limits for the solution plot, ensuring that:
      • All profile_lines are included,
      • The max_depth (toe) is included below,
      • Stress bars (length = max slice height * scale_frac) fit inside,
      • And an extra pad_fraction of breathing‐room is added.

    Parameters:
        data: dict containing at least 'profile_lines' and 'max_depth'
        df: pandas.DataFrame with slice data, must have 'y_lt' and 'y_lb' for stress‐bar sizing
        scale_frac: fraction of max slice height used when drawing stress bars
        pad_fraction: fraction of total range to pad above/below finally

    Returns:
        (y_min, y_max) suitable for ax.set_ylim(...)
    """
    import numpy as np

    y_vals = []

    # 1) collect all profile line elevations
    for line in data.get('profile_lines', []):
        if hasattr(line, "xy"):
            _, ys = line.xy
        else:
            _, ys = zip(*line)
        y_vals.extend(ys)

    # 2) explicitly include the deepest allowed depth
    if "max_depth" in data and data["max_depth"] is not None:
        y_vals.append(data["max_depth"])

    if not y_vals:
        return 0.0, 1.0

    y_min = min(y_vals)
    y_max = max(y_vals)

    # 3) ensure the largest stress bar will fit
    #    stress‐bar length = scale_frac * slice height
    heights = df["y_lt"] - df["y_lb"]
    if not heights.empty:
        max_bar = heights.max() * scale_frac
        y_min -= max_bar
        y_max += max_bar

    # 4) add a final small pad
    pad = (y_max - y_min) * pad_fraction
    return y_min - pad, y_max + pad

# ========== FOR PLOTTING INPUT DATA  =========

def plot_inputs(data, title="Slope Geometry and Inputs", width=12, height=6):
    """
    Creates a plot showing the slope geometry and input parameters.

    Parameters:
        data: Dictionary containing plot data
        title: Title for the plot
        width: Width of the plot in inches
        height: Height of the plot in inches

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))

    # Plot contents
    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'], data)
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

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    
    # Add distributed load to legend if present
    if data['dloads']:
        handler_class, dummy_line = get_dload_legend_handler()
        handles.append(dummy_line)
        labels.append('Distributed Load')
    
    ax.legend(
        handles=handles,
        labels=labels,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.12),
        ncol=2
    )

    ax.set_title(title)

    plt.tight_layout()
    plt.show()

# ========== Main Plotting Function =========

def plot_solution(data, df, failure_surface, results, width=12, height=7, slice_numbers=True):
    """
    Creates a plot showing the slope stability analysis solution.

    Parameters:
        data: Dictionary containing plot data
        df: DataFrame containing slice data
        failure_surface: Shapely LineString representing the failure surface
        results: Dictionary containing analysis results
        width: Width of the plot in inches
        height: Height of the plot in inches

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)

    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_slices(ax, df, fill=False)
    plot_failure_surface(ax, failure_surface)
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'], data)
    plot_tcrack_surface(ax, data['tcrack_surface'])
    if slice_numbers:
        plot_slice_numbers(ax, df)

    alpha = 0.3
    if results['method'] == 'spencer':
        FS = results['FS']
        thrust_line = compute_line_of_thrust(df, FS,  debug=False)
        plot_thrust_line(ax, thrust_line)

    plot_base_stresses(ax, df, alpha=alpha)

    import matplotlib.patches as mpatches
    normal_patch = mpatches.Patch(facecolor='none', edgecolor='green', hatch='.....', label="Eff Normal Stress (σ')")
    pore_patch = mpatches.Patch(color='blue', alpha=alpha, label='Pore Pressure (u)')

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    handles.extend([normal_patch, pore_patch])
    labels.extend(["Eff Normal Stress (σ')", 'Pore Pressure (u)'])
    
    # Add distributed load to legend if present
    if data['dloads']:
        handler_class, dummy_line = get_dload_legend_handler()
        handles.append(dummy_line)
        labels.append('Distributed Load')
    
    ax.legend(
        handles=handles,
        labels=labels,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.12),
        ncol=2
    )

    # Add vertical space below for the legend
    plt.subplots_adjust(bottom=0.2)
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
    elif method == 'janbu':
        fo = results['fo']
        title = f'Janbu-Corrected: FS = {fs:.3f}, fo = {fo:.2f}'
    elif method == 'corps_engineers':
        theta = results['theta']
        title = f'Corps Engineers: FS = {fs:.3f}, θ = {theta:.2f}°'
    elif method == 'lowe_karafiath':
        title = f'Lowe & Karafiath: FS = {fs:.3f}'
    ax.set_title(title)

    # zoom y‐axis to just cover the slope and depth, with a little breathing room (thrust line can be outside)
    ymin, ymax = compute_ylim(data, df, pad_fraction=0.05)
    ax.set_ylim(ymin, ymax)

    plt.tight_layout()
    plt.show()

# ========== Functions for Search Results =========

def plot_failure_surfaces(ax, fs_cache):
    """
    Plots all failure surfaces from the factor of safety cache.

    Parameters:
        ax: matplotlib Axes object
        fs_cache: List of dictionaries containing failure surface data and FS values

    Returns:
        None
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
    Plots the centers of circular failure surfaces.

    Parameters:
        ax: matplotlib Axes object
        fs_cache: List of dictionaries containing circle center data

    Returns:
        None
    """
    for result in fs_cache:
        ax.plot(result['Xo'], result['Yo'], 'ko', markersize=3, alpha=0.6)

def plot_search_path(ax, search_path):
    """
    Plots the search path used to find the critical failure surface.

    Parameters:
        ax: matplotlib Axes object
        search_path: List of dictionaries containing search path coordinates

    Returns:
        None
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
    Creates a plot showing the results of a circular failure surface search.

    Parameters:
        data: Dictionary containing plot data
        fs_cache: List of dictionaries containing failure surface data and FS values
        search_path: List of dictionaries containing search path coordinates
        highlight_fs: Boolean indicating whether to highlight the critical failure surface
        width: Width of the plot in inches
        height: Height of the plot in inches

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))

    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'], data)
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

def plot_noncircular_search_results(data, fs_cache, search_path=None, highlight_fs=True, width=12, height=7):
    """
    Creates a plot showing the results of a non-circular failure surface search.

    Parameters:
        data: Dictionary containing plot data
        fs_cache: List of dictionaries containing failure surface data and FS values
        search_path: List of dictionaries containing search path coordinates
        highlight_fs: Boolean indicating whether to highlight the critical failure surface
        width: Width of the plot in inches
        height: Height of the plot in inches

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))

    # Plot basic profile elements
    plot_profile_lines(ax, data['profile_lines'])
    plot_max_depth(ax, data['profile_lines'], data['max_depth'])
    plot_piezo_line(ax, data['piezo_line'])
    plot_dloads(ax, data['dloads'], data)
    plot_tcrack_surface(ax, data['tcrack_surface'])

    # Plot all failure surfaces from cache
    for i, result in reversed(list(enumerate(fs_cache))):
        surface = result['failure_surface']
        if surface is None or surface.is_empty:
            continue
        x, y = zip(*surface.coords)
        color = 'red' if i == 0 else 'gray'
        lw = 2 if i == 0 else 1
        ax.plot(x, y, color=color, linestyle='-', linewidth=lw, alpha=1.0 if i == 0 else 0.6)

    # Plot search path if provided
    if search_path:
        for i in range(len(search_path) - 1):
            start = search_path[i]
            end = search_path[i + 1]
            # For non-circular search, we need to plot the movement of each point
            start_points = np.array(start['points'])
            end_points = np.array(end['points'])
            
            # Plot arrows for each moving point
            for j in range(len(start_points)):
                dx = end_points[j, 0] - start_points[j, 0]
                dy = end_points[j, 1] - start_points[j, 1]
                if abs(dx) > 1e-6 or abs(dy) > 1e-6:  # Only plot if point moved
                    ax.arrow(start_points[j, 0], start_points[j, 1], dx, dy,
                            head_width=1, head_length=2, fc='green', ec='green',
                            length_includes_head=True, alpha=0.6)

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