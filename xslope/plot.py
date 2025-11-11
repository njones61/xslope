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

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.path import Path
from shapely.geometry import LineString

from .slice import generate_failure_surface

# Configure matplotlib for better text rendering
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 10
})

# Consistent color for materials (Tableau tab10)
def get_material_color(idx):
    tableau_colors = plt.get_cmap('tab10').colors  # 10 distinct colors
    return tableau_colors[idx % len(tableau_colors)]

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
        ax.plot(xs, ys, color=get_material_color(i), linewidth=1, label=f'Profile {i+1}')

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

def plot_slices(ax, slice_df, fill=True):
    """
    Plots the slices used in the analysis.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data
        fill: Boolean indicating whether to fill the slices with color

    Returns:
        None
    """
    if slice_df is not None:
        for _, row in slice_df.iterrows():
            if fill:
                xs = [row['x_l'], row['x_l'], row['x_r'], row['x_r'], row['x_l']]
                ys = [row['y_lb'], row['y_lt'], row['y_rt'], row['y_rb'], row['y_lb']]
                ax.plot(xs, ys, 'r-')
                ax.fill(xs, ys, color='red', alpha=0.1)
            else:
                ax.plot([row['x_l'], row['x_l']], [row['y_lb'], row['y_lt']], 'k-', linewidth=0.5)
                ax.plot([row['x_r'], row['x_r']], [row['y_rb'], row['y_rt']], 'k-', linewidth=0.5)

def plot_slice_numbers(ax, slice_df):
    """
    Plots the slice number in the middle of each slice at the middle height.
    Numbers are 1-indexed.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data

    Returns:
        None
    """
    if slice_df is not None:
        for _, row in slice_df.iterrows():
            # Calculate middle x-coordinate of the slice
            x_middle = row['x_c']
            
            # Calculate middle height of the slice
            y_middle = (row['y_cb'] + row['y_ct']) / 2
            
            # Plot the slice number (1-indexed)
            slice_number = int(row['slice #'])
            ax.text(x_middle, y_middle, str(slice_number), 
                   ha='center', va='center', fontsize=8, fontweight='bold',
                   bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.8))

def plot_piezo_line(ax, slope_data):
    """
    Plots the piezometric line(s) with markers at their midpoints.

    Parameters:
        ax: matplotlib Axes object
        data: Dictionary containing plot data with 'piezo_line' and optionally 'piezo_line2'

    Returns:
        None
    """
    
    def plot_single_piezo_line(ax, piezo_line, color, label):
        """Internal function to plot a single piezometric line"""
        if not piezo_line:
            return
            
        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, color=color, linewidth=2, label=label)
        
        # Find middle x-coordinate and corresponding y value
        x_min, x_max = min(piezo_xs), max(piezo_xs)
        mid_x = (x_min + x_max) / 2
        
        # Interpolate y value at mid_x
        from scipy.interpolate import interp1d
        if len(piezo_xs) > 1:
            f = interp1d(piezo_xs, piezo_ys, kind='linear', bounds_error=False, fill_value='extrapolate')
            mid_y = f(mid_x)
            ax.plot(mid_x, mid_y + 6, marker='v', color=color, markersize=8)
    
    # Plot both piezometric lines
    plot_single_piezo_line(ax, slope_data.get('piezo_line'), 'b', "Piezometric Line")
    plot_single_piezo_line(ax, slope_data.get('piezo_line2'), 'skyblue', "Piezometric Line 2")

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

def plot_dloads(ax, slope_data):
    """
    Plots distributed loads as arrows along the surface.
    """
    gamma_w = slope_data['gamma_water']
    ground_surface = slope_data['ground_surface']

    def plot_single_dload_set(ax, dloads, color, label):
        """Internal function to plot a single set of distributed loads"""
        if not dloads:
            return
            
        # find the max horizontal length of the ground surface
        max_horizontal_length_ground = 0
        for pt in ground_surface.coords:
            max_horizontal_length_ground = max(max_horizontal_length_ground, pt[0])

        arrow_spacing = max_horizontal_length_ground / 60

        # find the max dload value
        max_dload = 0
        for line in dloads:
            max_dload = max(max_dload, max(pt['Normal'] for pt in line))

        arrow_height = max_dload / gamma_w
        head_length = arrow_height / 12
        head_width = head_length * 0.8
        
        # Find the maximum load value for scaling
        max_load = 0
        for line in dloads:
            max_load = max(max_load, max(pt['Normal'] for pt in line))
        
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
                dx_abs = abs(x2 - x1)
                num_arrows = max(1, int(round(dx_abs / arrow_spacing)))
                if dx_abs == 0:
                    t_values = np.array([0.0, 1.0])
                else:
                    t_values = np.linspace(0, 1, num_arrows + 1)
                
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
                    
                    # Scale arrow height based on equivalent water depth
                    if max_load > 0:
                        water_depth = n / gamma_w
                        arrow_height = water_depth  # Direct water depth, not scaled relative to max
                    else:
                        arrow_height = 0
                    
                    # For very small arrows, just store surface point for connecting line
                    if arrow_height < 0.5:
                        top_xs.append(x)
                        top_ys.append(y)
                        continue
            
                    
                    # Calculate arrow start point (above surface)
                    arrow_start_x = x + perp_dx * arrow_height
                    arrow_start_y = y + perp_dy * arrow_height
                    
                    # Store points for connecting line
                    top_xs.append(arrow_start_x)
                    top_ys.append(arrow_start_y)
                    
                    # Draw arrow - extend all the way to surface point
                    arrow_length = np.sqrt((x - arrow_start_x)**2 + (y - arrow_start_y)**2)
                    if head_length > arrow_length:
                        # Draw a simple line without arrowhead
                        ax.plot([arrow_start_x, x], [arrow_start_y, y], 
                               color=color, linewidth=2, alpha=0.7)
                    else:
                        # Draw arrow with head
                        ax.arrow(arrow_start_x, arrow_start_y, 
                                x - arrow_start_x, y - arrow_start_y,
                                head_width=head_width, head_length=head_length, 
                                fc=color, ec=color, alpha=0.7,
                                length_includes_head=True)
                
                # Add end point if it's the last segment and load is zero
                if i == len(line) - 2 and n2 == 0:
                    top_xs.append(x2)
                    top_ys.append(y2)
                
                # Draw connecting line at arrow tops
                if top_xs:
                    ax.plot(top_xs, top_ys, color=color, linewidth=1.5, alpha=0.8)
            
            # Draw the surface line itself
            ax.plot(xs, ys, color=color, linewidth=1.5, alpha=0.8, label=label)
    
    dloads = slope_data['dloads']
    dloads2 = slope_data.get('dloads2', [])
    plot_single_dload_set(ax, dloads, 'purple', 'Distributed Load')
    plot_single_dload_set(ax, dloads2, 'orange', 'Distributed Load 2')

def plot_circles(ax, slope_data):
    """
    Plots starting circles with center markers and arrows.

    Parameters:
        ax (matplotlib axis): The plotting axis
        slope_data (dict): Slope data dictionary containing circles

    Returns:
        None
    """
    circles = slope_data['circles']
    for circle in circles:
        Xo = circle['Xo']
        Yo = circle['Yo']
        R = circle['R']
        # theta = np.linspace(0, 2 * np.pi, 100)
        # x_circle = Xo + R * np.cos(theta)
        # y_circle = Yo + R * np.sin(theta)
        # ax.plot(x_circle, y_circle, 'r--', label='Circle')

        # Plot the portion of the circle in the slope
        ground_surface = slope_data['ground_surface']
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
    if not non_circ or len(non_circ) == 0:
        return
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

    # Check if any materials have non-zero d and psi values
    has_d_psi = any(mat.get('d', 0) > 0 or mat.get('psi', 0) > 0 for mat in materials)

    # Check material options
    options = set(mat['option'] for mat in materials)

    # Decide column headers
    if options == {'mc'}:
        if has_d_psi:
            col_labels = ["Mat", "Name", "γ", "c", "φ", "d", "ψ"]
        else:
            col_labels = ["Mat", "Name", "γ", "c", "φ"]
    elif options == {'cp'}:
        if has_d_psi:
            col_labels = ["Mat", "Name", "γ", "cp", "rₑ", "d", "ψ"]
        else:
            col_labels = ["Mat", "Name", "γ", "cp", "rₑ"]
    else:
        if has_d_psi:
            col_labels = ["Mat", "Name", "γ", "c / cp", "φ / rₑ", "d", "ψ"]
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
            if has_d_psi:
                d = mat.get('d', 0)
                psi = mat.get('psi', 0)
                d_str = f"{d:.1f}" if d > 0 or psi > 0 else "-"
                psi_str = f"{psi:.1f}" if d > 0 or psi > 0 else "-"
                row = [idx+1, name, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}", d_str, psi_str]
            else:
                row = [idx+1, name, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}"]
        elif option == 'cp':
            cp = mat['cp']
            r_elev = mat['r_elev']
            if has_d_psi:
                d = mat.get('d', 0)
                psi = mat.get('psi', 0)
                d_str = f"{d:.1f}" if d > 0 or psi > 0 else "-"
                psi_str = f"{psi:.1f}" if d > 0 or psi > 0 else "-"
                row = [idx+1, name, f"{gamma:.1f}", f"{cp:.2f}", f"{r_elev:.1f}", d_str, psi_str]
            else:
                row = [idx+1, name, f"{gamma:.1f}", f"{cp:.2f}", f"{r_elev:.1f}"]
        else:
            if has_d_psi:
                d = mat.get('d', 0)
                psi = mat.get('psi', 0)
                d_str = f"{d:.1f}" if d > 0 or psi > 0 else "-"
                psi_str = f"{psi:.1f}" if d > 0 or psi > 0 else "-"
                row = [idx+1, name, f"{gamma:.1f}", "-", "-", d_str, psi_str]
            else:
                row = [idx+1, name, f"{gamma:.1f}", "-", "-"]
        table_data.append(row)

    # Adjust table width based on number of columns
    table_width = 0.25 if has_d_psi else 0.2

    # Add the table
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, table_width, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(8)

def plot_base_stresses(ax, slice_df, scale_frac=0.5, alpha=0.3):
    """
    Plots base normal stresses for each slice as bars.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data
        scale_frac: Fraction of plot height for bar scaling
        alpha: Transparency for bars

    Returns:
        None
    """
    u = slice_df['u'].values
    n_eff = slice_df['n_eff'].values
    dl = slice_df['dl'].values
    heights = slice_df['y_ct'] - slice_df['y_cb']
    max_ht = heights.max() if not heights.empty else 1.0
    max_bar_len = max_ht * scale_frac

    max_stress = np.max(np.abs(n_eff)) if len(n_eff) > 0 else 1.0
    max_u = np.max(u) if len(u) > 0 else 1.0

    for i, (index, row) in enumerate(slice_df.iterrows()):
        if i >= len(n_eff):
            break

        x1, y1 = row['x_l'], row['y_lb']
        x2, y2 = row['x_r'], row['y_rb']

        stress = n_eff[i]
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


def plot_thrust_line_from_df(ax, slice_df,
                            color: str = 'red',
                            linestyle: str = '--',
                            linewidth: float = 1,
                            label: str = 'Line of Thrust'):
    """
    Plots the line of thrust from the slice dataframe.

    Parameters:
        ax: matplotlib Axes object
        slice_df: DataFrame containing slice data with 'yt_l' and 'yt_r' columns
        color: Color of the line
        linestyle: Style of the line
        linewidth: Width of the line
        label: Label for the line in the legend

    Returns:
        None
    """
    # Check if required columns exist
    if 'yt_l' not in slice_df.columns or 'yt_r' not in slice_df.columns:
        return
    
    # Create thrust line coordinates from slice data
    thrust_xs = []
    thrust_ys = []
    
    for _, row in slice_df.iterrows():
        # Add left point of current slice
        thrust_xs.append(row['x_l'])
        thrust_ys.append(row['yt_l'])
        
        # Add right point of current slice (same as left point of next slice)
        thrust_xs.append(row['x_r'])
        thrust_ys.append(row['yt_r'])
    
    # Plot the thrust line
    ax.plot(thrust_xs, thrust_ys,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label)

def compute_ylim(data, slice_df, scale_frac=0.5, pad_fraction=0.1):
    """
    Computes y-limits for plotting based on slice data.

    Parameters:
        data: Input data
        slice_df: pandas.DataFrame with slice data, must have 'y_lt' and 'y_lb' for stress‐bar sizing
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
    heights = slice_df["y_lt"] - slice_df["y_lb"]
    if not heights.empty:
        max_bar = heights.max() * scale_frac
        y_min -= max_bar
        y_max += max_bar

    # 4) add a final small pad
    pad = (y_max - y_min) * pad_fraction
    return y_min - pad, y_max + pad

# ========== FOR PLOTTING INPUT DATA  =========

def plot_reinforcement_lines(ax, slope_data):
    """
    Plots the reinforcement lines from slope_data.
    
    Parameters:
        ax: matplotlib Axes object
        slope_data: Dictionary containing slope data with 'reinforce_lines' key
        
    Returns:
        None
    """
    if 'reinforce_lines' not in slope_data or not slope_data['reinforce_lines']:
        return
        
    tension_points_plotted = False  # Track if tension points have been added to legend
    
    for i, line in enumerate(slope_data['reinforce_lines']):
        # Extract x and y coordinates from the line points
        xs = [point['X'] for point in line]
        ys = [point['Y'] for point in line]
        
        # Plot the reinforcement line with a distinctive style
        ax.plot(xs, ys, color='darkgray', linewidth=3, linestyle='-', 
                alpha=0.8, label='Reinforcement Line' if i == 0 else "")
        
        # Add markers at each point to show tension values
        for j, point in enumerate(line):
            tension = point.get('T', 0.0)
            if tension > 0:
                # Use smaller marker size proportional to tension (normalized)
                max_tension = max(p.get('T', 0.0) for p in line)
                marker_size = 10 + 15 * (tension / max_tension) if max_tension > 0 else 10
                ax.scatter(point['X'], point['Y'], s=marker_size, 
                          color='red', alpha=0.7, zorder=5,
                          label='Tension Points' if not tension_points_plotted else "")
                tension_points_plotted = True


def plot_inputs(slope_data, title="Slope Geometry and Inputs", width=12, height=6, mat_table=True):
    """
    Creates a plot showing the slope geometry and input parameters.

    Parameters:
        slope_data: Dictionary containing plot data
        title: Title for the plot
        width: Width of the plot in inches
        height: Height of the plot in inches
        mat_table: Controls material table display. Can be:
            - True: Auto-position material table to avoid overlaps
            - False: Don't show material table
            - 'auto': Auto-position material table to avoid overlaps
            - String: Specific location for material table ('upper left', 'upper right', 'upper center',
                     'lower left', 'lower right', 'lower center', 'center left', 'center right', 'center')

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))

    # Plot contents
    plot_profile_lines(ax, slope_data['profile_lines'])
    plot_max_depth(ax, slope_data['profile_lines'], slope_data['max_depth'])
    plot_piezo_line(ax, slope_data)
    plot_dloads(ax, slope_data)
    plot_tcrack_surface(ax, slope_data['tcrack_surface'])
    plot_reinforcement_lines(ax, slope_data)

    if slope_data['circular']:
        plot_circles(ax, slope_data)
    elif slope_data.get('non_circ') and len(slope_data['non_circ']) > 0:
        plot_non_circ(ax, slope_data['non_circ'])

    # Handle material table display
    if mat_table:
        if isinstance(mat_table, str) and mat_table != 'auto':
            # Convert location string to xloc, yloc coordinates (inside plot area with margins)
            location_map = {
                'upper left': (0.05, 0.70),
                'upper right': (0.70, 0.70),
                'upper center': (0.35, 0.70),
                'lower left': (0.05, 0.05),
                'lower right': (0.70, 0.05),
                'lower center': (0.35, 0.05),
                'center left': (0.05, 0.35),
                'center right': (0.70, 0.35),
                'center': (0.35, 0.35)
            }
            if mat_table in location_map:
                xloc, yloc = location_map[mat_table]
                plot_material_table(ax, slope_data['materials'], xloc=xloc, yloc=yloc)
            else:
                # Default to upper right if invalid location
                plot_material_table(ax, slope_data['materials'], xloc=0.75, yloc=0.75)
        else:
            # Auto-position or default: find best location
            plot_elements_bounds = get_plot_elements_bounds(ax, slope_data)
            xloc, yloc = find_best_table_position(ax, slope_data['materials'], plot_elements_bounds)
            plot_material_table(ax, slope_data['materials'], xloc=xloc, yloc=yloc)

    ax.set_aspect('equal')  # ✅ Equal aspect
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    
    # Add distributed load to legend if present
    if slope_data['dloads']:
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

def plot_solution(slope_data, slice_df, failure_surface, results, width=12, height=7, slice_numbers=False):
    """
    Plots the full solution including slices, numbers, thrust line, and base stresses.

    Parameters:
        data: Input data
        slice_df: DataFrame containing slice data
        failure_surface: Failure surface geometry
        results: Solution results
        width: Width of the plot in inches
        height: Height of the plot in inches

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)

    plot_profile_lines(ax, slope_data['profile_lines'])
    plot_max_depth(ax, slope_data['profile_lines'], slope_data['max_depth'])
    plot_slices(ax, slice_df, fill=False)
    plot_failure_surface(ax, failure_surface)
    plot_piezo_line(ax, slope_data)
    plot_dloads(ax, slope_data)
    plot_tcrack_surface(ax, slope_data['tcrack_surface'])
    plot_reinforcement_lines(ax, slope_data)
    if slice_numbers:
        plot_slice_numbers(ax, slice_df)
    # plot_material_table(ax, data['materials'], xloc=0.75) # Adjust this so that it fits with the legend

    alpha = 0.3
    if results['method'] == 'spencer':
        plot_thrust_line_from_df(ax, slice_df)

    plot_base_stresses(ax, slice_df, alpha=alpha)

    import matplotlib.patches as mpatches
    normal_patch = mpatches.Patch(facecolor='none', edgecolor='green', hatch='.....', label="Eff Normal Stress (σ')")
    pore_patch = mpatches.Patch(color='blue', alpha=alpha, label='Pore Pressure (u)')

    # Get legend handles and labels
    handles, labels = ax.get_legend_handles_labels()
    handles.extend([normal_patch, pore_patch])
    labels.extend(["Eff Normal Stress (σ')", 'Pore Pressure (u)'])
    
    # Add distributed load to legend if present
    if slope_data['dloads']:
        handler_class, dummy_line = get_dload_legend_handler()
        handles.append(dummy_line)
        labels.append('Distributed Load')
    
    ax.legend(
        handles=handles,
        labels=labels,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),
        ncol=3
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
    ymin, ymax = compute_ylim(slope_data, slice_df, pad_fraction=0.05)
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

def plot_circular_search_results(slope_data, fs_cache, search_path=None, highlight_fs=True, width=12, height=7):
    """
    Creates a plot showing the results of a circular failure surface search.

    Parameters:
        slope_data: Dictionary containing plot data
        fs_cache: List of dictionaries containing failure surface data and FS values
        search_path: List of dictionaries containing search path coordinates
        highlight_fs: Boolean indicating whether to highlight the critical failure surface
        width: Width of the plot in inches
        height: Height of the plot in inches

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))

    plot_profile_lines(ax, slope_data['profile_lines'])
    plot_max_depth(ax, slope_data['profile_lines'], slope_data['max_depth'])
    plot_piezo_line(ax, slope_data)
    plot_dloads(ax, slope_data)
    plot_tcrack_surface(ax, slope_data['tcrack_surface'])

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

def plot_noncircular_search_results(slope_data, fs_cache, search_path=None, highlight_fs=True, width=12, height=7):
    """
    Creates a plot showing the results of a non-circular failure surface search.

    Parameters:
        slope_data: Dictionary containing plot data
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
    plot_profile_lines(ax, slope_data['profile_lines'])
    plot_max_depth(ax, slope_data['profile_lines'], slope_data['max_depth'])
    plot_piezo_line(ax, slope_data)
    plot_dloads(ax, slope_data)
    plot_tcrack_surface(ax, slope_data['tcrack_surface'])

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

def plot_reliability_results(slope_data, reliability_data, width=12, height=7):
    """
    Creates a plot showing the results of reliability analysis.
    
    Parameters:
        slope_data: Dictionary containing plot data
        reliability_data: Dictionary containing reliability analysis results
        width: Width of the plot in inches
        height: Height of the plot in inches
    
    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=(width, height))

    # Plot basic slope elements (same as other search functions)
    plot_profile_lines(ax, slope_data['profile_lines'])
    plot_max_depth(ax, slope_data['profile_lines'], slope_data['max_depth'])
    plot_piezo_line(ax, slope_data)
    plot_dloads(ax, slope_data)
    plot_tcrack_surface(ax, slope_data['tcrack_surface'])

    # Plot reliability-specific failure surfaces
    fs_cache = reliability_data['fs_cache']
    
    # Plot all failure surfaces
    for i, fs_data in enumerate(fs_cache):
        result = fs_data['result']
        name = fs_data['name']
        failure_surface = result['failure_surface']
        
        # Convert failure surface to coordinates
        if hasattr(failure_surface, 'coords'):
            coords = list(failure_surface.coords)
        else:
            coords = failure_surface
        
        x_coords = [pt[0] for pt in coords]
        y_coords = [pt[1] for pt in coords]
        
        # Color and styling based on surface type
        if name == "MLV":
            # Highlight the MLV (critical) surface in red
            ax.plot(x_coords, y_coords, color='red', linewidth=3, 
                   label=f'$F_{{MLV}}$ Surface (FS={result["FS"]:.3f})', zorder=10)
        else:
            # Other surfaces in different colors
            if '+' in name:
                color = 'blue'
                alpha = 0.7
                label = f'$F^+$ ({name}) (FS={result["FS"]:.3f})'
            else:  # '-' in name
                color = 'green'
                alpha = 0.7
                label = f'$F^-$ ({name}) (FS={result["FS"]:.3f})'
            
            ax.plot(x_coords, y_coords, color=color, linewidth=1.5, 
                   alpha=alpha, label=label, zorder=5)



    # Standard finalization
    ax.set_aspect('equal')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Title with reliability statistics using mathtext
    F_MLV = reliability_data['F_MLV']
    sigma_F = reliability_data['sigma_F']
    COV_F = reliability_data['COV_F']
    reliability = reliability_data['reliability']
    prob_failure = reliability_data['prob_failure']
    
    ax.set_title(f"Reliability Analysis Results\n"
                f"$F_{{MLV}}$ = {F_MLV:.3f}, $\\sigma_F$ = {sigma_F:.3f}, "
                f"$COV_F$ = {COV_F:.3f}\n"
                f"Reliability = {reliability*100:.2f}%, $P_f$ = {prob_failure*100:.2f}%")

    plt.tight_layout()
    plt.show()

def plot_mesh(mesh, materials=None, figsize=(14, 6), pad_frac=0.05, show_nodes=True, label_elements=False, label_nodes=False):
    """
    Plot the finite element mesh with material regions.
    
    Parameters:
        mesh: Mesh dictionary with 'nodes', 'elements', 'element_types', and 'element_materials' keys
        materials: Optional list of material dictionaries for legend labels
        figsize: Figure size tuple
        pad_frac: Fraction of mesh size to use for padding around plot
        show_nodes: If True, plot points at node locations
        label_elements: If True, label each element with its number at its centroid
        label_nodes: If True, label each node with its number
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.collections import PolyCollection
    import numpy as np
    
    nodes = mesh["nodes"]
    elements = mesh["elements"]
    element_types = mesh["element_types"]
    mat_ids = mesh["element_materials"]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Group elements by material ID
    material_elements = {}
    for i, (element, elem_type, mid) in enumerate(zip(elements, element_types, mat_ids)):
        if mid not in material_elements:
            material_elements[mid] = []
        
        # Only process 2D elements (skip 1D elements which have elem_type 2)
        if elem_type == 2:  # Skip 1D elements
            continue
        
        # Use corner nodes to define element boundary (no subdivision needed)
        if elem_type in [3, 6]:  # Triangular elements (linear or quadratic)
            element_coords = [nodes[element[0]], nodes[element[1]], nodes[element[2]]]
        elif elem_type in [4, 8, 9]:  # Quadrilateral elements (linear or quadratic)
            element_coords = [nodes[element[0]], nodes[element[1]], nodes[element[2]], nodes[element[3]]]
        else:
            continue  # Skip unknown element types
        
        material_elements[mid].append(element_coords)
    
    legend_elements = []
    
    # Plot 1D elements FIRST (bottom layer) if present in mesh
    if "elements_1d" in mesh and "element_types_1d" in mesh and "element_materials_1d" in mesh:
        elements_1d = mesh["elements_1d"]
        element_types_1d = mesh["element_types_1d"]
        mat_ids_1d = mesh["element_materials_1d"]
        
        # Group 1D elements by material ID
        material_lines = {}
        for i, (element_1d, elem_type_1d, mid_1d) in enumerate(zip(elements_1d, element_types_1d, mat_ids_1d)):
            if mid_1d not in material_lines:
                material_lines[mid_1d] = []
            
            # Get line coordinates based on actual number of nodes
            # elem_type_1d contains the number of nodes (2 for linear, 3 for quadratic)
            if elem_type_1d == 2:  # Linear 1D element (2 nodes)
                # Skip zero-padded elements
                if element_1d[1] != 0:  # Valid second node
                    line_coords = [nodes[element_1d[0]], nodes[element_1d[1]]]
                else:
                    continue  # Skip invalid element
            elif elem_type_1d == 3:  # Quadratic 1D element (3 nodes)
                # For visualization, connect all three nodes or just endpoints
                line_coords = [nodes[element_1d[0]], nodes[element_1d[1]], nodes[element_1d[2]]]
            else:
                continue  # Skip unknown 1D element types
            
            material_lines[mid_1d].append(line_coords)
        
        # Plot 1D elements with distinctive style
        for mid_1d, lines_list in material_lines.items():
            for line_coords in lines_list:
                xs = [coord[0] for coord in line_coords]
                ys = [coord[1] for coord in line_coords]
                ax.plot(xs, ys, color='red', linewidth=3, alpha=0.8, solid_capstyle='round')
        
        # Add 1D elements to legend
        if material_lines:
            legend_elements.append(plt.Line2D([0], [0], color='red', linewidth=3, 
                                            alpha=0.8, label='1D Elements'))
    
    # Plot 2D elements SECOND (middle layer)
    for mid, elements_list in material_elements.items():
        # Create polygon collection for this material
        poly_collection = PolyCollection(elements_list, 
                                       facecolor=get_material_color(mid),
                                       edgecolor='k',
                                       alpha=0.4,
                                       linewidth=0.5)
        ax.add_collection(poly_collection)
        
        # Add to legend
        if materials and mid <= len(materials) and materials[mid-1].get('name'):
            label = materials[mid-1]['name']  # Convert to 0-based indexing
        else:
            label = f'Material {mid}'
        
        legend_elements.append(Patch(facecolor=get_material_color(mid), 
                                   edgecolor='k', 
                                   alpha=0.4, 
                                   label=label))
    
    # Label 2D elements if requested
    if label_elements:
        for idx, (element, element_type) in enumerate(zip(elements, element_types)):
            # Calculate element centroid based on element type
            if element_type == 3:  # 3-node triangle
                element_coords = nodes[element[:3]]
            elif element_type == 6:  # 6-node triangle
                element_coords = nodes[element[:6]]
            elif element_type == 4:  # 4-node quad
                element_coords = nodes[element[:4]]
            elif element_type == 8:  # 8-node quad
                element_coords = nodes[element[:8]]
            elif element_type == 9:  # 9-node quad
                element_coords = nodes[element[:9]]
            else:
                continue  # Skip unknown element types
            
            centroid = np.mean(element_coords, axis=0)
            ax.text(centroid[0], centroid[1], str(idx+1),
                    ha='center', va='center', fontsize=6, color='black', alpha=0.7,
                    zorder=12)
    
    # Label 1D elements if requested (with different color)
    if label_elements and "elements_1d" in mesh:
        elements_1d = mesh["elements_1d"]
        element_types_1d = mesh["element_types_1d"]
        
        for idx, (element_1d, elem_type_1d) in enumerate(zip(elements_1d, element_types_1d)):
            # Skip zero-padded elements
            if elem_type_1d == 2 and element_1d[1] != 0:  # Linear 1D element
                # Calculate midpoint of line element
                coord1 = nodes[element_1d[0]]
                coord2 = nodes[element_1d[1]]
                midpoint = (coord1 + coord2) / 2
                ax.text(midpoint[0], midpoint[1], f"1D{idx+1}",
                        ha='center', va='center', fontsize=6, color='black', alpha=0.9,
                        zorder=13)
            elif elem_type_1d == 3 and element_1d[2] != 0:  # Quadratic 1D element
                # Use middle node as label position (if it exists)
                midpoint = nodes[element_1d[1]]
                ax.text(midpoint[0], midpoint[1], f"1D{idx+1}",
                        ha='center', va='center', fontsize=6, color='black', alpha=0.9,
                        zorder=13)
    
    # Plot nodes LAST (top layer) if requested
    if show_nodes:
        # Plot all nodes - if meshing is correct, all nodes should be used
        ax.plot(nodes[:, 0], nodes[:, 1], 'k.', markersize=2)
        # Add to legend
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                        markerfacecolor='k', markersize=6, 
                                        label=f'Nodes ({len(nodes)})', linestyle='None'))
    
    # Label nodes if requested
    if label_nodes:
        # Label all nodes
        for i, (x, y) in enumerate(nodes):
            ax.text(x + 0.5, y + 0.5, str(i+1), fontsize=6, color='blue', alpha=0.7,
                    ha='left', va='bottom', zorder=14)
    
    ax.set_aspect('equal')
    ax.set_title("Finite Element Mesh with Material Regions (Triangles and Quads)")
    
    # Add legend if we have materials
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=min(len(legend_elements), 4))

    # Add cushion
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    y_min, y_max = nodes[:, 1].min(), nodes[:, 1].max()
    x_pad = (x_max - x_min) * pad_frac
    y_pad = (y_max - y_min) * pad_frac
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    
    # Add extra cushion for legend space
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    plt.tight_layout()
    plt.show()


def plot_polygons(polygons, title="Material Zone Polygons"):
    """
    Plot all material zone polygons in a single figure.
    
    Parameters:
        polygons: List of polygon coordinate lists
        title: Plot title
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(12, 8))
    for i, polygon in enumerate(polygons):
        xs = [x for x, y in polygon]
        ys = [y for x, y in polygon]
        ax.fill(xs, ys, color=get_material_color(i), alpha=0.6, label=f'Material {i}')
        ax.plot(xs, ys, color=get_material_color(i), linewidth=1)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()


def plot_polygons_separately(polygons, title_prefix='Material Zone'):
    """
    Plot each polygon in a separate matplotlib frame (subplot), with vertices as round dots.
    
    Parameters:
        polygons: List of polygon coordinate lists
        title_prefix: Prefix for each subplot title
    """
    import matplotlib.pyplot as plt
    
    n = len(polygons)
    fig, axes = plt.subplots(n, 1, figsize=(8, 3 * n), squeeze=False)
    for i, polygon in enumerate(polygons):
        xs = [x for x, y in polygon]
        ys = [y for x, y in polygon]
        ax = axes[i, 0]
        ax.fill(xs, ys, color=get_material_color(i), alpha=0.6, label=f'Material {i}')
        ax.plot(xs, ys, color=get_material_color(i), linewidth=1)
        ax.scatter(xs, ys, color='k', s=30, marker='o', zorder=3, label='Vertices')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_title(f'{title_prefix} {i}')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        ax.legend()
    plt.tight_layout()
    plt.show()


def find_best_table_position(ax, materials, plot_elements_bounds):
    """
    Find the best position for the material table to avoid overlaps.
    
    Parameters:
        ax: matplotlib Axes object
        materials: List of materials to determine table size
        plot_elements_bounds: List of (x_min, x_max, y_min, y_max) for existing elements
        
    Returns:
        (xloc, yloc) coordinates for table placement
    """
    # Calculate table size based on number of materials and columns
    num_materials = len(materials)
    has_d_psi = any(mat.get('d', 0) > 0 or mat.get('psi', 0) > 0 for mat in materials)
    table_height = 0.05 + 0.025 * num_materials  # Height per row
    table_width = 0.25 if has_d_psi else 0.2
    
    # Define candidate positions (priority order) - with margins from borders
    candidates = [
        (0.05, 0.70),  # upper left
        (0.70, 0.70),  # upper right  
        (0.05, 0.05),  # lower left
        (0.70, 0.05),  # lower right
        (0.35, 0.70),  # upper center
        (0.35, 0.05),  # lower center
        (0.05, 0.35),  # center left
        (0.70, 0.35),  # center right
        (0.35, 0.35),  # center
    ]
    
    # Check each candidate position for overlaps
    for xloc, yloc in candidates:
        table_bounds = (xloc, xloc + table_width, yloc - table_height, yloc)
        
        # Check if table overlaps with any plot elements
        overlap = False
        for elem_bounds in plot_elements_bounds:
            elem_x_min, elem_x_max, elem_y_min, elem_y_max = elem_bounds
            table_x_min, table_x_max, table_y_min, table_y_max = table_bounds
            
            # Check for overlap
            if not (table_x_max < elem_x_min or table_x_min > elem_x_max or
                   table_y_max < elem_y_min or table_y_min > elem_y_max):
                overlap = True
                break
        
        if not overlap:
            return xloc, yloc
    
    # If all positions have overlap, return the first candidate
    return candidates[0]


def get_plot_elements_bounds(ax, slope_data):
    """
    Get bounding boxes of existing plot elements to avoid overlaps.
    
    Parameters:
        ax: matplotlib Axes object
        slope_data: Dictionary containing slope data
        
    Returns:
        List of (x_min, x_max, y_min, y_max) tuples for plot elements
    """
    bounds = []
    
    # Get axis limits
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    
    # Profile lines bounds
    if 'profile_lines' in slope_data:
        for line in slope_data['profile_lines']:
            if line:
                xs = [p[0] for p in line]
                ys = [p[1] for p in line]
                bounds.append((min(xs), max(xs), min(ys), max(ys)))
    
    # Distributed loads bounds
    if 'dloads' in slope_data and slope_data['dloads']:
        for dload_set in slope_data['dloads']:
            if dload_set:
                xs = [p['X'] for p in dload_set]
                ys = [p['Y'] for p in dload_set]
                bounds.append((min(xs), max(xs), min(ys), max(ys)))
    
    # Reinforcement lines bounds
    if 'reinforce_lines' in slope_data and slope_data['reinforce_lines']:
        for line in slope_data['reinforce_lines']:
            if line:
                xs = [p['X'] for p in line]
                ys = [p['Y'] for p in line]
                bounds.append((min(xs), max(xs), min(ys), max(ys)))
    
    return bounds