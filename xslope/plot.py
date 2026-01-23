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
import math

from .slice import generate_failure_surface

# Configure matplotlib for better text rendering
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 10
})

# Consistent color for materials (Tableau tab10)
def get_material_color(idx):
    cmap = plt.get_cmap("tab10")
    # Use the colormap callable rather than relying on the (typing-unknown) `.colors` attribute.
    return cmap(idx % getattr(cmap, "N", 10))

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


def plot_profile_lines(ax, profile_lines, materials=None, labels=False):
    """
    Plots the profile lines for each material in the slope.

    Parameters:
        ax: matplotlib Axes object
        profile_lines: List of profile line dicts, each with 'coords' and 'mat_id' keys
        materials: List of material dictionaries (optional, for color mapping)
        labels: If True, add index labels to each profile line (default: False)

    Returns:
        None
    """
    for i, line in enumerate(profile_lines):
        coords = line['coords']
        xs, ys = zip(*coords)
        
        # Get material index from mat_id (already 0-based)
        if materials and line.get('mat_id') is not None:
            mat_idx = line['mat_id']
            if 0 <= mat_idx < len(materials):
                color = get_material_color(mat_idx)
            else:
                # Fallback to index-based color if mat_id out of range
                color = get_material_color(i)
        else:
            # Fallback to index-based color if no materials or mat_id
            color = get_material_color(i)
        
        ax.plot(xs, ys, color=color, linewidth=1, label=f'Profile {i+1}')

        if labels:
            _add_profile_index_label(ax, coords, i + 1, color)


def _add_profile_index_label(ax, line, index, color):
    """
    Adds an index label to a profile line, positioned on a suitable segment.

    Parameters:
        ax: matplotlib Axes object
        line: List of (x, y) coordinates for the profile line
        index: The index number to display (1-based)
        color: Color for the label text

    Returns:
        None
    """
    if len(line) < 2:
        return

    # Build list of segments with their properties
    segments = []
    for j in range(len(line) - 1):
        x1, y1 = line[j]
        x2, y2 = line[j + 1]
        dx = x2 - x1
        dy = y2 - y1
        length = np.sqrt(dx**2 + dy**2)

        if length < 1e-9:
            continue

        # Calculate how horizontal the segment is (1.0 = perfectly horizontal)
        horizontalness = abs(dx) / length

        # Midpoint of the segment
        mid_x = (x1 + x2) / 2
        mid_y = (y1 + y2) / 2

        segments.append({
            'length': length,
            'horizontalness': horizontalness,
            'mid_x': mid_x,
            'mid_y': mid_y,
            'position': j  # segment index in the line
        })

    if not segments:
        return

    # Calculate the total line length and find segments in the middle third
    total_length = sum(s['length'] for s in segments)

    # Score segments: prefer longer, more horizontal segments near the middle
    # Avoid very short segments (less than 5% of total length)
    min_length_threshold = 0.05 * total_length

    n_segments = len(segments)
    best_segment = None
    best_score = -1

    for seg in segments:
        # Skip very short segments
        if seg['length'] < min_length_threshold:
            continue

        # Calculate how close to the middle this segment is (0-1 scale, 1 = center)
        position_ratio = (seg['position'] + 0.5) / n_segments
        middle_score = 1.0 - 2.0 * abs(position_ratio - 0.5)  # 1.0 at center, 0.0 at ends

        # Score: weight length, horizontalness, and middle position
        # Length is normalized by total length
        length_score = seg['length'] / total_length

        # Combined score: prioritize horizontalness, then middle position, then length
        score = (seg['horizontalness'] * 2.0 +
                 middle_score * 1.5 +
                 length_score * 1.0)

        if score > best_score:
            best_score = score
            best_segment = seg

    # Fallback: if no segment passed the threshold, use the longest one
    if best_segment is None:
        best_segment = max(segments, key=lambda s: s['length'])

    # Place the label at the midpoint of the chosen segment
    ax.text(
        best_segment['mid_x'],
        best_segment['mid_y'],
        str(index),
        fontsize=7,
        color=color,
        fontfamily='monospace',
        ha='center',
        va='center',
        bbox=dict(
            boxstyle='round,pad=0.3',
            facecolor='white',
            edgecolor=color,
            linewidth=0.8
        ),
        zorder=10
    )

def plot_max_depth(ax, profile_lines, max_depth):
    """
    Plots a horizontal line representing the maximum depth limit with hash marks.

    Parameters:
        ax: matplotlib Axes object
        profile_lines: List of profile line dicts, each with 'coords' key containing coordinate tuples
        max_depth: Maximum allowed depth for analysis

    Returns:
        None
    """
    if max_depth is None:
        return
    x_vals = [x for line in profile_lines for x, _ in line['coords']]
    x_min = min(x_vals)
    x_max = max(x_vals)
    ax.hlines(max_depth, x_min, x_max, colors='black', linewidth=1.5, label='Max Depth')

    x_diff = x_max - x_min
    spacing = x_diff / 100
    length = x_diff / 80

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
    
    def _plot_touching_v_marker(ax, x, y, color, markersize=8, extra_gap_points=0.0):
        """
        Place an inverted triangle marker so its tip visually touches the line at (x, y).
        We do this in display coordinates (points/pixels) so it scales consistently.
        """
        from matplotlib.markers import MarkerStyle
        from matplotlib.transforms import offset_copy

        # Compute the distance (in "marker units") from the marker origin to the tip.
        # Matplotlib scales marker vertices by `markersize` (in points) for Line2D.
        ms = MarkerStyle("v")
        path = ms.get_path().transformed(ms.get_transform())
        verts = np.asarray(path.vertices)
        min_y = float(verts[:, 1].min())  # tip is the lowest y
        tip_offset_points = (-min_y) * float(markersize) + float(extra_gap_points)

        # Offset the marker center upward in point units so the tip lands at (x, y).
        trans = offset_copy(ax.transData, fig=ax.figure, x=0.0, y=tip_offset_points, units="points")
        ax.plot([x], [y], marker="v", color=color, markersize=markersize, linestyle="None", transform=trans)

    def plot_single_piezo_line(ax, piezo_line, color, label):
        """Internal function to plot a single piezometric line"""
        if not piezo_line:
            return
            
        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, color=color, linewidth=2, label=label)
        
        # Find middle x-coordinate and corresponding y value
        if len(piezo_xs) > 1:
            # Sort by x to ensure monotonic input for interpolation
            pairs = sorted(zip(piezo_xs, piezo_ys), key=lambda p: p[0])
            sx, sy = zip(*pairs)
            x_min, x_max = min(sx), max(sx)
            mid_x = (x_min + x_max) / 2
            mid_y = float(np.interp(mid_x, sx, sy))
            # Slight negative gap so the marker visually "touches" the line (not floating above it)
            _plot_touching_v_marker(ax, mid_x, mid_y, color=color, markersize=8, extra_gap_points=2.0)
    
    # Plot both piezometric lines
    plot_single_piezo_line(ax, slope_data.get('piezo_line'), 'b', "Piezometric Line")
    plot_single_piezo_line(ax, slope_data.get('piezo_line2'), 'skyblue', "Piezometric Line 2")

def plot_seepage_bc_lines(ax, slope_data):
    """
    Plots seep boundary-condition lines for seep-only workflows.

    - Specified head geometry: solid dark blue, thicker than profile lines
    - Exit face geometry: solid red
    - Derived "water level" line for each specified head: y = h
      plotted as a lighter blue solid line with an inverted triangle marker
      (styled similarly to piezometric line markers).
    """
    def _plot_touching_v_marker(ax, x, y, color, markersize=8, extra_gap_points=2.0):
        """Place an inverted triangle so its tip visually sits on the line at (x, y)."""
        from matplotlib.markers import MarkerStyle
        from matplotlib.transforms import offset_copy

        ms = MarkerStyle("v")
        path = ms.get_path().transformed(ms.get_transform())
        verts = np.asarray(path.vertices)
        min_y = float(verts[:, 1].min())
        tip_offset_points = (-min_y) * float(markersize) + float(extra_gap_points)
        trans = offset_copy(ax.transData, fig=ax.figure, x=0.0, y=tip_offset_points, units="points")
        ax.plot([x], [y], marker="v", color=color, markersize=markersize, linestyle="None", transform=trans)

    # Geometry x-extent (used for vertical-head-line derived segment length / side)
    x_vals = []
    for line in slope_data.get("profile_lines", []):
        try:
            xs_line, _ = zip(*line['coords'])
            x_vals.extend(xs_line)
        except Exception:
            pass
    gs = slope_data.get("ground_surface", None)
    if not x_vals and gs is not None and hasattr(gs, "coords"):
        x_vals.extend([x for x, _ in gs.coords])
    x_min_geom = min(x_vals) if x_vals else 0.0
    x_max_geom = max(x_vals) if x_vals else 1.0
    geom_width = max(1e-9, x_max_geom - x_min_geom)

    seepage_bc = slope_data.get("seepage_bc") or {}
    specified_heads = seepage_bc.get("specified_heads") or []
    exit_face = seepage_bc.get("exit_face") or []

    # --- Specified head lines + derived water-level lines ---
    for i, sh in enumerate(specified_heads):
        coords = sh.get("coords") or []
        if len(coords) < 2:
            continue

        xs, ys = zip(*coords)
        ax.plot(
            xs,
            ys,
            color="darkblue",
            linewidth=3,
            linestyle="--",
            label="Specified Head Line" if i == 0 else "",
        )

        # Head values may be scalar (typical) or per-point array-like
        head_val = sh.get("head", None)
        if head_val is None:
            continue

        if isinstance(head_val, (list, tuple, np.ndarray)):
            if len(head_val) != len(coords):
                continue
            heads = [float(h) for h in head_val]
        else:
            try:
                head_scalar = float(head_val)
            except (TypeError, ValueError):
                continue
            heads = [head_scalar] * len(coords)

        # If specified-head geometry is vertical, draw a short horizontal derived line outside the boundary.
        # This avoids drawing a (nearly) vertical derived line that doesn't convey a water level.
        tol = 1e-6
        is_vertical = (max(xs) - min(xs)) <= tol
        if is_vertical:
            x0 = float(xs[0])
            y_head = float(heads[0])
            seg_len = 0.04 * geom_width
            gap = 0.01 * geom_width
            is_right = x0 >= 0.5 * (x_min_geom + x_max_geom)
            if is_right:
                wl_xs = [x0 + gap, x0 + gap + seg_len]
            else:
                wl_xs = [x0 - gap - seg_len, x0 - gap]
            wl_ys = [y_head, y_head]
        else:
            wl_xs = list(xs)
            wl_ys = heads

        ax.plot(
            wl_xs,
            wl_ys,
            color="lightskyblue",
            linewidth=2,
            linestyle="-",
            label="Specified Head Water Level" if i == 0 else "",
        )

        # Inverted triangle marker near the midpoint (like piezometric lines)
        if len(wl_xs) > 1:
            try:
                pairs = sorted(zip(wl_xs, wl_ys), key=lambda p: p[0])
                sx, sy = zip(*pairs)
                mid_x = 0.5 * (min(sx) + max(sx))
                mid_y = float(np.interp(mid_x, sx, sy))
                _plot_touching_v_marker(ax, mid_x, mid_y, color="lightskyblue", markersize=8, extra_gap_points=2.0)
            except Exception:
                # If interpolation fails for any reason, skip marker
                pass

    # --- Exit face line ---
    if len(exit_face) >= 2:
        ex_xs, ex_ys = zip(*exit_face)
        ax.plot(
            ex_xs,
            ex_ys,
            color="red",
            linewidth=3,
            linestyle="--",
            label="Exit Face",
        )

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
        if not isinstance(clipped_surface, LineString):
            clipped_surface = LineString(clipped_surface)
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

        # Draw arrow with pixel-based head size
        ax.annotate('',
                    xy=(Xo + dx * R, Yo + dy * R),  # arrow tip
                    xytext=(Xo, Yo),                 # arrow start
                    arrowprops=dict(
                        arrowstyle='-|>',
                        color='red',
                        lw=1.0,            # shaft width in points
                        mutation_scale=20  # head size in points
                    ))

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

def plot_lem_material_table(ax, materials, xloc=0.6, yloc=0.7):
    """
    Adds a limit equilibrium material properties table to the plot.

    Displays soil properties for limit equilibrium analysis including unit weight (γ),
    cohesion (c), friction angle (φ), and optionally dilation angle (d) and
    dilatancy angle (ψ). Supports both Mohr-Coulomb (mc) and constant-phi (cp) options.

    Parameters:
        ax: matplotlib Axes object to add the table to
        materials: List of material property dictionaries with keys:
            - 'name': Material name (str)
            - 'gamma': Unit weight (float)
            - 'option': Material model - 'mc' or 'cp' (str)
            - 'c': Cohesion for mc option (float)
            - 'phi': Friction angle for mc option (float)
            - 'cp': Constant phi for cp option (float)
            - 'r_elev': Reference elevation for cp option (float)
            - 'd': Dilation angle, optional (float)
            - 'psi': Dilatancy angle, optional (float)
        xloc: x-location of table bottom-left corner in axes coordinates (0-1, default: 0.6)
        yloc: y-location of table bottom-left corner in axes coordinates (0-1, default: 0.7)

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

    # Choose table height based on number of materials (uniform across table types)
    num_rows = max(1, len(materials))
    table_height = 0.06 + 0.035 * num_rows  # header + per-row estimate
    table_height = min(0.35, table_height)  # cap to avoid overflows for many rows

    # Add the table
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, table_width, table_height])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    # Auto layout based on content (shared method for all table types)
    auto_size_table_to_content(ax, table, col_labels, table_data, table_width, table_height)

def plot_seep_material_table(ax, seep_data, xloc=0.6, yloc=0.7):
    """
    Adds a seep material properties table to the plot.

    Displays hydraulic properties for seep analysis including hydraulic conductivities
    (k₁, k₂), anisotropy angle, and unsaturated flow parameters (kr₀, h₀).

    Parameters:
        ax: matplotlib Axes object to add the table to
        seep_data: Dictionary containing seep material properties with keys:
            - 'k1_by_mat': List of primary hydraulic conductivity values (float)
            - 'k2_by_mat': List of secondary hydraulic conductivity values (float)
            - 'angle_by_mat': List of anisotropy angles in degrees (float)
            - 'kr0_by_mat': List of relative permeability at residual saturation (float)
            - 'h0_by_mat': List of pressure head parameters (float)
            - 'material_names': List of material names (str), optional
        xloc: x-location of table bottom-left corner in axes coordinates (0-1, default: 0.6)
        yloc: y-location of table bottom-left corner in axes coordinates (0-1, default: 0.7)

    Returns:
        None
    """
    k1_by_mat = seep_data.get("k1_by_mat")
    k2_by_mat = seep_data.get("k2_by_mat")
    angle_by_mat = seep_data.get("angle_by_mat")
    kr0_by_mat = seep_data.get("kr0_by_mat")
    h0_by_mat = seep_data.get("h0_by_mat")
    material_names = seep_data.get("material_names", [])
    if k1_by_mat is None or len(k1_by_mat) == 0:
        return
    col_labels = ["Mat", "Name", "k₁", "k₂", "Angle", "kr₀", "h₀"]
    table_data = []
    for idx in range(len(k1_by_mat)):
        k1 = k1_by_mat[idx]
        k2 = k2_by_mat[idx] if k2_by_mat is not None else 0.0
        angle = angle_by_mat[idx] if angle_by_mat is not None else 0.0
        kr0 = kr0_by_mat[idx] if kr0_by_mat is not None else 0.0
        h0 = h0_by_mat[idx] if h0_by_mat is not None else 0.0
        material_name = material_names[idx] if idx < len(material_names) else f"Material {idx+1}"
        row = [idx + 1, material_name, f"{k1:.3f}", f"{k2:.3f}", f"{angle:.1f}", f"{kr0:.4f}", f"{h0:.2f}"]
        table_data.append(row)
    # Dimensions
    num_rows = max(1, len(k1_by_mat))
    table_width = 0.45
    table_height = 0.10 + 0.06 * num_rows
    table_height = min(0.50, table_height)
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, table_width, table_height])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    auto_size_table_to_content(ax, table, col_labels, table_data, table_width, table_height)

def plot_fem_material_table(ax, fem_data, xloc=0.6, yloc=0.7, width=0.6, height=None):
    """
    Adds a finite element material properties table to the plot.

    Displays material properties for FEM analysis including unit weight (γ), cohesion (c),
    friction angle (φ), Young's modulus (E), and Poisson's ratio (ν).

    Parameters:
        ax: matplotlib Axes object to add the table to
        fem_data: Dictionary containing FEM material properties with keys:
            - 'c_by_mat': List of cohesion values (float)
            - 'phi_by_mat': List of friction angle values in degrees (float)
            - 'E_by_mat': List of Young's modulus values (float)
            - 'nu_by_mat': List of Poisson's ratio values (float)
            - 'gamma_by_mat': List of unit weight values (float)
            - 'material_names': List of material names (str), optional
        xloc: x-location of table bottom-left corner in axes coordinates (0-1, default: 0.6)
        yloc: y-location of table bottom-left corner in axes coordinates (0-1, default: 0.7)
        width: Table width in axes coordinates (0-1, default: 0.6)
        height: Table height in axes coordinates (0-1, default: auto-calculated)

    Returns:
        None
    """
    c_by_mat = fem_data.get("c_by_mat")
    phi_by_mat = fem_data.get("phi_by_mat")
    E_by_mat = fem_data.get("E_by_mat")
    nu_by_mat = fem_data.get("nu_by_mat")
    gamma_by_mat = fem_data.get("gamma_by_mat")
    material_names = fem_data.get("material_names", [])
    if c_by_mat is None or len(c_by_mat) == 0:
        return
    col_labels = ["Mat", "Name", "γ", "c", "φ", "E", "ν"]
    table_data = []
    for idx in range(len(c_by_mat)):
        c = c_by_mat[idx]
        phi = phi_by_mat[idx] if phi_by_mat is not None else 0.0
        E = E_by_mat[idx] if E_by_mat is not None else 0.0
        nu = nu_by_mat[idx] if nu_by_mat is not None else 0.0
        gamma = gamma_by_mat[idx] if gamma_by_mat is not None else 0.0
        material_name = material_names[idx] if idx < len(material_names) else f"Material {idx+1}"
        row = [idx + 1, material_name, f"{gamma:.1f}", f"{c:.1f}", f"{phi:.1f}", f"{E:.0f}", f"{nu:.2f}"]
        table_data.append(row)
    if height is None:
        num_rows = max(1, len(c_by_mat))
        height = 0.06 + 0.035 * num_rows
        height = min(0.32, height)
    table = ax.table(cellText=table_data,
                     colLabels=col_labels,
                     loc='upper right',
                     colLoc='center',
                     cellLoc='center',
                     bbox=[xloc, yloc, width, height])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    auto_size_table_to_content(ax, table, col_labels, table_data, width, height)
def auto_size_table_to_content(ax, table, col_labels, table_data, table_width, table_height, min_row_frac=0.02, row_pad=1.35, col_min_frac=0.08, col_max_frac=0.15):
    """
    Automatically adjusts table column widths and row heights based on content.

    Measures text extents using the matplotlib renderer and sets column widths proportional
    to content while enforcing minimum and maximum constraints. The "Name" column gets
    more space (18-30%) while numeric columns are constrained to prevent excessive whitespace.
    Row heights are uniform and based on the tallest content in each row.

    Parameters:
        ax: matplotlib Axes object containing the table
        table: matplotlib Table object to be resized
        col_labels: List of column header labels (str)
        table_data: List of lists containing table cell data
        table_width: Total table width in axes coordinates (0-1)
        table_height: Total table height in axes coordinates (0-1)
        min_row_frac: Minimum row height as fraction of axes height (default: 0.02)
        row_pad: Padding factor applied to measured row heights (default: 1.35)
        col_min_frac: Minimum column width as fraction of table width for numeric columns (default: 0.08)
        col_max_frac: Maximum column width as fraction of table width for numeric columns (default: 0.15)

    Returns:
        None

    Notes:
        - The "Name" column is automatically left-aligned and gets 18-30% of table width
        - Numeric columns are center-aligned and constrained to 8-15% of table width
        - All three material table types (LEM, SEEP, FEM) use the same sizing parameters
    """
    # Force draw to get a valid renderer
    try:
        ax.figure.canvas.draw()
        renderer = ax.figure.canvas.get_renderer()
    except Exception:
        renderer = None

    ncols = len(col_labels)
    nrows = len(table_data) + 1  # include header
    # Measure text widths per column in pixels
    widths_px = [1.0] * ncols
    if renderer is not None:
        for c in range(ncols):
            max_w = 1.0
            for r in range(nrows):
                cell = table[(r, c)]
                text = cell.get_text()
                try:
                    bbox = text.get_window_extent(renderer=renderer)
                    max_w = max(max_w, bbox.width)
                except Exception:
                    pass
            widths_px[c] = max_w
    total_w = sum(widths_px) if sum(widths_px) > 0 else float(ncols)
    col_fracs = [w / total_w for w in widths_px]
    # Clamp extreme column widths to keep numeric columns from becoming too wide
    clamped = []
    for i, frac in enumerate(col_fracs):
        label = str(col_labels[i]).lower()
        min_frac = col_min_frac
        max_frac = col_max_frac
        if label == "name":
            min_frac = 0.18
            max_frac = 0.30
        clamped.append(min(max(frac, min_frac), max_frac))
    # Re-normalize to sum to 1.0
    s = sum(clamped)
    if s > 0:
        col_fracs = [c / s for c in clamped]
    # Compute per-row pixel heights based on text extents, convert to axes fraction
    axes_h_px = None
    if renderer is not None:
        try:
            axes_h_px = ax.get_window_extent(renderer=renderer).height
        except Exception:
            axes_h_px = None
    # Fallback axes height if needed (avoid division by zero)
    if not axes_h_px or axes_h_px <= 0:
        axes_h_px = 800.0  # arbitrary but reasonable default
    row_heights_frac = []
    for r in range(nrows):
        max_h_px = 1.0
        if renderer is not None:
            for c in range(ncols):
                try:
                    bbox = table[(r, c)].get_text().get_window_extent(renderer=renderer)
                    max_h_px = max(max_h_px, bbox.height)
                except Exception:
                    pass
        # padding factor to provide breathing room around text
        padded_px = max_h_px * row_pad
        # Convert to axes fraction with minimum clamp
        rh = max(padded_px / axes_h_px, min_row_frac)
        row_heights_frac.append(rh)

    # Apply column widths and per-row heights
    for r in range(nrows):
        for c in range(ncols):
            cell = table[(r, c)]
            cell.set_width(table_width * col_fracs[c])
            cell.set_height(row_heights_frac[r])
            # Left-align the "Name" column if present
            label = str(col_labels[c]).lower()
            if label == "name":
                cell.get_text().set_ha('left')

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
        coords = line['coords']
        if hasattr(coords, "xy"):
            _, ys = coords.xy
        else:
            _, ys = zip(*coords)
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


def plot_inputs(
    slope_data,
    title="Slope Geometry and Inputs",
    figsize=(12, 6),
    mat_table=False,
    save_png=False,
    dpi=300,
    mode="lem",
    tab_loc="top",
    legend_ncol="auto",
    legend_max_cols=6,
    legend_max_rows=4,
):
    """
    Creates a plot showing the slope geometry and input parameters.

    Parameters:
        slope_data: Dictionary containing plot data
        title: Title for the plot
        figsize: Tuple of (width, height) in inches for the plot
        mat_table: Controls material table display. Can be:
            - True: Use tab_loc for positioning (default)
            - False: Don't show material table
            - 'auto': Use tab_loc for positioning
            - String: Specific location from valid placements (see tab_loc)
        save_png: If True, save plot as PNG file (default: False)
        dpi: Resolution for saved PNG file (default: 300)
        mode: Which material properties table to display:
            - "lem": Limit equilibrium materials (γ, c, φ, optional d/ψ)
            - "seep": Seepage properties (k₁, k₂, Angle, kr₀, h₀)
            - "fem": FEM properties (γ, c, φ, E, ν)
        tab_loc: Table placement when mat_table is True or 'auto'. Valid options:
            - "upper left": Top-left corner of plot area
            - "upper right": Top-right corner of plot area
            - "upper center": Top-center of plot area
            - "lower left": Bottom-left corner of plot area
            - "lower right": Bottom-right corner of plot area
            - "lower center": Bottom-center of plot area
            - "center left": Middle-left of plot area
            - "center right": Middle-right of plot area
            - "center": Center of plot area
            - "top": Above plot area, horizontally centered
        legend_ncol: Legend column count. Use "auto" (default) to choose a value that
            keeps the legend from getting too tall, or pass an int to force a width.
        legend_max_cols: When legend_ncol="auto", cap the number of columns (default: 6).
        legend_max_rows: When legend_ncol="auto", try to keep legend rows <= this value
            by increasing columns (default: 4).

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Plot contents
    plot_profile_lines(ax, slope_data['profile_lines'], materials=slope_data.get('materials'), labels=True)
    plot_max_depth(ax, slope_data['profile_lines'], slope_data['max_depth'])
    plot_piezo_line(ax, slope_data)
    if mode == "seep":
        plot_seepage_bc_lines(ax, slope_data)
    plot_dloads(ax, slope_data)
    plot_tcrack_surface(ax, slope_data['tcrack_surface'])
    plot_reinforcement_lines(ax, slope_data)

    if slope_data['circular']:
        plot_circles(ax, slope_data)
    elif slope_data.get('non_circ') and len(slope_data['non_circ']) > 0:
        plot_non_circ(ax, slope_data['non_circ'])

    # Handle material table display
    if mat_table:
        # Helpers to adapt slope_data materials into formats expected by table functions
        def _build_seep_data():
            materials = slope_data.get('materials', [])
            return {
                "k1_by_mat": [m.get('k1', 0.0) for m in materials],
                "k2_by_mat": [m.get('k2', 0.0) for m in materials],
                "angle_by_mat": [m.get('alpha', 0.0) for m in materials],
                "kr0_by_mat": [m.get('kr0', 0.0) for m in materials],
                "h0_by_mat": [m.get('h0', 0.0) for m in materials],
                "material_names": [m.get('name', f"Material {i+1}") for i, m in enumerate(materials)],
            }

        def _build_fem_data():
            materials = slope_data.get('materials', [])
            return {
                "c_by_mat": [m.get('c', 0.0) for m in materials],
                "phi_by_mat": [m.get('phi', 0.0) for m in materials],
                "E_by_mat": [m.get('E', 0.0) for m in materials],
                "nu_by_mat": [m.get('nu', 0.0) for m in materials],
                "gamma_by_mat": [m.get('gamma', 0.0) for m in materials],
                "material_names": [m.get('name', f"Material {i+1}") for i, m in enumerate(materials)],
            }

        def _estimate_table_dims():
            """Estimate table dimensions based on mode and materials."""
            materials = slope_data.get('materials', [])
            num_rows = max(1, len(materials))

            if mode == "lem":
                has_d_psi = any(mat.get('d', 0) > 0 or mat.get('psi', 0) > 0 for mat in materials)
                width = 0.25 if has_d_psi else 0.2
                height = min(0.35, 0.06 + 0.035 * num_rows)
            elif mode == "fem":
                width = 0.60
                height = min(0.32, 0.06 + 0.035 * num_rows)
            elif mode == "seep":
                width = 0.45
                height = min(0.50, 0.10 + 0.06 * num_rows)
            else:
                raise ValueError(f"Unknown mode '{mode}'. Expected one of: 'lem', 'seep', 'fem'.")

            return width, height

        def _plot_table(ax, xloc, yloc):
            """Plot the appropriate material table based on mode."""
            if mode == "lem":
                plot_lem_material_table(ax, slope_data['materials'], xloc=xloc, yloc=yloc)
            elif mode == "seep":
                plot_seep_material_table(ax, _build_seep_data(), xloc=xloc, yloc=yloc)
            elif mode == "fem":
                width, height = _estimate_table_dims()
                plot_fem_material_table(ax, _build_fem_data(), xloc=xloc, yloc=yloc, width=width, height=height)

        def _calculate_position(location, margin=0.03):
            """Calculate xloc, yloc for a given location string."""
            width, height = _estimate_table_dims()

            # Location map with default coordinates
            position_map = {
                'upper left': (margin, max(0.0, 1.0 - margin - height)),
                'upper right': (max(0.0, 1.0 - width - margin), max(0.0, 1.0 - margin - height)),
                'upper center': (0.35, 0.70),
                'lower left': (0.05, 0.05),
                'lower right': (0.70, 0.05),
                'lower center': (0.35, 0.05),
                'center left': (0.05, 0.35),
                'center right': (0.70, 0.35),
                'center': (0.35, 0.35),
                'top': ((1.0 - width) / 2.0, 1.16)
            }

            return position_map.get(location, position_map['upper right'])

        # Determine which location to use
        placement = mat_table if isinstance(mat_table, str) and mat_table != 'auto' else tab_loc

        # Validate placement
        valid_placements = ['upper left', 'upper right', 'upper center', 'lower left',
                          'lower right', 'lower center', 'center left', 'center right', 'center', 'top']
        if placement not in valid_placements:
            raise ValueError(f"Unknown placement '{placement}'. Expected one of: {', '.join(valid_placements)}.")

        # Calculate position and plot table
        xloc, yloc = _calculate_position(placement)
        _plot_table(ax, xloc, yloc)

        # Adjust y-limits to prevent table overlap with plot data
        if placement in ("upper left", "upper right", "upper center"):
            _, height = _estimate_table_dims()
            margin = 0.03
            bottom_fraction = max(0.0, 1.0 - margin - height)

            y_min_curr, y_max_curr = ax.get_ylim()
            y_range = y_max_curr - y_min_curr
            if y_range > 0:
                elem_bounds = get_plot_elements_bounds(ax, slope_data)
                if elem_bounds:
                    y_top = max(b[3] for b in elem_bounds)
                    y_norm = (y_top - y_min_curr) / y_range
                    if y_norm >= bottom_fraction and bottom_fraction > 0:
                        y_max_new = y_min_curr + (y_top - y_min_curr) / bottom_fraction
                        ax.set_ylim(y_min_curr, y_max_new)

    ax.set_aspect('equal')  # ✅ Equal aspect

    # Add a bit of headroom so plotted lines/markers don't touch the top border
    y0, y1 = ax.get_ylim()
    if y1 > y0:
        pad = 0.05 * (y1 - y0)
        ax.set_ylim(y0, y1 + pad)
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
    
    # --- Legend layout ---
    # Historically this legend used ncol=2 and was anchored below the axes.
    # If there are many entries, that makes the legend tall and it can fall
    # off the bottom of the window. We auto-increase columns to cap row count,
    # and we also reserve bottom margin so the legend stays visible.
    n_items = len(labels)
    if legend_ncol == "auto":
        # Choose enough columns to keep row count <= legend_max_rows (as best we can),
        # but never exceed legend_max_cols.
        required_cols = max(1, math.ceil(n_items / max(1, int(legend_max_rows))))
        ncol = min(int(legend_max_cols), required_cols)
        # Keep at least 2 columns once there's more than one entry (matches prior look).
        if n_items > 1:
            ncol = max(2, ncol)
    else:
        ncol = max(1, int(legend_ncol))

    n_rows = max(1, math.ceil(n_items / max(1, ncol)))
    # Reserve a bit more space as the legend grows so it doesn't get clipped.
    bottom_margin = min(0.45, 0.10 + 0.04 * n_rows)

    ax.legend(
        handles=handles,
        labels=labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=ncol,
    )

    ax.set_title(title)

    plt.subplots_adjust(bottom=bottom_margin)
    plt.tight_layout()
    
    if save_png:
        filename = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '') + '.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()

# ========== Main Plotting Function =========

def plot_solution(slope_data, slice_df, failure_surface, results, figsize=(12, 7), slice_numbers=False, save_png=False, dpi=300):
    """
    Plots the full solution including slices, numbers, thrust line, and base stresses.

    Parameters:
        data: Input data
        slice_df: DataFrame containing slice data
        failure_surface: Failure surface geometry
        results: Solution results
        figsize: Tuple of (width, height) in inches for the plot

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(False)

    plot_profile_lines(ax, slope_data['profile_lines'], materials=slope_data.get('materials'))
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
    else:
        title = f'{method}: FS = {fs:.3f}'
    ax.set_title(title)

    # zoom y‐axis to just cover the slope and depth, with a little breathing room (thrust line can be outside)
    ymin, ymax = compute_ylim(slope_data, slice_df, pad_fraction=0.05)
    ax.set_ylim(ymin, ymax)

    plt.tight_layout()
    
    if save_png:
        filename = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '').replace('°', 'deg') + '.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
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

def plot_circular_search_results(slope_data, fs_cache, search_path=None, highlight_fs=True, figsize=(12, 7), save_png=False, dpi=300):
    """
    Creates a plot showing the results of a circular failure surface search.

    Parameters:
        slope_data: Dictionary containing plot data
        fs_cache: List of dictionaries containing failure surface data and FS values
        search_path: List of dictionaries containing search path coordinates
        highlight_fs: Boolean indicating whether to highlight the critical failure surface
        figsize: Tuple of (width, height) in inches for the plot

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=figsize)

    plot_profile_lines(ax, slope_data['profile_lines'], materials=slope_data.get('materials'))
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
    
    if save_png:
        filename = 'plot_circular_search_results.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()

def plot_noncircular_search_results(slope_data, fs_cache, search_path=None, highlight_fs=True, figsize=(12, 7), save_png=False, dpi=300):
    """
    Creates a plot showing the results of a non-circular failure surface search.

    Parameters:
        slope_data: Dictionary containing plot data
        fs_cache: List of dictionaries containing failure surface data and FS values
        search_path: List of dictionaries containing search path coordinates
        highlight_fs: Boolean indicating whether to highlight the critical failure surface
        figsize: Tuple of (width, height) in inches for the plot

    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Plot basic profile elements
    plot_profile_lines(ax, slope_data['profile_lines'], materials=slope_data.get('materials'))
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
    
    if save_png:
        filename = 'plot_noncircular_search_results.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()

def plot_reliability_results(slope_data, reliability_data, figsize=(12, 7), save_png=False, dpi=300):
    """
    Creates a plot showing the results of reliability analysis.
    
    Parameters:
        slope_data: Dictionary containing plot data
        reliability_data: Dictionary containing reliability analysis results
        figsize: Tuple of (width, height) in inches for the plot
    
    Returns:
        None
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Plot basic slope elements (same as other search functions)
    plot_profile_lines(ax, slope_data['profile_lines'], materials=slope_data.get('materials'))
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
    
    if save_png:
        filename = 'plot_reliability_results.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()

def plot_mesh(mesh, materials=None, figsize=(14, 6), pad_frac=0.05, show_nodes=True, label_elements=False, label_nodes=False, save_png=False, dpi=300):
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
    
    if save_png:
        filename = 'plot_mesh.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()


def plot_polygons(
    polygons,
    materials=None,
    nodes=False,
    legend=True,
    title="Material Zone Polygons",
    figsize=(10, 6),
    save_png=False,
    dpi=300,
):
    """
    Plot all material zone polygons in a single figure.
    
    Parameters:
        polygons: List of polygon coordinate lists or dicts with "coords"/"mat_id"
        materials: Optional list of material dicts (with key "name") or list of material
            name strings. If provided, the material name will be used in the legend.
        nodes: If True, plot each polygon vertex as a dot.
        legend: If True, show the legend.
        title: Plot title
        figsize: Matplotlib figure size tuple, e.g. (10, 6)
    """
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=figsize)
    for i, polygon in enumerate(polygons):
        coords = polygon.get("coords", []) if isinstance(polygon, dict) else polygon
        xs = [x for x, y in coords]
        ys = [y for x, y in coords]
        mat_idx = polygon.get("mat_id") if isinstance(polygon, dict) else i
        if mat_idx is None:
            mat_idx = i
        mat_name = None
        if materials is not None and 0 <= mat_idx < len(materials):
            item = materials[mat_idx]
            if isinstance(item, dict):
                mat_name = item.get("name", None)
            elif isinstance(item, str):
                mat_name = item
        label = mat_name if mat_name else f"Material {mat_idx}"
        ax.fill(xs, ys, color=get_material_color(mat_idx), alpha=0.6, label=label)
        ax.plot(xs, ys, color=get_material_color(mat_idx), linewidth=1)
        if nodes:
            # Avoid legend clutter by not adding a label here.
            ax.scatter(xs, ys, color='k', s=30, marker='o', zorder=3)
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title(title)
    if legend:
        ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    plt.tight_layout()
    
    if save_png:
        filename = 'plot_' + title.lower().replace(' ', '_').replace(':', '').replace(',', '') + '.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
    plt.show()


def plot_polygons_separately(polygons, materials=None, save_png=False, dpi=300):
    """
    Plot each polygon in a separate matplotlib frame (subplot), with vertices as round dots.
    
    Parameters:
        polygons: List of polygon coordinate lists or dicts with "coords"/"mat_id"
        materials: Optional list of material dicts (with key "name") or list of material
            name strings. If provided, the material name will be included in each subplot title.
    """
    import matplotlib.pyplot as plt
    
    n = len(polygons)
    fig, axes = plt.subplots(n, 1, figsize=(8, 3 * n), squeeze=False)
    for i, polygon in enumerate(polygons):
        coords = polygon.get("coords", []) if isinstance(polygon, dict) else polygon
        xs = [x for x, y in coords]
        ys = [y for x, y in coords]
        ax = axes[i, 0]
        mat_idx = polygon.get("mat_id") if isinstance(polygon, dict) else i
        if mat_idx is None:
            mat_idx = i
        ax.fill(xs, ys, color=get_material_color(mat_idx), alpha=0.6, label=f'Material {mat_idx}')
        ax.plot(xs, ys, color=get_material_color(mat_idx), linewidth=1)
        ax.scatter(xs, ys, color='k', s=30, marker='o', zorder=3, label='Vertices')
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        mat_name = None
        if materials is not None and 0 <= mat_idx < len(materials):
            item = materials[mat_idx]
            if isinstance(item, dict):
                mat_name = item.get("name", None)
            elif isinstance(item, str):
                mat_name = item
        if mat_name:
            ax.set_title(f'Material {mat_idx}: {mat_name}')
        else:
            ax.set_title(f'Material {mat_idx}')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
        # Intentionally no legend: these plots are typically used for debugging geometry,
        # and legends can obscure key vertices/edges.
    plt.tight_layout()
    
    if save_png:
        filename = 'plot_polygons_separately.png'
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
    
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
                coords = line['coords']
                xs = [p[0] for p in coords]
                ys = [p[1] for p in coords]
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