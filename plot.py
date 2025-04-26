import matplotlib.pyplot as plt
import numpy as np

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

def plot_circles(ax, circles):
    """
    Plots starting circles with center markers and arrows.

    Parameters:
        ax (matplotlib axis): The plotting axis
        circles (list of dicts): List of circles with 'Xo', 'Yo', 'R'

    Returns:
        None
    """
    for circle in circles:
        Xo = circle['Xo']
        Yo = circle['Yo']
        R = circle['R']
        theta = np.linspace(0, 2 * np.pi, 100)
        x_circle = Xo + R * np.cos(theta)
        y_circle = Yo + R * np.sin(theta)
        ax.plot(x_circle, y_circle, 'r--', label='Circle')

        # Center marker
        ax.plot(Xo, Yo, 'r+', markersize=10)

        # Arrow direction (fixed downward, angle = -90 degrees)
        angle = -np.pi / 2
        dx = np.cos(angle)
        dy = np.sin(angle)

        # Shorten arrow shaft slightly
        shaft_length = R - 5
        ax.arrow(Xo, Yo, dx * shaft_length, dy * shaft_length,
                 head_width=5, head_length=5, fc='red', ec='red')

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