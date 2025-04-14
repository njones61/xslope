import matplotlib.pyplot as plt
import numpy as np

def plot_slices(profile_lines, circle, df, piezo_line=None):
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot profile lines
    for i, line in enumerate(profile_lines):
        xs, ys = zip(*line)
        ax.plot(xs, ys, label=f'Material {i+1}')

    # Plot slip circle (bottom half only)
    Xo, Yo, D = circle['Xo'], circle['Yo'], circle['Depth']
    R = Yo - D
    theta = np.linspace(np.pi, 2 * np.pi, 500)
    circle_x = Xo + R * np.cos(theta)
    circle_y = Yo + R * np.sin(theta)
    ax.plot(circle_x, circle_y, 'k--', label="Slip Circle")

    # Plot slices
    for _, row in df.iterrows():
        ax.plot([row['xc'], row['xc']], [row['base_y'], row['top_y']], 'r-')
        # Add piezo point if available
        if row['py'] is not None:
            ax.plot(row['xb'], row['py'], 'bo')  # blue dot at piezo elevation
            ax.vlines(row['xb'], row['base_y'], row['py'], colors='blue', linestyles='dotted')

    # Plot piezometric line if provided
    if piezo_line:
        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, 'b-', label="Piezometric Line")

    ax.set_aspect('equal')
    ax.set_title("Soil Profile with Circular Slip Surface, Piezometric Line, and Slices")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()
