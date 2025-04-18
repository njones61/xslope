import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def plot_slices(profile_lines, circle, df, piezo_line=None, failure_surface=None, FS=None, dloads=None):
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot profile lines
    for i, line in enumerate(profile_lines):
        xs, ys = zip(*line)
        ax.plot(xs, ys, label=f'Material {i+1}')

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

    if FS is not None:
        ax.set_title(f"Factor of Safety = {FS:.3f}")

    plt.tight_layout()
    plt.show()