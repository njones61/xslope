import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plot_slices(profile_lines, circle, df, piezo_line=None, failure_surface=None):
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

        if not pd.isna(row['piezo_y']):
            ax.plot(row['x_c'], row['piezo_y'], 'bo')

    # Plot piezometric line
    if piezo_line:
        piezo_xs, piezo_ys = zip(*piezo_line)
        ax.plot(piezo_xs, piezo_ys, 'b-', label="Piezometric Line")

    ax.set_aspect('equal')
    ax.set_title("Soil Profile with Failure Surface and Slices")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()