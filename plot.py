import matplotlib.pyplot as plt
import numpy as np

def plot_slices(profile_lines, circle, df):
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot profile lines
    for i, line in enumerate(profile_lines):
        xs, ys = zip(*line)
        ax.plot(xs, ys, label=f'Material {i+1}')

    # Plot slip circle
    Xo, Yo, D = circle['Xo'], circle['Yo'], circle['Depth']
    R = Yo - D
    theta = np.linspace(np.pi, 2 * np.pi, 500)
    circle_x = Xo + R * np.cos(theta)
    circle_y = Yo + R * np.sin(theta)
    ax.plot(circle_x, circle_y, 'k--', label="Slip Circle")

    # Plot slices
    for _, row in df.iterrows():
        ax.plot([row['xc'], row['xc']], [row['base_y'], row['top_y']], 'r-')

    ax.set_aspect('equal')
    ax.set_title("Soil Profile with Circular Slip Surface and Slices")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()