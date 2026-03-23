"""Generate pile_diagram.png — conceptual diagram of piles stabilizing a slope."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch

fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# --- Slope geometry ---
# Ground surface: flat crest, 2H:1V slope, flat toe
x_ground = [0, 6, 6, 18, 24, 30]
y_ground = [12, 12, 12, 6, 6, 6]
# Smooth the slope transition
x_ground = np.array([0, 5, 7, 10, 13, 16, 19, 22, 25, 30])
y_ground = np.array([12, 12, 11.8, 10.5, 9.0, 7.5, 6.2, 6.0, 6.0, 6.0])

# Base of slope (bedrock/stable ground)
y_base = 0

# Fill the ground
ax.fill_between(x_ground, y_ground, y_base, color='#D2B48C', alpha=0.4, zorder=1)
ax.plot(x_ground, y_ground, 'k-', linewidth=2.0, zorder=5)
ax.plot([0, 30], [y_base, y_base], 'k-', linewidth=1.5, zorder=5)

# --- Circular failure surface ---
# Circle center and radius
Xo, Yo = 11.0, 22.0
R = 16.5
theta_arc = np.linspace(np.radians(225), np.radians(310), 200)
x_arc = Xo + R * np.cos(theta_arc)
y_arc = Yo + R * np.sin(theta_arc)

# Clip to within the ground
mask = np.ones_like(x_arc, dtype=bool)
y_ground_at_arc = np.interp(x_arc, x_ground, y_ground)
mask = (y_arc <= y_ground_at_arc) & (y_arc >= y_base) & (x_arc >= 0) & (x_arc <= 30)
x_arc_clipped = x_arc[mask]
y_arc_clipped = y_arc[mask]

ax.plot(x_arc_clipped, y_arc_clipped, 'r--', linewidth=2.0, zorder=6, label='Failure surface')

# --- Shade sliding mass ---
# Build polygon: ground surface from entry to exit, then failure arc back
x_entry = x_arc_clipped[0]
x_exit = x_arc_clipped[-1]
# Ground surface segment between entry and exit
gnd_mask = (x_ground >= x_entry) & (x_ground <= x_exit)
x_gnd_seg = x_ground[gnd_mask]
y_gnd_seg = y_ground[gnd_mask]
# Add entry and exit points
x_gnd_seg = np.concatenate([[x_entry], x_gnd_seg, [x_exit]])
y_gnd_seg = np.concatenate([[np.interp(x_entry, x_ground, y_ground)], y_gnd_seg,
                             [np.interp(x_exit, x_ground, y_ground)]])

# Polygon: ground forward, arc backward
x_poly = np.concatenate([x_gnd_seg, x_arc_clipped[::-1]])
y_poly = np.concatenate([y_gnd_seg, y_arc_clipped[::-1]])
ax.fill(x_poly, y_poly, color='#E8C88A', alpha=0.5, zorder=2)

# --- Piles ---
pile_x_positions = [12.5, 15.0, 17.5]
pile_width = 0.5
pile_color = '#555555'

for px in pile_x_positions:
    y_top = np.interp(px, x_ground, y_ground)
    y_bot = 1.5  # Embedded into stable ground below failure surface

    # Pile rectangle
    rect = patches.Rectangle((px - pile_width/2, y_bot), pile_width, y_top - y_bot,
                              linewidth=1.5, edgecolor='k', facecolor=pile_color,
                              zorder=8)
    ax.add_patch(rect)

    # Failure surface intersection
    # Find where the arc crosses this x
    y_fs = Yo - np.sqrt(max(R**2 - (px - Xo)**2, 0))
    ax.plot(px, y_fs, 'ro', markersize=6, zorder=10)

# --- Force arrows H at each pile ---
for px in pile_x_positions:
    y_fs = Yo - np.sqrt(max(R**2 - (px - Xo)**2, 0))
    arr = FancyArrowPatch((px - 0.3, y_fs), (px - 1.8, y_fs),
                           arrowstyle='->', mutation_scale=18,
                           color='blue', linewidth=2.5, zorder=11)
    ax.add_patch(arr)

# H label on middle pile only
px_mid = 15.0
y_fs_mid = Yo - np.sqrt(max(R**2 - (px_mid - Xo)**2, 0))
ax.text(px_mid - 1.0, y_fs_mid + 0.3, '$H$', fontsize=16, color='blue',
        va='bottom', ha='center', fontweight='bold', zorder=11)

# --- Movement arrows (sliding mass direction) ---
arrow_style = dict(arrowstyle='->', mutation_scale=14, color='#CC4400',
                   linewidth=1.8, zorder=9)
movement_pts = [(6, 10.5, 8, 9.5), (8.5, 8.8, 10.5, 7.8)]
for x1, y1, x2, y2 in movement_pts:
    arr = FancyArrowPatch((x1, y1), (x2, y2), **arrow_style)
    ax.add_patch(arr)

# --- Labels ---
# Place each label in clear open space, away from lines and piles.
# The sliding mass occupies roughly x=3-18, y=6-12 (between ground and arc).
# The stable zone is below the arc. Piles are at x=12.5, 15, 17.5.

ax.text(6, 12.4, 'Crest', fontsize=11, ha='center', va='bottom', zorder=12)
ax.text(20, 6.4, 'Toe', fontsize=11, ha='center', va='bottom', zorder=12)

# "Sliding mass" in open space above failure arc, between arrows and piles
ax.text(5, 9.0, 'Sliding mass', fontsize=12, ha='center', va='center',
        fontstyle='italic', color='#993300', zorder=12)

# "Stable zone" below the arc, left of piles where nothing overlaps
ax.text(6, 2.5, 'Stable zone', fontsize=12, ha='center', va='center',
        fontstyle='italic', color='#336633', zorder=12)

# "Failure surface" along the left portion of the arc
idx_fs = len(x_arc_clipped) // 6
ax.text(x_arc_clipped[idx_fs] - 0.3, y_arc_clipped[idx_fs] - 0.8,
        'Failure\nsurface', fontsize=10, ha='center', va='top', color='red', zorder=12)

# Pile label — short arrow from right side
ax.annotate('Piles', xy=(pile_x_positions[2] + 0.3, 4.5),
            xytext=(19.5, 4.5), fontsize=11, ha='left', va='center',
            arrowprops=dict(arrowstyle='->', color='black', lw=1.2),
            zorder=12)


# --- Formatting ---
ax.set_xlim(-0.5, 30.5)
ax.set_ylim(-0.5, 14)
ax.set_aspect('equal')
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.tick_params(labelsize=10)

# Hatch the bedrock/base
ax.fill_between([0, 30], [y_base, y_base], [-0.5, -0.5],
                hatch='///', facecolor='#E0E0E0', edgecolor='#999999',
                linewidth=0.5, zorder=0)

plt.tight_layout()
plt.savefig('docs/lem/images/pile_diagram.png', dpi=150, bbox_inches='tight',
            facecolor='white')
plt.close()
print("Saved pile_diagram.png")
