import numpy as np
import math
from xslope.global_config import non_circ

from xslope.fileio import load_slope_data, build_ground_surface
from xslope.plot import plot_circular_search_results, plot_noncircular_search_results, plot_solution, plot_reliability_results, plot_inputs
from xslope.solve import oms, bishop, spencer, janbu, corps_engineers, lowe_karafiath, solve_selected, solve_all
from xslope.search import circular_search, noncircular_search
from xslope.slice import generate_slices
from xslope.advanced import reliability as reliability_analysis

slope_data = load_slope_data("docs/lem/files/xslope_design2.xlsx")
plot_inputs(slope_data, mode='lem', save_png=True)

method = "spencer" # @param ["oms","bishop","janbu","corps_engineers","lowe_karafiath","spencer"]
num_slices = 30 # @param {"type":"integer"}
analysis_type = "single_surface" # @param ["single_surface","all_methods", "auto_search","reliability"]
surface_type = "circular" # @param ["circular","non_circular"]
save_png = True # @param {"type":"boolean"}
diagnostic = False # @param {"type":"boolean"}

beta1 = 25  # desired slope angle in degrees
beta2 = 35  # desired slope angle in degrees
design_fs = 1.4  # target factor of safety for design
toe_index = 1  # index of the toe point in the first profile line (zero-based index)
slope_index = 2  # index of the slope top point (toe_index+1 for right-facing, toe_index-1 for left-facing)

betas = np.linspace(beta1, beta2, num=10)  # generate slope angles between beta1 and beta2
fs_results = np.zeros_like(betas)
for i, beta in enumerate(betas):

  profile = slope_data['profile_lines'][0]['coords']
  x_toe, y_toe = profile[toe_index]
  x_top, y_top = profile[slope_index]

  # Calculate new x_top from beta: tan(beta) = (y_top - y_toe) / (x_top - x_toe)
  x_top_new = x_toe + (y_top - y_toe) / math.tan(math.radians(beta))
  slope_data['profile_lines'][0]['coords'][slope_index] = (x_top_new, y_top)
  slope_data['ground_surface'] = build_ground_surface(slope_data['profile_lines'])

  print(f"Slope angle: {beta:.1f}°, toe: ({x_toe}, {y_toe}), slope point: ({x_top_new:.2f}, {y_top})")

  fs_cache, converged, search_path, circle_cache = circular_search(slope_data, method, num_slices=num_slices, diagnostic=diagnostic)
  fs_results[i] = fs_cache[0]['FS']

print("Beta sweep finished. Results:")
for beta, fs in zip(betas, fs_results):
    print(f"  β = {beta:.1f}\u00b0: FS = {fs:.3f}")

# From the results, infer what the critical slope angle is to produce FS = design_fs.
fs_min, fs_max = fs_results.min(), fs_results.max()
if not (fs_min <= design_fs <= fs_max):
    print(f"Design FS={design_fs} is outside the computed range [{fs_min:.3f}, {fs_max:.3f}].")
    print(f"Adjust beta1/beta2 so the range brackets FS={design_fs} and re-run.")
    exit()

critical_slope_angle = float(np.interp(design_fs, fs_results[::-1], betas[::-1]))
print(f"Critical slope angle for FS={design_fs}: {critical_slope_angle:.2f}°")

# Redo the analysis at the critical slope angle to confirm FS is close to design_fs
profile = slope_data['profile_lines'][0]['coords']
x_toe, y_toe = profile[toe_index]
x_top, y_top = profile[slope_index]
x_top_new = x_toe + (y_top - y_toe) / math.tan(math.radians(critical_slope_angle))
slope_data['profile_lines'][0]['coords'][slope_index] = (x_top_new, y_top)
slope_data['ground_surface'] = build_ground_surface(slope_data['profile_lines'])
fs_cache, converged, search_path, circle_cache = circular_search(slope_data, method, num_slices=num_slices, diagnostic=diagnostic)
fs_at_critical = fs_cache[0]['FS']
print(f"FS at critical slope angle ({critical_slope_angle:.2f}°): {fs_at_critical:.3f}")  
critical_surface = fs_cache[0]
slice_df = critical_surface['slices']
failure_surface = critical_surface['failure_surface']
results = critical_surface['solver_result']
plot_solution(slope_data, slice_df, failure_surface, results, save_png=save_png)

# Plot FS vs slope angle
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(betas, fs_results, marker='o')

# Horizontal line at FS = design_fs
ax.axhline(y=design_fs, color='r', linestyle='--', linewidth=0.8, label=f'FS = {design_fs}')

# Vertical line at critical slope angle
ax.axvline(x=critical_slope_angle, color='gray', linestyle='--', linewidth=0.8,
           label=f'β = {critical_slope_angle:.1f}°')

# Mark the intersection point
ax.plot(critical_slope_angle, fs_at_critical, 's', color='r', markersize=8, zorder=5)

ax.set_xlabel('Slope Angle (degrees)')
ax.set_ylabel('Factor of Safety')
ax.set_title('Factor of Safety vs Slope Angle')
ax.legend()
ax.grid()
plt.tight_layout()
if save_png:
    plt.savefig('fs_vs_slope_angle.png')
plt.show()

