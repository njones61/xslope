# NOTE: This design driver only works with PROFILE-LINE based input templates.
# The slope-angle sweep edits profile_lines and regenerates the material polygons
# from them (see rebuild_geometry). Polygon-sheet inputs have no profile lines to
# reshape, so the slope angle cannot be varied — use a profile-based template.
#
# Two design studies are supported via `design_mode`:
#   "fs"          -- sweep the critical factor of safety vs slope angle and find the
#                    angle that yields a target FS (deterministic design).
#   "reliability" -- sweep the lognormal reliability index (beta) vs slope angle and
#                    find the angle that yields a target reliability index. Requires
#                    material standard deviations (mat sheet cols L-Q) in the input.
import numpy as np
import math
from xslope.global_config import non_circ

from shapely.geometry import Polygon
from xslope.fileio import load_slope_data, build_ground_surface_from_polygons
from xslope.mesh import build_polygons
from xslope.plot import plot_circular_search_results, plot_noncircular_search_results, plot_solution, plot_reliability_results, plot_inputs
from xslope.solve import oms, bishop, spencer, janbu, corps_engineers, lowe_karafiath, morgenstern_price, solve_selected, solve_all
from xslope.search import circular_search, noncircular_search
from xslope.slice import generate_slices
from xslope.advanced import reliability as reliability_analysis

slope_data = load_slope_data("docs/lem/files/xslope_acads_simple.xlsx")
plot_inputs(slope_data, mode='lem', save_png=True)

method = "bishop" # @param ["oms","bishop","janbu","corps_engineers","lowe_karafiath","spencer","morgenstern_price"]
# Bishop is used for the sweep: for simple circular surfaces it gives the same FS as
# Spencer but stays fast at steep angles. Spencer cascades into an expensive scipy
# fallback on steep/marginal circles, making each search 4-8x slower near beta2.
# Use Spencer only if you need non-circular surfaces or a force+moment method.
num_slices = 30 # @param {"type":"integer"}
design_mode = "reliability" # @param ["fs","reliability"]  -- sweep FS or reliability index vs slope angle
surface_type = "circular" # @param ["circular","non_circular"]
save_png = True # @param {"type":"boolean"}
diagnostic = False # @param {"type":"boolean"}

beta1 = 20  # start of slope-angle sweep, degrees
beta2 = 30  # end of slope-angle sweep, degrees
num_angles = 10  # number of slope angles to evaluate across [beta1, beta2]
toe_index = 1  # index of the toe point in the first profile line (zero-based index)
slope_index = 2  # index of the slope top point (toe_index+1 for right-facing, toe_index-1 for left-facing)

design_fs = 1.2            # target factor of safety   (design_mode="fs")
design_reliability = 0.75  # target reliability R = P(FS > 1)  (design_mode="reliability")
                           # The target must fall within the swept range. For the shipped
                           # sample (a weak c-phi soil), FS spans ~0.88-1.27 and R spans
                           # ~0.5-0.92 over beta=20-30°, so 1.2 and 0.75 bracket the sweep.


def rebuild_geometry(slope_data):
    """Rebuild the polygon-based geometry after editing profile_lines.

    Slice weights, layer heights, and base materials are all computed from
    slope_data['polygons'] (see generate_slices), NOT from 'ground_surface'.
    After moving a profile point we must regenerate the material polygons and
    re-derive the ground surface / domain polygon from them — exactly as
    load_slope_data does — or the weights stay pinned to the original geometry.
    """
    polys = [
        {'polygon': Polygon(p['coords']), 'mat_id': p['mat_id']}
        for p in build_polygons(slope_data={'profile_lines': slope_data['profile_lines'],
                                             'max_depth': slope_data.get('max_depth')})
    ]
    slope_data['polygons'] = polys
    ground_surface, domain_polygon = build_ground_surface_from_polygons(polys)
    slope_data['ground_surface'] = ground_surface
    slope_data['domain_polygon'] = domain_polygon


def run_search(slope_data):
    """Run the appropriate search for the selected surface_type and return the
    sorted fs_cache (lowest FS first), convergence flag, and search path."""
    if surface_type == "circular":
        fs_cache, converged, search_path, _circle_cache = circular_search(
            slope_data, method, num_slices=num_slices, diagnostic=diagnostic)
    else:
        fs_cache, converged, search_path = noncircular_search(
            slope_data, method, num_slices=num_slices, diagnostic=diagnostic)
    return fs_cache, converged, search_path


def evaluate(slope_data):
    """Evaluate the design metric for the current geometry.

    Returns (metric_value, payload) where payload carries whatever the final
    verification plot needs:
      - "fs":          metric = critical FS,         payload = critical fs_cache entry
      - "reliability": metric = reliability P(FS>1),  payload = reliability result dict
    """
    if design_mode == "fs":
        fs_cache, _converged, _search_path = run_search(slope_data)
        return fs_cache[0]['FS'], fs_cache[0]
    else:
        circular = (surface_type == "circular")
        success, result = reliability_analysis(
            slope_data, method, circular=circular, debug_level=0)
        if not success:
            raise RuntimeError(f"Reliability analysis failed: {result}")
        return result['reliability'], result


# Mode-specific configuration: target value, labels, and output filename.
if design_mode == "fs":
    target = design_fs
    metric_label = "Factor of Safety"
    metric_short = "FS"
    out_png = "fs_vs_slope_angle.png"
elif design_mode == "reliability":
    target = design_reliability
    metric_label = "Reliability,  P(FS > 1)"
    metric_short = "R"
    out_png = "reliability_vs_slope_angle.png"
else:
    raise ValueError(f"Unknown design_mode: {design_mode!r} (use 'fs' or 'reliability')")


def set_slope_angle(slope_data, beta_deg):
    """Move the slope-top point so the slope face makes angle beta_deg, then
    rebuild the polygon geometry. Returns the new (x_top, y_top)."""
    profile = slope_data['profile_lines'][0]['coords']
    x_toe, y_toe = profile[toe_index]
    _x_top, y_top = profile[slope_index]
    # tan(beta) = (y_top - y_toe) / (x_top - x_toe)
    x_top_new = x_toe + (y_top - y_toe) / math.tan(math.radians(beta_deg))
    slope_data['profile_lines'][0]['coords'][slope_index] = (x_top_new, y_top)
    rebuild_geometry(slope_data)
    return x_top_new, y_top


# ---------------------------------------------------------------------------
# Lightweight profiling. Transparent wrappers count solver/search/slice calls
# and accumulate solver time for THIS run only. They reassign in-memory module
# attributes (no xslope source files are modified) and are pass-throughs, so
# results are unchanged. Unlike cProfile this adds negligible overhead, so the
# reported wall time reflects the real run. A summary table prints at the end.
# ---------------------------------------------------------------------------
import time
import xslope.search as _search_mod
import xslope.solve as _solve_mod

_counts = {"reliability": 0, "search": 0, "generate_slices": 0, "solver": 0}
_solver_time = [0.0]

def _count(orig, key):
    def wrapped(*args, **kwargs):
        _counts[key] += 1
        return orig(*args, **kwargs)
    return wrapped

def _count_and_time_solver(orig):
    def wrapped(*args, **kwargs):
        _counts["solver"] += 1
        _t = time.perf_counter()
        try:
            return orig(*args, **kwargs)
        finally:
            _solver_time[0] += time.perf_counter() - _t
    return wrapped

# Search functions: wrap them in the search module so the calls made INSIDE
# reliability() (which does `from .search import ...` at call time) are counted,
# and rebind the module-level names that run_search() looks up.
_search_mod.circular_search = _count(_search_mod.circular_search, "search")
_search_mod.noncircular_search = _count(_search_mod.noncircular_search, "search")
circular_search = _search_mod.circular_search
noncircular_search = _search_mod.noncircular_search

# generate_slices: search.py holds its own binding, so wrap it there.
_search_mod.generate_slices = _count(_search_mod.generate_slices, "generate_slices")

# Solver: search.py resolves it via getattr(solve, method), so wrap that attribute.
setattr(_solve_mod, method, _count_and_time_solver(getattr(_solve_mod, method)))

# Reliability driver (counted via the name evaluate() looks up).
reliability_analysis = _count(reliability_analysis, "reliability")

_wall_start = time.perf_counter()


betas = np.linspace(beta1, beta2, num=num_angles)  # slope angles to evaluate
metric_results = np.zeros_like(betas)
payloads = [None] * len(betas)
for i, beta in enumerate(betas):
    x_top_new, y_top = set_slope_angle(slope_data, beta)
    metric_results[i], payloads[i] = evaluate(slope_data)
    if design_mode == "reliability":
        pf = payloads[i]['prob_failure']
        print(f"Slope angle: {beta:.1f}°, slope point x={x_top_new:.2f}  "
              f"R={metric_results[i]*100:.2f}%, Pf={pf*100:.2f}%")
    else:
        print(f"Slope angle: {beta:.1f}°, slope point x={x_top_new:.2f}  "
              f"FS={metric_results[i]:.3f}")

print("\nBeta sweep finished. Results:")
for beta, m in zip(betas, metric_results):
    print(f"  β_slope = {beta:.1f}°: {metric_short} = {m:.3f}")

# Both metrics decrease monotonically as the slope steepens, so interpolate the
# target with the arrays reversed (np.interp needs an ascending x table).
m_min, m_max = metric_results.min(), metric_results.max()
if not (m_min <= target <= m_max):
    print(f"\nTarget {metric_short}={target} is outside the computed range "
          f"[{m_min:.3f}, {m_max:.3f}].")
    print(f"Adjust beta1/beta2 so the range brackets {metric_short}={target} and re-run.")
    exit()

critical_slope_angle = float(np.interp(target, metric_results[::-1], betas[::-1]))
print(f"\nCritical slope angle for {metric_short}={target}: {critical_slope_angle:.2f}°")

# Re-run at the interpolated angle to confirm the metric is close to the target.
set_slope_angle(slope_data, critical_slope_angle)
metric_at_critical, payload = evaluate(slope_data)
print(f"{metric_short} at critical slope angle ({critical_slope_angle:.2f}°): "
      f"{metric_at_critical:.3f}")

if design_mode == "fs":
    plot_solution(slope_data, payload['slices'], payload['failure_surface'],
                  payload['solver_result'], save_png=save_png)
else:
    print(f"Probability of failure at critical angle: {payload['prob_failure']*100:.2f}%")
    plot_reliability_results(slope_data, payload, save_png=save_png)

# Plot the design metric vs slope angle
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(betas, metric_results, marker='o')

# Horizontal line at the design target
ax.axhline(y=target, color='r', linestyle='--', linewidth=0.8,
           label=f'{metric_short} = {target}')

# Vertical line at critical slope angle
ax.axvline(x=critical_slope_angle, color='gray', linestyle='--', linewidth=0.8,
           label=f'β = {critical_slope_angle:.1f}°')

# Mark the intersection point
ax.plot(critical_slope_angle, metric_at_critical, 's', color='r', markersize=8, zorder=5)

ax.set_xlabel('Slope Angle (degrees)')
ax.set_ylabel(metric_label)
ax.set_title(f'{metric_label} vs Slope Angle')
ax.legend()
ax.grid()
plt.tight_layout()
if save_png:
    plt.savefig(out_png, dpi=600)

# --- Profiling summary ---
_wall = time.perf_counter() - _wall_start
_search_name = "circular_search" if surface_type == "circular" else "noncircular_search"
print(f"\n===== main_design.py PROFILE ({method}, {design_mode} mode) =====")
print(f"{'Total wall time:':<34}{_wall:8.2f} s")
if design_mode == "reliability":
    print(f"{'reliability() calls:':<34}{_counts['reliability']:8d}")
print(f"{_search_name + ' calls:':<34}{_counts['search']:8d}")
print(f"{'generate_slices calls (circles):':<34}{_counts['generate_slices']:8d}")
print(f"{method + '() solver calls:':<34}{_counts['solver']:8d}   "
      f"(total {_solver_time[0]:.2f}s)")

plt.show()
