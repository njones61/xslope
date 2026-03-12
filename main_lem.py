from xslope.global_config import non_circ

from xslope.fileio import load_slope_data
from xslope.plot import plot_circular_search_results, plot_noncircular_search_results, plot_solution, plot_reliability_results, plot_inputs
from xslope.solve import oms, bishop, spencer, janbu, corps_engineers, lowe_karafiath, solve_selected, solve_all
from xslope.search import circular_search, noncircular_search
from xslope.slice import generate_slices
from xslope.advanced import reliability as reliability_analysis

slope_data = load_slope_data("docs/lem/files/xslope_earth_dam_rapid.xlsx")
plot_inputs(slope_data, mode='lem', save_png=True)

method = "spencer" # @param ["oms","bishop","janbu","corps_engineers","lowe_karafiath","spencer"]
num_slices = 30 # @param {"type":"integer"}
analysis_type = "auto_search" # @param ["single_surface","all_methods", "auto_search","reliability"]
surface_type = "circular" # @param ["circular","non_circular"]
rapid_drawdown = True # @param {"type":"boolean"}
reliability = False # @param {"type":"boolean"}
save_png = True # @param {"type":"boolean"}
diagnostic = False # @param {"type":"boolean"}

if analysis_type == 'single_surface': # analyze the specified failure surface
  circle = slope_data['circles'][0] if slope_data['circular'] else None
  non_circ = slope_data['non_circ'] if slope_data['non_circ'] else None
  success, result = generate_slices(slope_data, circle=circle, non_circ=non_circ, num_slices=num_slices)
  if success:
      slice_df, failure_surface = result
      results = solve_selected(method, slice_df, rapid=rapid_drawdown)
      plot_solution(slope_data, slice_df, failure_surface, results, save_png=save_png)
  else:
      print(result)
      exit()

elif analysis_type == "auto_search": # automated search for critical surface
  if surface_type == "circular":
    fs_cache, converged, search_path, circle_cache = circular_search(slope_data, method, rapid=rapid_drawdown, num_slices=num_slices, diagnostic=diagnostic)
    plot_circular_search_results(slope_data, fs_cache, search_path, circle_cache=circle_cache, save_png=save_png)
  else:
    fs_cache, converged, search_path = noncircular_search(slope_data, method, rapid=rapid_drawdown, diagnostic=diagnostic)
    plot_noncircular_search_results(slope_data, fs_cache, search_path, save_png=save_png)

  # Extract critical failure surface (lowest FS is first in sorted list)
  critical_surface = fs_cache[0]
  slice_df = critical_surface['slices']
  failure_surface = critical_surface['failure_surface']
  results = critical_surface['solver_result']
  if rapid_drawdown and results and 'stage1_FS' in results:
    print(f"=== RAPID DRAWDOWN SUMMARY (Critical Surface) ===")
    print(f"Stage 1 FS = {results['stage1_FS']:.4f}")
    print(f"Stage 2 FS = {results['stage2_FS']:.4f}")
    print(f"Stage 3 FS = {results['stage3_FS']:.4f}")
    print(f"Final rapid drawdown FS = {results['FS']:.4f}")
  plot_solution(slope_data, slice_df, failure_surface, results, save_png=save_png)

elif analysis_type == "reliability": # reliability analysis (supports both circular and non-circular)
  circular = (surface_type == "circular")
  success, result = reliability_analysis(slope_data, method, rapid=rapid_drawdown, circular=circular, debug_level=1)
  if success:
    plot_reliability_results(slope_data, result, save_png=save_png)
  else:
    print(f"Reliability analysis failed: {result}")
