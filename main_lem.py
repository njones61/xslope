from xslope.global_config import non_circ

from xslope.fileio import load_slope_data
from xslope.plot import plot_circular_search_results, plot_noncircular_search_results, plot_solution, plot_reliability_results, plot_inputs
from xslope.solve import oms, bishop, spencer, janbu, corps_engineers, lowe_karafiath, solve_selected, solve_all
from xslope.search import circular_search, noncircular_search
from xslope.slice import generate_slices
from xslope.advanced import reliability as reliability_analysis

slope_data = load_slope_data("inputs/slope/input_template_reliability6.xlsx")
plot_inputs(slope_data)

method = "bishop" # @param ["oms","bishop","janbu","corps_engineers","lowe_karafiath","spencer"]
num_slices = 20 # @param {"type":"integer"}
analysis_type = "single_surface" # @param ["single_surface","all_methods", "auto_search","reliability"]
surface_type = "circular" # @param ["circular","non_circular"]
rapid_drawdown = False # @param {"type":"boolean"}
reliability = False # @param {"type":"boolean"}
save_png = True # @param {"type":"boolean"}


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
    fs_cache, converged, search_path = circular_search(slope_data, method, rapid=rapid_drawdown, diagnostic=False)
    plot_circular_search_results(slope_data, fs_cache, search_path, save_png=save_png)
  else:
    fs_cache, converged, search_path = noncircular_search(slope_data, method, rapid=rapid_drawdown, diagnostic=False)
    plot_noncircular_search_results(slope_data, fs_cache, search_path, save_png=save_png)

  # Extract critical failure surface (lowest FS is first in sorted list)
  critical_surface = fs_cache[0]
  slice_df = critical_surface['slices']
  failure_surface = critical_surface['failure_surface']
  results = critical_surface['solver_result']
  plot_solution(slope_data, slice_df, failure_surface, results, save_png=save_png)

elif analysis_type == "reliability": # reliability analysis (supports both circular and non-circular)
  circular = (surface_type == "circular")
  success, result = reliability_analysis(slope_data, method, rapid=rapid_drawdown, circular=circular, debug_level=1)
  if success:
    plot_reliability_results(slope_data, result, save_png=save_png)
  else:
    print(f"Reliability analysis failed: {result}")
