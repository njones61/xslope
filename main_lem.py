from xslope.global_config import non_circ

from xslope.fileio import load_slope_data
from xslope.plot import plot_circular_search_results, plot_noncircular_search_results, plot_solution, plot_reliability_results, plot_inputs
from xslope.solve import oms, bishop, spencer, janbu, corps_engineers, lowe_karafiath, solve_selected, solve_all
from xslope.search import circular_search, noncircular_search
from xslope.slice import generate_slices
from xslope.advanced import reliability as reliability_analysis

slope_data = load_slope_data("docs/lem/files/xslope_piles.xlsx")
plot_inputs(slope_data, mode='lem', save_png=True)

method = "spencer" # @param ["oms","bishop","janbu","corps_engineers","lowe_karafiath","spencer"]
num_slices = 40 # @param {"type":"integer"}
analysis_type = "auto_search" # @param ["single_surface","all_methods", "auto_search","reliability"]
surface_type = "circular" # @param ["circular","non_circular"]
rapid_drawdown = False # @param {"type":"boolean"}
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
  if slice_df is not None and 'h_pile' in slice_df.columns:
    pile_slices = slice_df[slice_df['h_pile'] > 0]
    if not pile_slices.empty and slope_data.get('pile_lines'):
      for pile_info, (_, ps) in zip(slope_data['pile_lines'], pile_slices.iterrows()):
        if pile_info.get('H') is None and pile_info.get('D_pile') and pile_info.get('S'):
          from shapely.geometry import Point as Pt
          from xslope.ito_matsui import ito_matsui_coefficients, intersect_pile_with_materials, compute_ito_matsui_force
          D_p, S_p = pile_info['D_pile'], pile_info['S']
          D1 = S_p - D_p
          y_gnd = slope_data['ground_surface'].interpolate(
              slope_data['ground_surface'].project(Pt(ps['x_pile'], 1e6))).y
          depth = y_gnd - ps['y_pile']
          segments = intersect_pile_with_materials(
              ps['x_pile'], y_gnd, ps['y_pile'],
              slope_data['profile_lines'], slope_data['materials'])
          H_val, F_pile = compute_ito_matsui_force(D_p, S_p, segments)
          A1, A2 = ito_matsui_coefficients(D1, D_p, segments[0]['phi']) if segments else (0, 0)
          print(f'  === Ito & Matsui Summary (Pile at x={ps["x_pile"]:.2f}) ===')
          print(f'  Pile diameter (D)          = {D_p:.1f}')
          print(f'  Pile spacing (S)           = {S_p:.1f}')
          print(f'  Clear spacing (D1 = S - D) = {D1:.1f}')
          print(f'  Depth to failure surface   = {depth:.1f}')
          print(f'  Coefficients: A1 = {A1:.3f}, A2 = {A2:.3f}')
          print(f'  Force per pile (F_pile)    = {F_pile:.0f}')
          print(f'  Force per unit width (H)   = {H_val:.1f}')
  if rapid_drawdown and results and 'stage1_FS' in results:
    print(f"=== RAPID DRAWDOWN SUMMARY (Critical Surface) ===")
    print(f"Stage 1 FS = {results['stage1_FS']:.4f}")
    print(f"Stage 2 FS = {results['stage2_FS']:.4f}")
    print(f"Stage 3 FS = {results['stage3_FS']:.4f}")
    print(f"Final rapid drawdown FS = {results['FS']:.4f}")
  if results is None:
    print(f"[⚠️] No valid solution found (FS=9999 for all surfaces).")
    # Check if Ito & Matsui forces may be the cause
    from xslope.slice import _ito_matsui_warned
    if _ito_matsui_warned:
      print(f"    Ito & Matsui pile forces exceed driving forces for all trial surfaces.")
      print(f"    This may mean the piles are sufficient to stabilize the slope. Note that")
      print(f"    Ito & Matsui computes an upper bound on soil resistance and does not check")
      print(f"    pile structural capacity. The available resistance should be taken as the")
      print(f"    lesser of the Ito & Matsui soil resistance and the pile shear/moment")
      print(f"    capacity. If the structural capacity governs, specify H directly in the")
      print(f"    piles sheet as H = min(structural capacity, Ito & Matsui) / S.")
    else:
      print(f"    Check input parameters.")
  else:
    plot_solution(slope_data, slice_df, failure_surface, results, save_png=save_png)

elif analysis_type == "reliability": # reliability analysis (supports both circular and non-circular)
  circular = (surface_type == "circular")
  success, result = reliability_analysis(slope_data, method, rapid=rapid_drawdown, circular=circular, debug_level=1)
  if success:
    plot_reliability_results(slope_data, result, save_png=save_png)
  else:
    print(f"Reliability analysis failed: {result}")
