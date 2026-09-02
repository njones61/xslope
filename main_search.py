from xslope.fileio import load_slope_data
from xslope.plot import plot_circular_search_results, plot_noncircular_search_results
from xslope.search import circular_search, noncircular_search

slope_data = load_slope_data("docs/inputs/slope/xslope_reliability.xlsx")

# Circular search: an adaptive grid over trial centers, halving the grid spacing
# each time the best center stops moving. `method` can be any of 'oms', 'bishop',
# 'janbu', 'corps', 'lowe', 'spencer', 'mprice'.
fs_cache, converged, search_path, circle_cache = circular_search(slope_data, 'spencer', diagnostic=False)
plot_circular_search_results(slope_data, fs_cache, search_path, circle_cache=circle_cache)
print(f"Critical circular FS = {fs_cache[0]['FS']:.4f} (converged={converged})")

# Non-circular search: starts from the critical circle and moves the surface
# vertices one at a time.
# fs_cache, converged, search_path = noncircular_search(slope_data, 'spencer', diagnostic=False)
# plot_noncircular_search_results(slope_data, fs_cache, search_path)
