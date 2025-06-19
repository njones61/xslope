from global_config import non_circ

from fileio import load_globals
from plot import plot_circular_search_results, plot_noncircular_search_results
from solve import oms, bishop, spencer, janbu, corps_engineers, lowe_karafiath
from search import circular_search, noncircular_search


data = load_globals("docs/input_template_lface.xlsx")

# Run non-circular search
fs_cache, converged, search_path = noncircular_search(data, corps_engineers, diagnostic=False)
plot_noncircular_search_results(data, fs_cache, search_path)

# For circular search:
# fs_cache, converged, search_path = circular_search(data, oms, diagnostic=False)
# plot_circular_search_results(data, fs_cache, search_path)


# import cProfile
# import pstats
#
# cProfile.run('circular_search(data, oms, diagnostic=False)', 'profile_output')
#
# # Then view the results:
# p = pstats.Stats('profile_output')
# p.strip_dirs().sort_stats('cumtime').print_stats(30)  # Top 30 slowest functions by cumulative time