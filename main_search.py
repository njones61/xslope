from global_config import non_circ

from fileio import load_globals
from plot import plot_circular_search_results
from solve import oms, bishop, spencer, janbu_corrected
from search import circular_search


data = load_globals("docs/input_template.xlsx")

fs_cache, converged, search_path = circular_search(data, oms, diagnostic=True)

plot_circular_search_results(data, fs_cache, search_path)