from global_config import non_circ
from slice import build_ground_surface, generate_slices
from fileio import load_globals
from plot import plot_slope, plot_circular_search_results
from solve import oms, bishop, spencer, janbu_corrected
from search import circular_search


data = load_globals("docs/input_template.xlsx")

fs_cache, converged, search_path = circular_search(data, oms)

plot_circular_search_results(data, fs_cache, search_path)