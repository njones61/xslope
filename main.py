# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from xslope.slice import generate_slices
from xslope.fileio import load_slope_data
from xslope.plot import plot_circular_search_results, plot_inputs, plot_solution
from xslope.solve import solve_selected, solve_all
from xslope.search import circular_search

slope_data = load_slope_data("docs/inputs/slope/xslope_simple1.xlsx")

plot_inputs(slope_data, mode='lem')

# Solve the first circle on the input sheet with one method, then compare every
# method on the same slices.
circle = slope_data['circles'][0]
success, result = generate_slices(slope_data, circle=circle, num_slices=20)

if not success:
    print(result)
    exit()

slice_df, failure_surface = result

results = solve_selected('spencer', slice_df, rapid=False)
plot_solution(slope_data, slice_df, failure_surface, results)

solve_all(slice_df)

# Now search for the critical circle rather than solving the one on the sheet.
fs_cache, converged, search_path, circle_cache = circular_search(slope_data, 'spencer', diagnostic=False)
plot_circular_search_results(slope_data, fs_cache, search_path, circle_cache=circle_cache)
