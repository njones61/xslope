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

import xslope as xslope

from xslope.global_config import non_circ
from xslope.slice import generate_slices
from xslope.fileio import load_slope_data, load_data_from_pickle
from xslope.plot import plot_circular_search_results, plot_inputs, plot_solution
from xslope.solve import solve_selected, solve_all
from xslope.search import circular_search, noncircular_search


slope_data = load_slope_data("docs/inputs/slope/input_template_griffiths1_6.xlsx")

# plot_inputs(slope_data)

circle = slope_data['circles'][0] if slope_data['circular'] else None
non_circ = slope_data['non_circ'] if slope_data['non_circ'] else None

success, result = generate_slices(slope_data, circle=circle, non_circ=None, num_slices=20)

if success:
    slice_df, failure_surface = result
else:
    print(result)
    exit()

methods = ['oms', 'bishop', 'janbu', 'corps', 'lowe', 'spencer']
results = solve_selected(methods[5], slice_df, rapid=False)  # spencer = methods[5]
plot_solution(slope_data, slice_df, failure_surface, results)


solve_all(slice_df)


fs_cache, converged, search_path, circle_cache = circular_search(slope_data, methods[5], diagnostic=False)
plot_circular_search_results(slope_data, fs_cache, search_path)

