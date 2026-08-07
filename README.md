# xslope

**xslope** is an open-source Python package for geotechnical slope stability and
seepage analysis. It answers the question a slope engineer has to answer on every
embankment, cut, dam, and levee — *is this slope stable, and by how much?* — with
three analysis engines that share one input file:

- **Limit equilibrium (LEM)** — the method of slices with seven methods: Ordinary
  Method of Slices, Simplified Janbu, Bishop's Simplified, Corps of Engineers,
  Lowe & Karafiath, Spencer, and Morgenstern–Price. Circular and non-circular
  surfaces, automated critical-surface search, rapid drawdown, reinforcement,
  seismic loading, and Monte Carlo reliability analysis.
- **Finite element seepage** — saturated/unsaturated groundwater flow, steady-state
  or transient, on a mesh generated from the same slope geometry. It produces the
  pore pressure field used by the stability analyses, rather than assuming pore
  pressures from a piezometric line, and also serves as a standalone 2D groundwater
  flow solver.
- **Finite element slope stability** — elastic–perfectly plastic Mohr–Coulomb
  analysis with the Shear Strength Reduction Method, which lets the failure
  mechanism emerge instead of requiring an assumed failure surface.

Problems are defined in an Excel template, which keeps the workflow accessible to
practitioners while the analysis and plotting run in Python. A companion desktop
application, **XSLOPE Studio**, provides a point-and-click interface to the same
engine.

> **⚠️ Beta — under active development.** xslope is still in development mode:
> changes land daily, interfaces and input templates may shift between releases,
> and results should be independently verified before use in practice. We hope
> to formally issue version 1.0 soon. Feedback and issue reports are welcome in
> the meantime.

## Installation

xslope requires Python 3.9 or later and is published on
[PyPI](https://pypi.org/project/xslope/):

```bash
pip install xslope
```

That installs everything needed for limit equilibrium analysis. Optional extras add
the remaining capabilities:

```bash
pip install "xslope[fem]"       # seepage and finite element analysis (adds gmsh)
pip install "xslope[gui]"       # XSLOPE Studio, the desktop application
pip install "xslope[ai]"        # Studio's AI assistant
pip install "xslope[cad]"       # DXF import/export
pip install "xslope[gui,fem,ai]"  # everything
```

Studio is launched with the `xslope-studio` command. On Debian/Ubuntu Linux
(including Google Colab), gmsh needs system OpenGL libraries — run
`apt-get install -y libgl1 libglu1-mesa` before installing the `fem` extra.

## Example

Load a slope problem, build slices on its failure surface, solve for the factor of
safety, and plot the result:

```python
from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope.solve import solve_selected
from xslope.plot import plot_solution

slope_data = load_slope_data("docs/inputs/slope/xslope_simple1.xlsx")

circle = slope_data['circles'][0]
success, result = generate_slices(slope_data, circle=circle, num_slices=20)
slice_df, failure_surface = result

results = solve_selected('bishop', slice_df)
print(f"FS = {results['FS']:.3f}")          # FS = 1.279

plot_solution(slope_data, slice_df, failure_surface, results)
```

To search for the critical surface instead of analyzing a single one:

```python
from xslope.search import circular_search
from xslope.plot import plot_circular_search_results

fs_cache, converged, search_path, circle_cache = circular_search(slope_data, 'bishop')
print(f"minimum FS = {fs_cache[0]['FS']:.3f}")   # minimum FS = 1.215

plot_circular_search_results(slope_data, fs_cache, search_path)
```

Sample input files, including the one above, are in `docs/inputs/slope/`, and a
blank copy of the input template ships with the package
(`xslope.fileio.default_template_path()`).

## Documentation

Full documentation — theory, input template reference, worked examples, XSLOPE
Studio, and verification against published benchmarks — is at
**[xslope.org](https://xslope.org)**.

## Citation

If you use xslope in published work, please cite the archived release:

> Jones, N. L. (2026). *xslope* (Version 0.2.1) [Computer software].
> https://doi.org/10.5281/zenodo.21830232

Machine-readable metadata is in [CITATION.cff](CITATION.cff).

## License

This project is licensed under the Apache License, Version 2.0 - see the [LICENSE.txt](LICENSE.txt) file for details.

## Copyright

Copyright 2025 Norman L. Jones

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
