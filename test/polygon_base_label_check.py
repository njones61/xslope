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

"""A polygon model's base is labeled "Bottom boundary"; a profile model's keeps
"Max Depth".

A profile-line model's base is its Max Depth, an input the user types. A polygon
model has no Max Depth: its base is the bottom of its zones, drawn rather than
typed, and Tutorial LEM-6 says so -- while every plot of that model labeled the
hatched base line "Max Depth". Studio's editors and the docs call that line the
model's bottom boundary, and the plots now do too.

  A. **a polygon model** (LEM-6's sloping-bottom model): the inputs plot's legend
     names the base "Bottom boundary" and has no "Max Depth" entry;
  B. **a profile model** (LEM-1's simple embankment): the legend keeps
     "Max Depth" and has no "Bottom boundary" entry;
  C. **the Styles dialog** styles both lines with one row, and the row names both:
     "Max depth / bottom boundary" (skipped without PySide6).

Run directly:  PYTHONPATH=. python3 test/polygon_base_label_check.py
"""

import sys
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parent.parent
POLYGON = ROOT / 'docs' / 'lem' / 'files' / 'xslope_sloping_bottom.xlsx'
PROFILE = ROOT / 'docs' / 'lem' / 'files' / 'xslope_simple_embankment.xlsx'
STYLE_ROW = "Max depth / bottom boundary"


def _legend_labels(path):
    """The legend entries of the model's inputs plot."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from xslope.fileio import load_slope_data
    from xslope.plot import plot_inputs

    sd = load_slope_data(str(path))
    fig = plt.figure(figsize=(12, 7))
    try:
        plot_inputs(sd, fig=fig)
        labels = []
        for ax in fig.axes:
            leg = ax.get_legend()
            if leg is not None:
                labels += [t.get_text() for t in leg.get_texts()]
        for leg in fig.legends:
            labels += [t.get_text() for t in leg.get_texts()]
        return sd, labels
    finally:
        plt.close(fig)


def check_the_polygon_model(fails):
    """A: the polygon model's base is the Bottom boundary."""
    sd, labels = _legend_labels(POLYGON)
    if sd.get('profile_lines') or not sd.get('polygons'):
        fails.append(f"A: {POLYGON.name} is not a polygon model")
    if 'Bottom boundary' not in labels:
        fails.append(f"A: the polygon model's legend has no 'Bottom boundary': "
                     f"{labels}")
    if 'Max Depth' in labels:
        fails.append("A: the polygon model's legend still says 'Max Depth'")


def check_the_profile_model(fails):
    """B: the profile model's base keeps its input's name."""
    sd, labels = _legend_labels(PROFILE)
    if not sd.get('profile_lines'):
        fails.append(f"B: {PROFILE.name} is not a profile-line model")
    if 'Max Depth' not in labels:
        fails.append(f"B: the profile model's legend has no 'Max Depth': "
                     f"{labels}")
    if 'Bottom boundary' in labels:
        fails.append("B: the profile model's legend says 'Bottom boundary'")


def check_the_styles_row(fails):
    """C: the one Styles row for both lines names both."""
    try:
        from studio.styles_dialog import FEATURES
    except ImportError:
        return                              # no PySide6: no Studio to check
    rows = [label for label, key in FEATURES if key == 'max_depth']
    if rows != [STYLE_ROW]:
        fails.append(f"C: the Styles dialog's base-line row reads {rows}, not "
                     f"[{STYLE_ROW!r}]")


def run():
    fails = []
    for check in (check_the_polygon_model, check_the_profile_model,
                  check_the_styles_row):
        try:
            check(fails)
        except Exception as exc:                             # pragma: no cover
            fails.append(f"{check.__name__} raised {type(exc).__name__}: {exc}")
    return fails


if __name__ == '__main__':
    failures = run()
    for f in failures:
        print("FAIL", f)
    print("polygon_base_label:", "FAIL" if failures else "PASS")
    sys.exit(1 if failures else 0)
