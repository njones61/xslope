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

"""The FS versus time figure's legend covers none of the figure's own words.

The figure labels its own marks: the full-pool reference carries its value at the
right edge, the lowest point is ringed and annotated, and the drawdown band is
named at its top. Its legend was placed by matplotlib's ``loc='best'``, which
counts a covered label the same as one crossed line and cannot see the schedule
drawn on the twin axis; on Tutorial COMBO-3's comparison of the rapid drawdown
and single-stage curves it settled in the upper right, over ``full pool, 1.470``.

  A. **the three COMBO-3 figures**, drawn from the tables the tutorial publishes
     (Part 1's curve, Part 2's rapid drawdown curve with each instant's governing
     stage, and that curve with the single-stage curve beside it), at the
     function's own size and at a narrower and a wider one: the legend lies
     inside the figure and overlaps no text and no annotation arrow on the axes.
     Each is drawn twice, the second time with ``loc='best'`` blind to text, as
     matplotlib was before it learned to avoid it: that is how the shipped
     comparison figure came out with its legend on the full-pool label, and the
     placement must not depend on the matplotlib that draws it;
  B. **no free corner.** Where every position inside the axes would cover a
     label, the legend goes below the axes: under the tick labels and the axis
     label, overlapping neither, and inside the figure.

The tables are the page's own (docs/tutorials/combo03_fs_vs_time.md), so no
seepage is marched and no search is run; the check takes seconds.

Run directly:  PYTHONPATH=. python3 test/fs_vs_time_legend_check.py
"""

import contextlib
import sys
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

ROOT = Path(__file__).resolve().parent.parent
FILES = ROOT / 'docs' / 'tutorials' / 'files'
COMBO03 = FILES / 'xslope_earth_dam_fs_time.xlsx'
COMBO03R = FILES / 'xslope_johnson_fs_time.xlsx'

#: COMBO-3 Part 1's published curve: (t, FS, Xo, Yo, R).
PART1 = (
    (0.0, 1.5311, 103.00, 56.79, 56.79), (2.0, 1.5311, 103.00, 56.79, 56.79),
    (5.0, 1.5132, 5.20, 65.90, 65.88), (10.0, 1.4572, 5.20, 65.90, 65.88),
    (15.0, 1.4128, 5.20, 65.90, 65.88), (20.0, 1.3773, 6.25, 60.58, 60.58),
    (25.0, 1.3509, 7.00, 57.81, 57.75), (30.0, 1.3344, 7.00, 56.91, 56.91),
    (35.0, 1.3313, 7.00, 56.91, 56.91), (40.0, 1.3436, 7.00, 56.91, 56.91),
    (47.0, 1.3813, 7.00, 57.16, 57.16), (55.0, 1.4386, 7.00, 57.16, 57.16),
    (65.0, 1.4828, 7.00, 57.16, 57.16), (80.0, 1.5187, 7.00, 57.16, 57.16),
    (100.0, 1.5482, 103.00, 56.79, 56.79), (130.0, 1.5516, 103.00, 56.79, 56.79),
    (180.0, 1.5550, 103.00, 56.79, 56.79), (240.0, 1.5566, 103.00, 56.79, 56.79),
    (300.0, 1.5572, 103.00, 56.79, 56.79),
)
#: COMBO-3 Part 2's published drawdown table: (t, stage 1, stage 2, stage 3, FS,
#: governs, Xo, Yo, R). The drawdown sweep starts after stage 1 (t = 0), as the
#: tutorial's own run does, so it has no t = 0 row.
PART2 = (
    (5, 1.5242, 1.4697, None, 1.4697, 2, 242.70, 257.12, 167.06),
    (10, 1.5267, 1.3659, None, 1.3659, 2, 244.05, 254.46, 165.56),
    (15, 1.5349, 1.2886, None, 1.2886, 2, 246.78, 249.06, 162.07),
    (20, 1.5477, 1.2232, None, 1.2232, 2, 249.59, 244.90, 159.85),
    (25, 1.5535, 1.1690, None, 1.1690, 2, 251.00, 242.07, 157.77),
    (30, 1.5540, 1.1257, None, 1.1257, 2, 251.00, 242.07, 157.87),
    (35, 1.5556, 1.0919, None, 1.0919, 2, 248.17, 244.90, 161.64),
    (40, 1.5575, 1.0662, None, 1.0662, 2, 247.11, 244.54, 161.87),
    (45, 1.5606, 1.0483, None, 1.0483, 2, 245.35, 244.90, 162.91),
    (50, 1.5654, 1.0368, 1.0366, 1.0366, 3, 243.93, 244.90, 163.64),
    (60, 1.5610, 1.0699, None, 1.0699, 2, 245.35, 245.60, 163.62),
    (70, 1.5594, 1.0905, None, 1.0905, 2, 246.05, 246.31, 163.93),
    (80, 1.5608, 1.1062, None, 1.1062, 2, 245.31, 246.38, 164.30),
    (100, 1.5612, 1.1304, None, 1.1304, 2, 245.31, 244.96, 163.06),
    (130, 1.5612, 1.1550, None, 1.1550, 2, 245.31, 244.96, 163.06),
    (170, 1.5612, 1.1754, None, 1.1754, 2, 245.31, 244.96, 163.06),
    (220, 1.5612, 1.1890, None, 1.1890, 2, 245.31, 244.96, 163.06),
    (300, 1.5612, 1.1991, None, 1.1991, 2, 245.31, 244.96, 163.06),
    (400, 1.5612, 1.2044, None, 1.2044, 2, 245.31, 244.96, 163.06),
    (500, 1.5612, 1.2067, None, 1.2067, 2, 245.31, 244.96, 163.06),
)
#: The single-stage column of the page's comparison table, t = 5 to 500, with the
#: t = 0 instant the single-stage sweep also solves (the pool is full until day 5,
#: so it reads as t = 5 does: 1.5221 in the tutorial's run, the value its figure
#: plots at t = 0).
SINGLE_STAGE_T0 = 1.5221
SINGLE_STAGE = (1.5221, 1.4114, 1.3287, 1.2592, 1.2010, 1.1539, 1.1169, 1.0892,
                1.0698, 1.0578, 1.0965, 1.1218, 1.1419, 1.1743, 1.2091, 1.2404,
                1.2634, 1.2836, 1.2971, 1.3051)
#: The function's own figure size, and a narrower and a wider one.
SIZES = ((8.6, 5.0), (6.4, 4.8), (11.0, 4.6))


def _result(rows):
    """A ``fs_vs_time`` result from rows of dicts (``value``, ``fs`` and the rest)."""
    import pandas as pd
    df = pd.DataFrame(rows)
    df['param'] = 'time'
    df['method'] = 'spencer'
    df['success'] = df['fs'].notna()
    df['output'] = 'FS'
    df['output_label'] = 'Factor of Safety'
    ok = df.loc[df['success']]
    i = ok['fs'].astype(float).idxmin()
    return {'df': df, 'param': 'time', 'critical_time': float(df.at[i, 'value']),
            'min_fs': float(df.at[i, 'fs']), 'output': 'FS',
            'output_label': 'Factor of Safety'}


def _part1():
    return _result([{'value': t, 'fs': fs, 'Xo': x, 'Yo': y, 'R': r}
                    for t, fs, x, y, r in PART1])


def _part2():
    rows = []
    for t, s1, s2, s3, fs, g, x, y, r in PART2:
        rows.append({'value': float(t), 'fs': fs, 'Xo': x, 'Yo': y, 'R': r,
                     'stage1_FS': s1, 'stage2_FS': s2,
                     'stage3_FS': np.nan if s3 is None else s3, 'governs': g})
    return _result(rows)


def _single_stage():
    rows = [{'value': 0.0, 'fs': SINGLE_STAGE_T0}]
    rows += [{'value': float(t), 'fs': fs} for (t, *_rest), fs in zip(PART2, SINGLE_STAGE)]
    for r in rows:
        r.update(Xo=np.nan, Yo=np.nan, R=np.nan)
    return _result(rows)


def _covered(fig):
    """What the legend of ``fig``'s first axes covers: the texts and annotation
    arrows it overlaps, and whether it runs off the figure."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ax = fig.axes[0]
    leg = ax.get_legend()
    if leg is None:
        return None, ["no legend"]
    box = leg.get_window_extent(r)
    hits = []
    for a in fig.axes:
        for t in a.texts:
            if t.get_visible() and t.get_text() \
                    and box.overlaps(t.get_window_extent(r)):
                hits.append(repr(t.get_text()))
            arrow = getattr(t, 'arrow_patch', None)
            if arrow is not None and box.overlaps(arrow.get_window_extent(r)):
                hits.append(f"the arrow of {t.get_text()!r}")
    fb = fig.bbox
    if box.x0 < fb.x0 - 0.5 or box.x1 > fb.x1 + 0.5 or box.y0 < fb.y0 - 0.5 \
            or box.y1 > fb.y1 + 0.5:
        hits.append("the figure's edge")
    return box, hits


def check_the_combo03_figures(fails):
    """A: the legend covers no label on any of the three figures, at three sizes."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from xslope.fileio import load_slope_data
    from xslope.plot import plot_fs_vs_time

    for path in (COMBO03, COMBO03R):
        if not path.exists():
            fails.append(f"A: missing {path}")
            return
    sd1 = load_slope_data(str(COMBO03))
    sd2 = load_slope_data(str(COMBO03R))
    figures = (
        ("Part 1's curve", _part1, sd1, {}),
        ("Part 2's rapid drawdown curve", _part2, sd2, {}),
        ("Part 2's comparison with the single-stage curve", _part2, sd2,
         {'compare': _single_stage(), 'compare_label': 'single-stage, drained'}),
    )
    for blind in (False, True):
        with _text_blind_best(blind):
            for name, make, sd, kw in figures:
                for size in SIZES:
                    fig = plt.figure(figsize=size)
                    try:
                        plot_fs_vs_time(make(), slope_data=sd, fig=fig, **kw)
                        _box, hits = _covered(fig)
                    finally:
                        plt.close(fig)
                    if hits:
                        fails.append(
                            f"A: {name} at {size[0]:g} x {size[1]:g} in"
                            f"{' (loc=best blind to text)' if blind else ''}: "
                            f"the legend covers {', '.join(hits)}")


@contextlib.contextmanager
def _text_blind_best(on):
    """While ``on``, matplotlib's ``loc='best'`` does not see text, as it did not
    before matplotlib taught it to (the shipped COMBO-3 comparison figure of
    2026-11-13 was drawn that way: its legend sits on ``full pool, 1.470``). The
    figure's placement must not depend on which matplotlib draws it."""
    if not on:
        yield
        return
    import matplotlib.legend as mlegend
    from matplotlib.text import Text
    real = mlegend.Legend._auto_legend_data

    def blind(self, renderer):
        kids = self.parent._children
        self.parent._children = [a for a in kids if not isinstance(a, Text)]
        try:
            return real(self, renderer)
        finally:
            self.parent._children = kids

    mlegend.Legend._auto_legend_data = blind
    try:
        yield
    finally:
        mlegend.Legend._auto_legend_data = real


def check_no_free_corner(fails):
    """B: with every in-axes position covering a label, the legend goes below."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    try:
        from xslope.plot import _legend_clear_of_labels
    except ImportError as exc:
        fails.append(f"B: {exc}")
        return
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    try:
        for k in range(4):
            ax.plot([0, 1, 2], [k, k + 1, k], label=f"series {k + 1}")
        ax.set_xlabel("time (day)")
        # A label in every place a legend inside the axes could go.
        for x in (0.02, 0.5, 0.98):
            for y in (0.02, 0.5, 0.98):
                ax.annotate("a label " * 6, xy=(x, y), xycoords='axes fraction',
                            ha='center', va='center', fontsize=9)
        fig.tight_layout()
        h, l = ax.get_legend_handles_labels()
        _legend_clear_of_labels(ax, h, l, fontsize=9)
        box, hits = _covered(fig)
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        if box is not None:
            under = [t for t in list(ax.get_xticklabels()) + [ax.xaxis.get_label()]
                     if t.get_text() and box.overlaps(t.get_window_extent(r))]
            if under:
                hits.append("the x tick labels or the x-axis label")
            if box.y1 > ax.xaxis.get_label().get_window_extent(r).y0 + 0.5:
                hits.append("the space above the x-axis label (it is not below "
                            "the axes)")
        if hits:
            fails.append(f"B: with no free position inside the axes, the legend "
                         f"covers {', '.join(hits)}")
    finally:
        plt.close(fig)


def run():
    fails = []
    for check in (check_the_combo03_figures, check_no_free_corner):
        try:
            check(fails)
        except Exception as exc:                             # pragma: no cover
            fails.append(f"{check.__name__} raised {type(exc).__name__}: {exc}")
    return fails


if __name__ == '__main__':
    failures = run()
    for f in failures:
        print("FAIL", f)
    print("fs_vs_time_legend:", "FAIL" if failures else "PASS")
    sys.exit(1 if failures else 0)
