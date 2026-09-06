"""Render the figures for the Hassiotis et al. (1997) sample on docs/lem/samples.md.

Four figures into docs/lem/sample_images/, all through the production plot path,
framed to their content the way the rest of that page is framed (equal aspect,
the axes box at the model's true proportions, one uniform cushion):

  hassiotis_inputs.png      plot_inputs: the slope, the two pile stations drawn
                            together so one picture places both rows.
  hassiotis_results.png     plot_solution, unreinforced, Spencer.
  hassiotis_p1_results.png  plot_solution, pile row 13.7 m from the toe.
  hassiotis_p2_results.png  plot_solution, pile row 23.1 m from the toe.

The searches read the model's own search window (the entry/exit ranges and
tangent limit the two pile files declare), so the surface drawn is the surface
the page's locks are taken on rather than the deep bypass an unwindowed search
finds.  Method and slice count match the page's test tags.

Run from the repo root:  PYTHONPATH=. python3 benchmarks/make_hassiotis_figures.py
"""
import contextlib
import io
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from xslope.fileio import load_slope_data
from xslope.plot import plot_inputs, plot_solution
from xslope.search import circular_search, file_search_window

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
FILES = os.path.join(ROOT, 'docs', 'lem', 'files')
IMG = os.path.join(ROOT, 'docs', 'lem', 'sample_images')
METHOD = 'spencer'
NUM_SLICES = 50
DPI = 200


def capture(path, fn, *args, **kwargs):
    saved = []

    def _show(*a, **k):
        plt.gcf().savefig(path, dpi=DPI, bbox_inches='tight')
        saved.append(path)
        plt.close('all')

    orig = plt.show
    plt.show = _show
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn(*args, **kwargs)
        if not saved:
            plt.gcf().savefig(path, dpi=DPI, bbox_inches='tight')
            plt.close('all')
    finally:
        plt.show = orig
    print('wrote', os.path.relpath(path, ROOT))


def solution(stem, out):
    sd = load_slope_data(os.path.join(FILES, stem + '.xlsx'))
    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _c, _p, _cc = circular_search(
            sd, METHOD, num_slices=NUM_SLICES, diagnostic=False,
            **file_search_window(sd))
    crit = fs_cache[0]
    capture(os.path.join(IMG, out), plot_solution, sd, crit['slices'],
            crit['failure_surface'], crit['solver_result'], frame='content')
    print('   %s: Spencer FS = %.3f' % (stem, crit['solver_result']['FS']))


def main():
    # The inputs panel carries both pile stations at once: the two rows are the
    # same model with the row moved, and one picture places them relative to the
    # slope better than two identical pictures with a line in a different place.
    sd = load_slope_data(os.path.join(FILES, 'xslope_hassiotis_p1.xlsx'))
    sd2 = load_slope_data(os.path.join(FILES, 'xslope_hassiotis_p2.xlsx'))
    sd['pile_lines'] = sd['pile_lines'] + sd2['pile_lines']
    capture(os.path.join(IMG, 'hassiotis_inputs.png'), plot_inputs, sd,
            mode='lem', frame='content')

    solution('xslope_hassiotis', 'hassiotis_results.png')
    solution('xslope_hassiotis_p1', 'hassiotis_p1_results.png')
    solution('xslope_hassiotis_p2', 'hassiotis_p2_results.png')


if __name__ == '__main__':
    main()
