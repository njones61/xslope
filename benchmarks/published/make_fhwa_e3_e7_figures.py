"""Figures for the FHWA Examples E3 to E7 entries of docs/verification/published.md.

Two PNGs per example into docs/verification/images/:

  fhwa_e<n>.png            the input model -- soil zones, reinforcement layers,
                           surcharges
  fhwa_e<n>_solution.png   the critical Spencer surface the search returns

The input figure is ``plot_inputs`` on the committed model, so what it draws is
exactly what the locks read.  The solution figure re-runs the same search the
page's ``circular_search`` tag runs.  Both are framed to the model
(``frame="content"``) rather than to the whole domain: these walls sit inside an
opened-out foundation block, and a full-domain frame shrinks the wall to a
sliver.

    PYTHONPATH=. python3 benchmarks/published/make_fhwa_e3_e7_figures.py
    PYTHONPATH=. python3 benchmarks/published/make_fhwa_e3_e7_figures.py e3 e7
"""
import contextlib
import io
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, '..', '..'))
sys.path.insert(0, _REPO)

import matplotlib                                                    # noqa: E402
matplotlib.use('Agg')
import matplotlib.pyplot as plt                                      # noqa: E402

from xslope.fileio import load_slope_data                            # noqa: E402
from xslope.plot import plot_inputs, plot_solution                   # noqa: E402
from xslope.search import circular_search                            # noqa: E402

FILES = os.path.join(_REPO, 'docs', 'verification', 'files', 'published')
IMAGES = os.path.join(_REPO, 'docs', 'verification', 'images')
DPI = 150
METHOD = 'spencer'
NUM_SLICES = 30
EXAMPLES = ('e3', 'e4', 'e5', 'e6', 'e7')


def capture(path, fn, *args, **kwargs):
    """Run a production plot function and save whatever it draws.

    The plot functions end on ``plt.show()``; with a non-interactive backend that
    warns and throws the figure away, so ``show`` is briefly replaced by a save.
    """
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
    print('wrote', os.path.relpath(path, _REPO))


def make(tag):
    slope_data = load_slope_data(os.path.join(FILES, 'fhwa_%s.xlsx' % tag))
    capture(os.path.join(IMAGES, 'fhwa_%s.png' % tag), plot_inputs,
            slope_data, mode='lem', frame='content')

    with contextlib.redirect_stdout(io.StringIO()):
        fs_cache, _conv, _path, _cc = circular_search(
            slope_data, METHOD, num_slices=NUM_SLICES)
    crit = fs_cache[0]
    capture(os.path.join(IMAGES, 'fhwa_%s_solution.png' % tag), plot_solution,
            slope_data, crit['slices'], crit['failure_surface'],
            crit['solver_result'], frame='content')
    print('   %s: Spencer FS = %.4f' % (tag, crit['solver_result']['FS']))


def main(argv):
    tags = [a.lower() for a in argv] or list(EXAMPLES)
    os.makedirs(IMAGES, exist_ok=True)
    for tag in tags:
        make(tag)


if __name__ == '__main__':
    main(sys.argv[1:])
