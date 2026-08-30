"""Figure for the FHWA Example E1 entry of docs/verification/published.md.

One PNG into docs/verification/images/:

  fhwa_e1.png    the input model — the three soil zones, the eleven geogrid
                 layers, the broken backslope and the traffic surcharge

Nothing is solved here.  The figure is `plot_inputs` on the committed model, so
what it draws is exactly what the locks read.

    PYTHONPATH=. python3 benchmarks/published/make_fhwa_e1_figures.py
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
from xslope.plot import plot_inputs                                  # noqa: E402

MODEL = os.path.join(_REPO, 'docs', 'verification', 'files', 'published',
                     'fhwa_e1.xlsx')
OUT = os.path.join(_REPO, 'docs', 'verification', 'images', 'fhwa_e1.png')
DPI = 150


def main():
    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(MODEL)
        plot_inputs(sd, mode='lem')
    fig = plt.gcf()
    fig.savefig(OUT, dpi=DPI, bbox_inches='tight')
    plt.close('all')
    print('wrote', OUT)


if __name__ == '__main__':
    main()
