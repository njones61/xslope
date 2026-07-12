"""Render the side-by-side (inputs | solution) figures for the verification
corpus pages (docs/verification/rocscience.md).

For every built corpus file this renders plot_inputs (with coordinate labels,
verification-manual style) and a representative solved surface, then stitches
the two panels into one PNG in docs/verification/images/. Search problems are
re-solved here so the figure shows the critical surface, not the stored seed.

Run from the repo root:  python benchmarks/rocscience/make_figures.py [vpNNN ...]
"""

import io
import os
import sys
import warnings
import contextlib

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

from xslope.fileio import load_slope_data
from xslope.slice import generate_slices
from xslope import solve
from xslope.plot import plot_inputs, plot_solution

SRC = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'files', 'rocscience')

# Talbingo (vp005/vp006) has 26 vertices, half of them within a few metres at
# the crest - no readable label layout exists at panel size, so those two
# figures skip the coordinate labels (the manual's Table 5.2 carries them).
NO_COORD_LABELS = {'vp005', 'vp006'}
# Geometries with tight vertex clusters where the leader-line layout helps;
# everywhere else labels sit plainly beside their vertices (coord_arrows=False).
COORD_ARROWS = {'vp051'}
OUT = os.path.join(os.path.dirname(__file__), '..', '..', 'docs', 'verification', 'images')

# (file stem, surface kind, method) — kind: circle = the stored circle is the
# critical one (single_circle problems and searches whose critical circle was
# stored back); csearch/ncsearch = re-run the search; noncirc = the stored
# non-circular surface; rapid = rapid drawdown on the stored circle / search.
CASES = [
    ('xslope_acads_simple', 'csearch', 'bishop', 'lem'),
    ('xslope_arai_tagyo', 'csearch', 'bishop', 'lem'),
    ('xslope_acads_weak_layer', 'ncsearch', 'spencer', 'lem'),
    ('vp002', 'csearch', 'bishop'),
    ('vp003', 'csearch', 'spencer'),
    ('vp004', 'csearch', 'spencer'),
    ('vp005', 'circle', 'spencer'),
    ('vp006', 'circle', 'spencer'),
    ('vp008', 'noncirc', 'spencer'),
    ('vp009', 'ncsearch', 'spencer'),
    ('vp015', 'csearch', 'spencer'),
    ('vp016', 'csearch', 'spencer'),
    ('vp017', 'csearch', 'spencer'),
    ('vp018', 'ncsearch', 'spencer'),
    ('vp019', 'csearch', 'spencer'),
    ('vp020', 'csearch', 'spencer'),
    ('vp021a', 'circle', 'spencer'),
    ('vp021b', 'circle', 'spencer'),
    ('vp023', 'csearch', 'bishop'),
    ('vp024', 'csearch', 'bishop'),
    ('vp025', 'noncirc', 'spencer'),
    ('vp027', 'circle', 'spencer'),
    ('vp036', 'csearch', 'bishop'),
    ('vp041', 'csearch', 'bishop'),
    ('vp044a', 'csearch', 'spencer'),
    ('vp044b', 'csearch', 'spencer'),
    ('vp044c', 'csearch', 'spencer'),
    ('vp045a', 'csearch', 'spencer'),
    ('vp045b', 'csearch', 'spencer'),
    ('vp047', 'noncirc', 'janbu'),
    ('vp048', 'noncirc', 'janbu'),
    ('vp050', 'noncirc', 'janbu'),
    ('vp051', 'circle', 'spencer'),
    ('vp052a', 'csearch', 'spencer'),
    ('vp052b', 'csearch', 'spencer'),
    ('vp054a', 'circle', 'bishop'),
    ('vp054b', 'circle', 'bishop'),
    ('vp061a', 'csearch', 'spencer'),
    ('vp061b', 'csearch', 'spencer'),
    ('vp062a', 'csearch', 'spencer'),
    ('vp062b', 'csearch', 'spencer'),
    ('vp064', 'circle', 'spencer'),
    ('vp065', 'circle', 'bishop'),
    ('vp066', 'circle', 'spencer'),
    ('vp067', 'circle', 'spencer'),
    ('vp068', 'circle', 'bishop'),
    ('vp070a', 'csearch', 'spencer'),
    ('vp074', 'csearch', 'spencer'),
    ('vp078', 'csearch', 'bishop'),
    ('vp079', 'csearch', 'spencer'),
    ('vp080a', 'circle', 'spencer'),
    ('vp080b', 'circle', 'spencer'),
    ('vp081', 'csearch', 'spencer'),
    ('vp085a', 'circle', 'bishop'),
    ('vp085b', 'circle', 'bishop'),
    ('vp086', 'csearch', 'spencer'),
    ('vp087', 'circle', 'bishop'),
    ('vp088', 'circle', 'spencer'),
    ('vp089', 'circle', 'spencer'),
    ('vp090', 'circle', 'bishop'),
    ('vp091', 'circle', 'spencer'),
    ('vp092', 'circle', 'bishop'),
    ('vp093', 'circle', 'bishop'),
    ('vp094', 'circle', 'bishop'),
    ('vp096', 'rapid_circle', 'spencer'),
    ('vp097', 'rapid_csearch', 'spencer'),
    ('vp098', 'rapid_csearch', 'spencer'),
    ('vp099', 'rapid_csearch', 'spencer'),
    ('vp100', 'csearch', 'bishop'),
    ('vp101', 'csearch', 'bishop'),
]


def _solve_case(sd, kind, method):
    """Return (slice_df, failure_surface, results) for the representative surface."""
    fn = getattr(solve, method)
    if kind in ('circle', 'rapid_circle'):
        ok, res = generate_slices(sd, circle=sd['circles'][0], num_slices=60)
        if not ok:
            raise RuntimeError(res)
        df, fs = res
        if kind == 'rapid_circle':
            from xslope.advanced import rapid_drawdown
            s, r = rapid_drawdown(df, method, debug_level=0)
        else:
            s, r = fn(df)
        if not s:
            raise RuntimeError(r)
        return df, fs, r
    if kind == 'noncirc':
        ok, res = generate_slices(sd, non_circ=sd['non_circ'], num_slices=60)
        if not ok:
            raise RuntimeError(res)
        df, fs = res
        s, r = fn(df)
        if not s:
            raise RuntimeError(r)
        return df, fs, r
    if kind in ('csearch', 'rapid_csearch'):
        from xslope.search import circular_search
        cache, _, _, _ = circular_search(sd, method, num_slices=50,
                                         rapid=(kind == 'rapid_csearch'))
        best = cache[0]
        return best['slices'], best['failure_surface'], best['solver_result']
    if kind == 'ncsearch':
        from xslope.search import noncircular_search
        cache, _, _ = noncircular_search(sd, method, num_slices=50)
        best = cache[0]
        return best['slices'], best['failure_surface'], best['solver_result']
    raise ValueError(kind)


def make_figure(stem, kind, method, src=None, panel_size=(8.0, 5.0), dpi=150):
    base = SRC if src is None else os.path.join(os.path.dirname(__file__), '..', '..', 'docs', src, 'files')
    sd = load_slope_data(os.path.join(base, f'{stem}.xlsx'))
    with contextlib.redirect_stdout(io.StringIO()):
        df, fs, results = _solve_case(sd, kind, method)

    paths = []
    for which in ('inputs', 'solution'):
        fig = plt.figure(figsize=panel_size)
        if which == 'inputs':
            plot_inputs(sd, fig=fig, mat_table=False, show_title=True,
                        title=f'{stem} — inputs',
                        label_coordinates=stem not in NO_COORD_LABELS,
                        coord_arrows=stem in COORD_ARROWS,
                        coord_label_size=6)
        else:
            plot_solution(sd, df, fs, results, fig=fig, show_title=True)
        p = os.path.join(OUT, f'_{stem}_{which}.png')
        fig.savefig(p, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        paths.append(p)

    # stitch horizontally at equal height
    imgs = [Image.open(p) for p in paths]
    h = min(im.height for im in imgs)
    imgs = [im.resize((int(im.width * h / im.height), h)) for im in imgs]
    combo = Image.new('RGB', (sum(im.width for im in imgs) + 20, h), 'white')
    x = 0
    for im in imgs:
        combo.paste(im, (x, 0))
        x += im.width + 20
    out = os.path.join(OUT, f'{stem}.png')
    combo.save(out)
    for p in paths:
        os.remove(p)
    return out


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    only = set(sys.argv[1:])
    for case in CASES:
        stem, kind, method = case[:3]
        src = case[3] if len(case) > 3 else None
        if only and stem not in only:
            continue
        try:
            out = make_figure(stem, kind, method, src=src)
            print('ok ', stem)
        except Exception as e:
            print('FAIL', stem, ':', repr(e)[:100])
