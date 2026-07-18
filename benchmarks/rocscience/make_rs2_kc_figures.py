"""Render the side-by-side (inputs | critical seismic surface) figures for the
RS2 #68 critical-seismic-coefficient rows on docs/verification/rs2.md.

RS2 #68 is the one LEM problem on the RS2 page whose lock is a critical seismic
coefficient k꜀ (FS = 1), not a factor of safety, so it is not covered by the
`fem_ssrm` 4-panel renderer (make_rs2_figures.py). Each figure sets the seismic
coefficient to the locked k꜀ from the problem's `critical_kc` tag — the same value
the regression uses — re-runs the circular search, and shows the critical seismic
surface (with the pseudo-static arrow) beside the inputs.

Run from the repo root:
    python benchmarks/rocscience/make_rs2_kc_figures.py            # all
    python benchmarks/rocscience/make_rs2_kc_figures.py rs2_68b    # a subset (by file stem)
"""

import io
import os
import re
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
from xslope.search import circular_search
from xslope.plot import plot_inputs, plot_solution

ROOT = os.path.join(os.path.dirname(__file__), '..', '..')
RS2_MD = os.path.join(ROOT, 'docs', 'verification', 'rs2.md')
OUT = os.path.join(ROOT, 'docs', 'verification', 'images')

TAG_RE = re.compile(r'<!--\s*test:\s*(.*?)\s*-->')


def parse_kc_tags(path=RS2_MD):
    """One entry per input file: the Spencer `critical_kc` tag (the rigorous
    method), carrying the locked k꜀ and the search settings."""
    by_file = {}
    with open(path) as fh:
        for line in fh:
            m = TAG_RE.search(line)
            if not m:
                continue
            kv = {}
            for part in m.group(1).split(','):
                if '=' in part:
                    k, v = part.split('=', 1)
                    kv[k.strip()] = v.strip()
            if kv.get('type') != 'critical_kc':
                continue
            stem = os.path.basename(kv['file']).split('.')[0]
            # prefer Spencer; fall back to whatever method appears first
            if stem not in by_file or kv.get('method') == 'spencer':
                by_file[stem] = kv
    return by_file


def make_figure(stem, tag, panel_size=(8.0, 5.0), dpi=150):
    path = os.path.normpath(os.path.join(ROOT, 'docs', 'verification', tag['file']))
    sd = load_slope_data(path)
    sd['k_seismic'] = float(tag['expected_kc'])          # the locked k꜀ (FS = 1)
    method = tag.get('method', 'spencer')
    ns = int(tag.get('num_slices', 40))

    with contextlib.redirect_stdout(io.StringIO()):
        cache, _conv, _path, _cc = circular_search(sd, method, num_slices=ns)
    best = cache[0]
    df, surface, results = best['slices'], best['failure_surface'], best['solver_result']

    paths = []
    for which in ('inputs', 'solution'):
        fig = plt.figure(figsize=panel_size)
        if which == 'inputs':
            plot_inputs(sd, fig=fig, mat_table=False, show_title=True,
                        title=f'{stem} — inputs', label_coordinates=False)
        else:
            plot_solution(sd, df, surface, results, fig=fig, show_title=True)
        p = os.path.join(OUT, f'_{stem}_{which}.png')
        fig.savefig(p, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        paths.append(p)

    imgs = [Image.open(p) for p in paths]
    h = min(im.height for im in imgs)
    imgs = [im.resize((round(im.width * h / im.height), h), Image.LANCZOS) for im in imgs]
    combo = Image.new('RGB', (sum(im.width for im in imgs) + 20, h), 'white')
    x = 0
    for im in imgs:
        combo.paste(im, (x, 0))
        x += im.width + 20
    out = os.path.join(OUT, f'{stem}.png')
    combo.save(out)
    for p in paths:
        os.remove(p)
    return out, best['FS']


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    only = set(sys.argv[1:])
    tags = parse_kc_tags()
    for stem, tag in sorted(tags.items()):
        if only and stem not in only:
            continue
        out, fs = make_figure(stem, tag)
        print(f'ok   {stem}  k꜀={float(tag["expected_kc"]):.3f}  FS={fs:.3f}  '
              f'{os.path.basename(out)}', flush=True)
