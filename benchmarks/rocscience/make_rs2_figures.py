"""Render the side-by-side (FEM inputs | SSRM shear strain) figures for the RS2
SSR corpus page (docs/verification/rs2.md).

Two panels per problem, matching what the SSR corpus is actually about:

  left  — plot_fem_data: the mesh, coloured by material zone, with the
          displacement/boundary conditions drawn.
  right — plot_fem_results('shear_strain'): the viscoplastic max shear strain at
          the critical F. Colour fill only (no mesh lines) so the failure
          mechanism reads at panel size.

The cases are parsed straight out of the ``fem_ssrm`` test tags in rs2.md rather
than kept in a second list here. That is deliberate: the figure is then rendered
on the SAME mesh and F-bracket the regression lock uses, so a retagged mesh size
cannot silently leave a stale figure behind.

Each SSRM solve costs about a minute, so a full run is slow. Pass benchmark ids
to render a subset.

Run from the repo root:
    python benchmarks/rocscience/make_rs2_figures.py            # all
    python benchmarks/rocscience/make_rs2_figures.py RS2-30 RS2-31a
"""

import io
import os
import re
import sys
import time
import warnings
import contextlib

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

from xslope.fileio import load_slope_data
from xslope.fem import build_fem_data, solve_ssrm
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints)
from xslope.plot_fem import plot_fem_data, plot_fem_results

ROOT = os.path.join(os.path.dirname(__file__), '..', '..')
RS2_MD = os.path.join(ROOT, 'docs', 'verification', 'rs2.md')
OUT = os.path.join(ROOT, 'docs', 'verification', 'images')

TAG_RE = re.compile(r'<!--\s*test:\s*(.*?)\s*-->')


def parse_tags(path=RS2_MD):
    """Every fem_ssrm tag on the page, as a list of key->value dicts. The tag is
    the single source of truth for mesh size, element type and F-bracket."""
    cases = []
    with open(path) as fh:
        for line in fh:
            m = TAG_RE.search(line)
            if not m:
                continue
            kv = {}
            for part in m.group(1).split(','):
                if '=' not in part:
                    continue
                k, v = part.split('=', 1)
                kv[k.strip()] = v.strip()
            if kv.get('type') != 'fem_ssrm':
                continue
            cases.append(kv)
    return cases


def _build(tag):
    """Mesh + fem_data exactly as run_tests.run_fem_test does."""
    # tag paths are relative to docs/verification/
    path = os.path.normpath(os.path.join(ROOT, 'docs', 'verification', tag['file']))
    sd = load_slope_data(path)

    # seepage-coupled models must reuse the stored mesh (nodal seep_u)
    uses_seep = any(str(m.get('u', '')).strip().lower() == 'seep'
                    for m in sd.get('materials', []))
    if uses_seep and sd.get('mesh') is not None:
        mesh = sd['mesh']
    else:
        target = tag.get('target_size')
        if target is None:
            xs = [x for x, _ in sd['ground_surface'].coords]
            target = (max(xs) - min(xs)) / 100
        else:
            target = float(target)
        lines, _nr, _np = extract_constraint_line_geometry(sd)
        polys = get_material_polygons(sd, reinf_lines=lines)
        mesh = build_mesh_from_polygons(
            polys, target_size=target,
            element_type=tag.get('element_type', 'tri6'), lines=lines,
            point_constraints=extract_point_constraints(sd))
    return sd, build_fem_data(sd, mesh)


def make_figure(tag, panel_size=(8.0, 5.0), dpi=150):
    bench = tag.get('benchmark', os.path.basename(tag['file']).split('.')[0])
    sd, fem_data = _build(tag)

    with contextlib.redirect_stdout(io.StringIO()):
        sol = solve_ssrm(fem_data,
                         F_min=float(tag.get('f_min', 0.5)),
                         F_max=float(tag.get('f_max', 3.0)),
                         tolerance=float(tag.get('tolerance', 0.02)),
                         max_iterations=int(tag.get('max_iter', 4000)),
                         debug_level=0)
    if not sol.get('converged'):
        raise RuntimeError(f'SSRM did not converge: {sol.get("error")}')

    # solve_ssrm returns the bracket result; the FIELD solution (displacements,
    # viscoplastic strains) is the last converged solve, at the F just below
    # critical -- that is the developed mechanism worth plotting.
    field = sol.get('last_solution')
    if field is None:
        raise RuntimeError('SSRM returned no last_solution to plot')

    paths = []
    # left: mesh + materials + boundary conditions
    fig = plt.figure(figsize=panel_size)
    plot_fem_data(fem_data, fig=fig, show_title=True)
    p = os.path.join(OUT, f'_{bench}_data.png')
    fig.savefig(p, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    paths.append(p)

    # right: viscoplastic shear strain at the critical F, colour fill only
    fig = plt.figure(figsize=panel_size)
    plot_fem_results(fem_data, field, plot_type=['shear_strain'],
                     show_mesh=False, show_reinforcement=True,
                     fig=fig, show_title=True, show_legend=False)
    p = os.path.join(OUT, f'_{bench}_strain.png')
    fig.savefig(p, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    paths.append(p)

    imgs = [Image.open(p) for p in paths]
    h = min(im.height for im in imgs)
    imgs = [im.resize((int(im.width * h / im.height), h)) for im in imgs]
    combo = Image.new('RGB', (sum(im.width for im in imgs) + 20, h), 'white')
    x = 0
    for im in imgs:
        combo.paste(im, (x, 0))
        x += im.width + 20
    out = os.path.join(OUT, f'{bench}.png')
    combo.save(out)
    for p in paths:
        os.remove(p)
    return out, sol['FS']


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    only = set(sys.argv[1:])
    cases = parse_tags()
    print(f'{len(cases)} fem_ssrm tags on rs2.md')
    for tag in cases:
        bench = tag.get('benchmark', '?')
        if only and bench not in only:
            continue
        t0 = time.time()
        try:
            out, fs = make_figure(tag)
            print(f'ok   {bench:10s} FS={fs:.3f}  ({time.time()-t0:.0f}s)  '
                  f'{os.path.basename(out)}', flush=True)
        except Exception as e:
            print(f'FAIL {bench:10s} {type(e).__name__}: {e}', flush=True)
