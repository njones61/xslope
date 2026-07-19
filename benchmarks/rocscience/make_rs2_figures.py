"""Render the 4-panel FEM/SSRM figures for the RS2 SSR corpus page
(docs/verification/rs2.md).

Four panels per problem, arranged 2x2, giving the reader the inputs, the mesh,
and the two SSRM field solutions at a glance:

  upper-left  — plot_inputs(mode="fem"): the slope geometry with the FEM-relevant
                overlays (material zones, water table / piezo line, reinforcement,
                loads). This is the "what went in" panel.
  lower-left  — plot_fem_data: the mesh, coloured by material zone, with the
                displacement/boundary conditions drawn.
  upper-right — plot_fem_results('shear_strain'): the viscoplastic max shear
                strain at the critical F. Colour fill only (no mesh lines) so the
                failure mechanism reads at panel size.
  lower-right — plot_fem_results('displace_vector'): the viscoplastic displacement
                vectors at the critical F, over the mesh boundary outline — the
                kinematics of the same failure mechanism.

The two right panels share the same field solution (the last converged SSRM
solve, at the F just below critical), so the strain contours and the displacement
arrows are two views of one mechanism.

The cases are parsed straight out of the ``fem_ssrm`` test tags in rs2.md rather
than kept in a second list here. That is deliberate: the figure is then rendered
on the SAME mesh and F-bracket the regression lock uses, so a retagged mesh size
cannot silently leave a stale figure behind.

Each SSRM solve costs about a minute (more on a fine mesh), so a full run is slow.
Pass benchmark ids to render a subset.

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

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

from xslope.fileio import load_slope_data
from xslope.fem import build_fem_data, solve_ssrm
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints)
from xslope.plot import plot_inputs
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
        # Feature-aware mesh refinement, matching run_tests._refine_kwargs: a tag
        # carrying refine_factor=N (and optional refine_features=a;b) must build the
        # SAME refined mesh the lock was recorded on, or the figure would show the
        # unrefined coarse mesh and a different FS.
        refine_kw = {}
        if tag.get('refine_factor') is not None:
            refine_kw['refine_factor'] = float(tag['refine_factor'])
            feats = tag.get('refine_features')
            if feats:
                refine_kw['refine_features'] = [s.strip() for s in str(feats).split(';')
                                                if s.strip()]
        mesh = build_mesh_from_polygons(
            polys, target_size=target,
            element_type=tag.get('element_type', 'tri6'), lines=lines,
            point_constraints=extract_point_constraints(sd), **refine_kw)
    return sd, build_fem_data(sd, mesh)


def build_and_solve(tag):
    """Build the mesh, run the SSRM bracket, and return the pieces the figure
    needs: (sd, fem_data, field, FS).

    ``field`` is the FIELD solution (displacements, viscoplastic strains) from the
    last converged solve — the F just below critical, i.e. the developed
    mechanism worth plotting — not the bracket midpoint.
    """
    sd, fem_data = _build(tag)

    # SSR-exclusion material names (semicolon-separated within the tag value,
    # since tag key=value pairs are comma-split) — held at full strength.
    ssr_exclude = None
    if tag.get('ssr_exclude'):
        ssr_exclude = [s.strip() for s in str(tag['ssr_exclude']).split(';') if s.strip()]

    # SSR search-area polygon (RS2's "SSR Search Area"): flat x1;y1;x2;y2;... list.
    ssr_zone = None
    if tag.get('ssr_zone'):
        _z = [float(v) for v in str(tag['ssr_zone']).split(';') if v.strip() != '']
        ssr_zone = list(zip(_z[0::2], _z[1::2]))

    with contextlib.redirect_stdout(io.StringIO()):
        sol = solve_ssrm(fem_data,
                         F_min=float(tag.get('f_min', 0.5)),
                         F_max=float(tag.get('f_max', 3.0)),
                         tolerance=float(tag.get('tolerance', 0.02)),
                         max_iterations=int(tag.get('max_iter', 4000)),
                         ssr_exclude=ssr_exclude,
                         ssr_zone=ssr_zone,
                         debug_level=0)
    if not sol.get('converged'):
        raise RuntimeError(f'SSRM did not converge: {sol.get("error")}')

    field = sol.get('last_solution')
    if field is None:
        raise RuntimeError('SSRM returned no last_solution to plot')
    return sd, fem_data, field, sol['FS']


def _compose_2x2(paths, gutter_frac=0.02):
    """Stitch the four panel PNGs (keyed 'UL','UR','LL','LR') into a 2x2 grid.

    Each panel is rendered independently with its own tight_layout, so the sizing
    here is deliberately self-adjusting rather than hand-tuned:

      * every panel is scaled to a common column width (the narrowest panel's
        width, so nothing is upscaled and blurred) — this aligns the two columns;
      * each row's height is the taller of its two panels, and the shorter panel
        is vertically centred in the row — this aligns the two rows and keeps the
        left column (which carries a legend below the axes) flush with the right;
      * the gutter between panels is a fraction of the column width, so it tracks
        the render resolution instead of being a fixed pixel constant.
    """
    ims = {k: Image.open(p) for k, p in paths.items()}
    W = min(im.width for im in ims.values())

    def fit(im):
        if im.width == W:
            return im
        return im.resize((W, round(im.height * W / im.width)), Image.LANCZOS)

    ims = {k: fit(v) for k, v in ims.items()}
    row0 = max(ims['UL'].height, ims['UR'].height)
    row1 = max(ims['LL'].height, ims['LR'].height)
    gut = max(1, round(gutter_frac * W))

    combo = Image.new('RGB', (2 * W + gut, row0 + row1 + gut), 'white')

    def place(key, col, row_top, row_h):
        im = ims[key]
        x = col * (W + gut)
        y = row_top + (row_h - im.height) // 2
        combo.paste(im, (x, y))

    place('UL', 0, 0, row0)
    place('UR', 1, 0, row0)
    place('LL', 0, row0 + gut, row1)
    place('LR', 1, row0 + gut, row1)
    return combo


def render_figure(bench, sd, fem_data, field, out_dir=OUT, panel_size=(8.0, 5.0), dpi=150):
    """Render the four panels and compose them into <out_dir>/<bench>.png.

    Pure plotting: no solve happens here, so a cached (sd, fem_data, field) can be
    re-rendered cheaply while iterating on the layout.
    """
    os.makedirs(out_dir, exist_ok=True)
    tmp = {}

    def _panel(key, draw, figsize=panel_size):
        fig = plt.figure(figsize=figsize)
        draw(fig)
        p = os.path.join(out_dir, f'_{bench}_{key}.png')
        fig.savefig(p, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        tmp[key] = p

    # The three FEM panels box-adjust to a tight wide strip (cropped by the tight
    # bbox), but plot_inputs uses adjustable='datalim' and so pads the data range
    # out to whatever figure aspect it is given. Size its figure to the actual
    # slope aspect so it crops to the same strip as the others instead of leaving
    # a tall band of dead space above and below the geometry.
    nodes = fem_data['nodes']
    dw = float(nodes[:, 0].max() - nodes[:, 0].min())
    dh = float(nodes[:, 1].max() - nodes[:, 1].min())
    aspect = (dh / dw) if dw > 0 else 0.5
    ul_h = float(np.clip(panel_size[0] * aspect, 2.0, panel_size[1]))

    # upper-left: FEM inputs (geometry, material zones, water table, reinforcement)
    _panel('UL', lambda fig: plot_inputs(
        sd, mode='fem', title='FEM Inputs', fig=fig,
        show_title=True, show_legend=True), figsize=(panel_size[0], ul_h))
    # lower-left: mesh + material zones + boundary conditions
    _panel('LL', lambda fig: plot_fem_data(fem_data, fig=fig, show_title=True))
    # upper-right: viscoplastic shear strain at the critical F, colour fill only
    _panel('UR', lambda fig: plot_fem_results(
        fem_data, field, plot_type=['shear_strain'],
        show_mesh=False, show_reinforcement=True,
        fig=fig, show_title=True, show_legend=False))
    # lower-right: viscoplastic displacement vectors at the critical F
    _panel('LR', lambda fig: plot_fem_results(
        fem_data, field, plot_type=['displace_vector'],
        show_reinforcement=True,
        fig=fig, show_title=True, show_legend=False))

    combo = _compose_2x2(tmp)
    out = os.path.join(out_dir, f'{bench}.png')
    combo.save(out)
    for p in tmp.values():
        os.remove(p)
    return out


def make_figure(tag, panel_size=(8.0, 5.0), dpi=150):
    bench = tag.get('benchmark', os.path.basename(tag['file']).split('.')[0])
    sd, fem_data, field, fs = build_and_solve(tag)
    out = render_figure(bench, sd, fem_data, field, panel_size=panel_size, dpi=dpi)
    return out, fs


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
