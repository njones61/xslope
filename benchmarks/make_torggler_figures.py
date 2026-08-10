"""Render the mesh and SSRM results figures for the Torggler (2016) cases on
docs/verification/ssrm.md.

Two figures per case, both through the production plot path at its defaults:

  ``<prefix>_mesh.png``     — plot_fem_data: the mesh coloured by material zone
                              with the displacement boundary conditions and the
                              plate drawn.
  ``<prefix>_results.png``  — plot_fem_results at its default three panels,
                              rendered from the at-failure field solve_ssrm
                              captures, with the locked FS in the titles.

Like the Griffiths generator, the cases come from the page's own ``fem_ssrm``
tags rather than a second list here, so a figure is always drawn on the mesh and
F-bracket its lock runs, and a retagged mesh size cannot leave a stale figure
behind.

Each solve costs several minutes. Pass prefixes to render a subset.

Run from the repo root:
    PYTHONPATH=. python3 benchmarks/make_torggler_figures.py
    PYTHONPATH=. python3 benchmarks/make_torggler_figures.py torggler_3a_plate
"""
import contextlib
import io
import os
import sys
import time
import warnings

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import run_tests as RT
from xslope.fileio import load_slope_data
from xslope.fem import build_fem_data, solve_ssrm
import xslope.fem as _fem
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints,
                         extract_size_regions)
from xslope.plot_fem import plot_fem_data, plot_fem_results

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
SSRM_MD = os.path.join(ROOT, 'docs', 'verification', 'ssrm.md')
OUT = os.path.join(ROOT, 'docs', 'fem', 'images')
DPI = 200

# One mesh figure per section, drawn on the section's supported case: it is the
# unsupported mesh plus the plate, so it shows both models' discretization and
# where the plate sits in it.
MESH_FIGURE = {'xslope_torggler_3a_plate': 'torggler_3a',
               'xslope_torggler_3b_plate': 'torggler_3b'}


def figure_tags():
    """Every Torggler ``fem_ssrm`` tag on the ssrm page, in page order."""
    return [t for t in RT.parse_test_tags(SSRM_MD)
            if t.get('type') == 'fem_ssrm' and 'torggler' in t.get('file', '')]


def stem(tag):
    return os.path.splitext(os.path.basename(tag['file']))[0]


def render(tag):
    name = stem(tag)
    t0 = time.time()
    sd = load_slope_data(tag['file'])
    lines, _nr, _np = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            polys, target_size=tag.get('target_size'),
            element_type=tag['element_type'], lines=lines,
            point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd), **RT._refine_kwargs(tag))
        fem_data = build_fem_data(sd, mesh)

    if name in MESH_FIGURE:
        with contextlib.redirect_stdout(io.StringIO()):
            plot_fem_data(fem_data)
        path = os.path.join(OUT, f'{MESH_FIGURE[name]}_mesh.png')
        plt.gcf().savefig(path, dpi=DPI, bbox_inches='tight')
        plt.close('all')
        print(f'  wrote {os.path.relpath(path, ROOT)}  '
              f'({len(mesh["elements"])} elements)', flush=True)

    options = {'F_min': tag.get('f_min', 0.5), 'F_max': tag.get('f_max', 3.0),
               'tolerance': tag.get('tolerance', 0.05)}
    kwargs = {'max_iterations': int(tag['max_iter'])} if 'max_iter' in tag else {}
    with contextlib.redirect_stdout(io.StringIO()):
        with RT._force_fast_kernel(_fem, False):
            result = solve_ssrm(fem_data, debug_level=0,
                                capture_failure_state=True, **options, **kwargs)
    if not result.get('converged', False) or result.get('last_solution') is None:
        print(f'  {name}: SSRM did not converge — no figure written', flush=True)
        return
    with contextlib.redirect_stdout(io.StringIO()):
        plot_fem_results(fem_data, result['last_solution'], fs=result['FS'],
                         failure_solution=result.get('failure_solution'))
    prefix = name.replace('xslope_', '')
    path = os.path.join(OUT, f'{prefix}_results.png')
    plt.gcf().savefig(path, dpi=DPI, bbox_inches='tight')
    plt.close('all')
    print(f'  wrote {os.path.relpath(path, ROOT)}  FS={result["FS"]:.4f}  '
          f'({time.time() - t0:.0f}s)', flush=True)


def main(only=()):
    for tag in figure_tags():
        name = stem(tag)
        if only and not any(s in name for s in only):
            continue
        print(name, flush=True)
        render(tag)


if __name__ == '__main__':
    main(only=tuple(sys.argv[1:]))
