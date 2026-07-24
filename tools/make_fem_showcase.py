"""Render the FEM at-failure deformation showcase for the docs landing page
(docs/index.md, Finite Element Slope Stability section).

The figure is a single wide "Variant E" deformation panel: the at-failure
(unconverged) mechanism drawn as a crisp BLACK deformed mesh over the dashed
light outline of the undeformed extent (Griffiths & Lane convention), with an
honest stated scale and a slim "… FS = X" title. It shows the paper-famous
rotational slump of Griffiths & Lane (1999) Example 1 — the whole slope has
slumped: the crest drops below the dashed original crest line, the face bulges
out and down, the toe heaves — the failure mechanism the FEM/SSRM engine lets
emerge on its own, with no assumed surface.

Solve-free by construction. The at-failure field is read from the committed
``*_fem_failure_*`` sidecars next to the case xlsx; this script only rebuilds
the mesh (matching the recorded recipe) and imports them, so it never touches
the solver. The mesh recipe AND the reported FS are read straight from the
SSRM-1 benchmark tag in docs/verification/ssrm.md — the single source of truth
for the lock — so a retagged mesh or a re-recorded FS can never leave a stale
showcase behind.

Run from the repo root:
    python tools/make_fem_showcase.py
"""

import io
import os
import re
import sys
import contextlib
import warnings

warnings.filterwarnings('ignore')

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from xslope.fileio import load_slope_data
from xslope.fem import build_fem_data, import_fem_solution
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints)
from xslope.plot_fem import plot_fem_results

SSRM_MD = os.path.join(ROOT, 'docs', 'verification', 'ssrm.md')
OUT = os.path.join(ROOT, 'docs', 'fem', 'images', 'fem_failure_showcase.png')

# The showcase case: Griffiths & Lane (1999) Example 1, the classic homogeneous
# rotational slump, locked as benchmark SSRM-1.
BENCHMARK = 'SSRM-1'

TAG_RE = re.compile(r'<!--\s*test:\s*(.*?)\s*-->')


def find_tag(benchmark, path=SSRM_MD):
    """The fem_ssrm test tag carrying ``benchmark=<benchmark>`` on ssrm.md, as a
    key->value dict. The tag is the single source of truth for the mesh recipe
    (element_type, target_size) and the reported factor of safety (expected_fs),
    so the showcase renders on the SAME mesh and reports the SAME FS as the lock."""
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
            if kv.get('type') == 'fem_ssrm' and kv.get('benchmark') == benchmark:
                return kv
    raise RuntimeError(f'no fem_ssrm tag with benchmark={benchmark} on {path}')


def build_fem(xlsx, element_type, target_size):
    """Mesh + fem_data exactly as run_tests.run_fem_test builds them, so the node
    and element counts match the recorded sidecars (import_fem_solution validates
    the match and raises on any drift)."""
    sd = load_slope_data(xlsx)
    lines, _nr, _np = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    mesh = build_mesh_from_polygons(
        polys, target_size=target_size, element_type=element_type,
        lines=lines, point_constraints=extract_point_constraints(sd))
    with contextlib.redirect_stdout(io.StringIO()):
        fem_data = build_fem_data(sd, mesh)
    return sd, fem_data


def main():
    tag = find_tag(BENCHMARK)
    # tag file path is relative to docs/verification/
    xlsx = os.path.normpath(os.path.join(ROOT, 'docs', 'verification', tag['file']))
    element_type = tag.get('element_type', 'quad8')
    target_size = float(tag['target_size'])
    fs = float(tag['expected_fs'])
    stem = os.path.splitext(xlsx)[0]

    _sd, fem_data = build_fem(xlsx, element_type, target_size)
    with contextlib.redirect_stdout(io.StringIO()):
        solution = import_fem_solution(fem_data, stem)
    failure = solution.get('failure_solution')
    if failure is None:
        raise RuntimeError(
            f'no at-failure sidecar next to {os.path.basename(xlsx)} — expected '
            f'{os.path.basename(stem)}_fem_failure_*.csv')

    # Single wide "Variant E" deformation panel: at-failure field, BLACK deformed
    # mesh, dashed original outline, honest stated scale, slim "… FS = X" title.
    # fig=None lets plot_fem_results size the single-panel figure to the slope's own
    # (wide, short) data aspect, matching the wide landing-page neighbours; the outer
    # figsize width drives the render resolution together with dpi below.
    fig, _ax = plot_fem_results(
        fem_data, solution, plot_type=['deformation'],
        field_state='failure', failure_solution=failure, fs=fs,
        show_original='outline', deformed_color='k',
        figsize=(12, 7), show_legend=True)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=160, bbox_inches='tight')
    plt.close(fig)

    from PIL import Image
    size = Image.open(OUT).size
    print(f'wrote {OUT}  {size[0]}x{size[1]}px  '
          f'({os.path.basename(xlsx)}, {element_type} @ {target_size}, FS={fs})')


if __name__ == '__main__':
    main()
