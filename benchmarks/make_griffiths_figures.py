"""Render the mesh and SSRM results figures for the Griffiths & Lane (1999)
corpus page (docs/verification/ssrm.md).

Two figures per case, both through the production plot path at its defaults:

  ``<prefix>_mesh.png``     — plot_fem_data: the mesh coloured by material zone
                              with the displacement boundary conditions drawn.
  ``<prefix>_results.png``  — plot_fem_results at its default three panels
                              (deformation / shear_strain / displace_vector),
                              rendered from the AT-FAILURE field solve_ssrm
                              captures a margin beyond critical, so all three
                              panels are views of one developed mechanism.
                              Titles lead with the locked FS.

The cases are parsed straight out of the ``fem_ssrm`` test tags in ssrm.md rather
than kept in a second list here, so a figure is always rendered on the SAME mesh
and F-bracket the regression lock uses and a retagged mesh size cannot leave a
stale figure behind. Only the rows the page actually figures are rendered; a row
whose section carries no image is named in NO_FIGURE with its reason.

Every solve exports its converged + at-failure field sidecars next to the case
xlsx via export_fem_solution, so a later re-render is solve-free, together with
the MESH those fields were solved on ({stem}_mesh.json) and the run's own record
of how it reached its answer (bracket, tolerance, criterion, per-trial verdicts).
A model's mesh companion is discovered by name, so a case that wrote fields
without one was read back afterwards on whatever mesh happened to be beside it.

Solves run on the pure-NumPy REFERENCE kernel, the path that defines the locks,
so the FS printed in a title is the FS the suite verifies.

Each SSRM solve costs one to four minutes, so a full run is slow. Pass benchmark
prefixes to render a subset.

``--sweeps`` renders the three summary plots instead: the Example 3 and Example 5
factor-of-safety sweeps, which overlay the quad8 locks on a coarse-tri6 curve, and
the Example 1 displacement-vs-F catastrophe sweep. Both FS sweeps take every point
they plot from a tag on the same page, so a moved lock moves the figure.

Run from the repo root:
    python benchmarks/make_griffiths_figures.py                 # all case figures
    python benchmarks/make_griffiths_figures.py griffiths1
    python benchmarks/make_griffiths_figures.py --sweeps        # the summary plots
    python benchmarks/make_griffiths_figures.py --audit         # coverage only
"""

import io
import os
import re
import sys
import time
import warnings
import contextlib

warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import run_tests as RT
from xslope.fileio import load_slope_data
from xslope.fem import (build_fem_data, solve_ssrm, export_fem_solution,
                        ssrm_run_record)
import xslope.fem as _fem
from xslope.mesh import (get_material_polygons, build_mesh_from_polygons,
                         extract_constraint_line_geometry, extract_point_constraints,
                         extract_size_regions, export_mesh_to_json)
from xslope.plot_fem import plot_fem_data, plot_fem_results

ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
SSRM_MD = os.path.join(ROOT, 'docs', 'verification', 'ssrm.md')
OUT = os.path.join(ROOT, 'docs', 'fem', 'images')
DPI = 200

# One mesh figure per section, drawn on the section's first-listed case. The page
# states the geometry is identical across a section's stations ("geometry is
# identical for both cases; only the foundation strength differs"), so a second
# mesh panel would repeat the first picture.
MESH_FIGURE = {
    'xslope_griffiths1': 'griffiths1',
    'xslope_griffiths2': 'griffiths2',
    'xslope_griffiths3_r1': 'griffiths3',
    'xslope_griffiths4_r1': 'griffiths4',
}

# Example 5's mesh panel sits directly under its inputs panel, and that panel is
# drawn at the representative station L/H = 0.5 — the one the page names — so the
# mesh is drawn there too, at the same quad8 settings the Example 5 locks run.
# That station carries no tag of its own, so it is declared here.
MESH_ONLY = [('xslope_griffiths5_0p5', 'griffiths5', 'quad8', 3.5)]

# Rows that get no results figure of their own, each with the reason. Enforced
# both ways by --audit: an entry naming a row that is not tagged, or a row that
# HAS a PNG, is reported as a dead exemption exactly like a missing figure.
NO_FIGURE = {
    'xslope_griffiths5_1': 'the drained end of the sweep IS the Example 1 dry slope, '
                           'figured as griffiths1_results.png',
}

TAG_RE = re.compile(r'<!--\s*test:\s*(.*?)\s*-->')


def figure_tags():
    """Every GATED ``fem_ssrm`` tag on the ssrm page, in page order. The tag is the
    single source of truth for mesh size, element type and F-bracket.

    Gated means the tag names a ``benchmark`` — the rows the suite verifies against
    a published answer, which are exactly the rows the page figures. Selecting on
    the element type instead (quad8) dropped Example 6's full-reservoir station,
    whose lock is tri6: the page shows a results figure for it that no run of this
    script could redraw, and its sidecars aged in place. The untagged tri6 rows are
    sweep points, drawn by the sweep plots rather than figured on their own."""
    return [t for t in RT.parse_test_tags(SSRM_MD)
            if t.get('type') == 'fem_ssrm' and t.get('benchmark')]


def stem(tag):
    return os.path.splitext(os.path.basename(tag['file']))[0]


def render(tag):
    """Solve one case at its tag settings and write its figures + sidecars."""
    name = stem(tag)
    t0 = time.time()
    sd = load_slope_data(tag['file'])
    lines, _n_reinf, _n_pile = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            polys, target_size=tag.get('target_size'),
            element_type=tag['element_type'], lines=lines,
            point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd), **RT._refine_kwargs(tag))
        fem_data = build_fem_data(sd, mesh)

    # The mesh panel needs no solution, so it is written before the solve.
    if name in MESH_FIGURE:
        with contextlib.redirect_stdout(io.StringIO()):
            plot_fem_data(fem_data)
        path = os.path.join(OUT, f'{MESH_FIGURE[name]}_mesh.png')
        plt.gcf().savefig(path, dpi=DPI, bbox_inches='tight')
        plt.close('all')
        print(f'  wrote {os.path.relpath(path, ROOT)}  '
              f'({len(mesh["elements"])} elements)', flush=True)

    if name in NO_FIGURE:
        print(f'  {name}: no results figure ({NO_FIGURE[name]})', flush=True)
        return None

    kwargs = {'max_iterations': int(tag['max_iter'])} if 'max_iter' in tag else {}
    # The search options, in one dict, so the run record persists the same numbers
    # the solve was driven with rather than a second transcription of the tag.
    options = {'F_min': tag.get('f_min', 0.5), 'F_max': tag.get('f_max', 3.0),
               'tolerance': tag.get('tolerance', 0.05)}
    with contextlib.redirect_stdout(io.StringIO()):
        with RT._force_fast_kernel(_fem, False):
            result = solve_ssrm(fem_data, debug_level=0,
                                capture_failure_state=True, **options, **kwargs)
    if not result.get('converged', False):
        print(f'  {name}: SSRM did not converge — no figure written', flush=True)
        return None

    fs = result['FS']
    # The last CONVERGED field is the solution the panels default to; the at-failure
    # (unconverged) snapshot is what they actually render, so all three panels trace
    # the developed mechanism rather than the sub-critical settlement.
    field = result.get('last_solution')
    if field is None:
        print(f'  {name}: SSRM returned no last_solution — no figure written', flush=True)
        return None
    failure = result.get('failure_solution')
    with contextlib.redirect_stdout(io.StringIO()):
        plot_fem_results(fem_data, field, fs=fs, failure_solution=failure)
    prefix = name.replace('xslope_', '')
    path = os.path.join(OUT, f'{prefix}_results.png')
    plt.gcf().savefig(path, dpi=DPI, bbox_inches='tight')
    plt.close('all')

    with contextlib.redirect_stdout(io.StringIO()):
        base = os.path.splitext(tag['file'])[0]
        # What the RUN chose and what its trials found, from the run itself
        # (xslope.fem.ssrm_run_record) — the same record Studio's writer lays its
        # own facts over, so a reader of either file reads one schema. Over it go
        # the facts this script knows and the run does not: what was solved, on
        # what discretization, and which locked row it is.
        meta = ssrm_run_record(result, fem_data, options)
        meta.update({'analysis_type': 'ssrm', 'FS': fs,
                     'element_type': tag['element_type'],
                     'target_size': tag.get('target_size'),
                     'max_iter': int(tag['max_iter']) if 'max_iter' in tag else None,
                     'benchmark': tag.get('benchmark', '')})
        export_fem_solution(fem_data, field, base, meta=meta,
                            failure_solution=failure)
        # The mesh the fields were solved on, beside them. Without it a reload
        # discovers whatever {stem}_mesh.json happens to sit there — for two of
        # these cases, a mesh from an older run with a different node count —
        # and everything downstream reads a section that was never solved.
        export_mesh_to_json(mesh, f'{base}_mesh.json')

    print(f'  wrote {os.path.relpath(path, ROOT)}  FS = {fs:.4f}  '
          f'({len(mesh["elements"])} elements, {time.time() - t0:.0f}s)', flush=True)
    return fs


# ---------------------------------------------------------------- sweep figures
#
# Three summary plots on the page overlay the quad8 locks on a coarse-tri6 curve,
# so a moved lock moves the figure. Every point they plot comes from a fem_ssrm
# tag on the same page, except the three Example 5 stations the page reports but
# does not lock (L/H = 0.2, 0.5, 0.9); those are solved here at the same coarse
# tri6 settings their tagged neighbours use, declared in SWEEP5_UNTAGGED.

SWEEP3_RATIO = {'xslope_griffiths3_r1': 1.0, 'xslope_griffiths3_r0p8': 0.8,
                'xslope_griffiths3_r0p6': 0.6, 'xslope_griffiths3_r0p5': 0.5,
                'xslope_griffiths3_r0p4': 0.4, 'xslope_griffiths3_r0p2': 0.2}
SWEEP5_LH = {'xslope_griffiths5_m0p2': -0.2, 'xslope_griffiths5_0': 0.0,
             'xslope_griffiths5_0p2': 0.2, 'xslope_griffiths5_0p4': 0.4,
             'xslope_griffiths5_0p5': 0.5, 'xslope_griffiths5_0p7': 0.7,
             'xslope_griffiths5_0p9': 0.9, 'xslope_griffiths5_1': 1.0}
SWEEP5_UNTAGGED = [('xslope_griffiths5_0p2', 1.2, 2.0),
                   ('xslope_griffiths5_0p5', 1.0, 1.8),
                   ('xslope_griffiths5_0p9', 0.9, 1.8)]
# Example 3's ratio-1 station is locked on quad8 only; its coarse-tri6 companion is
# reported (the summary table's "1.45 tri6") but carries no tag of its own.
SWEEP3_UNTAGGED = [('xslope_griffiths3_r1', 1.0, 1.8)]
TRI6_COARSE = dict(element_type='tri6', target_size=6, max_iter=4000)


def _tagged(element_type, prefix):
    """{case stem: expected_fs} for every fem_ssrm tag of one element type whose
    file name starts with ``prefix``."""
    out = {}
    for t in RT.parse_test_tags(SSRM_MD):
        if t.get('type') != 'fem_ssrm' or t.get('element_type') != element_type:
            continue
        name = stem(t)
        if name.startswith(prefix):
            out[name] = t['expected_fs']
    return out


def _solve_tri6(name, f_min, f_max, tolerance=0.01):
    """One coarse-tri6 SSRM on the same mesh the tagged stations of this sweep use,
    for a station the page reports but does not lock.

    Solved to ``tolerance=0.01`` rather than the tags' 0.05. At 0.05 the bisection
    stops while the bracket is still 0.05 wide, so the reported midpoint depends on
    where the starting bracket happened to fall — the L/H = 0.9 station reads 1.328
    from [1.0, 1.7] and 1.336 from [0.9, 1.8], the same solve either way. A tagged
    station does not have that problem, because its tag pins the bracket alongside
    the tolerance; an untagged one has no such record, so it is solved tightly enough
    that the answer is a property of the model and not of the bracket."""
    path = os.path.join(ROOT, 'docs', 'fem', 'files', f'{name}.xlsx')
    sd = load_slope_data(path)
    lines, _a, _b = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            polys, target_size=TRI6_COARSE['target_size'],
            element_type=TRI6_COARSE['element_type'], lines=lines,
            point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd))
        fem_data = build_fem_data(sd, mesh)
        with RT._force_fast_kernel(_fem, False):
            res = solve_ssrm(fem_data, F_min=f_min, F_max=f_max, tolerance=tolerance,
                             debug_level=0,
                             max_iterations=TRI6_COARSE['max_iter'])
    if not res.get('converged'):
        raise RuntimeError(f'{name}: coarse tri6 SSRM did not converge')
    return res['FS']


def render_mesh_only():
    """Mesh panels for the stations declared in MESH_ONLY — the ones the page draws
    at a station that carries no tag of its own. No solve."""
    for name, prefix, element_type, target_size in MESH_ONLY:
        path_xlsx = os.path.join(ROOT, 'docs', 'fem', 'files', f'{name}.xlsx')
        sd = load_slope_data(path_xlsx)
        lines, _a, _b = extract_constraint_line_geometry(sd)
        polys = get_material_polygons(sd, reinf_lines=lines)
        with contextlib.redirect_stdout(io.StringIO()):
            mesh = build_mesh_from_polygons(
                polys, target_size=target_size, element_type=element_type,
                lines=lines, point_constraints=extract_point_constraints(sd),
                size_regions=extract_size_regions(sd))
            fem_data = build_fem_data(sd, mesh)
            plot_fem_data(fem_data)
        path = os.path.join(OUT, f'{prefix}_mesh.png')
        plt.gcf().savefig(path, dpi=DPI, bbox_inches='tight')
        plt.close('all')
        print(f'  wrote {os.path.relpath(path, ROOT)}  '
              f'({len(mesh["elements"])} elements, {name})', flush=True)


def _sweep_axes(title):
    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    return fig, ax


def sweep_griffiths3():
    tri = _tagged('tri6', 'xslope_griffiths3')
    quad = _tagged('quad8', 'xslope_griffiths3')
    for name, f_min, f_max in SWEEP3_UNTAGGED:
        tri[name] = _solve_tri6(name, f_min, f_max)
        print(f'  {name} (reported, not locked): tri6 FS = {tri[name]:.3f}', flush=True)
    pts = sorted((SWEEP3_RATIO[k], v) for k, v in tri.items() if k in SWEEP3_RATIO)
    qpts = sorted((SWEEP3_RATIO[k], v) for k, v in quad.items() if k in SWEEP3_RATIO)
    fig, ax = _sweep_axes('Griffiths & Lane (1999) Example 3 — thin weak layer\n'
                          r'$c_{u1}/\gamma H = 0.25$')
    ax.plot([p[0] for p in pts], [p[1] for p in pts], '-s', color='#1f5c8b',
            mfc='white', ms=8, lw=2, label='xslope SSRM (coarse tri6)')
    ax.plot([p[0] for p in qpts], [p[1] for p in qpts], 'o', color='#b02020',
            ms=11, label='xslope SSRM (gated quad8)')
    ax.plot([1.0], [1.47], '*', color='black', ms=20,
            label='Taylor (1937): FOS = 1.47')
    ax.axvline(0.6, color='0.6', ls=':', lw=1.5)
    ax.text(0.62, 1.87, 'transition\n' r'$c_{u2}/c_{u1} \approx 0.6$', color='0.45')
    ax.text(0.62, 1.55, 'base-circle plateau', color='#1f5c8b')
    ax.text(0.24, 0.80, 'weak-layer\nmechanism\n(falls ~linearly)', color='#1f5c8b')
    ax.set_xlabel(r'$c_{u2}/c_{u1}$')
    ax.set_ylabel('Factor of safety, FOS')
    ax.set_xlim(0.0, 1.12)
    ax.set_ylim(0.0, 2.0)
    ax.legend(loc='upper left')
    fig.tight_layout()
    path = os.path.join(OUT, 'griffiths3_sweep.png')
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print(f'  wrote {os.path.relpath(path, ROOT)}', flush=True)


def sweep_griffiths5():
    tri = _tagged('tri6', 'xslope_griffiths5')
    quad = _tagged('quad8', 'xslope_griffiths5')
    for name, f_min, f_max in SWEEP5_UNTAGGED:
        tri[name] = _solve_tri6(name, f_min, f_max)
        print(f'  {name} (reported, not locked): tri6 FS = {tri[name]:.3f}', flush=True)
    pts = sorted((SWEEP5_LH[k], v) for k, v in tri.items() if k in SWEEP5_LH)
    qpts = sorted((SWEEP5_LH[k], v) for k, v in quad.items() if k in SWEEP5_LH)
    fig, ax = _sweep_axes("Griffiths & Lane (1999) Example 5 — 'slow' drawdown sweep\n"
                          'XSLOPE SSRM vs Fig. 15')
    ax.plot([p[0] for p in pts], [p[1] for p in pts], '-s', color='#1f5c8b',
            mfc='white', ms=8, lw=2, label='XSLOPE SSRM (coarse tri6)')
    ax.plot([p[0] for p in qpts], [p[1] for p in qpts], 'o', color='#b02020',
            ms=11, label='XSLOPE SSRM gated locks (quad8)')
    ax.plot([0.0, 1.0], [1.85, 1.40], '*', color='#1f6b3a', ms=20,
            label='Published LE chart anchors')
    ax.annotate('Morgenstern (1963)\nF = 1.85', xy=(0.0, 1.86), xytext=(-0.18, 1.93),
                color='#1f6b3a',
                arrowprops=dict(arrowstyle='-', color='#1f6b3a'))
    ax.annotate('Bishop & Morgenstern (1960)\nFOS = 1.4', xy=(1.0, 1.41),
                xytext=(0.52, 1.58), color='#1f6b3a',
                arrowprops=dict(arrowstyle='-', color='#1f6b3a'))
    lo = min(pts, key=lambda p: p[1])
    ax.annotate(f'minimum FOS = {lo[1]:.2f}\nat L/H = {lo[0]:g}\n(paper: ~1.3 at 0.7)',
                xy=lo, xytext=(0.16, 1.24), color='#1f5c8b',
                arrowprops=dict(arrowstyle='-', color='#1f5c8b'))
    ax.set_xlabel('Drawdown ratio  L/H')
    ax.set_ylabel('Factor of safety, FOS')
    ax.set_xlim(-0.28, 1.12)
    ax.set_ylim(1.2, 2.0)
    ax.legend(loc='upper right')
    fig.tight_layout()
    path = os.path.join(OUT, 'griffiths5_sweep.png')
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print(f'  wrote {os.path.relpath(path, ROOT)}', flush=True)


# Example 1's displacement-vs-F catastrophe sweep: one solve_fem per trial F on the
# SSRM-1 mesh, plotting the maximum displacement. The F list spans the flat branch
# and the upturn Griffiths & Lane's own Fig. 2 reports.
SWEEP1_F = [1.0, 1.2, 1.3, 1.35, 1.4, 1.45, 1.5, 1.6]


def sweep_griffiths1():
    from xslope.fem import solve_fem
    import numpy as np
    tag = next(t for t in figure_tags() if stem(t) == 'xslope_griffiths1')
    sd = load_slope_data(tag['file'])
    lines, _a, _b = extract_constraint_line_geometry(sd)
    polys = get_material_polygons(sd, reinf_lines=lines)
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = build_mesh_from_polygons(
            polys, target_size=tag.get('target_size'), element_type=tag['element_type'],
            lines=lines, point_constraints=extract_point_constraints(sd),
            size_regions=extract_size_regions(sd))
        fem_data = build_fem_data(sd, mesh)
    fs = None
    with contextlib.redirect_stdout(io.StringIO()):
        with RT._force_fast_kernel(_fem, False):
            res = solve_ssrm(fem_data, F_min=tag['f_min'], F_max=tag['f_max'],
                             tolerance=tag['tolerance'], debug_level=0,
                             max_iterations=int(tag['max_iter']))
    fs = res['FS']
    dmax = []
    for f in SWEEP1_F:
        with contextlib.redirect_stdout(io.StringIO()):
            sol = solve_fem(fem_data, F=f, max_iterations=int(tag['max_iter']),
                            debug_level=0, fast_kernel=False)
        disp = np.asarray(sol['displacements']).reshape(-1, 2)
        d = float(np.max(np.hypot(disp[:, 0], disp[:, 1])))
        dmax.append(d)
        print(f'  F = {f:.2f}: max|u| = {d:.4f}', flush=True)
    fig, ax = _sweep_axes('Griffiths & Lane Example 1 — displacement vs strength '
                          'reduction factor')
    ax.plot(SWEEP1_F, dmax, '-o', color='#1f77b4', ms=9, lw=2.5)
    ax.axvline(fs, color='#d62728', lw=2)
    ax.axvline(1.40, color='#d62728', lw=2, ls='--')
    span = max(dmax) - min(dmax)
    ax.annotate(f'FS = {fs:.2f}\n(equilibrium criterion)',
                xy=(fs, min(dmax) + 0.78 * span), xytext=(1.06, min(dmax) + 0.82 * span),
                color='#d62728', arrowprops=dict(arrowstyle='->', color='#d62728'))
    ax.annotate('displacement upturn ≈ 1.40\n(G&L report 1.4)',
                xy=(1.40, min(dmax) + 0.46 * span),
                xytext=(1.07, min(dmax) + 0.50 * span), color='#d62728',
                arrowprops=dict(arrowstyle='->', color='#d62728'))
    ax.set_xlabel('Strength reduction factor  F')
    ax.set_ylabel('Maximum displacement (ft)')
    fig.tight_layout()
    path = os.path.join(OUT, 'griffiths1_sweep.png')
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    print(f'  wrote {os.path.relpath(path, ROOT)}', flush=True)
    print('  sweep readings: ' +
          ', '.join(f'{f:.2f}->{d:.4f}' for f, d in zip(SWEEP1_F, dmax)), flush=True)


SWEEPS = {'griffiths1': sweep_griffiths1, 'griffiths3': sweep_griffiths3,
          'griffiths5': sweep_griffiths5}


def audit(tags):
    """Report every tagged row with no rendered results PNG, and every NO_FIGURE
    entry that no longer describes a real figure-less row. Non-zero on either."""
    bad = []
    for t in tags:
        name = stem(t)
        png = os.path.join(OUT, f'{name.replace("xslope_", "")}_results.png')
        has = os.path.exists(png)
        if name in NO_FIGURE and has:
            bad.append(f'DEAD exemption: {name} has {os.path.basename(png)}')
        elif name not in NO_FIGURE and not has:
            bad.append(f'MISSING figure: {name} -> {os.path.basename(png)}')
    for name in NO_FIGURE:
        if name not in {stem(t) for t in tags}:
            bad.append(f'DEAD exemption: {name} carries no quad8 fem_ssrm tag')
    for line in bad:
        print(line)
    print(f'{len(tags)} tagged rows, {len(bad)} unexplained')
    return 1 if bad else 0


def main(argv):
    tags = figure_tags()
    if '--audit' in argv:
        return audit(tags)
    wanted = [a for a in argv if not a.startswith('-')]
    if '--mesh-only' in argv:
        render_mesh_only()
        return 0
    if '--sweeps' in argv:
        for name, fn in SWEEPS.items():
            if wanted and name not in wanted:
                continue
            print(f'{name} sweep', flush=True)
            fn()
        return 0
    if wanted:
        tags = [t for t in tags if any(w in stem(t) for w in wanted)]
    print(f'{len(tags)} case(s)')
    for t in tags:
        print(stem(t), flush=True)
        render(t)
    if not wanted or any('griffiths5' in w for w in wanted):
        render_mesh_only()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
