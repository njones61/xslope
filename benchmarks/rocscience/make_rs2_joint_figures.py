"""Render the 4-panel FEM/SSRM figures for the RS2 JOINT corpus page
(docs/verification/rs2_joints.md).

The panels, the domain framing, the legend packing, the sidecar export and the
uniform-panel audit are all :mod:`make_rs2_figures`' — this page's figures must
read exactly like the RS2 page's, and a second copy of that code would be a
second thing to keep in step. What is here is the difference: the page the tags
are read from, the rows that are reported without a tag, and the rows the page
documents as not yet built, which get no figure and say so.

Usage, the same as its sibling's:

    python benchmarks/rocscience/make_rs2_joint_figures.py            # every row
    python benchmarks/rocscience/make_rs2_joint_figures.py RJ-18      # one row
    python benchmarks/rocscience/make_rs2_joint_figures.py --from-sidecar
    python benchmarks/rocscience/make_rs2_joint_figures.py --audit
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))

import make_rs2_figures as F                                        # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAGE = os.path.join(ROOT, 'docs', 'verification', 'rs2_joints.md')

#: Rows the page reports without locking, so they carry no fem_ssrm tag and are
#: registered here instead. Settings are the ones the page states for the row.
_JOINT = dict(element_type='tri6', tolerance='0.02', f_min='0.5', f_max='3.0',
              max_iter='250000', tension_srf='false', k0='1')

EXTRA_CASES = [
    # Rows that bracket a factor but cannot lock, because at least one trial of
    # the bracket that defines it -- or of the refinement step that confirms it --
    # reaches the sweep budget without a verdict. They carry no tag, so they are
    # registered here at the settings the page states, and their figures read the
    # mechanism the same way a locked row's does.
    {**_JOINT, 'file': 'files/rocscience/joints/rj003.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-3'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj004.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-4'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj005.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-5'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj006.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-6'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj007.xlsx',
     'target_size': '12.0', 'benchmark': 'RJ-7'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj015.xlsx',
     'target_size': '2.0', 'benchmark': 'RJ-15'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj019.xlsx',
     'target_size': '3.0', 'benchmark': 'RJ-19'},
    # The termination family: problem 1's stepped base and Alejano's release
    # traces. Corpus size is the joint spacing, as everywhere else here — the
    # column width on problem 1, the bedding spacing on 9 to 14.
    {**_JOINT, 'file': 'files/rocscience/joints/rj001a.xlsx',
     'target_size': '10.0', 'benchmark': 'RJ-1a'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj001b.xlsx',
     'target_size': '10.0', 'benchmark': 'RJ-1b'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj001c.xlsx',
     'target_size': '10.0', 'benchmark': 'RJ-1c'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj001d.xlsx',
     'target_size': '10.0', 'benchmark': 'RJ-1d'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj009.xlsx',
     'target_size': '3.0', 'benchmark': 'RJ-9'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj010.xlsx',
     'target_size': '3.0', 'benchmark': 'RJ-10'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj011.xlsx',
     'target_size': '1.5', 'benchmark': 'RJ-11'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj012.xlsx',
     'target_size': '1.5', 'benchmark': 'RJ-12'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj013.xlsx',
     'target_size': '1.5', 'benchmark': 'RJ-13'},
    {**_JOINT, 'file': 'files/rocscience/joints/rj014.xlsx',
     'target_size': '1.5', 'benchmark': 'RJ-14'},
]

#: Rows the page documents as not yet built or blocked, each with the reason. A
#: row here gets no figure and the audit does not count it as missing.
NO_FIGURE = {
    'RJ-16': 'scores a tilt angle found by a gravity sweep, not a strength reduction',
    'RJ-17': "the second material's zone is an element-edge staircase in the vendor "
             'mesh rather than a boundary of its model',
    'RJ-21': 'a two-stage model whose opening is cut in stage 2, and the FEM has no '
             'staged excavation',
    'RJ-20': "the vendor's Voronoi network is 523 digitized traces with no block "
             'size or seed, so the input is not reproducible',
    'RJ-22': "exercises RS2's hyperbolic softening joint law, which the interface "
             'element does not have, and reports no factor of safety',
    'RJ-23': 'reports a stress-displacement curve rather than a factor of safety',
}


def parse_tags(path=PAGE):
    """Every fem_ssrm tag on the joint page."""
    return F.parse_tags(path)


def registered():
    """Every row this producer is responsible for."""
    return parse_tags() + EXTRA_CASES


def audit(out_dir=None, verbose=True):
    """Registered rows with no rendered PNG, and NO_FIGURE entries that are dead."""
    out_dir = out_dir or F.OUT
    missing, dead = [], []
    for tag in registered():
        bench = tag.get('benchmark', '?')
        png = os.path.join(out_dir, f'{bench}.png')
        if not os.path.exists(png):
            missing.append((bench, os.path.basename(png)))
    names = {t.get('benchmark') for t in registered()}
    for bench, why in NO_FIGURE.items():
        if bench in names:
            dead.append((bench, why))
    if verbose:
        print(f'{len(registered())} registered rows '
              f'({len(EXTRA_CASES)} reported-only); figures in '
              f'{os.path.normpath(out_dir)}')
        print(f'  registered, NO figure: {len(missing)}')
        for bench, f in missing:
            print(f'    {bench:22s} {f}')
        print(f'  expected no figure (named): {len(NO_FIGURE)}')
        for bench, why in NO_FIGURE.items():
            print(f'    {bench:22s} {why}')
        print(f'  DEAD exemptions: {len(dead)}')
        for bench, why in dead:
            print(f'    {bench:22s} {why}')
    return missing, dead


if __name__ == '__main__':
    os.makedirs(F.OUT, exist_ok=True)
    args = sys.argv[1:]
    if '--audit' in args:
        missing, dead = audit()
        sys.exit(1 if (missing or dead) else 0)
    from_sidecar = '--from-sidecar' in args
    only = set(a for a in args if not a.startswith('--'))
    cases = registered()
    print(f'{len(cases)} registered rows (fem_ssrm tags + EXTRA_CASES)'
          f"{'  (solve-free re-render from sidecars)' if from_sidecar else ''}")
    for tag in cases:
        bench = tag.get('benchmark', '?')
        if only and bench not in only:
            continue
        t0 = time.time()
        try:
            out, fs = (F.make_figure_from_sidecar(tag) if from_sidecar
                       else F.make_figure(tag))
            exp = tag.get('expected_fs')
            fsx = ('inputs-only' if tag.get('figure') == 'inputs'
                   else f'{fs:.3f}' if fs is not None else 'n/a')
            lock = (f'  lock={float(exp):.3f} d={fs-float(exp):+.3f}'
                    if exp and fs is not None else '')
            print(f'ok   {bench:10s} FS={fsx}{lock}  ({time.time()-t0:.0f}s)  '
                  f'{os.path.basename(out)}', flush=True)
        except Exception as e:
            import traceback
            print(f'FAIL {bench:10s} {type(e).__name__}: {e}', flush=True)
            traceback.print_exc()
