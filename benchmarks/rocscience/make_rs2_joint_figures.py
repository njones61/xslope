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
    python benchmarks/rocscience/make_rs2_joint_figures.py --capture-only RJ-18
        # redraw the at-failure picture from the row's stored search record
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))

import make_rs2_figures as F                                        # noqa: E402

# Every row on this page is jointed: its deformed-section panel draws each block
# under its own tint, so the bodies can be told apart.
F.COLOR_BLOCKS = True

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAGE = os.path.join(ROOT, 'docs', 'verification', 'rs2_joints.md')

#: The rows registered without a tag -- ``EXTRA_CASES`` (reported-only brackets)
#: and ``SWEEP_CASES`` (tilt sweeps) -- live in :mod:`rs2_joint_cases`, the one
#: table this producer and the corpus builder both read: the builder writes each
#: row's K0 into its workbook from the same entry this producer solves it with.
from rs2_joint_cases import EXTRA_CASES, SWEEP_CASES                 # noqa: E402


def parse_sweep_tags(path=PAGE):
    """Every ``fem_tilt`` tag on the joint page, as string key->value dicts."""
    cases = []
    with open(path) as fh:
        for line in fh:
            m = F.TAG_RE.search(line)
            if not m:
                continue
            kv = {}
            for part in m.group(1).split(','):
                if '=' in part:
                    k, v = part.split('=', 1)
                    kv[k.strip()] = v.strip()
            if kv.get('type') == 'fem_tilt':
                cases.append(kv)
    return cases


def sweep_cases():
    """Every sweep row, each ONCE: the page's fem_tilt tags, then any
    SWEEP_CASES entry whose row carries no tag yet."""
    tagged = parse_sweep_tags()
    have = {t.get('benchmark') for t in tagged}
    return tagged + [t for t in SWEEP_CASES if t.get('benchmark') not in have]


#: The iteration ceiling a jointed row's trials extend into, as a multiple of the
#: per-trial limit. solve_fem extends a trial still slowing at its limit up to
#: max_iterations_ceiling and takes max(ceiling, limit), so a row at 250 000 with
#: no ceiling of its own has a ceiling EQUAL to its limit and cannot extend at
#: all. Four times the limit lets a slow equilibrium finish; a trial holding
#: steady still stops at its limit, since the extension is conditional on the
#: movement dying away. A tag's own max_iter_ceiling wins.
CEILING_FACTOR = 4


def _with_ceiling(tag):
    """``tag`` with its ``max_iter_ceiling`` stated: the tag's own, else
    ``CEILING_FACTOR`` times its ``max_iter``."""
    max_iterations = int(float(tag.get('max_iter', 250000)))
    return {**tag, 'max_iter_ceiling': str(int(float(
        tag.get('max_iter_ceiling', CEILING_FACTOR * max_iterations))))}


def sweep_registered():
    """Every sweep row, by benchmark."""
    return {t['benchmark']: t for t in sweep_cases()}


def _tilt(k):
    """The tilt a coefficient stands for (run_tests.tilt_of, the one rule)."""
    sys.path.insert(0, F.ROOT)
    import run_tests as RT
    return RT.tilt_of(k)


def _solve_at(sd, mesh, k, tag):
    """One solve at FULL STRENGTH with the model pushed by ``k`` toward the face.

    F = 1 and no reduction: what the sweep asks of each coefficient is whether the
    model stands under it, which is the question the manual's tilt table asks and
    not the one a strength reduction answers. The solve itself is
    ``run_tests.solve_fem_tilt`` — the same call the suite's ``fem_tilt`` check
    makes — so the figure and the lock are one computation, reference kernel
    forced in both.
    """
    sys.path.insert(0, F.ROOT)
    import run_tests as RT
    return RT.solve_fem_tilt(sd, mesh, k, tag)


def make_sweep_figure(tag, dpi=150):
    """Render a sweep row's composite: the state at its first toppling coefficient.

    Two solves, each seconds: the last coefficient the model stands at and the
    first it goes at. The composite's right-hand panels draw the toppling state,
    and their titles carry the tilt both coefficients stand for instead of a
    factor of safety, because a factor of safety is not what this row measures.
    """
    bench = tag['benchmark']
    sd, _fd, _path, mesh = F._build(tag)
    _sd_s, _fd_s, standing = _solve_at(sd, mesh, tag['k_stand'], tag)
    sd_f, fem_f, toppling = _solve_at(sd, mesh, tag['k_fail'], tag)
    note = (f"at a tilt of {_tilt(tag['k_fail']):.2f}\u00b0, "
            f"the first it topples at; it stands at {_tilt(tag['k_stand']):.2f}\u00b0")
    # The inputs panel annotates the coefficient itself (plot.plot_inputs draws it
    # for mode='fem', arrow and all), so the panel says which way the model is
    # pushed and the titles say what tilt that is.
    out = F.render_figure(bench, sd_f, fem_f, standing, failure=toppling,
                          fs=None, out_dir=F.OUT, dpi=dpi, standing_note=note)
    return out, None


#: Rows the page documents as not yet built or blocked, each with the reason. A
#: row here gets no figure and the audit does not count it as missing.
NO_FIGURE = {
    'RJ-21': 'a two-stage model whose opening is cut in stage 2, and the FEM has no '
             'staged excavation',
    'RJ-22': "exercises RS2's hyperbolic softening joint law, which the interface "
             'element does not have, and reports no factor of safety',
    'RJ-23': 'reports a stress-displacement curve rather than a factor of safety',
}


def parse_tags(path=PAGE):
    """Every fem_ssrm tag on the joint page."""
    return F.parse_tags(path)


def registered():
    """Every row this producer is responsible for, each row ONCE.

    `EXTRA_CASES` exists for rows the page reports without locking, which carry no
    tag for `parse_tags` to find. A row that later LOCKS gains a tag and keeps its
    entry here until someone remembers to delete it, and the producer then solves
    it twice — the second solve overwriting the first's sidecars with the identical
    answer, at the cost of the whole bracket. The tag wins wherever both exist, so
    the registration stays right whether or not the list has been tidied.
    """
    tagged = parse_tags()
    have = {t.get('benchmark') for t in tagged}
    return [_with_ceiling(t) for t in
            tagged + [t for t in EXTRA_CASES if t.get('benchmark') not in have]]


def superseded():
    """EXTRA_CASES entries whose row now carries a tag — dead registrations."""
    have = {t.get('benchmark') for t in parse_tags()}
    return [t.get('benchmark') for t in EXTRA_CASES if t.get('benchmark') in have]


def audit(out_dir=None, verbose=True):
    """Registered rows with no rendered PNG, and NO_FIGURE entries that are dead."""
    out_dir = out_dir or F.OUT
    missing, dead = [], []
    for tag in list(registered()) + list(sweep_cases()):
        bench = tag.get('benchmark', '?')
        png = os.path.join(out_dir, f'{bench}.png')
        if not os.path.exists(png):
            missing.append((bench, os.path.basename(png)))
    names = {t.get('benchmark') for t in registered()} | set(sweep_registered())
    for bench, why in NO_FIGURE.items():
        if bench in names:
            dead.append((bench, why))
    stale = superseded()
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
        print(f'  EXTRA_CASES entries superseded by a tag: {len(stale)}')
        for bench in stale:
            print(f'    {bench:22s} the row locked; delete its EXTRA_CASES entry')
    return missing, dead + [(b, 'superseded by a tag') for b in stale]


if __name__ == '__main__':
    os.makedirs(F.OUT, exist_ok=True)
    args = sys.argv[1:]
    if '--audit' in args:
        missing, dead = audit()
        bad = F.audit_captures(registered())
        sys.exit(1 if (missing or dead or bad) else 0)
    from_sidecar = '--from-sidecar' in args
    capture_only = '--capture-only' in args
    if from_sidecar and capture_only:
        sys.exit('--from-sidecar and --capture-only are two different redraws; '
                 'pass one.')
    only = set(a for a in args if not a.startswith('--'))
    sweeps = sweep_cases()
    cases = list(registered()) + list(sweeps)
    print(f'{len(cases)} registered rows (fem_ssrm tags + EXTRA_CASES + '
          f'{len(sweeps)} sweep row(s))'
          f"{'  (solve-free re-render from sidecars)' if from_sidecar else ''}")
    for tag in cases:
        bench = tag.get('benchmark', '?')
        if only and bench not in only:
            continue
        t0 = time.time()
        try:
            if bench in sweep_registered():
                if from_sidecar or capture_only:
                    print(f'skip {bench:10s} a sweep row is re-rendered by '
                          f're-solving (seconds); it has no search record',
                          flush=True)
                    continue
                out, fs = make_sweep_figure(tag)
            else:
                out, fs = (F.make_figure_from_sidecar(tag) if from_sidecar
                           else F.make_figure_capture_only(tag) if capture_only
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
