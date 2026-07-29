"""Guard: the corpus builders must reproduce the checked-in corpus files.

DOCTRINE (project): the ``benchmarks/`` builders are authoritative for every
generated input file under ``docs/``.  A corpus ``.xlsx`` is never hand-patched;
it is regenerated.  That contract is only real if a rebuild actually reproduces
what is checked in, so this script rebuilds every declared target into a scratch
directory and compares it to the checked-in file FIELD BY FIELD through
``load_slope_data`` -- never byte-wise, because the xlsx container legitimately
differs run to run (zip mtimes, recalc flags, calcChain removal).

Usage
-----
    python3 benchmarks/verify_rebuild.py            # fast groups (~2 min)
    python3 benchmarks/verify_rebuild.py --all      # every group, incl. rs2
    python3 benchmarks/verify_rebuild.py --group lem --group rs2
    python3 benchmarks/verify_rebuild.py --list
    python3 benchmarks/verify_rebuild.py --json out.json

Exit code is nonzero if any target fails to rebuild or differs in any loaded
field.  Only ``rs2_seep`` is left out of the default set: those builders re-solve
a finite-element seepage problem to write their ``u='seep'`` sidecars, and one of
them (``vp102_transient``) is a 1500-hour transient solve that dominates the whole
run.  Run ``--all`` in a pre-release sweep, and expect it to take an hour.
"""

import argparse
import json
import math
import os
import shutil
import sys
import tempfile
import time
import warnings

warnings.filterwarnings('ignore')

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, 'benchmarks'))
sys.path.insert(0, os.path.join(_ROOT, 'benchmarks', 'rocscience'))


# ---------------------------------------------------------------------------
# Target groups.  Each entry: (module path, attribute holding the output dir,
# corpus dir relative to repo root, callable returning the builder list).
# ---------------------------------------------------------------------------

# The RS2 builders that only WRITE a file...
_RS2_FAST = ('rs2_56a rs2_56b rs2_57a rs2_57b rs2_58a rs2_58b hammah_hb1 '
             'rs2_60a rs2_60b rs2_60c rs2_31d rs2_61a rs2_59 rs2_63 '
             'rs2_66a rs2_66b rs2_66c rs2_66d rs2_66e '
             'rs2_62a rs2_62b rs2_62c rs2_65 rs2_51 '
             'rs2_67a rs2_67c rs2_67d '
             'rs2_68a rs2_68b rs2_68c '
             'rs2_64a rs2_64b rs2_64c rs2_64d rs2_64e rs2_64f '
             'rs2_64g rs2_64h rs2_64i rs2_64j rs2_64k rs2_64l '
             'rs2_64h_split rs2_64l_split rs2_29clay').split()

# ...and the ones that SOLVE a seepage problem to write their u='seep' sidecars.
# vp102_transient alone is a 1500-hour transient solve (tens of minutes); the rest
# are steady solves of a few minutes each.  Split out so the routine guard stays
# usable and the expensive set is an explicit opt-in.
_RS2_SEEP = 'rs2_67b rs2_67e rs2_67f rs2_28a rs2_28b rs2_28c vp102_transient'.split()


def _by_name(names):
    return lambda mod: [getattr(mod, n) for n in names]


def _ssrm_builders(m):
    """One callable per emitted file (the *_all() helpers build a whole station set)."""
    out = [m.build_griffiths1, m.build_griffiths2,
           m.build_griffiths4_r1, m.build_griffiths4_r2,
           m.build_griffiths6_full, m.build_griffiths6_dry]
    out += [(lambda r=r, t=t: m.build_griffiths3(r, t)) for r, t in m.GRIFFITHS3_STATIONS]
    out += [(lambda lh=lh, t=t: m.build_griffiths5(lh, t)) for lh, t in m.GRIFFITHS5_STATIONS]
    for fn in out:                       # __name__ is used to label a build failure
        if not hasattr(fn, '__name__'):
            fn.__name__ = 'griffiths_station'
    return out


GROUPS = {
    'lem': dict(module='benchmarks.build_lem', outattr='OUTDIR',
                corpus='docs/lem/files',
                builders=lambda m: [m.build_acads_simple,
                                    m.build_acads_weak_layer,
                                    m.build_arai_tagyo],
                slow=False),
    'problems': dict(module='build_problems', outattr='OUT',
                     corpus='docs/verification/files/rocscience',
                     builders=lambda m: list(m.BUILDERS),
                     slow=False),
    'gw': dict(module='build_groundwater', outattr='OUT',
               corpus='docs/verification/files/rocscience_gw',
               builders=lambda m: [getattr(m, n) for n in
                                   ('gw001 gw002 gw003 gw004 gw005 gw006a gw006b '
                                    'gw006c gw006e gw007 gw008 gw009a gw009b gw010 '
                                    'gw011 gw012 gw013').split()] + list(m._TRANSIENT),
               slow=False),
    'rs2': dict(module='build_rs2', outattr='OUT',
                corpus='docs/verification/files/rocscience',
                builders=_by_name(_RS2_FAST),
                slow=False),
    'rs2_seep': dict(module='build_rs2', outattr='OUT',
                     corpus='docs/verification/files/rocscience',
                     builders=_by_name(_RS2_SEEP),
                     slow=True),
    # vp038/vp046 do re-solve their seepage sidecars, but on small meshes (~15 s
    # for all four), so they stay in the default set.
    'vp038': dict(module='build_vp038', outattr='OUT',
                  corpus='docs/verification/files/rocscience',
                  builders=lambda m: list(m.BUILDERS),
                  slow=False),
    'vp046': dict(module='build_vp046', outattr='OUT',
                  corpus='docs/verification/files/rocscience',
                  builders=lambda m: list(m.BUILDERS),
                  slow=False),
    'ssrm': dict(module='benchmarks.build_ssrm', outattr='OUTDIR',
                 corpus='docs/fem/files',
                 builders=_ssrm_builders,
                 slow=False),
    'seep': dict(module='benchmarks.build_seep', outattr='OUTDIR',
                 corpus='docs/seep/files',
                 builders=lambda m: [m.build_confined_radial,
                                     lambda: m.build_sheetpile(0.5),
                                     lambda: m.build_sheetpile(0.75)],
                 slow=False),
}

DEFAULT_GROUPS = tuple(k for k, v in GROUPS.items() if not v['slow'])


# ---------------------------------------------------------------------------
# Field-level comparison of two loaded slope_data dicts.
# ---------------------------------------------------------------------------

# Derived geometry: rebuilt from the same primitives on both sides, so a diff
# here is always a shadow of a diff in the primitives that produced it.  Kept
# out of the report to avoid burying the real cause under shapely repr noise.
_DERIVED = {'ground_surface', 'domain_polygon', 'tcrack_surface', 'filename',
            'template_version'}

_TOL = 1e-9


def _norm(v):
    """Canonical, comparable form for any value load_slope_data can produce."""
    import numpy as np
    import pandas as pd
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return 'NaN'
    if isinstance(v, (np.floating, np.integer)):
        return _norm(v.item())
    if isinstance(v, np.ndarray):
        return [_norm(x) for x in v.tolist()]
    if isinstance(v, pd.DataFrame):
        return _norm(v.to_dict('records'))
    if isinstance(v, dict):
        return {str(k): _norm(x) for k, x in sorted(v.items(), key=lambda kv: str(kv[0]))}
    if isinstance(v, (list, tuple)):
        return [_norm(x) for x in v]
    if hasattr(v, 'wkt'):                      # shapely geometry
        return v.wkt
    return v


def _diff(a, b, path=''):
    """Yield (path, a, b) for every leaf mismatch between normalized values."""
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a:
                yield (f'{path}.{k}', '<missing>', b[k])
            elif k not in b:
                yield (f'{path}.{k}', a[k], '<missing>')
            else:
                yield from _diff(a[k], b[k], f'{path}.{k}')
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            yield (f'{path}[len]', len(a), len(b))
            return
        for i, (x, y) in enumerate(zip(a, b)):
            yield from _diff(x, y, f'{path}[{i}]')
        return
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) \
            and not isinstance(a, bool) and not isinstance(b, bool):
        if abs(a - b) > _TOL * max(1.0, abs(a), abs(b)):
            yield (path, a, b)
        return
    if a != b:
        yield (path, a, b)


def compare_files(corpus_path, rebuilt_path):
    """Field-level diff of two input files.  Returns list of (field, corpus, rebuilt)."""
    from xslope.fileio import load_slope_data
    a = load_slope_data(corpus_path)
    b = load_slope_data(rebuilt_path)
    keys = (set(a) | set(b)) - _DERIVED
    out = []
    for k in sorted(keys):
        out.extend(_diff(_norm(a.get(k)), _norm(b.get(k)), k))
    return out


# ---------------------------------------------------------------------------

# Companion files auto-discovered by load_slope_data next to the .xlsx.  A builder
# that does not regenerate them (most don't -- the seepage solve is a separate,
# expensive step) would otherwise read back as "lost the mesh", so the corpus copy
# is staged beside the rebuilt file.  A builder that DOES emit one keeps its own,
# which is then genuinely compared.
_COMPANIONS = ('_mesh.json', '_seep.csv', '_seep2.csv')


def _stage_companions(corpus_dir, outdir, fname):
    stem = os.path.splitext(fname)[0]
    for suffix in _COMPANIONS:
        src = os.path.join(corpus_dir, stem + suffix)
        dst = os.path.join(outdir, stem + suffix)
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copy(src, dst)


def run_group(name, spec, scratch, verbose=True):
    import importlib
    mod = importlib.import_module(spec['module'])
    outdir = os.path.join(scratch, name)
    os.makedirs(outdir, exist_ok=True)
    setattr(mod, spec['outattr'], outdir)
    corpus_dir = os.path.join(_ROOT, spec['corpus'])

    results = []
    for fn in spec['builders'](mod):
        try:
            produced = fn()
        except Exception as e:                                   # builder crash
            results.append({'group': name, 'builder': fn.__name__,
                            'status': 'BUILD_ERROR',
                            'error': f'{type(e).__name__}: {e}'})
            if verbose:
                print(f'  BUILD_ERROR {fn.__name__}: {e}')
            continue
        # A builder may emit a whole set from one solve (vp102_transient writes the
        # five drawdown snapshots off a single 1500 h march) and return the list of
        # names. Check every one of them: str(list) is not a filename, and treating
        # it as one silently reported the whole set as NO_OUTPUT.
        names = (list(produced) if isinstance(produced, (list, tuple))
                 else [produced])
        for produced_one in names:
            rec = {'group': name, 'builder': fn.__name__}
            fname = os.path.basename(str(produced_one))
            rec['file'] = fname
            corpus_path = os.path.join(corpus_dir, fname)
            rebuilt_path = os.path.join(outdir, fname)
            if not os.path.exists(rebuilt_path):
                rec.update(status='NO_OUTPUT')
            elif not os.path.exists(corpus_path):
                rec.update(status='NOT_IN_CORPUS')
            else:
                _stage_companions(corpus_dir, outdir, fname)
                try:
                    diffs = compare_files(corpus_path, rebuilt_path)
                except Exception as e:
                    rec.update(status='LOAD_ERROR', error=f'{type(e).__name__}: {e}')
                else:
                    rec['status'] = 'OK' if not diffs else 'DIFF'
                    if diffs:
                        rec['diffs'] = [{'field': f, 'corpus': c, 'rebuilt': r}
                                        for f, c, r in diffs]
            results.append(rec)
            if verbose and rec['status'] != 'OK':
                print(f"  {rec['status']:12s} {fname}"
                      + (f" ({len(rec.get('diffs', []))} fields)" if rec.get('diffs') else ''))
                for d in rec.get('diffs', [])[:12]:
                    print(f"      {d['field']}: corpus={d['corpus']!r} rebuilt={d['rebuilt']!r}")
    return results


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--group', action='append', default=[],
                    help='group to check (repeatable); default: '
                         + ', '.join(DEFAULT_GROUPS))
    ap.add_argument('--all', action='store_true', help='check every group')
    ap.add_argument('--list', action='store_true', help='list groups and exit')
    ap.add_argument('--json', help='write the full result record to this path')
    ap.add_argument('--keep', action='store_true',
                    help='keep the scratch rebuild directory')
    args = ap.parse_args(argv)

    if args.list:
        for k, v in GROUPS.items():
            print(f"{k:10s} {v['corpus']:40s}"
                  f"{'  (slow: re-solves seepage)' if v['slow'] else ''}")
        return 0

    names = list(GROUPS) if args.all else (args.group or list(DEFAULT_GROUPS))
    for n in names:
        if n not in GROUPS:
            ap.error(f'unknown group {n!r}; choose from {", ".join(GROUPS)}')

    scratch = tempfile.mkdtemp(prefix='xslope_rebuild_')
    t0 = time.time()
    results = []
    try:
        for n in names:
            print(f'[{n}]')
            results.extend(run_group(n, GROUPS[n], scratch))
    finally:
        if args.keep:
            print(f'scratch rebuild kept at {scratch}')
        else:
            shutil.rmtree(scratch, ignore_errors=True)
    elapsed = time.time() - t0

    bad = [r for r in results if r['status'] != 'OK']
    print(f'\n{len(results) - len(bad)}/{len(results)} targets reproduce '
          f'({elapsed:.1f}s)')
    if bad:
        print('FAILED:')
        for r in bad:
            print(f"  {r['status']:12s} {r.get('file', r['builder'])}")
    if args.json:
        with open(args.json, 'w') as f:
            json.dump({'elapsed_s': elapsed, 'groups': names,
                       'results': results}, f, indent=1, default=str)
        print(f'wrote {args.json}')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
