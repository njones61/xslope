"""Migrate every model input file in the repo to the v12 template.

For each model .xlsx (any workbook with a 'main' sheet, excluding the templates
themselves): load slope_data with the version-aware reader, rewrite it into a
fresh copy of the current (v12) template via save_slope_data_to_xlsx, reload
the result, and require every input category to compare equal (the same
equivalence check as run_tests.py's round-trip suite).

Default is a DRY RUN reporting per-file equivalence. --write overwrites the
originals in place (git is the safety net). After a --write run, re-run the
full test suite: the <!-- test: --> FS tags are the outer regression net.

Companion seep artifacts (_mesh.json, _seep.csv) are referenced by basename
convention, so in-place rewriting keeps them attached.

Usage:
    python migrate_to_v12.py            # dry run, report only
    python migrate_to_v12.py --write    # overwrite originals with v12 files
    python migrate_to_v12.py --only PATTERN   # limit to paths containing PATTERN
"""

import argparse
import glob
import os
import shutil
import sys
import tempfile
import warnings

warnings.filterwarnings('ignore')

import run_tests  # reuse the round-trip equivalence machinery
from xslope.fileio import load_slope_data, save_slope_data_to_xlsx, default_template_path

SKIP_BASENAMES = {'input_template.xlsx', 'input_template_v11.xlsx'}


def find_model_files(only=None):
    files = sorted(set(glob.glob('docs/**/*.xlsx', recursive=True)
                       + glob.glob('inputs/**/*.xlsx', recursive=True)))
    out = []
    for f in files:
        base = os.path.basename(f)
        if base.startswith('~$') or base in SKIP_BASENAMES:
            continue
        # archive/ folders hold historical pre-v8 files kept as artifacts; they
        # are not referenced by any test and are deliberately left unmigrated.
        if os.sep + 'archive' + os.sep in f:
            continue
        if only and only not in f:
            continue
        out.append(f)
    return out


def migrate_file(path, write=False):
    """Returns (status, detail): status in PASS/FAIL/SKIP."""
    try:
        d1 = load_slope_data(path)
    except Exception as e:
        return 'SKIP', f'load: {type(e).__name__}: {str(e)[:90]}'

    tmp = tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False).name
    try:
        try:
            save_slope_data_to_xlsx(d1, tmp)
            d2 = load_slope_data(tmp)
        except Exception as e:
            return 'FAIL', f'rewrite: {type(e).__name__}: {str(e)[:90]}'

        mismatches = []
        for k in run_tests.ROUNDTRIP_KEYS:
            mismatches += run_tests._roundtrip_diff(d1.get(k), d2.get(k), k)
        if not d1.get('profile_lines'):
            p1, p2 = d1.get('polygons') or [], d2.get('polygons') or []
            if len(p1) != len(p2):
                mismatches.append(f'polygons: len {len(p1)} vs {len(p2)}')
            else:
                for i, (a, b) in enumerate(zip(p1, p2)):
                    if a.get('mat_id') != b.get('mat_id'):
                        mismatches.append(f'polygons[{i}].mat_id')
                    if not a['polygon'].equals(b['polygon']):
                        mismatches.append(f'polygons[{i}].geom')
        if mismatches:
            return 'FAIL', '; '.join(mismatches[:4])

        if write:
            shutil.move(tmp, path)
            return 'PASS', 'migrated'
        return 'PASS', 'equivalent (dry run)'
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true', help='overwrite originals')
    ap.add_argument('--only', help='limit to paths containing this substring')
    args = ap.parse_args()

    files = find_model_files(args.only)
    counts = {'PASS': 0, 'FAIL': 0, 'SKIP': 0}
    failures = []
    for f in files:
        status, detail = migrate_file(f, write=args.write)
        counts[status] += 1
        flag = '' if status == 'PASS' else f'  <-- {detail}'
        print(f'{status:4}  {f}{flag}')
        if status == 'FAIL':
            failures.append((f, detail))

    print(f"\n{counts['PASS']} pass, {counts['FAIL']} fail, {counts['SKIP']} skip "
          f"of {len(files)} files ({'WROTE' if args.write else 'dry run'})")
    if failures:
        print('\nFailures:')
        for f, d in failures:
            print(f'  {f}: {d}')
        sys.exit(1)


if __name__ == '__main__':
    main()
