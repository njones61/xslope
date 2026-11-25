# Copyright 2025 Norman L. Jones
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""A git worktree tests like its main checkout.

Two things the suite reads belong to the machine's main checkout, not to the
branch a worktree holds: the row-timing record (``test/row_timings.json``, which
decides what the standard run holds back for the gate) and the private test repo
(a sibling folder of the main checkout). Worktrees live in a folder of their own,
so a runner that looked beside itself found neither. With no timing record a
fresh worktree's standard run held nothing back and ran the gate's eight-hour
rows; with no private repo every run made in a worktree, the 1.1.0 release
checks among them, silently skipped the private rows.

Built on throwaway folders laid out the way git lays them out:

  A. a plain clone (``.git`` is a folder) is its own main checkout;
  B. a worktree (``.git`` is a ``gitdir:`` file, its git folder has a
     ``commondir``) resolves to the main checkout, by an absolute and by a
     relative ``gitdir``, and without a ``commondir`` file;
  C. a ``.git`` file that is not a worktree pointer, or one that points at a
     folder with no ``test/``, leaves the clone as its own main checkout;
  D. the private repo is found beside the main checkout from a worktree, and
     not beside the worktree;
  E. this repo's own timing record path is in its main checkout.

Runs on its own (``python test/runner_home_check.py``) or as the suite's
``runner_home`` row.
"""

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _runner():
    spec = importlib.util.spec_from_file_location('_run_tests_home', REPO / 'run_tests.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _layout(base):
    """main/ (a clone with test/), main/.git/worktrees/wt, and wt/ beside it."""
    main = base / 'clones' / 'xslope'
    (main / 'test').mkdir(parents=True)
    gitdir = main / '.git' / 'worktrees' / 'wt'
    gitdir.mkdir(parents=True)
    (gitdir / 'commondir').write_text('../..\n')
    wt = base / 'xslope_worktrees' / 'wt'
    wt.mkdir(parents=True)
    (wt / '.git').write_text(f'gitdir: {gitdir}\n')
    return main, gitdir, wt


def run():
    rt = _runner()
    fails = []

    def same(got, want, label):
        if Path(got).resolve() != Path(want).resolve():
            fails.append(f"{label}: got {got}, expected {want}")

    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        main, gitdir, wt = _layout(base)

        # A
        same(rt._main_checkout(main), main, "A plain clone")
        # B
        same(rt._main_checkout(wt), main, "B worktree, absolute gitdir")
        (wt / '.git').write_text(f'gitdir: {os.path.relpath(gitdir, wt)}\n')
        same(rt._main_checkout(wt), main, "B worktree, relative gitdir")
        (gitdir / 'commondir').unlink()
        same(rt._main_checkout(wt), main, "B worktree, no commondir")
        (gitdir / 'commondir').write_text('../..\n')
        (wt / '.git').write_text(f'gitdir: {gitdir}\n')
        # C
        odd = base / 'odd'
        odd.mkdir()
        (odd / '.git').write_text('not a pointer\n')
        same(rt._main_checkout(odd), odd, "C .git file that is not a pointer")
        bare = base / 'bare'
        (bare / '.git' / 'worktrees' / 'x').mkdir(parents=True)
        lone = base / 'lone'
        lone.mkdir()
        (lone / '.git').write_text(f"gitdir: {bare / '.git' / 'worktrees' / 'x'}\n")
        same(rt._main_checkout(lone), lone, "C pointer to a checkout with no test/")

        # D
        saved = os.environ.pop('XSLOPE_PRIVATE_TESTS', None)
        try:
            (main.parent / 'xslope_private').mkdir()
            got = rt._private_dir(wt)
            if got is None:
                fails.append("D worktree: private repo beside the main checkout not found")
            else:
                same(got, main.parent / 'xslope_private', "D worktree private repo")
            (wt.parent / 'xslope_private').mkdir()
            (main.parent / 'xslope_private').rmdir()
            if rt._private_dir(wt) is not None:
                fails.append("D worktree: took a private repo beside the worktree, "
                             "not beside its main checkout")
        finally:
            if saved is not None:
                os.environ['XSLOPE_PRIVATE_TESTS'] = saved

    # E
    same(Path(rt.ROW_TIMINGS_PATH).parent.parent, rt._main_checkout(),
         "E this repo's timing record")
    return fails


if __name__ == '__main__':
    failures = run()
    for f in failures:
        print("FAIL", f)
    print("runner_home:", "FAIL" if failures else "PASS")
    sys.exit(1 if failures else 0)
