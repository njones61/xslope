"""The printed-precision rule bites: a tagged number on any page under docs/ is
checked at the precision the tag writes it, not only inside its tolerance, and a
tutorial whose prose no longer agrees with its tag fails the restatement row.

1.0.8's circle slicing moved factors of safety in the third decimal (0.918 ->
0.9186, 1.301 -> 1.305). Every tag guarding them carried a tolerance of 0.005 or
0.01, so every tag kept passing, and the pages printed the old numbers: the
tutorials for weeks (fixed in 1.1.0, when the rule covered them alone), and 162
values on 93 sample and verification tags until 1.1.1. run_tests.py requires a
row's run to ROUND to the tagged number at the tag's own precision
(``_as_written_miss``) on every page under docs/. This file proves that rule
fails what it must and passes what it must, on synthetic tags, solving nothing:

  1. THE RULE. A tag written 0.918 fails a run of 0.9186 (it rounds to 0.919)
     and passes 0.9184. A per-method tag keeps each method's written precision,
     trailing zeros included: fs_bishop=1.210 is three decimals, so 1.2106 fails
     although the float 1.21 reads as two. A number written in scientific
     notation is read to its significant figures. One element of a list-valued
     key (an instant of a curve, a head at a point) is read from its own text.

  2. THE SCOPE. Every page under docs/ is read this way: the same tag on a
     tutorial, a verification page and a sample page fails 0.9186. A page
     outside docs/ (the private course rows, the private repo's tests/ folder)
     keeps the tolerance alone and passes it.

  3. THE OPT-OUT. ``precision=tolerance`` with a ``precision_reason`` takes a tag
     back to the tolerance alone. Without the reason the row is refused before
     anything is solved; so is a misspelt value, and a reason with no opt-out.
     Every tag on every real page under docs/ must carry a well-formed opt-out or
     none, including rows the standard tier holds back and so never reaches.

  4. THE FAST KERNEL. A strength-reduction row's fast-kernel pass must reproduce
     the page's number at its written precision too: 1.730 is three decimals,
     so a kernel reading of 1.7334, inside the old two-decimal window of the
     float 1.73, falls through to the reference solve.

  5. THE PROSE. Once a tag is current, a page still printing the old number
     disagrees with it, and the restatement row fails on that: a method table
     printing Bishop 1.251 under fs_bishop=1.250 fails the row and names the
     line; the page printing 1.250, or 1.25, passes it.

Synthetic pages are written to a temporary folder; the docs are never written.
"""
import os
import sys
import tempfile
from pathlib import Path

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

FAILURES = []


def check(label, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'}  {label}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILURES.append(label + (f" ({detail})" if detail else ""))


TAGS = """
<!-- test: file=files/a.xlsx, type=circular_search, method=spencer, expected_fs=0.918, tolerance=0.005 -->
<!-- test: file=files/a.xlsx, type=circular_search, fs_spencer=0.918, fs_bishop=1.210, tolerance=0.005 -->
<!-- test: file=files/a.xlsx, type=seep, expected_flowrate=4.916e-07, tolerance=0.005 -->
<!-- test: file=files/a.xlsx, type=circular_search, method=spencer, expected_fs=0.918, tolerance=0.005, precision=tolerance -->
<!-- test: file=files/a.xlsx, type=circular_search, method=spencer, expected_fs=0.918, tolerance=0.005, precision=tolerance, precision_reason=the count differs between computers -->
<!-- test: file=files/a.xlsx, type=circular_search, method=spencer, expected_fs=0.918, tolerance=0.005, precision=tol, precision_reason=misspelt -->
<!-- test: file=files/a.xlsx, type=circular_search, method=spencer, expected_fs=0.918, tolerance=0.005, precision_reason=no opt-out -->
<!-- test: file=files/a.xlsx, type=fem_ssrm, expected_fs=1.730, tolerance=0.01 -->
"""


def _pages(RT, root):
    """The synthetic tags parsed from a tutorial page, a verification page, a
    sample page two folders down, and a private course page outside docs/, all
    under ``root``, with the rule's docs folder pointed at ``root``'s."""
    docs = Path(root) / 'docs'
    folders = (docs / 'tutorials', docs / 'verification', docs / 'lem' / 'files',
               Path(root) / 'private' / 'tests' / 'hw')
    for folder in folders:
        folder.mkdir(parents=True)
        (folder / 'selftest.md').write_text(TAGS)
    RT.DOCS_DIR = docs
    return tuple(RT.parse_test_tags(folder / 'selftest.md') for folder in folders)


def check_rule(RT, tut):
    """1. The rule, on a tutorial page."""
    single, spencer, bishop, seep = tut[0], tut[1], tut[2], tut[3]
    miss = RT._primary_as_written_miss(single, 0.918, 0.9186)
    check("0.918 written, 0.9186 computed: fails", bool(miss), miss or "passed")
    check("the message says what the page prints, the run gives, and it rounds to",
          bool(miss) and all(s in miss for s in ("0.918", "0.91860", "0.919")), miss)
    check("0.918 written, 0.9184 computed: passes",
          RT._primary_as_written_miss(single, 0.918, 0.9184) is None)
    check("fs_spencer=0.918 row: 0.9186 fails",
          bool(RT._primary_as_written_miss(spencer, 0.918, 0.9186)))
    miss = RT._primary_as_written_miss(bishop, 1.21, 1.2106)
    check("fs_bishop=1.210 row: 1.2106 fails (the trailing zero is a decimal)",
          bool(miss) and "fs_bishop=1.210" in miss, miss or "passed")
    check("fs_bishop=1.210 row: 1.2104 passes",
          RT._primary_as_written_miss(bishop, 1.21, 1.2104) is None)
    check("4.916e-07 written: 4.9166e-07 fails",
          bool(RT._primary_as_written_miss(seep, 4.916e-07, 4.9166e-07)))
    check("4.916e-07 written: 4.9163e-07 passes",
          RT._primary_as_written_miss(seep, 4.916e-07, 4.9163e-07) is None)
    check("one element of a list (1.0578 at t=50): 1.05786 fails",
          bool(RT._as_written_miss(single, 'expected (t=50)', 1.05786, written='1.0578')))
    check("one element of a list (1.0578 at t=50): 1.05784 passes",
          RT._as_written_miss(single, 'expected (t=50)', 1.05784, written='1.0578') is None)


def check_scope(RT, ver, sample, private):
    """2. Every page under docs/ is read at its written precision; a page outside
    docs/ (the private course rows) keeps the tolerance alone."""
    miss = RT._primary_as_written_miss(ver[0], 0.918, 0.9186)
    check("the same tag on a verification page fails 0.9186", bool(miss),
          miss or "passed")
    check("... and passes 0.9184",
          RT._primary_as_written_miss(ver[0], 0.918, 0.9184) is None)
    check("the same tag on a sample page two folders under docs/ fails 0.9186",
          bool(RT._primary_as_written_miss(sample[0], 0.918, 0.9186)))
    check("a fs_bishop=1.210 row on a verification page fails 1.2106",
          bool(RT._primary_as_written_miss(ver[2], 1.21, 1.2106)))
    check("the same tag on a private course page (outside docs/) passes 0.9186",
          RT._primary_as_written_miss(private[0], 0.918, 0.9186) is None)


def check_opt_out(RT, tut):
    """3. The opt-out, its reason, and the refusals."""
    no_reason, with_reason, misspelt, orphan = tut[4], tut[5], tut[6], tut[7]
    dispatched = []
    real_dispatch = RT._dispatch_test

    def fake_dispatch(test):
        dispatched.append(test)
        return 0.9186, None
    RT._dispatch_test = fake_dispatch
    try:
        computed, err, _ann = RT.run_test(no_reason)
        check("precision=tolerance without a reason is refused",
              computed is None and bool(err) and "precision_reason" in err, err)
        check("... before anything is solved", not dispatched)
        for name, test in (("a misspelt precision value", misspelt),
                           ("a reason with no opt-out", orphan)):
            computed, err, _ann = RT.run_test(test)
            check(f"{name} is refused", computed is None and bool(err), err)
        check("... none of them solved", not dispatched)
        computed, err, _ann = RT.run_test(with_reason)
        check("an opt-out with its reason runs", err is None and len(dispatched) == 1, err)
    finally:
        RT._dispatch_test = real_dispatch
    check("an opt-out with its reason passes 0.9186 inside its tolerance",
          RT._primary_as_written_miss(with_reason, 0.918, 0.9186) is None)


def check_real_pages(RT):
    """3. Every tag on every real page under docs/: a well-formed opt-out or none."""
    bad, n = [], 0
    for md in sorted((Path(_ROOT) / 'docs').rglob('*.md')):
        for t in RT.parse_test_tags(md):
            n += 1
            why = RT._precision_refusal(t)
            if why:
                bad.append(f"{md.relative_to(_ROOT)}: "
                           f"{t.get('benchmark') or t.get('type')}: {why}")
    check(f"every docs tag's precision opt-out is well formed ({n} rows)",
          not bad and n > 0, "; ".join(bad[:3]))


def check_fast_kernel(RT, tut):
    """4. The fast-kernel pass reads the page's written precision."""
    row = dict(tut[8], _default_tol=0.01)
    real = (RT.run_fem_test, RT._fast_kernel_available)
    RT.run_fem_test = lambda test, fast_kernel=None: (1.7334 if fast_kernel else 1.7302, None)
    RT._fast_kernel_available = lambda: True
    try:
        computed, err, annotation = RT._run_fem_ssrm_bracket(row)
    finally:
        RT.run_fem_test, RT._fast_kernel_available = real
    check("a kernel reading of 1.7334 against 1.730 falls through to the reference",
          err is None and computed == 1.7302 and annotation[0] == 'fallback',
          f"computed {computed}, {annotation}")


PROSE = """# Synthetic tutorial

## Results

Open `files/synthetic.xlsx` and run the circular search.

| OMS | Bishop | Spencer |
|---|---|---|
| 1.215 | {bishop} | 1.282 |

<!-- test: file=files/synthetic.xlsx, type=circular_search, num_slices=40, fs_oms=1.215, fs_bishop=1.250, fs_spencer=1.282 -->
"""


def check_prose(RT, root):
    """5. The restatement row fails on a page that disagrees with its own tag."""
    from tools.verification_checks import tutorials
    page = Path(root) / 'docs' / 'tutorials' / 'synthetic.md'
    page.parent.mkdir(parents=True, exist_ok=True)
    real = (tutorials.doc_locks, tutorials.tutorial_pages)
    # The row sweeps the real docs; point the sweep at this one page instead.
    tutorials.doc_locks = lambda repo=None: real[0](root)
    tutorials.tutorial_pages = lambda repo=None: [str(page)]
    try:
        verdicts = {}
        for bishop in ('1.251', '1.250', '1.25'):
            page.write_text(PROSE.format(bishop=bishop))
            verdicts[bishop] = RT.run_tutorial_restatements_test({})
    finally:
        tutorials.doc_locks, tutorials.tutorial_pages = real
    computed, err = verdicts['1.251']
    check("Bishop printed 1.251 under fs_bishop=1.250 fails the row",
          computed is None and bool(err) and "L9 prints 1.251" in err, err or "passed")
    for bishop in ('1.250', '1.25'):
        computed, err = verdicts[bishop]
        check(f"Bishop printed {bishop} passes the row", err is None, err)


def main():
    print("=" * 72)
    print("Printed-precision rule (every page under docs/)")
    print("=" * 72)
    import run_tests as RT
    saved = RT.DOCS_DIR
    with tempfile.TemporaryDirectory() as root:
        try:
            tut, ver, sample, private = _pages(RT, root)
            check_rule(RT, tut)
            check_scope(RT, ver, sample, private)
            check_opt_out(RT, tut)
            check_fast_kernel(RT, tut)
        finally:
            RT.DOCS_DIR = saved
    with tempfile.TemporaryDirectory() as root:
        check_prose(RT, root)
    check_real_pages(RT)
    print("\n" + "=" * 72)
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}): " + ", ".join(FAILURES))
        return 1
    print("All printed-precision checks passed.")
    return 0


def run():
    """Failures as a list, for run_tests.py."""
    del FAILURES[:]
    main()
    return list(FAILURES)


if __name__ == '__main__':
    sys.exit(main())
