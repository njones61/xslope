"""Standing check: every tutorial page walked as written (owner decision 591a).

The 1.1.0 tutorial pass found 26 of 31 pages that a first-time reader could not
follow to the end, and only a walk found them: a profile line one row short, a
mesh the page said to keep that Studio had dropped, a Log line the page quoted
from a run nobody could repeat. This check makes that walk standing. For each
page under ``docs/tutorials/`` it runs ``test/tutorial_walks/<page>.py``, which
does what the page tells the reader to do, in order:

* Studio steps run in the real Studio window, offscreen, with Studio's settings
  replaced by a scratch file (``tools/studio_settings.py``) and on copies of the
  page's files in a temporary folder;
* script and API steps run through the engine.

A page fails when a step cannot be done as written (a control or label the page
names is not there, a paste is refused, a dialog does not open or does not
behave as the page says), and when a number the page prints does not round, at
the precision it is printed, to what the run gives. The machinery is in
``test/tutorial_walks/harness.py``; each page module says what is and is not
scripted on its page.

Every page is walked in its own Python process, so one page's Studio state (or a
crash) cannot reach the next, and each page's wall time is reported. The full set
takes hours (strength reduction and transient seepage runs among them) and rides
the release gate; ``--page`` walks the named pages alone.

    python test/tutorial_walk_check.py --page lem01 --page seep04
    python test/tutorial_walk_check.py --list
    python test/tutorial_walk_check.py                # every page

``run(pages=None)`` returns the failure strings (the suite's contract).
"""
from __future__ import annotations

import argparse
import contextlib
import importlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
for _p in (_HERE, _REPO):
    if _p not in sys.path:
        sys.path.insert(0, _p)

TUTORIALS = os.path.join(_REPO, "docs", "tutorials")
WALKS = os.path.join(_HERE, "tutorial_walks")

#: How long one page may take before its walk is stopped and failed. A page
#: module may set its own ``TIMEOUT`` (seconds).
DEFAULT_TIMEOUT = 3 * 3600


def page_files():
    """Every tutorial page, by file name."""
    return sorted(f for f in os.listdir(TUTORIALS)
                  if f.endswith(".md") and f != "index.md")


def module_name(page):
    """``lem01_simple_embankment.md`` -> ``lem01``; a page without a numbered
    prefix keeps its whole stem (``building_models``)."""
    stem = os.path.splitext(page)[0]
    head = stem.split("_")[0]
    return head if any(ch.isdigit() for ch in head) else stem


def select(names):
    """The pages ``names`` pick (a module name or any part of a file name)."""
    pages = page_files()
    if not names:
        return pages
    out = []
    for n in names:
        hits = [p for p in pages if module_name(p) == n or n in p]
        if not hits:
            raise SystemExit("no tutorial page matches %r (pages: %s)"
                             % (n, ", ".join(module_name(p) for p in pages)))
        out += [h for h in hits if h not in out]
    return out


# --------------------------------------------------------------------------- #
# One page, in this process
# --------------------------------------------------------------------------- #
def walk_page(page, verbose=False, keep=False):
    """Walk one page here and now. Returns a dict: failures, notes, seconds,
    checked (the number of behaviors and numbers held to the page)."""
    from tutorial_walks import harness
    from tools import studio_settings

    t0 = time.time()
    name = module_name(page)
    path = os.path.join(WALKS, name + ".py")
    if not os.path.exists(path):
        return {"page": page, "failures": ["%s: no walk (test/tutorial_walks/%s.py)"
                                           % (page, name)],
                "notes": [], "seconds": 0.0, "checked": 0}
    workdir = tempfile.mkdtemp(prefix="xslope_walk_%s_" % name)
    walk = None
    try:
        with studio_settings.redirected(prefix="tutorial_walk"):
            mod = importlib.import_module("tutorial_walks." + name)
            if getattr(mod, "PAGE", None) != page:
                raise SystemExit("tutorial_walks/%s.py walks %r, not %r"
                                 % (name, getattr(mod, "PAGE", None), page))
            walk = harness.Walk(mod, workdir, verbose=verbose)
            try:
                mod.walk(walk)
            except harness._Abort:
                pass
            except Exception as exc:          # outside any step: the walk's own bug
                import traceback
                walk.failures.append("%s: the walk stopped outside a step: %s: %s"
                                     % (page, type(exc).__name__, exc))
                if verbose:
                    traceback.print_exc()
            finally:
                walk.close()
    finally:
        if keep:
            print("   work folder kept: %s" % workdir)
        else:
            shutil.rmtree(workdir, ignore_errors=True)
    return {"page": page, "failures": walk.failures if walk else [],
            "notes": walk.notes if walk else [], "seconds": time.time() - t0,
            "checked": walk.checked if walk else 0}


# --------------------------------------------------------------------------- #
# Every page, one process each
# --------------------------------------------------------------------------- #
def _timeout(page):
    try:
        mod = importlib.import_module("tutorial_walks." + module_name(page))
        return float(getattr(mod, "TIMEOUT", DEFAULT_TIMEOUT))
    except Exception:
        return DEFAULT_TIMEOUT


def walk_in_child(page, verbose=False):
    """Walk one page in a fresh Python process; its result as walk_page's."""
    fd, out = tempfile.mkstemp(prefix="xslope_walk_", suffix=".json")
    os.close(fd)
    cmd = [sys.executable, "-P", os.path.abspath(__file__), "--child", page,
           "--json", out] + (["--verbose"] if verbose else [])
    env = dict(os.environ)
    env["PYTHONPATH"] = _REPO + (os.pathsep + env["PYTHONPATH"]
                                 if env.get("PYTHONPATH") else "")
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("XSLOPE_NO_UPDATE_CHECK", "1")
    t0 = time.time()
    timeout = _timeout(page)
    try:
        proc = subprocess.run(cmd, env=env, timeout=timeout,
                              stdout=None if verbose else subprocess.PIPE,
                              stderr=subprocess.STDOUT, text=True)
        tail = (proc.stdout or "")[-1500:] if not verbose else ""
        code = proc.returncode
    except subprocess.TimeoutExpired as exc:
        tail = ((exc.stdout or b"").decode(errors="replace")
                if isinstance(exc.stdout, bytes) else (exc.stdout or ""))[-800:]
        os.unlink(out)
        return {"page": page, "seconds": time.time() - t0, "checked": 0, "notes": [],
                "failures": ["%s: the walk did not finish in %d s" % (page, timeout)]}
    try:
        with open(out, encoding="utf-8") as fh:
            res = json.load(fh)
    except Exception:
        res = {"page": page, "seconds": time.time() - t0, "checked": 0, "notes": [],
               "failures": ["%s: the walk crashed (exit %s): %s"
                            % (page, code, tail.strip().replace("\n", " | "))]}
    finally:
        with contextlib.suppress(OSError):
            os.unlink(out)
    return res


def run(pages=None, verbose=False, results=None):
    """Walk ``pages`` (default: every tutorial page), each in its own process.
    Returns the failure strings; ``results``, if a list, receives each page's
    result dict (failures, notes, seconds, checked)."""
    failures = []
    for page in select(pages):
        res = walk_in_child(page, verbose=verbose)
        if results is not None:
            results.append(res)
        failures += res["failures"]
        print("   %-34s %7.0f s  %4d checked  %s"
              % (page, res["seconds"], res.get("checked", 0),
                 "ok" if not res["failures"] else "%d FAILED" % len(res["failures"])),
              flush=True)
    return failures


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--page", action="append", default=[],
                    help="walk this page (module name such as lem01, or part of "
                         "its file name); repeatable; default every page")
    ap.add_argument("--list", action="store_true", help="list the pages and walks")
    ap.add_argument("--inline", action="store_true",
                    help="walk in this process (debugging; no per-page isolation)")
    ap.add_argument("--keep", action="store_true",
                    help="keep each page's work folder (with --inline)")
    ap.add_argument("--verbose", action="store_true", help="print every step and check")
    ap.add_argument("--times", help="write each page's result (seconds, checks, "
                                    "failures) to this JSON file")
    ap.add_argument("--child", help=argparse.SUPPRESS)
    ap.add_argument("--json", help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    if args.child:
        res = walk_page(args.child, verbose=args.verbose)
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(res, fh)
        sys.stdout.flush()
        os._exit(0)          # Qt teardown after a long walk can hang; the result is out

    if args.list:
        for p in page_files():
            m = module_name(p)
            have = os.path.exists(os.path.join(WALKS, m + ".py"))
            print("%-36s %s" % (p, ("tutorial_walks/%s.py" % m) if have else "NO WALK"))
        return 0

    import xslope
    print("xslope from %s" % xslope.__file__)
    results = []
    if args.inline:
        failures = []
        for page in select(args.page):
            res = walk_page(page, verbose=args.verbose, keep=args.keep)
            results.append(res)
            failures += res["failures"]
            print("   %-34s %7.0f s  %4d checked  %s"
                  % (page, res["seconds"], res["checked"],
                     "ok" if not res["failures"] else "%d FAILED" % len(res["failures"])))
    else:
        failures = run(args.page or None, verbose=args.verbose, results=results)
    for r in results:
        for n in r.get("notes", []):
            print("   note  %s: %s" % (r["page"], n))
    if args.times:
        with open(args.times, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=1)
    if failures:
        print("\n%d failure(s):" % len(failures))
        for f in failures:
            print("  - %s" % f)
    else:
        print("\nall %d page(s) walk as written" % len(results))
    rc = 1 if failures else 0
    if args.inline:
        # Studio ran in this process: a worker thread a page left behind (the
        # assistant's model-list refresh, say) aborts Qt's teardown at exit, so
        # the result is out and the process leaves without it.
        sys.stdout.flush()
        os._exit(rc)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
