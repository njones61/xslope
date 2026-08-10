"""Standing checks on project packaging: collecting a project into one ``.xslz``
file, and taking it apart again.

A project is a SET — the workbook plus the sidecars that basename-key to it — and the
package exists so that set can travel as one file. Every way that can go wrong is
quiet, which is why each is checked here rather than left to a hands-on session:

  A. THE SET IS THE RIGHT SET. Sidecars are collected by the ``{base}_*`` convention
     rather than a list of known suffixes, so a sidecar a future solver invents
     travels without anyone remembering to add it. Two workbooks in one folder, one
     name extending the other, is where attribution gets decided, and it has to be
     decided the way the LOADERS decide it — a file belongs to the workbook whose
     loader reads it by that exact name. Reading it as "the longest workbook name
     that prefixes the file" instead is wrong in both directions on real corpus
     pairs: ``vp091.xlsx``'s own FEM results are ``vp091_fem_nodes.csv`` and friends
     (``vp091_fem.xlsx``'s would be ``vp091_fem_fem_nodes.csv``), so longest-prefix
     packs ``vp091`` as a bare workbook and silently drops its entire solved FEM set;
     and read the other way, ``vp091_fem`` walks off with results that are not its.
     Meanwhile ``xslope_earth_dam1_vg_mesh.json`` really is the ``_vg`` project's
     mesh and must stay out of ``xslope_earth_dam1``. And because that attribution
     is decided by a set of suffixes, the solvers are read to check they write
     nothing that set has not heard of — a name drifting out of it is silent until
     two workbooks share a folder, and then it is a misattributed results file.
  B. WHAT WENT IN COMES OUT. Every packed file is restored byte for byte, and the
     model the extracted workbook loads to is the model the original loads to —
     mesh and pore-pressure field included, since those live in the sidecars and
     nowhere else.
  C. A SINGLE-FILE PROJECT PACKAGES TOO. Most samples are a workbook and nothing
     else. Zipping one file is uniform rather than clever: one rule covers every
     project, and a project that gains a sidecar tomorrow does not change type.
  D. THE COLLISION IS A REFUSAL, NOT A GUESS. The folder a package would unpack into
     may already hold the user's edits. The library cannot ask, so it raises — and
     the message names all three ways forward, because an error that states only
     the problem sends the reader to the source.
  E. ``overwrite=True`` AND ``dest=`` DO WHAT THEY SAY. The first replaces the
     package's own files and leaves everything else in the folder alone; the second
     extracts where it is told.
  F. A PACKAGE IS A FLAT ZIP OF ONE PROJECT. An archive with folders, with paths
     that would escape the destination, with a drive-relative name that would escape
     it on Windows only, or with two workbooks is refused by name rather than
     extracted.
  H. WHAT WILL NOT FIT IS REFUSED. Compression ratio is unbounded, so the size of a
     package says nothing about what it writes: a megabyte of one repeated byte
     expands past a gigabyte. The ceiling is enforced on the sizes the archive
     declares AND on the bytes that arrive, and a refusal leaves nothing behind.
  G. STUDIO OPENS AND EXPORTS ONE. Open reaches the normal open path with the
     EXTRACTED workbook, so recent files and the window title point at the loose
     .xlsx and never at the package; Export writes the whole set beside the project;
     and the destination dialog never silently reuses or overwrites a folder that is
     already there.

Skips its Studio leg cleanly when PySide6 is absent; the rest run either way.
"""
import contextlib
import io
import math
import os
import re
import shutil
import sys
import tempfile
import zipfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import xslope
from xslope.fileio import load_slope_data
from xslope.package import (PACKAGE_EXT, SIDECAR_SUFFIXES, package_contents,
                            project_files)

#: A real project WITH sidecars: a mesh and a solved steady seepage field its
#: materials actually read (``u = seep``), so a package that lost a sidecar loses a
#: visible half of the model rather than a file nobody opens.
WITH_SIDECARS = os.path.join(_REPO, "docs/seep/files/xslope_rface_SEEP_KEY.xlsx")
#: A project whose folder also holds ``xslope_earth_dam1_vg.xlsx`` — a separate
#: project whose name extends this one's — which is the over-collection hazard.
NEIGHBOURED = os.path.join(_REPO, "docs/seep/files/xslope_earth_dam1.xlsx")
#: The neighbouring project that must NOT be swept into the one above.
NEIGHBOUR = os.path.join(_REPO, "docs/seep/files/xslope_earth_dam1_vg.xlsx")
#: A project that is a workbook and nothing else.
SINGLE_FILE = os.path.join(_REPO, "docs/inputs/slope/xslope_simple1.xlsx")

#: The corpus is read-only here: every check copies what it needs into a temporary
#: folder and works there, so nothing in docs/ is ever written to.
_TMPDIRS = []


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _tmp():
    d = tempfile.mkdtemp(prefix="xslz_check_")
    _TMPDIRS.append(d)
    return d


def _copy_project(*models):
    """Copy each model's whole file set (by the same convention) into a fresh
    temporary folder, and return the folder. Copies, never moves: the corpus is a
    fixture."""
    folder = _tmp()
    for model in models:
        for path in project_files(model):
            shutil.copy2(path, folder)
    return folder


def _same_bytes(a, b):
    with open(a, "rb") as fa, open(b, "rb") as fb:
        return fa.read() == fb.read()


def _eq(a, b):
    """Scalar equality with a float tolerance; everything else compared verbatim."""
    if str(a) == str(b):
        return True
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return False
    return math.isclose(fa, fb, rel_tol=1e-9, abs_tol=1e-12)


def _diff(a, b, path=""):
    """Recursively compare two loaded models; return a list of mismatch paths."""
    import numpy as np

    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a_arr, b_arr = np.asarray(a), np.asarray(b)
        if a_arr.shape != b_arr.shape:
            return [f"{path}: shape {a_arr.shape} vs {b_arr.shape}"]
        if a_arr.dtype.kind in "fc" and not np.allclose(a_arr, b_arr, equal_nan=True):
            return [f"{path}: array values differ"]
        if a_arr.dtype.kind not in "fc" and not np.array_equal(a_arr, b_arr):
            return [f"{path}: array values differ"]
        return []
    if isinstance(a, dict):
        if not isinstance(b, dict):
            return [f"{path}: dict vs {type(b).__name__}"]
        out = []
        for k in set(a) | set(b):
            if k not in a or k not in b:
                out.append(f"{path}.{k}: missing on one side")
            else:
                out += _diff(a[k], b[k], f"{path}.{k}")
        return out
    if isinstance(a, (list, tuple)):
        if not isinstance(b, (list, tuple)) or len(a) != len(b):
            blen = len(b) if hasattr(b, "__len__") else "?"
            return [f"{path}: len {len(a)} vs {blen}"]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += _diff(x, y, f"{path}[{i}]")
        return out
    return [] if _eq(a, b) else [f"{path}: {a!r} != {b!r}"]


# ------------------------------------------------------------------ A. the file set
def test_file_set():
    """What a project's package holds, and — just as much — what it does not."""
    fails = []
    folder = _copy_project(NEIGHBOURED, NEIGHBOUR)
    book = os.path.join(folder, os.path.basename(NEIGHBOURED))
    names = [os.path.basename(p) for p in project_files(book)]

    if names[0] != os.path.basename(NEIGHBOURED):
        fails.append(f"the workbook is not the first file collected: {names}")
    for want in ("xslope_earth_dam1_mesh.json", "xslope_earth_dam1_seep.csv"):
        if want not in names:
            fails.append(f"{want} was not collected with its workbook: {names}")
    # The neighbouring project, whose name extends this one's, stays out — workbook
    # and sidecars both.
    for unwanted in ("xslope_earth_dam1_vg.xlsx", "xslope_earth_dam1_vg_mesh.json",
                     "xslope_earth_dam1_vg_seep.csv"):
        if unwanted in names:
            fails.append(f"{unwanted} belongs to another project but was collected")

    # The convention is a glob, not a list of known suffixes: a sidecar no solver has
    # invented yet still travels.
    invented = os.path.join(folder, "xslope_earth_dam1_future_results.csv")
    with open(invented, "w") as fh:
        fh.write("x\n")
    if os.path.basename(invented) not in [
            os.path.basename(p) for p in project_files(book)]:
        fails.append("a sidecar under the naming convention was not collected — the "
                     "convention has become a list of known suffixes")

    # The neighbour keeps its own set.
    vg = [os.path.basename(p)
          for p in project_files(os.path.join(folder, os.path.basename(NEIGHBOUR)))]
    if sorted(vg) != sorted(["xslope_earth_dam1_vg.xlsx", "xslope_earth_dam1_vg_mesh.json",
                             "xslope_earth_dam1_vg_seep.csv"]):
        fails.append(f"the neighbouring project's own set is wrong: {vg}")
    return fails


# ------------------------------------------- A2. attribution the loaders' way, both ways
def _fem_pair(folder):
    """The corpus shape that decides this: ``X.xlsx`` and ``X_fem.xlsx`` in one
    folder, where X's FEM results are ``X_fem_*`` and X_fem's own are ``X_fem_fem_*``
    (docs/verification/files/rocscience/vp091 and vp027 are both really like this).
    Synthesized rather than copied so the check states the shape it is about."""
    made = {}
    for name in ("model.xlsx", "model_fem.xlsx",
                 "model_fem_nodes.csv", "model_fem_meta.json",
                 "model_fem_failure_nodes.csv",
                 "model_fem_fem_nodes.csv", "model_fem_fem_meta.json",
                 "model_notes.txt"):
        path = os.path.join(folder, name)
        with open(path, "w") as fh:
            fh.write(name)
        made[name] = path
    return made


def test_attribution():
    """Which of two workbooks in one folder owns a results file.

    The rule is the loaders' rule: a file belongs to the workbook whose loader reads
    it by that exact name — its stem plus one of the suffixes the writers write. The
    tempting shorthand, "the longest workbook name that prefixes it", is wrong in
    both directions here, and wrong in a way that loses a solved FEM set without a
    word."""
    fails = []
    folder = _tmp()
    _fem_pair(folder)

    got = [os.path.basename(p)
           for p in project_files(os.path.join(folder, "model.xlsx"))]
    want = ["model.xlsx", "model_fem_failure_nodes.csv", "model_fem_meta.json",
            "model_fem_nodes.csv", "model_notes.txt"]
    if sorted(got) != sorted(want):
        fails.append(f"model.xlsx collected {sorted(got)}, expected {sorted(want)}")

    got_fem = [os.path.basename(p)
               for p in project_files(os.path.join(folder, "model_fem.xlsx"))]
    want_fem = ["model_fem.xlsx", "model_fem_fem_meta.json", "model_fem_fem_nodes.csv"]
    if sorted(got_fem) != sorted(want_fem):
        fails.append(f"model_fem.xlsx collected {sorted(got_fem)}, "
                     f"expected {sorted(want_fem)}")

    # LOCALITY: what a project packs must not depend on which neighbours happen to be
    # in the folder. Its own set, alone in a folder, is the same set.
    alone = _tmp()
    for name in want:
        shutil.copy2(os.path.join(folder, name), alone)
    solo = [os.path.basename(p)
            for p in project_files(os.path.join(alone, "model.xlsx"))]
    if sorted(solo) != sorted(want):
        fails.append(f"alone in a folder model.xlsx collects {sorted(solo)} — its set "
                     f"depends on its neighbours")

    # MUTATION: reinstate the longest-prefix rule and this must go red. A guard that
    # cannot fail is not a guard.
    import xslope.package as P

    original = P.companion_of
    try:
        P.companion_of = lambda name, stems: next(
            (s for s in stems if name.startswith(s + "_")), None)
        mutated = [os.path.basename(p)
                   for p in project_files(os.path.join(folder, "model.xlsx"))]
        if sorted(mutated) == sorted(want):
            fails.append("the longest-prefix rule passes this check — the check does "
                         "not actually pin the attribution rule")
    finally:
        P.companion_of = original

    # And the real corpus pair the defect was found on, if it is in this checkout.
    vp091 = os.path.join(_REPO, "docs/verification/files/rocscience/vp091.xlsx")
    if os.path.exists(vp091):
        real = [os.path.basename(p) for p in project_files(vp091)]
        fem_results = [n for n in real if n.startswith("vp091_fem")]
        if len(fem_results) != 8:
            fails.append(f"vp091.xlsx collected {len(fem_results)} of its 8 FEM "
                         f"result files: {real}")
        if "vp091_fem.xlsx" in real:
            fails.append("vp091's package carries the neighbouring workbook")
        # And it survives the round trip as something the FEM loader can still read.
        # dest= throughout: the corpus is read here, never written to.
        from xslope.fem import import_fem_meta

        before = import_fem_meta(os.path.splitext(vp091)[0])
        pkg = xslope.pack(vp091, dest=_tmp())
        book = xslope.unpack(pkg, dest=os.path.join(_tmp(), "vp091"))
        after = import_fem_meta(os.path.splitext(book)[0])
        if before is None:
            fails.append("the vp091 fixture no longer carries an FEM solution")
        elif after is None:
            fails.append("vp091's FEM solution did not survive the round trip — the "
                         "package arrived without the results")
        elif _diff(before, after, "fem_meta"):
            fails.append("vp091's FEM metadata came back changed: "
                         + "; ".join(_diff(before, after, "fem_meta")[:3]))
    return fails


# ------------------------------------------- A3. the writers spell what the set knows
#: The solver modules that write result files beside a workbook, and so decide what a
#: project's file set actually contains. They are SCANNED rather than imported because
#: neither spells its sidecar names as a module constant — xslope.seep builds them
#: inline (``f"{base}_tseep.csv"``, and a ``("_tseep.csv", "_tseep_meta.json")`` tuple
#: inside a function), and xslope.fem builds them from ``output_stem``. Reading the
#: source is the only way to see the names those modules really use.
SIDECAR_WRITERS = ("xslope/seep.py", "xslope/fem.py")

#: A suffix as it appears in those files: a quoted literal, or the tail of an f-string
#: (or a docstring) that follows a ``{base}`` / ``{stem}`` / ``{output_stem…}`` slot.
_SUFFIX_LITERAL = re.compile(r"""["'](_[a-z0-9_]+\.(?:csv|json))["']""")
_SUFFIX_INTERP = re.compile(
    r"""\{[A-Za-z_.]*(?:base|stem)[^}]*\}(_[a-z0-9_]+\.(?:csv|json))""")


def test_writer_suffixes():
    """The names the solvers write must be names the file set knows.

    ``SIDECAR_SUFFIXES`` decides which of two workbooks in one folder owns a results
    file. A solver that starts writing a name that set has never heard of does not
    break loudly: the file still travels (an unrecognized ``{base}_*`` sibling packs),
    and everything looks right until two workbooks share a folder and one name extends
    the other — and then it is attributed to the wrong project, silently, which is
    exactly the failure this module's A2 leg exists for. The corpus already holds that
    shape for transient seepage: ``xslope_johnson_res.xlsx`` sits beside
    ``xslope_johnson_res_tseep.xlsx``.

    So this leg reads the writers and asserts they spell nothing the set does not
    carry. It cannot be replaced by a comment in those files, and it fires on the
    commit that introduces the drift rather than on the corpus that trips over it."""
    fails = []
    found = {}
    for rel in SIDECAR_WRITERS:
        path = os.path.join(_REPO, rel)
        if not os.path.exists(path):
            fails.append(f"{rel} is missing — this guard is scanning nothing")
            continue
        with open(path, encoding="utf-8") as fh:
            src = fh.read()
        for suffix in set(_SUFFIX_LITERAL.findall(src)) | set(_SUFFIX_INTERP.findall(src)):
            found.setdefault(suffix, set()).add(rel)

    # A scan that finds nothing would pass forever. These two are the names the
    # transient writer has used since it was written; if they are gone, the scan is
    # looking at the wrong thing, not at a module that stopped writing sidecars.
    for anchor in ("_tseep.csv", "_tseep_meta.json"):
        if anchor not in found:
            fails.append(f"the scan of {', '.join(SIDECAR_WRITERS)} did not find "
                         f"{anchor} — the writers have been restructured and this "
                         f"guard is no longer reading them")

    known = set(SIDECAR_SUFFIXES)
    for suffix, where in sorted(found.items()):
        if suffix not in known:
            fails.append(
                f"{' and '.join(sorted(where))} writes '{{base}}{suffix}' beside a "
                f"workbook, but SIDECAR_SUFFIXES in xslope/package.py does not carry "
                f"it. Until it does, that file is attributed to the wrong project "
                f"wherever two workbooks share a folder and one name extends the "
                f"other (the corpus has that shape: xslope_johnson_res.xlsx beside "
                f"xslope_johnson_res_tseep.xlsx). Add '{suffix}' to SIDECAR_SUFFIXES.")
    return fails


# ------------------------------------------------------- B. the round trip, with sidecars
def test_round_trip():
    """Pack a project with sidecars, unpack it somewhere fresh, and ask whether it is
    the same project."""
    fails = []
    folder = _copy_project(WITH_SIDECARS)
    book = os.path.join(folder, os.path.basename(WITH_SIDECARS))

    pkg = xslope.pack(book)
    if pkg != os.path.splitext(book)[0] + PACKAGE_EXT:
        fails.append(f"the package was not written beside the workbook: {pkg}")
    if package_contents(pkg)[0] != os.path.basename(book):
        fails.append("the workbook is not the first entry in the archive")

    dest = os.path.join(_tmp(), "unpacked")
    book2 = xslope.unpack(pkg, dest=dest)
    if os.path.dirname(book2) != dest or not os.path.isfile(book2):
        fails.append(f"unpack did not return the extracted workbook: {book2}")

    packed = {os.path.basename(p) for p in project_files(book)}
    got = set(os.listdir(dest))
    if got != packed:
        fails.append(f"the extracted folder holds {sorted(got)}, packed {sorted(packed)}")
    for name in packed & got:
        if not _same_bytes(os.path.join(folder, name), os.path.join(dest, name)):
            fails.append(f"{name} came out of the package changed")

    before = _quiet(load_slope_data, book)
    after = _quiet(load_slope_data, book2)
    diffs = _diff(before, after, "model")
    if diffs:
        fails.append("the unpacked project loads to a different model: "
                     + "; ".join(diffs[:5]))
    # The mesh and the pore-pressure field live in the sidecars and nowhere else, so
    # they are the half of the model a package that dropped them would lose silently.
    if after.get("mesh") is None or not len(after["mesh"].get("nodes", [])):
        fails.append("the unpacked project carries no mesh")
    if after.get("seep_u") is None or not len(after["seep_u"]):
        fails.append("the unpacked project carries no pore-pressure field")
    return fails


# --------------------------------------------------------- C. the single-file project
def test_single_file():
    """A workbook with no sidecars packages and round-trips like any other project."""
    fails = []
    folder = _copy_project(SINGLE_FILE)
    book = os.path.join(folder, os.path.basename(SINGLE_FILE))
    if len(project_files(book)) != 1:
        fails.append("the single-file fixture has grown sidecars; pick another")

    pkg = xslope.pack(book)
    names = package_contents(pkg)
    if names != [os.path.basename(book)]:
        fails.append(f"a single-file package holds {names}")

    book2 = xslope.unpack(pkg, dest=os.path.join(_tmp(), "one"))
    if not _same_bytes(book, book2):
        fails.append("the single-file workbook came out of the package changed")
    diffs = _diff(_quiet(load_slope_data, book), _quiet(load_slope_data, book2), "model")
    if diffs:
        fails.append("the unpacked single-file project loads differently: "
                     + "; ".join(diffs[:5]))
    return fails


# ------------------------------------------------------------------ D. the collision
def test_collision():
    """The destination already exists. The library cannot ask, so it refuses — and
    says what to do about it."""
    fails = []
    folder = _copy_project(WITH_SIDECARS)
    book = os.path.join(folder, os.path.basename(WITH_SIDECARS))
    pkg = xslope.pack(book)

    dest = os.path.join(_tmp(), "twice")
    xslope.unpack(pkg, dest=dest)
    try:
        xslope.unpack(pkg, dest=dest)
    except FileExistsError as exc:
        message = str(exc)
    else:
        return ["unpacking over an existing folder was allowed"]

    for phrase, way_out in (
            (os.path.basename(book), "the workbook already extracted"),
            ("dest=", "extracting somewhere else"),
            ("overwrite=True", "writing over the folder")):
        if phrase not in message:
            fails.append(f"the collision message does not name {way_out} "
                         f"({phrase!r}): {message}")

    # load_slope_data on a package refuses the same way rather than loading whatever
    # happens to be in the folder already.
    stray = os.path.join(os.path.dirname(pkg),
                         os.path.splitext(os.path.basename(pkg))[0])
    os.makedirs(stray, exist_ok=True)
    try:
        _quiet(load_slope_data, pkg)
    except FileExistsError as exc:
        if "overwrite=True" not in str(exc):
            fails.append(f"load_slope_data's collision message is thinner: {exc}")
    else:
        fails.append("load_slope_data unpacked over an existing folder")

    # dest=/overwrite= are meaningless for a workbook, and saying so beats ignoring it.
    try:
        _quiet(load_slope_data, book, dest="/tmp/nowhere")
    except ValueError:
        pass
    except Exception as exc:
        fails.append(f"dest= on a workbook raised {exc!r}")
    else:
        fails.append("dest= was accepted on a workbook, where it does nothing")
    return fails


# --------------------------------------------------------- E. overwrite= and dest=
def test_overwrite_and_dest():
    """``overwrite=True`` replaces the package's own files and nothing else."""
    fails = []
    folder = _copy_project(WITH_SIDECARS)
    book = os.path.join(folder, os.path.basename(WITH_SIDECARS))
    pkg = xslope.pack(book)

    dest = os.path.join(_tmp(), "over")
    xslope.unpack(pkg, dest=dest)
    mesh = os.path.join(dest, "xslope_rface_SEEP_KEY_mesh.json")
    with open(mesh, "w") as fh:
        fh.write("clobbered")
    mine = os.path.join(dest, "my_notes.txt")
    with open(mine, "w") as fh:
        fh.write("keep me")

    book2 = xslope.unpack(pkg, dest=dest, overwrite=True)
    if not _same_bytes(os.path.join(folder, "xslope_rface_SEEP_KEY_mesh.json"), mesh):
        fails.append("overwrite=True did not restore a sidecar it had been asked to")
    if not os.path.exists(mine) or open(mine).read() != "keep me":
        fails.append("overwrite=True removed a file that was not the package's")
    if os.path.dirname(book2) != dest:
        fails.append("overwrite=True extracted somewhere other than dest")

    # The default destination is a folder named for the package, beside it.
    plain = os.path.join(_tmp(), "beside")
    os.makedirs(plain)
    shutil.copy2(pkg, plain)
    beside = xslope.unpack(os.path.join(plain, os.path.basename(pkg)))
    want = os.path.join(plain, os.path.splitext(os.path.basename(pkg))[0])
    if os.path.dirname(beside) != want:
        fails.append(f"the default destination is {os.path.dirname(beside)}, not {want}")

    # And load_slope_data("...xslz") does the whole thing.
    fresh = os.path.join(_tmp(), "loadme")
    os.makedirs(fresh)
    shutil.copy2(pkg, fresh)
    model = _quiet(load_slope_data, os.path.join(fresh, os.path.basename(pkg)))
    if model.get("mesh") is None:
        fails.append("load_slope_data on a package did not pick up the sidecars")

    # An existing package is left alone unless pack is told otherwise.
    try:
        xslope.pack(book)
    except FileExistsError:
        pass
    else:
        fails.append("pack wrote over an existing package without being asked")
    if xslope.pack(book, overwrite=True) != pkg:
        fails.append("pack(overwrite=True) wrote somewhere unexpected")

    # A dest that is not itself a package name is a FOLDER to write into, whether or
    # not it exists yet — otherwise a destination meant as a folder quietly becomes an
    # extensionless file.
    outbox = os.path.join(_tmp(), "outbox")
    into = xslope.pack(book, dest=outbox)
    if into != os.path.join(outbox, os.path.basename(pkg)):
        fails.append(f"pack(dest=<folder>) wrote {into}")
    named = xslope.pack(book, dest=os.path.join(outbox, "handover" + PACKAGE_EXT))
    if os.path.basename(named) != "handover" + PACKAGE_EXT:
        fails.append(f"pack(dest=<package path>) wrote {named}")
    return fails


# -------------------------------------------------------------- F. what is not a package
def test_not_a_package():
    """An archive that is not a project package is refused by name."""
    fails = []
    folder = _tmp()

    def archive(name, entries):
        path = os.path.join(folder, name)
        with zipfile.ZipFile(path, "w") as zf:
            for entry, data in entries:
                zf.writestr(entry, data)
        return path

    cases = [
        ("folders", [("project/model.xlsx", "x")]),
        ("escape", [("../model.xlsx", "x")]),
        ("two_books", [("a.xlsx", "x"), ("b.xlsx", "y")]),
        ("no_book", [("notes.txt", "x")]),
        # A drive-relative name. On POSIX "D:pwned.txt" is an ordinary file name and
        # every traversal test above passes it; on Windows it is a path on another
        # drive, and ntpath.join(dest, "D:pwned.txt") throws the destination away.
        # The refusal has to be made on the name, not on the platform.
        ("drive_letter", [("model.xlsx", "x"), ("D:pwned.txt", "x")]),
        ("unc_drive", [("model.xlsx", "x"), ("C:/Windows/system32/evil.dll", "x")]),
    ]
    for name, entries in cases:
        path = archive(name + PACKAGE_EXT, entries)
        try:
            package_contents(path)
        except ValueError as exc:
            if os.path.basename(path) not in str(exc):
                fails.append(f"the refusal of {name} does not name the file: {exc}")
        else:
            fails.append(f"{name} was accepted as a project package")
        # And the same refusal reaches anyone who tries to unpack or load it.
        try:
            xslope.unpack(path, dest=os.path.join(folder, name + "_out"))
        except ValueError:
            pass
        else:
            fails.append(f"unpack extracted {name}, which is not a project package")
    return fails


# ---------------------------------------------------------------- H. the expansion
def test_expansion():
    """A package that would not fit on disk is refused, twice, and writes nothing.

    Compression ratio is unbounded and invisible: a megabyte of archive holding one
    repeated byte expands past a gigabyte, and a package is exactly the thing a link,
    an email or a colleague hands a user. The ceiling is checked against the sizes the
    archive DECLARES and again against the bytes that actually arrive, because the
    declared size is whatever the person who built the archive typed.
    """
    from xslope.package import MAX_UNPACKED_BYTES, package_contents as _contents
    import xslope.package as P

    fails = []
    folder = _tmp()
    bomb = os.path.join(folder, "bomb" + PACKAGE_EXT)
    block = b"\0" * (1 << 23)
    chunks = int(MAX_UNPACKED_BYTES // len(block)) + 12      # comfortably over
    with zipfile.ZipFile(bomb, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("model.xlsx", b"PK\x03\x04")
        with zf.open("model_mesh.json", "w") as fh:
            for _ in range(chunks):
                fh.write(block)
    ratio = (len(block) * chunks) / os.path.getsize(bomb)
    if ratio < 100:
        fails.append(f"the fixture only compresses {ratio:.0f}:1, so it is not the "
                     f"hazard being checked")

    # 1. The declared sizes. Refused before a byte is extracted.
    out = os.path.join(folder, "bomb_out")
    try:
        _contents(bomb)
    except ValueError as exc:
        if "bomb" + PACKAGE_EXT not in str(exc):
            fails.append(f"the refusal does not name the package: {exc}")
    else:
        fails.append(f"a {ratio:.0f}:1 package was accepted by package_contents")
    try:
        xslope.unpack(bomb, dest=out)
    except ValueError:
        pass
    else:
        fails.append("unpack extracted a package that expands past the ceiling")
    if os.path.isdir(out) and os.listdir(out):
        fails.append(f"the refused package still wrote {os.listdir(out)}")

    # 2. The bytes themselves. With the declared sizes believed — a header can say
    #    anything — the count of what arrives is what stops it, and what it wrote is
    #    removed on the way out.
    out2 = os.path.join(folder, "bomb_out2")
    P.package_contents = lambda pkg: ["model.xlsx", "model_mesh.json"]
    try:
        xslope.unpack(bomb, dest=out2)
    except ValueError as exc:
        if "MB" not in str(exc):
            fails.append(f"the streamed refusal does not say how big it is: {exc}")
    else:
        fails.append("with the declared sizes believed, the bomb was extracted whole")
    finally:
        P.package_contents = _contents
    if os.path.isdir(out2) and os.listdir(out2):
        fails.append(f"the streamed refusal left {os.listdir(out2)} behind")

    # 3. A real project is nowhere near the ceiling — the guard is a bound, not a
    #    limit anyone meets.
    real = xslope.pack(os.path.join(_copy_project(WITH_SIDECARS),
                                    os.path.basename(WITH_SIDECARS)))
    if sum(i.file_size for i in zipfile.ZipFile(real).infolist()) > MAX_UNPACKED_BYTES / 8:
        fails.append("a real sample project is within a factor of eight of the "
                     "ceiling, which is set too low")
    _quiet(xslope.unpack, real, dest=os.path.join(folder, "real_out"))
    return fails


# ------------------------------------------------------------------------ G. Studio
def test_studio():
    """Open and Export, through the real window."""
    from PySide6.QtWidgets import QApplication, QDialog, QFileDialog, QMessageBox

    from studio.dialogs import UnpackPackageDialog
    from studio.main_window import MainWindow

    QApplication.instance() or QApplication([])
    fails = []
    folder = _copy_project(WITH_SIDECARS)
    book = os.path.join(folder, os.path.basename(WITH_SIDECARS))
    pkg = xslope.pack(book)

    # --- the destination dialog ------------------------------------------------
    dlg = UnpackPackageDialog(pkg)
    if dlg.dest.text() != os.path.splitext(pkg)[0]:
        fails.append(f"the dialog's default destination is {dlg.dest.text()}")
    if not dlg.btn_unpack.isVisibleTo(dlg) or dlg.btn_existing.isVisibleTo(dlg):
        fails.append("a free destination should offer Unpack and Open, and nothing else")
    taken = os.path.join(_tmp(), "taken")
    os.makedirs(taken)
    dlg.dest.setText(taken)
    if dlg.btn_unpack.isVisibleTo(dlg) or not dlg.btn_existing.isVisibleTo(dlg) \
            or not dlg.btn_fresh.isVisibleTo(dlg):
        fails.append("an existing destination must ask: open the existing copy, or "
                     "extract fresh")
    dlg._clicked(dlg.btn_fresh)
    fresh_dest, mode = dlg.chosen()
    if fresh_dest != taken + "-2" or mode != "extract":
        fails.append(f"Extract Fresh chose {fresh_dest} ({mode}), not {taken}-2")
    dlg.close()

    # --- File -> Open on a package --------------------------------------------
    mw = MainWindow()
    if not mw.act_export_pkg.text().startswith("Export"):
        fails.append("the File menu has no Export Project Package action")
    if mw.act_export_pkg.isEnabled():
        fails.append("Export Project Package is offered with no project open")

    recent = []
    mw._add_recent = recent.append       # never touch the user's settings from a test
    dest = os.path.join(_tmp(), "studio_open")
    _accept = UnpackPackageDialog.exec

    def _auto_exec(self):
        self.dest.setText(dest)
        self._mode = "extract"
        return QDialog.Accepted

    UnpackPackageDialog.exec = _auto_exec
    try:
        _quiet(mw.open_path, pkg)
    finally:
        UnpackPackageDialog.exec = _accept

    if not mw.doc.is_open:
        fails.append("opening a package did not open the project")
    elif os.path.normpath(mw.doc.path) != os.path.normpath(
            os.path.join(dest, os.path.basename(book))):
        fails.append(f"the project was opened from {mw.doc.path}, not the extracted "
                     f"workbook")
    if recent != [os.path.join(dest, os.path.basename(book))]:
        fails.append(f"recent files got {recent} — it must point at the loose "
                     f"workbook, never at the package")
    if os.path.basename(book) not in mw.windowTitle():
        fails.append(f"the window is titled {mw.windowTitle()!r} after opening a "
                     f"package")
    if mw.doc.slope_data.get("mesh") is None:
        fails.append("the project opened from a package carries no mesh")
    if not mw.act_export_pkg.isEnabled():
        fails.append("Export Project Package stayed disabled with a project open")

    # --- File -> Export Project Package ---------------------------------------
    # The save dialog is stubbed as an OBJECT, not as the static getSaveFileName:
    # the export builds a QFileDialog so it can set a default suffix, and a stub of
    # the static call would leave the real dialog to open and block forever.
    import studio.main_window as MW

    out = os.path.join(_tmp(), "exported" + PACKAGE_EXT)
    configured = {}

    class _StubSaveDialog:
        AcceptMode = QFileDialog.AcceptMode
        AcceptSave = QFileDialog.AcceptSave

        def __init__(self, parent=None, caption="", directory="", filter=""):
            configured["directory"] = directory
            configured["filter"] = filter

        def setAcceptMode(self, mode):
            configured["accept_mode"] = mode

        def setDefaultSuffix(self, suffix):
            configured["suffix"] = suffix

        def exec(self):
            return QDialog.Accepted

        def selectedFiles(self):
            return [configured.get("picked", out)]

    _dialog_cls = MW.QFileDialog
    _question = QMessageBox.question
    MW.QFileDialog = _StubSaveDialog
    QMessageBox.question = staticmethod(lambda *a, **k: QMessageBox.Cancel)
    try:
        _quiet(mw.export_package_dialog)
    finally:
        MW.QFileDialog = _dialog_cls
        QMessageBox.question = _question
    if not os.path.exists(out):
        fails.append("Export Project Package wrote nothing")
    else:
        got = sorted(package_contents(out))
        want = sorted(os.path.basename(p) for p in project_files(mw.doc.path))
        if got != want:
            fails.append(f"the exported package holds {got}, the project is {want}")
    # The extension is the DIALOG's default suffix, not something appended after it
    # returns: a name typed without one has to become "name.xslz" before the dialog
    # decides whether to ask about overwriting an existing file.
    if configured.get("suffix") != PACKAGE_EXT.lstrip("."):
        fails.append(f"the save dialog's default suffix is {configured.get('suffix')!r},"
                     f" so a typed name skips the overwrite confirmation")
    if configured.get("accept_mode") != QFileDialog.AcceptSave:
        fails.append("the export dialog is not in save mode, so it never asks about "
                     "overwriting")
    if not str(configured.get("directory", "")).endswith(PACKAGE_EXT):
        fails.append(f"the export dialog opens on {configured.get('directory')!r}, "
                     f"which does not name a package")

    # A project with a solution the session has not written out yet is saved first,
    # not packaged around: the question is asked, and Cancel stops the export.
    mw.doc.results["fem_solution"] = {"fem_data": None, "solution": None}
    os.remove(out)
    asked = []
    MW.QFileDialog = _StubSaveDialog
    QMessageBox.question = staticmethod(
        lambda *a, **k: (asked.append(a[1]), QMessageBox.Cancel)[1])
    try:
        _quiet(mw.export_package_dialog)
    finally:
        MW.QFileDialog = _dialog_cls
        QMessageBox.question = _question
    if not asked:
        fails.append("a run that has not been saved was packaged without a word")
    if os.path.exists(out):
        fails.append("cancelling the save prompt still exported a package")
    mw.doc.results.pop("fem_solution")

    mw.doc._dirty = False
    mw.close()
    return fails


CHECKS = [
    ("A. the file set is the right set", test_file_set),
    ("A2. attribution is the loaders' rule", test_attribution),
    ("A3. the writers spell what the set knows", test_writer_suffixes),
    ("B. what went in comes out", test_round_trip),
    ("C. a single-file project packages too", test_single_file),
    ("D. the collision is a refusal", test_collision),
    ("E. overwrite= and dest=", test_overwrite_and_dest),
    ("F. what is not a package", test_not_a_package),
    ("G. Studio opens and exports one", test_studio),
    ("H. what will not fit is refused", test_expansion),
]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("project packaging: PySide6 not installed — Studio check skipped.")
        checks = [c for c in CHECKS if c[1] is not test_studio]
    failures = []
    try:
        for name, fn in checks:
            try:
                fs = fn()
            except Exception as exc:
                import traceback
                traceback.print_exc()
                fs = [f"{name} raised: {exc!r}"]
            print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
            failures += fs
    finally:
        for d in _TMPDIRS:
            shutil.rmtree(d, ignore_errors=True)
        _TMPDIRS.clear()
    return failures


def main():
    print("project packaging (.xslz):")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll project-packaging checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication, QMessageBox
        _app = QApplication.instance() or QApplication([])
        QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
        QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
        QMessageBox.critical = staticmethod(lambda *a, **k: QMessageBox.Ok)
    except Exception:
        pass
    main()
