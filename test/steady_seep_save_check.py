"""Standing checks on two things a desktop session does between one run and the next:
attaching a seepage solution it just computed, and writing the file back out.

Both were reported from a first hands-on session with the frozen Windows build, and
neither is visible from a script — a script loads a file that already carries its
seepage sidecar, and it saves on a machine that happens to have a ``zip`` binary.

  A. THE STEADY FIELD IS PART OF THE MODEL. A stability run with ``u = seep`` reads
     ``slope_data['seep_u']``, and so does the gate that decides whether the run may
     start. Three things put a field there — the loader reading ``{base}_seep.csv``,
     the transient stagers reading one instant of a march, and a steady solve — and
     until now the third did not: it filled a results tab and a sidecar file, leaving
     the model itself carrying nothing. The gate then refused the run, correctly for
     what it could see, one second after the solution appeared on screen. A field
     belongs to the model the moment it is solved, not the moment it is reloaded.
  B. THE NO-SOLUTION REFUSAL STILL FIRES. Satisfying the gate from memory must not
     weaken it: a model with u = seep and no field of any kind is still refused, in
     the wording that rule is locked to, and a field is still no substitute for the
     mesh it would be read onto.
  C. THE FIELD IS THE ONE THAT WAS SOLVED. Applying it is refused when it was
     computed on a different mesh, it lands on the BC set it solved (1 -> seep_u,
     2 -> seep_u2, which is rapid drawdown's stage 2), and it clears any transient
     instant left behind — a steady field belongs to no moment of a march, and under
     automatic water loads that moment decides where the pool stands.
  D. SAVE NEEDS NOTHING BUT PYTHON. The .xlsx writer edits the workbook at the XML
     level and used to re-zip the archive by shelling out to ``zip``. There is no such
     executable on Windows: every Save there raised FileNotFoundError, surfaced as
     "[WinError 2] The system cannot find the file specified", while macOS and Linux
     worked because the tool ships with the OS. The check runs a real save with every
     subprocess entry point booby-trapped and PATH emptied, which is what a Windows
     machine looks like from inside this code.
  E. SAVE SURVIVES THE FROZEN LAYOUT. The template the writer copies is resolved from
     the xslope package directory, which in a PyInstaller bundle is a data directory
     inside the app. The check relocates the packaged template into a frozen-style
     tree and saves through it, so reader and writer are proven to agree on one
     resolution rather than on one that happens to work.
  F. THE ARCHIVE SURVIVES THE WRITE. Nothing is lost from the workbook, the stale
     calcChain is dropped with both of its references, fullCalcOnLoad is set, and a
     write that fails part way leaves the original file intact rather than truncated.
  G. NO HANDLE OUTLIVES A LOAD. Reading a workbook must not leave it open. POSIX
     hides a leaked handle completely -- an open file can still be replaced and
     unlinked -- so this is the half of the save story that only Windows can fail,
     and it fails it as "[WinError 32] ... being used by another process" on a save,
     a delete, or a user double-clicking their own file in Excel.

Skips its Studio leg cleanly when PySide6 is absent; A-F otherwise run either way.
"""
import contextlib
import copy
import io
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np

from xslope.fileio import (_round_cell_float, default_template_path,
                           load_slope_data, save_slope_data_to_xlsx,
                           write_cells_to_xlsx)
from xslope.preflight import preflight
from xslope.seep import apply_steady_stability_field

#: A limit-equilibrium model with a mesh, circles and u = seep throughout — the only
#: kind of model the missing-field gate has an opinion about. Its seepage BC set is
#: solvable in well under a second, so the Studio leg runs the REAL solver.
SEEP_LEM = os.path.join(_REPO, "docs/seep/files/xslope_rface_SEEP_KEY.xlsx")
#: A small, fast, seepage-free model for the save checks.
SIMPLE = os.path.join(_REPO, "docs/inputs/slope/xslope_simple1.xlsx")
#: The refusal wording ``seep_field.missing`` is locked to.
NO_FIELD = ("takes pore pressure from a seepage solution (u = seep), but this model "
            "carries no pore-pressure field")


def _quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def _seep_lem_model(with_mesh=True):
    """The sample with its stored field REMOVED — the state a user is in the moment a
    seepage solve finishes: a mesh, ``u = seep`` throughout, no ``{base}_seep.csv``."""
    sd = load_slope_data(SEEP_LEM)
    sd["seep_u"] = None
    if not with_mesh:
        sd["mesh"] = None
    return sd


def _steady_solution(n_nodes, value=100.0):
    """A solved-solution dict as ``run_seepage_analysis`` returns one. Only ``u`` is
    read by the staging path, so only ``u`` is arranged."""
    return {"u": np.full(n_nodes, float(value))}


# ------------------------------------------------ A + B. the gate, both directions
def test_gate():
    """A solved steady field satisfies the run gate; nothing else about it moves."""
    fails = []
    sel = {"surface": "circular"}

    # 1. No field of any kind: the ERROR stands, in the words it is locked to.
    sd = _seep_lem_model()
    errs = [f for f in preflight(sd, "lem", sel).errors
            if f.rule_id == "seep_field.missing"]
    if not errs:
        fails.append("a model with u = seep and no field of any kind was let through")
    elif NO_FIELD not in errs[0].message:
        fails.append("the no-solution refusal has been reworded: " + errs[0].message)

    # 2. The same model one steady solve later. This is the reported bug: the run was
    #    refused for a solution that had just been computed and was on screen.
    _quiet(apply_steady_stability_field, sd, _steady_solution(len(sd["mesh"]["nodes"])))
    report = preflight(sd, "lem", sel)
    if any(f.rule_id == "seep_field.missing" for f in report.errors):
        fails.append("a model carrying a freshly solved steady field was still "
                     "refused for a missing pore-pressure field")
    if any(f.rule_id == "seep_field.node_count_mismatch" for f in report.errors):
        fails.append("the applied field was read as belonging to another mesh")
    # A steady field names no instant, so nothing about a march may be claimed.
    if any(f.rule_id == "seep_field.transient_frame" for f in report.findings):
        fails.append("a steady field was reported as an instant of a transient march")

    # 3. A field is not a substitute for the MESH it would be read onto.
    sd = _seep_lem_model(with_mesh=False)
    sd["seep_u"] = np.full(10, 100.0)
    errs = [f for f in preflight(sd, "lem", sel).errors
            if f.rule_id == "seep_field.missing"]
    if not errs:
        fails.append("a stored field was allowed to stand in for a missing mesh")
    elif "carries no mesh" not in errs[0].message:
        fails.append("the missing-mesh refusal has been reworded: " + errs[0].message)
    return fails


# ------------------------------------------------------- C. what applying it means
def test_apply():
    """The staging entry point: which key, which mesh, and what it clears."""
    fails = []
    sd = _seep_lem_model()
    n = len(sd["mesh"]["nodes"])

    # BC set 1 is the model's pore-pressure field; set 2 is rapid drawdown's stage 2.
    _quiet(apply_steady_stability_field, sd, _steady_solution(n, 7.0), bc=1)
    _quiet(apply_steady_stability_field, sd, _steady_solution(n, 9.0), bc=2)
    if sd.get("seep_u") is None or float(np.asarray(sd["seep_u"])[0]) != 7.0:
        fails.append(f"BC set 1 did not land in seep_u: {sd.get('seep_u')}")
    if sd.get("seep_u2") is None or float(np.asarray(sd["seep_u2"])[0]) != 9.0:
        fails.append(f"BC set 2 did not land in seep_u2: {sd.get('seep_u2')}")

    # A steady field belongs to no instant. A stale seep_u_time would date these
    # pressures to a moment of a march they never came from — and under automatic
    # water loads that moment is where the pool is read.
    sd["seep_u_time"] = 47.0
    _quiet(apply_steady_stability_field, sd, _steady_solution(n, 7.0))
    if "seep_u_time" in sd:
        fails.append("a steady field kept a transient instant: "
                     f"seep_u_time = {sd['seep_u_time']}")

    # A field solved on another mesh is refused, and the model is left alone.
    before = np.asarray(sd["seep_u"]).copy()
    try:
        _quiet(apply_steady_stability_field, sd, _steady_solution(n + 3, 99.0))
        fails.append("a field with the wrong node count was accepted")
    except ValueError as exc:
        if "different mesh" not in str(exc):
            fails.append("the wrong-mesh refusal does not say so: " + str(exc))
    if not np.array_equal(np.asarray(sd["seep_u"]), before):
        fails.append("a refused field replaced the one already in the model")

    # No solution at all is a caller error, not a silent no-op.
    try:
        _quiet(apply_steady_stability_field, sd, None)
        fails.append("apply_steady_stability_field(None) was accepted")
    except ValueError:
        pass
    try:
        _quiet(apply_steady_stability_field, sd, _steady_solution(n), bc=3)
        fails.append("an unknown BC set was accepted")
    except ValueError:
        pass
    return fails


# ----------------------------------------------------- D. save needs nothing but Python
class _NoSubprocess:
    """Every way out of this process, closed — which is what Windows looks like to
    code that shells out to a POSIX tool. ``subprocess`` raises the exception Windows
    raises for a missing executable, and PATH is emptied so a lookup finds nothing."""

    _NAMES = ("run", "call", "check_call", "check_output", "Popen", "getoutput",
              "getstatusoutput")

    def __enter__(self):
        def boom(*a, **k):
            raise FileNotFoundError(2, "The system cannot find the file specified")

        self._saved = {n: getattr(subprocess, n) for n in self._NAMES
                       if hasattr(subprocess, n)}
        for n in self._saved:
            setattr(subprocess, n, boom)
        self._path = os.environ.get("PATH")
        os.environ["PATH"] = ""
        self._popen = os.popen
        self._system = os.system
        os.popen = boom
        os.system = boom
        return self

    def __exit__(self, *exc):
        for n, fn in self._saved.items():
            setattr(subprocess, n, fn)
        os.environ["PATH"] = self._path if self._path is not None else ""
        os.popen = self._popen
        os.system = self._system
        return False


def test_save_without_a_shell():
    """A full save with no external process available anywhere.

    This is bug (b) exactly: the writer re-zipped the workbook by running ``zip``,
    which does not exist on Windows, so Save raised FileNotFoundError there and
    Studio showed "[WinError 2] The system cannot find the file specified". Nothing
    on the save path may depend on a binary that is not part of what we ship.
    """
    fails = []
    sd = load_slope_data(SIMPLE)
    tmpdir = tempfile.mkdtemp(prefix="steady_seep_save_")
    try:
        out = os.path.join(tmpdir, "saved.xlsx")
        try:
            with _NoSubprocess():
                _quiet(save_slope_data_to_xlsx, sd, out)
        except FileNotFoundError as exc:
            fails.append(f"save shelled out to a missing executable: {exc!r}")
            return fails
        back = load_slope_data(out)
        if len(back["materials"]) != len(sd["materials"]):
            fails.append("the saved file lost materials")
        if abs(float(back["gamma_water"]) - float(sd["gamma_water"])) > 1e-9:
            fails.append("the saved file did not carry gamma_water")
        if not (back.get("circles") or back.get("non_circ")):
            fails.append("the saved file carries no failure surface")
        # In place, onto a file that already exists — the owner's actual action.
        try:
            with _NoSubprocess():
                _quiet(save_slope_data_to_xlsx, back, out)
        except FileNotFoundError as exc:
            fails.append(f"re-save shelled out to a missing executable: {exc!r}")
        # And the low-level writer on its own, on a workbook that HAS a calcChain.
        legacy = os.path.join(_REPO, "docs/inputs/input_template_v15.xlsx")
        if os.path.exists(legacy):
            cc = os.path.join(tmpdir, "cc.xlsx")
            shutil.copy(legacy, cc)
            try:
                with _NoSubprocess():
                    _quiet(write_cells_to_xlsx, cc, {"main": {"D10": 62.4}})
            except FileNotFoundError as exc:
                fails.append(f"write_cells_to_xlsx shelled out: {exc!r}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


# ------------------------------------------------------ E. the frozen resource layout
def test_frozen_template():
    """Save through the packaged template as a frozen bundle presents it.

    PyInstaller places ``xslope/resources`` as a DATA directory beside the frozen
    modules, so ``default_template_path`` resolves inside the app rather than inside
    site-packages. The reader (Studio's Save As, ``default_template_path``) and the
    writer must agree on that one resolution: the build's smoke test proves the file
    is READABLE there, which is not the same as proving a save works from it.
    """
    import xslope.fileio as fio

    fails = []
    tmpdir = tempfile.mkdtemp(prefix="steady_seep_frozen_")
    try:
        # A frozen-style tree: <_MEIPASS>/xslope/resources/input_template.xlsx, with
        # no repository anywhere above it.
        res = os.path.join(tmpdir, "_internal", "xslope", "resources")
        os.makedirs(res)
        shutil.copy(default_template_path(),
                    os.path.join(res, "input_template.xlsx"))
        frozen = os.path.join(res, "input_template.xlsx")

        saved = fio.default_template_path
        fio.default_template_path = lambda: frozen
        try:
            if fio.default_template_path() != frozen:
                fails.append("the frozen resolution did not take effect")
            sd = load_slope_data(SIMPLE)
            out = os.path.join(tmpdir, "frozen_save.xlsx")
            with _NoSubprocess():
                _quiet(save_slope_data_to_xlsx, sd, out)
            back = load_slope_data(out)
            if len(back["materials"]) != len(sd["materials"]):
                fails.append("the frozen-layout save lost materials")
        finally:
            fio.default_template_path = saved

        # The template it copied is untouched: a save must never write into the
        # bundle's own read-only resource.
        if os.path.getsize(frozen) != os.path.getsize(default_template_path()):
            fails.append("the save wrote into the packaged template itself")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


# ------------------------------------------------------------ F. the archive itself
def test_archive():
    """What the writer leaves behind: every part, the calcChain gone, and no debris."""
    fails = []
    tmpdir = tempfile.mkdtemp(prefix="steady_seep_zip_")
    try:
        legacy = os.path.join(_REPO, "docs/inputs/input_template_v15.xlsx")
        if not os.path.exists(legacy):
            return fails
        dst = os.path.join(tmpdir, "cc.xlsx")
        shutil.copy(legacy, dst)
        _quiet(write_cells_to_xlsx, dst, {"main": {"D10": 62.4}})
        with zipfile.ZipFile(legacy) as a, zipfile.ZipFile(dst) as b:
            na = [i.filename for i in a.infolist()]
            nb = [i.filename for i in b.infolist()]
            lost = set(na) - set(nb) - {"xl/calcChain.xml"}
            if lost:
                fails.append(f"the save dropped workbook parts: {sorted(lost)}")
            if set(nb) - set(na):
                fails.append(f"the save invented parts: {sorted(set(nb) - set(na))}")
            if "xl/calcChain.xml" in nb:
                fails.append("the stale calcChain survived the save")
            ct = b.read("[Content_Types].xml").decode("utf-8")
            rels = b.read("xl/_rels/workbook.xml.rels").decode("utf-8")
            if "calcChain" in ct or "calcChain" in rels:
                fails.append("calcChain was deleted but still referenced — Excel "
                             "will 'recover' the file and discard the edits")
            if 'fullCalcOnLoad="1"' not in b.read("xl/workbook.xml").decode("utf-8"):
                fails.append("fullCalcOnLoad was not set, so formulas stay stale")
            if b.testzip() is not None:
                fails.append(f"the written archive is corrupt at {b.testzip()}")
            # Untouched parts come through byte-for-byte, not re-encoded.
            untouched = [n for n in nb if n.startswith("xl/media/")
                         or n.endswith(".rels") and n != "xl/_rels/workbook.xml.rels"]
            for n in untouched:
                if a.read(n) != b.read(n):
                    fails.append(f"an untouched part was rewritten: {n}")
                    break
        # No temp file left beside the destination (the atomic-replace scratch).
        debris = [f for f in os.listdir(tmpdir) if f != "cc.xlsx"]
        if debris:
            fails.append(f"the save left debris beside the file: {debris}")

        # A write that fails part way leaves the ORIGINAL intact, not a truncated file.
        before = open(dst, "rb").read()
        import xslope.fileio as fio
        boom = fio._reset_view

        def explode(_xml):
            raise RuntimeError("simulated failure mid-write")

        fio._reset_view = explode
        try:
            _quiet(write_cells_to_xlsx, dst, {"main": {"D10": 1.0}})
            fails.append("a failing write reported success")
        except RuntimeError:
            pass
        finally:
            fio._reset_view = boom
        if open(dst, "rb").read() != before:
            fails.append("a failed write damaged the original file")
        debris = [f for f in os.listdir(tmpdir) if f != "cc.xlsx"]
        if debris:
            fails.append(f"a failed write left debris: {debris}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


# ------------------------------------------------------- G. no handle outlives a load
def _open_handles(substring):
    """Paths this process currently holds open that contain ``substring``.

    Returns ``None`` when the platform offers no way to ask, so a caller can say so
    rather than report a clean result it did not measure. Three ways are tried,
    because the one that works differs by platform and the check must not quietly
    become a no-op on the machine that matters:

      * ``psutil.Process().open_files()`` -- everywhere, when it is installed;
      * ``/proc/self/fd`` -- Linux, where the entries are readable symlinks;
      * ``lsof -p`` -- macOS, where ``/dev/fd`` entries are NOT resolvable symlinks
        and a readlink scan silently finds nothing at all. That false clean is how
        this leak survived a mac-side check in the first place.
    """
    try:
        import psutil
        return [f.path for f in psutil.Process().open_files() if substring in f.path]
    except Exception:
        pass
    if os.path.isdir("/proc/self/fd"):
        out = []
        for fd in os.listdir("/proc/self/fd"):
            try:
                target = os.readlink("/proc/self/fd/" + fd)
            except OSError:
                continue
            if substring in target:
                out.append(target)
        return out
    if shutil.which("lsof"):
        try:
            res = subprocess.run(["lsof", "-p", str(os.getpid())],
                                 capture_output=True, text=True, timeout=60)
        except Exception:
            return None
        return [line.split()[-1] for line in res.stdout.splitlines()
                if substring in line]
    return None


def test_no_leaked_handle():
    """Nothing may still hold the workbook open once it has been read.

    ``pd.ExcelFile(path)`` keeps the .xlsx open for as long as it lives and defines
    no ``__del__``, so the handle went back to the operating system only when the
    object graph behind it happened to be collected -- and ``load_slope_data`` both
    outlives that call and can raise on a dozen validation paths that keep the frame,
    and the parser, alive with the exception. A handle on the input file therefore
    outlived the call.

    POSIX hides this completely: an open file can still be renamed, replaced and
    unlinked, so the mac build and every mac test passed. Windows cannot -- the frozen
    build's self-test died in temp-directory cleanup with "[WinError 32] the process
    cannot access the file because it is being used by another process". For a user
    that is worse than a failed test: the file they just opened cannot be saved over,
    cannot be deleted, and cannot be opened in Excel while Studio is running.

    Both halves of the platform are asserted, so this is caught wherever it runs: the
    open-handle scan is the signal on macOS and Linux, and deleting the file
    immediately after loading it is the signal on Windows, where the delete is what
    actually fails.
    """
    fails = []
    tmpdir = tempfile.mkdtemp(prefix="steady_seep_handle_")
    try:
        probe = os.path.join(tmpdir, "handle_probe.xlsx")
        shutil.copy(SIMPLE, probe)

        sd = load_slope_data(probe)
        held = _open_handles("handle_probe.xlsx")
        if held is None:
            print("    (no open-handle detector on this platform; "
                  "the delete half still applies)")
        elif held:
            fails.append(f"load_slope_data left the workbook open: {held}")

        # The whole save path, the sequence the frozen self-test walks.
        out = os.path.join(tmpdir, "handle_saved.xlsx")
        save_slope_data_to_xlsx(sd, out)
        reread = load_slope_data(out)
        save_slope_data_to_xlsx(reread, out)
        load_slope_data(out)
        held = _open_handles("handle_saved.xlsx")
        if held:
            fails.append(f"save/load/save/load left handles open: {held}")

        # What Windows actually refuses. On POSIX this always succeeds, so it is the
        # delete on the OTHER platform that this line is standing guard over -- and
        # it is the exact operation TemporaryDirectory cleanup performs.
        for path in (probe, out):
            try:
                os.remove(path)
            except OSError as exc:
                fails.append(f"could not delete {os.path.basename(path)} right after "
                             f"loading it -- a handle is still open: {exc}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


# ------------------------------------------------------------- the Studio leg (A, live)
def test_studio():
    """The reported session, offscreen: open a seepage model with no saved solution,
    watch Run LEM refuse it, solve the seepage through the real runner, and watch the
    same dialog let the run start — with the field in the model, not only on screen.
    """
    from studio.main_window import MainWindow
    from studio.runners import SeepRunner
    from studio.dialogs import RunLemDialog
    from xslope.slice import generate_slices

    fails = []
    tmpdir = tempfile.mkdtemp(prefix="steady_seep_studio_")
    mw = None
    try:
        # The model as a user has it before the first solve: workbook + mesh, and
        # NO {base}_seep.csv. Nothing on disk can supply a pore-pressure field.
        stem = os.path.splitext(SEEP_LEM)[0]
        xlsx = os.path.join(tmpdir, "seeprun.xlsx")
        shutil.copy(SEEP_LEM, xlsx)
        shutil.copy(stem + "_mesh.json", os.path.join(tmpdir, "seeprun_mesh.json"))

        mw = MainWindow()
        _quiet(mw.open_path, xlsx)
        sd = mw.doc.slope_data
        if sd.get("seep_u") is not None:
            fails.append("the fixture already carries a seepage field — the check "
                         "would prove nothing")
            return fails
        if sd.get("mesh") is None:
            fails.append("the fixture lost its mesh on open")
            return fails

        # Before the solve the refusal is CORRECT, and stays the contract's wording.
        dlg = RunLemDialog(None, slope_data=sd)
        if not dlg.preflight.blocked:
            fails.append("Run LEM was offered on a model with no seepage solution")
        elif NO_FIELD not in dlg.preflight.block_reason():
            fails.append("the pre-solve refusal is not the missing-field one: "
                         + dlg.preflight.block_reason())
        if dlg._ok.isEnabled():
            fails.append("Run was pressable with no seepage solution")
        dlg.deleteLater()

        # The solve, through the real runner on the real solver.
        runner = SeepRunner(sd, {"bc": 1, "tol": 1e-4})
        got, err = {}, {}
        runner.succeeded.connect(lambda b: got.update(b))
        runner.failed.connect(lambda m: err.setdefault("msg", m))
        _quiet(runner.run)
        if err:
            fails.append(f"the steady runner failed: {err['msg']}")
            return fails
        if "solution" not in got:
            fails.append("the steady runner emitted no solution bundle")
            return fails
        _quiet(mw._on_seep_succeeded, got)

        # The reported bug: everything above happened and the model still carried
        # nothing, so the gate refused the run a second time.
        field = mw.doc.slope_data.get("seep_u")
        if field is None:
            fails.append("the solved steady field never reached the model "
                         "(slope_data['seep_u'] is still None)")
            return fails
        if len(field) != len(sd["mesh"]["nodes"]):
            fails.append(f"the attached field has {len(field)} values for "
                         f"{len(sd['mesh']['nodes'])} mesh nodes")

        dlg = RunLemDialog(None, slope_data=mw.doc.slope_data)
        if dlg.preflight.blocked:
            fails.append("Run LEM was refused after a seepage solve: "
                         + dlg.preflight.block_reason())
        if any(f.rule_id == "seep_field.missing"
               for f in dlg.preflight.report.errors):
            fails.append("the missing-field error still fires with a solved field "
                         "in the model")
        if not dlg._ok.isEnabled():
            fails.append("Run stayed disabled after a seepage solve: "
                         + dlg._ok.toolTip())
        dlg.deleteLater()

        # ...and the run it would start reads those pore pressures, not zero.
        run_sd = mw.doc.slope_data
        ok, res = _quiet(generate_slices, run_sd, circle=run_sd["circles"][0],
                         num_slices=20)
        if not ok:
            fails.append(f"the run did not get past the slicer: {res}")
        elif float(res[0]["u"].max()) <= 0.0:
            fails.append("every slice base read zero pore pressure after the solve")

        # A steady solve does not leave the document looking edited: like the mesh,
        # the field is a computed artifact with its own sidecar, not a user edit.
        if mw.doc.dirty:
            fails.append("a seepage solve marked the document dirty")
    finally:
        if mw is not None:
            mw.doc._dirty = False
            mw.close()
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


def test_small_magnitudes_survive_the_write():
    """A cell keeps ten significant digits however small the number is.

    The writer strips binary-repr noise on the way into a cell. Measuring that
    tail from the decimal point instead of from the number destroys real digits
    as the value shrinks, and seepage is where the small numbers live: a 2:1 rain
    boundary carries q_n = 1e-8 * cos(atan 0.5) = 8.94427191e-9 m/s, which a
    ten-DECIMAL-place rule writes as 8.9e-9 — two significant figures, a 0.6%
    error in the boundary, silently. Below 5e-11 the same rule writes 0, so a
    1e-13 m/s conductivity loads back as a hole in the model.

    Every value here is written into a real workbook through ``save_slope_data_to_xlsx``
    and read back through ``load_slope_data``, so the cell builders are covered along
    with the rounding helper: a value can be lost in either, and only the round trip
    sees both. The helper is checked first so a failure says which half broke.
    """
    fails = []
    # (label, value) — the rain-boundary normal flux, a conductivity below the old
    # floor, one either side of it, and the noise the rounding exists to remove.
    NOISE = 0.1 + 0.2
    cases = [("2:1 rain boundary q_n", 1e-8 * 2.0 / 5.0 ** 0.5),
             ("low-k lens conductivity", 1e-13),
             ("just under the old floor", 4.9e-11),
             ("nine significant digits", 1.234567891e-9),
             ("negative small", -8.94427191e-9),
             ("everyday magnitude", 0.3),
             ("repr noise", NOISE)]

    def _judge(where, label, v, got):
        if v == NOISE:                               # the one value meant to move
            if got != 0.3:
                fails.append(f"{where}: {label}: repr noise survived, {got!r}")
            return
        if got == 0.0:
            fails.append(f"{where}: {label}: {v!r} came back as zero")
            return
        rel = abs(got - v) / abs(v)
        if rel > 5e-10:                              # ten significant digits
            fails.append(f"{where}: {label}: {v!r} became {got!r}, relative error "
                         f"{rel:.2g} — significant digits lost")

    tmpdir = tempfile.mkdtemp(prefix="steady_seep_round_")
    try:
        for label, v in cases:
            _judge("helper", label, v, _round_cell_float(v))

        # The real path. gw006d is the donor because it carries flux boundary blocks,
        # which is where a rate of this magnitude actually lives; each value is written
        # into all of them, saved, and read back. A save writes some cells fresh and
        # edits others in place, so both cell builders are on this path.
        src = os.path.join(_REPO, "docs/verification/files/rocscience_gw/gw006d.xlsx")
        if not os.path.exists(src):
            fails.append("gw006d.xlsx is missing — the round trip did not run")
        else:
            donor = _quiet(load_slope_data, src)
            if not donor["seepage_bc"].get("specified_fluxes"):
                fails.append("the donor carries no flux block to write through")
            out = os.path.join(tmpdir, "roundtrip.xlsx")
            for label, v in cases:
                sd = copy.deepcopy(donor)
                for b in sd["seepage_bc"]["specified_fluxes"]:
                    b["flux"] = v
                _quiet(save_slope_data_to_xlsx, sd, out)
                back = _quiet(load_slope_data, out)
                got = [b["flux"] for b in back["seepage_bc"]["specified_fluxes"]]
                if len(got) != len(sd["seepage_bc"]["specified_fluxes"]):
                    fails.append(f"round trip: {label}: {len(got)} flux blocks came "
                                 f"back, {len(sd['seepage_bc']['specified_fluxes'])} "
                                 "went in")
                for r in got:
                    _judge("round trip", label, v, r)

            # ...and the shipped artifact still carries its exact rate.
            want = 1e-8 * 2.0 / 5.0 ** 0.5
            face = [b["flux"] for b in donor["seepage_bc"]["specified_fluxes"]
                    if b["flux"] < 9.9e-9]
            if not face:
                fails.append("gw006d carries no sloped-face flux block to check")
            for r in face:
                if abs(r - want) / want > 5e-10:
                    fails.append(f"gw006d's 2:1 face flux reads {r!r}, not {want!r} "
                                 "— the writer truncated it")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


CHECKS = [("the run gate, both directions", test_gate),
          ("what applying a steady field means", test_apply),
          ("save with no external process available", test_save_without_a_shell),
          ("save through the frozen resource layout", test_frozen_template),
          ("the written archive", test_archive),
          ("small magnitudes survive the write", test_small_magnitudes_survive_the_write),
          ("no handle outlives a load", test_no_leaked_handle),
          ("Studio: solve then run, offscreen", test_studio)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("steady seep + save: PySide6 not installed — Studio check skipped.")
        checks = [c for c in CHECKS if c[1] is not test_studio]
    failures = []
    for name, fn in checks:
        try:
            fs = fn()
        except Exception as exc:
            import traceback
            traceback.print_exc()
            fs = [f"{name} raised: {exc!r}"]
        print(f"  {name:44s} {'ok' if not fs else f'FAIL ({len(fs)})'}")
        failures += fs
    return failures


def main():
    print("steady seepage field + save path:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll steady-field and save-path checks passed.")


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
