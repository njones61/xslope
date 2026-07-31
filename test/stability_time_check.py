"""Standing checks on the transient stability time — WHICH instant of a transient
seepage march a stability run reads, and whether that choice survives the file.

A stability run against a time-dependent seepage field reports a factor of safety for
ONE instant. The instant is chosen in four places that must agree, and each of them is
a behaviour a wording change or a re-wired signal could break in silence:

  A. RESOLUTION — the precedence a run entry point applies: an explicit ``time=``
     beats the model's ``tseep`` ``stability_time``, which beats the engine default.
     The engine default is the LAST saved frame, and that is also what a blank
     ``stability_time`` means; the resolver has to name which of the three it used,
     because the run's log line is built from it.
  B. THE FILE — ``stability_time`` reads off the v22 tseep sheet, round-trips through
     the writer, and is ABSENT (never a stray number in the save_times column) on a
     pre-v22 destination. A model that never declared one comes back declaring none,
     which is the blank-vs-zero contract applied to a time.
  C. THE SCHEDULE — the march has to LAND on the instants a stability run will ask
     for, so ``stability_time`` joins ``stage_1``/``stage_2`` in the saved-frame
     schedule, and a re-march injects a requested time into the save times rather than
     interpolating a field between two frames.
  D. THE DIALOGS — the Run LEM / Run FEM seepage-time selector: it opens on the
     model's own answer, says in words which instant will be read (and that a blank
     stability time means the last frame), marks a free-entry time as needing a
     re-march, refuses a time outside the run, and hands back the stage times the
     rapid-drawdown path edits at its point of use.
  E. THE RUN GATE — the staging failure that used to be silent. Stage times with no
     saved frames are a real, reachable mistake; the run must be REFUSED and the model
     left exactly as it was, because the alternative is a factor of safety computed
     from the previous run's pore pressures with nothing on screen to say so.

No seepage is solved here: the solver has its own checks, and every behaviour above is
about selection, not about the field. The transient solutions are synthetic frame
ledgers, which is all the selection path reads.

Skips cleanly (exit 0) when PySide6 is not installed; the non-Studio checks (A-C) run
either way.
"""
import os
import shutil
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np

from xslope.fileio import load_slope_data, save_slope_data_to_xlsx
from xslope.seep import (apply_transient_stability_frame, has_transient_frame,
                         resolve_stability_time, select_transient_frame_u,
                         stage_transient_for_drawdown, _transient_saved_schedule)

#: A tseep-bearing model (duration 360, stages 0/47, save_times 2/15/30/47/80).
DAM = os.path.join(_REPO, "docs/seep/files/xslope_earth_dam_tseep.xlsx")
MASTER = os.path.join(_REPO, "docs/inputs/input_template.xlsx")


def _solution(times, n_nodes=6):
    """A synthetic transient solution: the frame ledger the selection path reads.

    Each frame's pore pressures are the frame's own time, so a test can tell WHICH
    frame was selected by looking at the field alone."""
    return {"times": [float(t) for t in times],
            "frames": [{"time": float(t), "u": np.full(n_nodes, float(t))}
                       for t in times]}


def _u_time(slope_data):
    """The time the pore pressures in ``slope_data`` came from, read off the field
    itself rather than off the bookkeeping key."""
    return float(np.asarray(slope_data["seep_u"])[0])


# ---------------------------------------------------------------- A. resolution
def test_resolution():
    fails = []
    sol = _solution([2.0, 15.0, 30.0, 47.0, 360.0])

    # No stability_time anywhere: the engine default, the LAST saved frame.
    sd = {"tseep": {"duration": 360.0}}
    t, src = resolve_stability_time(sd, sol)
    if (t, src) != (360.0, "default"):
        fails.append(f"blank stability_time resolved to ({t}, {src}), want (360.0, 'default')")

    # The file's value, which is what makes a scripted re-run reproduce a choice.
    sd = {"tseep": {"duration": 360.0, "stability_time": 30.0}}
    t, src = resolve_stability_time(sd, sol)
    if (t, src) != (30.0, "file"):
        fails.append(f"file stability_time resolved to ({t}, {src}), want (30.0, 'file')")

    # An explicit time overrides the file — the dialog's choice, for that run only.
    t, src = resolve_stability_time(sd, sol, time=15.0)
    if (t, src) != (15.0, "argument"):
        fails.append(f"explicit time resolved to ({t}, {src}), want (15.0, 'argument')")

    # select_transient_frame_u obeys the same order without being told.
    sd = {"tseep": {"duration": 360.0, "stability_time": 30.0}}
    select_transient_frame_u(sd, sol)
    if _u_time(sd) != 30.0:
        fails.append(f"select_transient_frame_u(time=None) read t={_u_time(sd)}, "
                     f"want the file's 30.0")
    sd = {"tseep": {"duration": 360.0}}
    select_transient_frame_u(sd, sol)
    if _u_time(sd) != 360.0:
        fails.append(f"select_transient_frame_u with a blank stability_time read "
                     f"t={_u_time(sd)}, want the last frame 360.0")
    # A model with no tseep sheet at all must not trip the resolver.
    sd = {}
    select_transient_frame_u(sd, sol)
    if _u_time(sd) != 360.0:
        fails.append("select_transient_frame_u failed on a model with no tseep sheet")

    # The entry surface: same resolution, plus a record of what it did.
    sd = {"tseep": {"duration": 360.0, "stability_time": 15.0}}
    info = apply_transient_stability_frame(sd, sol, verbose=False)
    if info["times"] != [15.0] or info["source"] != "file" or info["remarched"]:
        fails.append(f"apply_transient_stability_frame reported {info!r}")
    if _u_time(sd) != 15.0 or sd["seep_u_time"] != 15.0:
        fails.append("apply_transient_stability_frame read the wrong frame")

    # Rapid drawdown reads its two stage times and ignores any single-run time.
    sd = {"tseep": {"duration": 360.0, "stage_1": 15.0, "stage_2": 47.0,
                    "stability_time": 30.0}}
    info = apply_transient_stability_frame(sd, sol, rapid=True, verbose=False)
    if _u_time(sd) != 15.0 or float(np.asarray(sd["seep_u2"])[0]) != 47.0:
        fails.append("rapid staging did not read the stage_1/stage_2 frames")
    if info["source"] != "stages":
        fails.append(f"rapid staging reported source={info['source']!r}")

    # A time with no saved frame is refused by name, not silently snapped.
    sd = {"tseep": {"duration": 360.0}}
    try:
        apply_transient_stability_frame(sd, sol, time=100.0, verbose=False)
        fails.append("a time with no saved frame was accepted")
    except ValueError as e:
        if "100" not in str(e):
            fails.append(f"the refusal does not name the requested time: {e}")
    if has_transient_frame(sol, 100.0):
        fails.append("has_transient_frame said 100.0 is a saved frame")
    if not has_transient_frame(sol, 47.0):
        fails.append("has_transient_frame did not recognise a saved frame")
    return fails


# ------------------------------------------------------------------- B. the file
def test_file():
    fails = []
    sd = load_slope_data(DAM)
    ts = sd.get("tseep") or {}
    if "stability_time" not in ts:
        fails.append("the parsed tseep dict carries no 'stability_time' key")
    if ts.get("stability_time") is not None:
        fails.append(f"a model that declares no stability time read "
                     f"{ts['stability_time']!r} instead of None")

    tmpdir = tempfile.mkdtemp(prefix="stability_time_")
    try:
        # Round-trip a declared value through the v22 master.
        sd["tseep"] = {**ts, "stability_time": 30.0}
        out = os.path.join(tmpdir, "declared.xlsx")
        save_slope_data_to_xlsx(sd, out, template=MASTER)
        back = load_slope_data(out)["tseep"]
        if back.get("stability_time") != 30.0:
            fails.append(f"stability_time round-tripped as {back.get('stability_time')!r}, "
                         f"want 30.0")
        # The cell it lands in must not disturb its neighbours: the save-times list
        # sits directly under it and would read short if the row were off by one.
        if back.get("save_times") != ts.get("save_times"):
            fails.append(f"save_times changed across the round-trip: "
                         f"{ts.get('save_times')} -> {back.get('save_times')}")
        for key in ("duration", "save_interval", "stage_1", "stage_2"):
            if back.get(key) != ts.get(key):
                fails.append(f"tseep {key} changed across the round-trip: "
                             f"{ts.get(key)!r} -> {back.get(key)!r}")

        # A blank stays blank — never a stamped zero.
        sd["tseep"] = {**ts, "stability_time": None}
        out = os.path.join(tmpdir, "blank.xlsx")
        save_slope_data_to_xlsx(sd, out, template=MASTER)
        back = load_slope_data(out)["tseep"]
        if back.get("stability_time") is not None:
            fails.append(f"a blank stability_time came back as "
                         f"{back['stability_time']!r}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


# --------------------------------------------------------------- C. the schedule
def test_schedule():
    fails = []
    # The march has to land on the instant a stability run will read, or the run
    # cannot have it at all.
    tseep_data = {"duration": 100.0, "save_interval": 50.0, "save_times": [],
                  "stage_1": 10.0, "stage_2": 40.0, "stability_time": 33.0,
                  "breakpoints": []}
    sched, _ = _transient_saved_schedule(tseep_data, 100.0)
    for t in (10.0, 33.0, 40.0, 100.0):
        if t not in sched:
            fails.append(f"the saved-frame schedule omits t={t}: {sched}")

    # An unset stability time adds nothing.
    sched2, _ = _transient_saved_schedule({**tseep_data, "stability_time": None}, 100.0)
    if 33.0 in sched2:
        fails.append("a blank stability_time still put t=33 on the schedule")

    # A re-march is an INJECTION into the save schedule, not a blend of two frames:
    # the requested instant becomes a computed state.
    sched3, _ = _transient_saved_schedule(
        {**tseep_data, "stability_time": None, "save_times": [77.0]}, 100.0)
    if 77.0 not in sched3:
        fails.append(f"an injected save time is missing from the schedule: {sched3}")

    # remarch_for_times reaches the solver with the requested instants on the save
    # list and the model's own save times intact. The march itself is the solver's
    # business (and has its own checks), so it is stubbed out here.
    import xslope.seep as seep
    captured = {}
    _orig = seep.run_transient_seepage
    seep.run_transient_seepage = lambda sd, td, **kw: captured.setdefault("td", td)
    try:
        sd = load_slope_data(DAM)
        seep.remarch_for_times(None, sd, [100.0, 0.0, 1e9, None])
    finally:
        seep.run_transient_seepage = _orig
    got = sorted(captured.get("td", {}).get("save_times", []))
    if 100.0 not in got:
        fails.append(f"a re-march did not inject the requested time: {got}")
    for t in (sd["tseep"]["save_times"] or []):
        if t not in got:
            fails.append(f"a re-march dropped the model's own save time {t}: {got}")
    for t in (0.0, 1e9):
        if t in got:
            fails.append(f"a re-march kept t={t}, which no march can land on: {got}")
    if sd["tseep"]["save_times"] == got:
        fails.append("the re-march mutated the model's own save_times")
    return fails


# ----------------------------------------------------------------- D+E. Studio
def test_studio():
    fails = []
    from PySide6.QtWidgets import QMessageBox
    from studio.dialogs import RunFemDialog, RunLemDialog
    from studio.main_window import MainWindow

    sd = load_slope_data(DAM)
    sol = _solution([0.0, 2.0, 15.0, 30.0, 47.0, 80.0, 360.0])

    # --- the selector opens on the model's answer and says which one it is ---
    dlg = RunLemDialog(None, slope_data=sd, transient=sol, current_time=15.0)
    opts = dlg.options()
    if opts["seep_time"] is None:
        fails.append("Run LEM offered no seepage-time selector on a transient model")
    elif opts["seep_time"]["time"] != 360.0:
        fails.append(f"the selector defaulted to t={opts['seep_time']['time']}, want "
                     f"the last saved frame 360.0")
    note = dlg.seep_time.note.text()
    if "last saved frame" not in note.lower():
        fails.append("a blank stability time is not stated in the dialog: " + note)
    if opts["seep_time"]["remarch"]:
        fails.append("a saved frame was marked as needing a re-march")

    # "use current" is the named version of the play-bar coupling.
    dlg.seep_time.mode.setCurrentIndex(dlg.seep_time.mode.findData("current"))
    if dlg.options()["seep_time"]["time"] != 15.0:
        fails.append("'frame shown in the results viewer' did not read the viewer's time")

    # free entry: allowed, flagged as a re-march, and refused outside the run
    dlg.seep_time.mode.setCurrentIndex(dlg.seep_time.mode.findData("other"))
    dlg.seep_time.other.setText("100")
    o = dlg.options()["seep_time"]
    if not (o["time"] == 100.0 and o["remarch"]):
        fails.append(f"a free-entry time was not marked for a re-march: {o!r}")
    if dlg.seep_time.problem() is not None:
        fails.append("a valid free-entry time was refused")
    if "re-march" not in dlg.seep_time.note.text():
        fails.append("the re-march cost is not stated: " + dlg.seep_time.note.text())
    dlg.seep_time.other.setText("9999")
    if dlg.seep_time.problem() is None:
        fails.append("a time beyond the run duration was accepted")
    dlg.seep_time.other.setText("nonsense")
    if dlg.seep_time.problem() is None:
        fails.append("an unreadable seepage time was accepted")

    # The run gate: an unreadable seepage time refuses the run in its own words, and
    # a preflight ERROR still speaks first. (This fixture carries one, so the
    # preflight half is silenced to isolate the seepage half.)
    if not dlg._ok.isEnabled():
        if dlg._ok.toolTip() != dlg.preflight.block_reason():
            fails.append("a preflight error was not the reason Run gave")
    dlg.preflight._report = None                  # clean preflight, seepage time only
    dlg._sync_run()
    if dlg._ok.isEnabled():
        fails.append("Run stayed live with an unreadable seepage time")
    if "not a number" not in dlg._ok.toolTip():
        fails.append("Run gave no reason for refusing the seepage time: "
                     + dlg._ok.toolTip())
    dlg.seep_time.other.setText("30")
    if not dlg._ok.isEnabled():
        fails.append("Run stayed blocked after the seepage time was corrected")

    # --- the model's own stability time is what the dialog opens on ---
    sd_declared = load_slope_data(DAM)
    sd_declared["tseep"] = {**sd_declared["tseep"], "stability_time": 30.0}
    d2 = RunLemDialog(None, slope_data=sd_declared, transient=sol)
    if d2.options()["seep_time"]["time"] != 30.0:
        fails.append("the dialog did not open on the model's stability_time")
    d2.seep_time.write_back.setChecked(True)
    if not d2.options()["seep_time"]["write_back"]:
        fails.append("the write-back tick did not reach options()")

    # --- rapid drawdown swaps the single instant for the two stage times ---
    d3 = RunLemDialog(None, slope_data=sd, transient=sol)
    d3.rapid.setChecked(True)
    o3 = d3.options()
    if o3["seep_time"] is not None:
        fails.append("a rapid drawdown still carried a single-instant seepage time")
    if o3["stage_times"] != {"stage_1": 0.0, "stage_2": 47.0, "changed": False}:
        fails.append(f"the stage fields did not prefill from the file: {o3['stage_times']!r}")
    d3.stage_times.stage_2.setText("30")
    o3 = d3.options()
    if not (o3["stage_times"]["stage_2"] == 30.0 and o3["stage_times"]["changed"]):
        fails.append(f"an edited stage time was not reported: {o3['stage_times']!r}")
    d3.stage_times.stage_2.setText("0")
    if d3.stage_times.problem() is None:
        fails.append("stage 2 at or before stage 1 was accepted")
    d3.preflight._report = None
    d3._sync_run()
    if d3._ok.isEnabled():
        fails.append("Run stayed live with stage 2 at or before stage 1")
    d3.stage_times.stage_2.setText("")
    if d3.stage_times.problem() is None:
        fails.append("one blank stage time was accepted")
    d3.stage_times.stage_2.setText("47")

    # --- the FEM dialog carries the same selector ---
    f = RunFemDialog(None, slope_data=sd, transient=sol, current_time=15.0)
    if (f.options().get("seep_time") or {}).get("time") != 360.0:
        fails.append("Run FEM did not carry the seepage-time selector")

    # --- a steady model has no choice to offer ---
    steady = load_slope_data(os.path.join(_REPO, "docs/inputs/slope/xslope_dam.xlsx"))
    d4 = RunLemDialog(None, slope_data=steady)
    if d4.seep_time is not None or d4.options()["seep_time"] is not None:
        fails.append("a steady model was offered a seepage-time selector")

    # --- MainWindow: write-back, the run gate, and the re-march trigger ---
    _warned = []
    _orig_warning = QMessageBox.warning
    QMessageBox.warning = staticmethod(lambda *a, **k: (_warned.append(a[2] if len(a) > 2 else ""),
                                                        QMessageBox.Ok)[1])
    tmpdir = tempfile.mkdtemp(prefix="stability_time_mw_")
    mw = MainWindow()
    try:
        xlsx = os.path.join(tmpdir, "damT.xlsx")
        shutil.copy(DAM, xlsx)
        mw.open_path(xlsx)
        mw.doc.path = None          # no sidecar export from the test
        mw.doc.results["transient_seep"] = {"mode": "transient", "transient": sol,
                                            "frames": [], "seep_data": None}

        # write-back: the dialog's choice becomes the model's, and survives the file.
        mw._store_transient_times({"seep_time": {"time": 30.0, "write_back": True,
                                                 "remarch": False},
                                   "stage_times": None})
        if mw.doc.slope_data["tseep"].get("stability_time") != 30.0:
            fails.append("write-back did not reach the model")
        out = os.path.join(tmpdir, "saved.xlsx")
        save_slope_data_to_xlsx(mw.doc.slope_data, out, template=MASTER)
        if load_slope_data(out)["tseep"].get("stability_time") != 30.0:
            fails.append("the written-back stability time did not survive save/load")

        # without the tick, the file is untouched — the choice governs the run only.
        mw._store_transient_times({"seep_time": {"time": 15.0, "write_back": False,
                                                 "remarch": False},
                                   "stage_times": None})
        if mw.doc.slope_data["tseep"].get("stability_time") != 30.0:
            fails.append("an unticked seepage time was written to the model anyway")

        # stage times are model inputs edited at the point of use: always stored.
        mw._store_transient_times({"seep_time": None,
                                   "stage_times": {"stage_1": 2.0, "stage_2": 30.0,
                                                   "changed": True}})
        if (mw.doc.slope_data["tseep"].get("stage_1"),
                mw.doc.slope_data["tseep"].get("stage_2")) != (2.0, 30.0):
            fails.append("edited stage times did not reach the model")

        # the run gate reads the dialog's instant
        mw.doc.slope_data["seep_u"] = None
        if not mw._apply_transient_analysis_frame(rapid=False, time=47.0):
            fails.append("a valid seepage time was refused")
        if _u_time(mw.doc.slope_data) != 47.0:
            fails.append(f"the run read t={_u_time(mw.doc.slope_data)}, want 47.0")

        # ...and with no instant given, the model's stability_time.
        mw.doc.slope_data["seep_u"] = None
        mw._apply_transient_analysis_frame(rapid=False)
        if _u_time(mw.doc.slope_data) != 30.0:
            fails.append("the run did not fall back to the model's stability_time")

        # THE LOUD FAILURE. Stage times with no saved frames must refuse the run and
        # leave the model exactly as it was — the silent-wrong-pore-pressures bug.
        sentinel = np.full(6, -1.0)
        mw.doc.slope_data["seep_u"] = sentinel
        mw.doc.slope_data["seep_u2"] = sentinel
        mw.doc.slope_data["tseep"] = {**mw.doc.slope_data["tseep"],
                                      "stage_1": 2.0, "stage_2": 999.0}
        _warned.clear()
        undo_before = len(mw.doc.undo_labels())
        if mw._apply_transient_analysis_frame(rapid=True) is not False:
            fails.append("a failed rapid-drawdown staging let the run proceed")
        if not _warned:
            fails.append("a failed staging raised no dialog")
        elif "999" not in _warned[0]:
            fails.append(f"the staging failure does not name the bad time: {_warned[0]}")
        if float(np.asarray(mw.doc.slope_data["seep_u"])[0]) != -1.0:
            fails.append("a failed staging left new pore pressures behind")
        if len(mw.doc.undo_labels()) != undo_before:
            fails.append("a failed staging left a half-finished edit on the undo stack")

        # the re-march trigger fires only for instants this solution never saved
        mw.doc.slope_data["tseep"] = {**mw.doc.slope_data["tseep"],
                                      "stage_1": 2.0, "stage_2": 47.0}
        need = mw._remarch_times_needed({"seep_time": {"time": 30.0}, "stage_times": None})
        if need:
            fails.append(f"a saved frame was scheduled for a re-march: {need}")
        need = mw._remarch_times_needed({"seep_time": {"time": 100.0}, "stage_times": None})
        if need != [100.0]:
            fails.append(f"a free-entry time did not trigger a re-march: {need}")
        need = mw._remarch_times_needed(
            {"seep_time": None,
             "stage_times": {"stage_1": 2.0, "stage_2": 123.0, "changed": True}})
        if need != [123.0]:
            fails.append(f"an unsaved stage time did not trigger a re-march: {need}")
        # t = 0 is the initial condition and t > duration is unreachable; neither is
        # a re-march's business.
        need = mw._remarch_times_needed({"seep_time": {"time": 9999.0}, "stage_times": None})
        if need:
            fails.append(f"a time beyond the duration was sent to a re-march: {need}")
    finally:
        QMessageBox.warning = _orig_warning
        mw.doc._dirty = False
        mw.close()
        shutil.rmtree(tmpdir, ignore_errors=True)
    return fails


CHECKS = [("time resolution + entry surface", test_resolution),
          ("stability_time on the tseep sheet", test_file),
          ("saved-frame schedule + re-march injection", test_schedule),
          ("run dialogs + write-back + the staging gate", test_studio)]


def run():
    """Every check; returns a list of failure strings (empty = pass)."""
    try:
        import PySide6                                    # noqa: F401
        checks = CHECKS
    except Exception:
        print("stability time: PySide6 not installed — Studio checks skipped.")
        checks = CHECKS[:-1]
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
    print("transient stability time:")
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nAll transient stability-time checks passed.")


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    try:
        from PySide6.QtWidgets import QApplication
        _app = QApplication.instance() or QApplication([])
    except Exception:
        pass
    main()
