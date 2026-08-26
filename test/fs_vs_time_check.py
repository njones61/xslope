"""The FS-versus-time mode's row contract (``xslope.sensitivity.fs_vs_time``).

A factor-of-safety-versus-time curve is read for its SHAPE — where it dips, how
far, how fast it recovers — so the two ways it can mislead are both about rows
rather than about numbers:

  * **a dropped instant.** Lose the frame just after a drawdown and the curve's
    minimum moves to the next instant, which is a healthier slope than the model
    describes, with nothing in the record to say a point went missing. Every
    requested instant must come back as a row, whether or not it produced an
    answer.
  * **a reason that goes missing.** A row that failed and says nothing is worse
    than no row: it reads as a gap in the data rather than as a refusal, and the
    reader has no way to tell a frame that does not exist from a solver that did
    not converge.

Everything below pins one of those two, or the machinery that has to hold for
them to mean anything: the sentinel screen (an inadmissible search scores 9999
and would otherwise be recorded as a very safe slope), the refusal to blend two
frames into a field that solves nothing, one re-march for a whole set of
requested times rather than one per instant, the base-model gate running once
before the first solve, and the critical instant the run reports.

Two Studio sections follow the engine ones, because the mode is reached from the
Parametric dialog and a control that offers a run the model cannot make is its own
kind of dropped instant:

  * **the dialog** — the mode list carries the entry, its fields replace the
    parameter picker rather than sitting beside it, and it is disabled with a plain
    reason on each of the three models that cannot make a curve (no march, no
    material reading ``u = seep``, the seepage engine).
  * **the whole path, end to end** — the COMBO-3 tutorial model meshed, marched and
    swept through the dialog's own options and the Studio runner, asserting the
    published curve. This one solves real seepage and searches at every instant, so
    it costs about two minutes; the rest of this module is seconds.

Both skip cleanly when PySide6 is not installed. The engine sections are file-light
and fast: they load the shipped ACADS sample for a real geometry, build SYNTHETIC
frame ledgers, and solve no seepage at all.

Run directly:  PYTHONPATH=. python3 test/fs_vs_time_check.py
"""

import os
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings('ignore')

from xslope import sensitivity as S
from xslope import seep as SEEP
from xslope.fileio import load_slope_data

SAMPLE = Path(__file__).resolve().parent.parent / 'docs' / 'lem' / 'files' / \
    'xslope_acads_simple.xlsx'


def _model():
    """A real model with no seepage field: the frames below are synthetic, and the
    materials read their pore pressure some other way, so the LEM answer is the
    same at every instant. That is deliberate — this module checks the TABLE, and
    a constant curve makes an unexplained change in it unmistakable."""
    sd = load_slope_data(str(SAMPLE))
    # The shipped sample gained FEM sidecars (a *_mesh.json among them, commit
    # f87cc320) after this check was written, and the loader auto-attaches the
    # mesh. This fixture wants the bare LEM model — the synthetic 8-node frames
    # below are the whole field — so the mesh is dropped.
    sd.pop('mesh', None)
    return sd


def _frames(times, n=8):
    """A transient solution ledger of the shape seep.py hands out."""
    return {'times': [float(t) for t in times],
            'frames': [{'u': np.zeros(n), 'time': float(t)} for t in times]}


def _run(sd, sol, **kw):
    kw.setdefault('methods', ('bishop',))
    kw.setdefault('search', False)          # the stored circle: no search cost
    kw.setdefault('num_slices', 20)
    return S.fs_vs_time(sd, sol, **kw)


def check_every_instant_returns_a_row(failures):
    """No requested instant may go missing, present or absent from the march."""
    sd, sol = _model(), _frames([0.0, 10.0, 20.0, 30.0])
    ok, res = _run(sd, sol)
    if not ok:
        failures.append(f"the whole-march run was refused: {res}")
        return
    df = res['df']
    if len(df) != 4:
        failures.append(f"4 saved frames gave {len(df)} rows")
    got = [float(v) for v in df['value']]
    if got != [0.0, 10.0, 20.0, 30.0]:
        failures.append(f"rows are not the saved instants, in order: {got}")
    if list(res['times']) != [0.0, 10.0, 20.0, 30.0]:
        failures.append(f"result['times'] disagrees with the rows: {res['times']}")
    if not bool(df['success'].all()):
        failures.append(f"a saved frame failed: {df.loc[~df['success'], 'msg'].tolist()}")

    # An instant with no frame is a ROW, not a gap. This is the mutation that
    # matters most: a mode that skipped it would return 3 rows and a curve whose
    # minimum silently moved.
    ok, res = _run(sd, sol, times=[0.0, 15.0, 30.0])
    if not ok:
        failures.append(f"a run naming an unsaved instant was refused: {res}")
        return
    df = res['df']
    if len(df) != 3:
        failures.append(f"3 requested instants gave {len(df)} rows — an instant "
                        f"with no frame was DROPPED rather than reported")
    row = df.loc[(df['value'] - 15.0).abs() < 1e-9]
    if row.empty:
        failures.append("the instant with no saved frame is missing from the table")
    elif bool(row['success'].iloc[0]):
        failures.append("an instant with no saved frame was reported as a success")


def check_a_failed_row_carries_its_reason(failures):
    """Every unsuccessful row must say why, in words a reader can act on."""
    sd, sol = _model(), _frames([0.0, 10.0, 20.0])
    ok, res = _run(sd, sol, times=[0.0, 5.0, 20.0])
    if not ok:
        failures.append(f"refused: {res}")
        return
    df = res['df']
    bad = df.loc[~df['success']]
    if bad.empty:
        failures.append("the instant with no saved frame did not fail at all")
        return
    for _i, r in bad.iterrows():
        msg = str(r['msg'] or '').strip()
        if not msg:
            failures.append(f"the failed row at t = {r['value']:g} carries NO "
                            f"reason — a gap in the curve with nothing to explain it")
        elif 'saved' not in msg and 'frame' not in msg:
            failures.append(f"the failed row at t = {r['value']:g} does not say a "
                            f"frame is missing: {msg!r}")
    if res['n_failed'] != len(bad):
        failures.append(f"n_failed = {res['n_failed']} but {len(bad)} rows failed")
    # ...and the instants that DID have frames still carry results.
    good = df.loc[df['success']]
    if len(good) != 2:
        failures.append(f"one bad instant cost {3 - len(good)} good ones — a "
                        f"failure must not stop the run")


def check_no_field_is_interpolated(failures):
    """The refusal itself: a time between two frames is never served by blending
    them. Without a re-march it is a failed row, with one it is a computed state."""
    sd, sol = _model(), _frames([0.0, 10.0, 20.0])
    ok, res = _run(sd, sol, times=[10.0, 12.0])
    if not ok:
        failures.append(f"refused: {res}")
        return
    row = res['df'].loc[(res['df']['value'] - 12.0).abs() < 1e-9]
    if row.empty or bool(row['success'].iloc[0]):
        failures.append("t = 12 lies between two saved frames and was ANSWERED — "
                        "a field blended between frames is not a solution")
    if res['remarched']:
        failures.append("no seep_data was given, so no re-march was possible, but "
                        "the run reports one")


def check_one_remarch_for_the_whole_set(failures):
    """Times with no frames are served by ONE re-march carrying all of them, not
    by one march per instant — the difference between a minute and ten."""
    sd = _model()
    sol = _frames([0.0, 10.0, 20.0])
    calls = []

    def fake_remarch(seep_data, slope_data, times, **kw):
        calls.append(sorted(float(t) for t in times))
        return _frames(sorted({0.0, 10.0, 20.0} | {float(t) for t in times}))

    real = SEEP.remarch_for_times
    SEEP.remarch_for_times = fake_remarch
    try:
        ok, res = _run(sd, sol, times=[5.0, 10.0, 15.0],
                       seep_data={'stub': True}, remarch=True)
    finally:
        SEEP.remarch_for_times = real
    if not ok:
        failures.append(f"the re-march run was refused: {res}")
        return
    if len(calls) != 1:
        failures.append(f"{len(calls)} re-marches for 2 unsaved instants — the "
                        f"whole set must be served by one")
    elif calls[0] != [5.0, 15.0]:
        failures.append(f"the re-march was asked for {calls[0]}, not the two "
                        f"instants that had no frames")
    if not res['remarched']:
        failures.append("a re-march happened but the result does not say so")
    if not bool(res['df']['success'].all()):
        failures.append("an instant the re-march created still has no result")


def check_sentinel_and_negative_are_not_results(failures):
    """9999 is the search's no-admissible-surface flag and a negative factor of
    safety is not a result. Both must become failed rows carrying their reason."""
    sd, sol = _model(), _frames([0.0, 10.0])

    for value, want in ((9999.0, 'sentinel'), (-0.4, 'negative')):
        def fake_point(_sd, _mode, methods, *a, **kw):
            return [{'method': m, 'fs': value, 'success': True, 'msg': '',
                     'Xo': np.nan, 'Yo': np.nan, 'R': np.nan} for m in methods]
        real = S._run_point
        S._run_point = fake_point
        try:
            ok, res = _run(sd, sol)
        finally:
            S._run_point = real
        if not ok:
            failures.append(f"refused: {res}")
            continue
        df = res['df']
        if bool(df['success'].any()):
            failures.append(f"FS = {value:g} was recorded as a SUCCESS")
        if not all(want in str(m) for m in df['msg']):
            failures.append(f"FS = {value:g} was refused without saying it is a "
                            f"{want}: {list(df['msg'])}")
        if res['min_fs'] is not None:
            failures.append(f"FS = {value:g} reached the reported minimum")


def check_the_base_model_is_gated_once(failures):
    """A defect in the model is a defect in every instant: it is reported once, at
    the door, and no frame is ever solved."""
    sd, sol = _model(), _frames([0.0, 10.0, 20.0])
    solved = []

    def counting_point(*a, **kw):
        solved.append(1)
        return [{'method': 'bishop', 'fs': 1.0, 'success': True, 'msg': '',
                 'Xo': np.nan, 'Yo': np.nan, 'R': np.nan}]

    real_gate, real_point = S._sweep_gate, S._run_point
    S._sweep_gate = lambda *a, **kw: "the model cannot be analysed"
    S._run_point = counting_point
    try:
        ok, res = _run(sd, sol)
    finally:
        S._sweep_gate, S._run_point = real_gate, real_point
    if ok:
        failures.append("a refused base model still produced a curve")
    elif 'cannot be analysed' not in str(res):
        failures.append(f"the refusal lost the gate's own words: {res!r}")
    if solved:
        failures.append(f"{len(solved)} instants were solved AFTER the base model "
                        f"was refused")


def check_the_reported_minimum(failures):
    """critical_time / min_fs are the curve's own minimum over the rows that
    produced one, and a failed instant can never supply it."""
    sd, sol = _model(), _frames([0.0, 10.0, 20.0, 30.0])
    curve = {0.0: 1.60, 10.0: 0.92, 20.0: 1.05, 30.0: 1.20}

    def fake_point(sdx, _mode, methods, *a, **kw):
        fs = curve[float(sdx['seep_u_time'])]
        return [{'method': m, 'fs': fs, 'success': True, 'msg': '',
                 'Xo': np.nan, 'Yo': np.nan, 'R': np.nan} for m in methods]

    real = S._run_point
    S._run_point = fake_point
    try:
        ok, res = _run(sd, sol)
    finally:
        S._run_point = real
    if not ok:
        failures.append(f"refused: {res}")
        return
    if res['critical_time'] != 10.0:
        failures.append(f"critical time is {res['critical_time']} but the curve "
                        f"dips at t = 10")
    if abs(float(res['min_fs']) - 0.92) > 1e-9:
        failures.append(f"minimum FS is {res['min_fs']} but the curve dips to 0.92")
    if res['critical'].get('bishop', {}).get('time') != 10.0:
        failures.append(f"the per-method critical instant disagrees: {res['critical']}")


def check_the_refusals(failures):
    """The three arguments that cannot make a curve."""
    sd, sol = _model(), _frames([0.0, 10.0])
    ok, msg = S.fs_vs_time(sd, sol, mode='seep')
    if ok or 'seep' not in str(msg):
        failures.append(f"mode='seep' was not refused with its reason: {msg!r}")
    ok, msg = S.fs_vs_time(sd, None)
    if ok or 'transient' not in str(msg):
        failures.append(f"a missing transient solution was not refused: {msg!r}")
    ok, msg = S.fs_vs_time(sd, {'times': [], 'frames': []})
    if ok or 'frame' not in str(msg):
        failures.append(f"a march with no saved frames was not refused: {msg!r}")


# ---------------------------------------------------------------------------
# Studio: the Parametric dialog's fourth mode
# ---------------------------------------------------------------------------
COMBO03 = Path(__file__).resolve().parent.parent / 'docs' / 'tutorials' / \
    'files' / 'xslope_earth_dam_fs_time.xlsx'
#: The twelve instants COMBO-3's march saves, and the two factors of safety the
#: page quotes. Literals: the check never computes its own reference.
COMBO03_TIMES = (0.0, 2.0, 15.0, 30.0, 47.0, 60.0, 80.0, 120.0, 180.0, 240.0,
                 300.0, 360.0)
COMBO03_FULL_POOL_FS = 1.8308
COMBO03_MIN_FS, COMBO03_CRITICAL_TIME = 1.4962, 30.0
COMBO03_METHOD, COMBO03_SLICES, COMBO03_DIVISIONS = 'spencer', 40, 64


def _qt():
    """The QApplication, or None when Studio is not installed."""
    os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
    try:
        from PySide6.QtWidgets import QApplication
    except Exception:                                           # noqa: BLE001
        return None
    return QApplication.instance() or QApplication([])


def check_the_dialog_offers_the_mode(failures):
    """The mode is in the list, its fields replace the picker, and it refuses with
    a reason on each model that cannot make a curve."""
    if _qt() is None:
        return
    from studio.dialogs import SensitivityDialog

    dry = _model()                 # the shipped sample reads its u some other way
    sd = dict(dry, materials=[dict(m, u='seep') for m in dry['materials']])
    times = list(COMBO03_TIMES)
    live = SensitivityDialog(slope_data=sd, app_mode='lem',
                             transient={'times': times})
    idx = live.mode.findData('fs_vs_time')
    if idx < 0:
        failures.append("the Parametric dialog's mode list has no FS-versus-time "
                        "entry")
        return
    if live.mode.itemText(idx) != 'Factor of safety vs time':
        failures.append(f"the mode is labelled {live.mode.itemText(idx)!r}")
    if not live.mode.model().item(idx).isEnabled():
        failures.append("the mode is disabled on a model that carries a march and "
                        "reads u = seep")
    live.show()
    live.mode.setCurrentIndex(idx)
    if live.picker.isVisible():
        failures.append("the parameter picker is still shown in FS-versus-time "
                        "mode, where no input is substituted")
    if live.selected_times() != times:
        failures.append(f"the saved frames are not all ticked by default: "
                        f"{live.selected_times()}")
    opts = live.options()
    if opts.get('mode') != 'fs_vs_time' or opts.get('times') != times:
        failures.append(f"options() does not carry the mode and its times: {opts}")
    live._set_all_times(False)
    if live._ok.isEnabled():
        failures.append("Run stays enabled with no instant ticked")
    live.close()

    # The three refusals, each with its own sentence.
    for label, kw, want in (
            ("no transient solution", dict(app_mode='lem', transient=None),
             'transient seepage'),
            ("the seepage engine", dict(app_mode='seep',
                                        transient={'times': times}), 'LEM or FEM'),
    ):
        dlg = SensitivityDialog(slope_data=sd, **kw)
        i = dlg.mode.findData('fs_vs_time')
        item = dlg.mode.model().item(i)
        if item.isEnabled():
            failures.append(f"the mode is offered with {label}")
        elif want not in item.toolTip():
            failures.append(f"{label}: the reason does not say why — "
                            f"{item.toolTip()!r}")
        dlg.close()
    dlg = SensitivityDialog(slope_data=dry, app_mode='lem',
                            transient={'times': times})
    i = dlg.mode.findData('fs_vs_time')
    item = dlg.mode.model().item(i)
    if item.isEnabled():
        failures.append("the mode is offered on a model where no material reads "
                        "u = seep")
    elif 'u = seep' not in item.toolTip():
        failures.append(f"the u = seep reason does not name the column — "
                        f"{item.toolTip()!r}")
    dlg.close()


def check_the_combo03_curve(failures):
    """COMBO-3 end to end: mesh, march, and the dialog's own options through the
    Studio runner, against the two factors of safety the tutorial publishes."""
    if _qt() is None or not COMBO03.exists():
        return
    import contextlib
    import io

    from xslope.mesh import (build_mesh_from_polygons, extract_size_regions,
                             get_material_polygons)
    from xslope.seep import build_seep_data, build_tseep_data, run_transient_seepage
    from studio.dialogs import SensitivityDialog
    from studio.runners import SensitivityRunner

    with contextlib.redirect_stdout(io.StringIO()):
        sd = load_slope_data(str(COMBO03))
        xs = [x for x, _ in sd['ground_surface'].coords]
        sd['mesh'] = build_mesh_from_polygons(
            get_material_polygons(sd),
            (max(xs) - min(xs)) / COMBO03_DIVISIONS, 'tri3',
            size_regions=extract_size_regions(sd))
        solution = run_transient_seepage(
            build_seep_data(sd['mesh'], sd, seep_bc=1), build_tseep_data(sd),
            verbose=False)
    saved = [float(t) for t in solution['times']]
    if saved != list(COMBO03_TIMES):
        failures.append(f"the march saved {saved}, not the twelve frames COMBO-3 "
                        f"publishes")
        return

    dlg = SensitivityDialog(slope_data=sd, app_mode='lem',
                            defaults={'method': COMBO03_METHOD,
                                      'num_slices': COMBO03_SLICES},
                            transient=solution)
    dlg.mode.setCurrentIndex(dlg.mode.findData('fs_vs_time'))
    opts = dlg.options()
    dlg.close()

    runner = SensitivityRunner(sd, opts, transient=solution)
    got = {}
    runner.succeeded.connect(lambda b: got.update(b))
    runner.failed.connect(lambda m: got.setdefault('error', m))
    with contextlib.redirect_stdout(io.StringIO()):
        runner._run_fs_vs_time()
    if 'error' in got or not got:
        failures.append(f"the Studio FS-versus-time run failed: "
                        f"{got.get('error', 'no result')}")
        return
    df = got['df']
    if len(df) != len(COMBO03_TIMES):
        failures.append(f"the curve carries {len(df)} rows for "
                        f"{len(COMBO03_TIMES)} instants")
    for t, want in ((0.0, COMBO03_FULL_POOL_FS),
                    (COMBO03_CRITICAL_TIME, COMBO03_MIN_FS)):
        row = df.loc[(df['value'] - t).abs() < 1e-9]
        if row.empty or not bool(row['success'].iloc[0]):
            failures.append(f"t = {t:g}: no factor of safety in the curve")
            continue
        if abs(float(row['fs'].iloc[0]) - want) > 0.005:
            failures.append(f"t = {t:g}: COMBO-3 publishes {want:.4f}, the run "
                            f"reports {float(row['fs'].iloc[0]):.4f}")
    if got.get('critical_time') != COMBO03_CRITICAL_TIME:
        failures.append(f"the critical instant is {got.get('critical_time')}, not "
                        f"day {COMBO03_CRITICAL_TIME:g}")
    if got.get('min_fs') is None or abs(got['min_fs'] - COMBO03_MIN_FS) > 0.005:
        failures.append(f"the minimum is {got.get('min_fs')}, not "
                        f"{COMBO03_MIN_FS:.4f}")
    # The plot the result tab draws, on the real result.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from xslope.plot import plot_fs_vs_time
    fig = plot_fs_vs_time(got, slope_data=sd)
    if not fig.axes:
        failures.append("plot_fs_vs_time drew no axes")
    plt.close(fig)


def run():
    failures = []
    if not SAMPLE.exists():
        return [f"missing {SAMPLE}"]
    for fn in (check_every_instant_returns_a_row,
               check_a_failed_row_carries_its_reason,
               check_no_field_is_interpolated,
               check_one_remarch_for_the_whole_set,
               check_sentinel_and_negative_are_not_results,
               check_the_base_model_is_gated_once,
               check_the_reported_minimum,
               check_the_refusals,
               check_the_dialog_offers_the_mode,
               check_the_combo03_curve):
        try:
            fn(failures)
        except Exception as e:                                  # noqa: BLE001
            failures.append(f"{fn.__name__} raised {type(e).__name__}: {e}")
    return failures


def main():
    failures = run()
    if failures:
        print("\nFAILURES:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("\nEvery requested instant comes back as a row, every row that produced "
          "no factor of safety says why, and Studio's Parametric mode reproduces "
          "COMBO-3's published curve.")


if __name__ == '__main__':
    main()
