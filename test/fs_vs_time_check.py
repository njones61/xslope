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

File-light and fast: it loads the shipped ACADS sample for a real geometry,
builds SYNTHETIC frame ledgers, and solves no seepage at all.

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
               check_the_refusals):
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
    print("\nEvery requested instant comes back as a row, and every row that "
          "produced no factor of safety says why.")


if __name__ == '__main__':
    main()
