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

The **rapid-drawdown** mode adds a third way to mislead, and one section pins it:
a drawdown row that is not the drawdown. Each point of that curve is a full
three-stage analysis staged from the march's initial state to the instant, so the
row at the model's own stage-2 time must be the SAME number a direct
``run_lem_analysis(..., rapid=True)`` on the staged model reports — otherwise the
curve is a family of near-drawdowns nobody can trace back to a run they could
make by hand. The stage columns and the printed table are checked with it: a
reported value with no stage 2 and stage 3 beside it hides which stage governed,
which is the reading of a drawdown.

A third section pins the marks the drawn curve carries, which are read the same
way the rows are. The instants are colored by the face their critical circle came
out on, so a face read wrong moves the handover between the two slopes and a face
legend on a march that never left one slope invents a handover; the drawdown band
and the full-pool reference are the two marks the dip is sized against, so a band
where the pool never falls, or a reference taken after the fall began, misplaces
it. That section draws the tutorial's own published tables and needs no march.

Two Studio sections follow the engine ones, because the mode is reached from the
Parametric dialog and a control that offers a run the model cannot make is its own
kind of dropped instant:

  * **the dialog** — the mode list carries the entry, its fields replace the
    parameter picker rather than sitting beside it, and it is disabled with a plain
    reason on each of the three models that cannot make a curve (no march, no
    material reading ``u = seep``, the seepage engine). The drawdown box beside the
    frames is offered only where the model carries d and psi, and holds the
    Re-search toggle on while it is ticked.
  * **grid seeding** — the **Grid search** box on the FS-versus-time page reaches
    the engine as ``search_opts={'seed': 'grid'}``. The option decides which
    mechanism each instant reports, so a box that is ticked in the dialog and
    dropped on the way down is a curve drawn on the wrong family, silently.
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


# ---------------------------------------------------------------------------
# Rapid drawdown at each instant
# ---------------------------------------------------------------------------

def _rapid_model():
    """The ACADS sample given the two things a drawdown needs and nothing else:
    the d / psi strengths on its one material, and stage times.

    The frames are still synthetic and the material still reads ``u = none``, so
    the pore pressures are the same at both stages and the drawdown is a pure
    exercise of the STAGING — which is what this section is about. The factor of
    safety it produces is whatever the geometry gives; nothing here quotes one.
    """
    sd = _model()
    sd['materials'] = [dict(m, d=200.0, psi=15.0) for m in sd['materials']]
    sd['tseep'] = dict(sd.get('tseep') or {}, stage_1=0.0, stage_2=20.0)
    return sd


def _staged_by_hand(sd, sol, t):
    """The same model a rapid point stages, staged here instead — the reference the
    curve's row is measured against."""
    import copy

    from xslope.seep import stage_transient_for_drawdown
    from xslope.water import with_water_loads

    out = copy.copy(sd)
    out['materials'] = copy.deepcopy(sd['materials'])
    out['tseep'] = dict(sd.get('tseep') or {}, stage_1=0.0, stage_2=float(t))
    stage_transient_for_drawdown(out, sol)
    return with_water_loads(out)


def check_a_rapid_row_is_the_drawdown(failures):
    """The row at the model's stage-2 time is the same three-stage answer a direct
    ``run_lem_analysis(..., rapid=True)`` on the staged model gives — value, stages
    and circle. A curve nobody can reproduce by hand is not a record of anything."""
    import contextlib
    import io

    from xslope.search import run_lem_analysis

    sd, sol = _rapid_model(), _frames([0.0, 10.0, 20.0])
    with contextlib.redirect_stdout(io.StringIO()):
        ok, res = S.fs_vs_time(sd, sol, rapid=True, methods=('bishop',),
                               num_slices=20, check_inputs=False)
    if not ok:
        failures.append(f"the rapid-drawdown curve was refused: {res}")
        return
    if not res.get('rapid'):
        failures.append("a rapid run does not report itself as one")
    if res.get('stage_1_time') != 0.0:
        failures.append(f"stage 1 was read at t = {res.get('stage_1_time')}, not "
                        f"the march's initial state")
    df = res['df']
    for col in ('stage1_FS', 'stage2_FS', 'stage3_FS', 'stage3_run', 'governs'):
        if col not in df.columns:
            failures.append(f"a drawdown row carries no {col} — the reported value "
                            f"is the lower of stages 2 and 3 and nothing says which")
            return

    row = df.loc[(df['value'] - 20.0).abs() < 1e-9]
    if row.empty or not bool(row['success'].iloc[0]):
        failures.append("the stage-2 instant produced no drawdown: "
                        f"{'' if row.empty else row['msg'].iloc[0]}")
        return
    with contextlib.redirect_stdout(io.StringIO()):
        bundle = run_lem_analysis(_staged_by_hand(sd, sol, 20.0), 'bishop',
                                  analysis='auto_search', num_slices=20,
                                  rapid=True, announce=False)
    ref = bundle['results']
    got = row.iloc[0]
    for key, want in (('fs', ref['FS']), ('stage1_FS', ref['stage1_FS']),
                      ('stage2_FS', ref['stage2_FS']),
                      ('stage3_FS', ref['stage3_FS'])):
        if abs(float(got[key]) - float(want)) > 1e-6:
            failures.append(f"the t = 20 row's {key} is {float(got[key]):.4f}, but "
                            f"a direct rapid run on the staged model gives "
                            f"{float(want):.4f}")
    if bool(got['stage3_run']) != bool(ref.get('stage3_run')):
        failures.append("the row disagrees with the direct run about whether "
                        "stage 3 was required")
    want_gov = 3 if (ref.get('stage3_run')
                     and ref['stage3_FS'] < ref['stage2_FS']) else 2
    if int(got['governs']) != want_gov:
        failures.append(f"the row says stage {int(got['governs'])} governs; the "
                        f"stage values say {want_gov}")

    # Stage 1 is not a drawdown: a fall from the pool to itself.
    first = df.loc[(df['value'] - 0.0).abs() < 1e-9]
    if first.empty:
        failures.append("the stage-1 instant is missing from the table")
    elif bool(first['success'].iloc[0]):
        failures.append("a drawdown from t = 0 to t = 0 was answered")
    elif 'drawdown' not in str(first['msg'].iloc[0]):
        failures.append(f"the stage-1 instant's refusal does not say why: "
                        f"{first['msg'].iloc[0]!r}")


def check_the_table_is_printed(failures):
    """The run's record reaches the console once, with a column per stage, and says
    'not required' where stage 3 did not run rather than repeating stage 2's number
    as though it had."""
    import contextlib
    import io

    sd, sol = _rapid_model(), _frames([0.0, 10.0])
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        ok, res = S.fs_vs_time(sd, sol, rapid=True, methods=('bishop',),
                               num_slices=20, check_inputs=False)
    if not ok:
        failures.append(f"refused: {res}")
        return
    text = out.getvalue()
    if res.get('table_text') not in text:
        failures.append("result['table_text'] is not what was printed")
    header = res.get('table_text', '').splitlines()[:1]
    if not header:
        failures.append("the run printed no table")
        return
    if text.count(header[0]) != 1:
        failures.append(f"the header line appears {text.count(header[0])} times; "
                        f"the table is printed once, after the march")
    for col in ('stage 1', 'stage 2', 'stage 3', 'FS', 'governs'):
        if col not in header[0]:
            failures.append(f"the printed table has no {col!r} column: {header[0]!r}")
    rows = res.get('table')
    if not isinstance(rows, list) or len(rows) != len(res['df']):
        failures.append(f"result['table'] carries {rows and len(rows)} rows for "
                        f"{len(res['df'])} instants")
        return
    for r in rows:
        for key in ('time', 'fs', 'stage1_FS', 'stage2_FS', 'stage3_FS',
                    'stage3_run', 'governs', 'Xo', 'Yo', 'R', 'success', 'msg'):
            if key not in r:
                failures.append(f"a table row carries no {key!r}: {sorted(r)}")
                break
    solved = [r for r in rows if r['success']]
    if not solved:
        failures.append("no instant produced a drawdown at all")
        return

    # An instant where stage 3 did not run says so, rather than repeating stage 2's
    # number where a stage-3 answer belongs — which would read as a stage 3 that ran
    # and agreed. Written against the formatter directly so both cases are on the
    # page whichever way the solves happen to fall.
    import pandas as pd
    made_up = pd.DataFrame([
        {'value': 1.0, 'method': 'bishop', 'fs': 1.20, 'success': True, 'msg': '',
         'Xo': 10.0, 'Yo': 20.0, 'R': 5.0, 'stage1_FS': 1.50, 'stage2_FS': 1.20,
         'stage3_FS': 1.20, 'stage3_run': False, 'governs': 2},
        {'value': 2.0, 'method': 'bishop', 'fs': 1.05, 'success': True, 'msg': '',
         'Xo': 10.0, 'Yo': 20.0, 'R': 5.0, 'stage1_FS': 1.50, 'stage2_FS': 1.10,
         'stage3_FS': 1.05, 'stage3_run': True, 'governs': 3}])
    made_rows, made_text = S._fs_vs_time_table(made_up, {'time_unit': 'day'},
                                               rapid=True)
    lines = made_text.splitlines()
    if len(lines) != 3:
        failures.append(f"two instants made {len(lines) - 1} table rows")
    elif 'not required' not in lines[1]:
        failures.append(f"the instant where stage 3 did not run does not say so: "
                        f"{lines[1]!r}")
    elif '1.0500' not in lines[2]:
        failures.append(f"the instant where stage 3 DID run does not carry its "
                        f"value: {lines[2]!r}")
    if 't (day)' not in lines[0]:
        failures.append(f"the time column does not carry the model's unit: "
                        f"{lines[0]!r}")
    if len({len(ln) for ln in lines}) != 1:
        failures.append(f"the table's lines are {sorted({len(ln) for ln in lines})} "
                        f"characters wide — a column does not line up")
    if [r['governs'] for r in made_rows] != [2, 3]:
        failures.append(f"the returned rows lose the governing stage: {made_rows}")

    # ...and the single-stage table, which has no stages to print.
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        ok, res = _run(sd, sol)
    if not ok:
        failures.append(f"the single-stage run was refused: {res}")
        return
    text = out.getvalue()
    if 'stage 2' in text or 'governs' in text:
        failures.append("the single-stage table carries drawdown columns")
    if not res.get('table') or res.get('rapid'):
        failures.append("a single-stage run carries no table, or claims to be a "
                        "drawdown")
    if 'FS' not in text:
        failures.append(f"the single-stage table has no FS column: {text!r}")
    # Silenceable, without losing the rows. (The solvers still announce themselves;
    # what must be gone is the table.)
    quiet = io.StringIO()
    with contextlib.redirect_stdout(quiet):
        ok, res = _run(sd, sol, print_table=False)
    if not ok:
        failures.append(f"the silenced run was refused: {res}")
        return
    if not res.get('table'):
        failures.append("print_table=False dropped the rows as well as the print")
    elif res['table_text'].splitlines()[0] in quiet.getvalue():
        failures.append(f"print_table=False still printed the table: "
                        f"{quiet.getvalue()!r}")
    if 'versus time' in quiet.getvalue():
        failures.append("print_table=False still printed the table's title")


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
    # The three-stage procedure is a limit-equilibrium construction; SSRM has no
    # equivalent, so the combination is refused rather than quietly ignored.
    sd_fem = dict(sd, mesh={'nodes': [], 'elements': []})
    ok, msg = S.fs_vs_time(sd_fem, sol, mode='fem', rapid=True)
    if ok or 'limit-equilibrium' not in str(msg):
        failures.append(f"a rapid drawdown under mode='fem' was not refused with "
                        f"its reason: {msg!r}")


# ---------------------------------------------------------------------------
# Studio: the Parametric dialog's fourth mode
# ---------------------------------------------------------------------------
COMBO03 = Path(__file__).resolve().parent.parent / 'docs' / 'tutorials' / \
    'files' / 'xslope_earth_dam_fs_time.xlsx'
#: The nineteen instants COMBO-3's march saves, and the two factors of safety the
#: page quotes. Literals: the check never computes its own reference.
COMBO03_TIMES = (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 47.0,
                 55.0, 65.0, 80.0, 100.0, 130.0, 180.0, 240.0, 300.0)
#: COMBO-3 Part 1 publishes these under the run the page makes: Spencer, 40
#: slices, searched at every instant from the two circles the file carries, with
#: Grid search OFF. The full-pool answer is on the dam's DOWNSTREAM face and the
#: minimum on its upstream one, which is the reading the page is built around, so
#: a check that reached only one face would pin neither number.
COMBO03_FULL_POOL_FS = 1.5311
COMBO03_MIN_FS, COMBO03_CRITICAL_TIME = 1.3313, 35.0
COMBO03_METHOD, COMBO03_SLICES, COMBO03_DIVISIONS = 'spencer', 40, 64


# ---------------------------------------------------------------------------
# The curve's marks: the faces, the drawdown band, the full-pool reference
#
# A curve that hands over from one slope to the other is two mechanisms in
# sequence rather than one moving surface, and the plot says so by coloring each
# instant's marker with the face its critical circle came out on. Both ways that
# can mislead are pinned here:
#
#   * **a face read wrong** puts the handover at the wrong instant, which is the
#     one thing the coloring exists to place.
#   * **a face legend on a march that never left one slope** claims a handover
#     that did not happen — worse than no colors, because it reads as a finding.
#
# The drawdown band and the full-pool reference are the two marks the dip is
# measured against, so they are checked beside them: a band drawn where the pool
# never falls, or a full-pool line taken from an instant after the fall began,
# both misplace the dip they exist to size.
#
# The circles below are the ones the tutorial's own tables publish, as literals —
# this section never computes its own reference — so it needs no march and runs in
# a second. The real march reaches the same plot through the Studio section below.
# ---------------------------------------------------------------------------
#: COMBO-3 Part 2's model: the Johnson dam, whose drawdown curve stays on the
#: upstream face at every instant of the march.
COMBO03R = Path(__file__).resolve().parent.parent / 'docs' / 'tutorials' / \
    'files' / 'xslope_johnson_fs_time.xlsx'
#: Part 1's published curve, row for row: (t, FS, Xo, Yo, R, face). The dam's crest
#: runs from x = 51 to x = 59, so the face is the center against that span — the
#: reading the page's own **face** column states.
COMBO03_CURVE = (
    (0.0, 1.5311, 103.00, 56.79, 56.79, 'downstream'),
    (2.0, 1.5311, 103.00, 56.79, 56.79, 'downstream'),
    (5.0, 1.5132, 5.20, 65.90, 65.88, 'upstream'),
    (10.0, 1.4572, 5.20, 65.90, 65.88, 'upstream'),
    (15.0, 1.4128, 5.20, 65.90, 65.88, 'upstream'),
    (20.0, 1.3773, 6.25, 60.58, 60.58, 'upstream'),
    (25.0, 1.3509, 7.00, 57.81, 57.75, 'upstream'),
    (30.0, 1.3344, 7.00, 56.91, 56.91, 'upstream'),
    (35.0, 1.3313, 7.00, 56.91, 56.91, 'upstream'),
    (40.0, 1.3436, 7.00, 56.91, 56.91, 'upstream'),
    (47.0, 1.3813, 7.00, 57.16, 57.16, 'upstream'),
    (55.0, 1.4386, 7.00, 57.16, 57.16, 'upstream'),
    (65.0, 1.4828, 7.00, 57.16, 57.16, 'upstream'),
    (80.0, 1.5187, 7.00, 57.16, 57.16, 'upstream'),
    (100.0, 1.5482, 103.00, 56.79, 56.79, 'downstream'),
    (130.0, 1.5516, 103.00, 56.79, 56.79, 'downstream'),
    (180.0, 1.5550, 103.00, 56.79, 56.79, 'downstream'),
    (240.0, 1.5566, 103.00, 56.79, 56.79, 'downstream'),
    (300.0, 1.5572, 103.00, 56.79, 56.79, 'downstream'),
)
#: Part 2's published drawdown table, the same way. The Johnson crest runs from
#: x = 360 to x = 380 and every center lands near x = 245, so the whole curve is
#: upstream. The t = 0 row produced no result (a drawdown from stage 1 to itself)
#: and is carried as one, because a failed row must not be colored either.
COMBO03R_CURVE = (
    (0.0, None, None, None, None, None),
    (5.0, 1.4563, 242.74, 256.34, 165.90, 'upstream'),
    (10.0, 1.3496, 244.06, 253.05, 164.14, 'upstream'),
    (15.0, 1.2704, 246.78, 250.47, 163.54, 'upstream'),
    (20.0, 1.2036, 249.59, 244.90, 159.60, 'upstream'),
    (25.0, 1.1489, 249.59, 243.48, 159.38, 'upstream'),
    (30.0, 1.1051, 249.59, 243.48, 159.76, 'upstream'),
    (35.0, 1.0711, 248.17, 244.19, 161.02, 'upstream'),
    (40.0, 1.0457, 247.47, 244.90, 162.02, 'upstream'),
    (45.0, 1.0278, 244.64, 244.19, 162.68, 'upstream'),
    (50.0, 1.0157, 243.93, 244.90, 163.64, 'upstream'),
    (60.0, 1.0522, 245.35, 242.07, 160.48, 'upstream'),
    (70.0, 1.0743, 245.31, 246.38, 164.30, 'upstream'),
    (80.0, 1.0902, 245.31, 244.96, 163.08, 'upstream'),
    (100.0, 1.1159, 245.31, 243.53, 161.82, 'upstream'),
    (130.0, 1.1414, 245.31, 243.53, 161.82, 'upstream'),
    (170.0, 1.1612, 245.31, 243.53, 161.82, 'upstream'),
    (220.0, 1.1744, 245.31, 243.53, 161.82, 'upstream'),
    (300.0, 1.1848, 245.31, 243.53, 161.82, 'upstream'),
    (400.0, 1.1893, 245.31, 243.53, 161.82, 'upstream'),
    (500.0, 1.1912, 245.31, 243.53, 161.82, 'upstream'),
)
#: The face entries the legend carries, and the two marks that go with a falling
#: schedule, exactly as the plot writes them.
FACE_LABELS = ('critical on the upstream face', 'critical on the downstream face')


def _published_curve(rows, rapid=False):
    """A ``fs_vs_time`` result built from a published table, so the plot is drawn on
    the numbers the page states rather than on a march this check re-runs."""
    import pandas as pd

    recs = []
    for t, fs, Xo, Yo, R, _face in rows:
        rec = {'param': 'time', 'value': float(t), 'rel': np.nan,
               'is_base': False, 'analysis': 'lem', 'method': 'spencer',
               'fs': np.nan if fs is None else float(fs), 'success': fs is not None,
               'msg': '' if fs is not None else 'no result',
               'Xo': np.nan if Xo is None else float(Xo),
               'Yo': np.nan if Yo is None else float(Yo),
               'R': np.nan if R is None else float(R),
               'output': 'FS', 'output_label': 'Factor of Safety'}
        if rapid:
            rec.update({'stage1_FS': np.nan, 'stage2_FS': rec['fs'],
                        'stage3_FS': np.nan, 'stage3_run': False,
                        'governs': np.nan})
        recs.append(rec)
    df = pd.DataFrame(recs)
    ok = df.loc[df['success']]
    i = ok['fs'].astype(float).idxmin() if len(ok) else None
    return {'df': df, 'param': 'time',
            'times': [float(r[0]) for r in rows],
            'critical_time': None if i is None else float(df.at[i, 'value']),
            'min_fs': None if i is None else float(df.at[i, 'fs']),
            'rapid': bool(rapid), 'mode': 'lem', 'output': 'FS',
            'output_label': 'Factor of Safety'}


def _figure_marks(fig):
    """What a drawn curve says: its legend labels, its annotation texts, and how
    many shaded bands it carries.

    The band is the only patch this figure adds to any of its axes -- the curve,
    the schedule, the guides and the references are all lines -- so counting
    patches counts bands.
    """
    ax = fig.axes[0]
    leg = ax.get_legend()
    labels = [t.get_text() for t in leg.get_texts()] if leg else []
    texts = [t.get_text() for a in fig.axes for t in a.texts]
    bands = sum(len(a.patches) for a in fig.axes)
    return labels, texts, bands


def check_the_curve_colors_a_two_face_march(failures):
    """COMBO-3 Part 1's dam hands over from one slope to the other, so every
    instant is colored by its face and the legend names both."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from xslope.plot import _circle_face, _crest_span, plot_fs_vs_time

    if not COMBO03.exists():
        failures.append(f"missing {COMBO03}")
        return
    sd = load_slope_data(str(COMBO03))
    span = _crest_span(sd)
    if span is None or (round(span[0], 6), round(span[1], 6)) != (51.0, 59.0):
        failures.append(f"the dam's crest reads as {span}, not x = 51 to x = 59")
        return
    for t, _fs, Xo, Yo, R, face in COMBO03_CURVE:
        got = _circle_face(sd, span, Xo, Yo, R)
        if got != face:
            failures.append(f"t = {t:g}: a center at x = {Xo:g} reads as {got}, "
                            f"not the {face} face the page publishes")

    fig = plot_fs_vs_time(_published_curve(COMBO03_CURVE), slope_data=sd)
    labels, texts, bands = _figure_marks(fig)
    for want in FACE_LABELS:
        if want not in labels:
            failures.append(f"a two-face curve's legend does not say {want!r}: "
                            f"{labels}")
    if bands != 1:
        failures.append(f"the falling pool gave {bands} shaded band(s), not one")
    if 'drawdown' not in texts:
        failures.append(f"the shaded interval is not labeled 'drawdown': {texts}")
    full_pool = [s for s in texts if s.startswith('full pool,')]
    if full_pool != [f"full pool, {COMBO03_FULL_POOL_FS:.3f}"]:
        failures.append(f"the full-pool reference reads {full_pool}, not the "
                        f"{COMBO03_FULL_POOL_FS:.3f} the curve starts from")
    plt.close(fig)


def check_a_one_face_march_is_not_colored(failures):
    """COMBO-3 Part 2's dam never leaves the upstream face, so the curve is drawn
    plain — face colors there would claim a handover the run never made."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from xslope.plot import _circle_face, _crest_span, plot_fs_vs_time

    if not COMBO03R.exists():
        failures.append(f"missing {COMBO03R}")
        return
    sd = load_slope_data(str(COMBO03R))
    span = _crest_span(sd)
    if span is None or (round(span[0], 6), round(span[1], 6)) != (360.0, 380.0):
        failures.append(f"the Johnson crest reads as {span}, not x = 360 to "
                        f"x = 380")
        return
    for t, _fs, Xo, Yo, R, face in COMBO03R_CURVE:
        if Xo is None:
            continue
        got = _circle_face(sd, span, Xo, Yo, R)
        if got != face:
            failures.append(f"t = {t:g}: a center at x = {Xo:g} reads as {got}, "
                            f"not the {face} face")

    fig = plot_fs_vs_time(_published_curve(COMBO03R_CURVE, rapid=True),
                          slope_data=sd)
    labels, texts, bands = _figure_marks(fig)
    named = [s for s in labels if s in FACE_LABELS]
    if named:
        failures.append(f"a curve that stays on one face named a face in its "
                        f"legend: {named}")
    # The drawdown's own marks are still there: this dam's pool falls.
    if bands != 1 or 'drawdown' not in texts:
        failures.append(f"the drawdown band is missing from the rapid curve: "
                        f"{bands} band(s), {texts}")
    if not [s for s in texts if s.startswith('full pool,')]:
        failures.append(f"the rapid curve carries no full-pool reference: {texts}")
    plt.close(fig)


def check_a_schedule_that_never_falls_has_no_band(failures):
    """The band and the full-pool reference are read off the schedule, so a pool
    that only fills — or one held flat — gets neither."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from xslope.plot import _drawdown_interval, plot_fs_vs_time

    if not COMBO03.exists():
        failures.append(f"missing {COMBO03}")
        return
    sd = load_slope_data(str(COMBO03))
    if _drawdown_interval(sd['tseep']) != (2.0, 47.0):
        failures.append(f"the dam's fall reads as "
                        f"{_drawdown_interval(sd['tseep'])}, not day 2 to day 47")
    filling = dict(sd)
    ts = dict(sd['tseep'])
    ts['series'] = {'pool': [2.0, 2.0, 18.0]}          # the same schedule, rising
    filling['tseep'] = ts
    if _drawdown_interval(ts) is not None:
        failures.append(f"a rising pool reads as a drawdown: "
                        f"{_drawdown_interval(ts)}")
    fig = plot_fs_vs_time(_published_curve(COMBO03_CURVE), slope_data=filling)
    _labels, texts, bands = _figure_marks(fig)
    if bands:
        failures.append(f"a pool that never falls drew {bands} drawdown band(s)")
    if 'drawdown' in texts or [s for s in texts if s.startswith('full pool,')]:
        failures.append(f"a pool that never falls carried drawdown marks: {texts}")
    plt.close(fig)


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


def check_the_dialog_offers_the_drawdown(failures):
    """The drawdown box sits with the frames, is offered only where the model can
    make a drawdown, holds Re-search on while it is ticked, and reaches the runner."""
    if _qt() is None:
        return
    from studio.dialogs import SensitivityDialog

    dry = _model()
    seepy = dict(dry, materials=[dict(m, u='seep') for m in dry['materials']])
    times = list(COMBO03_TIMES)

    # No d / psi anywhere: every stage would read the same strengths.
    dlg = SensitivityDialog(slope_data=seepy, app_mode='lem',
                            transient={'times': times})
    dlg.mode.setCurrentIndex(dlg.mode.findData('fs_vs_time'))
    if dlg.rapid.isEnabled():
        failures.append("the drawdown box is offered on a model that carries "
                        "neither d nor \N{GREEK SMALL LETTER PSI}")
    elif 'd and' not in dlg.rapid.toolTip():
        failures.append(f"the greyed drawdown box does not name the columns: "
                        f"{dlg.rapid.toolTip()!r}")
    if dlg.rapid.isChecked():
        failures.append("the drawdown box is ticked on a model that cannot make one")
    dlg.close()

    # ...and with them.
    rapid_sd = dict(seepy, materials=[dict(m, d=200.0, psi=15.0)
                                      for m in seepy['materials']],
                    tseep={'stage_1': 0.0, 'stage_2': times[-1]})
    dlg = SensitivityDialog(slope_data=rapid_sd, app_mode='lem',
                            transient={'times': times})
    dlg.show()
    dlg.mode.setCurrentIndex(dlg.mode.findData('fs_vs_time'))
    if not dlg.rapid.isEnabled():
        failures.append("the drawdown box is greyed on a model carrying d and "
                        "\N{GREEK SMALL LETTER PSI}")
    if dlg.options().get('rapid'):
        failures.append("the drawdown is on by default; the single-stage curve is")
    dlg.rapid.setChecked(True)
    opts = dlg.options()
    if not opts.get('rapid'):
        failures.append(f"options() does not carry the drawdown: {opts}")
    if not opts.get('search'):
        failures.append("a drawdown run reached the runner with the search off, "
                        "which is not a run this mode makes")
    if dlg.search.isEnabled():
        failures.append("Re-search stays editable under a drawdown, where every "
                        "point is searched from the starting circle")
    if not dlg.search.isChecked():
        failures.append("Re-search was left off under a drawdown")
    if 'drawdown' not in dlg.note.text().lower():
        failures.append(f"the mode note does not say the points are drawdowns: "
                        f"{dlg.note.text()!r}")
    sel = dlg._selection()
    if sel.get('base') != 'rapid':
        failures.append(f"the model checks are not asked the drawdown's own "
                        f"questions: base = {sel.get('base')!r}")
    if [float(t) for t in (sel.get('seep_frame') or {}).get('times', [])] != \
            [0.0, float(times[-1])]:
        failures.append(f"the checks are not told the two staged instants: "
                        f"{sel.get('seep_frame')}")
    dlg.rapid.setChecked(False)
    if not dlg.search.isEnabled():
        failures.append("Re-search stayed greyed after the drawdown was unticked")
    dlg.close()


def check_the_dialog_offers_grid_seeding(failures):
    """**Grid search** on the FS-versus-time page reaches ``fs_vs_time`` as
    ``search_opts={'seed': 'grid'}``.

    The option decides WHICH mechanism each instant reports, so a box that is
    ticked in the dialog and dropped on the way to the engine is a curve drawn on
    the wrong family with nothing in the record to say so. This runs the whole path
    the dialog uses -- ``options()`` into ``SensitivityRunner._run_fs_vs_time`` --
    and reads the keyword ``fs_vs_time`` was actually called with.
    """
    if _qt() is None:
        return
    from studio.dialogs import SensitivityDialog
    from studio import runners as R

    dry = _model()
    seepy = dict(dry, materials=[dict(m, u='seep') for m in dry['materials']])
    times = list(COMBO03_TIMES)
    sol = _frames(times)

    dlg = SensitivityDialog(slope_data=seepy, app_mode='lem',
                            transient={'times': times})
    dlg.mode.setCurrentIndex(dlg.mode.findData('fs_vs_time'))
    if not dlg.grid_seed.isEnabled():
        failures.append("Grid search is greyed on a circular LEM curve, where it "
                        "is the option the mode most often needs")
    if dlg.options().get('grid_seed'):
        failures.append("Grid search is on by default; the seeded search is")
    if 'grid of circle centers' not in dlg.grid_seed.toolTip():
        failures.append(f"Grid search does not carry Run LEM's own tooltip: "
                        f"{dlg.grid_seed.toolTip()!r}")

    # Off: no seeding keyword at all, so the engine's default is untouched.
    seen = {}

    def _spy(*a, **kw):
        seen.clear()
        seen.update(kw)
        return False, 'stopped by the check'

    # The runner imports the engine function inside the call, so the module
    # attribute is what it will reach for.
    import xslope.sensitivity as XS
    orig = XS.fs_vs_time
    XS.fs_vs_time = _spy
    try:
        for want, tick in (({}, False), ({'seed': 'grid'}, True)):
            dlg.grid_seed.setChecked(tick)
            opts = dlg.options()
            if bool(opts.get('grid_seed')) != tick:
                failures.append(f"options() carries grid_seed="
                                f"{opts.get('grid_seed')!r} with the box {tick}")
            runner = R.SensitivityRunner(seepy, opts, transient=sol)
            runner.failed.connect(lambda *_: None)
            import contextlib
            import io
            with contextlib.redirect_stdout(io.StringIO()):
                runner._run_fs_vs_time()
            got = seen.get('search_opts')
            if tick and got != want:
                failures.append(f"Grid search ticked reached the engine as "
                                f"search_opts={got!r}, not {want!r}")
            if not tick and got not in (None, {}):
                failures.append(f"Grid search unticked still sent "
                                f"search_opts={got!r}")
    finally:
        XS.fs_vs_time = orig
    dlg.close()

    # A non-circular model has no circle centers to sweep, so the box goes away.
    noncirc = dict(seepy, circular=False)
    dlg = SensitivityDialog(slope_data=noncirc, app_mode='lem',
                            transient={'times': times})
    dlg.mode.setCurrentIndex(dlg.mode.findData('fs_vs_time'))
    if dlg.grid_seed.isEnabled() or dlg.grid_seed.isChecked():
        failures.append("Grid search is live on a non-circular model, which has "
                        "no circle centers to sweep")
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
        failures.append(f"the march saved {saved}, not the nineteen frames "
                        f"COMBO-3 publishes")
        return

    dlg = SensitivityDialog(slope_data=sd, app_mode='lem',
                            defaults={'method': COMBO03_METHOD,
                                      'num_slices': COMBO03_SLICES},
                            transient=solution)
    dlg.mode.setCurrentIndex(dlg.mode.findData('fs_vs_time'))
    # Grid search stays OFF, as the page runs it: the two published numbers are on
    # opposite faces of the dam and both are reached from the circles sheet alone,
    # so a run that needed the grid seed to find them would not be the page's run.
    opts = dlg.options()
    dlg.close()
    if opts.get('grid_seed'):
        failures.append('the dialog carried Grid search into a run the page makes '
                        'with it off')

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
    # The marks the page's figure shows, on the real result rather than on the
    # published table the section above draws: this dam changes face, its pool
    # falls, and the curve starts at full pool.
    labels, texts, bands = _figure_marks(fig)
    for want in FACE_LABELS:
        if want not in labels:
            failures.append(f"the searched curve's legend does not say {want!r}: "
                            f"{labels}")
    if bands != 1 or 'drawdown' not in texts:
        failures.append(f"the searched curve lost its drawdown band: {bands} "
                        f"band(s), {texts}")
    if f"full pool, {COMBO03_FULL_POOL_FS:.3f}" not in texts:
        failures.append(f"the searched curve's full-pool reference is not "
                        f"{COMBO03_FULL_POOL_FS:.3f}: {texts}")
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
               check_a_rapid_row_is_the_drawdown,
               check_the_table_is_printed,
               check_the_refusals,
               check_the_curve_colors_a_two_face_march,
               check_a_one_face_march_is_not_colored,
               check_a_schedule_that_never_falls_has_no_band,
               check_the_dialog_offers_the_mode,
               check_the_dialog_offers_the_drawdown,
               check_the_dialog_offers_grid_seeding,
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
