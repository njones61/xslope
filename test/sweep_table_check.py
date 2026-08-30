"""The Table sub-tab every Parametric result view carries (``studio.sweep_table``).

A Parametric plot is a picture of a table the run already computed, and the table
is the half that gets quoted: read off a value, check it against a hand
calculation, paste it into a report. So each result view carries its numbers
beside its figure, and this module pins the three ways that can go wrong:

  * **a short table.** A grid with fewer rows than the run produced is the
    dangerous one, because nothing on screen says a row is missing: drop the
    instant a drawdown curve dips at and the table reports a healthier slope than
    the model describes. Every mode's row count is checked against its own result
    — the instants of a march, the points of a design sweep, the bars of a
    tornado, the parameters of a reliability run.
  * **a CSV that is not the grid.** Save CSV… is what leaves the program, so a
    file that drops a column, re-orders rows or re-formats a number would be
    quoted in place of the run. Every mode's file is read back and compared cell
    for cell with what the grid shows.
  * **a table that disagrees with its plot.** The two sit side by side under one
    tab, so a face column that names the other slope, a drawdown row with no stage
    values behind its reported factor of safety, or a failed instant carried as a
    blank rather than as a reason, all read as findings. Those are checked against
    COMBO-3's own published tables.

The fixtures are the tutorial's published curves and small synthetic sweeps, as
literals — this check never computes its own reference and solves nothing, so it
runs in seconds. Skips cleanly (exit 0) when PySide6 is not installed.
"""

import csv
import os
import sys
import tempfile
from pathlib import Path

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

try:
    from PySide6.QtWidgets import QApplication, QMessageBox
except Exception:                       # engine-only install — nothing to test
    print("sweep table check: PySide6 not installed — skipped.")
    sys.exit(0)

_app = QApplication.instance() or QApplication([])
QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.information = staticmethod(lambda *a, **k: QMessageBox.Ok)
QMessageBox.critical = staticmethod(lambda *a, **k: QMessageBox.Ok)

from xslope.fileio import load_slope_data
from xslope.sensitivity import tornado_from_sweeps
from studio.main_window import MainWindow
from studio import sweep_table

#: COMBO-3 Part 1's model and its published curve, row for row: (t, FS, Xo, Yo, R,
#: face). The dam's crest runs from x = 51 to x = 59, so the face is the center
#: against that span — the column the tutorial's own table states.
COMBO03 = Path(_REPO) / 'docs' / 'tutorials' / 'files' / \
    'xslope_earth_dam_fs_time.xlsx'
COMBO03_CURVE = (
    (0.0, 1.5311, 103.00, 56.79, 56.79, 'downstream'),
    (5.0, 1.5132, 5.20, 65.90, 65.88, 'upstream'),
    (20.0, 1.3773, 6.25, 60.58, 60.58, 'upstream'),
    (35.0, 1.3313, 7.00, 56.91, 56.91, 'upstream'),
    (80.0, 1.5187, 7.00, 57.16, 57.16, 'upstream'),
    (300.0, 1.5572, 103.00, 56.79, 56.79, 'downstream'),
)
#: Part 2's drawdown march, with an instant that produced no result (a drawdown
#: from stage 1 to itself) — the row that must carry its reason rather than a gap.
COMBO03R_CURVE = (
    (0.0, None, None, None, None),
    (25.0, 1.1489, 249.59, 243.48, 159.38),
    (50.0, 1.0157, 243.93, 244.90, 163.64),
    (500.0, 1.1912, 245.31, 243.53, 161.82),
)

_SHARED = {}


def _dam():
    """COMBO-3 Part 1's model, loaded once — the geometry the face column is read
    against."""
    if 'dam' not in _SHARED:
        _SHARED['dam'] = load_slope_data(str(COMBO03))
    return _SHARED['dam']


# ---------------------------------------------------------------------------
# fixtures: one bundle per Parametric mode, shaped as its runner emits it
# ---------------------------------------------------------------------------

def _march(rows, rapid=False):
    """An ``fs_vs_time`` bundle built from a published table, so the grid is
    checked against the numbers the page states rather than a march re-run here."""
    from xslope.sensitivity import _fs_vs_time_table
    recs = []
    for row in rows:
        t, fs = row[0], row[1]
        Xo, Yo, R = row[2], row[3], row[4]
        rec = {'param': 'time', 'value': float(t), 'rel': np.nan, 'is_base': False,
               'analysis': 'lem', 'method': 'spencer',
               'fs': np.nan if fs is None else float(fs),
               'success': fs is not None,
               'msg': '' if fs is not None else 'no result',
               'Xo': np.nan if Xo is None else float(Xo),
               'Yo': np.nan if Yo is None else float(Yo),
               'R': np.nan if R is None else float(R),
               'output': 'FS', 'output_label': 'Factor of Safety'}
        if rapid:
            rec.update({'stage1_FS': np.nan if fs is None else float(fs) + 0.3,
                        'stage2_FS': rec['fs'], 'stage3_FS': np.nan,
                        'stage3_run': False,
                        'governs': np.nan if fs is None else 2})
        recs.append(rec)
    df = pd.DataFrame(recs)
    table, _text = _fs_vs_time_table(df, _dam(), rapid)
    ok = df.loc[df['success']]
    i = ok['fs'].astype(float).idxmin() if len(ok) else None
    return {'kind': 'fs_vs_time', 'df': df, 'table': table, 'rapid': bool(rapid),
            'param': 'time', 'times': [float(r[0]) for r in rows],
            'critical_time': None if i is None else float(df.at[i, 'value']),
            'min_fs': None if i is None else float(df.at[i, 'fs']),
            'method': 'spencer', 'mode': 'lem', 'output': 'FS',
            'output_label': 'Factor of Safety'}


def _sweep_df(param, base_value, values, fs_at):
    rows = [{'param': param, 'value': base_value, 'rel': 0.0, 'is_base': True,
             'analysis': 'lem', 'method': 'spencer', 'fs': fs_at(base_value),
             'success': True, 'msg': '', 'Xo': 60.0, 'Yo': 90.0, 'R': 55.0}]
    for v in values:
        rows.append({'param': param, 'value': v,
                     'rel': (v - base_value) / base_value, 'is_base': False,
                     'analysis': 'lem', 'method': 'spencer', 'fs': fs_at(v),
                     'success': True, 'msg': '', 'Xo': 60.0, 'Yo': 90.0,
                     'R': 55.0})
    df = pd.DataFrame(rows)
    df['output'] = 'FS'
    df['output_label'] = 'Factor of Safety'
    return df


def _sweeps():
    return {'mat:Clay:c': _sweep_df('mat:Clay:c', 20.0, [10.0, 20.0, 30.0],
                                    lambda v: 1.0 + 0.02 * v),
            'mat:Clay:phi': _sweep_df('mat:Clay:phi', 30.0, [25.0, 30.0, 35.0],
                                      lambda v: 0.9 + 0.015 * v)}


def _bundles():
    """One bundle per Parametric mode, keyed by the name this check reports it
    under, with the row count its own result implies."""
    sweeps = _sweeps()
    tornado = tornado_from_sweeps(sweeps, base_fs=1.40, method='spencer')
    design_df = _sweep_df('mat:Clay:c', 20.0, [10.0, 15.0, 20.0, 25.0, 30.0],
                          lambda v: 1.0 + 0.02 * v)
    scaled = {'bars': [{'param': 'mat:Clay:c', 'base_value': 20.0, 'fs_base': 1.4,
                        'dfdp': 0.02, 'elasticity': 0.286, 'per_1pct': 0.004,
                        'per_sigma': 0.04, 'sigma': 2.0, 'sign': 1},
                       {'param': 'mat:Clay:phi', 'base_value': 30.0,
                        'fs_base': 1.4, 'dfdp': 0.015, 'elasticity': 0.321,
                        'per_1pct': 0.0045, 'per_sigma': 0.03, 'sigma': 2.0,
                        'sign': 1}],
              'skipped': [], 'fs_base': 1.4, 'method': 'spencer', 'mode': 'lem',
              'rel_step': 0.01, 'output': 'FS', 'has_sigma': True}
    variance = {'bars': [{'param': 'mat:Clay:c', 'label': 'Clay · c',
                          'variance': 0.0064, 'delta_F': 0.16, 'sigma': 2.0,
                          'pct': 64.0, 'cumulative': 64.0},
                         {'param': 'mat:Clay:phi', 'label': 'Clay · phi',
                          'variance': 0.0036, 'delta_F': 0.12, 'sigma': 2.0,
                          'pct': 36.0, 'cumulative': 100.0}],
                'sigma_F': 0.1, 'F_MLV': 1.4, 'COV_F': 0.071, 'method': 'spencer'}
    rank = {'bars': [{'param': 'mat:Clay:c', 'label': 'Clay · c', 'rho': 0.62,
                      'sign': 1, 'sigma': 2.0},
                     {'param': 'mat:Clay:phi', 'label': 'Clay · phi',
                      'rho': -0.31, 'sign': -1, 'sigma': 2.0}],
            'method': 'bishop', 'n_valid': 4980, 'n_samples': 5000,
            'distribution': 'normal'}
    return {
        'tornado': ({'kind': 'sensitivity', 'plot_type': 'tornado',
                     'sweeps': sweeps, 'tornado': tornado, 'method': 'spencer'},
                    len(sweeps)),
        'spider': ({'kind': 'sensitivity', 'plot_type': 'spider',
                    'sweeps': sweeps, 'tornado': tornado, 'method': 'spencer'},
                   sum(len(df) for df in sweeps.values())),
        'scaled': ({'kind': 'sensitivity', 'plot_type': 'scaled',
                    'scaled': scaled, 'scaling': 'elasticity',
                    'method': 'spencer'}, len(scaled['bars'])),
        'variance': ({'kind': 'sensitivity', 'plot_type': 'variance',
                      'variance': variance, 'method': 'spencer'},
                     len(variance['bars'])),
        'rank': ({'kind': 'sensitivity', 'plot_type': 'rank', 'rank': rank,
                  'method': 'bishop'}, len(rank['bars'])),
        'design': ({'kind': 'design', 'study': 'design', 'df': design_df,
                    'param': 'mat:Clay:c', 'target_fs': 1.5, 'crossing': 25.0,
                    'bracketed': True, 'output': 'FS',
                    'message': 'FS = 1.5 at mat:Clay:c = 25'}, len(design_df)),
        'back_analysis': ({'kind': 'design', 'study': 'back_analysis',
                           'df': design_df, 'param': 'mat:Clay:c',
                           'target_fs': 1.0, 'crossing': 12.0, 'bracketed': True,
                           'output': 'FS', 'message': 'FS = 1.0 at 12'},
                          len(design_df)),
        'fs_vs_time': (_march(COMBO03_CURVE), len(COMBO03_CURVE)),
        'drawdown_vs_time': (_march(COMBO03R_CURVE, rapid=True),
                             len(COMBO03R_CURVE)),
    }


#: The Taylor-series and Monte Carlo reliability results, as their engines return
#: them — the two "reliability summary" views.
RELIABILITY_TAYLOR = {
    'param_info': [
        {'material_id': 1, 'material_name': 'Clay', 'param': 'c', 'mlv': 20.0,
         'std': 2.0, 'F_plus': 1.48, 'F_minus': 1.32, 'delta_F': 0.16},
        {'material_id': 1, 'material_name': 'Clay', 'param': 'phi', 'mlv': 30.0,
         'std': 2.0, 'F_plus': 1.46, 'F_minus': 1.34, 'delta_F': 0.12},
    ],
    'F_MLV': 1.40, 'sigma_F': 0.10, 'COV_F': 0.0714, 'beta_ln': 3.42,
    'reliability': 0.9997, 'prob_failure': 0.0003, 'fs_cache': None,
}
RELIABILITY_MC = {
    'mean_FS': 1.402, 'sigma_F': 0.101, 'COV_F': 0.072, 'beta_normal': 3.98,
    'beta_ln': 4.21, 'pf_empirical': 0.0002, 'pf_normal': 0.00003,
    'pf_lognormal': 0.00001, 'n_valid': 4980, 'n_samples': 5000,
    'distribution': 'normal', 'fs_samples': None,
}


def _window():
    """A Studio window carrying every Parametric result, each shown through the
    same ``_show_*`` the runner calls on success."""
    if 'win' in _SHARED:
        return _SHARED['win']
    win = MainWindow()
    win.doc.slope_data = _dam()
    bundles = _bundles()
    views = {}
    for name in ('tornado', 'spider', 'scaled', 'variance', 'rank'):
        win.doc.results['sensitivity'] = bundles[name][0]
        win._show_sensitivity()
        views[name] = win.sens_canvas
        # Each plot type re-uses the one Sensitivity tab, so its table is read
        # while it is the one showing.
        views[name + '_rows'] = win.sens_canvas.table_view.rows()
        views[name + '_headers'] = win.sens_canvas.table_view.headers()
    for name in ('design', 'back_analysis'):
        win.doc.results['design'] = bundles[name][0]
        win._show_design()
        views[name] = win.design_canvas
        views[name + '_rows'] = win.design_canvas.table_view.rows()
        views[name + '_headers'] = win.design_canvas.table_view.headers()
    for name in ('fs_vs_time', 'drawdown_vs_time'):
        win.doc.results['fs_vs_time'] = bundles[name][0]
        win._show_fs_vs_time()
        views[name] = win.fs_time_canvas
        views[name + '_rows'] = win.fs_time_canvas.table_view.rows()
        views[name + '_headers'] = win.fs_time_canvas.table_view.headers()
    win.doc.results['reliability'] = {'engine': 'taylor', 'app_mode': 'lem',
                                      'reliability': RELIABILITY_TAYLOR}
    win._show_reliability(RELIABILITY_TAYLOR)
    views['reliability_taylor'] = win.reliability_canvas
    views['reliability_taylor_rows'] = win.reliability_canvas.table_view.rows()
    views['reliability_taylor_headers'] = \
        win.reliability_canvas.table_view.headers()
    win.doc.results['reliability'] = {'engine': 'mc', 'app_mode': 'lem',
                                      'reliability': RELIABILITY_MC}
    win._show_reliability_histogram()
    views['reliability_mc'] = win.reliability_hist_canvas
    views['reliability_mc_rows'] = win.reliability_hist_canvas.table_view.rows()
    views['reliability_mc_headers'] = \
        win.reliability_hist_canvas.table_view.headers()
    # The click-through curve, reached the way a user reaches it: a double-click on
    # a tornado bar. It needs the drawn axes, so the tornado is rendered first.
    win.doc.results['sensitivity'] = bundles['tornado'][0]
    win._show_sensitivity()
    win.sens_canvas.canvas.render_now()
    win._on_tornado_pick(0.0, 0.0, 0.1)
    views['sensitivity_curve'] = win.sens_curve_canvas
    if win.sens_curve_canvas is not None:
        views['sensitivity_curve_rows'] = win.sens_curve_canvas.table_view.rows()
        views['sensitivity_curve_headers'] = \
            win.sens_curve_canvas.table_view.headers()
    _SHARED['win'] = win
    _SHARED['views'] = views
    _SHARED['bundles'] = bundles
    return win


def _views():
    _window()
    return _SHARED['views']


# ---------------------------------------------------------------------------
# A. every Parametric result view carries the two sub-tabs
# ---------------------------------------------------------------------------

def check_every_view_has_a_table_sub_tab(failures):
    """A Parametric result view is its plot and its table, in that order.

    Plot is tab 0 and the one a fresh result opens on — the figure is what a sweep
    is run to see. A mode whose view carries only a figure leaves its numbers
    reachable nowhere but the Log pane, which is the state this whole feature
    replaces.
    """
    views = _views()
    for name in ('tornado', 'spider', 'scaled', 'variance', 'rank', 'design',
                 'back_analysis', 'fs_vs_time', 'drawdown_vs_time',
                 'reliability_taylor', 'reliability_mc', 'sensitivity_curve'):
        view = views.get(name)
        if view is None:
            failures.append(f"{name}: no result view was built")
            continue
        if not isinstance(view, sweep_table.SweepResultView):
            failures.append(f"{name}: result view is {type(view).__name__}, "
                            "not the shared plot/table view")
            continue
        labels = [view.tabs.tabText(i) for i in range(view.tabs.count())]
        if labels != ["Plot", "Table"]:
            failures.append(f"{name}: sub-tabs are {labels}, expected "
                            "['Plot', 'Table']")
        if view.tabs.currentIndex() != 0:
            failures.append(f"{name}: opens on sub-tab "
                            f"{view.tabs.currentIndex()}, not the plot")
        if not hasattr(view.table_view, "save_csv"):
            failures.append(f"{name}: the Table sub-tab carries no Save CSV…")


# ---------------------------------------------------------------------------
# B. the grid holds every row the result holds
# ---------------------------------------------------------------------------

def check_the_row_count_is_the_results_own(failures):
    """One row per thing the run produced — instant, sweep point, bar, parameter.

    A grid shorter than its result is the failure with no symptom: nothing on
    screen says a row went missing, so the shortest table always reads as the
    whole answer.
    """
    views, bundles = _views(), _SHARED['bundles']
    for name, (_bundle, expected) in bundles.items():
        rows = views.get(name + '_rows')
        if rows is None:
            failures.append(f"{name}: no table rows were captured")
            continue
        if len(rows) != expected:
            failures.append(f"{name}: {len(rows)} table row(s), expected "
                            f"{expected}")
    pairs = (('reliability_taylor', len(RELIABILITY_TAYLOR['param_info'])),
             ('reliability_mc', 11),
             ('sensitivity_curve', len(_sweeps()['mat:Clay:c'])))
    for name, expected in pairs:
        rows = views.get(name + '_rows')
        if rows is None:
            failures.append(f"{name}: no table rows were captured")
        elif len(rows) != expected:
            failures.append(f"{name}: {len(rows)} table row(s), expected "
                            f"{expected}")
    # Every row is as wide as the header — a ragged grid loses a column silently
    # on the way to the file.
    for name in list(bundles) + ['reliability_taylor', 'reliability_mc',
                                 'sensitivity_curve']:
        headers = views.get(name + '_headers') or []
        for i, row in enumerate(views.get(name + '_rows') or []):
            if len(row) != len(headers):
                failures.append(f"{name}: row {i} has {len(row)} cell(s) under "
                                f"{len(headers)} header(s)")
                break


# ---------------------------------------------------------------------------
# C. Save CSV… writes the grid, and nothing else
# ---------------------------------------------------------------------------

def check_the_csv_round_trips(failures):
    """The file that leaves the program is the grid that was on screen.

    Save CSV… is how these numbers get quoted, so a file that drops a column,
    re-orders rows or re-formats a value would be read in place of the run. Every
    mode's file is written, read back and compared cell for cell.
    """
    views = _views()
    names = list(_SHARED['bundles']) + ['reliability_taylor', 'reliability_mc',
                                        'sensitivity_curve']
    with tempfile.TemporaryDirectory() as tmp:
        for name in names:
            view = views.get(name)
            if view is None:
                continue
            headers = views[name + '_headers']
            rows = views[name + '_rows']
            path = os.path.join(tmp, f"{name}.csv")
            # The view's own writer, the one Save CSV… calls once the dialog has
            # named a file.
            view.table_view.set_table(headers, rows, path)
            view.table_view.write_csv(path)
            with open(path, newline='', encoding='utf-8') as fh:
                back = list(csv.reader(fh))
            if not back:
                failures.append(f"{name}: the CSV is empty")
                continue
            if back[0] != headers:
                failures.append(f"{name}: CSV header {back[0]} != grid header "
                                f"{headers}")
            if back[1:] != rows:
                bad = next((i for i, (a, b) in enumerate(zip(back[1:], rows))
                            if a != b), None)
                failures.append(
                    f"{name}: CSV rows differ from the grid"
                    + (f" (first at row {bad}: {back[1:][bad]} != {rows[bad]})"
                       if bad is not None
                       else f" ({len(back) - 1} row(s) written, {len(rows)} shown)"))


def check_the_offered_csv_name(failures):
    """Save CSV… opens on ``<model>_<mode>.csv`` beside the project, so a sweep's
    numbers land next to the file they were computed from rather than in whatever
    directory the dialog last visited."""
    win = _window()
    win.doc.path = os.path.join(os.sep, "models", "johnson.xlsx")
    for stem, expected in (("fs_vs_time", "johnson_fs_vs_time.csv"),
                           ("tornado", "johnson_tornado.csv"),
                           ("reliability_mc", "johnson_reliability_mc.csv")):
        got = win._sweep_csv_path(stem)
        if os.path.basename(got) != expected:
            failures.append(f"{stem}: offered name {os.path.basename(got)!r}, "
                            f"expected {expected!r}")
        if os.path.dirname(got) != os.path.join(os.sep, "models"):
            failures.append(f"{stem}: offered {got!r}, not beside the project")
    win.doc.path = None
    if win._sweep_csv_path("tornado") != "tornado.csv":
        failures.append("an unsaved project should offer a bare name, got "
                        f"{win._sweep_csv_path('tornado')!r}")


# ---------------------------------------------------------------------------
# D. the grid says what the figure beside it says
# ---------------------------------------------------------------------------

def check_the_march_carries_its_faces(failures):
    """COMBO-3 Part 1's dam hands over from one slope to the other, and the figure
    colors each instant by the face its critical circle came out on. The table
    states the same face in words — read from the same classifier — so a reader
    quoting the grid places the handover where the plot does."""
    views = _views()
    headers = views.get('fs_vs_time_headers') or []
    rows = views.get('fs_vs_time_rows') or []
    if 'face' not in headers:
        failures.append(f"the march's table has no face column: {headers}")
        return
    j = headers.index('face')
    for row, published in zip(rows, COMBO03_CURVE):
        if row[j] != published[5]:
            failures.append(f"t = {published[0]:g}: table says face {row[j]!r}, "
                            f"the published curve says {published[5]!r}")
    # …and the circle it was found on, beside it.
    for col in ('Xo', 'Yo', 'R'):
        if col not in headers:
            failures.append(f"the march's table has no {col} column: {headers}")


def check_a_failed_instant_carries_its_reason(failures):
    """An instant that produced nothing is a row saying why, never a gap. A blank
    row reads as a frame that does not exist rather than as a solve that did not
    land, and the two call for opposite things from the reader."""
    views = _views()
    headers = views.get('drawdown_vs_time_headers') or []
    rows = views.get('drawdown_vs_time_rows') or []
    if 'note' not in headers:
        failures.append(f"the drawdown table has no note column: {headers}")
        return
    j_note, j_fs = headers.index('note'), headers.index('FS')
    first = rows[0] if rows else []
    if not first:
        failures.append("the drawdown table has no rows")
        return
    if first[j_fs] != sweep_table.DASH:
        failures.append(f"the failed instant reports FS {first[j_fs]!r} rather "
                        "than no number")
    if not first[j_note].strip():
        failures.append("the failed instant's row carries no reason")


def check_a_drawdown_row_shows_its_stages(failures):
    """A rapid-drawdown row reports the lower of stages 2 and 3, so the stage
    values and the stage that governed sit beside it: a reported number with no
    stages behind it hides which analysis produced it, which is the reading of a
    drawdown curve."""
    views = _views()
    headers = views.get('drawdown_vs_time_headers') or []
    for col in ('stage 1', 'stage 2', 'stage 3', 'governs'):
        if col not in headers:
            failures.append(f"the drawdown table has no {col!r} column: {headers}")
    rows = views.get('drawdown_vs_time_rows') or []
    if 'stage 3' in headers and len(rows) > 1:
        j = headers.index('stage 3')
        # Stage 3 is run only where stage 2 did not already govern; the row says so
        # rather than leaving a blank the reader has to account for.
        if rows[1][j] != 'not required':
            failures.append("a solved drawdown row whose stage 3 never ran should "
                            f"say 'not required', got {rows[1][j]!r}")
    # A single-stage march carries none of those columns — they would be four
    # columns of dashes claiming an analysis the run did not make.
    plain = views.get('fs_vs_time_headers') or []
    for col in ('stage 1', 'stage 2', 'stage 3', 'governs'):
        if col in plain:
            failures.append(f"a single-stage march's table carries {col!r}")


def check_the_tornado_table_is_its_bars(failures):
    """The tornado's table is the bars as numbers: each parameter's low and high
    bound with the factor of safety at each, widest swing first — the order the
    bars are stacked in. A table sorted some other way than the plot is a second
    reading of the same run."""
    views = _views()
    headers = views.get('tornado_headers') or []
    rows = views.get('tornado_rows') or []
    for col in ('parameter', 'low value', 'high value', 'swing'):
        if col not in headers:
            failures.append(f"the tornado table has no {col!r} column: {headers}")
            return
    j = headers.index('swing')
    swings = [float(r[j]) for r in rows]
    if swings != sorted(swings, reverse=True):
        failures.append(f"the tornado table's swings are not widest-first: "
                        f"{swings}")
    j_lo, j_hi = headers.index('FS at low'), headers.index('FS at high')
    sweeps = _sweeps()
    for row in rows:
        df = sweeps.get(row[0])
        if df is None:
            failures.append(f"the tornado table names {row[0]!r}, which was not "
                            "swept")
            continue
        ok = df.loc[~df['is_base'] & df['success']].sort_values('value')
        for cell, want in ((row[j_lo], ok['fs'].iloc[0]),
                           (row[j_hi], ok['fs'].iloc[-1])):
            if abs(float(cell) - float(want)) > 5e-5:
                failures.append(f"{row[0]}: table says {cell}, the sweep says "
                                f"{want:.4f}")


def run():
    failures = []
    if not COMBO03.exists():
        return [f"missing {COMBO03}"]
    for fn in (check_every_view_has_a_table_sub_tab,
               check_the_row_count_is_the_results_own,
               check_the_csv_round_trips,
               check_the_offered_csv_name,
               check_the_march_carries_its_faces,
               check_a_failed_instant_carries_its_reason,
               check_a_drawdown_row_shows_its_stages,
               check_the_tornado_table_is_its_bars):
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
    print("\nEvery Parametric result view carries its numbers beside its plot, "
          "row for row, and Save CSV… writes the grid exactly as shown.")


if __name__ == '__main__':
    main()
