"""The numbers behind a Parametric result plot, as a grid — and the CSV of it.

Every Parametric result view (a design or back-analysis curve, the five sweep
plots, a reliability summary, a factor-of-safety-versus-time march) is a picture
of a table the run already computed. A picture is read for its shape; the table is
what gets quoted, checked against a hand calculation and pasted into a report. So
each of those views carries two sub-tabs — **Plot**, which stays the default, and
**Table** — and the Table tab is the same widget everywhere, so a reader who has
used one has used all of them.

The module is two halves:

  * :class:`SweepTableView`, the grid plus its **Save CSV…** button, and
    :class:`SweepResultView`, the plot/table sub-tab container that wraps a result
    canvas and forwards to it (so the main window's tab, display-panel and
    clear machinery is unchanged);
  * the builders — one per result kind — that turn a runner's bundle into
    ``(headers, rows, stem)``, where ``rows`` are already-formatted strings. They
    are plain functions over the bundles, so what a tab shows can be checked
    without a screen.

Each builder reads the SAME numbers its plot reads (the tornado's low/high pair
comes from the grouping ``plot_tornado`` uses, the march's rows are the ones
``fs_vs_time`` already printed), so the grid can never disagree with the figure
beside it.
"""

from __future__ import annotations

import csv
import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView, QFileDialog, QHBoxLayout, QHeaderView, QMessageBox,
    QTableWidget, QTableWidgetItem, QTabWidget, QToolButton, QVBoxLayout, QWidget,
)

#: What a missing number reads as, in the grid and in the CSV — the same em dash
#: the engine's own console tables use, so the two match cell for cell.
DASH = "—"


# ---------------------------------------------------------------------------
# formatting
# ---------------------------------------------------------------------------

def _f(v, fmt="{:.4f}"):
    """``v`` formatted, or the em dash where there is no number."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return DASH
    if x != x or x in (float("inf"), float("-inf")):
        return DASH
    return fmt.format(x)


def _g(v):
    return _f(v, "{:g}")


def _n(v):
    return _f(v, "{:.2f}")


def _sig(v):
    return _f(v, "{:.4g}")


# ---------------------------------------------------------------------------
# builders — one per result kind, each returning (headers, rows, stem)
# ---------------------------------------------------------------------------

def fs_vs_time_table(bundle, slope_data=None):
    """The march's own record: one row per evaluated instant.

    Reads ``bundle['table']`` — the row list ``fs_vs_time`` built and printed — so
    the grid and the console table carry identical numbers. The face each instant's
    critical circle came out on is added beside the circle, because it is the one
    thing the figure shows that the printed table does not: a curve that hands over
    from one slope to the other is two mechanisms in sequence, and the row says
    which.
    """
    rapid = bool(bundle.get("rapid"))
    rows_in = list(bundle.get("table") or [])
    if not rows_in and bundle.get("df") is not None:
        # A bundle that carries only its DataFrame (an older stored result) gets the
        # rows built the same way the run builds them, rather than a second reading
        # of the frame that could drift from it.
        from xslope.sensitivity import _fs_vs_time_table
        rows_in, _text = _fs_vs_time_table(bundle["df"], slope_data, rapid)
    unit = str((slope_data or {}).get("time_unit") or "").strip()
    faces = _faces_by_row(rows_in, slope_data)

    multi = len({r.get("method") for r in rows_in}) > 1
    has_circle = any(r.get("Xo") is not None and r.get("R") is not None
                     for r in rows_in)
    has_face = any(faces)
    has_note = any(not r.get("success") for r in rows_in)

    headers = [f"t ({unit})" if unit else "t"]
    if multi:
        headers.append("method")
    if rapid:
        headers += ["stage 1", "stage 2", "stage 3"]
    headers.append("FS")
    if rapid:
        headers.append("governs")
    if has_circle:
        headers += ["Xo", "Yo", "R"]
    if has_face:
        headers.append("face")
    if has_note:
        headers.append("note")

    rows = []
    for r, face in zip(rows_in, faces):
        cells = [_g(r.get("time"))]
        if multi:
            cells.append(str(r.get("method") or ""))
        if rapid:
            cells.append(_f(r.get("stage1_FS")))
            cells.append(_f(r.get("stage2_FS")))
            # Stage 3 is run only where stage 2 did not already govern; saying so
            # is the answer, not a blank the reader has to account for.
            cells.append("not required" if (r.get("success")
                                            and not r.get("stage3_run"))
                         else _f(r.get("stage3_FS")))
        cells.append(_f(r.get("fs")))
        if rapid:
            cells.append(str(r["governs"]) if r.get("governs") else DASH)
        if has_circle:
            cells += [_n(r.get("Xo")), _n(r.get("Yo")), _n(r.get("R"))]
        if has_face:
            cells.append(face or DASH)
        if has_note:
            cells.append("" if r.get("success") else (r.get("msg") or "no result"))
        rows.append(cells)
    stem = "drawdown_vs_time" if rapid else "fs_vs_time"
    return headers, rows, stem


def _faces_by_row(rows, slope_data):
    """The embankment face each row's critical circle sits on, or None where the
    geometry does not say. Uses the plotter's own classifier, so the column and the
    figure's marker colors can only ever agree."""
    if not slope_data:
        return [None] * len(rows)
    try:
        from xslope.plot import _circle_face, _crest_span
        span = _crest_span(slope_data)
    except Exception:                                          # noqa: BLE001
        return [None] * len(rows)
    if span is None:
        return [None] * len(rows)
    out = []
    for r in rows:
        if not r.get("success"):
            out.append(None)
            continue
        try:
            out.append(_circle_face(slope_data, span, r.get("Xo"), r.get("Yo"),
                                    r.get("R")))
        except Exception:                                      # noqa: BLE001
            out.append(None)
    return out


def design_table(bundle):
    """A design or back-analysis sweep: the parameter value at every point and the
    output it produced, in the order the sweep solved them."""
    df = bundle.get("df")
    param = str(bundle.get("param") or "value")
    output = str(bundle.get("output") or "FS")
    return _sweep_df_table(df, param, output,
                           stem=("back_analysis"
                                 if bundle.get("study") == "back_analysis"
                                 else "design"))


def curve_table(param, df, output="FS"):
    """One parameter's factor-of-safety-versus-value curve — the tornado bar's
    click-through view."""
    return _sweep_df_table(df, str(param), str(output), stem="sensitivity_curve")


def _sweep_df_table(df, param, output, stem):
    """The shared shape of the two views drawn from a single ``sensitivity()``
    DataFrame: value, output, the circle it was found on, and the reason where a
    point produced nothing."""
    headers = [param, output]
    if df is None or not len(df):
        return headers, [], stem
    has_circle = bool({"Xo", "R"} <= set(df.columns)) and bool(
        df["Xo"].notna().any() and df["R"].notna().any())
    multi = "method" in df.columns and int(df["method"].nunique()) > 1
    has_note = "success" in df.columns and bool((~df["success"].astype(bool)).any())
    if multi:
        headers.insert(1, "method")
    if has_circle:
        headers += ["Xo", "Yo", "R"]
    if has_note:
        headers.append("note")

    rows = []
    for _i, r in df.iterrows():
        ok = bool(r.get("success", True))
        cells = [_sig(r.get("value"))]
        if multi:
            cells.append(str(r.get("method") or ""))
        cells.append(_f(r.get("fs")) if ok else DASH)
        if has_circle:
            cells += [_n(r.get("Xo")), _n(r.get("Yo")), _n(r.get("R"))]
        if has_note:
            cells.append("" if ok else (str(r.get("msg") or "") or "no result"))
        rows.append(cells)
    return headers, rows, stem


def tornado_table(bundle):
    """The tornado's bars as numbers: each parameter's low and high bound, the
    output at each, and the swing between them — widest first, the order the bars
    are stacked in."""
    result = bundle.get("tornado") or {}
    df = result.get("df")
    output = str(result.get("output") or "FS")
    headers = ["parameter", "low value", "high value",
               f"{output} at low", f"{output} at high", "swing"]
    bars = []
    if df is not None and len(df):
        for param, g in df.loc[~df["is_base"]].groupby("param"):
            ok = g.loc[g["success"]].sort_values("value")
            if ok.empty:
                continue
            bars.append((str(param), ok["value"].iloc[0], ok["value"].iloc[-1],
                         ok["fs"].iloc[0], ok["fs"].iloc[-1],
                         abs(float(ok["fs"].iloc[-1]) - float(ok["fs"].iloc[0]))))
    bars.sort(key=lambda b: b[5], reverse=True)
    rows = [[b[0], _sig(b[1]), _sig(b[2]), _f(b[3]), _f(b[4]), _f(b[5])]
            for b in bars]
    return headers, rows, "tornado"


def spider_table(bundle):
    """Every point of every parameter's curve — the spider plot is one line per
    parameter over a shared relative-change axis, so its table is one row per
    solve."""
    sweeps = bundle.get("sweeps") or {}
    headers = ["parameter", "value", "change (%)", "FS", "note"]
    rows = []
    for param, df in sweeps.items():
        if df is None:
            continue
        for _i, r in df.iterrows():
            ok = bool(r.get("success", True))
            rel = r.get("rel")
            rows.append([
                str(param), _sig(r.get("value")),
                _f(None if rel is None else float(rel) * 100.0, "{:.1f}"),
                _f(r.get("fs")) if ok else DASH,
                "" if ok else (str(r.get("msg") or "") or "no result"),
            ])
    return headers, rows, "spider"


def scaled_table(bundle):
    """The scaled-sensitivity bars: the derivative at the base case, in each of the
    scalings that make different-unit parameters comparable."""
    result = bundle.get("scaled") or {}
    bars = list(result.get("bars") or [])
    has_sigma = bool(result.get("has_sigma"))
    headers = ["parameter", "base value", "dF/dp", "elasticity", "ΔF per 1%"]
    if has_sigma:
        headers += ["σ", "ΔF per σ"]
    rows = []
    for b in bars:
        cells = [str(b.get("param") or ""), _sig(b.get("base_value")),
                 _sig(b.get("dfdp")), _f(b.get("elasticity")),
                 _f(b.get("per_1pct"))]
        if has_sigma:
            cells += [_sig(b.get("sigma")), _f(b.get("per_sigma"))]
        rows.append(cells)
    return headers, rows, "scaled_sensitivity"


def variance_table(bundle):
    """The variance Pareto: each uncertain parameter's share of Var(FS), and the
    running total the Pareto line traces."""
    result = bundle.get("variance") or {}
    headers = ["parameter", "σ", "ΔF", "variance", "% of Var(FS)", "cumulative %"]
    rows = [[str(b.get("label") or b.get("param") or ""), _sig(b.get("sigma")),
             _f(b.get("delta_F")), _sig(b.get("variance")),
             _f(b.get("pct"), "{:.2f}"), _f(b.get("cumulative"), "{:.2f}")]
            for b in (result.get("bars") or [])]
    return headers, rows, "variance_contribution"


def rank_table(bundle):
    """The Monte Carlo rank correlations: how strongly each sampled input drives
    the factor of safety across the whole sampled distribution."""
    result = bundle.get("rank") or {}
    headers = ["parameter", "σ", "ρ (Spearman)"]
    rows = [[str(b.get("label") or b.get("param") or ""), _sig(b.get("sigma")),
             _f(b.get("rho"), "{:.3f}")]
            for b in (result.get("bars") or [])]
    return headers, rows, "rank_correlation"


#: The reported quantities of a sampling reliability run, in reading order.
_MC_SUMMARY = [
    ("mean FS", "mean_FS", "{:.4f}"),
    ("σ_F", "sigma_F", "{:.4f}"),
    ("COV_F", "COV_F", "{:.4f}"),
    ("β (normal)", "beta_normal", "{:.4f}"),
    ("β (lognormal)", "beta_ln", "{:.4f}"),
    ("Pf empirical (%)", "pf_empirical", "{:.4f}"),
    ("Pf normal (%)", "pf_normal", "{:.4f}"),
    ("Pf lognormal (%)", "pf_lognormal", "{:.4f}"),
]


def reliability_table(rel, engine=None):
    """A reliability run's own table.

    The Taylor series works parameter by parameter, so its table is the
    per-parameter one: the most-likely value, the one-sigma pair either side of it,
    the factor of safety at each, and the swing ΔF that the variance is built from.
    A sampling run has no such pair — it has a distribution — so its table is the
    reported statistics, the numbers the histogram is a picture of.
    """
    rel = rel or {}
    if engine in ("mc", "rs"):
        rows = []
        for label, key, fmt in _MC_SUMMARY:
            v = rel.get(key)
            if v is None:
                continue
            rows.append([label, _f(100.0 * float(v) if label.startswith("Pf") else v,
                                   fmt)])
        for label, key in (("valid samples", "n_valid"),
                           ("samples requested", "n_samples"),
                           ("distribution", "distribution")):
            if rel.get(key) is not None:
                rows.append([label, str(rel[key])])
        stem = "reliability_rs" if engine == "rs" else "reliability_mc"
        return ["quantity", "value"], rows, stem

    headers = ["parameter", "MLV", "σ", "MLV+σ", "MLV−σ", "F+", "F−", "ΔF"]
    rows = []
    for p in (rel.get("param_info") or []):
        mlv, std = p.get("mlv"), p.get("std")
        try:
            up = float(mlv) + float(std)
            dn = float(mlv) - float(std)
        except (TypeError, ValueError):
            up = dn = None
        rows.append([
            f"{p.get('material_name', '')} · {p.get('param', '')}".strip(" ·"),
            _sig(mlv), _sig(std), _sig(up), _sig(dn),
            _f(p.get("F_plus")), _f(p.get("F_minus")), _f(p.get("delta_F")),
        ])
    return headers, rows, "reliability_taylor"


def parametric_table(bundle, slope_data=None):
    """``(headers, rows, stem)`` for any Parametric result bundle — the one entry
    point the main window and the check both use."""
    bundle = bundle or {}
    kind = bundle.get("kind")
    if kind == "fs_vs_time":
        return fs_vs_time_table(bundle, slope_data)
    if kind == "design":
        return design_table(bundle)
    plot_type = bundle.get("plot_type", "tornado")
    return {"scaled": scaled_table, "spider": spider_table,
            "variance": variance_table,
            "rank": rank_table}.get(plot_type, tornado_table)(bundle)


# ---------------------------------------------------------------------------
# widgets
# ---------------------------------------------------------------------------

class SweepTableView(QWidget):
    """A result's numbers as a read-only grid, with **Save CSV…** beside them.

    The grid holds already-formatted strings — the same text the engine's console
    tables print — so what is on screen, what lands in the CSV and what the run
    reported are one set of numbers. Rows keep the order the run produced them
    (sorting is off): a march is read in time and a Pareto in rank, and a
    re-ordered grid would no longer be the table the plot is drawn from.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._headers = []
        self._rows = []
        self._default_path = ""

        self.table = QTableWidget(self)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectItems)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        # Each column sizes to its own numbers, and the last one keeps that size
        # rather than stretching: a text column (a face, a reason) would otherwise
        # take the whole spare width of the window and push the numbers together.
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(False)

        bar = QHBoxLayout()
        bar.addStretch(1)
        self._save_btn = QToolButton()
        self._save_btn.setText("Save CSV…")
        self._save_btn.setToolTip("Write these numbers to a comma-separated file")
        self._save_btn.clicked.connect(self.save_csv)
        bar.addWidget(self._save_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addLayout(bar)
        layout.addWidget(self.table, 1)

    # --- content ---------------------------------------------------------
    def set_table(self, headers, rows, default_path=""):
        """Fill the grid. ``default_path`` is what Save CSV… offers first."""
        self._headers = [str(h) for h in headers]
        self._rows = [[str(c) for c in row] for row in rows]
        self._default_path = default_path or ""
        self.table.clear()
        self.table.setColumnCount(len(self._headers))
        self.table.setRowCount(len(self._rows))
        self.table.setHorizontalHeaderLabels(self._headers)
        for i, row in enumerate(self._rows):
            for j, cell in enumerate(row):
                item = QTableWidgetItem(cell)
                # Numbers right, words left — a column of factors of safety is read
                # down its decimal point.
                if _is_number(cell):
                    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.table.setItem(i, j, item)
        # A freshly filled grid has nothing selected — Qt would otherwise leave the
        # first cell looking picked, which reads as a value the run singled out.
        self.table.setCurrentCell(-1, -1)

    def headers(self):
        return list(self._headers)

    def rows(self):
        return [list(r) for r in self._rows]

    def row_count(self):
        return len(self._rows)

    # --- export ----------------------------------------------------------
    def save_csv(self, _checked=False):
        path, _sel = QFileDialog.getSaveFileName(
            self, "Save table", self._default_path, "CSV file (*.csv)")
        if not path:
            return
        if not os.path.splitext(path)[1]:
            path += ".csv"
        try:
            self.write_csv(path)
        except Exception:                                      # noqa: BLE001
            import traceback
            traceback.print_exc()
            QMessageBox.warning(self, "Save table failed",
                                "Could not write the file — see the Log pane.")
            return
        win = self.window()
        if hasattr(win, "statusBar"):
            win.statusBar().showMessage(f"Saved {os.path.basename(path)}")

    def write_csv(self, path):
        """Write the grid to ``path`` exactly as shown — header row, then rows."""
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(self._headers)
            writer.writerows(self._rows)


def _is_number(text):
    try:
        float(str(text).replace(",", ""))
    except (TypeError, ValueError):
        return False
    return True


class SweepResultView(QWidget):
    """A Parametric result view: its plot and its table as two sub-tabs.

    Wraps a result canvas without changing it. **Plot** is the default sub-tab —
    the figure is what a sweep is run to see — and **Table** sits beside it. Every
    attribute the main window reaches for on a result view (``render_tornado``,
    ``ensure_fitted``, ``picked``, ``_main_axes`` …) is forwarded to the wrapped
    canvas, so the tab, display-panel and clear machinery needs no special case:
    the view is added to ``view_tabs`` and rendered exactly as the bare canvas was.
    """

    def __init__(self, canvas, parent=None):
        super().__init__(parent)
        self._canvas = canvas
        canvas.setParent(self)
        self.table_view = SweepTableView(self)
        self.tabs = QTabWidget(self)
        # A document-style sub-strip, visually subordinate to the result tabs above
        # it; the plot is tab 0, so a fresh result opens on its figure.
        self.tabs.setDocumentMode(True)
        self.tabs.addTab(canvas, "Plot")
        self.tabs.addTab(self.table_view, "Table")
        self.tabs.currentChanged.connect(self._on_sub_tab)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.tabs)

    @property
    def canvas(self):
        return self._canvas

    def set_table(self, headers, rows, default_path=""):
        self.table_view.set_table(headers, rows, default_path)

    def show_plot(self):
        self.tabs.setCurrentIndex(0)

    def _on_sub_tab(self, index):
        # The canvas rasterizes to the viewport size, so coming back from the table
        # (where it had no size) needs the same re-fit a result tab gets.
        if index == 0:
            self._canvas.ensure_fitted()

    def __getattr__(self, name):
        # Only reached when normal lookup failed, so it never shadows this class's
        # own attributes. Guarded on _canvas: during __init__ (before it is set) an
        # attribute miss must stay a miss rather than recursing.
        try:
            canvas = object.__getattribute__(self, "_canvas")
        except AttributeError:
            raise AttributeError(name) from None
        return getattr(canvas, name)
